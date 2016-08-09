#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <time.h>
#include <fstream>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <deque>
#include <utility>
#include <cmath>
#include <time.h>
#include <pthread.h>

using namespace std;

//if BOUNDED_HISTORY and HISTORY_BOUND are defined then user's histories will only contain the last HISTORY_BOUND products bought
//these symbols are currently defined in g++ command line.
//#define BOUNDED_HISTORY
//#define HISTORY_BOUND 100

//Model parameters
double alpha; //probability of following a recommendation
double beta; //beta (see the paper). beta=1000 means "infinity"
int num_products; //number of types of products
double advantage; //initial "advantage" of the first product users' personal preferences (see the paper).

//Social graph
char* graph_filename; //Filename of the input graph
long N=-1; //number of vertices
long* in_degrees = NULL;
long** in_neighbors = NULL;

char* backgrounds_filename; //Filename where backgrounds are saved
//Matrix of personal distributions.
//backgrounds[i][j] is the probability that user i buys product j when he does not follow a recommendation.
double** backgrounds = NULL; 
long* owned_products_no = NULL; //owned_products_no[y] is the overall number of products in the history of user i.
long** owned_products = NULL; //matrix of bought products. owned_products[i][j] is the number of products of type j bought by user i.

#ifdef BOUNDED_HISTORY
	deque<int>** history; //history of each user. history[i] is a dequeue containing all the products in the history of users i, in order of purchase.
	long* bought_products_no; //bought_products_no[i] is the overall number of products bought by user i during the whole simulation.
	long** bought_products; //bought_products[i][j] is the overall number of products of type j bought by user i during the whole simulation.
#endif

//Number of threads to use for the simulation
int num_threads;

//Random number generators (one main rng + one for each thread)
int main_seed;
gsl_rng *main_rng;
gsl_rng **rngs = NULL;

//Synchronization objects
pthread_mutex_t* mutex_vertex; //One mutex for each vertex/user. Access to user's purchases is in mutual exclusion

pthread_mutex_t mutex_iteration; //Mutex for updating iteration count
pthread_mutex_t* mutex_waiting_writer; //Mutexes and condition variable to sync thread with stats writer
pthread_cond_t* conds_waiting_writer_set;
pthread_cond_t* conds_waiting_writer_unset;
bool* waiting_writer;

bool shutdown=false; //Has the simulation ended? Setting this to true causes the threads to terminate
long long iteration=0; //Number of steps of the simulation done so far = number of purchases
long long next_write=0; //Iteration count of next stats write

long start_time; //Timestamp of program startup

static void* simulate(void *extra);
void write_stats();


// Samples a product from a probility distribution
// Input: a random number generator and a vector of num_products elements representing
// a probability distribution on [0, num_products-1]
// Output: a product p, such that  0 < p < num_products-1 and P(p=k) = probs[k]
int sample_product(gsl_rng* rng, double* probs)
{
	double unif = gsl_rng_uniform(rng);
	int product = 0;
	while(unif > probs[product] && product<num_products-1)
		unif -= probs[product++];

	return product;
}

//Writes a log message to stderr. Format and ... have the same meaning as in printf.
void log(const char *format, ...)
{
    va_list args;
    va_start(args, format);

    fprintf(stderr, "[%ld] ", (long)(time(NULL)-start_time));
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");

    va_end(args);
}

//Parses the command line arguments and exits on failure.
void parse_args(int argc, char** argv)
{
	if(argc<9)
	{
		fprintf(stderr, "Usage: %s graph alpha beta num_threads seed background_file num_products advantage\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	graph_filename = argv[1];
	alpha = atof(argv[2]);
	beta =  atof(argv[3]);
	num_products = atoi(argv[7]);

	num_threads = atoi(argv[4]);
	main_seed = atoi(argv[5]);
	advantage = atof(argv[8]);

	backgrounds_filename = argv[6];
}

void init()
{
	//Increase the stdout buffer (improves performance)
	static char stdoutbuf[5 * 1024 * 1024]; //5MiB
	setvbuf( stdout , stdoutbuf , _IOFBF , sizeof(stdoutbuf) );

	start_time = (long)time(NULL);
}

void load_graph(char* filename)
{
	FILE* input = fopen(filename, "r");

	log("Counting vertices...");

	char buf[2048];
	long u, v;
	while (fgets(buf, 2048, input) != NULL)
	{
		if (buf[0] == '#')
			continue;

		sscanf(buf, " %ld %ld ", &u, &v);

		N = (N > u) ? N : u;
		N = (N > v) ? N : v;
	}
	rewind(input);

	log("Computing degrees.");

	in_degrees = new long[N + 1];
	in_neighbors = new long *[N + 1];
	long *in_cnts = new long[N + 1];

	for (long i = 0; i <= N; i++)
		in_cnts[i] = in_degrees[i] = 0;

	while (fgets(buf, 2048, input) != NULL)
	{
		if (buf[0] == '#')
			continue;

		sscanf(buf, " %ld %ld ", &u, &v);

		in_degrees[u]++; //reversing edges
	}
	rewind(input);

	for (long i = 0; i <= N; i++)
		in_neighbors[i] = new long[in_degrees[i]];

	log("Loading graph.");

	while (fgets(buf, 2048, input) != NULL)
	{
		if (buf[0] == '#')
			continue;

		sscanf(buf, " %ld %ld ", &u, &v);

		in_neighbors[u][in_cnts[u]++] = v; //reversing edges
	}
	fclose(input);
	delete[] in_cnts;
}

//Creates and seeds the random number generators
void setup_rngs()
{
	gsl_rng_env_setup();
	const gsl_rng_type *T = gsl_rng_default;
	main_rng = gsl_rng_alloc(T);
	gsl_rng_set(main_rng, main_seed+42);

	rngs = new gsl_rng*[num_threads];

	for(int k=0; k<num_threads; k++)
	{
		rngs[k] = gsl_rng_alloc(T);
		gsl_rng_set(rngs[k], main_seed+64+100*k);
	}
}

//Compares two double. Helper function for setup_backgrounds(...)
int cmp(const void* a, const void* b)
{
	double A = *((double*)a);
	double B = *((double*)b);

	if(A<B)
		return -1;
	if(A>B)
		return 1;
	return 0;
}

//Sets the backgrounds distribution for each users by slicing the unit interval in prodcuts-1+advantage pieces.
//The first product gets "advantage" pieces. The other produdcts get one piece each.
void setup_backgrounds(gsl_rng* rng)
{
	backgrounds = new double*[N + 1];

	int advantage100 = (int) (round(advantage * 100));
	int nsplits = (num_products - 1) * 100 + advantage100;
	double splits[nsplits];

	for (long i = 0; i <= N; i++)
	{
		splits[nsplits - 1] = 1;
		for (int j = 0; j < nsplits - 1; j++)
			splits[j] = gsl_rng_uniform(rng);

		qsort(splits, nsplits, sizeof(double), cmp);

		backgrounds[i] = new double[num_products];
		backgrounds[i][0] = splits[advantage100 - 1];
		for (int j = 1; j < num_products; j++)
			backgrounds[i][j] = splits[j * 100 + advantage100 - 1]	- splits[(j - 1) * 100 + advantage100 - 1];
	}
}

//Writes the users' background distribution to filename
void write_backgrounds(char* filename)
{
	FILE* bgfile = fopen(filename, "w");
	for (long i = 0; i <= N; i++)
	{
		for (int j = 0; j < num_products; j++)
			fprintf(bgfile, "%.10f ", backgrounds[i][j]);

		fprintf(bgfile, "\n");
	}
	fclose(bgfile);
}

void setup_products(gsl_rng* rng)
{

#ifdef BOUNDED_HISTORY
	history = new deque<int>*[N+1];
	bought_products_no = new long[N+1];
	bought_products = new long*[N+1];

	for(int i=0; i<=N; i++)
	{
		bought_products[i] = new long[num_products];
		history[i] = new deque<int>();
	}
#endif

	//Setup products
	owned_products_no = new long[N+1];
	owned_products = new long*[N+1];
	for(int i=0; i<=N; i++)
	{
		owned_products_no[i] = 1;
		owned_products[i]=new long[num_products];

		//Each user initially owns one prodcuts, chosen according to his personal preference
		int p = sample_product(rng, backgrounds[i]);
		for(int j=0; j<num_products; j++)
		{
			owned_products[i][j]= (j==p)?1:0;
		}

#ifdef BOUNDED_HISTORY
		bought_products_no[i] = owned_products_no[i];

		for(int j=0; j<num_products; j++)
		{
			bought_products[i][j] = owned_products[i][j];
			for(int k=0; k<bought_products[i][j]; k++)
				history[i]->push_back(j);
		}
#endif
	}
}

//Creates the mutex and condition variables used to synchronize threads
void setup_thread_sync_objects()
{
	mutex_vertex = new pthread_mutex_t[N + 1];
	for (int i = 0; i <= N; i++)
		pthread_mutex_init(&mutex_vertex[i], NULL);

	mutex_waiting_writer = new pthread_mutex_t[num_threads];
	conds_waiting_writer_set = new pthread_cond_t[num_threads];
	conds_waiting_writer_unset = new pthread_cond_t[num_threads];
	waiting_writer = new bool[num_threads];

	for (int k = 0; k < num_threads; k++)
	{
		pthread_mutex_init(&mutex_waiting_writer[k], NULL);
		pthread_cond_init(&conds_waiting_writer_set[k], NULL);
		pthread_cond_init(&conds_waiting_writer_unset[k], NULL);
		waiting_writer[k] = false;
	}
}

int main(int argc, char** argv)
{
	init();
	parse_args(argc, argv);

	log("Alpha is %f. Beta is %f", alpha, beta);
	log("Main seed is: %d", main_seed);
	log("Using %d thread(s)", num_threads);

	setup_rngs();

	load_graph(graph_filename);

	setup_backgrounds(main_rng);
	write_backgrounds(backgrounds_filename);

	setup_products(main_rng);

	setup_thread_sync_objects();

	log("Starting simulation");

	//Start threads
	shutdown=false;
	int writes=0;
	next_write = (writes+1)*100*N;
	//write_stats();

	pthread_t threads[num_threads];
	for(int k=0; k<num_threads; k++)
	{
		int* tid = new int;
		*tid = k;
		pthread_create(&threads[k], NULL, &simulate, tid);
	}


	while(!shutdown)
	{
		for(int k=0; k<num_threads; k++)
		{
			pthread_mutex_lock(&mutex_waiting_writer[k]);
			while( !waiting_writer[k] )
				pthread_cond_wait(&conds_waiting_writer_set[k], &mutex_waiting_writer[k]);
			pthread_mutex_unlock(&mutex_waiting_writer[k]);
		}


		write_stats();

		writes++;
		next_write = (writes+1)*100*N;

		if(writes==30)
			shutdown=true;

		for(int k=0; k<num_threads; k++)
		{
			pthread_mutex_lock(&mutex_waiting_writer[k]);

			waiting_writer[k] = false;
			pthread_cond_signal(&conds_waiting_writer_unset[k]);

			pthread_mutex_unlock(&mutex_waiting_writer[k]);
		}

	}


	for(int k=0; k<num_threads; k++)
	{
		pthread_join(threads[k], NULL);
	}

	fflush(stdout);
	fflush(stderr);

	return EXIT_SUCCESS;
}

//Writes (to stdout) the current state of the simulation, i.e., the number of products of each type owned by each user.
//If the simulation is using a bounded history, the overall number of products of each type ever bought by each user is also reported.
void write_stats()
{
	log("Iteration %lld (%f%%)", iteration, ((double)iteration)/(N*100*30)*100 );

	printf("#Iteration %lld\n", iteration);
	for(long i=0; i<=N; i++)
	{
		printf("%ld ", owned_products_no[i]);
		for(int j=0; j<num_products; j++)
			printf("%ld ", owned_products[i][j]);

#ifdef BOUNDED_HISTORY
		printf("%ld ", bought_products_no[i]);
		for(int j=0; j<num_products; j++)
			printf("%ld ", bought_products[i][j]);
#endif

		printf("\n");
	}
	fflush(stdout);
	fflush(stderr);
}

//Modifies a vector of probabilities "prob" so that the probability of an element i is proportional to prob[i]^beta
//beta = 1000 means "infinity": the item of maximizing probability will get a new probability of 1
//if there are k maxima, each of them will get a probability of 1/k
void apply_beta(double* prob)
{
	if (beta != 1 && beta < 1000)
	{
		double total = 0;
		for (int j = 0; j < num_products; j++)
		{
			prob[j] = pow(prob[j], beta);
			total += prob[j];
		}

		//normalize
		for (int j = 0; j < num_products; j++)
			prob[j] /= total;
	}
	else if (beta >= 1000)
	{
		double max = prob[0];
		int nmax = 1;
		for (int j = 1; j < num_products; j++) {
			if (max < prob[j]) {
				max = prob[j];
				nmax = 1;
			} else if (max == prob[j]) {
				nmax++;
			}
		}

		for (int j = 0; j < num_products; j++)
			prob[j] = (prob[j] == max) ? (1.0 / nmax) : 0;
	}
}

//This is the entry point function for the threads
static void* simulate(void *extra)
{
	int tid = *((int*)extra); //unique thread id from 0 to num_threads-1

	//Main loop
	int my_iteration=0; //number of iterations done by this thread.  Every 1000 iterations we try to sync with the main thread.
	double prob[num_products]; //probability vector of a recommendation
	double* p; //probability vector to choose a product from
	while(true)
	{
		long u = gsl_rng_uniform_int(rngs[tid], N+1); //pick a random user u
		p=backgrounds[u]; //u will buy according to his background distribution unless he passes the following test

		if(in_degrees[u]>0 && gsl_rng_uniform(rngs[tid]) <= alpha ) //with probability alpha u looks for a recommendation
		{
			long v = in_neighbors[u][ gsl_rng_uniform_int(rngs[tid], in_degrees[u]) ]; //v is the friend of u who will provide the recommendation

			//u will buy according to v's history
			p = prob;
			pthread_mutex_lock(&mutex_vertex[v]);
			for(int j=0; j<num_products; j++)
				prob[j] = ((double)owned_products[v][j])/owned_products_no[v];
			pthread_mutex_unlock(&mutex_vertex[v]);

			apply_beta(prob); //modify the probability distribution according to beta
		}

		int product = sample_product(rngs[tid], p); //pick a product from the distribution

		//update history and owned products
		pthread_mutex_lock(&mutex_vertex[u]);
		owned_products_no[u]++;
		owned_products[u][product]++;

#ifdef BOUNDED_HISTORY
		bought_products_no[u]++;
		bought_products[u][product]++;

		history[u]->push_back( product );

		if(history[u]->size() > HISTORY_BOUND)
		{
			owned_products_no[u]--;
			int top= history[u]->front();
				owned_products[u][top]--;

			history[u]->pop_front();
		}
#endif
		pthread_mutex_unlock(&mutex_vertex[u]);

		my_iteration++;
		if(my_iteration<1000)
			continue;


		//Time to try to sync with me main thread
		pthread_mutex_lock(&mutex_iteration);

		iteration+=my_iteration;
		bool wait_for_writer = (iteration > next_write); //should we stop to let the main thread collect the stats?

		pthread_mutex_unlock(&mutex_iteration);

		if(wait_for_writer)
		{
			pthread_mutex_lock(&mutex_waiting_writer[tid]);

			waiting_writer[tid]=true; //signal that we are waiting
			pthread_cond_signal(&conds_waiting_writer_set[tid]);

			while( waiting_writer[tid] ) //wait for the main thread to finish writing
				pthread_cond_wait(&conds_waiting_writer_unset[tid], &mutex_waiting_writer[tid]);

			pthread_mutex_unlock(&mutex_waiting_writer[tid]);

			if(shutdown) //has the simulation ended?
				return NULL;
		}

		my_iteration=0; //perform another 1000 iterations
	}

	return NULL;
}
