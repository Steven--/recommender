#include <iostream>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/IterativeLinearSolvers>
#include <random>
#include <cmath>
#include <unistd.h>
#include <string.h>
#include <iomanip>

#include <gsl/gsl_rng.h>
#include "rtnorm.hpp"

//If DIRECTED is not defined, then each edge (u,v) in the input graph will be tread as the two edges (u,v) and (v,u)
//#define DIRECTED

using namespace std;
using namespace Eigen;

typedef Eigen::Triplet<double> T;

//Writes the vector to v a file.
//the file name is composed of the base name of the path stored in "basename" plus the given "suffix"
void write_vector(VectorXd v, const char* basename, const char* suffix)
{
	char filename[strlen(basename)+strlen(suffix)+1];

	char *name = NULL;
	char *c = (char*)basename;
	while(c!=NULL && *c!='\0')
	{
		name=c;
		c=strstr(c, "/");
		if(c!=NULL)
		c++;
	}

	strcpy(filename, name);
	strcat(filename, suffix);

	ofstream file;
	file << fixed << setprecision(10);
	file.open(filename);
	file << v << endl;
	file.close();
}

int main(int argc, char** argv)
{

	//Parse parameters
	if(argc<4 || argc>5)
	{
		cout << "Usage: " << argv[0] << " graph alpha seed [background]" << endl;
		return EXIT_FAILURE;
	}
	bool use_bf=(argc==5); //should we use the background (personal distribution) vector given as a parameter?
	
	double alpha = atof(argv[2]);
	int seed = atoi(argv[3]);


	//Setup and seed the random number generator
	gsl_rng_env_setup();                          // Read variable environnement
	const gsl_rng_type* type = gsl_rng_default;   // Default algorithm 'twister'
	gsl_rng *gsl_gen = gsl_rng_alloc (type);      // Rand generator allocation
	gsl_rng_set(gsl_gen, 42+seed);



	//Create a list of elements for the adjacency matrix: T(i,j,x) means that the element of coordinates (i,j) has value x
	vector<T> tripletList;
	FILE* input = fopen(argv[1], "r");

	char buf[2048];

	int u,v;
	int N=0;
	while(fgets(buf, 2048, input)!=NULL)
	{
		if(buf[0]=='#')
				continue;

		sscanf(buf, " %d %d ", &u, &v);

		N = (N>u)?N:u;
		N = (N>v)?N:v;

#ifndef DIRECTED
		tripletList.push_back(T(u,v,1));
#endif
		tripletList.push_back(T(v,u,1)); //v influences u, we transpose the graph
	}
	N++;
	fclose(input);
	input=NULL;


//Should we add an hub vertex that influences everyone?
#ifdef ADD_HUB
	for(int i=0; i<N; i++) //the hub is vertex N
	{
		tripletList.push_back(T(N,i,1)); //N influences i
	}
	N++;
#endif

	cout << "Created list" << endl;

	//Create the actual Adjacency matrix
	Eigen::SparseMatrix<double> Adj(N,N);
	Adj.setFromTriplets(tripletList.begin(), tripletList.end());
	tripletList.clear();

	cout << "Created adjacency matrix" << endl;
	
	VectorXd ones = Eigen::VectorXd::Constant(N, 1);

	//Compute indegree and outdegrees
	VectorXd outdeg = Adj*ones;
	VectorXd indeg = (ones.transpose()*Adj).transpose();

	//Write the computed out/in-degrees to a file
	if( !use_bf )
	{
		write_vector(indeg, argv[1], ".ideg");
		write_vector(outdeg, argv[1], ".odeg");
	}


	cout << "Computed degrees" << endl;


	//Create a diagonal matrix with the inverse outdegree of vertices. This is used to normalize the adjacency matrix.
	DiagonalMatrix<double, Dynamic> inv_deg(N);
	for(int i = 0; i<N; i++)
	{
			if(indeg(i)==0)
					inv_deg.diagonal()(i)=0;
			else
					inv_deg.diagonal()(i) = 1.0/indeg(i);
	}


	//Create the normalized transposed adjacency matrix
	SparseMatrix<double> M = inv_deg *Adj.transpose(); // (A * inv_deg)^T
	Adj.resize(0,0); //destroy A

	cout << "Created M" << endl;

	//Initialize the vector A. A(i) must be 0 if user i has no in-neighbors.
	DiagonalMatrix<double, Dynamic> A(N);
	A.diagonal() = VectorXd::Constant(N, alpha);

	for(int i = 0; i<N; i++)
	{
			if(indeg(i)==0)
					A.diagonal()(i)=0;
	}

	cout << "Computing gamma" << endl;

	//Use the iterative relation (see the paper) to approximate the infliuence vector "gamma"
	RowVectorXd gamma = RowVectorXd::Constant(N, 1);
	RowVectorXd x = RowVectorXd::Constant(N, 1);
	for(int i=0; i<50; i++)
	{
			x=x*A*M;
			gamma += x;
			cerr << x.maxCoeff() << " " << x.minCoeff() << flush << endl;
	}

	gamma = gamma - gamma*A;

	//Write the influence vector to a file
	if( !use_bf )
		write_vector(gamma.transpose(), argv[1], ".gamma");

    //We solve the following system in order to compute the equilibrium vector x for different vectors b
    // x = \hat M b
    // x = (I-AM)^-1 (I-A) b
    // (I-AM)x = b - Ab

    SparseMatrix<double> I(N,N);
    I.setIdentity();

    //SparseLU< SparseMatrix<double> > solver;
    BiCGSTAB< SparseMatrix<double> > solver;
    solver.compute(I-A*M);
    if(solver.info()!=Success)
    {
      cerr << "Decomposition failed" << endl;
      return EXIT_FAILURE;
    }

    cerr << "Decomposition done" << endl;
 
  
    //if we have a input background file, we just solve the system for that background
	if( use_bf )
	{
		VectorXd bf(N);

		printf("Using background file\n");
		FILE* input = fopen(argv[4], "r");

		double d;
		int i=0;
		while(fgets(buf, 2048, input)!=NULL)
		{
			if(buf[0]=='#')
				continue;

			sscanf(buf, " %lf ", &d);
			bf(i)=d;
			i++;
		}
		if(i!=N)
		{
			printf("Wrong background file length\n");
			return EXIT_FAILURE;
		}

        VectorXd yf = solver.solve(bf-A*bf); //Solve the system
        if(solver.info()!=Success)
        {
           cerr << "Solve failed" << endl;
           return EXIT_FAILURE;
        }
        
        write_vector(yf, argv[1], ".yf"); //write the solution to a file
	}
  

	//if we are not given a background file as an input, then we generate seval background vectors (from different distributions)
	if(!use_bf)
	{
		VectorXd bn(N); //normal
		VectorXd bu(N); //uniform
		VectorXd bp(N); //power-law
		VectorXd be(N); //exponential
		for(int i=0; i<N; i++)
		{
			//Power-law in [0.0.1,1]
			double r = gsl_rng_uniform(gsl_gen);
			double L=0.01;
			double H=1;
			double a=0.01;
			bp(i) = pow( (r*(pow(L,a) - pow(H,a)) + pow(H,a))/(pow(L*H, a))  , -1.0/a );

			//Uniform in [0,1]
			bu(i) = gsl_rng_uniform_pos(gsl_gen);

			//Truncated normal in [0,1]
			bn(i)=rtnorm(gsl_gen, 0, 1, 0.5, 0.5/3).first;

			//Truncated exponential in [0,1]
			double lambda=0.5;
			r=gsl_rng_uniform_pos(gsl_gen)*(1-exp(-1.0/lambda)); //1.0 is the truncation
			be(i) = -log(1-r)*lambda;
		}


//The hub always has a background of 1 (as he only recommends a single product)
#ifdef ADD_HUB
		bp(N-1) = 1;
		be(N-1) = 1;
		bu(N-1) = 1;
		bn(N-1) = 1;
#endif
 

	//Solve!

    VectorXd yn = solver.solve(bn-A*bn);
    if(solver.info()!=Success) {
      cerr << "Solve failed" << endl;
      return EXIT_FAILURE;
    }

    VectorXd yu = solver.solve(bu-A*bu);
    if(solver.info()!=Success) {
      cerr << "Solve failed" << endl;
      return EXIT_FAILURE;
    }

    VectorXd yp = solver.solve(bp-A*bp);
    if(solver.info()!=Success) {
      cerr << "Solve failed" << endl;
      return EXIT_FAILURE;
    }

    VectorXd ye = solver.solve(be-A*be);
    if(solver.info()!=Success) {
      cerr << "Solve failed" << endl;
      return EXIT_FAILURE;
    }


    //Write the backgrounds and the solutions to files
    write_vector(bu, argv[1], ".bu");
    write_vector(bn, argv[1], ".bn");
    write_vector(bp, argv[1], ".bp");
    write_vector(be, argv[1], ".be");
    write_vector(yu, argv[1], ".yu");
    write_vector(yn, argv[1], ".yn");
    write_vector(yp, argv[1], ".yp");
    write_vector(ye, argv[1], ".ye");
}

    return EXIT_SUCCESS;
}
