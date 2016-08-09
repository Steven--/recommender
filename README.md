
=== Graphs File Format ===

Graph vertices must have consecutive integers IDs, from 0 to N-1, where N is the number of vertices.
Each line of the graph file contains two integers separated by a space.
A line of the form: "A B" (where A and B are integers), means that there is a directed edge from vertex with ID A to vertex with ID B.
If a graph is undirected both "A B" and "B A" must appear in the input file.
Duplicate edges are not allowed.
Lines starting with an hash character (#) are ignored.

=== Compiling the tools ===

In order to compile source you need a C++ compiler, the GNU make program, and the GNU Scientific Library.
Once the dependencies are met it suffices to run "make" in the directory where the provided sources have been extracted.


=== equilibrium and equilibrium_hub ===

The "equilibrium" tool takes a graph as an input, transposes it, and computes the following vectors:

- The outdegree vector: the element in i-th position is the outdegree of the i-th vertex, where i starts from 0.
- The indegree vector.
- The influence vector (gamma).

Moreover it generates four personal distribution vectors where each element is from the following distributions:
- A normal of mean 0.5 and standard deviation 1/6.
- A uniform in [0,1]
- A power-law in [0.01, 1] with exponent 0.01
- An exponential of mean 1/2.

For all the distribution with a support wider than [0,1], the truncated version is used instead.
For each of the above personal distribution vectors the corresponding equilibrium vector is computed.

The program is invoked as follows:

./equilibrium graph alpha seed

where graph is the input graph (see Graph File Format), 0 <= alpha < 1 is the parameter of the model, and seed is any integer value that is used for seeding the random number generator.
Example:
./equilibrium path/to/my_graph.txt 0.5 42

The output vectors are written to the current working directory. Their name is the filename of the input graph (without its leading directories) plus a suffix that depends on the output vector. These suffixes are:

.odeg for the outdegree vector.
.ideg for the indegree vector.
.gamma for the gamma vector.

.bn for the personal distribution vector sampled from a normal distribution.
.bu for the personal distribution vector sampled from a uniform distribution.
.bp for the personal distribution vector sampled from a power-law distribution.
.be for the personal distribution vector sampled from a exponential distribution.

.yX for the equilibrium vector corresponding to the personal distribution vector bX, where X is one of n,u,p,e.

Output vectors contain one floating point number per line. The number on the i-th line corresponds to the value associated to vertex i (starting from 0).


Finally, this tool can also be used to compute an equilibrium vector for a given graph and background vector. It is invoked as follows:
./equilibrium graph alpha seed background

Example:
./equilibrium path/to/my_graph.txt 0.5 42 my_background

In this case the seed is ignores, and only one file is produced, with suffix .yf.


The tool name "equilibrium_hub" acts exactly as "equilibrium" but an additional vertex is added after the graph has been loaded (and transposed).
This vertex has no incoming edge and has an outgoing edge to each other vertex in the graph.


=== simulate and simulate-bounded-100 ===

The "simulate" tool simulates the recommendation process on a given graph and for a given set of parameters as described in the following.
The simulation continues for 100000*N purchases. Every 100*N purchases (and hence for a total of 100 times) the current state of the simulation is printed on the standard output.
Each output consists of an header line "#Iteration X" where is is the number of purchases completed since the beginning of the process, followed by a line for each user in the network.
Each of these line contains K+1 space-separated integers, where K is the number of products, the first integer is the total number of products bought by the corresponding user, while the i-th integer (for 2<=i<=K+1) is the number of products of type i-1 bought by that user.

The personal preference distribution for each user is chosen by splitting the unit interval into K+advantage-1 pieces. The preference of the first product corresponds to the length of the first "advantage" intervals. The preference of the i-th product for i>=2 is the length of the (advantage+i-1)-th interval. I.e., the first product is "advantage" times more likely to be bought than each other product, on average.

This tool can be invoked as follows:

./simulate graph alpha beta num_threads seed background_file num_products advantage

where:

graph: is the path to the file containing the graph to be used. This graph is transposed so that an edge (u,v) in the original graph means that v is influencing the purchases of u.
alpha: is the parameter of the model (0 <= alpha < 1)
beta: is a parameter of our model (refer to the paper for a description of its meaning). A value of 1000 means infinity.
num_threads: is the number of threads to use for the simulation. This can be set to the number of available processors/cores to speed up the simulation.
seed: is an integer value used to seed the random number generator.
background_file: is a file where the generated personal distributions for the users will be saved. It will contain N lines, one for each user. The i-th line has K space-separated floating-point values, the j-th of these values is the probability that the i-th user will select the j-th product when he is not following a recommendation. 
num_products: is the number of different product types to use in the simulation.
advantage: this parameter influences the generation of the personal distribution, as described above.

Example:

./simulate path/to/my_graph.txt 0.2 1 4 42 personal_distributions.txt 4 1


The "simulate-bounded-100" tool acts similarly to simulate and it is invoked in the same way, except that users histories are limited to the last 100 purchases.
Each line of the output of this tool has 2*(K+1) additional columns.
The first K+1 are analogous to the ones described before except that they only report the number of products in the users' histories. The latter K+1 columns report the *overall* number of bought products since the beginning of the simulation.

Example:

./simulate-bounded-100 path/to/my_graph.txt 0.2 1 4 42 personal_distributions.txt 4 1
