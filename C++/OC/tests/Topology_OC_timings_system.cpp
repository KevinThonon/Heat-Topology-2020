// Project Mathematical Engineering: Topology Optimization
// 2019-2020
// Names: Kevin Thonon, Phillipe Dewart, Urbaan Lemmens


#include "Topology_project_OC.hpp"
#include <armadillo>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <math.h>
#include <chrono>

using namespace std;
using namespace arma;
using namespace top;
using namespace chrono; 


int main(int argc, char *argv[]) {

	for (int i = 0; i < 101; i+=10){

	// Start of the clock for timing of the program.

	auto start = high_resolution_clock::now();

	// Command line arguments: gridsize, initial guess and penalty 

	int N = atoi(argv[1]);
	N = N + i;
	cout<<"N = "<<N<<endl;
 	double pctmetal = atof(argv[2]);
	double penal = atof(argv[3]);

	// Pre-allocation of vectors and matrices and declaring some variables.

	vec a(N*N);
	mat k(N,N);
	vec f(pow(N+1,2));
	sp_mat K(((N+1)*(N+1)), ((N+1)*(N+1)));
	vec u((N+1)*(N+1));


	int iterations = 0;





	// fill "a" with the initial guess from the command line arguments

	a.fill(pctmetal);

	// Compute the RL of K*u = f

	f = RL(N);

	// Start of the iterations.

	while (iterations < 5) {

	cout<<"iteration = "<<iterations<<endl;


	
	// Creating and solving the system K*u=f

	k = create_k(a, N, penal);
	K = K_mat(k, N);
	u = spsolve(K,f,"lapack");


	// Incrementing the iteration number.

	iterations += 1;
	}

	// A clock to compute the timing of this program in seconds.

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start)/5;
	cout << duration.count() << endl;
	} 
}