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

	// Start of the clock for timing of the program.

	auto start = high_resolution_clock::now();

	// Command line arguments: gridsize, initial guess and penalty 

	int N = atoi(argv[1]);
	double pctmetal = atof(argv[2]);
	double penal = atof(argv[3]);

	// Pre-allocation of vectors and matrices and declaring some variables.

	vec a(N*N);
	mat k(N,N);
	vec f(pow(N+1,2));
	sp_mat K(((N+1)*(N+1)), ((N+1)*(N+1)));
	vec u((N+1)*(N+1));
	vec lambda((N+1)*(N+1));
	vec dcda_check(N*N);
	vec a_old(N*N);
	vec change_v(N*N);

	int iterations = 0;
	double change = 1.0;
	double rmin = 2.0; 

	// Choose which method you want to use for the computation of the gradients.

	//vec dcda_a(N*N); // arithmetic analytic adjoint method
	vec dcda_h(N*N); // harmonic analytic adjoint method
	//vec dcda_f(N*N); // finite differences method


	// fill "a" with the initial guess from the command line arguments

	a.fill(pctmetal);

	// Compute the RL of K*u = f

	f = RL(N);

	// Start of the iterations.

	while (iterations < 33) {

	cout<<"iteration = "<<iterations<<endl;

	
	// after some iterations we update the penalty factor.

	if (iterations == 25){
		penal = 4.0;
	}
	
	if (iterations == 30){
		penal = 5.0;
	}

	
	// Creating and solving the system K*u=f

	k = create_k(a, N, penal);
	K = K_mat(k, N);
	u = spsolve(K,f,"lapack");

	// Choose one of the four cost functions.

	double cost = objective_function1(u, N);
	cout<<"cost = "<<cost<<endl;
	lambda = lambda1(u, K, N);
	
	/*
	double cost = objective_function1w(u, N);
	cout<<"cost = "<<cost<<endl;
	vec lambda = lambda1w(K, u, N);
	*/

	/*
	double cost = objective_function2(u, N);
	cout<<"cost = "<<cost<<endl;
	vec lambda = lambda2(K, N);
	*/

	/*
	double cost = objective_function2w(u, N);
	cout<<"cost = "<<cost<<endl;
	vec lambda = lambda2w(K, N);
	*/

	// Choose which method is used for calculating the gradients.
	// The function "check" is used so no checkerboard pattern occurs.

	/*
	dcda_f = dcda_fd(u, f, a, N, penal);
	dcda_check = check(N, rmin, a, dcda_f);
	*/
	
	/*
	dcda_a = dcda_arit(lambda, u, a, N, penal);
	dcda_check = check(N, rmin, a, dcda_a);
	*/

	dcda_h = dcda_harm(lambda, u, a, k, N, penal);
	dcda_check = check(N, rmin, a, dcda_h);

	// Save the old vector with the percentages of metal.
	
	a_old = a;

	// Optimization method.
	// gridsize, percentages of metal per cell, max avg. amount of metal, gradients
	// N^2 + 1 constraints: 
	// N^2 constraints: 0 <= a(i,j) <= 1 (N constraints)
	// 1 constraint: sum(a(i,j))/N^2 <= 0.4

	a = OC(N, a, 0.4, dcda_check);

	// Writing solutions to txt files.
	// If you would like to choose the path you need to
	// change the path in the .hpp file in these 4 functions.

	/*
	metalToTxtFile(a, iterations, N);
	temperatureToTxtFile(u, iterations, N);
	gradientToTxtFile(dcda_check, iterations, N);
	costToTxtFile(cost, iterations, N);
	*/
	
	// calculating the change between the 2 iterations of the percentages of metal.

	double change_n = 0.0;

	for (int i = 0; i < N*N ; i++){
		change_v(i) = a(i) - a_old(i);  
		change_n = max(change_v(i),change_n);
	} 

	change = change_n;

	cout<<"change ="<<change<<endl;

	// Incrementing the iteration number.

	iterations += 1;
	}

	// A clock to compute the timing of this program in seconds.

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);
	cout << duration.count() << endl; 
}