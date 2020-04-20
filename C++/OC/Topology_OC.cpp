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
	auto start = high_resolution_clock::now(); 

	int N = atoi(argv[1]);
	double pctmetal = atof(argv[2]);
	double penal = atof(argv[3]);

	vec a(N*N);
	a.fill(pctmetal);
	mat k(N,N);
	vec f(pow(N+1,2));
	sp_mat K(((N+1)*(N+1)), ((N+1)*(N+1)));
	vec u((N+1)*(N+1));
	vec lambda((N+1)*(N+1));
	//vec dcda_a(N*N); // arithmetic
	vec dcda_h(N*N); // harmonic
	//vec dcda_f(N*N); // finite differences
	vec dcda_check(N*N);
	vec a_old(N*N);
	vec change_v(N*N);
	int iterations = 0;
	double change = 1.0; 

	while (iterations < 33) {

	cout<<"iteration = "<<iterations<<endl;

	
	if (iterations == 25){
		penal = 4.0;
	}
	
	if (iterations == 30){
		penal = 5.0;
	}

	
	k = create_k(a, N, penal);
	K = K_mat(k, N);
	f = RL(N);
	u = spsolve(K,f,"lapack");

  
	//temperatureToTxtFile(u, iterations, N);

	double cost = objective_function1(u, N);
	std::cout<<"cost = "<<cost<<std::endl;
	lambda = lambda1(u, K, N);
	
	//double cost = objective_function1w(u, N);
	//std::cout<<"cost = "<<cost<<std::endl;
	//vec lambda = lambda1w(K, u, N);

	//double cost = objective_function2(u, N);
	//vec lambda = lambda2(K, N);

	//double cost = objective_function2w(u, N);
	//vec lambda = lambda2w(K, N);

	//dcda_f = dcda_fd(u, f, a, N, penal);
	//dcda_a = dcda_arit(lambda, u, a, N, penal);
	dcda_h = dcda_harm(lambda, u, a, k, N, penal);

	double rmin = 2.0;
	//dcda_check = check(N, rmin, a, dcda_f);
	//dcda_check = check(N, rmin, a, dcda_a);
	dcda_check = check(N, rmin, a, dcda_h);
	

	gradientToTxtFile(dcda_check, iterations, N);
	costToTxtFile(cost, iterations, N);

	a_old = a;

	a = OC(N, a, 0.4, dcda_check);

	metalToTxtFile(a, iterations, N);
	
	double change_n = 0.0;

	for (int i = 0; i < N*N ; i++){
		change_v(i) = a(i) - a_old(i);  
		change_n = max(change_v(i),change_n);
	} 

	change = change_n;
	cout<<"change ="<<change<<endl;
	iterations += 1;
	}
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);
	cout << duration.count() << endl; 
}