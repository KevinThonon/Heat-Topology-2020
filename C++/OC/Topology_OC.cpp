#include "Topology_project_OC.hpp"
#include <armadillo>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <math.h>
#include <nlopt.h>

using namespace std;
using namespace arma;
using namespace top;


int main(int argc, char *argv[]) {

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
	vec dcda_a(N*N); // arithmetic
	vec dcda_h(N*N); // harmonic
	vec dcda_f(N*N); // finite differences
	vec difference(N*N);
	vec dcda_check(N*N);
	vec a_old(N*N);
	vec change_v(N*N);
	int iterations = 0;
	double change = 1.0; 

	while (change > 0.01) {

	cout<<"iteration = "<<iterations<<endl;

	if (iterations == 5){
		penal = 2.0;
	}

	if (iterations == 10){
		penal = 3.0;
	}
 
	if (iterations == 15){
		penal = 4.0;
	} 

	if (iterations == 20){
		penal = 5.0;
	}

	if (iterations == 25){
		penal = 10.0;
	}
 
	if (iterations == 25){
		penal = 9.0;
	} 
	
	k = create_k(a, N, penal);
	K = K_mat(k, N);
	f = RL(N);
	u = spsolve(K,f,"lapack");

  
	temperatureToTxtFile(u, iterations, N);

	//double cost = objective_function1(u, N);
	lambda = lambda1(u, K, N);

	//double cost = objective_function2(u, N);
	//vec lambda = lambda2(K, N);

	//double cost = objective_function3(u, N);
	//vec lambda = lambda3(K, N);

	//vec dcda_f = dcda_fd(u, f, a, N, penal);
	//dcda_a = dcda_arit(lambda, u, a, N, penal);
	vec dcda_h = dcda_harm(lambda, u, a, k, N, penal);

	double rmin = 2.0;
	dcda_check = check(N, rmin, a, dcda_h);

	gradientToTxtFile(dcda_check, iterations, N);

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
}