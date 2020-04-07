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


int main() {

	int N = 100;
	int iterations = 0;

	vec a(N*N);
	for (int i = 0; i < N*N; ++i) {
		a(i) = 0.3; 
	}

	double change = 1.0; 
	double penal = 2.0;

	while (change > 0.01) {

	std::cout<<"iteration = "<<iterations<<std::endl;

	
	if (iterations == 5){
		penal = 3.0;
	}

	if (iterations == 15){
		penal = 4.0;
	}

	if (iterations == 20){
		penal = 6.0;
	}
 
	if (iterations == 25){
		penal = 7.0;
	}
	
	mat k = top::create_k(a, N, penal);
	vec rl = top::RL(N);
	sp_mat ll = top::K_mat(k, N);
	vec u = spsolve(ll,rl,"lapack");

  
  	string t = "temperature_";
  	t += to_string(iterations);
  	t += ".txt";
	string path_temperature = "/Users/Urban/Documents/GitHub/Heat-Topology-2020/C++/OC/solutions_OC/temperature/";
	path_temperature += t;

	ofstream temperature_file;
        temperature_file.open(path_temperature);
        for (int i = 0; i < (N+1)*(N+1); ++i) {
    		temperature_file <<u(i)<<std::endl;
    	}
    	temperature_file.close();

	//double cost = top::objective_function1(u, N);
	vec lambda = top::lambda1(u, ll, N);

	//double cost = top::objective_function2(u, N);
	//vec lambda = top::lambda2(ll, N);

	//double cost = top::objective_function3(u, N);
	//vec lambda = top::lambda3(ll, N);

	// vec dcda = top::dcda_fd(u, rl, a, N, penal);
	//vec dcda = top::dcda_arit(lambda, u, a, N, penal);
	vec dcda = top::dcda_harm(lambda, u, a, k, N, penal);

	double rmin = 2.0;
	vec dcda_check = top::check(N, rmin, a, dcda);
	vec a_old = a;

	a = top::OC(N, a, 0.4, dcda_check);
	
	vec change_v (N*N);
	double change_n = 0.0;

	for (int i = 0; i < N*N ; i++){
		change_v(i) = a(i) - a_old(i);  
		change_n = max(change_v(i),change_n);
	} 

	change = change_n;
	std::cout<<"change ="<<change<<std::endl;

  	string g = "gradient_";
  	g += to_string(iterations);
  	g += ".txt";
	string path_gradient = "/Users/Urban/Documents/GitHub/Heat-Topology-2020/C++/OC/solutions_OC/gradient/";
	path_gradient += g;

	ofstream gradient_file;
        gradient_file.open(path_gradient);
        for (int i = 0; i < N*N; ++i) {
    		gradient_file <<dcda(i)<<std::endl;
    	}
    	gradient_file.close();

  	string m = "metal_";
  	m += to_string(iterations);
  	m += ".txt";
	string path_metal = "/Users/Urban/Documents/GitHub/Heat-Topology-2020/C++/OC/solutions_OC/metal/";
	path_metal += m;

	ofstream metal_file;
        metal_file.open(path_metal);
        for (int i = 0; i < N*N; ++i) {
    		metal_file <<a[i]<<std::endl;
    	}
    	metal_file.close();
	iterations += 1;
	if (change <= 0.1){
		penal += 1.0;
	}
	}
	


}