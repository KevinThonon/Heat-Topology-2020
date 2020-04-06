#include "vector.hpp"
#include "matrix.hpp"
#include "Topology_project.hpp"
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

	int N = 50;
	int iterations = 0;

	vec a(N*N);
	for (int i = 0; i < N*N; ++i) {
		a(i) = 0.3; 
	}

	double change = 1.0; 

	while (change > 0.02) {
	std::cout<<"iteration = "<<iterations<<std::endl;

	mat k = top::create_k(a, N);
	vec rl = top::RL(N);
	sp_mat ll = top::K_mat(k, N);
	vec u = spsolve(ll,rl,"lapack");

  
  	string t = "temperature_";
  	t += to_string(iterations);
  	t += ".txt";
	string path_temperature = "/Users/Urban/Documents/GitHub/Heat-Topology-2020/temperature/";
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

	//vec dcda = top::dcda(lambda, u, a, N);
	vec dcda = top::dcda_harm(lambda, u, a, k, N);

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

  	string g = "gradient_";
  	g += to_string(iterations);
  	g += ".txt";
	string path_gradient = "/Users/Urban/Documents/GitHub/Heat-Topology-2020/gradient/";
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
	string path_metal = "/Users/Urban/Documents/GitHub/Heat-Topology-2020/metal/";
	path_metal += m;

	ofstream metal_file;
        metal_file.open(path_metal);
        for (int i = 0; i < N*N; ++i) {
    		metal_file <<a[i]<<std::endl;
    	}
    	metal_file.close();
	iterations += 1;
	}
	


}