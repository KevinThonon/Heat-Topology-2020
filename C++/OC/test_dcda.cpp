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


int main(int argc, char *argv[]) {

	int N = atoi(argv[1]);
	double pctmetal = atof(argv[2]);
	double penal = atof(argv[3]);

	vec a(N*N);
	for (int i = 0; i < N*N; ++i) {
		a(i) = pctmetal; 
	}
	
	mat k = top::create_k(a, N, penal);
	vec rl = top::RL(N);
	sp_mat ll = top::K_mat(k, N);
	vec u = spsolve(ll,rl,"lapack");

	double cost = top::objective_function1(u, N);
	vec lambda = top::lambda1(u, ll, N);

	//double cost = top::objective_function2(u, N);
	//vec lambda = top::lambda2(ll, N);

	//double cost = top::objective_function3(u, N);
	//vec lambda = top::lambda3(ll, N);

	vec dcda_fd = top::dcda_fd(u, rl, a, N, penal);
	vec dcda_arit = top::dcda_arit(lambda, u, a, N, penal);
	vec dcda_harm = top::dcda_harm(lambda, u, a, k, N, penal);

	vec difference(N*N);

	for (int i = 0; i < N*N ; i++){
		difference(i) = (dcda_harm(i)-dcda_arit(i))/cost;
	}
	
	std::cout<<difference<<std::endl;
	//std::cout<<" "<<std::endl;
	//std::cout<<dcda_fd<<std::endl;
	//std::cout<<" "<<std::endl;
	//std::cout<<dcda_harm<<std::endl;

}