// Project Mathematical Engineering: Topology Optimization
// 2019-2020
// Names: Kevin Thonon, Phillipe Dewart, Urbaan Lemmens


#include "Topology_project_MMA.hpp"
#include <armadillo>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <math.h>
#include <nlopt.h>

using namespace std;
using namespace arma;

double myfunc(unsigned n, const double *a, double *grad, void *data){

	std::cout<<"iteration = "<<*(int *)data<<std::endl;
	*(int *)data += 1;
	
	double penal = 2.0;
	int d = *(int *)data;

	/*
	if (d > 10){
		penal = 4.0;
	}


	if (d > 10){
		penal = 4.0;
	}
	*/
	

	int N = sqrt(n);

	mat k = top::create_k(a, N, penal);
	vec rl = top::RL(N);
	sp_mat ll = top::K_mat(k, N);
	vec u = spsolve(ll,rl,"lapack");

  
  	string t = "temperature_";
  	t += to_string(*(int *)data);
  	t += ".txt";
	string path_temperature = "/Users/Urban/Documents/GitHub/Heat-Topology-2020/C++/MMA/solutions_MMA/temperature/";
	path_temperature += t;

	ofstream temperature_file;
        temperature_file.open(path_temperature);
        for (int i = 0; i < (N+1)*(N+1); ++i) {
    		temperature_file <<u(i)<<std::endl;
    	}
    	temperature_file.close();

	double cost = top::objective_function1(u, N);
	std::cout<<"cost ="<<cost<<std::endl;
	vec lambda = top::lambda1(u, ll, N);

	//double cost = top::objective_function2(u, N);
	//vec lambda = top::lambda2(ll, N);

	//double cost = top::objective_function3(u, N);
	//vec lambda = top::lambda3(ll, N);

	//vec dcda_a = top::dcda_arit(lambda, u, a, N, penal);
	//vec dcda_h = top::dcda_harm(lambda, u, a, k, N, penal);
	vec dcda_f = top::dcda_fd(u, rl, a, N, penal);

	double rmin = 2.0;
	vec dcda_check = top::check(N, rmin, a, dcda_f);

  	string g = "gradient_";
  	g += to_string(*(int *)data);
  	g += ".txt";
	string path_gradient = "/Users/Urban/Documents/GitHub/Heat-Topology-2020/C++/MMA/solutions_MMA/gradient/";
	path_gradient += g;

	ofstream gradient_file;
        gradient_file.open(path_gradient);
        for (int i = 0; i < n; ++i) {
    		gradient_file <<dcda_check(i)<<std::endl;
    	}
    	gradient_file.close();

  	string m = "metal_";
  	m += to_string(*(int *)data);
  	m += ".txt";
	string path_metal = "/Users/Urban/Documents/GitHub/Heat-Topology-2020/C++/MMA/solutions_MMA/metal/";
	path_metal += m;

	ofstream metal_file;
        metal_file.open(path_metal);
        for (int i = 0; i < n; ++i) {
    		metal_file <<a[i]<<std::endl;
    	}
    	metal_file.close();


	if (grad) {
	for (int i = 0; i < n; ++i) {
        	grad[i] = dcda_check(i);
	}
    }

    return cost;
}

// Inequality constraint: sum(a(i,j))/N^2 - 0.4 <= 0

double myconstraint(unsigned n, const double *a, double *grad, void *data){
	double sum = 0.0;
	for (int i = 0; i < n; ++i) {
		sum += a[i];
	}
	double average_pctmetal = sum/n;
	std::cout<<"average pct of metal = "<<average_pctmetal<<" (absolute error = "<<average_pctmetal-0.4<<")"<<std::endl;
 
    	return average_pctmetal - 0.4;
}




int main() {

int N = 40;
int iterations = 0;

// MMA - method

// Contstraints for the percentage of metal: 0.0 <= a(i,j) <= 1.0

// lowerbounds for the parameters a(i,j) (percentage of metal): 0.0 <= a(i,j) 
double lb[N*N];
for (int i = 0; i < N*N; ++i) {
	lb[i] = 0.0; /* lower bounds */
}

// upperbounds for the parameters a(i,j) (percentage of metal): a(i,j) <= 1.0
double ub[N*N];
for (int i = 0; i < N*N; ++i) {
	ub[i] = 1.0; /* upper bounds */
}

// create object

nlopt_opt opt;

// create algorithm with dimension N^2 (unsigned n = N*N)

opt = nlopt_create(NLOPT_LD_MMA, N*N); /* algorithm and dimensionality */

//set lower and upper bounds (0 <= a(i,j) <= 1)

nlopt_set_lower_bounds(opt, lb);
nlopt_set_upper_bounds(opt, ub);

// create minimization problem with myfunc (function is written in the top of this file)
// a reference to pctmetal is given as data so that 
// not every time in myfunc a new matrix pctmetal is made

nlopt_set_min_objective(opt, myfunc, &iterations);

// Add the inequality constraint to the problem with the parameters: opt (object), 
// myconstraint (the sum(a(i,j)/N^2) <= 0.4 constraint), 
// NULL (no data is given) and tolerance (1*10^-8)

nlopt_add_inequality_constraint(opt, myconstraint, NULL, 1e-14);

// stopping criteria
//nlopt_set_xtol_rel(opt, 1e-8);
nlopt_set_maxeval(opt,50);

// Initial guess

double a[N*N];  // some initial guess: average percentage of metal in one element is 0.4
for (int i = 0; i < N*N; ++i) {
	a[i] = 0.01; 
}

// value of the objective function during the iterations.
double minf;

if (nlopt_optimize(opt, a, &minf) < 0) {
    printf("nlopt failed!\n");
}
else {
    std::cout<<"optimization stopped because of the stopping criteria"<<std::endl;
    return 0;
}

nlopt_destroy(opt);

}