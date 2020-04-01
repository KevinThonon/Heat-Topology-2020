#include "vector.hpp"
#include "matrix.hpp"
#include "Topology_project.hpp"
#include <armadillo>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <nlopt.h>

using namespace std;
using namespace arma;

double myfunc(unsigned n, const double *a, double *grad, void *data){

	std::cout<<"runs"<<std::endl;

	int N = sqrt(n);

	mat k = top::create_k(a, N);
	vec rl = top::RL(N);
	mat ll = top::K_mat(k, N);

	vec u = solve(ll,rl);

	ofstream myfile;
        myfile.open ("temperature.txt");
        for (int i = 0; i < n; ++i) {
    		myfile <<u(i)<<std::endl;
    	}
    	myfile.close();

	double cost = top::objective_function1(u, N);
	vec lambda = top::lambda1(u, ll, N);

	//double cost = top::objective_function2(u, N);
	//vec lambda = top::lambda2(ll, N);

	//double cost = top::objective_function3(u, N);
	//vec lambda = top::lambda3(ll, N);

	vec dcda = top::dcda(lambda, u, a, N);

	ofstream myfile1;
        myfile1.open ("gradient.txt");
        for (int i = 0; i < n; ++i) {
    		myfile1 <<dcda(i)<<std::endl;
    	}
    	myfile1.close();

	if (grad) {
	for (int i = 0; i < n; ++i) {
        	grad[i] = dcda(i);
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
	std::cout<<"constraint ="<<average_pctmetal - 0.4<<std::endl;
 
    	return average_pctmetal - 0.4;
}




int main() {

int N = 10;

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

nlopt_set_min_objective(opt, myfunc, NULL);

// Add the inequality constraint to the problem with the parameters: opt (object), 
// myconstraint (the sum(a(i,j)/N^2) <= 0.4 constraint), 
// NULL (no data is given) and tolerance (1*10^-8)

nlopt_add_inequality_constraint(opt, myconstraint, NULL, 1e-8);

// stopping criteria
nlopt_set_xtol_rel(opt, 1e-4);

// Initial guess

double a[N*N];  // some initial guess: average percentage of metal in one element is 0.4
for (int i = 0; i < N*N; ++i) {
	a[i] = 0.2; 
}

// value of the objective function during the iterations.
double minf;

if (nlopt_optimize(opt, a, &minf) < 0) {
    printf("nlopt failed!\n");
}
else {
    ofstream myfile;
    myfile.open ("metal.txt");
    for (int i = 0; i < N*N; ++i) {
    	myfile <<a[i]<<std::endl;
    }
    myfile.close();
    return 0;
}

nlopt_destroy(opt);

}