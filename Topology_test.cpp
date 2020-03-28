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


typedef struct {
  mat pctmetal;
} *my_func_data;

double myfunc(unsigned n, const double *a, double *grad, void *data){

	int N = 10;

	// transfer *data (a pointer) to the real matrix pctmetal 
	// (algorithm is way faster because not every time a new matrix is being build every
	// iteration)

	my_func_data &temp = (my_func_data &) data;
	mat pctmetal = temp->pctmetal;


	// put the percentages of metal in the pctmetal matrix. 
	// (initial guess a is 0.4 for every element)

	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			pctmetal(i,j) = a[i*N + j];
		}
	}

	mat k = top::create_k(pctmetal, N);
	vec rl = top::RL(N);
	mat ll = top::K_mat(k, N);

	vec u = solve(ll,rl);
	double cost = top::objective_function(u, N);

	vec lambda = top::lambda(u, ll, N);
	vec dcda = top::dcda(lambda, u, pctmetal, N);

	if (grad) {
	for (int i = 0; i < dcda.size(); ++i) {
        	grad[i] = dcda(i);
	}
    }
    return cost;
}

// Inequality constraint: sum(a(i,j))/N^2 - 0.4 <= 0

double myconstraint(unsigned n, const double *a, double *grad, void *data){
	int N = 10;
	double sum = 0.0;
	for (int i = 0; i < N*N; ++i) {
		sum += a[i];
	}
	double average_pctmetal = sum/pow(N,2);
	std::cout<<average_pctmetal<<std::endl;
 
    	return average_pctmetal - 0.4;
}




int main() {

int N = 10;
mat pctmetal(N,N);

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

// create algorithm with dimension N^2

opt = nlopt_create(NLOPT_LD_MMA, N*N); /* algorithm and dimensionality */

//set lower and upper bounds (0 <= a(i,j) <= 1)

nlopt_set_lower_bounds(opt, lb);
nlopt_set_upper_bounds(opt, ub);

// create minimization problem with myfunc (function is written in the top of this file)
// a reference to pctmetal is given as data so that 
// not every time in myfunc a new matrix pctmetal is made

nlopt_set_min_objective(opt, myfunc, &pctmetal);

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
    for (int i = 0; i < N*N; ++i) {
    	std::cout<<a[i]<<std::endl;
    }
}

nlopt_destroy(opt);

}