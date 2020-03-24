#include "vector.hpp"
#include "matrix.hpp"
#include "Topology_project.hpp"
#include <armadillo>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;
using namespace arma;


int main(int argc, char* argv[]) {

int N = 10;

mat pctmetal(N,N);
pctmetal.fill(0.4);

mat k = top::create_k(pctmetal, N);
//mat bigmat = top::bigkmat(k, N);
vec rl = top::RL(N);
//mat ll = top::LL(bigmat, N);
mat ll = top::K_mat(k, N);

vec sol = solve(ll,rl);
double cost = top::objective_function(sol, N);

vec lambda = top::lambda(sol, ll, N);
mat dcda = top::dcda(lambda, sol, pctmetal, N);

ofstream myfile;
myfile.open ("results.txt");
myfile <<sol;
myfile.close();
return cost;



}
