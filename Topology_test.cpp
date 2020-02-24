#include "vector.hpp"
#include "matrix.hpp"
#include "Topology_project.hpp"
#include <armadillo>

using namespace std;
using namespace arma;

int main(int argc, char* argv[]) {

int N = 10;

mat pctmetal(N,N);
pctmetal.fill(0.4);

mat k = tws::create_k(pctmetal, N);
mat bigmat = tws::bigkmat(k, N);
vec rl = tws::RL(N);
mat ll = tws::LL(bigmat, N);

vec sol = solve(ll,rl);

std::cout<<sol<<std::endl;

}
