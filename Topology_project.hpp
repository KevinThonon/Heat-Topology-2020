#ifndef top_Topology_project_hpp
#define top_Topology_project_hpp
#include "vector.hpp"
#include "vector_expressions.hpp"
#include "matrix.hpp"
#include <cassert>
#include <iostream>
#include <cmath>
#include <typeinfo>
#include <algorithm>
#include <random>
#include <chrono>
#include <iterator>
#include <armadillo>
#include <nlopt.h>
#include <math.h>

using namespace std;
using namespace arma;

namespace top {

// variabelen definiëren

double penal = 1.0;

mat K_0() {

mat K_0(4,4);

K_0(0,0) = 2.0/3.0;
K_0(0,1) = -1.0/6.0;
K_0(0,2) = -1.0/3.0;
K_0(0,3) = -1.0/6.0;
K_0(1,0) = -1.0/6.0;
K_0(1,1) = 2.0/3.0;
K_0(1,2) = -1.0/6.0;
K_0(1,3) = -1.0/3.0;
K_0(2,0) = -1.0/3.0;
K_0(2,1) = -1.0/6.0;
K_0(2,2) = 2.0/3.0;
K_0(2,3) = -1.0/6.0;
K_0(3,0) = -1.0/6.0;
K_0(3,1) = -1.0/3.0;
K_0(3,2) = -1.0/6.0;
K_0(3,3) = 2.0/3.0;

return K_0;

}






mat create_k(mat pctmetal, int N) {
	mat k(N,N);
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			k(i,j) = (65-0.2)*pow(pctmetal(i,j),penal) + 0.2;
		}
	}
	return k;
}

mat bigkmat(mat k, int N){
	mat bigmat(2*N+1, 2*N+1);
	bigmat.fill(0.0);
	for (int i = 1; i < 2*N+1; i+=2){
		for (int j = 1; j < 2*N+1; j+=2){
			bigmat(i,j) = k((i-1)/2, (j-1)/2);
		}
	}
	
	for (int i = 1; i < 2*N; i+=2){
		for (int j = 2; j < 2*N; j+=2){
			bigmat(i,j) = (k((i-1)/2, j/2) + k((i-1)/2, j/2 - 1))/2;
		}
	}
	
	for (int i = 2; i < 2*N; i+=2){
		for (int j = 1; j < 2*N; j+=2){
			bigmat(i,j) = (k(i/2, (j-1)/2) + k(i/2 - 1, (j-1)/2))/2;
		}
	}
	
	
	for (int i = 1; i < 2*N; i+=2){
		
		bigmat(i, 0) = k((i-1)/2, 0);
		bigmat(0, i) = k(0, (i-1)/2);
		bigmat(i, 2*N) = k((i-1)/2, N-1);
		bigmat(2*N, i) = k(N-1, (i-1)/2);
		
	}
	return bigmat;
}

vec RL(int N){
	double q = 2.0/(0.01*0.01*0.001);
	double h = 0.01/N;
	double qhp = - q * pow(h,2.0)/4.0;
	double qn = - q * pow(h,2.0)/2.0;
	double qq = - q * pow(h,2.0);
	vec rl(pow(N+1,2));
	rl.fill(qq);

	
	rl(0) = qhp;
	rl(N) = qhp;
	rl(pow(N,2)+N) = qhp; 
	rl(pow(N,2)+2*N) = qhp;
	
	for (int i = 0.3*N; i < 0.7*N + 1; i++){
		rl(i) = 293.0;
		rl(pow(N+1,2)-i-1) = 293.0;
	}
	
	for (int i = N+1; i < N+1; i+=(N+1)){
		rl(i) = qn;
		rl(i + N) = qn;
	}
	
	for (int i = 1; i < 0.3*N; i++){
		rl(i) = qn;
		rl(N-i) = qn;
		rl(pow(N, 2) + N + i) = qn;
		rl(pow(N+1,2) - 1 - i) = qn;
	}
	
	return rl;
}

mat LL(mat bigkmat, int N){
	mat ll(((N+1)*(N+1)), ((N+1)*(N+1)));
	ll.fill(0.0);
	for (int i = 0.3*N; i < 0.7*N + 1; i++){
		ll(i,i) = 1.0;
		ll(((N+1)*(N+1))-i-1, ((N+1)*(N+1))-i-1) = 1.0;
	}
	
	//Linksboven hoekpunt
	ll(0,0) = -0.5*(bigkmat(1,0)+bigkmat(0,1));
	ll(0,1) = 0.5*bigkmat(1,0);
	ll(0,N+1) = 0.5*bigkmat(0,1);
	
	//Linksonder hoekpunt
	ll(N,0) = -0.5*(bigkmat(2*N-1,0)+bigkmat(2*N-1,1));
	//std::cout<<N<<std::endl;
	ll(N,N-1) = 0.5*bigkmat(2*N-1,0);
	ll(N,2*N+1) = 0.5*bigkmat(2*N,1);
	
	//Rechtsboven hoekpunt
	ll(N*(N+1),N*(N+1)) = -0.5*(bigkmat(0,2*N-1)+bigkmat(1,2*N));
	ll(N*(N+1),N*(N+1)+1) = 0.5*bigkmat(1,2*N);
	ll(N*(N+1),N*N-1) = 0.5*bigkmat(0,2*N-1);
	
	//Rechtsonder hoekpunt
	ll(N*(N+2),N*(N+2)) = -0.5*(bigkmat(2*N-1,2*N)+bigkmat(2*N,2*N-1));
	ll(N*(N+2),N*(N+1)-1) = 0.5*bigkmat(2*N,2*N-1);
	ll(N*(N+2),N*(N+2)-1) = 0.5*bigkmat(2*N-1,2*N);
	
	//Neumann links
	for (int i = 1; i < 0.3*N; i++){
		ll(i,i) = -0.5*bigkmat(2*i-1,0) -0.5*bigkmat(2*i+1,0) - bigkmat(2*i,1);
		ll(i,i-1) = 0.5*bigkmat(2*i-1,0);
		ll(i,i+1) = 0.5*bigkmat(2*i+1,0);
		ll(i,i+N+1) = bigkmat(2*i,1);
		
		ll(N-i,N-i) = -0.5*bigkmat(2*(N-i)-1,0) -0.5*bigkmat(2*(N-i)+1,0) - bigkmat(2*(N-i),1);
		ll(N-i,N-i-1) = 0.5*bigkmat(2*(N-i)-1,0);
		ll(N-i,N-i+1) = 0.5*bigkmat(2*(N-i)+1,0);
		ll(N-i,N-i+N+1) = bigkmat(2*(N-i),1);
	}

	//Neumann rechts
	for (int i = 1; i < 0.3*N; i++){
		ll(N*(N+1)+i,N*(N+1)+i) = -0.5*bigkmat(2*i-1,2*N) -0.5*bigkmat(2*i+1,2*N) - bigkmat(2*i,2*N-1);
		ll(N*(N+1)+i,N*(N+1)+i-1) = 0.5*bigkmat(2*i-1,2*N);
		ll(N*(N+1)+i,N*(N+1)+i+1) = 0.5*bigkmat(2*i+1,2*N);
		ll(N*(N+1)+i,N*(N+1)+i-(N+1)) = bigkmat(2*i,2*N-1);
		
		ll(N*(N+2)-i,N*(N+2)-i) = -0.5*bigkmat(2*(N-i)-1,2*N) -0.5*bigkmat(2*(N-i)+1,2*N) - bigkmat(2*(N-i),2*N-1);
		ll(N*(N+2)-i,N*(N+2)-i-1) = 0.5*bigkmat(2*(N-i)-1,2*N);
		ll(N*(N+2)-i,N*(N+2)-i+1) = 0.5*bigkmat(2*(N-i)+1,2*N);
		ll(N*(N+2)-i,N*(N+2)-i-(N+1)) = bigkmat(2*(N-i),2*N-1);
	}
	
	//Neumann boven
	for (int i = 1; i < N; i++){
		ll(i*(N+1),i*(N+1)) = -0.5*bigkmat(0,2*i-1) -0.5*bigkmat(0,2*i+1) - bigkmat(1,2*i);
		ll(i*(N+1),i*(N+1)+1) = bigkmat(1,2*i);
		ll(i*(N+1),(i-1)*(N+1)) = 0.5*bigkmat(0,2*i-1);
		ll(i*(N+1),(i+1)*(N+1)) = 0.5*bigkmat(0,2*i+1);
	}
	
	//Neumann onder
	for (int i = 1; i < N; i++){
		ll(i*(N+1)+N,i*(N+1)+N) = -0.5*bigkmat(2*N,2*i-1) -0.5*bigkmat(2*N,2*i+1) - bigkmat(2*N-1,2*i);
		ll(i*(N+1)+N,i*(N+1)+N-1) = bigkmat(2*N-1,2*i);
		ll(i*(N+1)+N,(i-1)*(N+1)+N) = 0.5*bigkmat(2*N,2*i-1);
		ll(i*(N+1)+N,(i+1)*(N+1)+N) = 0.5*bigkmat(2*N,2*i+1);
	}
	
	//Gewone kolommen
	for (int j = 1; j < N; j++){
		for (int i = 1; i < N; i++){
			ll(j*(N+1)+i,j*(N+1)+i) = -(bigkmat(2*i-1,2*j) + bigkmat(2*i+1,2*j) + bigkmat(2*i,2*j-1) + bigkmat(2*i,2*j+1));
			ll(j*(N+1)+i,j*(N+1)+i-1) = bigkmat(2*i-1,2*j);
			ll(j*(N+1)+i,j*(N+1)+i+1) = bigkmat(2*i+1,2*j);
			ll(j*(N+1)+i,(j-1)*(N+1)+i) = bigkmat(2*i,2*j-1);
			ll(j*(N+1)+i,(j+1)*(N+1)+i) = bigkmat(2*i,2*j+1);
		}
	}
	
	return ll;
}

double objective_function(vec T, mat K){
	double cost = dot(0.5*T.t(),K*T);
	return cost;
}


mat check(int N, double rmin, mat x, mat dc){

mat dcn = zeros<mat>(N,N);

for (int i = 1; i < N-2; i++) {
  for (int j = 1; j < N-2; j++) {
    double sum = 0.0; 
    for (int k = std::max(i - floor(rmin),2.0); k < std::min(i + floor(rmin),N-2.0); k++) {
      for (int l = std::max(j - floor(rmin),2.0); l< std::min(j + floor(rmin),N-2.0); l++) {
        double fac = rmin - sqrt(pow((i-k),2) + pow((j-l),2));
        sum = sum + std::max(0.0,fac);
        dcn(j,i) = dcn(j,i) + std::max(0.0,fac) * dc(l,k) * x(l,k);
      }
    }
    dcn(j,i) = dcn(j,i) / (x(j,i) * sum);
  }
}
return dcn;
}






}
#endif


