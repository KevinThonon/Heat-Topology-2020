#ifndef tws_Topology_project_hpp
#define tws_Topology_project_hpp
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

namespace tws {



int N = 10;

tws::matrix<double> create_pctmetal() {
	tws::matrix<double> pctmetal(N,N,0.0);
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			pctmetal(i,j) = 0.4;
		}
	}
	return pctmetal;
}

tws::matrix<double> create_k(tws::matrix<double> pctmetal) {
	tws::matrix<double> k(N,N,0.0);
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			k(i,j) = (65-0.2)*pctmetal(i,j) + 0.2;
		}
	}
	return k;
}

int pos(int i, int j) {
	int ij;
	ij = i+N*j;
	return ij;
}

double meank(double k1, double k2) {
	double m;
	m = (k1+k2)/2; //arithmetic
    	//m = 2/(1/k1+1/k2); //harmonic
	return m;
}


tws::vector<double> xi(int i, int j, tws::matrix<double> k) {
	
	tws::vector<double> result(4,0.0);

    	if ( i > 1 ) {
        	result(0) = meank(k(i,j),k(i-1,j));
    	}

    	if ( i < N ) {
       		result(1) = meank(k(i,j),k(i+1,j));
    	}

    	if ( j > 1 ) {
        	result(2) = meank(k(i,j),k(i,j-1));
    	}

    	if ( j < N ) {
        	result(3) = meank(k(i,j),k(i,j+1));
   	} 

	return result;



}

tws::matrix<double> create_K(tws::matrix<double> k) {
	tws::matrix<double> K(N*N,N*N,0.0);
	tws::vector<double> xi_list(4,0.0);
	for (int i = 1; i < N - 1; i++) {
    		for (int j = 1; j < N - 1; j++) {
			xi_list = xi(i, j, k);
            		K(pos(i,j),pos(i-1,j)) = -xi_list(0);
                	K(pos(i,j),pos(i+1,j)) = -xi_list(1);
                	K(pos(i,j),pos(i,j-1)) = -xi_list(2);
                	K(pos(i,j),pos(i,j+1)) = -xi_list(3);
            		K(pos(i,j),pos(i,j)) = xi_list(0)+xi_list(1)+xi_list(2)+xi_list(3);
		}
	}
	return K;
}






double k(int i, int j) {
	double y;
    	//double pctmetal = 0.4;
    	double pctmetal = 0.1;

    	if ( (j<=0.2*N) || (j > 0.8*N)) {
      	  	pctmetal = 0.9;
 		}

   	y = (65-0.2)*(pctmetal*pctmetal*pctmetal) + 0.2;
    
    	if ( (i<1) || (j<1) || (i>N) || (j>N) ) {
		y = 0;
		}

	return y;
}





}

#endif
