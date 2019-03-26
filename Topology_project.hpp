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
//#include <armadillo>


namespace tws {

// variabelen definiëren

int N = 10;

tws::matrix<double> x(N,N,0.4);
double penal = 3;

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
	int position;
	position = i+N*j;
	return position;
}

double meank(double k1, double k2) {
	double m;
	m = (k1+k2)/2; //arithmetic
    	//m = 2/(1/k1+1/k2); //harmonic
	return m;
}

tws::vector<int> create_dir_list() {
	int length_vector_dir = 2*(round(0.7*N) - round(0.3*N));
	int length_side = round(0.7*N) - round(0.3*N);
	tws::vector<int> Dir_positions(length_vector_dir,0.0);

	int iter = 0;

	for (int i = 0; i < length_side; i++) {
		Dir_positions(i) = round(0.3*N) + iter;
		iter += 1;
	}
	
	iter = 0;

	for (int i = length_side; i < length_vector_dir; i++) {
		Dir_positions(i) = pos(0,N-1) + round(0.3*N) + iter;
		iter += 1;
	}
	

	return Dir_positions;
}

tws::matrix<double> create_K(tws::matrix<double> k, tws::vector<int> dir_pos_list, tws::vector<double> f) {
	tws::matrix<double> K(N*N,N*N,0.0);
	bool is_Dir_pos = false;

	for (int i = 0; i < N; i++) {
    		for (int j = 0; j < N; j++) {
			int position = pos(i,j);
			if  ( std::find(std::begin(dir_pos_list), std::end(dir_pos_list), position) != std::end(dir_pos_list)) {
				is_Dir_pos = true;
				K(pos(i,j),pos(i,j)) = 1;
				f(pos(i,j)) = 293;
			}
			
			if (! is_Dir_pos) {
				std::cout<<pos(i,j)<<std::endl;
				double xin = 0, xis = 0, xiw = 0, xie = 0;
				if ( i > 1 ) {
					xin = meank(k(i,j),k(i-1,j));
            				K(pos(i,j),pos(i-1,j)) = -xin;
				}
				if ( i < N ) {
       					xis = meank(k(i,j),k(i+1,j));
					K(pos(i,j),pos(i+1,j)) = -xis;
    				}
    				if ( j > 1 ) {
        				xiw = meank(k(i,j),k(i,j-1));
					K(pos(i,j),pos(i,j-1)) = -xiw;
    				}
				if ( j < N ) {
        				xie = meank(k(i,j),k(i,j+1));
					K(pos(i,j),pos(i,j+1)) = -xie;
   				}
            			K(pos(i,j),pos(i,j)) = xin+xis+xiw+xie;
			}
			is_Dir_pos = false;
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

// Objective function and sensitivity analysis

tws::matrix<double> dc(tws::vector<double> U, tws::matrix<double> x) {
	double c = 0.0;
	tws::vector<double> Ue(4, 0.0);
	tws::matrix<double> dc(N-2, N-2, 0.0);
	tws::matrix<double> Stiff(4,4,0.);
	Stiff(0,0) = 2/3;
	Stiff(0,1) = -1/6;
	Stiff(0,2) = -1/3;
	Stiff(0,3) = -1/6;
	Stiff(1,0) = -1/6;
	Stiff(1,1) = 2/3;
	Stiff(1,2) = -1/6;
	Stiff(1,3) = -1/3;
	Stiff(2,0) = -1/3;
	Stiff(2,1) = -1/6;
	Stiff(2,2) = 2/3;
	Stiff(2,3) = -1/6;
	Stiff(3,0) = -1/6;
	Stiff(3,1) = -1/3;
	Stiff(3,2) = -1/6;
	Stiff(3,3) = 2/3;
	for (int ely = 1; ely < N-1; ely++) {
		for (int elx = 1; elx < N-1; elx++) {
			int n1 = N*(elx-1) + ely;
			int n2 = N*elx + ely;
			int n3 = N*(elx+1) + ely;
			Ue(0) = U(n1);
			Ue(1) = U(n2-1);
			Ue(2) = U(n3);
			Ue(3) = U(n2+1);
			c = c + (0.001 + 0.999*pow(x(ely,elx),penal))*inner_product(Ue, Stiff*Ue);
			dc(ely, elx) = -0.999*penal*pow(x(ely, elx), penal-1)*inner_product(Ue, Stiff*Ue);
		}
	}
	return dc;
}
}
#endif

