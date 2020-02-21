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

using namespace std;

namespace tws {

// variabelen definiëren

int N;
double penal = 3.0;
double q = 2.0/(0.01*0.01*0.001);
double h = 0.01/N;
double qhp = - q * pow(h,2.0)/4.0;
double qn = - q * pow(h,2.0)/2.0;
double qq = - q * pow(h,2.0);

tws::matrix<double> pctmetal(N,N,0.4);

tws::matrix<double> create_k(tws::matrix<double> pctmetal) {
	tws::matrix<double> k(N,N,0.0);
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			k(i,j) = (65-0.2)*pctmetal(i,j) + 0.2;
		}
	}
	return k;
}

tws::matrix<double> bigkmat(tws::matrix<double> k){
	tws::matrix<double> bigmat(2*N+1, 2*N+1, 0.0);
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
		bigmat(i, 2*N) = k((i-1)/2, 2*N);
		bigmat(2*N, i) = k(2*N, (i-1)/2);
	}
	return bigmat;
}

tws::vector<double> RL(){
	tws::vector<double> rl(pow(N+1,2),qq);
	
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

tws::matrix<double> LL(tws::matrix<double> bigkmat){
	tws::matrix<double> ll(pow(N+1,2), pow(N+1,2), 0.0);
	for (int i = 0.3*N; i < 0.7*N + 1; i++){
		ll(i,i) = 1.0;
		ll(pow(N+1,2)-i-1, pow(N+1,2)-i-1) = 1.0;
	}
	
	//Linksboven hoekpunt
	ll(0,0) = -0.5*(bigkmat(1,0)+bigkmat(0,1));
	ll(0,1) = 0.5*bigkmat(1,0);
	ll(0,N+1) = 0.5*bigkmat(0,1);
	
	//Linksonder hoekpunt
	ll(N,0) = -0.5*(bigkmat(2*N-1,0)+bigkmat(2*N-1,1));
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
		for (int i = 1; i <N; i++){
			ll(j*(N+1)+i,j*(N+1)+i) = -(bigkmat(2*i-1,2*j) + bigkmat(2*i+1,2*j) + bigkmat(2*i,2*j-1) + bigkmat(2*i,2*j+1));
			ll(j*(N+1)+i,j*(N+1)+i-1) = bigkmat(2*i-1,2*j);
			ll(j*(N+1)+i,j*(N+1)+i+1) = bigkmat(2*i+1,2*j);
			ll(j*(N+1)+i,(j-1)*(N+1)+i) = bigkmat(2*i,2*j-1);
			ll(j*(N+1)+i,(j+1)*(N+1)+i) = bigkmat(2*i,2*j+1);
		}
	}
	
	return ll;
}
	









/* int pos(int i, int j) {
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

tws::matrix<double> create_K(tws::matrix<double> k, tws::vector<int> dir_pos_list) {
	tws::matrix<double> K(N*N,N*N,0.0);
	bool is_Dir_pos = false;

	for (int i = 0; i < N; i++) {
    		for (int j = 0; j < N; j++) {
			int position = pos(i,j);
			if  ( std::find(std::begin(dir_pos_list), std::end(dir_pos_list), position) != std::end(dir_pos_list)) {
				is_Dir_pos = true;
				K(pos(i,j),pos(i,j)) = 1;
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

tws::vector<double> create_f(tws::vector<int> dir_pos_list, tws::vector<double> f) {
	for (int i = 0; i < N; i++) {
    		for (int j = 0; j < N; j++) {
			int position = pos(i,j);
			if  ( std::find(std::begin(dir_pos_list), std::end(dir_pos_list), position) != std::end(dir_pos_list)) {
				f(pos(i,j)) = 293;
			}
		}
	}
	return f;
}

double k(int i, int j) {
	double y;
    	double pctmetal = 0.4;
    	//double pctmetal = 0.1;

    	if ( (j<=0.2*N) || (j > 0.8*N)) {
      	  	pctmetal = 0.9;
 		}

   	y = (65-0.2)*(pctmetal*pctmetal*pctmetal) + 0.2;
    
    	if ( (i<1) || (j<1) || (i>N) || (j>N) ) {
		y = 0;
		}

	return y;
}



vector<double> gauss(tws::matrix<double> K, tws::vector<double> f) {
	int n = N*N;
	tws::matrix<double> A(N*N,N*N+1);
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			A(i,j) = K(i,j);
		}
	}
	for (int i = 0; i < n; i++){
		A(i,n) = f(i);
	}

    for (int i=0; i<n; i++) {
        // Search for maximum in this column
        double maxEl = abs(A(i,i));
        int maxRow = i;
        for (int k=i+1; k<n; k++) {
            if (abs(A(k,i)) > maxEl) {
                maxEl = abs(A(k,i));
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        for (int k=i; k<n+1;k++) {
            double tmp = A(maxRow, k);
            A(maxRow, k) = A(i,k);
            A(i,k) = tmp;
        }

        // Make all rows below this one 0 in current column
        for (int k=i+1; k<n; k++) {
            double c = -A(k,i)/A(i,i);
            for (int j=i; j<n+1; j++) {
                if (i==j) {
                    A(k,j) = 0;
                } else {
                    A(k,j) += c * A(i,j);
                }
            }
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    vector<double> x(n);
    for (int i=n-1; i>=0; i--) {
        x(i) = A(i,n)/A(i,i);
        for (int k=i-1;k>=0; k--) {
            A(k,n) -= A(k,i) * x(i);
        }
    }
    return x;
}

// Objective function and sensitivity analysis

tws::matrix<double> dcc(tws::vector<double> U, tws::matrix<double> x) {
	double c = 0.0;
	tws::vector<double> Ue(4, 0.0);
	tws::matrix<double> dc(N, N, 0.0);
	tws::matrix<double> Stiff(4,4,0.0);
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

//Optimality criteria

tws::matrix<double> xnew(tws::matrix<double> x, tws::matrix<double> dc){
	double l1 = 0.0;
	double l2 = 100000.0;
	double move = 0.2;
	tws::matrix<double> xnew = x;
	while (l2-l1 > 0.0001){
		double sum = 0.0;
		double lmid = 0.5 * (l2 + l1);
		for (int i = 0; i < N; i++){
			for (int j = 0; j < N; j++){
				xnew(i,j) = std::max(0.001, std::max(x(i,j) - move, std::min(1.0, std::min(x(i,j) + move, x(i,j) * std::sqrt(-dc(i,j)/lmid)))));
				sum = sum + xnew(i,j);
			}
		}
		if (sum > 0.4*N*N){
			l1 = lmid;
		}
		else{
			l2 = lmid;
		}
	}
	return xnew;
}


//Mesh-independency filter
tws::matrix<double> check(tws::matrix<double> x, tws::matrix<double> dc){
	tws::matrix<double> dcn(N, N, 0.0);
	for (int i = 1; i < N-1; i++){
		for (int j = 1; j < N-1; j++){
			double sum = 0.0;
			for (int k = std::max(i-1,2); k < std::min(i+1, N-1); k++){
				for (int l = std::max(j-1,2); l < std::min(j+1, N-1); l++){
					double fac = 1.2 - std::sqrt(pow((i-k),2)+pow((j-l),2));
					sum = sum + std::max(0.0,fac);
					dcn(j,i) = dcn(j,i) + std::max(0.0,fac)*dc(l,k)*x(l,k);
				}
			}
			dcn(j,i) = dcn(j,i)/(x(j,i)*sum);
		}
	}
	return dcn;
}

tws::matrix<double> solution(int N){
	int loop = 0;
	double change = 1.0;
	tws::matrix<double> K(N*N,N*N,0.0);
	double h = 0.01/N; 
	tws::vector<double> f(N*N,2*h*h*pow(10.0,6.0));
	tws::matrix<double> pctmetal(N,N,0.0);
	tws::matrix<double> k(N,N,0.0);
	tws::vector<int> list_dir(2*(round(0.7*N) - round(0.3*N)),0.0);
	pctmetal = tws::create_pctmetal();
	while (change > 0.01){
		loop = loop + 1;
		tws::matrix<double> xold = x;
		k = tws::create_k(pctmetal);
		//list_dir = tws::create_dir_list();
		//K = tws::create_K(k, list_dir);
		//f = tws::create_f(list_dir, f);
		//tws::vector<double> U = gauss(K, f);
		//tws::matrix<double> dc = dcc(U, x);
		//dc = check(x, dc);
		//x = xnew(x, dc);
		//double currmax = 0.0;
		//for (int i = 0; i < N; i++){
			//for (int j = 0; j < N; j++){
				//if (std::abs(x(i,j)-xold(i,j)) > currmax){
					//currmax = std::abs(x(i,j)-xold(i,j));
				//}
			//}
		//}
		//change = currmax;
	}
	return k;
} */

}
#endif


