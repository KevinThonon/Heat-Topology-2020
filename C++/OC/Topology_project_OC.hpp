#ifndef top_Topology_project_OC_hpp
#define top_Topology_project_OC_hpp
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


// Dichtheid k in elk element met SIMP methode. Input: percentage metaal in elk element. Output: Dichtheid k in elk element
mat create_k(vec a, int N, double penal) {
	mat k(N,N);
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			k(i,j) = (65-0.2)*pow(a(i + N*j),penal) + 0.2;
		}
	}
	return k;
}

// Grote matrix met ook k-waarden op randen van elementen (arithmisch gemiddelde)
mat bigkmat(mat k, int N){
	mat bigmat(2*N+1, 2*N+1);
	bigmat.fill(0.0);
	// Midden van elk element
	for (int i = 1; i < 2*N+1; i+=2){
		for (int j = 1; j < 2*N+1; j+=2){
			bigmat(i,j) = k((i-1)/2, (j-1)/2);
		}
	}
	
	// Randen van middelste elementen
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
	
	// k-waarden op de randen van het grid
	for (int i = 1; i < 2*N; i+=2){
		
		bigmat(i, 0) = k((i-1)/2, 0);
		bigmat(0, i) = k(0, (i-1)/2);
		bigmat(i, 2*N) = k((i-1)/2, N-1);
		bigmat(2*N, i) = k(N-1, (i-1)/2);
		
	}
	return bigmat;
}

// Rechterlid van Ku = f aanmaken. Output is een vector die het grid kolom per kolom bevat
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
	
	for (int i = N+1; i < N*(N+1); i+=(N+1)){
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

/* // Matrix K aanmaken
mat LL(mat& bigkmat, int N){
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
} */

double mean(double a, double b){
	//return (a+b)/2.0;	//arithmetic
	return 2.0*(a*b)/(a+b);		//harmonic
}



/* // Matrix K aanmaken zonder bigkmat
mat K_mat(mat& k, int N){
	//Aanmaken matrix en vullen met nullen
	mat ll(((N+1)*(N+1)), ((N+1)*(N+1)));
	ll.fill(0.0);
	//Dirichlett boundary conditions
	for (int i = 0.3*N; i < 0.7*N + 1; i++){
		ll(i,i) = 1.0;
		ll(((N+1)*(N+1))-i-1, ((N+1)*(N+1))-i-1) = 1.0;
	}
	
	//Linksboven hoekpunt
	ll(0,0) = -0.5*(k(0,0)+k(0,0));
	ll(0,1) = 0.5*k(0,0);
	ll(0,N+1) = 0.5*k(0,0);
	
	//Linksonder hoekpunt
	ll(N,N) = -0.5*(k(N-1,0)+k(N-1,0)); //index ll(N,0) aangepast
	ll(N,N-1) = 0.5*k(N-1,0);
	ll(N,2*N+1) = 0.5*k(N-1,0);
	
	//Rechtsboven hoekpunt
	ll(N*(N+1),N*(N+1)) = -0.5*(k(0,N-1)+k(0,N-1));
	ll(N*(N+1),N*(N+1)+1) = 0.5*k(0,N-1);
	ll(N*(N+1),(N-1)*(N+1)) = 0.5*k(0,N-1);
	
	//Rechtsonder hoekpunt
	ll(N*(N+2),N*(N+2)) = -0.5*(k(N-1,N-1)+k(N-1,N-1));
	ll(N*(N+2),N*(N+1)-1) = 0.5*k(N-1,N-1);
	ll(N*(N+2),N*(N+2)-1) = 0.5*k(N-1,N-1);
	
	//Neumann links
	for (int i = 1; i < 0.3*N; i++){
		ll(i,i) = -0.5*k(i-1,0) -0.5*k(i,0) - ((k(i-1,0) + k(i,0))/2.0);
		ll(i,i-1) = 0.5*k(i-1,0);
		ll(i,i+1) = 0.5*k(i,0);
		ll(i,i+N+1) = (k(i-1,0) + k(i,0))/2.0;
		
		ll(N-i,N-i) = -0.5*k(N-i,0) -0.5*k(N-1-i,0) - ((k(N-i,0) + k(N-1-i,0))/2.0);
		ll(N-i,N-i-1) = 0.5*k(N-1-i,0);
		ll(N-i,N-i+1) = 0.5*k(N-i,0);
		ll(N-i,N-i+N+1) = (k(N-i,0) + k(N-1-i,0))/2.0;
	}

	//Neumann rechts
	for (int i = 1; i < 0.3*N; i++){
		ll(N*(N+1)+i,N*(N+1)+i) = -0.5*k(i-1,N-1) -0.5*k(i,N-1) - ((k(i-1,N-1) + k(i,N-1))/2.0);
		ll(N*(N+1)+i,N*(N+1)+i-1) = 0.5*k(i-1,N-1);
		ll(N*(N+1)+i,N*(N+1)+i+1) = 0.5*k(i,N-1);
		ll(N*(N+1)+i,N*(N+1)+i-(N+1)) = (k(i-1,N-1) + k(i,N-1))/2.0;
		
		ll(N*(N+2)-i,N*(N+2)-i) = -0.5*k(N-i,N-1) -0.5*k(N-1-i,N-1) - ((k(N-i,N-1) + k(N-1-i,N-1))/2.0);
		ll(N*(N+2)-i,N*(N+2)-i-1) = 0.5*k(N-1-i,N-1);
		ll(N*(N+2)-i,N*(N+2)-i+1) = 0.5*k(N-i,N-1);
		ll(N*(N+2)-i,N*(N+2)-i-(N+1)) = (k(N-i,N-1) + k(N-1-i,N-1))/2.0;
	}
	

	for (int i = 1; i < N; i++){
	//Neumann boven
		ll(i*(N+1),i*(N+1)) = -0.5*k(0,i-1) -0.5*k(0,i) - ((k(0,i-1) + k(0,i))/2.0);
		ll(i*(N+1),i*(N+1)+1) = (k(0,i-1) + k(0,i))/2.0;
		ll(i*(N+1),(i-1)*(N+1)) = 0.5*k(0,i-1);
		ll(i*(N+1),(i+1)*(N+1)) = 0.5*k(0,i);
	
	//Neumann onder
		ll(i*(N+1)+N,i*(N+1)+N) = -0.5*k(N-1,i-1) -0.5*k(N-1,i) - ((k(N-1,i-1) + k(N-1,i))/2.0);
		ll(i*(N+1)+N,i*(N+1)+N-1) = (k(N-1,i-1) + k(N-1,i))/2.0;
		ll(i*(N+1)+N,(i-1)*(N+1)+N) = 0.5*k(N-1,i-1);
		ll(i*(N+1)+N,(i+1)*(N+1)+N) = 0.5*k(N-1,i);
	}
	
	//Gewone kolommen
	for (int j = 1; j < N; j++){
		for (int i = 1; i < N; i++){
			ll(j*(N+1)+i,j*(N+1)+i) = -((k(i-1,j-1)+k(i-1,j))/2.0 + (k(i,j-1)+k(i,j))/2.0 + (k(i-1,j-1)+k(i,j-1))/2.0 + (k(i-1,j)+k(i,j))/2.0);
			ll(j*(N+1)+i,j*(N+1)+i-1) = (k(i-1,j-1)+k(i-1,j))/2.0;
			ll(j*(N+1)+i,j*(N+1)+i+1) = (k(i,j-1)+k(i,j))/2.0;
			ll(j*(N+1)+i,(j-1)*(N+1)+i) = (k(i-1,j-1)+k(i,j-1))/2.0;
			ll(j*(N+1)+i,(j+1)*(N+1)+i) = (k(i-1,j)+k(i,j))/2.0;
		}
	}
	
	return ll;
} */

// Matrix K aanmaken zonder bigkmat met arithmetic/harmonic mean voor de k-waarden
sp_mat K_mat(mat& k, int N){
	//Aanmaken matrix en vullen met nullen
	sp_mat ll(((N+1)*(N+1)), ((N+1)*(N+1)));
	//Dirichlett boundary conditions
	for (int i = 0.3*N; i < 0.7*N + 1; i++){
		ll(i,i) = 1.0;
		ll(((N+1)*(N+1))-i-1, ((N+1)*(N+1))-i-1) = 1.0;
	}
	
	//Linksboven hoekpunt
	ll(0,0) = -1.0*k(0,0);
	ll(0,1) = 0.5*k(0,0);
	ll(0,N+1) = 0.5*k(0,0);
	
	//Linksonder hoekpunt
	ll(N,N) = -1.0*k(N-1,0); //index ll(N,0) aangepast
	ll(N,N-1) = 0.5*k(N-1,0);
	ll(N,2*N+1) = 0.5*k(N-1,0);
	
	//Rechtsboven hoekpunt
	ll(N*(N+1),N*(N+1)) = -1.0*k(0,N-1);
	ll(N*(N+1),N*(N+1)+1) = 0.5*k(0,N-1);
	ll(N*(N+1),(N-1)*(N+1)) = 0.5*k(0,N-1);
	
	//Rechtsonder hoekpunt
	ll(N*(N+2),N*(N+2)) = -1.0*k(N-1,N-1);
	ll(N*(N+2),N*(N+1)-1) = 0.5*k(N-1,N-1);
	ll(N*(N+2),N*(N+2)-1) = 0.5*k(N-1,N-1);
	
	//Neumann links
	for (int i = 1; i < 0.3*N; i++){
		ll(i,i) = -0.5*k(i-1,0) -0.5*k(i,0) - 1.0*mean(k(i-1,0),k(i,0));
		ll(i,i-1) = 0.5*k(i-1,0);
		ll(i,i+1) = 0.5*k(i,0);
		ll(i,i+N+1) = 1.0*mean(k(i-1,0),k(i,0));
		
		ll(N-i,N-i) = -0.5*k(N-i,0) -0.5*k(N-1-i,0) - 1.0*mean(k(N-i,0),k(N-1-i,0));
		ll(N-i,N-i-1) = 0.5*k(N-1-i,0);
		ll(N-i,N-i+1) = 0.5*k(N-i,0);
		ll(N-i,N-i+N+1) = 1.0*mean(k(N-i,0),k(N-1-i,0));
	}

	//Neumann rechts
	for (int i = 1; i < 0.3*N; i++){
		ll(N*(N+1)+i,N*(N+1)+i) = -0.5*k(i-1,N-1) -0.5*k(i,N-1) - 1.0*mean(k(i-1,N-1),k(i,N-1));
		ll(N*(N+1)+i,N*(N+1)+i-1) = 0.5*k(i-1,N-1);
		ll(N*(N+1)+i,N*(N+1)+i+1) = 0.5*k(i,N-1);
		ll(N*(N+1)+i,N*(N+1)+i-(N+1)) = 1.0*mean(k(i-1,N-1),k(i,N-1));
		
		ll(N*(N+2)-i,N*(N+2)-i) = -0.5*k(N-i,N-1) -0.5*k(N-1-i,N-1) - 1.0*mean(k(N-i,N-1),k(N-1-i,N-1));
		ll(N*(N+2)-i,N*(N+2)-i-1) = 0.5*k(N-1-i,N-1);
		ll(N*(N+2)-i,N*(N+2)-i+1) = 0.5*k(N-i,N-1);
		ll(N*(N+2)-i,N*(N+2)-i-(N+1)) = 1.0*mean(k(N-i,N-1),k(N-1-i,N-1));
	}
	

	for (int i = 1; i < N; i++){
	//Neumann boven
		ll(i*(N+1),i*(N+1)) = -0.5*k(0,i-1) -0.5*k(0,i) - 1.0*mean(k(0,i-1),k(0,i));
		ll(i*(N+1),i*(N+1)+1) = 1.0*mean(k(0,i-1),k(0,i));
		ll(i*(N+1),(i-1)*(N+1)) = 0.5*k(0,i-1);
		ll(i*(N+1),(i+1)*(N+1)) = 0.5*k(0,i);
	
	//Neumann onder
		ll(i*(N+1)+N,i*(N+1)+N) = -0.5*k(N-1,i-1) -0.5*k(N-1,i) - 1.0*mean(k(N-1,i-1),k(N-1,i));
		ll(i*(N+1)+N,i*(N+1)+N-1) = 1.0*mean(k(N-1,i-1),k(N-1,i));
		ll(i*(N+1)+N,(i-1)*(N+1)+N) = 0.5*k(N-1,i-1);
		ll(i*(N+1)+N,(i+1)*(N+1)+N) = 0.5*k(N-1,i);
	}
	
	//Gewone kolommen
	for (int j = 1; j < N; j++){
		for (int i = 1; i < N; i++){
			ll(j*(N+1)+i,j*(N+1)+i) = -1.0*(mean(k(i-1,j-1),k(i-1,j)) + mean(k(i,j-1),k(i,j)) + mean(k(i-1,j-1),k(i,j-1)) + mean(k(i-1,j),k(i,j)));
			ll(j*(N+1)+i,j*(N+1)+i-1) = mean(k(i-1,j-1),k(i-1,j));
			ll(j*(N+1)+i,j*(N+1)+i+1) = mean(k(i,j-1),k(i,j));
			ll(j*(N+1)+i,(j-1)*(N+1)+i) = mean(k(i-1,j-1),k(i,j-1));
			ll(j*(N+1)+i,(j+1)*(N+1)+i) = mean(k(i-1,j),k(i,j));
		}
	}
	
	return ll;
}

// Kostfunctie = u^T * u / #elementen
double objective_function1(vec& T, int N){
	double cost = dot(T,T)/(N*N);
	return cost;
}

// Kostfunctie = u / #elementen
// Deze kostfunctie is het gemiddelde van alle temperatuur-punten.
double objective_function2(vec& T, int N){
	vec w(pow((N+1),2));
	w.fill(1.0);
	double cost = dot(w,T)/(pow((N+1),2));
	return cost;
}

// Kostfunctie = w^T * u
// Deze kostfunctie is het gewogen gemiddelde van alle temperatuur-punten.
// w is het gewicht van elk temperatuur-punt. De som van alle gewichten is 1. Een middenpunt heeft gewicht 1/N² = h_x * h_y. Een randpunt heeft hiervan de helft (vakje is half zo klein).
// Een hoekpunt heeft gewicht h_x * h_y / 4 (vakje is een kwart van een middelste vakje).
double objective_function3(vec& T, int N){
	double A = pow((1.0/N),2);
	vec w(pow((N+1),2));
	w.fill(A);
	w(0) = A/4.0;
	w(N) = A/4.0;
	w(pow(N,2)+N) = A/4.0;
	w(pow(N,2)+2*N) = A/4.0;
	for (int i = 1; i < N; i++){
		w(i) = A/2.0;
		w(pow(N,2)+2*N-i) = A/2.0;
		w(i*(N+1)) = A/2.0;
		w(i*(N+1)+N) = A/2.0;
	}
	double cost = dot(w,T);
	return cost;
}

// Lambda is de oplossing van K^T * lambda = -dg/du = -2*u / #elementen
vec lambda1(vec& T, sp_mat& K, int N){
	vec lambda = spsolve(K.t(), (-2.0/(N*N))*T, "lapack"); //Veranderd van 2.0 naar -2.0
	return lambda;
}

// Lambda is de oplossing van K^T * lambda = -dg/du = - eenheidsvector / #elementen
vec lambda2(sp_mat& K, int N){
	vec w(pow((N+1),2));
	w.fill(-1.0);
	vec lambda = spsolve(K.t(), w, "lapack");
	return lambda;
}

// Lambda is de oplossing van K^T * lambda = -dg/du = -w
vec lambda3(sp_mat& K, int N){
	double A = pow((1.0/N),2);
	vec w(pow((N+1),2));
	w.fill(A);
	w(0) = A/4.0;
	w(N) = A/4.0;
	w(pow(N,2)+N) = A/4.0;
	w(pow(N,2)+2*N) = A/4.0;
	for (int i = 1; i < N; i++){
		w(i) = A/2.0;
		w(pow(N,2)+2*N-i) = A/2.0;
		w(i*(N+1)) = A/2.0;
		w(i*(N+1)+N) = A/2.0;
	}
	vec lambda = spsolve(K.t(), -1.0*w,"lapack");
	return lambda;
}


//Finite difference
vec dcda_fd(vec& T, vec& RL, vec& a, int N, double penal){
	vec dcda(N*N);
	dcda.fill(0.0);
	
	//cost function using original a
	double cost_original = objective_function1(T, N);
	
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			//Change element a_i,j a little bit (reduce percentage of metal in a_i,j with 0.01) 
			a(i + j*N) -= 0.01;
			//With updated a_i,j calculate k again, then K, then u
			mat k = create_k(a, N, penal);
			sp_mat ll = K_mat(k, N);
			vec u = spsolve(ll,RL,"lapack");
			//Use this updated u in the cost function
			double cost_update = objective_function1(u, N);
			//Use finite difference between original cost function and updated cost function where the difference between original and updated a_i,j = 0.01
			dcda(i + j*N) = (cost_original - cost_update)/0.01;
			a(i + j*N) += 0.01;
		}
	}
	
	return dcda;

}
			

// dc/da is een vector van gradienten die nodig is in de optimalisatie stap
// Deze dcda moet gebruikt worden indien met arithmetic mean wordt gewerkt
vec dcda_arit(vec lambda, vec T, vec& a, int N, double penal){
	//Initialiseren dc/da en opvullen met nullen
	vec dcda(N*N);
	dcda.fill(0.0);
	//Initialiseren dc/dk en opvullen met nullen
	vec dcdk(N*N);
	dcdk.fill(0.0);
	
	//dc/da element per element opvullen
	for (int i = 1; i < N-1; i++){
		for (int j = 1; j < N-1; j++){
			//Wat berekend moet worden is dc/da_i,j = dk/da * (lambda^T * dK/dk_i,j * u). Aangezien dK/dk_i,j slechts 4 niet-nul rijen zal bevatten zal niet eerst dK/dk_i,j aangemaakt worden om dan te vermenigvuldigen met u
			//maar zal dK/dk_i,j * u direct aangemaakt worden. 
			vec dKdk_u((N+1)*(N+1));
			dKdk_u.fill(0.0);
			dKdk_u(i + j*(N+1)) = -1.0*T(i + j*(N+1)) + 0.5*T(i+1 + j*(N+1)) + 0.5*T(i + (j+1)*(N+1));
			dKdk_u(i+1 + j*(N+1)) = -1.0*T(i+1 + j*(N+1)) + 0.5*T(i + j*(N+1)) + 0.5*T(i+1 + (j+1)*(N+1));
			dKdk_u(i + (j+1)*(N+1)) = -1.0*T(i + (j+1)*(N+1)) + 0.5*T(i + j*(N+1)) + 0.5*T(i+1 + (j+1)*(N+1));
			dKdk_u(i+1 + (j+1)*(N+1)) = -1.0*T(i+1 + (j+1)*(N+1)) + 0.5*T(i+1 + j*(N+1)) + 0.5*T(i + (j+1)*(N+1));
			//Vermenigvuldiging met lambda^T om tot dc/dk te komen
			dcdk(i + j*N) = dot(lambda, dKdk_u);
			//Vermenigvuldiging met dk/da om tot dc/da te komen
			dcda(i + j*N) = penal*(65.0-0.2)*pow(a(i + N*j),penal-1)*dcdk(i + j*N);
		}
	}
	
	//Linksboven
	vec dKdk_ulb((N+1)*(N+1));
	dKdk_ulb.fill(0.0);
	dKdk_ulb(0) = -1.0*T(0) + 0.5*T(1) + 0.5*T(N+1);
	dKdk_ulb(1) = -1.0*T(1) + 0.5*T(0) + 0.5*T(N+2);
	dKdk_ulb(N+1) = -1.0*T(N+1) + 0.5*T(0) + 0.5*T(N+2);
	dKdk_ulb(N+2) = -1.0*T(N+2) + 0.5*T(1) + 0.5*T(N+1);
	
	dcdk(0) = dot(lambda, dKdk_ulb);
	
	dcda(0) = penal*(65.0-0.2)*pow(a(0),penal-1)*dcdk(0);
	
	//Linksonder
	vec dKdk_ulo((N+1)*(N+1));
	dKdk_ulo.fill(0.0);
	dKdk_ulo(N) = -1.0*T(N) + 0.5*T(N-1) + 0.5*T(2*N+1);
	dKdk_ulo(N-1) = -1.0*T(N-1) + 0.5*T(N) + 0.5*T(2*N);
	dKdk_ulo(2*N+1) = -1.0*T(2*N+1) + 0.5*T(N) + 0.5*T(2*N);
	dKdk_ulo(2*N) = -1.0*T(2*N) + 0.5*T(N-1) + 0.5*T(2*N+1);
	
	dcdk(N-1) = dot(lambda, dKdk_ulo);
	
	dcda(N-1) = penal*(65.0-0.2)*pow(a(N-1),penal-1)*dcdk(N-1);
	
	//Rechtsboven
	vec dKdk_urb((N+1)*(N+1));
	dKdk_urb.fill(0.0);
	dKdk_urb(N*N+N) = -1.0*T(N*N+N) + 0.5*T(N*N+N+1) + 0.5*T(N*N-1);
	dKdk_urb(N*N-1) = -1.0*T(N*N-1) + 0.5*T(N*N+N) + 0.5*T(N*N);
	dKdk_urb(N*N+N+1) = -1.0*T(N*N+N+1) + 0.5*T(N*N+N) + 0.5*T(N*N);
	dKdk_urb(N*N) = -1.0*T(N*N) + 0.5*T(N*N-1) + 0.5*T(N*N+N+1);
	
	dcdk((N-1)*N) = dot(lambda, dKdk_urb);
	
	dcda((N-1)*N) = penal*(65.0-0.2)*pow(a((N-1)*N),penal-1)*dcdk((N-1)*N);
	
	//Rechtsonder
	vec dKdk_uro((N+1)*(N+1));
	dKdk_uro.fill(0.0);
	dKdk_uro(N*N + 2*N) = -1.0*T(N*N + 2*N) + 0.5*T(N*N + 2*N - 1) + 0.5*T(N*N + N - 1);
	dKdk_uro(N*N + 2*N - 1) = -1.0*T(N*N + 2*N - 1) + 0.5*T(N*N + 2*N) + 0.5*T(N*N + N - 2);
	dKdk_uro(N*N + N - 1) = -1.0*T(N*N + N - 1) + 0.5*T(N*N + 2*N) + 0.5*T(N*N + N - 2);
	dKdk_uro(N*N + N - 2) = -1.0*T(N*N + N - 2) + 0.5*T(N*N + 2*N - 1) + 0.5*T(N*N + N - 1);
	
	dcdk(N*N-1) = dot(lambda, dKdk_uro);
	
	dcda(N*N-1) = penal*(65.0-0.2)*pow(a(N*N-1),penal-1)*dcdk(N*N-1);
	
	//Bovenrand en onderrand
	for (int j = 1; j < N-1; j++){
		vec dKdk_ub((N+1)*(N+1));
		dKdk_ub.fill(0.0);
		dKdk_ub(j*(N+1)) = -1.0*T(j*(N+1)) + 0.5*T(1 + j*(N+1)) + 0.5*T((j+1)*(N+1));
		dKdk_ub(1 + j*(N+1)) = -1.0*T(1 + j*(N+1)) + 0.5*T(j*(N+1)) + 0.5*T(1 + (j+1)*(N+1));
		dKdk_ub((j+1)*(N+1)) = -1.0*T((j+1)*(N+1)) + 0.5*T(j*(N+1)) + 0.5*T(1 + (j+1)*(N+1));
		dKdk_ub(1 + (j+1)*(N+1)) = -1.0*T(1 + (j+1)*(N+1)) + 0.5*T(1 + j*(N+1)) + 0.5*T((j+1)*(N+1));
	
		dcdk(j*N) = dot(lambda, dKdk_ub);
	
		dcda(j*N) = penal*(65.0-0.2)*pow(a(j*N),penal-1)*dcdk(j*N);
		
		
		vec dKdk_uo((N+1)*(N+1));
		dKdk_uo.fill(0.0);
		dKdk_uo(j*(N+1)) = -1.0*T(N-1 + j*(N+1)) + 0.5*T(N + j*(N+1)) + 0.5*T(N-1 + (j+1)*(N+1));
		dKdk_uo(1 + j*(N+1)) = -1.0*T(N + j*(N+1)) + 0.5*T(N-1 + j*(N+1)) + 0.5*T(N + (j+1)*(N+1));
		dKdk_uo((j+1)*(N+1)) = -1.0*T(N-1 + (j+1)*(N+1)) + 0.5*T(N-1 + j*(N+1)) + 0.5*T(N + (j+1)*(N+1));
		dKdk_uo(1 + (j+1)*(N+1)) = -1.0*T(N + (j+1)*(N+1)) + 0.5*T(N + j*(N+1)) + 0.5*T(N-1 + (j+1)*(N+1));
	
		dcdk(N-1 + j*N) = dot(lambda, dKdk_uo);
	
		dcda(N-1 + j*N) = penal*(65.0-0.2)*pow(a(N-1 + N*j),penal-1)*dcdk(N-1 + j*N);
	}
	
	//Neamann links
	for (int i = 1; i < 0.3*N - 1; i++){
		vec dKdk_unlb((N+1)*(N+1));
		dKdk_unlb.fill(0.0);
		dKdk_unlb(i) = -1.0*T(i) + 0.5*T(i+1) + 0.5*T(i + (N+1));
		dKdk_unlb(i+1) = -1.0*T(i+1) + 0.5*T(i) + 0.5*T(i+1 + (N+1));
		dKdk_unlb(i + (N+1)) = -1.0*T(i + (N+1)) + 0.5*T(i) + 0.5*T(i+1 + (N+1));
		dKdk_unlb(i+1 + (N+1)) = -1.0*T(i+1 + (N+1)) + 0.5*T(i+1) + 0.5*T(i + (N+1));

		dcdk(i) = dot(lambda, dKdk_unlb);

		dcda(i) = penal*(65.0-0.2)*pow(a(i),penal-1)*dcdk(i);
		
			
		vec dKdk_unlo((N+1)*(N+1));
		dKdk_unlo.fill(0.0);
		dKdk_unlo(N-i) = -1.0*T(N-i) + 0.5*T(N-i-1) + 0.5*T(N-i + (N+1));
		dKdk_unlo(N-i-1) = -1.0*T(N-i-1) + 0.5*T(N-i) + 0.5*T(N-i-1 + (N+1));
		dKdk_unlo(N-i + (N+1)) = -1.0*T(N-i + (N+1)) + 0.5*T(N-i) + 0.5*T(N-i-1 + (N+1));
		dKdk_unlo(N-i-1 + (N+1)) = -1.0*T(N-i-1 + (N+1)) + 0.5*T(N-i-1) + 0.5*T(N-i + (N+1));

		dcdk(N-i-1) = dot(lambda, dKdk_unlo);

		dcda(N-i-1) = penal*(65.0-0.2)*pow(a(N-i-1),penal-1)*dcdk(N-i-1);
	}
	
	//Neamann rechts
	for (int i = 1; i < 0.3*N - 1; i++){
		vec dKdk_unrb((N+1)*(N+1));
		dKdk_unrb.fill(0.0);
		dKdk_unrb(i + N*(N+1)) = -1.0*T(i + N*(N+1)) + 0.5*T(i+1 + N*(N+1)) + 0.5*T(i + (N-1)*(N+1));
		dKdk_unrb(i+1 + N*(N+1)) = -1.0*T(i+1 + N*(N+1)) + 0.5*T(i + N*(N+1)) + 0.5*T(i+1 + (N-1)*(N+1));
		dKdk_unrb(i + (N-1)*(N+1)) = -1.0*T(i + (N-1)*(N+1)) + 0.5*T(i + N*(N+1)) + 0.5*T(i+1 + (N-1)*(N+1));
		dKdk_unrb(i+1 + (N-1)*(N+1)) = -1.0*T(i+1 + (N-1)*(N+1)) + 0.5*T(i + (N-1)*(N+1)) + 0.5*T(i+1 + N*(N+1));

		dcdk(N*N-N+i) = dot(lambda, dKdk_unrb);

		dcda(N*N-N+i) = penal*(65.0-0.2)*pow(a(N*N-N+i),penal-1)*dcdk(N*N-N+i);
		
		
		vec dKdk_unro((N+1)*(N+1));
		dKdk_unro.fill(0.0);
		dKdk_unro(N-i + N*(N+1)) = -1.0*T(N-i + N*(N+1)) + 0.5*T(N-i-1 + N*(N+1)) + 0.5*T(N-i + (N-1)*(N+1));
		dKdk_unro(N-i-1 + N*(N+1)) = -1.0*T(N-i-1 + N*(N+1)) + 0.5*T(N-i + N*(N+1)) + 0.5*T(N-i-1 + (N-1)*(N+1));
		dKdk_unro(N-i + (N-1)*(N+1)) = -1.0*T(N-i + (N-1)*(N+1)) + 0.5*T(N-i + N*(N+1)) + 0.5*T(N-i-1 + (N-1)*(N+1));
		dKdk_unro(N-i-1 + (N-1)*(N+1)) = -1.0*T(N-i-1 + (N-1)*(N+1)) + 0.5*T(N-i + (N-1)*(N+1)) + 0.5*T(N-i-1 + N*(N+1));

		dcdk(N*N-1-i) = dot(lambda, dKdk_unro);

		dcda(N*N-1-i) = penal*(65.0-0.2)*pow(a(N*N-1-i),penal-1)*dcdk(N*N-1-i);
	}
	
	//Vakje met 1 hoekpunt van 293K, linkerkant, boven dirichlett boundary condition
	vec dKdk_uhlb((N+1)*(N+1));
	dKdk_uhlb.fill(0.0);
	dKdk_uhlb(0.3*N - 1) = -1.0*T(0.3*N - 1) + 0.5*T(0.3*N) + 0.5*T(0.3*N - 1 + (N+1));
	//dKdk_uhlb(0.3*N) = -1.0*T(0.3*N) + 0.5*T(0.3*N - 1) + 0.5*T(0.3*N + (N+1));
	dKdk_uhlb(0.3*N - 1 + (N+1)) = -1.0*T(0.3*N - 1 + (N+1)) + 0.5*T(0.3*N - 1) + 0.5*T(0.3*N + (N+1));
	dKdk_uhlb(0.3*N + (N+1)) = -1.0*T(0.3*N + (N+1)) + 0.5*T(0.3*N) + 0.5*T(0.3*N - 1 + (N+1));

	dcdk(0.3*N - 1) = dot(lambda, dKdk_uhlb);

	dcda(0.3*N - 1) = penal*(65.0-0.2)*pow(a(0.3*N - 1),penal-1)*dcdk(0.3*N - 1);
	
	//Vakje met 1 hoekpunt van 293K, linkerkant, onder dirichlett boundary condition
	vec dKdk_uhlo((N+1)*(N+1));
	dKdk_uhlo.fill(0.0);
	//dKdk_uhlo(0.7*N) = -1.0*T(0.7*N) + 0.5*T(0.7*N + 1) + 0.5*T(0.7*N + (N+1));
	dKdk_uhlo(0.7*N + 1) = -1.0*T(0.7*N + 1) + 0.5*T(0.7*N) + 0.5*T(0.7*N + 1 + (N+1));
	dKdk_uhlo(0.7*N + (N+1)) = -1.0*T(0.7*N + (N+1)) + 0.5*T(0.7*N) + 0.5*T(0.7*N + 1 + (N+1));
	dKdk_uhlo(0.7*N + 1 + (N+1)) = -1.0*T(0.7*N + 1 + (N+1)) + 0.5*T(0.7*N + 1) + 0.5*T(0.7*N + (N+1));

	dcdk(0.7*N) = dot(lambda, dKdk_uhlo);

	dcda(0.7*N) = penal*(65.0-0.2)*pow(a(0.7*N),penal-1)*dcdk(0.7*N);
	
	//Vakjes direct naast de linkse dirichlett boundary condition. Deze hebben 2 hoekpunten van 293K
	for (int i = 0.3*N; i < 0.7*N; i++){
		vec dKdk_udl((N+1)*(N+1));
		dKdk_udl.fill(0.0);
		//dKdk_udl(i) = -1.0*T(i) + 0.5*T(i+1) + 0.5*T(i + (N+1));
		//dKdk_udl(i+1) = -1.0*T(i+1) + 0.5*T(i) + 0.5*T(i+1 + (N+1));
		dKdk_udl(i + (N+1)) = -1.0*T(i + (N+1)) + 0.5*T(i) + 0.5*T(i+1 + (N+1));
		dKdk_udl(i+1 + (N+1)) = -1.0*T(i+1 + (N+1)) + 0.5*T(i+1) + 0.5*T(i + (N+1));

		dcdk(i) = dot(lambda, dKdk_udl);

		dcda(i) = penal*(65.0-0.2)*pow(a(i),penal-1)*dcdk(i);
	}
	
	//Vakje met 1 hoekpunt van 293K, rechterkant, boven dirichlett boundary condition
	vec dKdk_uhrb((N+1)*(N+1));
	dKdk_uhrb.fill(0.0);
	dKdk_uhrb(N*N-1 + 0.3*N - 1) = -1.0*T(N*N-1 + 0.3*N - 1) + 0.5*T(N*N-1 + 0.3*N) + 0.5*T(N*N+N + 0.3*N - 1);
	dKdk_uhrb(N*N-1 + 0.3*N) = -1.0*T(N*N-1 + 0.3*N) + 0.5*T(N*N-1 + 0.3*N - 1) + 0.5*T(N*N+N + 0.3*N);
	dKdk_uhrb(N*N+N + 0.3*N - 1) = -1.0*T(N*N+N + 0.3*N - 1) + 0.5*T(N*N-1 + 0.3*N - 1) + 0.5*T(N*N+N + 0.3*N);
	//dKdk_uhrb(N*N+N + 0.3*N) = -1.0*T(N*N+N + 0.3*N) + 0.5*T(N*N-1 + 0.3*N) + 0.5*T(N*N+N + 0.3*N - 1);

	dcdk(N*(N-1) + 0.3*N - 1) = dot(lambda, dKdk_uhrb);

	dcda(N*(N-1) + 0.3*N - 1) = penal*(65.0-0.2)*pow(a(N*(N-1) + 0.3*N - 1),penal-1)*dcdk(N*(N-1) + 0.3*N - 1);
	
	//Vakje met 1 hoekpunt van 293K, rechterkant, onder dirichlett boundary condition
	vec dKdk_uhro((N+1)*(N+1));
	dKdk_uhro.fill(0.0);
	dKdk_uhro(N*N-1 + 0.7*N) = -1.0*T(N*N-1 + 0.7*N) + 0.5*T(N*N-1 + 0.7*N + 1) + 0.5*T(N*N+N + 0.7*N);
	dKdk_uhro(N*N-1 + 0.7*N + 1) = -1.0*T(N*N-1 + 0.7*N + 1) + 0.5*T(N*N-1 + 0.7*N) + 0.5*T(N*N+N + 0.7*N + 1);
	//dKdk_uhro(N*N+N + 0.7*N) = -1.0*T(N*N+N + 0.7*N) + 0.5*T(N*N-1 + 0.7*N) + 0.5*T(N*N+N + 0.7*N + 1);
	dKdk_uhro(N*N+N + 0.7*N + 1) = -1.0*T(N*N+N + 0.7*N + 1) + 0.5*T(N*N-1 + 0.7*N + 1) + 0.5*T(N*N+N + 0.7*N);

	dcdk(N*(N-1) + 0.7*N) = dot(lambda, dKdk_uhro);

	dcda(N*(N-1) + 0.7*N) = penal*(65.0-0.2)*pow(a(N*(N-1) + 0.7*N),penal-1)*dcdk(N*(N-1) + 0.7*N);
	
	//Vakjes direct naast de rechtse dirichlett boundary condition. Deze hebben 2 hoekpunten van 293K
	for (int i = 0.3*N; i < 0.7*N; i++){
		vec dKdk_udr((N+1)*(N+1));
		dKdk_udr.fill(0.0);
		dKdk_udr(N*N-1 + i) = -1.0*T(N*N-1 + i) + 0.5*T(N*N-1 + i+1) + 0.5*T(N*N+N + i);
		dKdk_udr(N*N-1 + i+1) = -1.0*T(N*N-1 + i+1) + 0.5*T(N*N-1 + i) + 0.5*T(N*N+N + i+1);
		//dKdk_udr(N*N+N + i) = -1.0*T(N*N+N + i) + 0.5*T(N*N-1 + i) + 0.5*T(N*N+N + i+1);
		//dKdk_udr(N*N+N + i+1) = -1.0*T(N*N+N + i+1) + 0.5*T(N*N-1 + i+1) + 0.5*T(N*N+N + i);

		dcdk(N*(N-1) + i) = dot(lambda, dKdk_udr);

		dcda(N*(N-1) + i) = penal*(65.0-0.2)*pow(a(N*(N-1) + i),penal-1)*dcdk(N*(N-1) + i);
	}
	
	return dcda;
	
}

// dc/da is een vector van gradienten die nodig zijn in de optimalisatie stap
// Deze dcda moet gebruikt worden indien met harmonic mean wordt gewerkt
vec dcda_harm(vec& lambda, vec& T, vec& a, mat& k, int N, double penal){
	//Initialiseren dc/da en opvullen met nullen
	vec dcda(N*N);
	dcda.fill(0.0);
	//Initialiseren dc/dk en opvullen met nullen
	vec dcdk(N*N);
	dcdk.fill(0.0);
	
	//dc/da element per element opvullen
	for (int i = 1; i < N-1; i++){
		for (int j = 1; j < N-1; j++){
			//Wat berekend moet worden is dc/da_i,j = dk/da * (lambda^T * dK/dk_i,j * u). Aangezien dK/dk_i,j slechts 4 niet-nul rijen zal bevatten zal niet eerst dK/dk_i,j aangemaakt worden om dan te vermenigvuldigen met u
			//maar zal dK/dk_i,j * u direct aangemaakt worden. 
			vec dKdk_u((N+1)*(N+1));
			dKdk_u.fill(0.0);
			dKdk_u(i + j*(N+1)) = (-2.0*(pow(k(i,j-1),2)/pow((k(i,j-1) + k(i,j)),2) + pow(k(i-1,j),2)/pow((k(i-1,j) + k(i,j)),2)))*T(i + j*(N+1)) + (2.0*(pow(k(i,j-1),2)/pow((k(i,j-1) + k(i,j)),2)))*T(i+1 + j*(N+1)) + (2.0*pow(k(i-1,j),2)/pow((k(i-1,j) + k(i,j)),2))*T(i + (j+1)*(N+1));
			dKdk_u(i+1 + j*(N+1)) = (-2.0*(pow(k(i,j-1),2)/pow((k(i,j-1) + k(i,j)),2) + pow(k(i+1,j),2)/pow((k(i+1,j) + k(i,j)),2)))*T(i+1 + j*(N+1)) + (2.0*(pow(k(i,j-1),2)/pow((k(i,j-1) + k(i,j)),2)))*T(i + j*(N+1)) + (2.0*pow(k(i+1,j),2)/pow((k(i+1,j) + k(i,j)),2))*T(i+1 + (j+1)*(N+1));
			dKdk_u(i + (j+1)*(N+1)) = (-2.0*(pow(k(i-1,j),2)/pow((k(i-1,j) + k(i,j)),2) + pow(k(i,j+1),2)/pow((k(i,j+1) + k(i,j)),2)))*T(i + (j+1)*(N+1)) + (2.0*(pow(k(i-1,j),2)/pow((k(i-1,j) + k(i,j)),2)))*T(i + j*(N+1)) + (2.0*pow(k(i,j+1),2)/pow((k(i,j+1) + k(i,j)),2))*T(i+1 + (j+1)*(N+1));
			dKdk_u(i+1 + (j+1)*(N+1)) = (-2.0*(pow(k(i+1,j),2)/pow((k(i+1,j) + k(i,j)),2) + pow(k(i,j+1),2)/pow((k(i,j+1) + k(i,j)),2)))*T(i+1 + (j+1)*(N+1)) + (2.0*(pow(k(i+1,j),2)/pow((k(i+1,j) + k(i,j)),2)))*T(i+1 + j*(N+1)) + (2.0*(pow(k(i,j+1),2)/pow((k(i,j+1) + k(i,j)),2)))*T(i + (j+1)*(N+1));
			//Vermenigvuldiging met lambda^T om tot dc/dk te komen
			dcdk(i + j*N) = dot(lambda, dKdk_u);
			//Vermenigvuldiging met dk/da om tot dc/da te komen
			dcda(i + j*N) = penal*(65.0-0.2)*pow(a(i + N*j),penal-1)*dcdk(i + j*N);
		}
	}
	
	//Bovenrand
	for (int j = 1; j < N-1; j++){
		vec dKdk_u((N+1)*(N+1));
		dKdk_u.fill(0.0);
		dKdk_u(j*(N+1)) = (-2.0*(pow(k(0,j-1),2)/pow((k(0,j-1) + k(0,j)),2)) - 1.0)*T(j*(N+1)) + (2.0*(pow(k(0,j-1),2)/pow((k(0,j-1) + k(0,j)),2)))*T(1 + j*(N+1)) + 1.0*T((j+1)*(N+1));
		dKdk_u(1 + j*(N+1)) = (-2.0*(pow(k(0,j-1),2)/pow((k(0,j-1) + k(0,j)),2) + pow(k(0+1,j),2)/pow((k(0+1,j) + k(0,j)),2)))*T(1 + j*(N+1)) + (2.0*(pow(k(0,j-1),2)/pow((k(0,j-1) + k(0,j)),2)))*T(j*(N+1)) + (2.0*pow(k(0+1,j),2)/pow((k(0+1,j) + k(0,j)),2))*T(1 + (j+1)*(N+1));
		dKdk_u((j+1)*(N+1)) = (-1.0 - 2.0*(pow(k(0,j+1),2)/pow((k(0,j+1) + k(0,j)),2)))*T((j+1)*(N+1)) + 1.0*T(j*(N+1)) + (2.0*pow(k(0,j+1),2)/pow((k(0,j+1) + k(0,j)),2))*T(1 + (j+1)*(N+1));
		dKdk_u(1 + (j+1)*(N+1)) = (-2.0*(pow(k(0+1,j),2)/pow((k(0+1,j) + k(0,j)),2) + pow(k(0,j+1),2)/pow((k(0,j+1) + k(0,j)),2)))*T(1 + (j+1)*(N+1)) + (2.0*(pow(k(0+1,j),2)/pow((k(0+1,j) + k(0,j)),2)))*T(1 + j*(N+1)) + (2.0*(pow(k(0,j+1),2)/pow((k(0,j+1) + k(0,j)),2)))*T((j+1)*(N+1));

		dcdk(j*N) = dot(lambda, dKdk_u);

		dcda(j*N) = penal*(65.0-0.2)*pow(a(N*j),penal-1)*dcdk(j*N);
	}
	
	//Onderste rand
	for (int j = 1; j < N-1; j++){
		vec dKdk_u((N+1)*(N+1));
		dKdk_u.fill(0.0);
		dKdk_u(N-1 + j*(N+1)) = (-2.0*(pow(k(N-1,j-1),2)/pow((k(N-1,j-1) + k(N-1,j)),2) + pow(k(N-2,j),2)/pow((k(N-2,j) + k(N-1,j)),2)))*T(N-1 + j*(N+1)) + (2.0*(pow(k(N-1,j-1),2)/pow((k(N-1,j-1) + k(N-1,j)),2)))*T(N + j*(N+1)) + (2.0*pow(k(N-2,j),2)/pow((k(N-2,j) + k(N-1,j)),2))*T(N-1 + (j+1)*(N+1));
		dKdk_u(N + j*(N+1)) = (-2.0*(pow(k(N-1,j-1),2)/pow((k(N-1,j-1) + k(N-1,j)),2)) - 1.0)*T(N + j*(N+1)) + (2.0*(pow(k(N-1,j-1),2)/pow((k(N-1,j-1) + k(N-1,j)),2)))*T(N-1 + j*(N+1)) + 1.0*T(N + (j+1)*(N+1));
		dKdk_u(N-1 + (j+1)*(N+1)) = (-2.0*(pow(k(N-2,j),2)/pow((k(N-2,j) + k(N-1,j)),2) + pow(k(N-1,j+1),2)/pow((k(N-1,j+1) + k(N-1,j)),2)))*T(N-1 + (j+1)*(N+1)) + (2.0*(pow(k(N-2,j),2)/pow((k(N-2,j) + k(N-1,j)),2)))*T(N-1 + j*(N+1)) + (2.0*pow(k(N-1,j+1),2)/pow((k(N-1,j+1) + k(N-1,j)),2))*T(N + (j+1)*(N+1));
		dKdk_u(N + (j+1)*(N+1)) = (-1.0 - 2.0*(pow(k(N-1,j+1),2)/pow((k(N-1,j+1) + k(N-1,j)),2)))*T(N + (j+1)*(N+1)) + 1.0*T(N + j*(N+1)) + (2.0*(pow(k(N-1,j+1),2)/pow((k(N-1,j+1) + k(N-1,j)),2)))*T(N-1 + (j+1)*(N+1));

		dcdk(N-1 + j*N) = dot(lambda, dKdk_u);

		dcda(N-1 + j*N) = penal*(65.0-0.2)*pow(a(N-1 + N*j),penal-1)*dcdk(N-1 + j*N);
	}
	
	//Linksboven hoekpunt
	vec dKdk_ulb((N+1)*(N+1));
	dKdk_ulb.fill(0.0);
	dKdk_ulb(0) = -2.0*T(0) + 1.0*T(1) + 1.0*T(N+1);
	dKdk_ulb(1) = (-1.0 - 2.0*(pow(k(1,0),2)/pow((k(1,0) + k(0,0)),2)))*T(1) + 1.0*T(0) + (2.0*pow(k(1,0),2)/pow((k(1,0) + k(0,0)),2))*T(N+2);
	dKdk_ulb(N+1) = (-1.0 - 2.0*(pow(k(0,1),2)/pow((k(0,1) + k(0,0)),2)))*T(N+1) + 1.0*T(0) + (2.0*pow(k(0,1),2)/pow((k(0,1) + k(0,0)),2))*T(N+2);
	dKdk_ulb(N+2) = (-2.0*(pow(k(1,0),2)/pow((k(1,0) + k(0,0)),2) + pow(k(0,1),2)/pow((k(0,1) + k(0,0)),2)))*T(N+2) + (2.0*(pow(k(1,0),2)/pow((k(1,0) + k(0,0)),2)))*T(1) + (2.0*(pow(k(0,1),2)/pow((k(0,1) + k(0,0)),2)))*T(N+1);

	dcdk(0) = dot(lambda, dKdk_ulb);

	dcda(0) = penal*(65.0-0.2)*pow(a(0),penal-1)*dcdk(0);
	
	//Linksonder hoekpunt
	vec dKdk_ulo((N+1)*(N+1));
	dKdk_ulo.fill(0.0);
	dKdk_ulo(N-1) = (-1.0 - 2.0*(pow(k(N-2,0),2)/pow((k(N-2,0) + k(N-1,0)),2)))*T(N-1) + 1.0*T(N) + (2.0*pow(k(N-2,0),2)/pow((k(N-2,0) + k(N-1,0)),2))*T(2*N);
	dKdk_ulo(N) = -2.0*T(N) + 1.0*T(N-1) + 1.0*T(2*N+1);
	dKdk_ulo(2*N) = (-2.0*(pow(k(N-2,0),2)/pow((k(N-2,0) + k(N-1,0)),2) + pow(k(N-1,1),2)/pow((k(N-1,1) + k(N-1,0)),2)))*T(2*N) + (2.0*(pow(k(N-2,0),2)/pow((k(N-2,0) + k(N-1,0)),2)))*T(N-1) + (2.0*pow(k(N-1,1),2)/pow((k(N-1,1) + k(N-1,0)),2))*T(2*N+1);
	dKdk_ulo(2*N+1) = (-1.0 - 2.0*(pow(k(N-1,1),2)/pow((k(N-1,1) + k(N-1,0)),2)))*T(2*N+1) + 1.0*T(N) + (2.0*(pow(k(N-1,1),2)/pow((k(N-1,1) + k(N-1,0)),2)))*T(2*N);

	dcdk(N-1) = dot(lambda, dKdk_ulo);

	dcda(N-1) = penal*(65.0-0.2)*pow(a(N-1),penal-1)*dcdk(N-1);
	
	//Rechtsboven hoekpunt
	vec dKdk_urb((N+1)*(N+1));
	dKdk_urb.fill(0.0);
	dKdk_urb(N*N - 1) = (-2.0*(pow(k(0,N-2),2)/pow((k(0,N-2) + k(0,N-1)),2)) - 1.0)*T(N*N - 1) + (2.0*(pow(k(0,N-2),2)/pow((k(0,N-2) + k(0,N-1)),2)))*T(N*N) + 1.0*T(N*(N+1));
	dKdk_urb(N*N) = (-2.0*(pow(k(0,N-2),2)/pow((k(0,N-2) + k(0,N-1)),2) + pow(k(1,N-1),2)/pow((k(1,N-1) + k(0,N-1)),2)))*T(N*N) + (2.0*(pow(k(0,N-2),2)/pow((k(0,N-2) + k(0,N-1)),2)))*T(N*N - 1) + (2.0*pow(k(1,N-1),2)/pow((k(1,N-1) + k(0,N-1)),2))*T(N*(N+1) + 1);
	dKdk_urb(N*(N+1)) = -2.0*T(N*(N+1)) + 1.0*T(N*N - 1) + 1.0*T(N*(N+1) + 1);
	dKdk_urb(N*(N+1) + 1) = (-2.0*(pow(k(1,N-1),2)/pow((k(1,N-1) + k(0,N-1)),2)) - 1.0)*T(N*(N+1) + 1) + (2.0*(pow(k(1,N-1),2)/pow((k(1,N-1) + k(0,N-1)),2)))*T(N*N) + 1.0*T(N*(N+1));
	
	dcdk((N-1)*N) = dot(lambda, dKdk_urb);

	dcda((N-1)*N) = penal*(65.0-0.2)*pow(a((N-1)*N),penal-1)*dcdk((N-1)*N);
	
	//Rechtsonder hoekpunt
	vec dKdk_uro((N+1)*(N+1));
	dKdk_uro.fill(0.0);
	dKdk_uro(N*N + N - 2) = (-2.0*(pow(k(N-1,N-2),2)/pow((k(N-1,N-2) + k(N-1,N-1)),2) + pow(k(N-2,N-1),2)/pow((k(N-2,N-1) + k(N-1,N-1)),2)))*T(N*N + N - 2) + (2.0*(pow(k(N-1,N-2),2)/pow((k(N-1,N-2) + k(N-1,N-1)),2)))*T(N*N + N - 1) + (2.0*pow(k(N-2,N-1),2)/pow((k(N-2,N-1) + k(N-1,N-1)),2))*T(N*N + 2*N - 1);
	dKdk_uro(N*N + N - 1) = (-2.0*(pow(k(N-1,N-2),2)/pow((k(N-1,N-2) + k(N-1,N-1)),2)) - 1.0)*T(N*N + N - 1) + (2.0*(pow(k(N-1,N-2),2)/pow((k(N-1,N-2) + k(N-1,N-1)),2)))*T(N*N + N - 2) + 1.0*T(N*N + 2*N);
	dKdk_uro(N*N + 2*N - 1) = (-2.0*(pow(k(N-2,N-1),2)/pow((k(N-2,N-1) + k(N-1,N-1)),2)) - 1.0)*T(N*N + 2*N - 1) + (2.0*(pow(k(N-2,N-1),2)/pow((k(N-2,N-1) + k(N-1,N-1)),2)))*T(N*N + N - 2) + 1.0*T(N*N + 2*N);
	dKdk_uro(N*N + 2*N) = -2.0*T(N*N + 2*N) + 1.0*T(N*N + N - 1) + 1.0*T(N*N + 2*N - 1);
	
	dcdk(N*N - 1) = dot(lambda, dKdk_uro);
	
	dcda(N*N - 1) = penal*(65.0-0.2)*pow(a(N*N - 1),penal-1)*dcdk(N*N - 1);
	
	//To check vanaf hier
	
	//Neumann links, boven dirichlett
	for (int i = 1; i < 0.3*N - 1; i++){
		vec dKdk_u((N+1)*(N+1));
		dKdk_u.fill(0.0);
		dKdk_u(i) = (-1.0 - 2.0*(pow(k(i-1,0),2)/pow((k(i-1,0) + k(i,0)),2)))*T(i) + 1.0*T(i+1) + (2.0*pow(k(i-1,0),2)/pow((k(i-1,0) + k(i,0)),2))*T(i + (N+1));
		dKdk_u(i+1) = (-1.0 - 2.0*(pow(k(i+1,0),2)/pow((k(i+1,0) + k(i,0)),2)))*T(i+1) + 1.0*T(i) + (2.0*pow(k(i+1,0),2)/pow((k(i+1,0) + k(i,0)),2))*T(i+1 + (N+1));
		dKdk_u(i + (N+1)) = (-2.0*(pow(k(i-1,0),2)/pow((k(i-1,0) + k(i,0)),2) + pow(k(i,0+1),2)/pow((k(i,0+1) + k(i,0)),2)))*T(i + (N+1)) + (2.0*(pow(k(i-1,0),2)/pow((k(i-1,0) + k(i,0)),2)))*T(i) + (2.0*pow(k(i,0+1),2)/pow((k(i,0+1) + k(i,0)),2))*T(i+1 + (N+1));
		dKdk_u(i+1 + (N+1)) = (-2.0*(pow(k(i+1,0),2)/pow((k(i+1,0) + k(i,0)),2) + pow(k(i,0+1),2)/pow((k(i,0+1) + k(i,0)),2)))*T(i+1 + (N+1)) + (2.0*(pow(k(i+1,0),2)/pow((k(i+1,0) + k(i,0)),2)))*T(i+1) + (2.0*(pow(k(i,0+1),2)/pow((k(i,0+1) + k(i,0)),2)))*T(i + (N+1));

		dcdk(i) = dot(lambda, dKdk_u);

		dcda(i) = penal*(65.0-0.2)*pow(a(i),penal-1)*dcdk(i);
	}
	
	//Neumann links, onder dirichlett
	for (int i = 0.7*N + 1; i < N - 1; i++){
		vec dKdk_u((N+1)*(N+1));
		dKdk_u.fill(0.0);
		dKdk_u(i) = (-1.0 - 2.0*(pow(k(i-1,0),2)/pow((k(i-1,0) + k(i,0)),2)))*T(i) + 1.0*T(i+1) + (2.0*pow(k(i-1,0),2)/pow((k(i-1,0) + k(i,0)),2))*T(i + (N+1));
		dKdk_u(i+1) = (-1.0 - 2.0*(pow(k(i+1,0),2)/pow((k(i+1,0) + k(i,0)),2)))*T(i+1) + 1.0*T(i) + (2.0*pow(k(i+1,0),2)/pow((k(i+1,0) + k(i,0)),2))*T(i+1 + (N+1));
		dKdk_u(i + (N+1)) = (-2.0*(pow(k(i-1,0),2)/pow((k(i-1,0) + k(i,0)),2) + pow(k(i,0+1),2)/pow((k(i,0+1) + k(i,0)),2)))*T(i + (N+1)) + (2.0*(pow(k(i-1,0),2)/pow((k(i-1,0) + k(i,0)),2)))*T(i) + (2.0*pow(k(i,0+1),2)/pow((k(i,0+1) + k(i,0)),2))*T(i+1 + (N+1));
		dKdk_u(i+1 + (N+1)) = (-2.0*(pow(k(i+1,0),2)/pow((k(i+1,0) + k(i,0)),2) + pow(k(i,0+1),2)/pow((k(i,0+1) + k(i,0)),2)))*T(i+1 + (N+1)) + (2.0*(pow(k(i+1,0),2)/pow((k(i+1,0) + k(i,0)),2)))*T(i+1) + (2.0*(pow(k(i,0+1),2)/pow((k(i,0+1) + k(i,0)),2)))*T(i + (N+1));

		dcdk(i) = dot(lambda, dKdk_u);

		dcda(i) = penal*(65.0-0.2)*pow(a(i),penal-1)*dcdk(i);
	}
	
	//Neumann rechts, boven dirichlett
	for (int i = 1; i < 0.3*N - 1; i++){
		vec dKdk_u((N+1)*(N+1));
		dKdk_u.fill(0.0);
		dKdk_u(i + (N-1)*(N+1)) = (-2.0*(pow(k(i,N-2),2)/pow((k(i,N-2) + k(i,N-1)),2) + pow(k(i-1,N-1),2)/pow((k(i-1,N-1) + k(i,N-1)),2)))*T(i + (N-1)*(N+1)) + (2.0*(pow(k(i,N-2),2)/pow((k(i,N-2) + k(i,N-1)),2)))*T(i+1 + (N-1)*(N+1)) + (2.0*pow(k(i-1,N-1),2)/pow((k(i-1,N-1) + k(i,N-1)),2))*T(i + (N)*(N+1));
		dKdk_u(i+1 + (N-1)*(N+1)) = (-2.0*(pow(k(i,N-2),2)/pow((k(i,N-2) + k(i,N-1)),2) + pow(k(i+1,N-1),2)/pow((k(i+1,N-1) + k(i,N-1)),2)))*T(i+1 + (N-1)*(N+1)) + (2.0*(pow(k(i,N-2),2)/pow((k(i,N-2) + k(i,N-1)),2)))*T(i + (N-1)*(N+1)) + (2.0*pow(k(i+1,N-1),2)/pow((k(i+1,N-1) + k(i,N-1)),2))*T(i+1 + (N)*(N+1));
		dKdk_u(i + (N)*(N+1)) = (-2.0*(pow(k(i-1,N-1),2)/pow((k(i-1,N-1) + k(i,N-1)),2)) - 1.0)*T(i + (N)*(N+1)) + (2.0*(pow(k(i-1,N-1),2)/pow((k(i-1,N-1) + k(i,N-1)),2)))*T(i + (N-1)*(N+1)) + 1.0*T(i+1 + (N)*(N+1));
		dKdk_u(i+1 + (N)*(N+1)) = (-2.0*(pow(k(i+1,N-1),2)/pow((k(i+1,N-1) + k(i,N-1)),2)) - 1.0)*T(i+1 + (N)*(N+1)) + (2.0*(pow(k(i+1,N-1),2)/pow((k(i+1,N-1) + k(i,N-1)),2)))*T(i+1 + (N-1)*(N+1)) + 1.0*T(i + (N)*(N+1));

		dcdk(i + (N-1)*N) = dot(lambda, dKdk_u);

		dcda(i + (N-1)*N) = penal*(65.0-0.2)*pow(a(i + N*(N-1)),penal-1)*dcdk(i + (N-1)*N);
	}
	
	//Neumann rechts, onder dirichlett
	for (int i = 0.7*N + 1; i < N - 1; i++){
		vec dKdk_u((N+1)*(N+1));
		dKdk_u.fill(0.0);
		dKdk_u(i + (N-1)*(N+1)) = (-2.0*(pow(k(i,N-2),2)/pow((k(i,N-2) + k(i,N-1)),2) + pow(k(i-1,N-1),2)/pow((k(i-1,N-1) + k(i,N-1)),2)))*T(i + (N-1)*(N+1)) + (2.0*(pow(k(i,N-2),2)/pow((k(i,N-2) + k(i,N-1)),2)))*T(i+1 + (N-1)*(N+1)) + (2.0*pow(k(i-1,N-1),2)/pow((k(i-1,N-1) + k(i,N-1)),2))*T(i + (N)*(N+1));
		dKdk_u(i+1 + (N-1)*(N+1)) = (-2.0*(pow(k(i,N-2),2)/pow((k(i,N-2) + k(i,N-1)),2) + pow(k(i+1,N-1),2)/pow((k(i+1,N-1) + k(i,N-1)),2)))*T(i+1 + (N-1)*(N+1)) + (2.0*(pow(k(i,N-2),2)/pow((k(i,N-2) + k(i,N-1)),2)))*T(i + (N-1)*(N+1)) + (2.0*pow(k(i+1,N-1),2)/pow((k(i+1,N-1) + k(i,N-1)),2))*T(i+1 + (N)*(N+1));
		dKdk_u(i + (N)*(N+1)) = (-2.0*(pow(k(i-1,N-1),2)/pow((k(i-1,N-1) + k(i,N-1)),2)) - 1.0)*T(i + (N)*(N+1)) + (2.0*(pow(k(i-1,N-1),2)/pow((k(i-1,N-1) + k(i,N-1)),2)))*T(i + (N-1)*(N+1)) + 1.0*T(i+1 + (N)*(N+1));
		dKdk_u(i+1 + (N)*(N+1)) = (-2.0*(pow(k(i+1,N-1),2)/pow((k(i+1,N-1) + k(i,N-1)),2)) - 1.0)*T(i+1 + (N)*(N+1)) + (2.0*(pow(k(i+1,N-1),2)/pow((k(i+1,N-1) + k(i,N-1)),2)))*T(i+1 + (N-1)*(N+1)) + 1.0*T(i + (N)*(N+1));

		dcdk(i + (N-1)*N) = dot(lambda, dKdk_u);

		dcda(i + (N-1)*N) = penal*(65.0-0.2)*pow(a(i + N*(N-1)),penal-1)*dcdk(i + (N-1)*N);
	}
	
	//links dirichlett
	for (int i = 0.3*N; i < 0.7*N; i++){
		vec dKdk_u((N+1)*(N+1));
		dKdk_u.fill(0.0);
		//dKdk_u(i) = (-1.0 - 2.0*(pow(k(i-1,0),2)/pow((k(i-1,0) + k(i,0)),2)))*T(i) + 1.0*T(i+1) + (2.0*pow(k(i-1,0),2)/pow((k(i-1,0) + k(i,0)),2))*T(i + (N+1));
		//dKdk_u(i+1) = (-1.0 - 2.0*(pow(k(i+1,0),2)/pow((k(i+1,0) + k(i,0)),2)))*T(i+1) + 1.0*T(i) + (2.0*pow(k(i+1,0),2)/pow((k(i+1,0) + k(i,0)),2))*T(i+1 + (N+1));
		dKdk_u(i + (N+1)) = (-2.0*(pow(k(i-1,0),2)/pow((k(i-1,0) + k(i,0)),2) + pow(k(i,0+1),2)/pow((k(i,0+1) + k(i,0)),2)))*T(i + (N+1)) + (2.0*(pow(k(i-1,0),2)/pow((k(i-1,0) + k(i,0)),2)))*T(i) + (2.0*pow(k(i,0+1),2)/pow((k(i,0+1) + k(i,0)),2))*T(i+1 + (N+1));
		dKdk_u(i+1 + (N+1)) = (-2.0*(pow(k(i+1,0),2)/pow((k(i+1,0) + k(i,0)),2) + pow(k(i,0+1),2)/pow((k(i,0+1) + k(i,0)),2)))*T(i+1 + (N+1)) + (2.0*(pow(k(i+1,0),2)/pow((k(i+1,0) + k(i,0)),2)))*T(i+1) + (2.0*(pow(k(i,0+1),2)/pow((k(i,0+1) + k(i,0)),2)))*T(i + (N+1));

		dcdk(i) = dot(lambda, dKdk_u);

		dcda(i) = penal*(65.0-0.2)*pow(a(i),penal-1)*dcdk(i);
	}
	
	//rechts dirichlett
	for (int i = 0.3*N; i < 0.7*N; i++){
		vec dKdk_u((N+1)*(N+1));
		dKdk_u.fill(0.0);
		dKdk_u(i + (N-1)*(N+1)) = (-2.0*(pow(k(i,N-2),2)/pow((k(i,N-2) + k(i,N-1)),2) + pow(k(i-1,N-1),2)/pow((k(i-1,N-1) + k(i,N-1)),2)))*T(i + (N-1)*(N+1)) + (2.0*(pow(k(i,N-2),2)/pow((k(i,N-2) + k(i,N-1)),2)))*T(i+1 + (N-1)*(N+1)) + (2.0*pow(k(i-1,N-1),2)/pow((k(i-1,N-1) + k(i,N-1)),2))*T(i + (N)*(N+1));
		dKdk_u(i+1 + (N-1)*(N+1)) = (-2.0*(pow(k(i,N-2),2)/pow((k(i,N-2) + k(i,N-1)),2) + pow(k(i+1,N-1),2)/pow((k(i+1,N-1) + k(i,N-1)),2)))*T(i+1 + (N-1)*(N+1)) + (2.0*(pow(k(i,N-2),2)/pow((k(i,N-2) + k(i,N-1)),2)))*T(i + (N-1)*(N+1)) + (2.0*pow(k(i+1,N-1),2)/pow((k(i+1,N-1) + k(i,N-1)),2))*T(i+1 + (N)*(N+1));
		//dKdk_u(i + (N)*(N+1)) = (-2.0*(pow(k(i-1,N-1),2)/pow((k(i-1,N-1) + k(i,N-1)),2)) - 1.0)*T(i + (N)*(N+1)) + (2.0*(pow(k(i-1,N-1),2)/pow((k(i-1,N-1) + k(i,N-1)),2)))*T(i + (N-1)*(N+1)) + 1.0*T(i+1 + (N)*(N+1));
		//dKdk_u(i+1 + (N)*(N+1)) = (-2.0*(pow(k(i+1,N-1),2)/pow((k(i+1,N-1) + k(i,N-1)),2)) - 1.0)*T(i+1 + (N)*(N+1)) + (2.0*(pow(k(i+1,N-1),2)/pow((k(i+1,N-1) + k(i,N-1)),2)))*T(i+1 + (N-1)*(N+1)) + 1.0*T(i + (N)*(N+1));

		dcdk(i + (N-1)*N) = dot(lambda, dKdk_u);

		dcda(i + (N-1)*N) = penal*(65.0-0.2)*pow(a(i + N*(N-1)),penal-1)*dcdk(i + (N-1)*N);
	}
	
	//Vanaf hier nog heel goed checken
	
	//Vakje met 1 hoekpunt van 293K, linkerkant, boven dirichlett boundary condition
	vec dKdk_hlb((N+1)*(N+1));
	dKdk_hlb.fill(0.0);
	dKdk_hlb(i) = (-1.0 - 2.0*(pow(k(i-1,0),2)/pow((k(i-1,0) + k(i,0)),2)))*T(i) + 1.0*T(i+1) + (2.0*pow(k(i-1,0),2)/pow((k(i-1,0) + k(i,0)),2))*T(i + (N+1));
	dKdk_hlb(i+1) = (-1.0 - 2.0*(pow(k(i+1,0),2)/pow((k(i+1,0) + k(i,0)),2)))*T(i+1) + 1.0*T(i) + (2.0*pow(k(i+1,0),2)/pow((k(i+1,0) + k(i,0)),2))*T(i+1 + (N+1));
	dKdk_hlb(i + (N+1)) = (-2.0*(pow(k(i-1,0),2)/pow((k(i-1,0) + k(i,0)),2) + pow(k(i,0+1),2)/pow((k(i,0+1) + k(i,0)),2)))*T(i + (N+1)) + (2.0*(pow(k(i-1,0),2)/pow((k(i-1,0) + k(i,0)),2)))*T(i) + (2.0*pow(k(i,0+1),2)/pow((k(i,0+1) + k(i,0)),2))*T(i+1 + (N+1));
	dKdk_hlb(i+1 + (N+1)) = (-2.0*(pow(k(i+1,0),2)/pow((k(i+1,0) + k(i,0)),2) + pow(k(i,0+1),2)/pow((k(i,0+1) + k(i,0)),2)))*T(i+1 + (N+1)) + (2.0*(pow(k(i+1,0),2)/pow((k(i+1,0) + k(i,0)),2)))*T(i+1) + (2.0*(pow(k(i,0+1),2)/pow((k(i,0+1) + k(i,0)),2)))*T(i + (N+1));

	dcdk(i) = dot(lambda, dKdk_hlb);

	dcda(i) = penal*(65.0-0.2)*pow(a(i),penal-1)*dcdk(i);
	
	//Vakje met 1 hoekpunt van 293K, linkerkant, onder dirichlett boundary condition
	vec dKdk_hlo((N+1)*(N+1));
	dKdk_hlo.fill(0.0);
	dKdk_hlo(i) = (-1.0 - 2.0*(pow(k(i-1,0),2)/pow((k(i-1,0) + k(i,0)),2)))*T(i) + 1.0*T(i+1) + (2.0*pow(k(i-1,0),2)/pow((k(i-1,0) + k(i,0)),2))*T(i + (N+1));
	dKdk_hlo(i+1) = (-1.0 - 2.0*(pow(k(i+1,0),2)/pow((k(i+1,0) + k(i,0)),2)))*T(i+1) + 1.0*T(i) + (2.0*pow(k(i+1,0),2)/pow((k(i+1,0) + k(i,0)),2))*T(i+1 + (N+1));
	dKdk_hlo(i + (N+1)) = (-2.0*(pow(k(i-1,0),2)/pow((k(i-1,0) + k(i,0)),2) + pow(k(i,0+1),2)/pow((k(i,0+1) + k(i,0)),2)))*T(i + (N+1)) + (2.0*(pow(k(i-1,0),2)/pow((k(i-1,0) + k(i,0)),2)))*T(i) + (2.0*pow(k(i,0+1),2)/pow((k(i,0+1) + k(i,0)),2))*T(i+1 + (N+1));
	dKdk_hlo(i+1 + (N+1)) = (-2.0*(pow(k(i+1,0),2)/pow((k(i+1,0) + k(i,0)),2) + pow(k(i,0+1),2)/pow((k(i,0+1) + k(i,0)),2)))*T(i+1 + (N+1)) + (2.0*(pow(k(i+1,0),2)/pow((k(i+1,0) + k(i,0)),2)))*T(i+1) + (2.0*(pow(k(i,0+1),2)/pow((k(i,0+1) + k(i,0)),2)))*T(i + (N+1));

	dcdk(i) = dot(lambda, dKdk_hlo);

	dcda(i) = penal*(65.0-0.2)*pow(a(i),penal-1)*dcdk(i);
	
	/*
	//Dirichlett hoekpunten van dcda_arit als basis
	
	//Vakje met 1 hoekpunt van 293K, linkerkant, boven dirichlett boundary condition
	vec dKdk_uhlb((N+1)*(N+1));
	dKdk_uhlb.fill(0.0);
	dKdk_uhlb(0.3*N - 1) = -1.0*T(0.3*N - 1) + 0.5*T(0.3*N) + 0.5*T(0.3*N - 1 + (N+1));
	//dKdk_uhlb(0.3*N) = -1.0*T(0.3*N) + 0.5*T(0.3*N - 1) + 0.5*T(0.3*N + (N+1));
	dKdk_uhlb(0.3*N - 1 + (N+1)) = -1.0*T(0.3*N - 1 + (N+1)) + 0.5*T(0.3*N - 1) + 0.5*T(0.3*N + (N+1));
	dKdk_uhlb(0.3*N + (N+1)) = -1.0*T(0.3*N + (N+1)) + 0.5*T(0.3*N) + 0.5*T(0.3*N - 1 + (N+1));

	dcdk(0.3*N - 1) = dot(lambda, dKdk_uhlb);

	dcda(0.3*N - 1) = penal*(65.0-0.2)*pow(a(0.3*N - 1),penal-1)*dcdk(0.3*N - 1);
	
	//Vakje met 1 hoekpunt van 293K, linkerkant, onder dirichlett boundary condition
	vec dKdk_uhlo((N+1)*(N+1));
	dKdk_uhlo.fill(0.0);
	//dKdk_uhlo(0.7*N) = -1.0*T(0.7*N) + 0.5*T(0.7*N + 1) + 0.5*T(0.7*N + (N+1));
	dKdk_uhlo(0.7*N + 1) = -1.0*T(0.7*N + 1) + 0.5*T(0.7*N) + 0.5*T(0.7*N + 1 + (N+1));
	dKdk_uhlo(0.7*N + (N+1)) = -1.0*T(0.7*N + (N+1)) + 0.5*T(0.7*N) + 0.5*T(0.7*N + 1 + (N+1));
	dKdk_uhlo(0.7*N + 1 + (N+1)) = -1.0*T(0.7*N + 1 + (N+1)) + 0.5*T(0.7*N + 1) + 0.5*T(0.7*N + (N+1));

	dcdk(0.7*N) = dot(lambda, dKdk_uhlo);

	dcda(0.7*N) = penal*(65.0-0.2)*pow(a(0.7*N),penal-1)*dcdk(0.7*N);
	
	//Vakje met 1 hoekpunt van 293K, rechterkant, boven dirichlett boundary condition
	vec dKdk_uhrb((N+1)*(N+1));
	dKdk_uhrb.fill(0.0);
	dKdk_uhrb(N*N-1 + 0.3*N - 1) = -1.0*T(N*N-1 + 0.3*N - 1) + 0.5*T(N*N-1 + 0.3*N) + 0.5*T(N*N+N + 0.3*N - 1);
	dKdk_uhrb(N*N-1 + 0.3*N) = -1.0*T(N*N-1 + 0.3*N) + 0.5*T(N*N-1 + 0.3*N - 1) + 0.5*T(N*N+N + 0.3*N);
	dKdk_uhrb(N*N+N + 0.3*N - 1) = -1.0*T(N*N+N + 0.3*N - 1) + 0.5*T(N*N-1 + 0.3*N - 1) + 0.5*T(N*N+N + 0.3*N);
	//dKdk_uhrb(N*N+N + 0.3*N) = -1.0*T(N*N+N + 0.3*N) + 0.5*T(N*N-1 + 0.3*N) + 0.5*T(N*N+N + 0.3*N - 1);

	dcdk(N*(N-1) + 0.3*N - 1) = dot(lambda, dKdk_uhrb);

	dcda(N*(N-1) + 0.3*N - 1) = penal*(65.0-0.2)*pow(a(N*(N-1) + 0.3*N - 1),penal-1)*dcdk(N*(N-1) + 0.3*N - 1);
	
	//Vakje met 1 hoekpunt van 293K, rechterkant, onder dirichlett boundary condition
	vec dKdk_uhro((N+1)*(N+1));
	dKdk_uhro.fill(0.0);
	dKdk_uhro(N*N-1 + 0.7*N) = -1.0*T(N*N-1 + 0.7*N) + 0.5*T(N*N-1 + 0.7*N + 1) + 0.5*T(N*N+N + 0.7*N);
	dKdk_uhro(N*N-1 + 0.7*N + 1) = -1.0*T(N*N-1 + 0.7*N + 1) + 0.5*T(N*N-1 + 0.7*N) + 0.5*T(N*N+N + 0.7*N + 1);
	//dKdk_uhro(N*N+N + 0.7*N) = -1.0*T(N*N+N + 0.7*N) + 0.5*T(N*N-1 + 0.7*N) + 0.5*T(N*N+N + 0.7*N + 1);
	dKdk_uhro(N*N+N + 0.7*N + 1) = -1.0*T(N*N+N + 0.7*N + 1) + 0.5*T(N*N-1 + 0.7*N + 1) + 0.5*T(N*N+N + 0.7*N);

	dcdk(N*(N-1) + 0.7*N) = dot(lambda, dKdk_uhro);

	dcda(N*(N-1) + 0.7*N) = penal*(65.0-0.2)*pow(a(N*(N-1) + 0.7*N),penal-1)*dcdk(N*(N-1) + 0.7*N);
	*/

	
	return dcda;
	
}
	
	
	


vec check(int N, double rmin, vec x, mat dc){

vec dcn = zeros<vec>(N*N);

for (int i = 0; i < N; i++) {
  for (int j = 0; j < N; j++) {
    double sum = 0.0; 
    for (int k = std::max(i - floor(rmin),0.0); k < std::min(i + floor(rmin),(double)N); k++) {
      for (int l = std::max(j - floor(rmin),0.0); l< std::min(j + floor(rmin),(double)N); l++) {
        double fac = rmin - sqrt(pow((i-k),2) + pow((j-l),2));
        sum = sum + std::max(0.0,fac);
        dcn(i*N + j) = dcn(i*N + j) + std::max(0.0,fac) * dc(N*k+l) * x(N*k+l);
      }
    }
    dcn(i*N + j) = dcn(i*N + j) / (x(i*N + j) * sum);
  }
}
return dcn;
}

vec OC(int N, vec x, double volfrac, vec dc){

double l1 = 0.0;
double l2 = 100000.0;
double move = 0.2;
vec xnew(N*N);

while ((l2-l1)>1e-4){
	double lmid = 0.5*(l2+l1);

	for (int i = 0; i < N*N; i++){
    		xnew(i) = max(0.001, max(x(i)-move,min(1.0,min(x(i)+move, x(i)*sqrt(-dc(i)/lmid)))));
	}	

	double sum_x = 0.0;

	for (int i = 0; i < N*N; i++){
		sum_x += xnew(i);
	}
 
   	if (sum_x - volfrac * (N*N) > 0){
   	     l1 = lmid;
   	}
   	else
   	{
   	    l2 = lmid;
   	}
}
return xnew;
}

}
#endif


