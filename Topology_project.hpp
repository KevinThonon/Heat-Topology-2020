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

double penal = 3.0;

// Dichtheid k in elk element met SIMP methode. Input: percentage metaal in elk element. Output: Dichtheid k in elk element
mat create_k(mat pctmetal, int N) {
	mat k(N,N);
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			k(i,j) = (65-0.2)*pow(pctmetal(i,j),penal) + 0.2;
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

/* // Matrix K aanmaken
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
} */

// Matrix K aanmaken zonder bigkmat
mat K_mat(mat k, int N){
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
}

// Kostfunctie = u^T * u / #elementen
double objective_function1(vec T, int N){
	double cost = dot(T,T)/(N*N);
	return cost;
}

// Kostfunctie = u / #elementen
// Deze kostfunctie is het gemiddelde van alle temperatuur-punten.
double objective_function2(vec T, int N){
	vec w(pow((N+1),2);
	w.fill(1.0);
	double cost = dot(w,T)/(pow((N+1),2));
	return cost;
}

// Kostfunctie = w^T * u / #elementen
// Deze kostfunctie is het gewogen gemiddelde van alle temperatuur-punten.
// w is het gewicht van elk temperatuur-punt. De som van alle gewichten is 1. Een middenpunt heeft gewicht 1/N² = h_x * h_y. Een randpunt heeft hiervan de helft (vakje is half zo klein).
// Een hoekpunt heeft gewicht h_x * h_y / 4 (vakje is een kwart van een middelste vakje).
double objective_function3(vec T, int N){
	double A = pow((1.0/N),2);
	vec w(pow((N+1),2);
	w.fill(A);
	w(0) = A/4.0;
	w(N) = A/4.0;
	w(pow(N,2)+N) = A/4.0;
	w(pow(N,2)+2*N) = A/4.0;
	for (int i = 1, i < N, i++){
		w(i) = A/2.0;
		w(pow(N,2)+2*N-i) = A/2.0;
		w(i*(N+1)) = A/2.0;
		w(i*(N+1)+N) = A/2.0;
	}
	double cost = dot(w,T);
	return cost;
}

// Lambda is de oplossing van K^T * lambda = -dg/du = -2*u / #elementen
vec lambda1(vec T, mat K, int N){
	vec lambda = solve(K.t(), (-2.0/(N*N))*T); //Veranderd van 2.0 naar -2.0
	return lambda;
}

// Lambda is de oplossing van K^T * lambda = -dg/du = - eenheidsvector / #elementen
vec lambda2(mat K, int N){
	vec w(pow((N+1),2);
	w.fill(-1.0);
	vec lambda = solve(K.t(), w);
	return lambda;
}

// Lambda is de oplossing van K^T * lambda = -dg/du = -w
vec lambda3(mat K, int N){
	double A = pow((1.0/N),2);
	vec w(pow((N+1),2);
	w.fill(A);
	w(0) = A/4.0;
	w(N) = A/4.0;
	w(pow(N,2)+N) = A/4.0;
	w(pow(N,2)+2*N) = A/4.0;
	for (int i = 1, i < N, i++){
		w(i) = A/2.0;
		w(pow(N,2)+2*N-i) = A/2.0;
		w(i*(N+1)) = A/2.0;
		w(i*(N+1)+N) = A/2.0;
	}
	vec lambda = solve(K.t(), -1.0*w);
	return lambda;
}



// dc/da is een matrix/vector van gradienten die nodig zijn in de optimalisatie stap
// Afhankelijk of matrix of vector nodig is in optimalisatie stap, moet mat of vec gecomment worden
//mat dcda(vec lambda, vec T, mat pctmetal, int N){
vec dcda(vec lambda, vec T, mat pctmetal, int N){
	//Initialiseren dc/da en opvullen met nullen
	//mat dcda(N,N);
	vec dcda(N*N);
	dcda.fill(0.0);
	//Initialiseren dc/dk en opvullen met nullen
	vec dcdk(N*N);
	dcdk.fill(0.0);
	
	//dc/da element per element opvullen
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			//Wat berekend moet worden is dc/da_i,j = dk/da * (lambda^T * dK/dk_i,j * u). Aangezien dK/dk_i,j slechts 4 niet-nul rijen zal bevatten zal niet eerst dK/dk_i,j aangemaakt worden om dan te vermenigvuldigen met u
			//maar zal dK/dk_i,j * u direct aangemaakt worden. 
			vec dKdk_u((N+1)*(N+1));
			dKdk_u.fill(0.0);
			dKdk_u(i + j*(N+1)) = -1.0*T(i + j*(N+1)) + 0.5*T(i+1 + j*(N+1)) + 0.5*T(i + (j+1)*(N+1));
			dKdk_u(i+1 + j*(N+1)) = -1.0*T(i+1 + j*(N+1)) + 0.5*T(i + j*(N+1)) + 0.5*T(i+1 + (j+1)*(N+1));
			dKdk_u(i + (j+1)*(N+1)) = -1.0*T(i + (j+1)*(N+1)) + 0.5*T(i + j*(N+1)) + 0.5*T(i+1 + (j+1)*(N+1));
			dKdk_u(i+1 + (j+1)*(N+1)) = -1.0*T(i+1 + (j+1)*(N+1)) + 0.5*T(i+1 + j*(N+1)) + 0.5*T(i + (j+1)*(N+1));
			//Vermenigvuldiging met lambda^T om tot dc/dk te komen
			//std::cout<<lambda.size()<<std::endl;
			//std::cout<<dKdk_u.size()<<std::endl;
			dcdk(i + j*N) = dot(lambda, dKdk_u);
			//Vermenigvuldiging met dk/da om tot dc/da te komen
			//dcda(i,j) = penal*(65.0-0.2)*pow(pctmetal(i,j),penal-1)*dcdk(i + j*N);
			dcda(i + j*(N)) = penal*(65.0-0.2)*pow(pctmetal(i,j),penal-1)*dcdk(i + j*(N));
		}
	}
	
	
	return dcda;
	
}
	
	
	

/* // Matrix dK/dk aanmaken
mat dKdk(int N){
	mat ll(((N+1)*(N+1)), ((N+1)*(N+1)));
	ll.fill(0.0);
	for (int i = 0.3*N; i < 0.7*N + 1; i++){
		ll(i,i) = 1.0;
		ll(((N+1)*(N+1))-i-1, ((N+1)*(N+1))-i-1) = 1.0;
	}
	
	//Linksboven hoekpunt
	ll(0,0) = -1.0;
	ll(0,1) = 0.5;
	ll(0,N+1) = 0.5;
	
	//Linksonder hoekpunt
	ll(N,0) = -1.0;
	ll(N,N-1) = 0.5;
	ll(N,2*N+1) = 0.5;
	
	//Rechtsboven hoekpunt
	ll(N*(N+1),N*(N+1)) = -1.0;
	ll(N*(N+1),N*(N+1)+1) = 0.5;
	ll(N*(N+1),N*N-1) = 0.5;
	
	//Rechtsonder hoekpunt
	ll(N*(N+2),N*(N+2)) = -1.0;
	ll(N*(N+2),N*(N+1)-1) = 0.5;
	ll(N*(N+2),N*(N+2)-1) = 0.5;
	
	//Neumann links
	for (int i = 1; i < 0.3*N; i++){
		ll(i,i) = -2.0;
		ll(i,i-1) = 0.5;
		ll(i,i+1) = 0.5;
		ll(i,i+N+1) = 1.0;
		
		ll(N-i,N-i) = -2.0;
		ll(N-i,N-i-1) = 0.5;
		ll(N-i,N-i+1) = 0.5;
		ll(N-i,N-i+N+1) = 1.0;
	}

	//Neumann rechts
	for (int i = 1; i < 0.3*N; i++){
		ll(N*(N+1)+i,N*(N+1)+i) = -2.0;
		ll(N*(N+1)+i,N*(N+1)+i-1) = 0.5;
		ll(N*(N+1)+i,N*(N+1)+i+1) = 0.5;
		ll(N*(N+1)+i,N*(N+1)+i-(N+1)) = 1.0;
		
		ll(N*(N+2)-i,N*(N+2)-i) = -2.0;
		ll(N*(N+2)-i,N*(N+2)-i-1) = 0.5;
		ll(N*(N+2)-i,N*(N+2)-i+1) = 0.5;
		ll(N*(N+2)-i,N*(N+2)-i-(N+1)) = 1.0;
	}
	
	//Neumann boven
	for (int i = 1; i < N; i++){
		ll(i*(N+1),i*(N+1)) = -2.0;
		ll(i*(N+1),i*(N+1)+1) = 1.0;
		ll(i*(N+1),(i-1)*(N+1)) = 0.5;
		ll(i*(N+1),(i+1)*(N+1)) = 0.5;
	}
	
	//Neumann onder
	for (int i = 1; i < N; i++){
		ll(i*(N+1)+N,i*(N+1)+N) = -2.0;
		ll(i*(N+1)+N,i*(N+1)+N-1) = 1.0;
		ll(i*(N+1)+N,(i-1)*(N+1)+N) = 0.5;
		ll(i*(N+1)+N,(i+1)*(N+1)+N) = 0.5;
	}
	
	//Gewone kolommen
	for (int j = 1; j < N; j++){
		for (int i = 1; i < N; i++){
			ll(j*(N+1)+i,j*(N+1)+i) = -4.0;
			ll(j*(N+1)+i,j*(N+1)+i-1) = 1.0;
			ll(j*(N+1)+i,j*(N+1)+i+1) = 1.0;
			ll(j*(N+1)+i,(j-1)*(N+1)+i) = 1.0;
			ll(j*(N+1)+i,(j+1)*(N+1)+i) = 1.0;
		}
	}
	
	return ll;
} */




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


