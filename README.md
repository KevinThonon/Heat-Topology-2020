Year: 2019-2020  
Project: heat  
Author 1: Kevin Thonon  
Author 2: Philippe Dewart  
Author 3: Urbaan Lemmens  
Needed libraries:
* Armadillo, http://arma.sourceforge.net, C++ 

Instructions:
* To compile the code:
* g++ -Wall -std=c++17 -o test Topology_OC.cpp -larmadillo
* To run the code and test you need to give some command line arguments: 
* The size of the grid (int), Initial guess for % metal for every cell (double), penalty (double)
* An example: 50, 0.3 and 3.0
* ./test 50 0.3 3.0
Can we use and modify your code for demonstrations concerning the course? Yes
