#ifndef tws_matrix_hpp
#define tws_matrix_hpp
#include <cassert>
#include <iostream>
#include <cmath>
#include <typeinfo>
#include <algorithm>
#include <random>
#include <chrono>
#include <iterator>
#include "vector.hpp"
#include "vector_expressions.hpp"

namespace tws {

template <typename T>
class matrix {
public:
    typedef T   value_type ;
    typedef int size_type ;
    
public:
    matrix( size_type n, size_type m )
    : size_n( n )
    , size_m( m )
    , data_( new value_type[size_n * size_m] )
    {}
    
    matrix( size_type n, size_type m, value_type val )
    : size_n( n )
    , size_m( m )
    , data_( new value_type[size_n * size_m] ) 
    {
        std::fill_n(data_, n * m, val); 
    }

    ~matrix()
    { delete [] data_ ; }
    
public: // Copy
    matrix( matrix const& that )
    : size_n( that.size_n )
    , size_m( that.size_m )
    , data_( new value_type[size_n * size_m] ) 
    {
        //calls operator=
        (*this) = that ;
    }
    
    matrix& operator=( matrix const& that ) {
      assert( that.size() == size() ) ;
      std::copy( that.data_, that.data_+(size_n * size_m), data_ ) ;
      return *this ;
    }
    
    template <typename M>
    matrix& operator=( M const& that ) {
        for (size_type i=0; i<size_n; ++i) {
		        for (size_type j=0; j<size_m; ++j) {
            			data_[i * size_m + j] = that(i , j) ;
        		}
        }
        return *this ;
    }

public:// Access
    value_type operator() ( size_type i , size_type j) const {
        assert( i>=0 ) ;
        assert( i<size_n ) ;
        assert( j>=0 ) ;
        assert( j<size_m ) ;
        return data_[i * size_m + j] ;
    }
    
    value_type& operator() ( size_type i , size_type j) {

        return data_[i * size_m + j] ;
    }
    
    size_type size() const {
        return size_n*size_m ;
    }

    size_type sizen() const {
        return size_n ;
    }

    size_type sizem() const {
        return size_m ;
    }

    inline value_type* begin(){
        return data_;
    }
    inline value_type* end(){
        return data_+(size_n*size_m);
    }
    inline const value_type* cbegin() const{
        return data_;
    }
    inline const value_type* cend() const{
        return data_+(size_n*size_m);
    }

    
public: // Fortran binding:
    typedef value_type* const& pointer ;
    pointer ptr() const {
        return data_ ;
    }

public: // Functions
	int num_columns() {return size_m;}

	int num_rows() {return size_n;}
    
private:
    size_type   size_n, size_m ;
    value_type* data_ ;
};

template <class T>
struct is_matrix : public std::false_type{};


template <class T>
struct is_matrix<tws::matrix<T> > : public std::true_type{};


template <typename Matrix>
class transpose{
	public:
		typedef typename Matrix::size_type size_type;
		typedef typename Matrix::value_type value_type;

	public:
		transpose( Matrix const& m )
		:m_( m )
		{}

        	size_type sizem() const {
            	return m_.sizen() ;
        	}

        	size_type sizen() const {
            	return m_.sizem() ;
        	}

	public:

		value_type operator()(int i, int j) const {
			return m_(j,i);
		}

	private:
		tws::matrix<double> const& m_;
};

template <class Matrix>
struct is_matrix<transpose<Matrix> > {
    static bool const value = tws::is_matrix<Matrix>::value ;
};

template <typename Matrix>
typename std::enable_if< tws::is_matrix<Matrix>::value, transpose<Matrix> >::type trans(Matrix const& mat) {
    return transpose<Matrix>(mat) ;
}

template <typename Mat, typename Vec>
    class matvecmul {
         static_assert(tws::is_matrix<Mat>::value,"matvecmul requires first argument to have succesfull is_matrix");
         static_assert(tws::is_vector<Vec>::value,"matvecmul requires second argument to have succesfull is_vector");
        public:
        typedef typename Vec::size_type size_type ;
        typedef decltype( typename Mat::value_type() * typename Vec::value_type() ) value_type ;

        public:
        matvecmul( Mat const& mat, Vec const& vec )
        : mat_( mat )
        , vec_( vec )
        {}

        size_type size() const {
            return mat_.sizen() ;
        }

        value_type operator()( size_type i ) const {
        assert( i>=0 ) ;
        assert( i<size() ) ;
	value_type temp = 0;
	for (size_type j=0; j<vec_.size(); ++j) {
		temp = temp + mat_(i,j)*vec_(j);
	}
        return temp ;
        }

        private:
        Mat const& mat_ ;
        Vec const& vec_ ;
    } ; //matvecmul

  template <class T,class T2>
  struct is_vector<matvecmul<T,T2> > {
    static bool const value = tws::is_vector<T2>::value && tws::is_matrix<T>::value ;
  };

    template <typename Mat, typename Vec>
    typename std::enable_if< tws::is_vector<Vec>::value && tws::is_matrix<Mat>::value, matvecmul<Mat,Vec> >::type operator*(Mat const& mat, Vec const& vec ) {
        return matvecmul<Mat,Vec>(mat,vec) ;
    }//operator*



// Functor SRDA1 is the implementation of task 3.

struct SRDA1 {

	SRDA1(tws::matrix<double> const& X, double const& beta)
	: X_(X), beta_(beta)
	{} 

	void operator()(tws::vector<double> const& x, tws::vector<double> & b) const{
		b = trans(X_)*(X_*x) + beta_*x;
	}

	tws::matrix<double> X_;
	double beta_;
};

// Functor SRDA2 is the implementation of task 4.

struct SRDA2 {

	SRDA2(tws::matrix<double> const& X, double const& beta)
	: X_(X), beta_(beta)
	{} 

	void operator()(tws::vector<double> const& x, tws::vector<double> & b) const{
		tws::vector<double> t(X_.sizen());
		t = X_*x;
		b = trans(X_)*t + beta_*x;
	}

	tws::matrix<double> X_;
	double beta_;
};

// Functor SRDA is the improved implementation of the SRDA functor of the 4th homework. Now the vector x is not copied anymore but a reference is used.

class SRDA {
public:
	typedef tws::vector<double>  vec ;
	typedef tws::matrix<double>  mat ;
    
public:
	SRDA(mat const& X , double const& beta)
	: X_(X), beta_(beta)
	{}

	void operator()(vec const& x, vec & b) const {
	b = transpose(X_)*(X_*x) + beta_*x;
	}

private:
	tws::matrix<double> X_;
	double beta_;

};

// Functor laplace is the improved implementation of the laplace functor of the 4th homework. Now the matrix A is not created anymore.
struct laplace {

	laplace(int N)
	: N_(N)
	{} 

	tws::vector<double> operator()(vector<double> const& x, vector<double> & y) const{
	y(0) = 2 * x(0) - x(1);
	y(N_-1) = 2 * x(N_-1) - x(N_-2);
	for (int j=1; j<N_-1; ++j) {
		y(j) = 2 * x(j) - x(j-1) - x(j+1);
	}

	return y;
	}

	int N_;
};

// Functor toeplitz is the improved implementation of the toeplitz functor of the 4th homework. Now the matrix B is not created anymore.
struct toeplitz {

	toeplitz(int N)
	: N_(N)
	{} 

	tws::vector<double> operator()(vector<double> const& x, vector<double> & y) const{
	tws::vector<double> v(2 * N_ - 1);
	for (int i=0; i<v.size(); ++i) {
		v(i) = rand() % 20 + 1 ;  // random integers in the range [1,20]
	}
	for (int r=0; r<N_; ++r){
		double temp = 0 ;
		for (int s=0; s<N_; ++s){
			temp = temp + v(N_-r-1+s)*x(s);
		}
		y(r) = temp;	
	}

	return y;
	}

	int N_;
};
  
// Functor hankel is the improved implementation of the hankel functor of the 4th homework. Now the matrix C is not created anymore.
struct hankel {

	hankel(int N)
	: N_(N)
	{} 

	tws::vector<double> operator()(vector<double> const& x, vector<double> & y) const{
	tws::vector<double> v(2 * N_ - 1);
	for (int i=0; i<v.size(); ++i) {
		v(i) = rand() % 20 + 1 ;  // random integers in the range [1,20]
	}
	for (int r=0; r<N_; ++r){
		double temp = 0 ;
		for (int s=0; s<N_; ++s){
			temp = temp + v(r+s)*x(s);
		}
		y(r) = temp;	
	}

	return y;
	}

	int N_;
};

}
#endif
