#ifndef tws_matrix_hpp
#define tws_matrix_hpp
#include "vector.hpp"
#include "vector_expressions.hpp"
#include <cassert>
#include <iostream>
#include <cmath>
#include <typeinfo>
#include <algorithm>
#include <random>
#include <chrono>
#include <iterator>

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
    , data_( new value_type[size_n*size_m] )
    {}
    
    matrix( size_type n, size_type m, value_type val )
    : size_n( n )
    , size_m( m )
    , data_( new value_type[size_n*size_m] )
    {
	std::fill_n( data_, n*m, val ); 

    }

    ~matrix()
    { delete [] data_ ; }
    
public: // Copy
    matrix( matrix const& that )
    : size_n( that.size_n )
    , size_m( that.size_m )
    , data_( new value_type[size_n*size_m] )
    {
        //calls operator=
        (*this) = that ;
    }
    
    matrix& operator=( matrix const& that ) {
      assert( that.size() == size() ) ;
      std::copy( that.data_, that.data_+(size_n*size_m), data_ ) ;
      return *this ;
    }
    
    template <typename M>
    matrix& operator=( M const& that ) {
        for (size_type i=0; i<size_n; ++i) {
		for (size_type j=0; j<size_m; ++j) {
            		data_[i*size_m+j] = that(i,j) ;
		}
        }
        return *this ;
    }


public:// Access
    value_type operator() ( size_type i, size_type j ) const {
        //assert( i>=0 ) ;
        //assert( i<size_n ) ;
        //assert( j>=0 ) ;
        //assert( j<size_m ) ;
        return data_[i*size_m+j] ;
    }
    
    value_type& operator() ( size_type i, size_type j  ) {
        //assert( i>=0 ) ;
        //assert( i<size_n ) ;
        //assert( j>=0 ) ;
        //assert( j<size_m ) ;
        return data_[i*size_m+j] ;
    }
    
    size_type size() const {
        return size_n*size_m ;
    }
    size_type size_row() const {
        return size_n ;
    }
    size_type size_column() const {
        return size_m ;
    }
    inline value_type* begin(){
        return data_;
    }
    inline value_type* end(){
        return data_+size_;
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

public: // Number of colums and rows
    int num_columns() {return size_m ;}
    int num_rows() {return size_n ;}


private:
    size_type   size_n, size_m, size_ ;
    value_type* data_;
};

template <class T>
struct is_matrix : public std::false_type{};

template <class T>
struct is_matrix<tws::matrix<T> > : public std::true_type{};


// task 1: transpose of matrix using an expression template.
template <typename Matrix>
class transpose_matrix {
	static_assert(is_matrix<Matrix>::value,"transpose requires matrix");
public:
	typedef typename Matrix::size_type size_type ;
	typedef typename Matrix::value_type value_type ;

public:
	transpose_matrix( Matrix const& m1)
	: m_(m1)
	{}

	double operator() (int i, int j) const
	{return m_(j,i);}

        size_type size_row() const {
            return m_.size_column() ;
        }
        size_type size_column() const {
            return m_.size_row() ;
        }

private:
	Matrix const& m_;

};

template <class Matrix>
struct is_matrix<transpose_matrix<Matrix> > {
	static bool const value = tws::is_matrix<Matrix>::value ;
};

template <typename Matrix>
typename std::enable_if< tws::is_matrix<Matrix>::value, transpose_matrix<Matrix> >::type transpose(Matrix const& m1){
return transpose_matrix<Matrix>(m1);
}

// task 2: matrix vector multiply
template <typename M1, typename V2>
class matrix_vector_multiply {
	static_assert(tws::is_matrix<M1>::value,"matrix_vector_multiply requires first argument to have succesfull is_matrix");
        static_assert(tws::is_vector<V2>::value,"matrix_vector_multiply requires second argument to have succesfull is_vector");
        public:
        typedef typename V2::size_type size_type ;
        typedef decltype(typename M1::value_type() * typename V2::value_type() ) value_type ;

        public:
        matrix_vector_multiply( M1 const& m1, V2 const& v2 )
        : m1_( m1 )
        , v2_( v2 )
        {
        }

        size_type size() const {
            return m1_.size_row() ;
        }

        value_type operator()( size_type i ) const {
            assert( i>=0 ) ;
            //assert( i<m1_.size_row() ) ;
            value_type sum = 0;
	    for (size_type j = 0; j<v2_.size(); j++){
            	sum = sum + m1_(i,j)*v2_(j);
            }
	    return sum;
	}


        private:
        M1 const& m1_ ;
        V2 const& v2_ ;
    	} ; //matrix_vector_multiply


template <class T,class T2>
struct is_vector<matrix_vector_multiply<T,T2> > {
	static bool const value = tws::is_vector<T2>::value && tws::is_matrix<T>::value ;
};


template <typename M1, typename V2>
typename std::enable_if< tws::is_vector<V2>::value && tws::is_matrix<M1>::value, matrix_vector_multiply<M1,V2> >::type operator*(M1 const& m1, V2 const& v2 ) {
        return matrix_vector_multiply<M1,V2>(m1,v2) ;
}//operator*


// task 3: matrix-vector product of y = transpose(X_)*(X_*x) + beta_*x
struct SRDA {
  SRDA(tws::matrix<double> const& X, double const& beta)
  : X_(X), beta_(beta)
  {}

  void operator()(tws::vector<double> const& x, tws::vector<double> &y) const {
  	y = transpose(X_)*(X_*x) + beta_*x;
  }

  tws::matrix<double> X_;
  double beta_;

  };


// task 4
struct SRDA1 {
  SRDA1(tws::matrix<double> const& X , double const& beta)
  : X_(X), beta_(beta)
  {}

  void operator()(tws::vector<double> const& x, tws::vector<double> &y) const {
	tws::vector<double> t(X_.size_row());
	t = X_*x;
  	y = transpose(X_)*t + beta_*x;
  }

  tws::matrix<double> X_;
  double beta_;

  };

// Improved code of assignment 4 (these functors were first located in functors.hpp but transferred to matrix.hpp now)
// There's no matrix created anymore 
// The vectors x and y are not copied
struct laplacian {

	laplacian(int n)
	: n_(n)
	{} 

	tws::vector<double> operator()(vector<double> const& x, vector<double> &y) const{
		assert(x.size() == n_);

		y(0) = 2*x(0) - x(1);
		y(n_-1) = - x(n_-2) + 2*x(n_-1);
		for (int i=1; i<n_-1; ++i) {
			y(i) = - x(i-1) + 2*x(i) - x(i+1);
		}
		return y;
	}

	int n_;
};

// Improved code of assignment 4 (these functors were first located in functors.hpp but transferred to matrix.hpp now)
// There's no matrix created anymore 
// The vectors x and y are not copied
struct toeplitz {

	toeplitz(int n)
	: n_(n)
	{} 

	tws::vector<double> operator()(vector<double> const& x, vector<double> &y) const{
		assert(x.size() == n_);
		assert(y.size() == n_);
		tws::vector<double> vec(2*n_-1);

		for (int i=0; i<vec.size(); ++i) {
			vec(i) = rand() % 30 + 1;
		}

		for (int i=0; i<n_; i++) {
			double sum = 0.0;
			for (int j=0; j<n_; j++) {
				sum = sum + vec(n_-1+j-i)*x(j);	
			}
			y(i) = sum;
		}
	return y;
	}

	int n_;
};

// Improved code of assignment 4 (these functors were first located in functors.hpp but transferred to matrix.hpp now)
// There's no matrix created anymore 
// The vectors x and y are not copied
struct hankel {
    

	hankel(int n)
	: n_(n)
	{} 

	tws::vector<double> operator()(vector<double> const& x, vector<double> &y) const{
		assert(x.size() == n_);
		assert(y.size() == n_);
		tws::vector<double> vec(2*n_-1);
		

		for (int i=0; i<vec.size(); ++i) {
			vec(i) = rand() % 30 + 1;
		}

		for (int i=0; i<n_; i++) {
			double sum = 0.0;
			for (int j=0; j<n_; j++) {
				sum = sum + vec(i+j)*x(j);	
			}
			y(i) = sum;
		}
	return y;
	}

	int n_;
};



}

#endif
