#ifndef tws_matrixexpr_hpp
#define tws_matrixexpr_hpp
#include <cassert>
#include <iostream>
#include <cmath>
#include <typeinfo>
#include <algorithm>
#include <random>
#include <chrono>
#include <iterator>
#include "vector.hpp"
#include "matrix.hpp"
#include "tws_util.hpp"
#include "vector_expressions.hpp"
#include "matrix_expressions.hpp"

namespace tws {
    template <typename M>
	class transpose{
	static_assert(is_matr<M>::value,"transpose requires matrix");
	public:
	        typedef typename M::size_type size_type ;
			typedef typename M::value_type value_type ;
	public: 
		transpose(M const& matr): 
		matr_(matr){
			
		}
		
		
		value_type operator()( size_type i, size_type j ) const {
            assert( i>=0 ) ;
			assert( j>=0 ) ;
            assert( i<matr_.num_columns() ) ;
			assert( j<matr_.num_rows() ) ;

            return matr_(j,i);
        }
	size_type num_rows() const{
		return matr_.num_columns();
	}
	size_type num_columns() const{
		return matr_.num_rows();
	}
	private:
	M const& matr_;
	};
	
	template <typename M, typename V>
    class matrix_vector_multiply {
         static_assert(is_matr<M>::value);
         static_assert(is_vector<V>::value);
        public:
        typedef typename V::size_type size_type ;
        typedef decltype( typename M::value_type() * typename V::value_type() ) value_type ;

        public:
        matrix_vector_multiply( M const& matrix, V const& vec )
        : matrix_( matrix )
        , vector_( vec )
        {
			assert(vec.size() == matrix.num_columns());
        }

        size_type size() const {
            return matrix_.num_rows() ;
        }

        value_type operator()( size_type i ) const {
            assert( i>=0 ) ;
            assert( i<matrix_.num_rows()) ;
			value_type sum = 0;
			for (size_type j=0; j<vector_.size(); j++){sum = sum + vector_(j)*matrix_(i,j);};
            return sum ;
        }

        private:
        M const& matrix_ ;
        V const& vector_ ;
    } ; //matrix_vector_multiply

	



  template <class M>
  struct is_matr<transpose<M>> {
    static bool const value = is_matr<M>::value;
  };


  template <class M,class V>
  struct is_vector<matrix_vector_multiply<M,V>> {
    static bool const value = is_matr<M>::value && is_vector<V>::value;
  };


template <typename M, typename V>
typename std::enable_if< is_matr<M>::value && is_vector<V>::value , matrix_vector_multiply<M,V> >::type  operator*(M const& matrix, V const& vector ) {
        return matrix_vector_multiply<M,V>(matrix,vector) ;
    }//operator*
}//namespace tws


#endif