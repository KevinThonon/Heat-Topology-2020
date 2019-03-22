#ifndef tws_vector_expressions_hpp
#define tws_vector_expressions_hpp
#include"vector.hpp"
#include<type_traits>
namespace tws{
       
    template <typename V1, typename V2>
    class vector_sum {
         static_assert(is_vector<V1>::value,"vector_sum requires first argument to have succesfull is_vector");
         static_assert(is_vector<V2>::value,"vector_sum requires second argument to have succesfull is_vector");
        public:
        typedef typename V1::size_type size_type ;
        typedef decltype( typename V1::value_type() + typename V2::value_type() ) value_type ;

        public:
        vector_sum( V1 const& v1, V2 const& v2 )
        : v1_( v1 )
        , v2_( v2 )
        {
            assert( v1.size()==v2.size() ) ;
        }

        size_type size() const {
            return v1_.size() ;
        }

        value_type operator()( size_type i ) const {
            assert( i>=0 ) ;
            assert( i<size() ) ;
            return v1_(i) + v2_(i) ;
        }

        private:
        V1 const& v1_ ;
        V2 const& v2_ ;
    } ; //vector_sum

  template <class T,class T2>
  struct is_vector<vector_sum<T,T2> > {
    static bool const value = is_vector<T>::value && is_vector<T2>::value;
  };

    template <typename V1, typename V2>
    typename std::enable_if< is_vector<V2>::value && is_vector<V1>::value, vector_sum<V1,V2> >::type operator+(V1 const& v1, V2 const& v2 ) {
        return vector_sum<V1,V2>(v1,v2) ;
    } //operator+


    template <typename V1, typename V2>
    class vector_minus {
        public:
         static_assert(is_vector<V1>::value,"vector_minus requires first argument to have succesfull is_vector");
         static_assert(is_vector<V2>::value,"vector_minus requires second argument to have succesfull is_vector");
        typedef typename V1::size_type size_type ;
        typedef decltype( typename V1::value_type() + typename V2::value_type() ) value_type ;

        public:
        vector_minus( V1 const& v1, V2 const& v2 )
        : v1_( v1 )
        , v2_( v2 )
        {
            assert( v1.size()==v2.size() ) ;
        }

        size_type size() const {
            return v1_.size() ;
        }
    
        value_type operator()( size_type i ) const {
            assert( i>=0 ) ;
            assert( i<size() ) ;
            return v1_(i) - v2_(i) ;
        }

        private:
        V1 const& v1_ ;
        V2 const& v2_ ;
    } ; //vector_minus

  template <class T,class T2>
  struct is_vector<vector_minus<T,T2> > {
    static bool const value = is_vector<T>::value && is_vector<T2>::value;
  };

    template <typename V1, typename V2>
     typename std::enable_if< is_vector<V2>::value && is_vector<V1>::value, vector_minus<V1,V2> >::type operator-(V1 const& v1, V2 const& v2 ) {
        return vector_minus<V1,V2>(v1,v2) ;
    }//operator-

    template <typename S1, typename V2>
    class scalar_vector_multiply {
         static_assert(std::is_arithmetic<S1>::value,"scalar_vector_multiply requires first argument to have succesfull is_arithmetic");
         static_assert(is_vector<V2>::value,"scalar_vector_multiply requires second argument to have succesfull is_vector");
        public:
        typedef typename V2::size_type size_type ;
        typedef decltype( S1() * typename V2::value_type() ) value_type ;

        public:
        scalar_vector_multiply( S1 const& s1, V2 const& v2 )
        : s1_( s1 )
        , v2_( v2 )
        {
        }

        size_type size() const {
            return v2_.size() ;
        }

        value_type operator()( size_type i ) const {
            assert( i>=0 ) ;
            assert( i<size() ) ;
            return s1_*v2_(i) ;
        }

        private:
        S1 const& s1_ ;
        V2 const& v2_ ;
    } ; //scalar_vector_multiply

  template <class T,class T2>
  struct is_vector<scalar_vector_multiply<T,T2> > {
    static bool const value = is_vector<T2>::value && std::is_arithmetic<T>::value ;
  };

    template <typename S1, typename V2>
    typename std::enable_if< is_vector<V2>::value && std::is_arithmetic<S1>::value, scalar_vector_multiply<S1,V2> >::type operator*(S1 const& s1, V2 const& v2 ) {
        return scalar_vector_multiply<S1,V2>(s1,v2) ;
    }//operator*

    template <typename S1, typename V2>
    typename std::enable_if< is_vector<V2>::value && std::is_arithmetic<S1>::value , scalar_vector_multiply<S1,V2> >::type  operator*(V2 const& v2, S1 const& s1 ) {
        return scalar_vector_multiply<S1,V2>(s1,v2) ;
    }//operator*

} // namespace tws

#endif
