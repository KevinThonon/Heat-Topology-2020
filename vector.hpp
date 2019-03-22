#ifndef tws_vector_hpp
#define tws_vector_hpp
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
class vector {
public:
    typedef T   value_type ;
    typedef int size_type ;
    
public:
    vector( size_type size )
    : size_( size )
    , data_( new value_type[size_] )
    {}
    
    vector( size_type size, value_type val )
    : size_( size )
    , data_( new value_type[size_] )
    {
        std::fill_n(data_, size, val); 
    }

    ~vector()
    { delete [] data_ ; }
    
public: // Copy
    vector( vector const& that )
    : size_( that.size_ )
    , data_( new value_type[size_] )
    {
        //calls operator=
        (*this) = that ;
    }
    
    vector& operator=( vector const& that ) {
      assert( that.size() == size() ) ;
      std::copy( that.data_, that.data_+size_, data_ ) ;
      return *this ;
    }
    
    template <typename V>
    vector& operator=( V const& that ) {
      assert( that.size() == size() ) ;
        for (size_type i=0; i<size_; ++i) {
            data_[i] = that(i) ;
        }
        return *this ;
    }

public:// Access
    value_type operator() ( size_type i ) const {
        assert( i>=0 ) ;
        assert( i<size_ ) ;
        return data_[i] ;
    }
    
    value_type& operator() ( size_type i ) {
        assert( i>=0 ) ;
        assert( i<size_ ) ;
        return data_[i] ;
    }
    
    size_type size() const {
        return size_ ;
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
        return data_+size_;
    }
public: //arithmetic operations
     template <typename Vector>
      inline vector& operator-=(Vector const& v ) {
         assert(this->size()==v.size());
         for (size_type i=0; i<size(); ++i) { data_[i] = data_[i]-v(i) ; }
         return *this;
      }

      template <typename Vector>
      inline vector& operator+=(Vector const& v ) { 
         assert(this->size()==v.size());
         for (size_type i=0; i<size(); ++i) { data_[i] = data_[i]+v(i) ; }
         return *this;
      }

      void randomize(value_type min, value_type max, int seed=1){
         #ifdef NDEBUG
            seed = std::chrono::system_clock::now().time_since_epoch().count();
         #endif
         auto engine = std::default_random_engine(seed);
         std::uniform_real_distribution<value_type> distribution(min,max);
         for (size_type i=0; i<size(); ++i) { data_[i] =distribution(engine) ; }
      }      
    
public: // Fortran binding:
    typedef value_type* const& pointer ;
    pointer ptr() const {
        return data_ ;
    }
    
private:
    size_type   size_ ;
    value_type* data_ ;
};

template <class T>
struct is_vector : public std::false_type{};


template <class T>
struct is_vector<tws::vector<T> > : public std::true_type{};

template <typename T>
std::ostream& operator<<( std::ostream& os, vector<T> const& v ) {
    os << "[" << v.size() << "](" ;
    for (typename vector<T>::size_type i=0; i< v.size()-1; ++i) {
        os << v(i) << "," ;
    }
        os << v(v.size()-1)<<")" ;
    return os ;
}

template<typename V1, typename V2>
inline decltype(auto) inner_product(V1 const& v1, V2 const& v2){
    assert(v1.size()==v2.size());
    assert(v1.size()>0);
    decltype(v1(0)*v2(0)) sum=0; 
    for (typename V1::size_type i=0;i<v1.size();i++) sum+=v1(i)*v2(i);
    return sum;
}

template<typename V1>
inline decltype(auto) norm_2(V1 const& v1){
    return std::sqrt(inner_product(v1,v1));
}


}//namespace tws

#endif
