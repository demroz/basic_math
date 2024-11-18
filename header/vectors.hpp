#ifndef __VECTORS_HPP__
#define __VECTORS_HPP__
#include <vector>
#include <cmath>
#include <iostream>
#include <complex>
#include <limits>
#include <assert.h>
#include "utils.hpp"
#include "rapidcsv.h"

template <class T>
class Vector
{
    private:
        unsigned int _dim;
        std::vector<T> _v;
    
    public:

        // constructors
        Vector();
        Vector(unsigned int dim);
        Vector(const Vector<T>& v);
        ~Vector();
        void read_csv(std::string filename);
        
        // int dim() const {return _v.size();};

        // utility operations
        Vector<T>& zero();
        Vector<T>& set(const T* v, unsigned int dim);
        Vector<T>& set(std::vector<T> v);
        unsigned int dim() {return _dim;};
        unsigned int dim() const {return _dim;};

        // math
        const double norm2();
        const double p_norm(unsigned int p);
        const T sum() const;
        const T dot(const Vector<T>& v) const;
        const Vector<T> cross(const Vector<T>& v) const;
        const Vector<T> abs() const;

        // operators

        // access 
        const T operator[] (const int i) const;
        T& operator[] (const int i);

        // equality
        bool operator == (const Vector<T>& v) const;
        bool operator != (const Vector<T>& v) const;

        // mathops
        Vector<T> operator* (const T& a) const;
        T operator* (const Vector<T>& v) const;
        Vector<T> operator+(const Vector<T> v) const;
        Vector<T> operator-(const Vector<T> v) const;
        Vector<T> operator-() const;

        Vector<T>& operator+=(const Vector<T>& v);
        Vector<T>& operator-=(const Vector<T>& v);
        
        Vector<T>& operator+=(const T s);
        Vector<T>& operator-=(const T s);
        Vector<T>& operator*=(const T s);
        Vector<T>& operator/=(const T s);

        Vector<T> operator < (const T s);            // compare each element with s, return vector of 1 or 0 based on test
        Vector<T> operator > (const T s);

        Vector<T> operator < (const Vector<T>& v);  // element-wise less than comparion, return vector of 1 or 0 based on test
        Vector<T> operator > (const Vector<T>& v);  // element-wise greater than comparion, return vector of 1 or 0 based on test

        Vector<T> Conj();
};
template <class T>
Vector<T>::~Vector()
{

}
template <class T>
inline void Vector<T>::read_csv(std::string filename)
{
    rapidcsv::Document doc(filename,rapidcsv::LabelParams(-1,-1));
    std::vector<std::string> col0 = doc.GetColumn<std::string>(0);
    std::vector<std::string> row0 = doc.GetRow<std::string>(0);

    int _rows, _cols, _numel;
    _rows = col0.size();
    _cols = row0.size();
    _numel = _rows*_cols;
    _dim = _numel;
    for(int i =0; i<_numel; i++)
    {
        _v.push_back(0);
    }
    std::cout<<(*this).dim()<<std::endl;
    std::cout<<" I AM HERE REEEEEEEEEEE "<<std::endl;
    std::cout<<"col " << col0.size()<<std::endl;
    std::cout<<row0.size()<<std::endl;
    int k = 0;
    for(int i = 0; i < _rows; i++)
    {
        std::cout<<"this is numrows " << _rows <<std::endl;
        std::vector<std::string> row = doc.GetRow<std::string>(i);
        for(int j = 0; j < _cols; j++)
        {
            std::cout<<"this is numcols " << _cols <<std::endl;
            std::string value = row[j];
            //std::string value = doc.GetCell<std::string>(i,j);
            T c;
            
            if ( value.find("(") != std::string::npos )
            {
                std::string s = convert_complex_string_to_streamstring(value);
                std::istringstream is(s);
                is >> c;

            }
            else
            {
                std::string s = value;
                std::istringstream is(s);
                is >> c;
            }
            std::cout<<k<< "  " << c<<std::endl;
            (*this)[k] = c;
            k++;
        }
    }
}
template <class T>
inline Vector<T>::Vector()
{
}
template <class T>
Vector<T>::Vector(unsigned int dim)
{
    _dim = dim;
    for (int i = 0; i < dim; i++ )
    {
        _v.push_back(0);
    }
}

template <class T>
Vector<T>::Vector(const Vector<T>& v)
{
    _dim = v.dim();
    _v.clear();
    for(int i = 0; i < _dim; i++)
    {
        _v.push_back(v[i]);
    }   
}

template <class T>
Vector<T>& Vector<T>::zero()
{
    for(int i = 0; i < _dim; i++ )
    {
        _v[i] = 0;
    }
    return *this;
}

template <class T>
Vector<T>& Vector<T>::set(const T* v, unsigned int dim)
{
    _dim = dim;
    _v.clear();
    for(int i = 0; i <_dim; i++)
    {
        _v.push_back(0);
    }
    for(int i = 0; i < _dim; i++)
    {
        _v[i] = v[i];
    }
    return *this;
}
template<class T>
Vector<T>& Vector<T>::set(std::vector<T> v)
{
    _dim = v.size();
    _v.clear();
    for(int i = 0; i <_dim; i++)
    {
        _v.push_back(0);
    }
    for(int i = 0; i < _dim; i++)
    {
        _v[i] = v[i];
    }
    return *this;
}

// math
template<class T>
inline const double Vector<T>::norm2()
{
    double ret = 0.;

    for(int i=0; i < _dim; i++)
    {
        ret += std::abs(_v[i])*std::abs(_v[i]);
    }
    return std::sqrt(ret);
}

template<class T>
inline const double Vector<T>::p_norm(unsigned int p) 
{
    double ret = 0.;

    for(int i=0; i < _dim; i++)
    {
        ret += std::pow(std::abs(_v[i]),p);
    }
    return std::pow(ret,1./p);

}

template<class T>
inline const T Vector<T>::sum() const
{
    T ret = 0.;

    for(int i=0; i < _dim; i++)
    {
        ret += _v[i];
    }
    return ret;
}

template<class T>
inline const T Vector<T>::dot(const Vector<T>& v) const
{
    T ret = 0.;

    for(int i=0; i < _dim; i++)
    {
        ret += _v[i]*v[i];
    }
    return ret;
}

template<class T>
inline const Vector<T> Vector<T>::cross(const Vector<T>& v) const
{
    assert(_dim == 3);

    Vector<T> ret(3);
    ret[0] = (_v[1] * v[2]) - (_v[2] * v[1]);
    ret[1] = (_v[2] * v[0]) - (_v[0] * v[2]);
    ret[2] = (_v[0] * v[1]) - (_v[1] * v[0]);
    return ret;

}

template<class T>
inline const Vector<T> Vector<T>::abs() const
{
    Vector<T> ret(_dim);
    for(int i =0; i < _dim; i++)
    {
        ret[i] = std::abs(_v[i]);
    }
    return ret;
}

// operator overload
template<class T>
const T Vector<T>::operator[] (const int i) const
{
    return _v[i];
}

template<class T>
T& Vector<T>::operator[] (const int i)
{
    return _v[i];
}
template<class T>
bool Vector<T>::operator == (const Vector<T>& v) const
{
    bool ret = true;
    for(int i = 0; i < _dim; i++)
    {
        ret = ret && (_v[i] == v[i]);
    }
    return ret;
}

template<class T>
bool Vector<T>::operator != (const Vector<T>& v) const
{
    bool ret = true;
    for(int i = 0; i < _dim; i++)
    {
        ret = ret && (_v[i] != v[i]);
    }
    return ret;
}

template<class T>
Vector<T> Vector<T>::operator*(const T& a) const
{
    Vector rtn(dim());
    for(int i = 0; i < dim(); i++)
    {
        rtn[i] = _v[i] * a;
    }
    return rtn;
}
template<class T>
T Vector<T>::operator* (const Vector<T>& v) const
{
    T rtn = 0;
    for (int i = 0; i < dim(); i++)
    {
        rtn += _v[i] * v[i];
    }
    
    return rtn;
}
template<class T>
Vector<T> Vector<T>::operator+(const Vector<T> v) const
{
    Vector<T> rtn(dim());
    for(int i =0; i < dim(); i++)
    {
        rtn[i] = _v[i] + v[i];
    }
    return rtn;
}
template<class T>
Vector<T> Vector<T>::operator-(const Vector<T> v) const
{
    Vector<T> rtn(dim());
    for(int i = 0; i < dim(); i++)
    {
        rtn[i] = _v[i] - v[i];
    }
    return rtn;
}
template<class T>
Vector<T> Vector<T>::operator-() const
{
    Vector<T> rtn(dim());

    for(int i=0; i<dim();i++)
    {
        rtn[i] = -_v[i];
    }
    return rtn;
}

template<class T>
Vector<T>& Vector<T>::operator+=(const Vector<T>& v)
{
    for(int i=0;i<dim();i++)
    {
        _v[i] += v[i];
    }
    return *this;
}
template<class T>
Vector<T>& Vector<T>::operator-=(const Vector<T>& v)
{
    for(int i=0;i<dim();i++)
    {
        _v[i] -= v[i];
    }
    return *this;
}
template<class T>
Vector<T>& Vector<T>::operator+=(const T s)
{
    for(int i=0;i<dim();i++)
    {
        _v[i] += s;
    }
    return *this;
}
template<class T>
Vector<T>& Vector<T>::operator-=(const T s)
{
    for(int i=0;i<dim();i++)
    {
        _v[i] -= s;
    }
    return *this;
}
template<class T>
Vector<T>& Vector<T>::operator*=(const T s)
{
    for(int i=0;i<dim();i++)
    {
        _v[i] *= s;
    }
    return *this;
}
template<class T>
Vector<T>& Vector<T>::operator/=(const T s)
{
    for(int i=0;i<dim();i++)
    {
        _v[i] /= s;
    }
    return *this;
}

template<class T>
inline Vector<T> Vector<T>::operator < (const T s)
{
    bool ret = true;
    for(int i= 0; i< dim(); i++)
    {
        ret &= _v[i] < s;
    }
    return ret;
}

template<class T>
inline Vector<T> Vector<T>::operator > (const T s)
{
    bool ret = true;
    for(int i= 0; i< dim(); i++)
    {
        ret &= _v[i] < s;
    }
    return ret;
}

template<class T>
inline Vector<T> Vector<T>::operator < (const Vector<T>& v)
{
    bool ret = true;
    for(int i= 0; i< dim(); i++)
    {
        ret &= _v[i] < v[i];
    }
    return ret;
}

template<class T>
inline Vector<T> Vector<T>::operator > (const Vector<T>& v)
{
    bool ret = true;
    for(int i= 0; i< dim(); i++)
    {
        ret &= _v[i] < v[i];
    }
    return ret;
}

template <class T>
inline Vector<T> Vector<T>::Conj()
{
    Vector<T> ret( dim() );
    for(int i = 0; i < dim(); i++)
    {
        ret[i] = conj((*this)[i]);
    }
    return ret;
}

template<class T>
inline Vector<T> operator*(const T& a, const Vector<T> v)
{
  return v*a;
}

/*
template<class T>
inline bool AlmostEqual(const Vector<T> v1, const Vector<T> v2, double TOL)
{
  Vector<T> tolerance_vector(v1.dim());
  tolerance_vector = v1-v2;
  tolerance_vector.abs();
  bool ret = true;
  for(int i =0; i < v1.dim(); i++)
  {
    ret &= tolerance_vector < TOL;
  }
  return ret;
}*/
template <class T>
inline std::ostream& operator<<(std::ostream& os, const Vector<T>& v)
{
  os << "[";
  for (int i=0; i<v.dim(); i++)
    {
      if (i!=0) os << ", ";
      os << v[i];
    }
  os << "]";
  return os;
}



#endif