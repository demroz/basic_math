#ifndef __MATRIX_HPP__
#define __MATRIX_HPP__
#include "vectors.hpp"
#include <initializer_list>
#include "utils.hpp"
#include "rapidcsv.h"
#include <sstream>
#include<complex>
#include<regex>
#include <algorithm>

template <class T>
class Matrix
{
    private:
        unsigned int _rows;
        unsigned int _cols;
        unsigned int _numel;

        T *_data;
        unsigned int ravel(unsigned int i, unsigned int j);
    
    public:
        Matrix();
        Matrix(unsigned int rows, unsigned int cols);

        // copy constructor
        Matrix(std::vector<std::vector<T>> m);
        Matrix(const Matrix<T>& m);
        // assignment
        Matrix(unsigned int rows, unsigned int cols, const std::initializer_list<T>& list);
        Matrix(const std::initializer_list<std::initializer_list<T>>& list);
        ~Matrix();
        void read_csv(std::string filename);

        // get dimensions
        unsigned int numRows() const;
        unsigned int numCols() const;

        // access op
        const T& operator()(const int& i, const int& j) const;
        T& operator()(const int& i,const int& j);

        // utility operations
        Matrix<T>& zero(); // zero matrix
        Matrix<T>& set(std::vector<std::vector<T>> m); // set from std::vector<std::vector>>

        // assignment operator
        Matrix<T>& operator=(const Matrix& matrix);
        Matrix<T>& operator=(const std::initializer_list<T>& list);
        Matrix<T>& operator=(const std::initializer_list<std::initializer_list<T>>& list);
        
        // bool operators
        bool operator == (const Matrix<T>& RHS) const;
        bool operator != (const Matrix<T>& RHS) const;
        bool almost_equal(const Matrix<T>& RHS, double TOL) const;
        

        /*
         * math
        */

       // matrix - scalar
        Matrix<T> operator* (const T& a) const;
        Matrix<T> operator+ (const T& a) const;
        Matrix<T> operator- (const T& a) const;
        Matrix<T> operator/ (const T& a) const;

        Matrix<T> operator* (const T& a);
        Matrix<T> operator+ (const T& a);
        Matrix<T> operator- (const T& a);
        Matrix<T> operator/ (const T& a);

        Matrix<T>& operator+= (const T& a);
        Matrix<T>& operator-= (const T& a);
        Matrix<T>& operator*= (const T& a);
        Matrix<T>& operator/= (const T& a);

        // matrix - vector
        Vector<T> operator* (const Vector<T>& v) const;

        // matrix/matrix
        Matrix<T> operator*(const Matrix<T>& RHS) const;
        Matrix<T> operator+(const Matrix<T>& RHS) const;
        Matrix<T> operator-(const Matrix<T>& RHS) const;

        
        Matrix<T> Conj();
        Matrix<T> Transpose();
        Matrix<T> H();

        // elementary ops
        Vector<T> getCol(unsigned int col);
        Vector<T> getRow(unsigned int row);
        Vector<T> getCol(unsigned int col) const ;
        Vector<T> getRow(unsigned int row) const ;

        Matrix<T> getSubBlock(unsigned int ibegin, unsigned int iend, unsigned int jbegin, unsigned int jend);
        Matrix<T> getSubBlock(unsigned int ibegin, unsigned int iend, unsigned int jbegin, unsigned int jend) const;

        void setSubBlock(unsigned int ibegin, unsigned int iend,
                        unsigned int jbegin, unsigned int jend,
                        const Matrix<T>& block);


        //
        void setCol(unsigned int col, const Vector<T>& v);
        void setRow(unsigned int row, const Vector<T>& v);

        void swapRows(unsigned int i, unsigned int j);
        void swapCols(unsigned int i, unsigned int j);
};

// base constructor
template <class T>
inline Matrix<T>::Matrix()
{
    _rows = 0;
    _cols = 0;
    _numel = 0;
}

// constructor with memory allocation
template <class T>
inline Matrix<T>::Matrix(unsigned int rows, unsigned int cols)
{
    _rows = rows;
    _cols = cols;
    _numel = rows*cols;
    _data = new T[_numel];
    for(int i = 0; i < _numel; i++)
    {
        _data[i] = 0;
    }

}

template <class T>
inline unsigned int Matrix<T>::ravel(unsigned int i, unsigned int j)
{
    return i*(numCols()) + j;
}

// set matrix data equal to vales in std::vector
template <class T>
inline Matrix<T>::Matrix(std::vector<std::vector<T>> m)
{
    _rows = m.size();
    _cols = m[0].size();
    _numel = _rows*_cols;

    int k = 0;
    for(int i =0; i < _rows; i++)
    {
        for(int j=0; j< _cols; j++)
        {
            _data[k] = m[i][j];
            k+=1;
        }
    }
}

template <class T>
inline Matrix<T>::Matrix(const Matrix<T> &m)
{
    _rows = m.numRows();
    _cols = m.numCols();
    _numel = _rows*_cols;
    _data = new T[_numel];
    
    int k = 0;
    for(int i = 0; i < numRows(); i++)
    {
        for(int j = 0; j < numCols(); j++)
        {
            _data[k] = m(i,j);
            k++;
        }
    }
    
}

template <class T>
inline Matrix<T>::Matrix(unsigned int rows, unsigned int cols, const std::initializer_list<T> &list) : Matrix<T>()
{   
    _rows = rows;
    _cols = cols;
    _numel = _rows*_cols;
    _data = new T[_numel];
    int i = 0;
    for(T value : list)
    {
        _data[i] = value;
        i++;
    }
 
}

template <class T>
inline Matrix<T>::Matrix(const std::initializer_list<std::initializer_list<T>> &list)
{
    _rows = list.size();
    _cols = list.begin()[0].size();
    _numel = _rows*_cols;
    _data = new T[_numel];
    
    int k = 0;
    for( std::initializer_list<T> row : list)
    {
        for( T value : row)
        {
            _data[k] = value;
            k++;
        }
    }

}

template <class T>
inline Matrix<T>::~Matrix()
{
    if (_numel != 0)
    {
        delete[] _data;
    }
}


template <class T>
inline void Matrix<T>::read_csv(std::string filename)
{
    rapidcsv::Document doc(filename,rapidcsv::LabelParams(-1,-1));
    std::vector<std::string> col0 = doc.GetColumn<std::string>(0);
    std::vector<std::string> row0 = doc.GetRow<std::string>(0);

    _rows = col0.size();
    _cols = row0.size();
    _numel = _rows*_cols;
    _data = new T[_numel];
    std::cout<<"col " << col0.size()<<std::endl;
    std::cout<<row0.size()<<std::endl;
    int k = 0;
    for(int i = 0; i < numRows(); i++)
    {
        std::cout<<"this is numrows " << numRows() <<std::endl;
        std::vector<std::string> row = doc.GetRow<std::string>(i);
        for(int j = 0; j < numCols(); j++)
        {
            std::cout<<"this is numcols " << numCols() <<std::endl;
            std::string value = row[j];
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
            std::cout<<(*this).numRows() << " " << (*this).numCols() << i<<","<<j<<" i am ceeee "<<c<<std::endl;
            (*this)(i,j) = c;
        }
    }
}

template<class T>
unsigned int Matrix<T>::numRows() const
{
    return _rows;
}
template<class T>
unsigned int Matrix<T>::numCols() const
{
    return _cols;
}

template <class T>
inline const T &Matrix<T>::operator()(const int &i, const int &j) const
{
    unsigned int idx;
    idx = i*(numCols()) + j;
    return _data[idx];
    
}

template <class T>
inline T &Matrix<T>::operator()(const int &i, const int &j)
{
    unsigned int idx;
    idx = i*(numCols()) + j;
    return _data[idx];
}


template <class T>
inline Matrix<T> &Matrix<T>::operator=(const Matrix &matrix)
{
    this->~Matrix();
    _rows = matrix.numRows();
    _cols = matrix.numCols();
    _numel = _rows*_cols;

    _data = new T[_rows*_cols];
    
    for(int i = 0; i < _numel; i++)
    {
        _data[i] = matrix._data[i];
    }
    return *this;
}

template <class T>
inline Matrix<T> &Matrix<T>::operator=(const std::initializer_list<T> &list)
{
    int i = 0;
    for(T value : list)
    {
        _data[i] = value;
        i++;
    }
}

template <class T>
inline Matrix<T>& Matrix<T>::operator=(const std::initializer_list<std::initializer_list<T>>& list)
{
    _rows = list.size();
    _cols = list.begin()[0].size();
    if ( _numel != 0 )
    {
        delete[] _data;
    }
    _numel = _rows*_cols;
    _data = new T[_numel];

    int k = 0;
    for (std::initializer_list<T> row : list)
    {
        for (T value : row)
        {
            _data[k] = value;
            k++;
        }
    }
}
template <class T>
inline Matrix<T> &Matrix<T>::zero()
{
    int k = 0;
    for(int i=0; i < _numel; i++)
    {
        _data[k] = 0;
        k++;
    }
}

template <class T>
inline Matrix<T> &Matrix<T>::set(std::vector<std::vector<T>> m)
{
    _rows = m.size();
    _cols = m[0].size();
    if ( _numel != 0 )
    {
        delete[] _data;
    }
    _numel = _rows*_cols;
    _data = new T[_numel];

    int k = 0;
    for(int i =0; i < _rows; i++)
    {
        for(int j=0; j< _cols; j++)
        {
            _data[k] = m[i][j];
            k+=1;
        }
    }
}

template <class T>
inline bool Matrix<T>::operator==(const Matrix<T> &RHS) const
{
    if ( (_rows != RHS.numRows() ) || (_cols != RHS.numCols() ) )
    {
        return false;
    }
    int k = 0;
    for(int i = 0; i < _rows; i++)
    {
        for(int j = 0; j < _cols; j++)
        {
            if (_data[k] != RHS(i,j))
            {
                return false;
            }
            k++;
        }
    }
    return true;
}

template <class T>
inline bool Matrix<T>::operator!=(const Matrix<T> &RHS) const
{
    if ( (_rows != RHS.numRows() ) || (_cols != RHS.numCols() ) )
    {
        return true;
    }
    
    if ( (*this) == RHS )
    {
        return false;
    }
    return false;
}

template <class T>
inline bool Matrix<T>::almost_equal(const Matrix<T> &RHS, double TOL) const
{
    for(int i = 0; i < _rows; i++)
    {
        for(int j = 0; j < _cols; j++)
        {
            if(std::abs((*this)(i,j) - RHS(i,j)) > TOL)
            {
                return false;
            }
        }
    }
    return true;
}

template <class T>
inline Matrix<T> Matrix<T>::operator*(const T &a) const
{
    Matrix<T> rtn(_rows,_cols);
    for(int i=0; i < _rows; i++)
    {
        for(int j=0;j<_cols;j++)
        {
            rtn(i,j) = (*this)(i,j) * a;
        }
    }
    return rtn;
}

template <class T>
inline Matrix<T> Matrix<T>::operator*(const T &a)
{
    Matrix<T> rtn(_rows,_cols);
    for(int i=0; i < _rows; i++)
    {
        for(int j=0;j<_cols;j++)
        {
            rtn(i,j) = (*this)(i,j) * a;
        }
    }
    return rtn;
}

template <class T>
inline std::ostream &operator<<(std::ostream &os, const Matrix<T> &m_RHS)
{
  os << "[";
  for (int i=0; i<m_RHS.numRows(); i++)
    {
      for(int j=0; j<m_RHS.numCols(); j++)
	{
	  if(j!=0) os << ", ";
	  os << m_RHS(i,j);
	}
      if(i!=m_RHS.numRows()-1) os << ";\n ";
    }
  os << "]";
  return os;	
  

}

template <class T>
inline Matrix<T> eye(int rows, int cols)
{

  Matrix<T> rtn(rows,cols);
  rtn.zero();
  unsigned int idx = std::min(rows, cols);
  T one = 1;
  for(int i = 0; i<idx;i++)
    {
        rtn(i,i) = one;
    }
  return rtn;
}


template <class T>
inline Matrix<T> Matrix<T>::operator+(const T &a) const
{
    Matrix<T> rtn(_rows,_cols);
    for(int i=0; i < _rows; i++)
    {
        for(int j=0;j<_cols;j++)
        {
            rtn(i,j) = (*this)(i,j) + a;
        }
    }
    return rtn;
}

template <class T>
inline Matrix<T> Matrix<T>::operator+(const T &a)
{
    Matrix<T> rtn(_rows,_cols);
    for(int i=0; i < _rows; i++)
    {
        for(int j=0;j<_cols;j++)
        {
            rtn(i,j) = (*this)(i,j) + a;
        }
    }
    return rtn;
}
template <class T>
inline Matrix<T> Matrix<T>::operator-(const T &a) const
{
    Matrix<T> rtn(_rows,_cols);
    for(int i=0; i < _rows; i++)
    {
        for(int j=0;j<_cols;j++)
        {
            rtn(i,j) = (*this)(i,j) - a;
        }
    }
    return rtn;
}

template <class T>
inline Matrix<T> Matrix<T>::operator-(const T &a)
{
    Matrix<T> rtn(_rows,_cols);
    for(int i=0; i < _rows; i++)
    {
        for(int j=0;j<_cols;j++)
        {
            rtn(i,j) = (*this)(i,j) - a;
        }
    }
    return rtn;
}

template <class T>
inline Matrix<T> Matrix<T>::operator/(const T &a) const
{
    Matrix<T> rtn(_rows,_cols);
    for(int i=0; i < _rows; i++)
    {
        for(int j=0;j<_cols;j++)
        {
            rtn(i,j) = (*this)(i,j) / a;
        }
    }
    return rtn;
}

template <class T>
inline Matrix<T> Matrix<T>::operator/(const T &a)
{
    Matrix<T> rtn(_rows,_cols);
    for(int i=0; i < _rows; i++)
    {
        for(int j=0;j<_cols;j++)
        {
            rtn(i,j) = (*this)(i,j) / a;
        }
    }
    return rtn;
}


template <class T>
inline Matrix<T>& Matrix<T>::operator+=(const T &a)
{
    int k = 0;
    for(int i=0; i < _rows; i++)
    {
        for(int j=0;j<_cols;j++)
        {
            _data[k] += a;
            k+=1;
        }
    }
    
    return *this;
}
template <class T>
inline Matrix<T>& Matrix<T>::operator*=(const T &a)
{
    int k = 0;
    for(int i=0; i < _rows; i++)
    {
        for(int j=0;j<_cols;j++)
        {
            _data[k] = (*this)(i,j) * a;
            k+=1;
        }
    }
    return *this;
}
template <class T>
inline Matrix<T>& Matrix<T>::operator-=(const T &a)
{
    int k = 0;
    for(int i=0; i < _rows; i++)
    {
        for(int j=0;j<_cols;j++)
        {
           _data[k] = (*this)(i,j) - a;
           k++;
        }
    }
    return *this;
}
template <class T>
inline Matrix<T>& Matrix<T>::operator/=(const T &a)
{
    int k = 0;
    for(int i=0; i < _rows; i++)
    {
        for(int j=0;j<_cols;j++)
        {
            _data[k] = (*this)(i,j) / a;
            k++;
        }
    }
    return *this;
}
template <class T>
inline Vector<T> Matrix<T>::operator*(const Vector<T> &v) const
{
    Vector<T> rtn(_rows);
    for(int i = 0; i < _rows; i++)
    {
        for(int j=0; j<_cols; j++)
        {
            rtn[i] += (*this)(i,j)* v[j];
        }
    }
    return rtn;
}

template <class T>
inline Matrix<T> Matrix<T>::operator*(const Matrix<T> &RHS) const
{
    int n, m, p;
    n = _rows;
    m = _cols;
    p = RHS.numCols();

    int idx;
    Matrix<T> C(n,p);
    C.zero();
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < p; j++)
        {
            T sum = 0;
            for(int k = 0; k < m; k++)
            {
                sum += (*this)(i,k) * RHS(k,j);
            }
            C(i,j) = sum;
        }
    }
    return C;
}

template <class T>
inline Matrix<T> Matrix<T>::operator+(const Matrix<T> &RHS) const
{
    Matrix rtn(_rows, _cols);
    for(int i =0; i < _rows; i++)
    {
        for(int j =0; j < _cols; j++)
        {
            rtn(i,j) = (*this)(i,j) + RHS(i,j);
        }
    }
    return rtn;
}

template <class T>
inline Matrix<T> Matrix<T>::operator-(const Matrix<T> &RHS) const
{
    Matrix rtn(_rows, _cols);
    for(int i =0; i < _rows; i++)
    {
        for(int j =0; j < _cols; j++)
        {
            rtn(i,j) = (*this)(i,j) - RHS(i,j);
        }
    }
    return rtn;
}


template <class T>
inline Matrix<T> operator*(const T &a, const Matrix<T> m)
{
    return m*a;
}

template<class T>
inline Matrix<T> Matrix<T>::Transpose()
{
    Matrix<T> ret(_cols, _rows);
    for(int i = 0; i < _rows; i++)
    {
        for (int j = 0; j < _cols; j++ )
        {
            ret(j,i) = (*this)(i,j);
        }
    }
    return ret;
}


template<class T>
inline Matrix<T> Matrix<T>::Conj()
{
    Matrix<T> ret(_rows, _cols);
    for(int i = 0; i < _rows; i++)
    {
        for (int j = 0; j < _cols; j++ )
        {
            ret(i,j) = conj((*this)(i,j));
        }
    }
    return ret;
}


template <class T>
inline Matrix<T> Matrix<T>::H()
{
    Matrix<T> ret(_cols, _rows);
    for(int i = 0; i < _rows; i++)
    {
        for (int j = 0; j < _cols; j++ )
        {
            ret(j,i) = conj((*this)(i,j));
        }
    }
    return ret;
}

template <class T>
inline Vector<T> Matrix<T>::getCol(unsigned int col)
{
    Vector<T> rtn(_rows);
    for (int i = 0; i < _rows; i++)
    {
        rtn[i] = (*this)(i,col);
    }
    return rtn;
}

template <class T>
inline Vector<T> Matrix<T>::getRow(unsigned int row)
{
    Vector<T> rtn(_cols);
    for (int i = 0; i < _cols; i++)
    {
        rtn[i] = (*this)(row,i);
    }
    return rtn;
}

template <class T>
inline Vector<T> Matrix<T>::getCol(unsigned int col) const 
{
    Vector<T> rtn(_rows);
    for (int i = 0; i < _rows; i++)
    {
        rtn[i] = (*this)(i,col);
    }
    return rtn;
}

template <class T>
inline Vector<T> Matrix<T>::getRow(unsigned int row) const
{
    Vector<T> rtn(_cols);
    for (int i = 0; i < _cols; i++)
    {
        rtn[i] = (*this)(row,i);
    }
    return rtn;
}

template <class T>
inline Matrix<T> Matrix<T>::getSubBlock(unsigned int ibegin, unsigned int iend, unsigned int jbegin, unsigned int jend)
{
    Matrix<T> rtn(iend-ibegin, jend-jbegin);
    for(int i = 0; i< iend - ibegin; i++)
    {
        for(int j = 0; j < jend- jbegin; j++)
        {
            rtn(i,j) = (*this)(ibegin + i,jbegin + j);
        }
    }
    return rtn;
}

template <class T>
inline Matrix<T> Matrix<T>::getSubBlock(unsigned int ibegin, unsigned int iend, unsigned int jbegin, unsigned int jend) const
{
    Matrix<T> rtn(iend-ibegin, jend-jbegin);
    for(int i = 0; i< iend - ibegin; i++)
    {
        for(int j = 0; j < jend - jbegin; j++)
        {
            rtn(i,j) = (*this)(ibegin + i,jbegin + j);
        }
    }
    return rtn;
}

template <class T>
inline void Matrix<T>::setSubBlock(unsigned int ibegin, unsigned int iend, unsigned int jbegin, unsigned int jend, const Matrix<T> &block)
{
    for(int i = 0; i< iend - ibegin; i++)
    {
        for(int j = 0; j < jend - jbegin; j++)
        {
            (*this)(ibegin + i,jbegin + j) = block(i,j);
        }
    }
}

template <class T>
inline void Matrix<T>::setCol(unsigned int col, const Vector<T>& v)
{
    for (int i = 0; i < _rows; i++)
    {
        (*this)(i, col) = v[i];
    }
    
}

template <class T>
inline void Matrix<T>::setRow(unsigned int row, const Vector<T>& v)
{
    for (int i = 0; i < _cols; i++)
    {
        (*this)(row, i) = v[i];
    }
    
}
template <class T>
inline void Matrix<T>::swapRows(unsigned int i, unsigned int j)
{
    std::cout<<_rows << " " <<_cols<<std::endl;
    for(int k = 0; k < _cols; k++)
    {
        T value_i, value_j;
        value_i = (*this)(i, k);
        value_j = (*this)(j, k);
        (*this)(i, k) = value_j;
        (*this)(j, k) = value_i;
    }

}

template <class T>
inline void Matrix<T>::swapCols(unsigned int i, unsigned int j)
{
    for(int k = 0; k < _rows; k++)
    {
        T value_i, value_j;
        value_i = (*this)(k,i);
        value_j = (*this)(k,j);
        (*this)(k,i) = value_j;
        (*this)(k,j) = value_i;
    }
}

template <class T>
inline bool AlmostEqual(const Matrix<T>& A, const Matrix<T>& B, double TOL)
{
    int m,n,p,q;
    m = A.numRows();
    n = A.numCols();
    p = B.numRows();
    q = B.numCols();

    if( (m != p) || (n != q))
    {
        return false;
    }

    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j < n; j++)
        {
            double norm;
            norm = std::abs(A(i,j) - B(i,j));
            if (norm > TOL)
            {
                return false;
            }
        }
    }
    return true;
}

#endif