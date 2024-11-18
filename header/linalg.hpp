#ifndef __LINALG_HPP__
#define __LINALG_HPP__
#include "matrix.hpp"
#include "vectors.hpp"
#include <complex>

#include <type_traits>
#include <typeinfo>
#ifndef _MSC_VER
#   include <cxxabi.h>
#endif
#include <memory>
#include <string>
#include <cstdlib>
namespace LA
{
    template<class T>
    Vector<T> proj(Vector<T>& v, Vector<T>& u)
    {
        Vector<T> ret(u.dim());
        ret = v.Conj()*u / (u.Conj()*u) * u;
        return ret;
    }

    template<class T>
    Matrix<T> GramSchmidt(const Matrix<T>& V)
    {
        Vector<T> uk(V.numRows());
        Vector<T> vk(V.numRows());

        Matrix<T> U(V.numRows(), V.numCols());

        // v_0
        vk = V.getCol(0);
        // U(0) = v_0
        U.setCol(0, vk);

        for(int i = 1; i<V.numCols(); i++)
        {
            vk = V.getCol(i);
            uk = vk;
            for(int j = 0; j < i; j++)
            {
                Vector<T> uj(V.numCols());
                uj = U.getCol(j);
                uk = uk - proj(vk, uj);
            }
            U.setCol(i,uk);
        }
        return U;
    }

    template<class T>
    std::vector<Matrix<T>> HouseHolderQR(const Matrix<T>& A)
    {
  
        int m = A.numRows();
        int n = A.numCols();
        Matrix<T> R(m,n);
        R = A;
        
        Matrix<T> Q(m,m);
        Q = eye<T>(m,m);
        T TWO = 2;
        int maxiter = 0;
        if ( m > n)
        {
            maxiter = n;
        }
        else
        {
            maxiter = n - 1;
        }
        for(int k = 0; k < maxiter ; k++)
        {
            Matrix<T> x = R.getSubBlock(k,m,k,k+1);
            Matrix<T> vk(m-k,1);
            vk = x;
            vk(0,0) = computeHouseHolderAlphaPhase(x(0,0)) * std::sqrt( (x.H() * x)(0,0) ) + x(0,0);

            vk = vk / std::sqrt( (vk.H() * vk)(0,0) );

            Matrix<T> Qk(m,m);
            Qk = eye<T>(m,m);
            Qk.setSubBlock(k,m,k,m,eye<T>(m-k,m-k) - TWO * (vk * vk.H()));
            Q = Qk*Q;
            R.setSubBlock(k,m,k,n, 
                        R.getSubBlock(k,m,k,n) 
                        - TWO * ( vk * (vk.H() * R.getSubBlock(k,m,k,n)) ) );
    
        }
        std::vector<Matrix<T>> ret;
        ret.push_back(Q.H());
        ret.push_back(R);
        return ret;
    }

    template<class T>
    Matrix<T> solveUpperTriangular(Matrix<T> R, Matrix<T> b)
    {
        int m,n,p,q;
        m = R.numRows();
        n = R.numCols();
        p = b.numRows();
        q = b.numCols();
        std::cout<<m<<" "<<n<<" "<<p<<" "<<q<<std::endl;
        if(m != p)
        {
            std::cout<< "matrix shapes must match" << std::endl;
        }
        Matrix<T> rtn(n,1);
        rtn.zero();
        for(int i = n-1; i >= 0; i--)
        {
            T s = 0;
            for(int j = i+1; j < n; j++) s = s+R(i,j)*rtn(j,0);
            rtn(i,0) = (b(i,0) - s)/R(i,i);
        } 
        return rtn;

    }

}
#endif