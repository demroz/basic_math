#ifndef __LINEAR_REGRESSION__
#define __LINEAR_REGRESSION__
#include "matrix.hpp"
#include "vectors.hpp"
#include <complex>
#include "linalg.hpp"
#include "rapidcsv.h"

namespace LR
{
    template<typename T>
    std::vector<T> read_least_squares_data_from_csv(
        std::string filename, 
        std::vector<std::string> xVarNames,
        std::string yVarNames
    )
    {
        Matrix<T> X;
        Matrix<T> Y;
        
        rapidcsv::Document doc(filename);
        std::vector<T> ydata = doc.GetColumn<T>(yVarNames);
        std::vector<std::vector<T>> y;
        y.push_back(ydata);
        std::vector<std::vector<T>> x;
        std::vector<T> ones;
        
        for(int i = 0; i < ydata.size(); i++)
        {
            ones.push_back(1);
        }
        x.push_back(ones);
        
        for(int i = 0; i<xVarNames.size(); i++)
        {
            std::vector<T> col = doc.GetColumn<T>(xVarNames[i]);
            x.push_back(col);            
        }
        X.set(x);
        X = X.Transpose();
        
        Y.set(y);
        Y = Y.Transpose();
        std::vector<Matrix<T>> qr = LA::HouseHolderQR(X);
        Matrix<T> Q = qr[0];
        Matrix<T> R = qr[1];
        
        Matrix<T> sol;
        sol = LA::solveUpperTriangular(R, Q.H()*Y);
        std::vector<T> ret;
        for(int i = 0; i < sol.numRows(); i++)
        {
            ret.push_back(sol(i,0));
        }

        Matrix<T> F(Y.numRows(), 1);
        F = X*sol;
        Matrix<T> error(Y.numRows(),1);
        error = F - Y;
        std::vector<double> ss_res, ss_tol;
        T mean_y = 0;
        for(int i = 0; i < Y.numRows(); i++)
        {
            mean_y += Y(i,0);
        }
        mean_y /= Y.numRows();
        /* print solution*/
        std::cout<< "--------------- Regression Completed --------------- \n";

        return ret;
        
    }
    std::vector<double> ndLSQ(Matrix<double>& x, Matrix<double>&y)
    {

    }
    
} // namespace LR


#endif