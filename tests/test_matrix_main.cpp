#include<iostream>
#include "../header/matrix.hpp"

int main()
{
    Matrix<double> M(2,2);
    M = {5,2,3,4};
    std::cout<<M<<std::endl;
    M = {{1,2,3},{4,5,6},{7,8,9}};
    std::cout<<M<<std::endl;
    Matrix<double> A(2,2,{1.,2.,3.,4.});
    std::cout<<A<<std::endl;
    Matrix<double> B({{1,2,3},{4,5,6},{7,8,9}});
    std::cout<<B<<std::endl;
}