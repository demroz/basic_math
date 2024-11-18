#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "../header/linalg.hpp"
#include<cmath>

TEST_CASE("proj", "[proj]")
{
    
    Vector<std::complex<double>> a(3), b(3), c(3), d(3);
    a[0] = std::complex<double>{0,0};
    a[1] = std::complex<double>{1,1};
    a[2] = std::complex<double>{1,1};

    b[0] = std::complex<double>{1,1};
    b[1] = std::complex<double>{0,0};
    b[2] = std::complex<double>{1,1};
    
    c[0] = std::complex<double>{0.5,0.5};
    c[1] = std::complex<double>{0.0,0.0};
    c[2] = std::complex<double>{0.5,0.5};
    std::cout<<LA::proj(a,b)<<std::endl;
    d = LA::proj(a,b);
    REQUIRE(d==c);
    
}

TEST_CASE("Gram Schmidt", "[GS]")
{
    Matrix<double> S(2,2), SGS(2,2);
    S = {{3,2},{1,2}};
    Matrix<double> U(2,2);
    U = {{3.,-2./5.},{1.,6./5.}};
    SGS = LA::GramSchmidt(S);
    bool aeq;
    aeq = U.almost_equal(SGS,1e-10);
    REQUIRE(aeq);
    
}

TEST_CASE("HouseholderQR", "[QR]")
{
    Matrix<std::complex<double>> A(3,1);
    // A = {{12,-51,4},{6,167,-68},{-4,24,-41}};
    A = {{1.,std::complex<double>(1,std::sqrt(5))},
        {2.,1.},
        {1,-std::sqrt(5)}};
    
    std::cout<<A<<std::endl;
    LA::HouseHolderQR(A);
    Matrix<std::complex<double>> B(3,1);
    B = {{1.,0,0},
        {0,1.,0},
        {0,0,1.}};
    std::cout<< " this is b "<<B<<std::endl;
    std::vector<Matrix<std::complex<double>>> qrvec;
    qrvec = LA::HouseHolderQR(A);
    Matrix<std::complex<double>> Q(qrvec[0]), R(qrvec[1]);
    std::cout<<Q<<"\n";
    std::cout<<R<<"\n";
    std::cout<<Q*R<<"\n";
    REQUIRE( AlmostEqual(Q*R, A, 1e-10) );

    std::string filename;
    filename = "/home/demroz/Documents/code/basic_math/tests/A.txt";
    Matrix<std::complex<double>> W;
    W.read_csv(filename);
    Matrix<std::complex<double>> b(3,1);
    b(0,0) = std::complex<double>{1,0};
    b(1,0) = std::complex<double>{1,0};
    b(2,0) = std::complex<double>{1,0};
    std::cout<< W << " "<<b<<" "<< LA::solveUpperTriangular(W,b); 
    std::cout<< W * LA::solveUpperTriangular(W,b); 


}