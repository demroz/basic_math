#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "../header/matrix.hpp"
#include "../header/vectors.hpp"
#include<cmath>

Matrix<double> m(5,5);
std::vector<std::vector<double>> vals = { {1,2,3,4,5},
                      {6,7,8,9,10},
                      {11,12,13,14,15},
                      {16,17,18,19,20},
                      {21,22,23,24,25}};
TEST_CASE( "Constructor with array given values", "[constructor]" ) 
{
    m.set(vals);

    for(int i = 0; i < 5; i++)
    {
        for(int j = 0; j < 5; j++)
        {
            REQUIRE( m(i,j) == vals[i][j]);
        }
    }

    m.zero();
    
    for(int i = 0; i < 5; i++)
    {
        for(int j = 0; j < 5; j++)
        {
            REQUIRE( m(i,j) == 0);
        }
    }
    
    
    Matrix<double> Q(5,11);
    REQUIRE(Q.numRows() == 5);
    REQUIRE(Q.numCols() == 11);

    Matrix<double> z;
    z.set(vals);


}


TEST_CASE( "Identity", "[eye]" ) 
{
    Matrix<double> I(3,5);
    I = eye<double>(3,5);
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 5; j++)
        {
            if (i == j)
            {
                REQUIRE(I(i,j) == 1);
            }
            else
            {
                REQUIRE(I(i,j) == 0);
            }
        }
    }

}

TEST_CASE( "Matrix/scalar", "[scalar operators]" ) 
{
    Matrix<double> m(5,5);
    std::vector<std::vector<double>> vals = { {1,2,3,4,5},
                      {6,7,8,9,10},
                      {11,12,13,14,15},
                      {16,17,18,19,20},
                      {21,22,23,24,25}};
    m.set(vals);

    Matrix<double> b(5,5);

    // *
    b = m*5;
    for(int i =0; i < 5; i++)
    {
        for(int j=0; j < 5; j++)
        {
            REQUIRE(b(i,j) == m(i,j) * 5);
        }
    }
    
    // +
    b = m+5;
    for(int i =0; i < 5; i++)
    {
        for(int j=0; j < 5; j++)
        {
            REQUIRE(b(i,j) == m(i,j) + 5);
        }
    }
    
    // -
    b = m-5;
    for(int i =0; i < 5; i++)
    {
        for(int j=0; j < 5; j++)
        {
            REQUIRE(b(i,j) == m(i,j) - 5);
        }
    }
    
    
    // /
    b = m/5;
    for(int i =0; i < 5; i++)
    {
        for(int j=0; j < 5; j++)
        {
            REQUIRE(b(i,j) == m(i,j) / 5);
        }
    }
    
    b = m - 1;
    b += 1;
    // += 
    
    b = m - 1;
    b += 1;
    for(int i =0; i < 5; i++)
    {
        for(int j=0; j < 5; j++)
        {
            REQUIRE(b(i,j) == m(i,j) );
        }
    }
    
    // -= 
    b = m + 1;
    b -= 1;
    for(int i =0; i < 5; i++)
    {
        for(int j=0; j < 5; j++)
        {
            REQUIRE(b(i,j) == m(i,j) );
        }
    }
    
    // *=
    b = m;
    b *= 2;
    for(int i =0; i < 5; i++)
    {
        for(int j=0; j < 5; j++)
        {
            REQUIRE(b(i,j) == m(i,j) * 2 );
        }
    }

    // /=
    std::cout<<b<<std::endl;
    b = m;
    std::cout<<b<<std::endl;
    b /= 2;
    std::cout<<b<<std::endl;
    std::cout<<m<<std::endl;
    for(int i =0; i < 5; i++)
    {
        for(int j=0; j < 5; j++)
        {
            REQUIRE(b(i,j) == m(i,j) / 2 );
        }
    }
    
}

TEST_CASE( "Matrix/vector", "[mmv]" ) 
{
    Matrix<double> m(5,5);
    std::vector<std::vector<double>> vals = { 
                      {1,2,3,4,5},
                      {6,7,8,9,10},
                      {11,12,13,14,15},
                      {16,17,18,19,20},
                      {21,22,23,24,25}};
    m.set(vals);

    Vector<double> v(5);
    std::vector<double> vector_vals = {1,0,5,8,15};

    v.set(vector_vals);

    Vector<double> res(5);
    std::vector<double> res_vals = {123,268,413,558,703};
    res.set(res_vals);
    REQUIRE( res == m*v);
    
    

}

TEST_CASE( "Matrix/matrix", "[dmm]" ) 
{  
    Matrix<double> m(5,5);
    std::vector<std::vector<double>> vals = { 
                      {1,2,3,4,5},
                      {6,7,8,9,10},
                      {11,12,13,14,15},
                      {16,17,18,19,20},
                      {21,22,23,24,25}};
    m.set(vals);
    Matrix<double> o(5,5);
    o = m + 1;

    std::vector<std::vector<double>> result = {
        {230,245,260,275,290},
        {530,570,610,650,690},
        {830,895,960,1025,1090},
        {1130,1220,1310,1400,1490},
        {1430,1545,1660,1775,1890}
    };
    std::cout<<m<<std::endl;
    std::cout<<o<<std::endl;
    Matrix<double> res(5,5);
    res.set(result);
    REQUIRE(res == m*o);
    REQUIRE( o - 1 == m);


}
TEST_CASE( "transpose", "[transpose]" )
{
    
    Matrix<double> A(2,3);
    std::cout<<"A"<<std::endl;
    std::cout<<A.numRows() << " " << A.numCols() <<std::endl;

    std::vector<std::vector<double>> vals = { 
                      {1,2,3},
                      {6,7,8}};
    A.set(vals);
    std::cout<<"A"<<std::endl;
    std::cout<<A.numRows() << " " << A.numCols() <<std::endl;
    std::cout<<A<<std::endl;
    Matrix<double> B = A.Transpose();
    std::cout<<"B"<<std::endl;
    std::cout<<B.numRows() << " " << B.numCols() <<std::endl;
    std::cout<<B<<std::endl;
    // std::cout<<"A"<<std::endl;
    // std::cout<<A<<std::endl;
    // std::cout<<"B = A.T"<<std::endl;
    // std::cout<<B<<std::endl;
    
    
    REQUIRE(B.Transpose() == A);
    REQUIRE(A.Conj() == A);

    Matrix<std::complex<double>> C(2,3);
    std::vector<std::vector<std::complex<double>>> cvals = { 
                      {{1,1},{2,2},{3,3}},
                      {{6,6},{7,7},{8,8}} };
    C.set(cvals);
    Matrix<std::complex<double>> D(2,3);
    std::cout<<C.numRows() << " " << C.numCols() << std::endl;
    D = C.Conj();
    std::cout<<D.numRows() << " " << D.numCols() << std::endl;
    for(int i  = 0; i < 2; i++ )
    {
        for(int j = 0; j < 3; j++)
        {
            std::complex<double> x;
            x = std::conj(cvals[i][j]);
            REQUIRE(D(i,j) == x);
        }
    }
    D = C.H();
    for(int i  = 0; i < 3; i++ )
    {
        for(int j = 0; j < 2; j++)
        {
            std::complex<double> x;
            x = std::conj(cvals[j][i]);
            REQUIRE(D(i,j) == x);
        }
    }

} 
TEST_CASE( "elementary ops", "[elementary ops]" )
{
    Matrix<std::complex<double>> C(2,3);
    std::vector<std::vector<std::complex<double>>> cvals = { 
                      {{1,1},{2,2},{3,3}},
                      {{6,6},{7,7},{8,8}} };
    C.set(cvals);
    std::vector<std::complex<double>> col_v_vals = {
        {1,1}, {6,6}
    };
    std::vector<std::complex<double>> row_v_vals = {
        {1,1},{2,2},{3,3}
    };
    Vector<std::complex<double>> CV(2);
    Vector<std::complex<double>> CR(3);
    CV.set(col_v_vals);
    CR.set(row_v_vals);
    REQUIRE(CV == C.getCol(0));
    REQUIRE(CR == C.getRow(0));

    CV[0] = {5,5};
    CV[1] = {11,11};

    Matrix<std::complex<double>> D(2,3);
    std::vector<std::vector<std::complex<double>>> dvals = { 
                      {{2,2},{3,3},{1,1}},
                      {{7,7},{8,8},{6,6}} };
    D.set(dvals);
    std::cout<<D<<std::endl;
    D.swapCols(0,1);
    std::cout<<D<<std::endl;
    D.swapRows(0,1);
    std::cout<<D<<std::endl;
    D.set(dvals);
    D.swapRows(0,1);
    std::vector<std::vector<std::complex<double>>> evals = { 
                      {{7,7},{8,8},{6,6}},
                      {{2,2},{3,3},{1,1}}};
    Matrix<std::complex<double>> E(2,3);
    E.set(evals);
    REQUIRE( D == E );
    
}
TEST_CASE("init", "[init]")
{
    Matrix<std::complex<double>> M(2,2,{{1.514,1},{2,2},{3,3},{4,4}});
    Matrix<std::complex<double>> A({ {{1.514,1},{2,2}},
                                    {{3,3},{4,4}}});
    Matrix<std::complex<double>> B(2,2);
    B = {std::complex<double>{1.514,1},
        std::complex<double>{2,2},
        std::complex<double>{3,3},
        std::complex<double>{4,4}};
    REQUIRE(A==M);
    REQUIRE(B==M);
}

TEST_CASE("extract subblocks", "[subblock]")
{
    Matrix<double> A(5,5);
    int k = 0;
    for(int i = 0; i<5; i++)
    {
        for(int j =0; j<5; j++)
        {
            A(i,j) = k;
            k++;
        }
    }
    Matrix<double> B(2,2);
    B = A.getSubBlock(1,3,1,3);
    std::cout<<A<<std::endl;
    std::cout<<B<<std::endl;
    REQUIRE(B == A.getSubBlock(1,3,1,3));
    Matrix<double> C(2,2);
    A.setSubBlock(1,3,1,3,C);
    std::cout<<A<<std::endl;
    std::cout<<A.getSubBlock(1,2,1,3);
}
TEST_CASE("read csv", "[csv]")
{
    std::string filename;
    filename = "/home/demroz/Documents/code/basic_math/tests/A.txt";
    Matrix<std::complex<double>> A;
    A.read_csv(filename);
    std::cout<<A<<std::endl;
    Vector<std::complex<double>> v;
    v.read_csv(filename);
    std::cout<<v<<std::endl;
       

}






