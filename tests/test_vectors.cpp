#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "../header/vectors.hpp"
#include<cmath>

Vector<double> v(5);
double vals[5] = {1,2,3,4,5};
TEST_CASE( "Constructor with array given values", "[constructor]" ) {
    v.set(vals,5);

    for(int i = 0; i < 5; i++)
    {
        REQUIRE( v[i] == vals[i]);
    }

}

TEST_CASE( "Constructor with vector", "[constructor]" ) 
{
    std::vector<double> vals_vec;
    for(int i = 0; i < 5; i++)
    {
        vals_vec.push_back(i);
    }
    v.set(vals_vec);

    for(int i = 0; i < 5; i++)
    {
        REQUIRE( v[i] == vals_vec[i]);
    }

}

TEST_CASE( "set to zero", "[zero]")
{
    v.zero();
    for(int i = 0; i < 5; i++)
    {
        REQUIRE( v[i] == 0.);
    }
}

TEST_CASE( "norm2", "[norm2]")
{
    Vector<double> v1(5);
    for(int i=0;i<5;i++)
    {
        v1[i] = 1.;
    }

    double expected_norm = std::sqrt(5);

    REQUIRE ( v1.norm2() == expected_norm);
}


TEST_CASE( "norm2 complex", "[norm2]")
{
    Vector<std::complex<double>> v1(5);
    for(int i=0;i<5;i++)
    {
        v1[i] = std::complex(1./std::sqrt(2),1/std::sqrt(2));
    }

    double expected_norm = std::sqrt(5);

    REQUIRE ( std::abs(v1.norm2() - expected_norm) < 1e-10);
}

TEST_CASE(" lp norm ", "[lp]")
{
    Vector<std::complex<double>> v1(5);
    for(int i=0;i<5;i++)
    {
        v1[i] = std::complex(1./std::sqrt(2),1/std::sqrt(2));
    }

    double expected_norm = std::sqrt(5);
    REQUIRE ( std::abs(v1.p_norm(2) - expected_norm) < 1e-10);
}


TEST_CASE(" sum ", "[sum]")
{
    Vector<double> v1(5);
    for(int i=0;i<5;i++)
    {
        v1[i] = i;
    }

    double expected_sum = 10;
    REQUIRE ( std::abs(v1.sum() - expected_sum) < 1e-10);
}

TEST_CASE(" dot ", "[dot]")
{
    Vector<double> v1(3);
    Vector<double> v2(3);
    for(int i=0;i<3;i++)
    {
        v1[i] = i;
        v2[i] = i;
    }

    double dot = 5;
    REQUIRE ( v1.dot(v2) == dot );
}

TEST_CASE(" cross ", "[cross]")
{
    Vector<double> v1(3);
    Vector<double> v2(3);
    
    v1[0] = 1;
    v1[1] = 2;
    v1[2] = 3;

    v2[0] = 1;
    v2[1] = 5;
    v2[2] = 7;
    
    Vector<double> v3(3);
    v3[0] = -1;
    v3[1] = -4;
    v3[2] = 3;
    REQUIRE ( v1.cross(v2) == v3 );
}

TEST_CASE(" abs ", "[abs]")
{
    Vector<double> v1(3);
    Vector<double> v2(3);
    
    v1[0] = -1;
    v1[1] = -2;
    v1[2] = -3;

    v2[0] = 1;
    v2[1] = 2;
    v2[2] = 3;
    
    REQUIRE ( v1.abs() == v2 );
}

TEST_CASE(" math ", "[math]")
{
    Vector<double> v1(3);
    Vector<double> v2(3);
    
    v1[0] = 1;
    v1[1] = 2;
    v1[2] = 3;

    v2[0] = 1;
    v2[1] = 5;
    v2[2] = 7;
    
    // addition
    Vector<double> plus(3);
    plus[0] = 2;
    plus[1] = 7;
    plus[2] = 10;

    // subtraction
    Vector<double> minus(3);
    minus[0] = 0;
    minus[1] = -3;
    minus[2] = -4;

    
    REQUIRE ( v1+v2 == plus );

    double dot;
    // dot product via *
    dot = 1*1+2*5+3*7;

    REQUIRE( v1*v2 == dot);
    
    // scalar by vector
    Vector<double> scalar_mult(3);
    scalar_mult[0] = 2;
    scalar_mult[1] = 4;
    scalar_mult[2] = 6;
    
    REQUIRE( (v1 * 2.0) == scalar_mult);
    REQUIRE( (2.0 * v1) == scalar_mult);

    // element wise negative
    Vector<double> neg(3);
    neg[0] = -1;
    neg[1] = -5;
    neg[2] = -7;

    REQUIRE(neg == -v2);
    Vector<double> pe1(3),pe2(3),pe3(3);
    pe1[0] = 1;
    pe1[1] = 2;
    pe1[2] = 3;

    pe2[0] = -1;
    pe2[1] = -2;
    pe2[2] = -3;

    pe3[0] = 0;
    pe3[1] = 0;
    pe3[2] = 0;


    pe1 += pe2;
    REQUIRE(pe1 == pe3);

    pe1 -= pe2;
    REQUIRE(pe1 == -pe2);

}
