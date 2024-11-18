#include <iostream>
#include "../header/vectors.hpp"

int main()
{

    Vector<double> p(3);
    p[0] = 1;
    p[1] = 2;
    p[2] = 3;


    Vector<double> q(3);
    q[0] = 1;
    q[1] = 2;
    q[2] = 3;
    
    p += q;

    std::cout << p << std::endl;
}