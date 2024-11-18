#ifndef __UTILS_HPP__
#define __UTILS_HPP__
#include <vector>
#include <cmath>
#include <iostream>
#include <complex>
#include <limits>
#include <assert.h>
double conj(double x);
float conj(float x);
int conj(int x);
std::complex<double> conj(std::complex<double> x);
std::complex<float> conj(std::complex<float> x);
std::complex<int> conj(std::complex<int> x);

double arg(double x);
float arg(float x);
int arg(int x);
std::complex<double> arg(std::complex<double> x);
std::complex<float> arg(std::complex<float> x);
std::complex<int> arg(std::complex<int> x);

double computeHouseHolderAlphaPhase(double x);
float computeHouseHolderAlphaPhase(float x);
int computeHouseHolderAlphaPhase(int x);
std::complex<double> computeHouseHolderAlphaPhase(std::complex<double> x);
std::complex<float> computeHouseHolderAlphaPhase(std::complex<float> x);
std::complex<int> computeHouseHolderAlphaPhase(std::complex<int> x);
std::string convert_complex_string_to_streamstring(std::string s);

std::string convert_complex_string_to_streamstring(std::string s)
{
    bool isScientificNotation = false;
    isScientificNotation = (s.find("e") != std::string::npos);
    std::string ret;
    if(isScientificNotation)
    {
        int first_e_index = s.find_first_of("e");
        std::string real, imag;
        real = s.substr(0,first_e_index+2);
        imag = s.substr(first_e_index + 2, s.size() - first_e_index - 2);
        std::string imaginary_front;
        size_t pos = imag.find_first_of("+-");
        
        imaginary_front = imag.substr(0,pos);
        real += imaginary_front;
        std::cout<<imag<<std::endl;
        imag.erase(0,pos);
        std::cout<<imag<<std::endl;
        if (imag[0] == '+')
        {
            imag[0] = ',';
        }
        else
        {
            imag = "," + imag;
        }
        ret = real+imag;
    }
    else
    {
        size_t pos = s.find_last_of("+-");
        std::cout<<pos<<" " <<s[pos] << std::endl;
        if (s[pos] == '+')
        {
            s[pos] = ',';
        }
        else
        {
            s.insert(pos,",");
        }
        std::cout<<s<<std::endl;
        ret = s;
    }
    ret.erase(std::remove(ret.begin(), ret.end(), 'j'), ret.end());
    return ret;
}

/*
**********************************************
* ***** overloaded conjugate functions *******
**********************************************
*/
double conj(double x)
{
    return x;
}

float conj(float x)
{
    return x;
}

int conj(int x)
{
    return x;
}


std::complex<double> conj(std::complex<double> x)
{
    return std::conj(x);
}

std::complex<float> conj(std::complex<float> x)
{
    return std::conj(x);
}

std::complex<int> conj(std::complex<int> x)
{
    return std::conj(x);
}

double arg(double x)
{
    return 0;
}
float arg(float x)
{
    return 0;
}
int arg(int x)
{
    return 0;
}
std::complex<double> arg(std::complex<double> x)
{
    return std::arg(x);
}
std::complex<float> arg(std::complex<float> x)
{
    return std::arg(x);
}
std::complex<int> arg(std::complex<int> x)
{
    return std::arg(x);
}

inline double computeHouseHolderAlphaPhase(double x)
{
    if (x >= 0)
    {
        return 1.0;
    }
    return -1.0;
}

inline float computeHouseHolderAlphaPhase(float x)
{
    if (x >= 0)
    {
        return 1.0f;
    }
    return -1.0f;
}

inline int computeHouseHolderAlphaPhase(int x)
{
    if (x >= 0)
    {
        return 1;
    }
    return -1;
}

std::complex<double> computeHouseHolderAlphaPhase(std::complex<double> x)
{
    return 1.0*std::exp( std::complex<double>{0,1} * arg(x));
}
std::complex<float> computeHouseHolderAlphaPhase(std::complex<float> x)
{
    return 1.0f*std::exp( std::complex<float>{0,1} * arg(x));
}
std::complex<int> computeHouseHolderAlphaPhase(std::complex<int> x)
{
    return 1*std::exp( std::complex<int>{0,1} * arg(x));
}


#endif