#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "../header/linalg.hpp"
#include "../header/linearRegression.hpp"
#include <cmath>


TEST_CASE("readcsv", "readcsv")
{
    std::string filename;
    filename = "/home/demroz/Documents/code/basic_math/tests/test_matrices/regression_test_data.csv";
    std::vector<std::string> xVars;
    std::string yVars;
    xVars.push_back("x1");
    xVars.push_back("x2");
    yVars = "y";
    LR::read_least_squares_data_from_csv<double>(filename,xVars, yVars);
}