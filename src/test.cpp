#pragma once
#include <Eigen/Dense>
#include "ndarray.cpp"
#include <iostream>

typedef Eigen::Array<int,3,1> int3;


int main(){
    int shape[] { 10, 3};
    int_2 f(shape);
    f[0][0] = -1;
    std::cout << f[0][0] << std::endl;
    for (auto r: f.range<int3>())
    {
        std::cout << r*r << std::endl;
    }
    return f[0][0];
}