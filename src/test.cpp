#include <Eigen/Dense>


typedef Eigen::Array<float,3,3> float33;	//used for triangle geometry computations


int main(){
    float33 f;
    f[0][0] = -1;
    return f[0][0];
}