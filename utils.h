#ifndef ANCF_MODEL_C___UTILS_H
#define ANCF_MODEL_C___UTILS_H
#include<vector>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;


class utils {
public:
    static void gauss_points(int, vector<double>&, vector<double>&);
    static void shape_fun(double, double, int, Eigen::Vector<double, 4>&);
    static Matrix3d tilde(Eigen::Vector3d);
    static Matrix<double, 12, 12> element_mass_matrix(double, double);
};


#endif //ANCF_MODEL_C___UTILS_H
