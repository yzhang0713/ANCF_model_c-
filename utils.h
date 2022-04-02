#ifndef ANCF_MODEL_C___UTILS_H
#define ANCF_MODEL_C___UTILS_H
#include <vector>
#include <Eigen/Dense>
#include "beam.h"

using namespace std;
using namespace Eigen;


class utils {
public:
    static void gauss_points(int, vector<double>&, vector<double>&);
    static void shape_fun(double, double, int, Eigen::Vector<double, 4>&);
    static Matrix3d tilde(Eigen::Vector3d);
    static Matrix<double, 12, 12> element_mass_matrix(double, double);
    static vector<Vector3d> get_points_from_beam(beam*, int, int, int);
    static vector<Vector3d> get_points_from_beam_element(beam*, int, vector<double>);
};


#endif //ANCF_MODEL_C___UTILS_H
