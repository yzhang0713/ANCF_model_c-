#include "utils.h"
#include <math.h>
#include <iostream>

// This method will return the gauss points for numerical integration
void utils::gauss_points(int n, vector<double> &x, vector<double> &weight) {
    if (n < 2 || n > 5) {
        std::cout << "Gauss points are only defined here in range 2 to 5.\n";
        return;
    }
    switch (n) {
        case 2: {
            x[0] = - sqrt(1.0/3.0);
            x[1] = - x[0];
            weight[0] = 1.0;
            weight[1] = 1.0;
            break;
        }
        case 3: {
            x[0] = - sqrt(3.0/5.0);
            x[1] = 0.0;
            x[2] = - x[0];
            weight[0] = 5.0/9.0;
            weight[1] = 8.0/9.0;
            weight[2] = weight[0];
            break;
        }
        case 4: {
            double value1 = 3.0/7.0;
            double value2 = sqrt(120.0)/35.0;
            x[0] = - sqrt(value1 + value2);
            x[1] = - sqrt(value1 - value2);
            x[2] = - x[1];
            x[3] = - x[0];
            value1 = 5.0/3.0/sqrt(120.0);
            weight[0] = 0.5 - value1;
            weight[1] = 0.5 + value1;
            weight[2] = weight[1];
            weight[3] = weight[0];
            break;
        }
        case 5: {
            double value1 = 2.0*sqrt(10.0/7.0);
            x[0] = - sqrt(5.0 + value1) / 3.0;
            x[1] = - sqrt(5.0 - value1) / 3.0;
            x[2] = 0.0;
            x[3] = - x[1];
            x[4] = - x[0];
            value1 = 13.0/sqrt(70.0);
            weight[0] = (322.0 - value1) / 900.0;
            weight[1] = (322.0 + value1) / 900.0;
            weight[2] = 128.0/225.0;
            weight[3] = weight[1];
            weight[4] = weight[0];
        }
    }
}

// Generate the shape function matrix, and also the derivatives of shape function based on
// third input of function
//      der = 0, return S
//      der = 1, return S_x
//      der = 2, return S_xx
void utils::shape_fun(double x, double L, int der, Vector4d &S) {
    double xi = x / L;
    switch (der) {
        case 0: {
            S(0) = 1.0 - 3.0 * xi * xi + 2.0 * xi * xi * xi;
            S(1) = L * (xi - 2.0 * xi * xi + xi * xi * xi);
            S(2) = 3.0 * xi * xi - 2.0 * xi * xi * xi;
            S(3) = L * (- xi * xi + xi * xi * xi);
        }
        case 1: {
            S(0) = (- 6.0 * xi + 6.0 * xi * xi) / L;
            S(1) = 1.0 - 4.0 * xi + 3.0 * xi * xi;
            S(2) = (6.0 * xi - 6.0 * xi * xi) / L;
            S(3) = -2.0 * xi + 3.0 * xi * xi;
        }
        case 2: {
            S(0) = (- 6.0 + 12.0 * xi) / (L * L);
            S(1) = (- 4.0 + 6.0 * xi) / L;
            S(2) = (6.0 - 12.0 * xi) / (L * L);
            S(3) = (-2.0 + 6.0 * xi) / L;
        }
    }
}

// Provide the matrix form of vector to help compute the cross product
Matrix3d utils::tilde(Vector3d vec) {
    Matrix3d vec_mat;
    vec_mat << 0.0, -vec(2), vec(1), vec(2), 0.0, -vec(0), -vec(1), vec(0), 0.0;
    return vec_mat;
}