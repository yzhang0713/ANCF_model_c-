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

// Provide the element level mass matrix
Matrix<double, 12, 12> utils::element_mass_matrix(double L, double factor) {
    double m11 = factor * 13.0 / 35.0;
    double m12 = factor * L * 11.0 / 210.0;
    double m13 = factor * 9.0 / 70.0;
    double m14 = - factor * L * 13.0 / 420.0;
    double m22 = factor * L * L / 105.0;
    double m23 = - m14;
    double m24 = - factor * L * L / 140.0;
    double m33 = m11;
    double m34 = - m12;
    double m44 = m22;
    Matrix<double, 12, 12> mass_matrix;
    mass_matrix.setZero();
    for (int i = 0; i < 3; i++) {
        mass_matrix(i, i) = m11;
        mass_matrix(i+3, i) = m12;
        mass_matrix(i+6,i) = m13;
        mass_matrix(i+9,i) = m14;
    }
    for (int i = 3; i < 6; i++) {
        mass_matrix(i-3,i) = m12;
        mass_matrix(i,i) = m22;
        mass_matrix(i+3,i) = m23;
        mass_matrix(i+6,i) = m24;
    }
    for (int i = 6; i < 9; i++) {
        mass_matrix(i-6,i) = m13;
        mass_matrix(i-3,i) = m23;
        mass_matrix(i,i) = m33;
        mass_matrix(i+3,i) = m34;
    }
    for (int i = 9; i < 12; i++) {
        mass_matrix(i-9,i) = m14;
        mass_matrix(i-6,i) = m24;
        mass_matrix(i-3,i) = m34;
        mass_matrix(i,i) = m44;
    }
    return mass_matrix;
}

vector<Vector3d> get_points_from_beam(beam* b, int start_element, int n_element, int n_point) {
    vector<Vector3d> points;
    if (n_point <= 1) {
        return points;
    }
    double l_element = b->get_length() / ((double) b->get_nelement());
    VectorXd pos = b->get_position();
    for (int i = 0; i < n_point; i++) {
        int k = start_element + (n_element*i)/(n_point-1);
        if ((n_element*i)%(n_point-1) == 0) {
            points.push_back(pos.segment(6*(k-1),3));
        } else {
            double xi = ((double) (n_element*i)) / ((double) (n_point-1)) - k;
            Vector<double, 4> S;
            utils::shape_fun(xi*l_element, l_element, 0, S);
            Matrix<double, 3, 4> disp_rs = pos.segment(6*(k-1), 12).reshaped(3, 4);
            points.push_back(disp_rs*S);
        }
    }
    return points;
}

vector<Vector3d> utils::get_points_from_beam_element(beam* b, int i_element, vector<double> positions) {

}