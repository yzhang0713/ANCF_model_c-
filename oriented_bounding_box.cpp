#include "oriented_bounding_box.h"
#include <iostream>

oriented_bounding_box::oriented_bounding_box() {
    for (int i = 0; i < 3; i++) {
        normal_dir.push_back(Vector3d::Zero());
        range_dir.push_back(Vector2d::Zero());
    }
}

void oriented_bounding_box::set_obb_info(vector<Vector3d> points) {
    int debug = 0;
    // Number of points
    double n_point = (double) points.size();
    if (debug)
        cout << "Number of points: " << n_point << endl;
    // Mean position
    mean_point.setZero();
    for (Vector3d point : points) {
        if (debug)
            cout << "Current point: " << point.transpose() << endl;
        mean_point += point;
    }
    mean_point /= n_point;
    if (debug)
        cout << "Mean point: " << mean_point.transpose() << endl;
    // Normal directions
    double cov_11 = 0.0;
    double cov_33 = 0.0;
    double cov_13 = 0.0;
    for (Vector3d point : points) {
        double diff_1 = point(0) - mean_point(0);
        double diff_3 = point(2) - mean_point(2);
        cov_11 += (diff_1 * diff_1);
        cov_13 += (diff_1 * diff_3);
        cov_33 += (diff_3 * diff_3);
    }
    cov_11 /= n_point;
    cov_13 /= n_point;
    cov_33 /= n_point;
    if (debug)
        cout << "Covariance done: " << cov_11 << ", " << cov_33 << ", " << cov_13 << endl;
    double factor = sqrt((cov_33-cov_11)*(cov_33-cov_11)+4.0*cov_13*cov_13);
    double factor_1 = ((cov_33-cov_11) + factor) / 2.0;
    double factor_2 = factor_1 - factor;
    double factor_3 = sqrt(factor_1 * factor_1 + cov_13 * cov_13);
    double factor_4 = sqrt(factor_2 * factor_2 + cov_13 * cov_13);
    if (abs(cov_13) <= 1.0e-5) {
        if (debug)
            cout << "case 1" << endl;
        normal_dir[0] << 1.0, 0.0, 0.0;
        normal_dir[2] << 0.0, 0.0, 1.0;
    } else {
        normal_dir[0] << cov_13/factor_3, 0.0, factor_1/factor_3;
        normal_dir[2] << cov_13/factor_4, 0.0, factor_2/factor_4;
    }
    normal_dir[1].setZero();
    if (debug) {
        cout << "normal direction 1: " << normal_dir[0].transpose() << endl;
        cout << "normal direction 1: " << normal_dir[1].transpose() << endl;
    }
    // Range of box over normal directions
    range_dir[0].setZero();
    range_dir[1].setZero();
    range_dir[2].setZero();
    for (Vector3d point : points) {
        double dot_prod = normal_dir[0].dot(point-mean_point);
        range_dir[0](0) = min(range_dir[0](0), dot_prod);
        range_dir[0](1) = max(range_dir[0](1), dot_prod);
        dot_prod = normal_dir[2].dot(point-mean_point);
        range_dir[2](0) = min(range_dir[2](0), dot_prod);
        range_dir[2](1) = max(range_dir[2](1), dot_prod);
    }
}

int check_disjoint(const oriented_bounding_box & obb1, const oriented_bounding_box & obb2, int base_dir, double thick) {
    // Get base vector, base range, and distance between mean points
    Vector3d base_vector = obb1.normal_dir[base_dir];
    if (base_vector.norm() < 1.0e-5) {
        return 0;
    }
    Vector2d base_range = obb1.range_dir[base_dir];
    Vector3d mean_vector = obb2.mean_point - obb1.mean_point;
    double distance = base_vector.dot(mean_vector);
    double r_base = (distance > 0) ? base_range(1) : -base_range(0);
    double r = 0.0;
    if (distance > 0) {
        r_base = base_range(1);
        for (int i = 0; i < 3; i++) {
            double dot_prod = base_vector.dot(obb2.normal_dir[i]);
            r -= (dot_prod * ((dot_prod > 0) ? obb2.range_dir[i](0) : obb2.range_dir[i](1)));
        }
    } else {
        distance *= -1.0;
        r_base = -base_range(0);
        for (int i = 0; i < 3; i++) {
            double dot_prod = base_vector.dot(obb2.normal_dir[i]);
            r += (dot_prod * ((dot_prod > 0) ? obb2.range_dir[i](1) : obb2.range_dir[i](0)));
        }
    }
    return (distance - r_base - r - thick) > 0;
}

int oriented_bounding_box::is_obb_disjoint(const oriented_bounding_box & other, double thick) {
    for (int i = 0; i < 3; i++) {
        if (check_disjoint(*this, other, i, thick)) {
            return 1;
        }
        if (check_disjoint(other, *this, i, thick)) {
            return 1;
        }
    }
    return 0;
}