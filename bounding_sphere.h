#ifndef ANCF_MODEL_C___BOUNDING_SPHERE_H
#define ANCF_MODEL_C___BOUNDING_SPHERE_H
#include <Eigen/Dense>
#include <vector>

using namespace Eigen;
using namespace std;

class bounding_sphere {
private:
    Vector3d center_point;
    double radius;
    friend int check_disjoint(const bounding_sphere &, const bounding_sphere &);
public:
    bounding_sphere() = default;
    ~bounding_sphere() = default;
    void set_center_point(Vector3d c_point) {this->center_point = c_point;};
    void set_radius(double r) {this->radius = r;};
    int is_bs_disjoint(const bounding_sphere &);
};


#endif //ANCF_MODEL_C___BOUNDING_SPHERE_H
