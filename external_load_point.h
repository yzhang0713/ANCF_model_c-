#ifndef ANCF_MODEL_C___EXTERNAL_LOAD_POINT_H
#define ANCF_MODEL_C___EXTERNAL_LOAD_POINT_H
#include <Eigen/Dense>

using namespace Eigen;

class external_load_point {
private:
    double position;
    Vector3d force;
public:
    external_load_point() = default;
    ~external_load_point() = default;
    void set_position(double pos) {position = pos;};
    void set_force(Vector3d f) {force = f;};
    double get_position() {return position;};
    Vector3d get_force() {return force;};
};

#endif //ANCF_MODEL_C___EXTERNAL_LOAD_POINT_H
