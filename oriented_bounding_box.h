#ifndef ANCF_MODEL_C___ORIENTED_BOUNDING_BOX_H
#define ANCF_MODEL_C___ORIENTED_BOUNDING_BOX_H
#include <Eigen/Dense>
#include <vector>

using namespace Eigen;
using namespace std;

class oriented_bounding_box {
private:
    Vector3d mean_point;
    vector<Vector3d> normal_dir;
    vector<Vector2d> range_dir;
    friend int check_disjoint(const oriented_bounding_box &, const oriented_bounding_box &, int, double);
public:
    oriented_bounding_box();
    ~oriented_bounding_box() = default;
    void set_obb_info(vector<Vector3d>);
    int is_obb_disjoint(const oriented_bounding_box &, double);
};


#endif //ANCF_MODEL_C___ORIENTED_BOUNDING_BOX_H
