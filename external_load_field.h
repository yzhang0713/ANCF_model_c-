#ifndef ANCF_MODEL_C___EXTERNAL_LOAD_FIELD_H
#define ANCF_MODEL_C___EXTERNAL_LOAD_FIELD_H
#include <Eigen/Dense>
#include <vector>
#include "external_load_point.h"

using namespace Eigen;
using namespace std;


class external_load_field {
private:
    vector<external_load_point> forces;
public:
    external_load_field() = default;
    ~external_load_field() = default;
    vector<external_load_point> get_forces() {return forces;};
    void set_forces(vector<external_load_point> fs) {forces = fs;};
};




#endif //ANCF_MODEL_C___EXTERNAL_LOAD_FIELD_H
