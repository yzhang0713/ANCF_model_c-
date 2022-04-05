#ifndef ANCF_MODEL_C___FORCE_ENGINE_H
#define ANCF_MODEL_C___FORCE_ENGINE_H
#include <Eigen/Dense>
#include <vector>
#include <map>
#include <math.h>
#include "beam.h"
#include "fluid_field.h"
#include "utils.h"
#include "oriented_bounding_box.h"

using namespace Eigen;
using namespace std;

class force_engine {
private:
    Vector3d gravity{0.0, 0.0, 9.8};
    Vector<double, 12> axial_force_integrand(Vector<double, 12>, double, double);
    Vector<double, 12> flexural_force_integrand(Vector<double, 12>, double, double);
public:
    force_engine() = default;
    ~force_engine() = default;
    void gravity_force(beam *);
    void elastic_force(beam *);
    void add_constraint_load(beam *);
    void distributed_load(beam *, fluid_field *, double);
    void point_load_element_level(beam *, int, beam *, int);
    void point_load(beam*, beam*);
    void point_load_segment_level(beam*, int, int, beam*, int, int);
};


#endif //ANCF_MODEL_C___FORCE_ENGINE_H
