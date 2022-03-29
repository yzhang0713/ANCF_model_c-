#ifndef ANCF_MODEL_C___FORCE_ENGINE_H
#define ANCF_MODEL_C___FORCE_ENGINE_H
#include <Eigen/Dense>
#include <vector>
#include <math.h>
#include "beam.h"
#include "fluid_field.h"
#include "utils.h"

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
    void add_constraint_load(int, int, VectorXd &);
//    Eigen::MatrixXd point_load(Eigen::MatrixXd, std::vector<beam>);
    void distributed_load(beam *, fluid_field *, double);
};


#endif //ANCF_MODEL_C___FORCE_ENGINE_H
