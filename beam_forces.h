#ifndef ANCF_MODEL_C___BEAM_FORCES_H
#define ANCF_MODEL_C___BEAM_FORCES_H
#include <Eigen/Dense>

using namespace Eigen;


class beam_forces {
private:
    VectorXd Q_gravity;         // Gravity force
    VectorXd Q_elastic;         // Elastic force
    VectorXd Q_dist;            // Distributed external force
    VectorXd Q_point;           // Concentrated external force
public:
    beam_forces() = default;
    ~beam_forces() = default;
    void set_Q_gravity(VectorXd Qg);
    VectorXd get_Q_gravity();
    void set_Q_elastic(VectorXd Qe);
    VectorXd get_Q_elastic();
    void set_Q_dist(VectorXd Qd);
    VectorXd get_Q_dist();
    void set_Q_point(VectorXd Qp);
    VectorXd get_Q_point();
};


#endif //ANCF_MODEL_C___BEAM_FORCES_H
