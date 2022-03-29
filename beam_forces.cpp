#include "beam_forces.h"

void beam_forces::set_Q_gravity(VectorXd Qg) {
    Q_gravity = Qg;
}

VectorXd beam_forces::get_Q_gravity() {
    return Q_gravity;
}

void beam_forces::set_Q_elastic(VectorXd Qe) {
    Q_elastic = Qe;
}

VectorXd beam_forces::get_Q_elastic() {
    return Q_elastic;
}

void beam_forces::set_Q_dist(VectorXd Qd) {
    Q_dist = Qd;
}

VectorXd beam_forces::get_Q_dist() {
    return Q_dist;
}

void beam_forces::set_Q_point(VectorXd Qp) {
    Q_point = Qp;
}

VectorXd beam_forces::get_Q_point() {
    return Q_point;
}