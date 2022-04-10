#ifndef ANCF_MODEL_C___SYSTEM_ENGINE_H
#define ANCF_MODEL_C___SYSTEM_ENGINE_H
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream>
#include "beam.h"
#include "particle.h"
#include "force_engine.h"
#include "utils.h"
#include "time_march.h"
#include "fluid_field.h"
#include <map>

using namespace Eigen;
using namespace std;

class system_engine {
private:
    vector<beam*> beams;
    vector<particle*> particles;
    force_engine* f_engine;
    fluid_field* f_field;
    double t;
    double h;
    map<int, vector<VectorXd>> beam_pos_record;
    map<int, vector<VectorXd>> beam_vel_record;
    map<int, vector<VectorXd>> beam_acc_record;

public:
    system_engine();
    ~system_engine() = default;
    force_engine * get_force_engine() {return f_engine;};
    fluid_field * get_fluid_field() {return f_field;};
    void read_beams();
    void test_read_beams();
    void read_particles();
    void test_read_particles();
    double get_cur_time() {return t;};
    void update_time() {t += h;};
    map<int, vector<VectorXd>> & get_beam_pos_record() {return beam_pos_record;};
    map<int, vector<VectorXd>> & get_beam_vel_record() {return beam_vel_record;};
    map<int, vector<VectorXd>> & get_beam_acc_record() {return beam_acc_record;};
    void read_fluid_field();
    void test_read_fluid_field();
    void update_beam_gravity_forces();
    void update_beam_elastic_forces();
    void update_beam_dist_forces();
    void reset_beam_point_forces();
    void update_beam_point_forces();
    void update_beam_total_forces();
    void update_beam_forces_with_constraint();
    VectorXd solving_linear_system_of_beam(beam *);
    void store_beam_information();
    void update_beams();
    void check_output_folder();
    void write_to_file();
};


#endif //ANCF_MODEL_C___SYSTEM_ENGINE_H
