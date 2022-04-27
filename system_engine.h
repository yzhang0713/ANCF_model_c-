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
#include "external_load_field.h"
#include "external_load_point.h"
#include <map>

using namespace Eigen;
using namespace std;

class system_engine {
private:
    vector<beam*> beams;
    vector<particle*> particles;
    force_engine* f_engine;
    fluid_field* f_field;
    external_load_field* el_field;
    beam_builder * b_builder;
    particle_builder * p_builder;
    fluid_field_builder * ff_builder;
    double t;
    double h;
    map<int, vector<VectorXd>> beam_pos_record;
    map<int, vector<VectorXd>> beam_vel_record;
    map<int, vector<VectorXd>> beam_acc_record;
    map<int, vector<VectorXd>> beam_elastic_record;
    map<int, vector<VectorXd>> beam_gravity_record;
    map<int, vector<VectorXd>> beam_damping_record;
    map<int, vector<VectorXd>> beam_external_record;
    map<int, vector<VectorXd>> beam_distribute_record;
    map<int, vector<VectorXd>> beam_point_record;
    vector<double> time_record;
    int distributed_force_switch;
    int point_force_switch;
    int external_load_switch;
    int damping_force_switch;
    int gravity_force_switch;
    int debug = 0;
public:
    system_engine();
    ~system_engine();
    force_engine * get_force_engine() {return f_engine;};
    fluid_field * get_fluid_field() {return f_field;};
    void initialize_time();
    void set_current_time(double);
    void set_time_step(double);
    void read_beams();
    void test_read_beams();
    void read_particles();
    void test_read_particles();
    void read_external_load_field();
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
    void initialize_beam_point_forces();
    void reset_beam_point_forces();
    void debug_beam_point_forces();
    void update_beam_point_forces();
    void update_beam_external_forces();
    void update_beam_damping_forces();
    void update_beam_total_forces();
    void update_beam_forces_with_constraint();
    VectorXd solving_linear_system_of_beam(beam *);
    void store_beam_information();
    void update_beams();
    void check_output_folder();
    void write_to_file();
    int check_write();
    void load_switches();
    void beam_dist_force_off();
    void beam_point_force_off();
    void beam_external_force_off();
    void beam_damping_force_off();
    void beam_gravity_force_off();
    int get_distributed_force_switch() {return distributed_force_switch;};
    int get_point_force_switch() {return point_force_switch;};
    int get_external_force_switch() {return external_load_switch;};
    int get_damping_force_switch() {return damping_force_switch;};
    int get_gravity_force_switch() {return gravity_force_switch;};
    void remove_outputs();
};


#endif //ANCF_MODEL_C___SYSTEM_ENGINE_H
