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

using namespace Eigen;
using namespace std;

class system_engine {
private:
    vector<beam*> beams;
    vector<particle*> particles;
    force_engine* f_engine;
    fluid_field* f_field;
    double t;
public:
    system_engine();
    ~system_engine() = default;
    void read_beams();
    void test_read_beams();
    void read_particles();
    void test_read_particles();
    void read_fluid_field();
    void test_read_fluid_field();
    void update_beam_gravity_forces();
    void update_beam_elastic_forces();
    void update_beam_dist_forces();
};


#endif //ANCF_MODEL_C___SYSTEM_ENGINE_H
