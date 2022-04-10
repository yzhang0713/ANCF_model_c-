#ifndef ANCF_MODEL_C___TIME_MARCH_H
#define ANCF_MODEL_C___TIME_MARCH_H
#include <Eigen/Dense>
#include <map>
#include <vector>
#include "beam.h"

using namespace Eigen;
using namespace std;

class system_engine;

class time_march {
public:
    static void range_kutta(beam *, double, double, system_engine *);
};


#endif //ANCF_MODEL_C___TIME_MARCH_H
