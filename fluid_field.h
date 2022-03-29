#ifndef ANCF_MODEL_C___FLUID_FIELD_H
#define ANCF_MODEL_C___FLUID_FIELD_H
#include <iostream>
#include <Eigen/Dense>
#include <math.h>

using namespace std;
using namespace Eigen;

class fluid_field {
private:
    double Uv;
    double Cv;
    double kv;
    double zv;
    double rhof;
    double Cm;
    double Cd;
    double tml;
    friend class fluid_field_builder;
public:
    fluid_field() = default;
    ~fluid_field() = default;
    double get_rhof() {return rhof;};
    double get_Cm() {return Cm;};
    double get_Cd() {return Cd;};
    Vector3d get_velocity(Vector3d, double);
    Vector3d get_acceleration(Vector3d, double);
    friend ostream& operator<<(ostream&, const fluid_field&);
};

class fluid_field_builder {
private:
    fluid_field * ff;
public:
    fluid_field_builder();
    ~fluid_field_builder();
    void Reset();
    void set_Uv(double);
    void set_Cv(double);
    void set_kv(double);
    void set_zv(double);
    void set_rhof(double);
    void set_Cm(double);
    void set_Cd(double);
    void set_tml(double);
    fluid_field * get_fluid_field();
};

#endif //ANCF_MODEL_C___FLUID_FIELD_H
