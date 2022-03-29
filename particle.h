#ifndef ANCF_MODEL_C___PARTICLE_H
#define ANCF_MODEL_C___PARTICLE_H
#include <Eigen/Dense>
#include <iostream>

using namespace std;
using namespace Eigen;

class particle {
private:
    double E = 4.5e10;
    double rho = 2.65e3;
    double r = 0.01;
    double volume;
    Vector3d position;
    Vector3d velocity;
    Vector3d acceleration;
    friend class particle_builder;
public:
    particle() = default;
    ~particle() = default;
    void set_volume();
    double get_E() {return E;}
    double get_rho() {return rho;}
    double get_r() {return r;}
    double get_volume() {return volume;}
    Vector3d get_position() {return position;}
    Vector3d get_velocity() {return velocity;}
    Vector3d get_acceleration() {return acceleration;}
    friend ostream& operator<<(ostream&, particle&);
};

class particle_builder {
private:
    particle* p;
public:
    particle_builder();
    ~particle_builder();
    void Reset();
    void set_E(double);
    void set_rho(double);
    void set_r(double);
    void set_position(Vector3d);
    void set_velocity();
    void set_velocity(Vector3d);
    void set_acceleration();
    void set_acceleration(Vector3d);
    particle* get_particle();
};


#endif //ANCF_MODEL_C___PARTICLE_H
