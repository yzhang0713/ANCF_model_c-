#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#include "particle.h"
#include <math.h>

void particle::set_volume() {
    volume = 4.0 / 3.0 * M_PI * r * r * r;
}

ostream& operator<<(ostream& os, particle& p) {
    os << "Particle details: " << endl;
    os << "    E - " << p.E << endl;
    os << "    rho - " << p.rho << endl;
    os << "    r - " << p.r << endl;
    return os;
}

particle_builder::particle_builder() {
    this->Reset();
}

particle_builder::~particle_builder() {
    delete p;
}

void particle_builder::Reset() {
    this->p = new particle();
}

void particle_builder::set_E(double E) {
    this->p->E = E;
}

void particle_builder::set_rho(double rho) {
    this->p->rho = rho;
}

void particle_builder::set_r(double r) {
    this->p->r = r;
    this->p->set_volume();
}

void particle_builder::set_position(Vector3d position) {
    this->p->position = position;
}

void particle_builder::set_velocity() {
    this->p->velocity << 0.0, 0.0, 0.0;
}

void particle_builder::set_velocity(Vector3d velocity) {
    this->p->velocity = velocity;
}

void particle_builder::set_acceleration() {
    this->p->acceleration << 0.0, 0.0, 0.0;
}

void particle_builder::set_acceleration(Vector3d acceleration) {
    this->p->acceleration = acceleration;
}

particle *particle_builder::get_particle() {
    particle* result = this->p;
    this->Reset();
    return result;
}
