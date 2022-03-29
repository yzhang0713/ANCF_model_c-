#include "fluid_field.h"

Vector3d fluid_field::get_velocity(Vector3d pos, double t) {
    Vector3d v_f = Vector3d::Zero();
    if (pos(2) > (zv - 0.5 * tml) && (pos(0) - Cv * t <= 0.0)) {
        v_f(0) = Uv * cos(kv * M_PI * (pos(0) - Cv * t)) * sin(kv * M_PI * (pos(2) - zv));
        v_f(2) = -Uv * sin(kv * M_PI * (pos(0) - Cv * t)) * cos(kv * M_PI * (pos(2) - zv));
    }
    if (pos(0) - Cv * t > 0.0) {
        v_f(0) = 0.05;
    }
    return v_f;
}

Vector3d fluid_field::get_acceleration(Vector3d pos, double t) {
    Vector3d a_f = Vector3d::Zero();
    if (pos(2) > (zv - 0.5 * tml) && (pos(0) - Cv * t <= 0.0)) {
        a_f(0) = Cv * Uv * sin(kv * M_PI * (pos(0) - Cv * t)) * sin(kv * M_PI * (pos(2) - zv));
        a_f(2) = Cv * Uv * cos(kv * M_PI * (pos(0) - Cv * t)) * cos(kv * M_PI * (pos(2) - zv));
    }
    return a_f;
}

ostream& operator<<(ostream& os, const fluid_field& ff) {
    os << "Fluid field details: " << endl;
    os << "    Uv - " << ff.Uv << endl;
    os << "    Cv - " << ff.Cv << endl;
    os << "    kv - " << ff.kv << endl;
    os << "    zv - " << ff.zv << endl;
    os << "    rhof - " << ff.rhof << endl;
    os << "    Cm - " << ff.Cm << endl;
    os << "    Cd - " << ff.Cd << endl;
    return os;
}

fluid_field_builder::fluid_field_builder() {
    this->Reset();
}

fluid_field_builder::~fluid_field_builder() {
    delete ff;
}

void fluid_field_builder::Reset() {
    this->ff = new fluid_field();
}

void fluid_field_builder::set_Uv(double Uv) {
    this->ff->Uv = Uv;
}

void fluid_field_builder::set_Cv(double Cv) {
    this->ff->Cv = Cv;
}

void fluid_field_builder::set_kv(double kv) {
    this->ff->kv = kv;
}

void fluid_field_builder::set_zv(double zv) {
    this->ff->zv = zv;
}

void fluid_field_builder::set_rhof(double rhof) {
    this->ff->rhof = rhof;
}

void fluid_field_builder::set_Cm(double Cm) {
    this->ff->Cm = Cm;
}

void fluid_field_builder::set_Cd(double Cd) {
    this->ff->Cd = Cd;
}

void fluid_field_builder::set_tml(double tml) {
    this->ff->tml = tml;
}

fluid_field * fluid_field_builder::get_fluid_field() {
    fluid_field * result = this->ff;
    this->Reset();
    return result;
}