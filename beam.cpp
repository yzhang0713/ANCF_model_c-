#include "beam.h"

void beam::set_area() {
    area = thick * width;
}

void beam::set_inertia() {
    inertia = width * thick * thick * thick / 12.0;
}

void beam::set_ndof() {
    ndof = (nelement + 1) * 6;
}

ostream& operator<<(ostream& os, const beam& b) {
    os << "Beam details: " << endl;
    os << "    E - " << b.E << endl;
    os << "    rho - " << b.rho << endl;
    os << "    length - " << b.length << endl;
    os << "    thick - " << b.thick << endl;
    os << "    width - " << b.width << endl;
    os << "    element - " << b.nelement << endl;
    return os;
}

beam_builder::beam_builder() {
    this->Reset();
}

beam_builder::~beam_builder() {
    delete b;
}

void beam_builder::Reset() {
    this->b = new beam();
}

void beam_builder::set_E(double E) {
    this->b->E = E;
}

void beam_builder::set_rho(double rho) {
    this->b->rho = rho;
}

void beam_builder::set_length(double length) {
    this->b->length = length;
}

void beam_builder::set_thick(double thick) {
    this->b->thick = thick;
}

void beam_builder::set_width(double width) {
    this->b->width = width;
}

void beam_builder::set_nelement(int nelement) {
    this->b->nelement = nelement;
    this->b->set_ndof();
}

void beam_builder::set_botCnstr(int botCnstr) {
    this->b->botCnstr = botCnstr;
}

void beam_builder::set_topCnstr(int topCnstr) {
    this->b->topCnstr = topCnstr;
}

void beam_builder::set_position(VectorXd position) {
    this->b->position = position;
}

void beam_builder::set_velocity() {
    this->b->velocity = VectorXd::Zero(this->b->ndof);
}

void beam_builder::set_velocity(VectorXd velocity) {
    this->b->velocity = velocity;
}

void beam_builder::set_acceleration() {
    this->b->acceleration = VectorXd::Zero(this->b->ndof);
}

void beam_builder::set_acceleration(VectorXd acceleration) {
    this->b->acceleration = acceleration;
}

beam *beam_builder::get_beam() {
    this->b->set_area();
    this->b->set_inertia();
    this->b->forces = new beam_forces();
    beam* result = this->b;
    this->Reset();
    return result;
}