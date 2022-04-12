#include "beam.h"

int beam::counter = 1;

beam::beam(const beam & b) {
    E = b.E;
    nu = b.nu;
    rho = b.rho;
    length = b.length;
    thick = b.thick;
    width = b.width;
    area = b.area;
    inertia = b.inertia;
    nelement = b.nelement;
    ndof = b.ndof;
    botCnstr = b.botCnstr;
    topCnstr = b.topCnstr;
    position = b.position;
    velocity = b.velocity;
    mass_LU = b.mass_LU;
    forces = new beam_forces();
}

void beam::set_beam_id() {
    beam_id = counter;
    counter++;
}

void beam::set_area() {
    area = thick * width;
}

void beam::set_inertia() {
    inertia = width * thick * thick * thick / 12.0;
}

void beam::set_ndof() {
    ndof = (nelement + 1) * 6;
}

void beam::set_mass_matrix() {
    // Get element mass matrix
    double l_element = length / nelement;
    double factor = rho * area * l_element;
    Matrix<double, 12, 12> mass_element = utils::element_mass_matrix(l_element, factor);
//    cout << "mass element done" << endl;
    // Assemble to get the global mass matrix of the beam
    mass_matrix = MatrixXd::Zero(ndof, ndof);
//    cout << "ndof " << ndof << endl;
//    cout << "nelement " << nelement << endl;
    mass_matrix.setZero();
    for (int i = 0; i < nelement; i++) {
        mass_matrix.block(6*i,6*i,6,6) += mass_element.topLeftCorner(6,6);
        mass_matrix.block(6*(i+1),6*i,6,6) = mass_element.bottomLeftCorner(6,6);
        mass_matrix.block(6*i,6*(i+1),6,6) = mass_element.topRightCorner(6,6);
        mass_matrix.block(6*(i+1),6*(i+1),6,6) += mass_element.bottomRightCorner(6,6);
    }
    // Adjust mass matrix according to constraint
    if (botCnstr != 0) {
        mass_matrix.topRows(3*botCnstr).setZero();
        for (int i = 0; i < 3*botCnstr; i++) {
            mass_matrix(i,i) = 1.0;
        }
    }
    if (topCnstr != 0) {
        mass_matrix.bottomRows(3*topCnstr).setZero();
        for (int i = 0; i < 3*topCnstr; i++) {
            mass_matrix(ndof-i,ndof-i) = 1.0;
        }
    }
//    cout << "mass matrix done" << endl;
}

void beam::set_mass_LU() {
    mass_LU = mass_matrix.fullPivLu();
}

void beam::set_total_force() {
    forces->set_Q_total(forces->get_Q_gravity() + forces->get_Q_dist()
                        + forces->get_Q_point() + forces->get_Q_external()
                        - forces->get_Q_elastic() + forces->get_Q_damping());
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

void beam_builder::set_nu(double nu) {
    this->b->nu = nu;
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
    cout << "beam velocity: " << this->b->velocity << endl;
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
    this->b->set_beam_id();
    this->b->set_area();
    this->b->set_inertia();
    this->b->set_mass_matrix();
    this->b->set_mass_LU();
    this->b->forces = new beam_forces();
    beam* result = this->b;
    this->Reset();
    return result;
}