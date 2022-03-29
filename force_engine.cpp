#include "force_engine.h"

void force_engine::gravity_force(beam * b) {
    // Beam parameters
    int ndof = b->get_ndof();
    int nelement = b->get_nelement();
    double lelement = b->get_length() / ((double) nelement);
    double rho = b->get_rho();
    double area = b->get_area();

    // Initialize gravity load vector
    VectorXd Qg = VectorXd::Zero(ndof);

    // Local gravity vector
    Vector<double, 12> Qg_local;
    double factors[4];
    factors[0] = rho * area * lelement / 2.0;
    factors[1] = rho * area * lelement * lelement / 12.0;
    factors[2] = factors[0];
    factors[3] = -factors[1];
    for (int i = 0; i < 4; i++) {
        Qg_local.middleRows(3*i, 3) = factors[i] * gravity;
    }

    // Construct global gravity load vector
    for (int i = 0; i < nelement; i++) {
        Qg.middleRows(6*i, 12) += Qg_local;
    }

    // Set beam gravity force
    b->set_gravity_force(Qg);
}

void force_engine::elastic_force(beam * b) {
    // Beam parameters
    int ndof = b->get_ndof();
    int nelement = b->get_nelement();
    double lelement = b->get_length() / ((double) nelement);
    double E = b->get_E();
    double area = b->get_area();
    double inertia = b->get_inertia();

    // Get gauss points of 3 and 5
    vector<double> x3(3, 0.0);
    vector<double> weight3(3, 0.0);
    utils::gauss_points(3, x3, weight3);

    vector<double> x5(5, 0.0);
    vector<double> weight5(5, 0.0);
    utils::gauss_points(5, x5, weight5);

    // Initialize elastic load vector
    VectorXd Qe = VectorXd::Zero(ndof);

    // Initialize local elastic forces
    Vector<double, 12> Qn;
    Vector<double, 12> Qf;

    // Current beam displacement
    VectorXd disp_cur = b->get_position();

    // Loop over all beam elements to compute elastic force
    for (int i = 0; i < nelement; i++) {
        // Calculate elastic force due to axial deformation
        Qn = VectorXd::Zero(12);
        for (int k = 0; k < 5; k++) {
            double x = (1.0 + x5[k]) * lelement / 2.0;
            Qn += weight5[k] * axial_force_integrand(disp_cur.middleRows(6*i, 12), x, lelement);
        }
        Qn *= (E * area * lelement / 2.0);

        // Calculate elastic force due to flexural deformation
        Qf = VectorXd::Zero(12);
        for (int k = 0; k < 3; k++) {
            double x = (1.0 + x3[k]) * lelement / 2.0;
            Qf += weight3[k] * flexural_force_integrand(disp_cur.middleRows(6*i, 12), x, lelement);
        }
        Qf *= (E * inertia * lelement / 2.0);

        // Add the two sources of elastic force
        Qe.middleRows(6*i, 12) += (Qn + Qf);
    }

    // Set elastic force
    b->set_elastic_force(Qe);
}

Vector<double, 12> force_engine::axial_force_integrand(Vector<double, 12> disp, double x, double lelement) {
    // First derivative of shape function with respect to x
    Vector<double, 4> Sx;
    utils::shape_fun(x, lelement, 1, Sx);

    // First derivative of position vector r with respect to x
    Vector3d rx;
    Matrix<double, 3, 4> disp_rs;
    disp_rs = disp.reshaped(3,4);
    rx = disp_rs * Sx;

    // Axial strain of element
    double epsilon;
    epsilon = (rx.dot(rx) - 1.0) / 2.0;

    // Partial derivative of axial strain to base vector e
    Vector<double, 12> epsilon_e;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 3; j++) {
            epsilon_e(3*i + j) = Sx(i) * rx(j);
        }
    }

    // Integrand
    Vector<double, 12> en_int;
    en_int = epsilon * epsilon_e;
    return en_int;
}

Vector<double, 12> force_engine::flexural_force_integrand(Vector<double, 12> disp, double x, double lelement) {
    // First and second derivative of shape function with respect to x
    Vector<double, 4> Sx, Sxx;
    utils::shape_fun(x, lelement, 1, Sx);
    utils::shape_fun(x, lelement, 2, Sxx);

    // Shape function matrix
    Matrix<double, 3, 12> Sx_mat, Sxx_mat;
    Sx_mat = MatrixXd::Zero(3, 12);
    Sxx_mat = MatrixXd::Zero(3, 12);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 3; j++) {
            Sx_mat(j, 3*i+j) = Sx(i);
            Sxx_mat(j, 3*i+j) = Sxx(i);
        }
    }

    // First and second derivative of position vector r with respect to x
    Matrix<double, 3, 4> disp_rs;
    disp_rs = disp.reshaped(3, 4);
    Vector3d rx, rxx;
    rx = disp_rs * Sx;
    rxx = disp_rs * Sxx;
    double rx_mag2 = rx.dot(rx);

    // Calculate v = rx x rxx
    Vector3d v;
    v = rx.cross(rxx);
    double v_mag2 = v.dot(v);

    // Evaluate the kernel of integral
    Vector<double, 12> ef_int;
    Matrix<double, 3, 12> temp;
    temp = utils::tilde(rx) * Sxx_mat - utils::tilde(rxx) * Sx_mat;
    ef_int = temp.transpose() * v / (rx_mag2 * rx_mag2 * rx_mag2)
            - 3.0 * v_mag2 / (rx_mag2 * rx_mag2 * rx_mag2 * rx_mag2) * (Sx_mat.transpose() * rx);
    return ef_int;
}

void force_engine::add_constraint_load(int botCnstr, int topCnstr, Eigen::VectorXd & Q) {
    // Update force based on constraint
    switch (botCnstr) {
        case 1:
            Q.topRows(3) = Eigen::Vector3d::Zero();
            break;
        case 2:
            Q.topRows(6) = Eigen::Vector<double, 6>::Zero();
            break;
        default:
            break;
    }

    switch (topCnstr) {
        case 1:
            Q.bottomRows(3) = Eigen::Vector3d::Zero();
            break;
        case 2:
            Q.bottomRows(6) = Eigen::Vector<double, 6>::Zero();
            break;
        default:
            break;
    }
}

void force_engine::distributed_load(beam * b, fluid_field * ff, double t) {
    // Gauss integration points
    int np = 5;
    vector<double> x5(np, 0.0);
    vector<double> weight5(np, 0.0);
    utils::gauss_points(np, x5, weight5);

    // Beam properties
    double l_element = b->get_length() / ((double) b->get_nelement());
    double width = b->get_width();

    // Fluid field properties
    double rhof = ff->get_rhof();
    double Cd = ff->get_Cd();
    double Cm = ff->get_Cm();

    // Shape function matrix
    vector<Matrix<double, 3, 12>> S_mat, Sx_mat;
    for (int i = 0; i < np; i++) {
        double x = (1.0 + x5[i]) / 2.0;
        Vector4d S = Vector4d::Zero();
        utils::shape_fun(x*l_element, l_element, 0, S);
        Vector4d Sx = Vector4d::Zero();
        utils::shape_fun(x*l_element, l_element, 1, Sx);
        S_mat[i] = MatrixXd::Zero(3,12);
        Sx_mat[i] = MatrixXd::Zero(3, 12);
        for (int ii = 0; ii < 4; ii++) {
            for (int jj = 0; jj < 3; jj++) {
                S_mat[i](jj, 3*(ii-1)+jj) = S(ii);
                Sx_mat[i](jj, 3*(ii-1)+jj) = Sx(ii);
            }
        }
    }

    // Initialize external force vector
    VectorXd Qd = VectorXd::Zero(b->get_ndof());

    // Get displacement vector
    VectorXd disp = b->get_position();
    VectorXd vel = b->get_velocity();

    // Loop over all beam elements
    for (int ie = 0; ie < b->get_nelement(); ie++) {
        int istart = 6 * ie;
        // Line force of element
        VectorXd Qd_ele = VectorXd::Zero(12);
        // Loop over all element points
        for (int i = 0; i < np; i++) {
            // Get position and velocity vectors at current point
            Vector3d r_c = (S_mat[i])*(disp.segment(istart,12));
            Vector3d v_c = (S_mat[i])*(vel.segment(istart,12));
            Vector3d rx_c = (Sx_mat[i])*(disp.segment(istart,12));
            // Fluid velocity
            Vector3d v_f = ff->get_velocity(r_c, t);
            Vector3d a_f = ff->get_acceleration(r_c, t);
            // Relative velocity
            Vector3d v_rel = v_f - v_c;
            Vector3d v_reln = v_rel - v_rel.dot(rx_c) / rx_c.squaredNorm() * rx_c;
            Vector3d fl_c = 0.5 * rhof * Cd * width * v_reln.norm() * v_reln;
            // Added mass
            Vector3d a_fn = a_f - a_f.dot(rx_c) / rx_c.squaredNorm() * rx_c;
            Vector3d f_am = rhof * Cm * M_PI * width * width / 4.0 * a_fn;
            Vector3d f_tot = fl_c + f_am;
            if (t <= 1.0) {
                f_tot *= t;
            }
            // Compute line force of element using Gauss integration
            Qd_ele += (weight5[i] * (S_mat[i]).transpose() * f_tot);
        }
        // Rescale distributed force for each element
        Qd_ele *= (l_element / 2.0);
        // Construct the global force vector
        Qd.segment(istart, 12) += Qd_ele;
    }

    // Set distributed force
    b->set_dist_force(Qd);
}