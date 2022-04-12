#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
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

void force_engine::add_constraint_load(beam * b) {
    // Update force based on constraint
    switch (b->get_botCnstr()) {
        case 1:
            b->get_total_force().topRows(3) = Eigen::Vector3d::Zero();
            break;
        case 2:
            b->get_total_force().topRows(6) = Eigen::Vector<double, 6>::Zero();
            break;
        default:
            break;
    }

    switch (b->get_topCnstr()) {
        case 1:
            b->get_total_force().bottomRows(3) = Eigen::Vector3d::Zero();
            break;
        case 2:
            b->get_total_force().bottomRows(6) = Eigen::Vector<double, 6>::Zero();
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

//    cout << "gather information done" << endl;

    // Shape function matrix
    vector<Matrix<double, 3, 12>> S_mat, Sx_mat;
//    cout << "S_mat size " << S_mat.size() << endl;
//    cout << "Sx_mat size " << Sx_mat.size() << endl;
    for (int i = 0; i < np; i++) {
        double x = (1.0 + x5[i]) / 2.0;
        Vector4d S = Vector4d::Zero();
        utils::shape_fun(x*l_element, l_element, 0, S);
        Vector4d Sx = Vector4d::Zero();
        utils::shape_fun(x*l_element, l_element, 1, Sx);
//        cout << "shape fun done" << endl;
        MatrixXd S_mat_local = MatrixXd::Zero(3, 12);
        MatrixXd Sx_mat_local = MatrixXd::Zero(3, 12);
        for (int ii = 0; ii < 4; ii++) {
            for (int jj = 0; jj < 3; jj++) {
                S_mat_local(jj, 3*ii+jj) = S(ii);
                Sx_mat_local(jj, 3*ii+jj) = Sx(ii);
            }
        }
//        cout << "shape matrix done" << endl;
        S_mat.push_back(S_mat_local);
        Sx_mat.push_back(Sx_mat_local);
    }

//    cout << "shape function matrix done" << endl;

    // Initialize external force vector
    VectorXd Qd = VectorXd::Zero(b->get_ndof());

    // Get displacement vector
    VectorXd disp = b->get_position();
    VectorXd vel = b->get_velocity();

    // Loop over all beam elements
    for (int ie = 0; ie < b->get_nelement(); ie++) {
        int istart = 6 * ie;
//        cout << "istart: " << istart << endl;
        // Line force of element
        VectorXd Qd_ele = VectorXd::Zero(12);
        // Loop over all element points
        for (int i = 0; i < np; i++) {
            // Get position and velocity vectors at current point
            Vector3d r_c = (S_mat[i])*(disp.segment(istart,12));
            Vector3d v_c = (S_mat[i])*(vel.segment(istart,12));
            Vector3d rx_c = (Sx_mat[i])*(disp.segment(istart,12));
//            cout << "position and velocity vector" << endl;
            // Fluid velocity
            Vector3d v_f = ff->get_velocity(r_c, t);
            Vector3d a_f = ff->get_acceleration(r_c, t);
//            cout << "fluid velocity and acceleration" << endl;
//            cout << v_f << endl;
//            cout << a_f << endl;
            // Relative velocity
            Vector3d v_rel = v_f - v_c;
            Vector3d v_reln = v_rel - v_rel.dot(rx_c) / rx_c.squaredNorm() * rx_c;
            Vector3d fl_c = 0.5 * rhof * Cd * width * v_reln.norm() * v_reln;
//            cout << "relative velocity" << endl;
//            cout << fl_c << endl;
            // Added mass
            Vector3d a_fn = a_f - a_f.dot(rx_c) / rx_c.squaredNorm() * rx_c;
            Vector3d f_am = rhof * Cm * M_PI * width * width / 4.0 * a_fn;
            Vector3d f_tot = fl_c + f_am;
            if (t <= 1.0) {
                f_tot *= t;
            }
            // Compute line force of element using Gauss integration
            Qd_ele += (weight5[i] * (S_mat[i]).transpose() * f_tot);
//            cout << "adding to element done" << endl;
        }
        // Rescale distributed force for each element
        Qd_ele *= (l_element / 2.0);
        // Construct the global force vector
        Qd.segment(istart, 12) += Qd_ele;
    }

    // Set distributed force
    b->set_dist_force(Qd);
}

void force_engine::point_load_element_level(beam * b1, int ielement_1, beam * b2, int ielement_2) {
    // Initialize point load
    Vector<double, 12> Qp1, Qp2;
    Qp1.setZero();
    Qp2.setZero();

    // Contact parameters
    double stiffness_factor = 1.0e-3;
    double E_star = stiffness_factor / ((2.0*(1.0-pow(b1->get_nu(),2))/b1->get_E()) +
            (2.0*(1.0-pow(b2->get_nu(),2)/b2->get_E())));
    double R = 1.0 / (2.0/b1->get_thick() + 2.0/b2->get_thick());
    double K = 4.0 / 3.0 * E_star * sqrt(R);

    // The particle level contact resolution will be done in three levels
    // Using bounding sphere to do level 1 and level 2
    // Using particle contact to do level 3
    double l_element_1 = b1->get_length() / ((double) b1->get_nelement());
    int n_particle_1 = (int) (l_element_1 / b1->get_thick());
    Vector3i n_centers_1 = utils::get_level_3_centers(n_particle_1);
    double l_element_2 = b2->get_length() / ((double) b2->get_nelement());
    int n_particle_2 = (int) (l_element_2 / b2->get_thick());
    Vector3i n_centers_2 = utils::get_level_3_centers(n_particle_2);

    // First level of overlap check
    vector<double> level_1_pos_1;
    for (int i = 0; i < n_centers_1(0); i++) {
        level_1_pos_1.push_back(((double)i+0.5)/((double)n_centers_1(0)));
    }
    double level_1_rad_1 = l_element_1 / ((double) n_centers_1(0));
    vector<bounding_sphere> level_1_bss_1 = utils::get_bs_from_beam_element(b1, ielement_1, level_1_pos_1, level_1_rad_1);
    vector<double> level_1_pos_2;
    for (int i = 0; i < n_centers_2(0); i++) {
        level_1_pos_2.push_back(((double)i+0.5)/((double)n_centers_2(0)));
    }
    double level_1_rad_2 = l_element_2 / ((double) n_centers_2(0));
    vector<bounding_sphere> level_1_bss_2 = utils::get_bs_from_beam_element(b2, ielement_2, level_1_pos_2, level_1_rad_2);

    // Overlap record
    map<int, vector<int>> level_1_overlap;
    for (int i = 0; i < n_centers_1(0); i++) {
        bounding_sphere bs_1 = level_1_bss_1[i];
        for (int j = 0; j < n_centers_2(0); j++) {
            bounding_sphere bs_2 = level_1_bss_2[j];
            if (!bs_1.is_bs_disjoint(bs_2)) {
                level_1_overlap[i].push_back(j);
            }
        }
    }

    // Second level of overlap check
    map<int, vector<int>> level_2_overlap;
    for (auto const& x : level_1_overlap) {
        int n1 = x.first;
        vector<int> n2_list = x.second;
        vector<double> level_2_pos_1;
        int add_1 = n1 * n_centers_1(1);
        int base_1 = n_centers_1(1) * n_centers_1(0);
        for (int i = 0; i < n_centers_1(1); i++) {
            level_2_pos_1.push_back(((double)i+0.5+(double)add_1)/((double) base_1));
        }
        double level_2_rad_1 = l_element_1 / ((double) base_1);
        vector<bounding_sphere> level_2_bss_1 = utils::get_bs_from_beam_element(b1, ielement_1, level_2_pos_1, level_2_rad_1);
        for (auto n2 : n2_list) {
            vector<double> level_2_pos_2;
            int add_2 = n2 * n_centers_2(1);
            int base_2 = n_centers_2(1) * n_centers_2(0);
            for (int i = 0; i < n_centers_2(1); i++) {
                level_2_pos_2.push_back(((double)i+0.5+(double)add_2)/((double) base_2));
            }
            double level_2_rad_2 = l_element_2 / ((double) base_2);
            vector<bounding_sphere> level_2_bss_2 = utils::get_bs_from_beam_element(b2, ielement_2, level_2_pos_2, level_2_rad_2);

            // Check overlap
            for (int i = 0; i < level_2_bss_1.size(); i++) {
                bounding_sphere bs_1 = level_2_bss_1[i];
                for (int j = 0; j < level_2_bss_2.size(); j++) {
                    bounding_sphere bs_2 = level_2_bss_2[j];
                    if (!bs_1.is_bs_disjoint(bs_2)) {
                        level_2_overlap[i+add_1].push_back(j+add_2);
                    }
                }
            }
        }
    }

    // Third level (particle level) of overlap check
    for (auto const& x  : level_2_overlap) {
        int n1 = x.first;
        vector<int> n2_list = x.second;
        vector<double> level_3_pos_1;
        int add_1 = n1 * n_centers_1(2);
        int base_1 = n_centers_1.prod();
        for (int i = 0; i < n_centers_1(2); i++) {
            level_3_pos_1.push_back(((double)i+0.5+(double)add_1)/((double) base_1));
        }
        double level_3_rad_1 = b1->get_thick()/2.0;
        vector<bounding_sphere> level_3_bss_1 = utils::get_bs_from_beam_element(b1, ielement_1, level_3_pos_1, level_3_rad_1);
        for (auto n2 : n2_list) {
            vector<double> level_3_pos_2;
            int add_2 = n2 * n_centers_2(2);
            int base_2 = n_centers_2.prod();
            for (int i = 0; i < n_centers_2(2); i++) {
                level_3_pos_2.push_back(((double)i+0.5+(double)add_2)/((double) base_2));
            }
            double level_3_rad_2 = b2->get_thick()/2.0;
            vector<bounding_sphere> level_3_bss_2 = utils::get_bs_from_beam_element(b2, ielement_2, level_3_pos_2, level_3_rad_2);

            // Check overlap and compute contact force
            for (int i = 0; i < level_3_bss_1.size(); i++) {
                bounding_sphere bs_1 = level_3_bss_1[i];
                for (int j = 0; j < level_3_bss_2.size(); j++) {
                    bounding_sphere bs_2 = level_3_bss_2[j];
                    if (!bs_1.is_bs_disjoint(bs_2)) {
                        double distance = bs_1.get_distance(bs_2);
                        double f_mag = K * pow(b1->get_thick()/2.0+b2->get_thick()/2.0-distance,1.5);
                        Vector3d direction = get_direction(bs_1, bs_2);
                        Matrix<double, 3, 12> S_mat = utils::get_shape_matrix(level_3_pos_1[i]*l_element_1, l_element_1, 0);
                        Qp1 += (f_mag / distance * (S_mat.transpose())*direction);
                        S_mat = utils::get_shape_matrix(level_3_pos_2[j]*l_element_2, l_element_2, 0);
                        Qp2 += (-f_mag / distance * (S_mat.transpose())*direction);
                    }
                }
            }
        }
    }

    // Update the point load of beam
    b1->get_point_force().segment(6*(ielement_1-1),12) += Qp1;
    b2->get_point_force().segment(6*(ielement_2-1),12) += Qp2;
}

void force_engine::point_load(beam * b1, beam * b2) {
    // Firstly, do beam level contact detection

    // Start with bounding sphere check to save computation
    bounding_sphere bs_1 {};
    bs_1.set_center_point(utils::get_points_from_beam(b1, 1, b1->get_nelement(), 3)[1]);
    bs_1.set_radius(b1->get_length()/2.0);
    bounding_sphere bs_2 {};
    bs_2.set_center_point(utils::get_points_from_beam(b2, 1, b2->get_nelement(), 3)[1]);
    bs_2.set_radius(b2->get_length()/2.0);
    if (bs_1.is_bs_disjoint(bs_2)) {
        return;
    }

    // Then check with oriented bounding box
    oriented_bounding_box obb_1 {};
    obb_1.set_obb_info(utils::get_points_from_beam(b1, 1, b1->get_nelement(), 5));
    oriented_bounding_box obb_2 {};
    obb_2.set_obb_info(utils::get_points_from_beam(b2, 1, b2->get_nelement(), 5));
    // In case of no disjoint, no contact force
    if (obb_1.is_obb_disjoint(obb_2, b1->get_thick()/2.0+b2->get_thick()/2.0)) {
        return;
    }

    // Otherwise, call segment collision detection recursively
    point_load_segment_level(b1, 1, b1->get_nelement(), b2, 1, b2->get_nelement());
}

void force_engine::point_load_segment_level(beam * b1, int start_element_1, int n_element_1,
                                            beam * b2, int start_element_2, int n_element_2) {
    // Base case, do element level check
    if (n_element_1 == 1 && n_element_2 == 1) {
        point_load_element_level(b1, start_element_1, b2, start_element_2);
        return;
    }

    double gap = b1->get_thick()/2.0 + b2->get_thick()/2.0;

    // In case n_element_1 = 1, only need to split n_element_2
    if (n_element_1 == 1) {
        // For beam 1
        oriented_bounding_box obb_1 {};
        obb_1.set_obb_info(utils::get_points_from_beam(b1, start_element_1, 1, 5));
        // For beam 2
        int n_element_21 = n_element_2/2;
        int start_element_22 = start_element_2 + n_element_21;
        oriented_bounding_box obb_21 {};
        obb_21.set_obb_info(utils::get_points_from_beam(b2, start_element_2, n_element_21, 5));
        oriented_bounding_box obb_22 {};
        obb_22.set_obb_info(utils::get_points_from_beam(b2, start_element_22, n_element_2-n_element_21, 5));
        // Check obb interaction
        if (!obb_1.is_obb_disjoint(obb_21, gap)) {
            point_load_segment_level(b1, start_element_1, n_element_1, b2, start_element_2, n_element_21);
        }
        if (!obb_1.is_obb_disjoint(obb_22, gap)) {
            point_load_segment_level(b1, start_element_1, n_element_1, b2, start_element_22, n_element_2-n_element_21);
        }
        return;
    }

    // In case n_element_2 = 1, only need to split n_element_1
    if (n_element_2 == 1) {
        // For beam 1
        int n_element_11 = n_element_1/2;
        int start_element_12 = start_element_1 + n_element_11;
        // Get obbs for each piece
        oriented_bounding_box obb_11 {};
        obb_11.set_obb_info(utils::get_points_from_beam(b1, start_element_1, n_element_11, 5));
        oriented_bounding_box obb_12 {};
        obb_12.set_obb_info(utils::get_points_from_beam(b1, start_element_12, n_element_1-n_element_11, 5));
        // For beam 2
        oriented_bounding_box obb_2 {};
        obb_2.set_obb_info(utils::get_points_from_beam(b2, start_element_2, 1, 5));
        // Check obb interaction
        if (!obb_2.is_obb_disjoint(obb_11, gap)) {
            point_load_segment_level(b1, start_element_1, n_element_11, b2, start_element_2, n_element_2);
        }
        if (!obb_2.is_obb_disjoint(obb_12, gap)) {
            point_load_segment_level(b1, start_element_12, n_element_1-n_element_11, b2, start_element_2, n_element_2);
        }
        return;
    }


    // Otherwise, split both segments into two pieces, and check collision between pieces
    // Beam 1
    int n_element_11 = n_element_1/2;
    int start_element_12 = start_element_1 + n_element_11;
    // Get obbs for each piece
    oriented_bounding_box obb_11 {};
    obb_11.set_obb_info(utils::get_points_from_beam(b1, start_element_1, n_element_11, 5));
    oriented_bounding_box obb_12 {};
    obb_12.set_obb_info(utils::get_points_from_beam(b1, start_element_12, n_element_1-n_element_11, 5));
    // Beam 2
    int n_element_21 = n_element_2/2;
    int start_element_22 = start_element_2 + n_element_21;
    oriented_bounding_box obb_21 {};
    obb_21.set_obb_info(utils::get_points_from_beam(b2, start_element_2, n_element_21, 5));
    oriented_bounding_box obb_22 {};
    obb_22.set_obb_info(utils::get_points_from_beam(b2, start_element_22, n_element_2-n_element_21, 5));
    // Check obb interaction
    if (!obb_11.is_obb_disjoint(obb_21, gap)) {
        point_load_segment_level(b1, start_element_1, n_element_11, b2, start_element_2, n_element_21);
    }
    if (!obb_11.is_obb_disjoint(obb_22, gap)) {
        point_load_segment_level(b1, start_element_1, n_element_11, b2, start_element_22, n_element_2-n_element_21);
    }
    if (!obb_12.is_obb_disjoint(obb_21, gap)) {
        point_load_segment_level(b1, start_element_12, n_element_1-n_element_11, b2, start_element_2, n_element_21);
    }
    if (!obb_12.is_obb_disjoint(obb_22, gap)) {
        point_load_segment_level(b1, start_element_12, n_element_1-n_element_11, b2, start_element_22, n_element_2-n_element_21);
    }
}

void force_engine::external_load(beam * b, external_load_field * el_field) {
    b->set_external_force(VectorXd::Zero(b->get_ndof()));
    cout << "initialize external force done" << endl;
    double l_element = b->get_length() / ((double) b->get_nelement());
    for (external_load_point el_point : el_field->get_forces()) {
        double pos = el_point.get_position();
        Vector3d force = el_point.get_force();
        int i_element = (int) (pos * b->get_nelement());
        cout << "i element: " << i_element << endl;
        if (i_element == b->get_nelement()) {
            // Meaning force at tip of beam
            Matrix<double, 3, 12> S_mat = utils::get_shape_matrix(l_element, l_element, 0);
            cout << "shape matrix " << S_mat << endl;
            b->get_external_force().tail(12) += ((S_mat.transpose())*force);
            cout << "assign force done" << endl;
        } else {
            // Internal points
            Matrix<double, 3, 12> S_mat = utils::get_shape_matrix((pos*b->get_nelement()-i_element)*l_element, l_element, 0);
            b->get_external_force().segment(6*i_element,12) += ((S_mat.transpose())*force);
        }
    }
}

void force_engine::damping_load(beam * b) {
    b->set_damping_force(-0.10 * b->get_velocity());
}