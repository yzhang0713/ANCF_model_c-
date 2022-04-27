#include <iostream>
#include "utils.h"
#include "system_engine.h"
#include <Eigen/Dense>

using Eigen::MatrixXd;

int main() {
//    int n = 4;
//    vector<double> x(n, 0.0);
//    vector<double> weight(n, 0.0);
//    utils::gauss_points(n, x, weight);
//    for (int i = 0; i < n; i++) {
//        cout << i << ": " << x[i] << ", " << weight[i] << endl;
//    }
//    MatrixXd m(2,2);
//    m(0,0) = 3;
//    m(1,0) = 2.5;
//    m(0,1) = -1;
//    m(1,1) = m(1,0) + m(0,1);
//    std::cout << m << std::endl;
//    Eigen::VectorXd v1(5);
//    v1(0) = 1.0;
//    v1(1) = 2.0;
//    v1(2) = 3.0;
//    v1(3) = 4.0;
//    v1(4) = 5.0;
//    Eigen::Vector<double, 3> v2(2.0, 3.0, 4.0);
//    v1.middleRows(1, 3) += v2;
//    std::cout << v1 << std::endl;

    int debug = 0;
    int show_time = 1;

    // Initialize system engine
    system_engine s{};

    // Remove output directory
    s.remove_outputs();

    // Read beams and particles from files
    s.read_beams();
    if (debug)
        cout << "Read beam done" << endl;
    s.test_read_beams();
    s.read_particles();

    // Test files reading
    s.test_read_particles();

    // Read fluid field
    s.read_fluid_field();

    // Test fluid field reading
    s.test_read_fluid_field();
    if (debug)
        cout << "test fluid field done" << endl;

    // Read external load field
    s.read_external_load_field();

    // Load force switches
    s.load_switches();

    // Initialize all forces
    if (s.get_gravity_force_switch()) {
        s.update_beam_gravity_forces();
    } else {
        s.beam_gravity_force_off();
    }

    if (debug)
        cout << "update gravity force done" << endl;

    s.update_beam_elastic_forces();
    if (debug)
        cout << "update elastic force done" << endl;
    if (s.get_distributed_force_switch()) {
        s.update_beam_dist_forces();
    } else {
        s.beam_dist_force_off();
    }
    if (debug)
        cout << "update distributed force done" << endl;
//    s.reset_beam_point_forces();
    s.initialize_beam_point_forces();
//    s.debug_beam_point_forces();
    if (debug)
        cout << "reset point force done" << endl;
    if (s.get_point_force_switch()) {
        s.update_beam_point_forces();
    } else {
        s.beam_point_force_off();
    }
    if (debug)
        cout << "update point force done" << endl;
    if (s.get_external_force_switch()) {
        s.update_beam_external_forces();
    } else {
        s.beam_external_force_off();
    }
    if (debug)
        cout << "update external force done" << endl;
    if (s.get_damping_force_switch()) {
        s.update_beam_damping_forces();
    } else {
        s.beam_damping_force_off();
    }
    if (debug)
        cout << "update damping force done" << endl;
    s.update_beam_total_forces();
    if (debug)
        cout << "update total force done" << endl;
    s.update_beam_forces_with_constraint();
    if (debug)
        cout << "update constraint done" << endl;

    if (debug)
        cout << "initialize beam forces done" << endl;

    // Start time marching
    s.initialize_time();
    s.set_time_step(0.00001);
    double t_end = 20.0;
    int update_count = 0;
    while (s.get_cur_time() < t_end) {
        // Store current beam information
//        s.store_beam_information();
        if (update_count == 100) {
            s.store_beam_information();
            update_count = 0;
        }
        // Moving beam to next position
        s.update_beams();
        // Update time
        s.update_time();
        // Update point load based on new position
        s.reset_beam_point_forces();
//        s.debug_beam_point_forces();
        if (s.get_point_force_switch()) {
            s.update_beam_point_forces();
        }
        if (s.get_external_force_switch()) {
            s.update_beam_external_forces();
        }
        if (s.get_damping_force_switch()) {
            s.update_beam_damping_forces();
        }
        if (s.check_write()) {
            if (debug)
                cout << "write to file at time " << s.get_cur_time() << endl;
            s.write_to_file();
        }
        if (show_time)
            cout << "current time " << s.get_cur_time() << endl;
        update_count++;
    }
    s.write_to_file();

    return 0;
}
