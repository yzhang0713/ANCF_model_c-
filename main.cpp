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


    // Initialize system engine
    system_engine s{};

    // Read beams and particles from files
    s.read_beams();
    s.read_particles();

    // Test files reading
    s.test_read_beams();
    s.test_read_particles();

    // Read fluid field
    s.read_fluid_field();

    // Test fluid field reading
    s.test_read_fluid_field();

    // Initialize all forces
    s.update_beam_gravity_forces();
    s.update_beam_elastic_forces();
    s.update_beam_dist_forces();
    s.reset_beam_point_forces();
    s.update_beam_point_forces();
    s.update_beam_total_forces();
    s.update_beam_forces_with_constraint();

    // Start time marching
    double t_end = 100.0;
    while (s.get_cur_time() < t_end) {
        // Store current beam information
        s.store_beam_information();
        // Moving beam to next position
        s.update_beams();
        // Update time
        s.update_time();
        // Update point load based on new position
        s.reset_beam_point_forces();
        s.update_beam_point_forces();
    }

    return 0;
}
