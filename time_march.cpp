#include "time_march.h"
#include "system_engine.h"

void time_march::range_kutta(beam * b, double t, double h, system_engine * s_engine) {
    // Computer the status at 'xi_1' step
    VectorXd pos_xi_1 = h * b->get_velocity();
    VectorXd vel_xi_1 = h * b->get_acceleration();
    VectorXd q_zero = VectorXd::Zero(b->get_ndof());
    beam b_xi_1(*b);
    b_xi_1.update_position(b->get_position() + 0.5*pos_xi_1);
    b_xi_1.update_velocity(b->get_velocity() + 0.5*vel_xi_1);
    b_xi_1.set_gravity_force(b->get_gravity_force());
    b_xi_1.set_point_force(b->get_point_force());
    b_xi_1.set_external_force(b->get_external_force());
    s_engine->get_force_engine()->elastic_force(&b_xi_1);
    if (s_engine->get_distributed_force_switch()) {
        s_engine->get_force_engine()->distributed_load(&b_xi_1, s_engine->get_fluid_field(), t+0.5*h);
    } else {
        b_xi_1.set_dist_force(q_zero);
    }
    if (s_engine->get_damping_force_switch()) {
        s_engine->get_force_engine()->damping_load(&b_xi_1);
    } else {
        b_xi_1.set_damping_force(q_zero);
    }
    b_xi_1.set_total_force();
    s_engine->get_force_engine()->add_constraint_load(&b_xi_1);
    VectorXd acc_xi_1 = s_engine->solving_linear_system_of_beam(&b_xi_1);

    // Compute the status at 'xi_2' step
    VectorXd pos_xi_2 = h * (b->get_velocity() + 0.5 * vel_xi_1);
    VectorXd vel_xi_2 = h * acc_xi_1;
    beam b_xi_2(*b);
    b_xi_2.update_position(b->get_position() + 0.5*pos_xi_2);
    b_xi_2.update_velocity(b->get_velocity() + 0.5*vel_xi_2);
    b_xi_2.set_gravity_force(b->get_gravity_force());
    b_xi_2.set_point_force(b->get_point_force());
    b_xi_2.set_external_force(b->get_external_force());
    s_engine->get_force_engine()->elastic_force(&b_xi_2);
    if (s_engine->get_distributed_force_switch()) {
        s_engine->get_force_engine()->distributed_load(&b_xi_2, s_engine->get_fluid_field(), t+0.5*h);
    } else {
        b_xi_2.set_dist_force(q_zero);
    }
    if (s_engine->get_damping_force_switch()) {
        s_engine->get_force_engine()->damping_load(&b_xi_2);
    } else {
        b_xi_2.set_damping_force(q_zero);
    }
    b_xi_2.set_total_force();
    s_engine->get_force_engine()->add_constraint_load(&b_xi_2);
    VectorXd acc_xi_2 = s_engine->solving_linear_system_of_beam(&b_xi_2);

    // Compute the status at 'xi_3' step
    VectorXd pos_xi_3 = h * (b->get_velocity() + 0.5 * vel_xi_2);
    VectorXd vel_xi_3 = h * acc_xi_2;
    beam b_xi_3(*b);
    b_xi_3.update_position(b->get_position() + pos_xi_3);
    b_xi_3.update_velocity(b->get_velocity() + vel_xi_3);
    b_xi_3.set_gravity_force(b->get_gravity_force());
    b_xi_3.set_point_force(b->get_point_force());
    b_xi_3.set_external_force(b->get_external_force());
    s_engine->get_force_engine()->elastic_force(&b_xi_3);
    if (s_engine->get_distributed_force_switch()) {
        s_engine->get_force_engine()->distributed_load(&b_xi_3, s_engine->get_fluid_field(), t+h);
    } else {
        b_xi_3.set_dist_force(q_zero);
    }
    if (s_engine->get_damping_force_switch()) {
        s_engine->get_force_engine()->damping_load(&b_xi_3);
    } else {
        b_xi_3.set_damping_force(q_zero);
    }
    b_xi_3.set_total_force();
    s_engine->get_force_engine()->add_constraint_load(&b_xi_3);
    VectorXd acc_xi_3 = s_engine->solving_linear_system_of_beam(&b_xi_3);

    // Compute the status at 'xi_4' step
    VectorXd pos_xi_4 = h * (b->get_velocity() + vel_xi_3);
    VectorXd vel_xi_4 = h * acc_xi_3;

    // Compute the status at next time step
    VectorXd pos = b->get_position() + (pos_xi_1 + 2.0*pos_xi_2 + 2.0*pos_xi_3 + pos_xi_4)/6.0;
    VectorXd vel = b->get_velocity() + (vel_xi_1 + 2.0*vel_xi_2 + 2.0*vel_xi_3 + vel_xi_4)/6.0;
    b->update_position(pos);
    b->update_velocity(vel);
    s_engine->get_force_engine()->elastic_force(b);
    if (s_engine->get_distributed_force_switch()) {
        s_engine->get_force_engine()->distributed_load(b, s_engine->get_fluid_field(), t+h);
    }
    if (s_engine->get_damping_force_switch()) {
        s_engine->get_force_engine()->damping_load(b);
    }
    b->set_total_force();
    s_engine->get_force_engine()->add_constraint_load(b);
    VectorXd acc = s_engine->solving_linear_system_of_beam(b);
    b->update_acceleration(acc);
}