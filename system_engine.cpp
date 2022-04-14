#include "system_engine.h"
#include <filesystem>

namespace fs = filesystem;

system_engine::system_engine() {
    f_engine = new force_engine();
    el_field = new external_load_field();
    t = 0.0;
}

void system_engine::initialize_time() {
    t = 0.0;
}

void system_engine::set_current_time(double time) {
    t = time;
}

void system_engine::set_time_step(double time_step) {
    h = time_step;
}

void system_engine::read_beams() {
    ifstream beam_file;
    beam_file.open("./inputs/input_beam");
    ifstream beam_init_position;
    beam_init_position.open("./inputs/input_beam_position");
    ifstream beam_init_velocity;
    beam_init_velocity.open("./inputs/input_beam_velocity");
    int def_vel = 0;
    if (beam_init_velocity) {
        def_vel = 1;
    }
    ifstream beam_init_acceleration;
    beam_init_acceleration.open("./inputs/input_beam_acceleration");
    int def_acc = 0;
    if (beam_init_acceleration) {
        def_acc = 1;
    }
    string beam_line;
    string beam_pos;
    string beam_vel;
    string beam_acc;
    while (getline(beam_file, beam_line)) {
        if (debug) {
            cout << "value in beam file: " << beam_line << endl;
        }
        istringstream iss(beam_line);
        double E, nu, rho, length, thick, width;
        int nelement, botCnstr, topCnstr;
        if (!(iss >> E >> nu >> rho >> length >> thick >> width >> nelement >> botCnstr >> topCnstr)) {
            cout << "Error reading beam file." << endl;
            break;
        }
        beam_builder * b_builder = new beam_builder();
        b_builder->set_E(E);
        b_builder->set_nu(nu);
        b_builder->set_rho(rho);
        b_builder->set_length(length);
        b_builder->set_thick(thick);
        b_builder->set_width(width);
        b_builder->set_nelement(nelement);
        b_builder->set_botCnstr(botCnstr);
        b_builder->set_topCnstr(topCnstr);
        int ndof = 6 * (nelement + 1);
        getline(beam_init_position, beam_pos);
        istringstream iss_pos(beam_pos);
        VectorXd pos(ndof);
        for (int i = 0; i < ndof; i++) {
            double cur_pos;
            iss_pos >> cur_pos;
            pos(i) = cur_pos;
        }
        b_builder->set_position(pos);
        if (debug)
            cout << "set position done" << endl;
        if (!def_vel) {
            b_builder->set_velocity();
        } else {
            getline(beam_init_velocity, beam_vel);
            istringstream iss_vel(beam_vel);
            VectorXd vel(ndof);
            for (int i = 0; i < ndof; i++) {
                double cur_vel;
                iss_vel >> cur_vel;
                vel(i) = cur_vel;
            }
            b_builder->set_velocity(vel);
        }
        if (debug)
            cout << "set velocity done" << endl;
        if (!def_acc) {
            b_builder->set_acceleration();
        } else {
            getline(beam_init_acceleration, beam_acc);
            istringstream iss_acc(beam_acc);
            VectorXd acc(ndof);
            for (int i = 0; i < ndof; i++) {
                double cur_acc;
                iss_acc >> cur_acc;
                acc(i) = cur_acc;
            }
            b_builder->set_acceleration(acc);
        }
        if (debug)
            cout << "set acceleration done" << endl;
        beam * b = b_builder->get_beam();
        beams.push_back(b);
    }
    beam_file.close();
    beam_init_position.close();
    if (def_vel) {
        beam_init_velocity.close();
    }
    if (def_acc) {
        beam_init_acceleration.close();
    }
    if (debug) {
        cout << "number of beams: " << beams.size() << endl;
    }
}

void system_engine::test_read_beams() {
    for (beam * b : beams) {
        cout << *b;
    }
}

void system_engine::read_particles() {
    ifstream par_file;
    par_file.open("./inputs/input_particle");
    ifstream par_init_position;
    par_init_position.open("./inputs/input_particle_position");
    ifstream par_init_velocity;
    par_init_velocity.open("./inputs/input_particle_velocity");
    int def_vel = 0;
    if (par_init_velocity) {
        def_vel = 1;
    }
    ifstream par_init_acceleration;
    par_init_acceleration.open("./inputs/input_particle_acceleration");
    int def_acc = 0;
    if (par_init_acceleration) {
        def_acc = 1;
    }
    string par_line;
    string par_pos;
    string par_vel;
    string par_acc;
    while (getline(par_file, par_line)) {
        istringstream iss(par_line);
        double E, rho, r;
        if (!(iss >> E >> rho >> r)) {
            cout << "Error reading particle file." << endl;
            break;
        }
        particle_builder * p_builder = new particle_builder();
        p_builder->set_E(E);
        p_builder->set_rho(rho);
        p_builder->set_r(r);
        getline(par_init_position, par_pos);
        istringstream iss_pos(par_pos);
        Vector3d pos;
        double pos1, pos2, pos3;
        iss_pos >> pos1 >> pos2 >> pos3;
        pos << pos1, pos2, pos3;
        p_builder->set_position(pos);
        if (!def_vel) {
            p_builder->set_velocity();
        } else {
            getline(par_init_velocity, par_vel);
            istringstream iss_vel(par_vel);
            Vector3d vel;
            double vel1, vel2, vel3;
            iss_vel >> vel1 >> vel2 >> vel3;
            vel << vel1, vel2, vel3;
            p_builder->set_velocity(vel);
        }
        if (!def_acc) {
            p_builder->set_acceleration();
        } else {
            getline(par_init_acceleration, par_acc);
            istringstream iss_acc(par_acc);
            Vector3d acc;
            double acc1, acc2, acc3;
            iss_acc >> acc1 >> acc2 >> acc3;
            acc << acc1, acc2, acc3;
            p_builder->set_acceleration(acc);
        }
        particle * p = p_builder->get_particle();
        particles.push_back(p);
    }
    par_file.close();
    par_init_position.close();
    if (def_vel) {
        par_init_velocity.close();
    }
    if (def_acc) {
        par_init_acceleration.close();
    }
}

void system_engine::test_read_particles() {
    for (particle * p : particles) {
        cout << *p;
    }
}

void system_engine::read_external_load_field() {
    ifstream el_file;
    el_file.open("./inputs/input_external_load_field");
    string el_line;
    vector<external_load_point> el_points;
    while (getline(el_file, el_line)) {
        istringstream iss(el_line);
        double pos;
        double fx;
        double fy;
        double fz;
        if (!(iss >> pos >> fx >> fy >> fz)) {
            cout << "Error reading external load field file." << endl;
        }
        external_load_point el_point {};
        el_point.set_position(pos);
        Vector3d force;
        force << fx, fy, fz;
        el_point.set_force(force);
        el_points.push_back(el_point);
    }
    el_field->set_forces(el_points);
    el_file.close();
}

void system_engine::read_fluid_field() {
    ifstream ff_file;
    ff_file.open("./inputs/input_fluid_field");
    string ff_line;
    fluid_field_builder * ff_builder = new fluid_field_builder();
    while (getline(ff_file, ff_line)) {
        istringstream iss(ff_line);
        string field_name;
        double field_value;
        if (!(iss >> field_name >> field_value)) {
            cout << "Error reading fluid field file." << endl;
        }
        if (field_name == "Uv") {
            ff_builder->set_Uv(field_value);
        } else if (field_name == "Cv") {
            ff_builder->set_Cv(field_value);
        } else if (field_name == "kv") {
            ff_builder->set_kv(field_value);
        } else if (field_name == "zv") {
            ff_builder->set_zv(field_value);
        } else if (field_name == "rhof") {
            ff_builder->set_rhof(field_value);
        } else if (field_name == "Cm") {
            ff_builder->set_Cm(field_value);
        } else if (field_name == "Cd") {
            ff_builder->set_Cd(field_value);
        } else if (field_name == "tml") {
            ff_builder->set_tml(field_value);
        }
    }
    f_field = ff_builder->get_fluid_field();
    ff_file.close();
}

void system_engine::test_read_fluid_field() {
    cout << *f_field;
}

void system_engine::update_beam_gravity_forces() {
    for (beam * b : beams) {
        f_engine->gravity_force(b);
        if (debug) {
            cout << "gravity force" << b->get_gravity_force().transpose() << endl;
        }
    }
}

void system_engine::update_beam_elastic_forces() {
    for (beam * b : beams) {
        f_engine->elastic_force(b);
        if (debug) {
            cout << "elastic force" << b->get_elastic_force().transpose() << endl;
        }
    }
}

void system_engine::update_beam_dist_forces() {
    for (beam * b : beams) {
        f_engine->distributed_load(b, f_field, t);
        if (debug) {
            cout << "distributed force" << b->get_dist_force().transpose() << endl;
        }
    }
}

void system_engine::reset_beam_point_forces() {
    for (beam * b : beams) {
        b->set_point_force(VectorXd::Zero(b->get_ndof()));
        if (debug) {
            cout << "reset point forces" << b->get_point_force().transpose() << endl;
        }
    }
}

void system_engine::update_beam_point_forces() {
    for (int i = 0; i < beams.size()-1; i++) {
        for (int j = i+1; j < beams.size(); j++) {
            f_engine->point_load(beams[i], beams[j]);
        }
    }
}

void system_engine::update_beam_external_forces() {
    for (beam * b : beams) {
        f_engine->external_load(b, el_field);
        if (debug) {
            cout << "external force" << b->get_external_force().transpose() << endl;
        }
    }
}

void system_engine::update_beam_damping_forces() {
    for (beam * b : beams) {
        f_engine->damping_load(b);
        if (debug) {
            cout << "damping force" << b->get_damping_force().transpose() << endl;
        }
    }
}

void system_engine::update_beam_total_forces() {
    for (beam * b : beams) {
        b->set_total_force();
        if (debug) {
            cout << "total force" << b->get_total_force().transpose() << endl;
        }
    }
}

void system_engine::update_beam_forces_with_constraint() {
    for (beam * b : beams) {
        f_engine->add_constraint_load(b);
    }
}

VectorXd system_engine::solving_linear_system_of_beam(beam * b) {
    return (b->get_mass_LU().solve(b->get_total_force()));
}

void system_engine::store_beam_information() {
    for (beam * b : beams) {
//        cout << "record size before store: " << beam_pos_record.size() << endl;
        beam_pos_record[b->get_beam_id()].push_back(b->get_position());
        beam_vel_record[b->get_beam_id()].push_back(b->get_velocity());
        beam_acc_record[b->get_beam_id()].push_back(b->get_acceleration());
        VectorXd qe(b->get_elastic_force());
        beam_elastic_record[b->get_beam_id()].push_back(qe);
        VectorXd qg(b->get_gravity_force());
        beam_gravity_record[b->get_beam_id()].push_back(qg);
        VectorXd qd(b->get_damping_force());
        beam_damping_record[b->get_beam_id()].push_back(qd);
        VectorXd qex(b->get_external_force());
        beam_external_record[b->get_beam_id()].push_back(qex);
        VectorXd qdist(b->get_dist_force());
        beam_distribute_record[b->get_beam_id()].push_back(qdist);
        VectorXd qp(b->get_point_force());
        beam_point_record[b->get_beam_id()].push_back(qp);
        time_record.push_back(t);
//        cout << "record size after store: " << beam_pos_record.size() << endl;
    }
}

void system_engine::update_beams() {
    for (beam * b : beams) {
        time_march::range_kutta(b, t, h, this);
    }
}

void system_engine::write_to_file() {
//    cout << "write to file" << endl;
    check_output_folder();
    // Write outputs to each beam
    for (auto const& x : beam_pos_record) {
        int beam_id = x.first;
        string pos_name("./outputs/" + to_string(beam_id) + "/position");
        ofstream pos_file;
        pos_file.open(pos_name, ios_base::app);
        for (int i = 0; i < time_record.size(); i++) {
            pos_file << time_record[i] << " " << beam_pos_record[beam_id][i].transpose() << endl;
        }
        pos_file.close();
        string vel_name("./outputs/" + to_string(beam_id) + "/velocity");
        ofstream vel_file;
        vel_file.open(vel_name, ios_base::app);
        for (int i = 0; i < time_record.size(); i++) {
            vel_file << time_record[i] << " " << beam_vel_record[beam_id][i].transpose() << endl;
        }
        vel_file.close();
        string acc_name("./outputs/" + to_string(beam_id) + "/acceleration");
        ofstream acc_file;
        acc_file.open(acc_name, ios_base::app);
        for (int i = 0; i < time_record.size(); i++) {
            acc_file << time_record[i] << " " << beam_acc_record[beam_id][i].transpose() << endl;
        }
        acc_file.close();
        string elastic_name("./outputs/" + to_string(beam_id) + "/elastic");
        ofstream elastic_file;
        elastic_file.open(elastic_name, ios_base::app);
        for (int i = 0; i < time_record.size(); i++) {
            elastic_file << time_record[i] << " " << beam_elastic_record[beam_id][i].transpose() << endl;
        }
        elastic_file.close();
        string gravity_name("./outputs/" + to_string(beam_id) + "/gravity");
        ofstream gravity_file;
        gravity_file.open(gravity_name, ios_base::app);
        for (int i = 0; i < time_record.size(); i++) {
            gravity_file << time_record[i] << " " << beam_gravity_record[beam_id][i].transpose() << endl;
        }
        gravity_file.close();
        string damping_name("./outputs/" + to_string(beam_id) + "/damping");
        ofstream damping_file;
        damping_file.open(damping_name, ios_base::app);
        for (int i = 0; i < time_record.size(); i++) {
            damping_file << time_record[i] << " " << beam_damping_record[beam_id][i].transpose() << endl;
        }
        damping_file.close();
        string external_name("./outputs/" + to_string(beam_id) + "/external");
        ofstream external_file;
        external_file.open(external_name, ios_base::app);
        for (int i = 0; i < time_record.size(); i++) {
            external_file << time_record[i] << " " << beam_external_record[beam_id][i].transpose() << endl;
        }
        external_file.close();
        string distribution_name("./outputs/" + to_string(beam_id) + "/distribution");
        ofstream distribution_file;
        distribution_file.open(distribution_name, ios_base::app);
        for (int i = 0; i < time_record.size(); i++) {
            distribution_file << time_record[i] << " " << beam_distribute_record[beam_id][i].transpose() << endl;
        }
        distribution_file.close();
        string point_name("./outputs/" + to_string(beam_id) + "/point");
        ofstream point_file;
        point_file.open(point_name, ios_base::app);
        for (int i = 0; i < time_record.size(); i++) {
            point_file << time_record[i] << " " << beam_point_record[beam_id][i].transpose() << endl;
        }
        point_file.close();
    }
    beam_pos_record.clear();
    beam_vel_record.clear();
    beam_acc_record.clear();
    beam_elastic_record.clear();
    beam_gravity_record.clear();
    beam_external_record.clear();
    beam_damping_record.clear();
    beam_distribute_record.clear();
    beam_point_record.clear();
    time_record.clear();
}

void system_engine::check_output_folder() {

    if (!fs::exists("./outputs")) {
        // Create output folder
        fs::create_directories("./outputs");
    }
    for (auto const& x : beam_pos_record) {
        int beam_id = x.first;
        if (!fs::exists("./outputs/" + to_string(beam_id))) {
            fs::create_directories("./outputs/" + to_string(beam_id));
            // Create file for beam position, velocity, and acceleration
            ofstream pos_file("./outputs/" + to_string(beam_id) + "/position");
            pos_file.close();
            ofstream vel_file("./outputs/" + to_string(beam_id) + "/velocity");
            vel_file.close();
            ofstream acc_file("./outputs/" + to_string(beam_id) + "/acceleration");
            acc_file.close();
            ofstream elastic_file("./outputs/" + to_string(beam_id) + "/elastic");
            elastic_file.close();
            ofstream gravity_file("./outputs/" + to_string(beam_id) + "/gravity");
            gravity_file.close();
            ofstream damping_file("./outputs/" + to_string(beam_id) + "/damping");
            damping_file.close();
            ofstream external_file("./outputs/" + to_string(beam_id) + "/external");
            external_file.close();
            ofstream distribution_file("./outputs/" + to_string(beam_id) + "/distribution");
            distribution_file.close();
            ofstream point_file("./outputs/" + to_string(beam_id) + "/point");
            point_file.close();
        }
    }
}

int system_engine::check_write() {
    return time_record.size() >= 100;
}

void system_engine::load_switches() {
    ifstream switch_file;
    switch_file.open("./inputs/input_switches");
    string switch_line;
    while (getline(switch_file, switch_line)) {
        istringstream iss(switch_line);
        string switch_name;
        int switch_value;
        if (!(iss >> switch_name >> switch_value)) {
            cout << "Error reading switch file." << endl;
        }
        if (switch_name == "distributed") {
            distributed_force_switch = switch_value;
        } else if (switch_name == "point") {
            point_force_switch = switch_value;
        } else if (switch_name == "external") {
            external_load_switch = switch_value;
        } else if (switch_name == "damping") {
            damping_force_switch = switch_value;
        } else if (switch_name == "gravity") {
            gravity_force_switch = switch_value;
        }
    }
    switch_file.close();
}

void system_engine::beam_dist_force_off() {
    for (beam * b : beams) {
        VectorXd q = VectorXd::Zero(b->get_ndof());
        b->set_dist_force(q);
    }
}

void system_engine::beam_point_force_off() {
    for (beam * b : beams) {
        VectorXd q = VectorXd::Zero(b->get_ndof());
        b->set_point_force(q);
    }
}

void system_engine::beam_external_force_off() {
    for (beam * b : beams) {
        VectorXd q = VectorXd::Zero(b->get_ndof());
        b->set_external_force(q);
    }
}

void system_engine::beam_damping_force_off() {
    for (beam * b : beams) {
        VectorXd q = VectorXd::Zero(b->get_ndof());
        b->set_damping_force(q);
    }
}

void system_engine::beam_gravity_force_off() {
    for (beam * b : beams) {
        VectorXd q = VectorXd::Zero(b->get_ndof());
        b->set_gravity_force(q);
    }
}

void system_engine::remove_outputs() {
    if (!fs::exists("./outputs")) {
        // Remove output folder
        fs::remove_all("./outputs");
    }
}