#include "system_engine.h"

system_engine::system_engine() {
    f_engine = new force_engine();
    t = 0.0;
}

void system_engine::read_beams() {
    ifstream beam_file;
    beam_file.open("input_beam");
    ifstream beam_init_position;
    beam_init_position.open("input_beam_position");
    ifstream beam_init_velocity;
    beam_init_velocity.open("input_beam_velocity");
    int def_vel = 0;
    if (beam_init_velocity) {
        def_vel = 1;
    }
    ifstream beam_init_acceleration;
    beam_init_acceleration.open("input_beam_acceleration");
    int def_acc = 0;
    if (beam_init_acceleration) {
        def_acc = 1;
    }
    string beam_line;
    string beam_pos;
    string beam_vel;
    string beam_acc;
    while (getline(beam_file, beam_line)) {
        istringstream iss(beam_line);
        double E, rho, length, thick, width;
        int nelement, botCnstr, topCnstr;
        if (!(iss >> E >> rho >> length >> thick >> width >> nelement >> botCnstr >> topCnstr)) {
            cout << "Error reading beam file." << endl;
        }
        beam_builder * b_builder = new beam_builder();
        b_builder->set_E(E);
        b_builder->set_rho(rho);
        b_builder->set_length(length);
        b_builder->set_thick(thick);
        b_builder->set_width(width);
        b_builder->set_nelement(nelement);
        b_builder->set_botCnstr(botCnstr);
        b_builder->set_topCnstr(topCnstr);
        int ndof = b_builder->get_beam()->get_ndof();
        getline(beam_init_position, beam_pos);
        istringstream iss_pos(beam_pos);
        VectorXd pos(ndof);
        for (int i = 0; i < ndof; i++) {
            double cur_pos;
            iss_pos >> cur_pos;
            pos(i) = cur_pos;
        }
        b_builder->set_position(pos);
        if (def_vel) {
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
        if (def_acc) {
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
        beam * b = b_builder->get_beam();
        beams.push_back(b);
    }
    beam_file.close();
    beam_init_position.close();
    if (!def_vel) {
        beam_init_velocity.close();
    }
    if (!def_acc) {
        beam_init_acceleration.close();
    }
}

void system_engine::test_read_beams() {
    for (beam * b : beams) {
        cout << *b;
    }
}

void system_engine::read_particles() {
    ifstream par_file;
    par_file.open("input_particle");
    ifstream par_init_position;
    par_init_position.open("input_particle_position");
    ifstream par_init_velocity;
    par_init_velocity.open("input_particle_velocity");
    int def_vel = 0;
    if (par_init_velocity) {
        def_vel = 1;
    }
    ifstream par_init_acceleration;
    par_init_acceleration.open("input_particle_acceleration");
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
            cout << "Error reading beam file." << endl;
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
        if (def_vel) {
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
        if (def_acc) {
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
    if (!def_vel) {
        par_init_velocity.close();
    }
    if (!def_acc) {
        par_init_acceleration.close();
    }
}

void system_engine::test_read_particles() {
    for (particle * p : particles) {
        cout << *p;
    }
}

void system_engine::read_fluid_field() {
    ifstream ff_file;
    ff_file.open("input_fluid_field");
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
    }
}

void system_engine::update_beam_elastic_forces() {
    for (beam * b : beams) {
        f_engine->elastic_force(b);
    }
}

void system_engine::update_beam_dist_forces() {
    for (beam * b : beams) {
        f_engine->distributed_load(b, f_field, t);
    }
}