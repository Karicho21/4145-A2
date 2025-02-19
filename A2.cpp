#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>
#include <iomanip>

struct Particle {
    double mass;
    double x, y, z;
    double vx, vy, vz;
    double fx, fy, fz;
};

const double G = 6.67430e-11; // Gravitational constant
double dt; // Time step
int num_steps; // Number of time steps
int dump_interval; // Output state interval
int num_particles; // Number of particles

std::vector<Particle> initialize_particles(int num_particles) {
    std::vector<Particle> particles(num_particles);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> pos_dist(-1.0, 1.0);
    std::uniform_real_distribution<> vel_dist(-0.1, 0.1);
    std::uniform_real_distribution<> mass_dist(1e24, 1e26);
    
    for (auto &p : particles) {
        p.mass = mass_dist(gen);
        p.x = pos_dist(gen);
        p.y = pos_dist(gen);
        p.z = pos_dist(gen);
        p.vx = vel_dist(gen);
        p.vy = vel_dist(gen);
        p.vz = vel_dist(gen);
        p.fx = p.fy = p.fz = 0.0;
    }
    return particles;
}

void compute_forces(std::vector<Particle> &particles) {
    for (auto &p : particles) {
        p.fx = p.fy = p.fz = 0.0; // Reset forces
    }
    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            double dx = particles[j].x - particles[i].x;
            double dy = particles[j].y - particles[i].y;
            double dz = particles[j].z - particles[i].z;
            double dist_sq = dx * dx + dy * dy + dz * dz + 1e-10;
            double dist = std::sqrt(dist_sq);
            double force = (G * particles[i].mass * particles[j].mass) / dist_sq;
            double fx = force * dx / dist;
            double fy = force * dy / dist;
            double fz = force * dz / dist;
            particles[i].fx += fx;
            particles[i].fy += fy;
            particles[i].fz += fz;
            particles[j].fx -= fx;
            particles[j].fy -= fy;
            particles[j].fz -= fz;
        }
    }
}

void update_particles(std::vector<Particle> &particles) {
    for (auto &p : particles) {
        p.vx += (p.fx / p.mass) * dt;
        p.vy += (p.fy / p.mass) * dt;
        p.vz += (p.fz / p.mass) * dt;
        p.x += p.vx * dt;
        p.y += p.vy * dt;
        p.z += p.vz * dt;
    }
}

void output_state(const std::vector<Particle> &particles, std::ofstream &outfile) {
    outfile << particles.size();
    for (const auto &p : particles) {
        outfile << "\t" << std::scientific << std::setprecision(10)
                << p.mass << "\t" << p.x << "\t" << p.y << "\t" << p.z
                << "\t" << p.vx << "\t" << p.vy << "\t" << p.vz
                << "\t" << p.fx << "\t" << p.fy << "\t" << p.fz;
    }
    outfile << "\n";
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <num_particles> <dt> <num_steps> <dump_interval>" << std::endl;
        return 1;
    }
    
    num_particles = std::stoi(argv[1]);
    dt = std::stod(argv[2]);
    num_steps = std::stoi(argv[3]);
    dump_interval = std::stoi(argv[4]);
    
    std::vector<Particle> particles = initialize_particles(num_particles);
    std::ofstream outfile("solar.tsv");
    for (int step = 0; step < num_steps; ++step) {
        compute_forces(particles);
        update_particles(particles);
        if (step % dump_interval == 0) {
            output_state(particles, outfile);
        }
    }
    outfile.close();
    std::cout << "Simulation completed. Output saved to solar.tsv" << std::endl;
    return 0;
}
