#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>
#include <iomanip>
#include <chrono>

struct Ele {
    double mass;
    double x, y, z;
    double vx, vy, vz;
    double fx, fy, fz;
};

const double G = 6.67430e-11; // Gravitational constant
double dt; // Time step
int nSteps; // Number of time steps
int dInter; // Output state interval
int nEle; // Number of particles

std::vector<Ele> ini_ele(int nEle) {
    std::vector<Ele> prtce(nEle);  
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> pos_dist(-1.0, 1.0);
    std::uniform_real_distribution<> vel_dist(-0.1, 0.1);
    std::uniform_real_distribution<> mass_dist(1e24, 1e26);
    
    for (int i = 0; i < nEle; i++) {
        Ele a;
        a.mass = mass_dist(gen);
        a.x = pos_dist(gen);
        a.y = pos_dist(gen);
        a.z = pos_dist(gen);
        a.vx = vel_dist(gen);
        a.vy = vel_dist(gen);
        a.vz = vel_dist(gen);
        a.fx = a.fy = a.fz = 0.0;
        prtce[i] = a;  
    }
    return prtce;  
}

void compf(std::vector<Ele> &particles) {
    for (auto &a : particles) {
        a.fx = a.fy = a.fz = 0.0; // Reset forces
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

void update(std::vector<Ele> &particles) {
    for (auto &a : particles) {
        a.vx += (a.fx / a.mass) * dt;
        a.vy += (a.fy / a.mass) * dt;
        a.vz += (a.fz / a.mass) * dt;
        a.x += a.vx * dt;
        a.y += a.vy * dt;
        a.z += a.vz * dt;
    }
}

void output(const std::vector<Ele> &particles, std::ofstream &outfile) {
    outfile << particles.size();
    for (const auto a : particles) {
        outfile << "\t" << std::scientific << std::setprecision(10)
                << a.mass << "\t" << a.x << "\t" << a.y << "\t" << a.z
                << "\t" << a.vx << "\t" << a.vy << "\t" << a.vz
                << "\t" << a.fx << "\t" << a.fy << "\t" << a.fz;
    }
    outfile << "\n";
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <nEle> <dt> <nStep> <dInter>" << std::endl;
        return 1;
    }
    
    nEle = std::stoi(argv[1]);
    dt = std::stod(argv[2]);
    nSteps = std::stoi(argv[3]);
    dInter = std::stoi(argv[4]);
    
    std::vector<Ele> particles = ini_ele(nEle);
    std::ofstream outfile("solar.tsv");

    auto start = std::chrono::high_resolution_clock::now();  //measuring starts  here

    for (int step = 0; step < nSteps; ++step) {
        compf(particles);
        update(particles);
        if (step % dInter == 0) {
            output(particles, outfile);
        }
    }
   
    auto end = std::chrono::high_resolution_clock::now(); //and  ends here
    
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Run time was: " << time.count() << " milliseconds." << std::endl;

    return 0;
}
