#include "Particle.hh"
#include "MultiplicativeCongruentialGenerator.hh"
#include <cmath>
#include <string>
#include <map>

#ifndef PI
#define PI 3.14159265358979323846
#endif

// particle constructor
Particle::Particle(double x, double y, double z, double E, double mu, double phi, std::vector<double> omega,
                   MultiplicativeCongruentialGenerator rng) : x(x), y(y), z(z), E(E), mu(mu), phi(phi), omega(omega), rng(rng) {}

// get particle depth in medium
double Particle::getX() {
    return x;
}

// sample collisional pdf to determine where the next collision will occur, given a medium M
void Particle::computeNextCollisionPosition(
	std::map<std::string, std::map<double, double>> M) {
	
	double R = rng.getNextRandomNumber();
    double current_R = rng.getRandomNumber();
    std::cout << "Current random number: " << current_R << std::endl;
	double collision_distance = -log(R/M["xs_t"][E])/M["xs_t"][E];
    double new_x = x + omega[0]*collision_distance;
    double new_y = y + omega[1]*collision_distance;
    double new_z = z + omega[2]*collision_distance;
    x = new_x;
    y = new_y;
    z = new_z;
}

// sample non-absorption probability to determine whether an interaction will be absorption or scattering
bool Particle::isAbsorbed(
    std::map<std::string, std::map<double, double>> M) {
    
    double R = rng.getNextRandomNumber();
    return R <= M["xs_a"][E]/M["xs_t"][E];
};

// sample angular distribution (polar angle and azimuthal angle) to determine new direction for particle in directional coordinate system
void Particle::computeNextAngle(
    std::map<std::string, std::map<double, double>> M) {
    
    // pdf(mu) = 1 + mu^2
    // not normalized
    // cdf(mu) = mu + (1/3)*mu^3 = R
    // mu =  ??? not invertible
    double R_polar = rng.getNextRandomNumber();
    double new_mu = R_polar; // in directional coordinate system

    // azimuthal angle uniformly distributed in 2pi
    double R_azimuthal = rng.getNextRandomNumber();
    double new_phi = 2.0*PI*R_azimuthal; // in directional coordinate system

    // compute directional cosines
    double sin_theta = pow(1.0 - pow(new_mu, 2), 0.5);
    double new_omega_x = omega[0]*new_mu +
                         (omega[2]*omega[0]*sin_theta*cos(new_phi)/pow(1 - pow(omega[2], 2), 0.5)) +
                         (omega[1]*sin_theta*sin(new_phi)/pow(1 - pow(omega[2], 2), 0.5));
    double new_omega_y = omega[1]*new_mu +
                         (omega[2]*omega[1]*sin_theta*cos(new_phi)/pow(1 - pow(omega[2], 2), 0.5)) +
                         (omega[0]*sin_theta*sin(new_phi)/pow(1 - pow(omega[2], 2), 0.5));
    double new_omega_z = omega[2]*new_mu - sin_theta*cos(new_phi)*pow(1 - pow(omega[2], 2), 0.5);
    // convert to fixed coordinate system
    omega = {new_omega_x, new_omega_y, new_omega_z};
}


// simulate particle's track in a medium M
int Particle::simulate(
	std::map<std::string, std::map<double, double>> M) {

    std::cout << "Simulating particle track..." << std::endl;

    bool firstCollision = true;
    while (x >= 0.0 && x <= 2.0) {
        // sample collisional pdf to determine distance to next interaction
        computeNextCollisionPosition(M);
        std::cout << "New position: (" << x << ", " << y << ", " << z << ")" << std::endl;

        // test whether interaction is outside the slab
        if (x < 0.0) {
            return -1; // backscatter
        } else if (x > 2.0) {
            if (firstCollision) {
                return 1; // transmitted (uncollided)
            } else {
                return 2; // transmitted (scattered)
            }
        }
        
        // sample non-absorption probability to determine whether particle is absorbed or scattered
        if (isAbsorbed(M)) {
            return 0; // stopped in the medium
        }

        //TODO: compute loss of energy

        // sample angular pdf to determine change in direction
        computeNextAngle(M);

        firstCollision = false; // faster than checking
    }
}