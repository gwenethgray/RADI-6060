#include "Particle.hh"
#include "MultiplicativeCongruentialGenerator.hh"
#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <map>

// parameters for the generator
long long multiplier = pow(7, 5);
long long modulus = pow(2, 31) - 1;
long long initialSeed = modulus / 2;

MultiplicativeCongruentialGenerator rng(multiplier, modulus, initialSeed);

// cross sections for energy in keV
std::map<double, double> xs_a{
	{50.0, 2.725*pow(10, -2)},
	{100.0, 2.763*pow(10, -3)}
};
std::map<double, double> xs_s{
	{50.0, 1.803*pow(10, -1)},
	{100.0, 1.626*pow(10, -1)}
};
std::map<double, double> xs_t{
	{50.0, xs_a[50.0] + xs_s[50.0]},
	{100.0, xs_a[100.0] + xs_s[100.0]}
};
// medium with above cross sections
std::map<std::string, std::map<double, double>> medium{
	{"xs_a", xs_a},
	{"xs_s", xs_s},
	{"xs_t", xs_t}
};

int main() {
    int N_histories = 10;
    std::cout << "Simulating tracks of " << N_histories << " particles in water..." << std::endl;

    int N_backscattered = 0;
    int N_absorbed = 0;
    int N_uncollided = 0;
    int N_scattered = 0;

    double initial_x = 0.0;
    double initial_y = 0.0;
    double initial_z = 0.0;
    double initial_E = 100.0; // keV
    double initial_mu = 1.0;
    double initial_phi = 0.0; // radians
    std::vector<double> initial_omega = {1.0, 0.0, 0.0};

    for (int i = 0; i < N_histories; i++) {
        Particle p(initial_x, initial_y, initial_z, initial_E, initial_mu, initial_phi, initial_omega, rng);
        std::cout << "Simulating track of particle " << i << "..." << std::endl;
        while (p.getX() >= 0.0 && p.getX() <= 2.0) {
            int result = p.simulate(medium);
            if (result == -1) { N_backscattered++; }
            else if (result == 0) { N_absorbed++; }
            else if (result == 1) { N_uncollided++; }
            else if (result == 2) { N_scattered++; }
        }
    }

    std::cout << "# particles backscattered: " << N_backscattered << std::endl
              << "# particles absorbed: " << N_absorbed << std::endl
              << "# particles transmitted without colliding: " << N_uncollided << std::endl
              << "# particles transmitted after colliding: " << N_scattered << std::endl;
}