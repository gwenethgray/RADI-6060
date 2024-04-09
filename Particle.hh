#include <map>
#include <string>
#include "MultiplicativeCongruentialGenerator.hh"

#ifndef PARTICLE
#define PARTICLE

class Particle {
private:
	// positional coordinates
	double x;
	double y;
	double z;
	// energy
	double E;
	// polar angle
	double mu;
	// azimuthal angle
	double phi;

	// directional cosines
	std::vector<double> omega;

    // random number generator
    MultiplicativeCongruentialGenerator rng;

public:
    // constructor
	Particle(double x, double y, double z, double E, double mu, double phi, std::vector<double> omega, MultiplicativeCongruentialGenerator rng);

	// get particle depth in medium
	double getX();

	// sample collisional pdf to determine where the next collision will occur, given a medium M
	void computeNextCollisionPosition(
		std::map<std::string, std::map<double, double>> M
	);

    // sample non-absorption probability to determine whether an interaction will be absorption or scattering
    bool isAbsorbed(
		std::map<std::string, std::map<double, double>> M
	);

	// sample angular distribution (polar angle and azimuthal angle) to determine new direction for particle in directional coordinate system
	void computeNextAngle(
		std::map<std::string, std::map<double, double>> M
	);

    // simulate particle's track in a medium M
	/*	Output		Result Type
		-1			backscattered
		0			absorbed in medium
		1			transmitted (uncollided)
		2			transmitted (scattered)
	 */
    int simulate(
		std::map<std::string, std::map<double, double>> M
    );
};

#endif