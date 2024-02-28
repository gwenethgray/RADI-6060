#include <iostream>
#include <vector>
#include "MultiplicativeCongruentialGenerator.hh"
#include <cassert>

#ifndef PI
#define PI 3.14159265358979323846
#endif

#ifndef NUCLIDE
#define NUCLIDE

// a nuclear species with its own energy dependent cross sections
class Nuclide {
public:
    // cross sections with respect to energy (must be same length)
    // assume cross sections are independent of direction
    std::vector<double> absorptionCrossSection;
    std::vector<double> scatterCrossSection;
    std::vector<double> totalCrossSection;
    std::vector<double> energy;
    // probability of scatter event being elastic vs inelastic
    double p_elastic;

    // nuclide constructor
    Nuclide(std::vector<double> AXS, std::vector<double> SXS, std::vector<double> energy, double p_elastic);
};

#endif

// nuclide constructor
Nuclide::Nuclide(std::vector<double> AXS, std::vector<double> SXS, std::vector<double> energy, double p_elastic) {
    // confirm that the cross sections are the same length as the values in the energy range
    assert(absorptionCrossSection.size() == scatterCrossSection.size() == energy.size());
    // compute total cross section by summing absorption and scatter cross sections
    totalCrossSection = {};
    for (int i=0; i < AXS.size(); i++) {
        totalCrossSection.push_back(AXS[i] + SXS[i]);
    }
}


#ifndef MATERIAL
#define MATERIAL

// a material consisting of one or more nuclides in relative proportions by mass
class Material {
public:
    // list of component nuclides each with their own respective scatter and absorption xs
    std::vector<Nuclide> nuclides;
    std::vector<double> massFractions;
    // composite cross sections with respect to energy (must be same length)
    // assume cross sections are independent of direction
    std::vector<double> compositeAbsorptionCrossSection;
    std::vector<double> compositeScatterCrossSection;
    std::vector<double> compositeTotalCrossSection;
    std::vector<double> energy;

    // constructor
    Material(std::vector<Nuclide> nuclides, std::vector<double> massFractions);
};

#endif

// material constructor
Material::Material(std::vector<Nuclide> nuclides, std::vector<double> massFractions) {
    // confirm that nuclides and mass fractions are 1-to-1
    assert(nuclides.size() == massFractions.size());
    // confirm the nuclides all have the same number of energy values
    int numberOfEnergies = nuclides[0].energy.size();
    for (Nuclide nuclide : nuclides) {
        assert(nuclide.energy.size() == numberOfEnergies);
    }
    // compute composite cross sections, assuming the nuclides all have cross section vectors that are the same length
    compositeAbsorptionCrossSection = {};
    compositeScatterCrossSection = {};
    compositeTotalCrossSection = {};
    for (int i=0; i < numberOfEnergies; i++) {
        double AXS = 0;
        double SXS = 0;
        double TXS = 0;
        for (int j=0; j < nuclides.size(); j++) {
            AXS += nuclides[j].absorptionCrossSection[i] * massFractions[j];
            SXS += nuclides[j].scatterCrossSection[i] * massFractions[j];
            TXS += nuclides[j].totalCrossSection[i] * massFractions[j];
        }
    }
}


#ifndef RANDOM_WALK
#define RANDOM_WALK

class RandomWalk {
private:
    // positional coordinates
    double x;
    double y;
    double z;
    // energy
    double E;
    // directional coordinates
    double theta;
    double phi;
    // random number generator
    MultiplicativeCongruentialGenerator rng;
    // assume homogeneous medium
    Material medium;
public:
    // constructor
    RandomWalk(
        double initialX, double initialY, double initialZ,
        double initialEnergy,
        double initialPolarAngle, double initialAzimuthalAngle,
        MultiplicativeCongruentialGenerator rng, Material medium);
    // initialize (r0,E0,omega0) by sampling fQ (source pdf)

    // GET PARAMETERS

    // get position vector
    std::vector<double> getPosition() const;
    // get energy
    double getEnergy() const;
    // get direction vector
    std::vector<double> getDirection() const;

    // determine which energy bin in the medium's xs data the particle's current energy falls within
    int currentEnergyBinNumber() const;

    // COMPUTE NEW PARAMETERS

    // sample ratio of absorption xs to total xs to determine whether the particle is absorbed in an interaction
    bool isAbsorbed();
    // determine which material the interaction occurs in
    int materialOfInteraction();
    // determine the type of scatter reaction (elastic vs inelastic) based on the probability for material m
    std::string typeOfInteraction(int m);
    // sample single differential for material m to get new scatter cosine and azimuth (assume azimuthal isotropy)
    // return (theta,phi)
    std::vector<double> newDirection(int m);
    // sample fT to get distance s to next collision site
    double distanceToNextCollision();
};

#endif

// random walk constructor
RandomWalk::RandomWalk(
    double initialX, double initialY, double initialZ,
    double initialEnergy,
    double initialPolarAngle, double initialAzimuthalAngle,
    MultiplicativeCongruentialGenerator rng, Material medium) {}

// get current position of particle as (x,y,z) coordinates
std::vector<double> RandomWalk::getPosition() const {
    return { x, y, z };
}
// get current energy of particle
double RandomWalk::getEnergy() const {
    return E;
}
// get current direction of particle as (theta,phi) coordinates
std::vector<double> RandomWalk::getDirection() const {
    return { theta, phi };
}

// determine which energy bin in the medium's xs data the particle's current energy falls within
int RandomWalk::currentEnergyBinNumber() const {
    int bin = 0;
    for (int i=0; i < medium.energy.size() - 1; i++) {
        if (E < medium.energy[0]) {
            // current energy is lower than minimum of the medium's xs data energy range; use lowest bin
            break;
        } else if (E >= medium.energy[i] && E < medium.energy[i+1]) {
            // found the energy bin the current energy is in
            bin = i;
            break;
        } else if (i == medium.energy.size() - 2) {
            // current energy is higher than maximum of the medium's xs data energy range; use highest bin
            bin = medium.energy.size();
            break;
        }
    }
    return bin;
}

// determine whether the particle is absorbed in the current interaction
bool RandomWalk::isAbsorbed() {
    // check the bin
    int bin = currentEnergyBinNumber();
    // compute probability of absorption in the medium at the particle's current energy
    double p_absorb = medium.compositeAbsorptionCrossSection[bin] / medium.compositeTotalCrossSection[bin];
    // generate a random number between 0 and 1
    double xi1 = rng.getNextRandomNumber();
    // test whether the random number is greater than the probability of absorption at the current energy
    return xi1 < p_absorb;
}

// determine which material the interaction occurs in
int RandomWalk::materialOfInteraction() {
    // check the bin
    int bin = currentEnergyBinNumber();
    // generate a random number between 0 and 1
    double xi2 = rng.getNextRandomNumber();
    // nuclides are assigned to bins of width equal to the fraction of the composite scatter xs that they contribute
    // accumulate the sum of consecutive bin widths from first to last until one is found that contains the random number
    double upperBound = 0;
    for (int m=0; m < medium.nuclides.size(); m++) {
        // compute this nuclide's contributed fraction of the composite scatter xs and add it to the upper bound
        upperBound += medium.nuclides[m].scatterCrossSection[bin] / medium.compositeScatterCrossSection[bin];
        if (xi2 < upperBound) {
            return m;
        }
    }
}

// determine the type of scatter reaction (elastic vs inelastic) based on the probability for material m
std::string RandomWalk::typeOfInteraction(int m) {
    // generate a random number between 0 and 1
    double xi3 = rng.getNextRandomNumber();
    if (xi3 <= medium.nuclides[m].p_elastic) {
        return "Elastic";
    } else {
        return "Inelastic";
    }
}

// sample single differential for material m to get new scatter cosine and azimuth (assume azimuthal isotropy)
// return (theta,phi)
std::vector<double> RandomWalk::newDirection(int m) {
    // generate a random number between 0 and 1
    double xi4 = rng.getNextRandomNumber();
    //TODO compute single diff, which is marginal pdf of the double diff scatter probability function
    //TODO then set xi4 equal to the CDF of this (with factor of 2Pi to account for azimuthal isotropy)
    double mu; // = CDF evaluated at xi4
    // generate another random number between 0 and 1
    double xi5 = rng.getNextRandomNumber();
    // map the random number to a number between 0 and 2Pi to sample the uniform distribution of azimuthal angles
    double phi = 2*PI*xi5;
    return { mu, phi };
    //TODO step 5 in the Analog Monte Carlo notes has some equations for calculating the new direction in the fixed coordinate system from these values
}

//TODO sample new energy

// sample fT to get distance s to next collision site
double RandomWalk::distanceToNextCollision() {
    // generate a random number between 0 and 1

}





class MultiplicativeCongruentialGenerator {
private:
    long long seed; // current random number sequences becomes seed for next
    long long I1; // multiplier
    long long M; // modulus
    double randomNumber; // current random number between 0 and 1
    std::vector<long long> seeds; // cumulative set of random number sequences between 0 and M
    std::vector<double> randomNumbers; // cumulative set of random numbers between 0 and 1

public:
    MultiplicativeCongruentialGenerator(long long multiplier, long long modulus, long long initialSeed);

    // compute the next random number in the sequence and store it
    void computeNext();

    // compute and return the next random number
    double MultiplicativeCongruentialGenerator::getNextRandomNumber();

    // get current seed
    long long getSeed() const;

    // get current random number
    double getRandomNumber() const;

    // get all stored seeds
    const std::vector<long long>& getSeeds() const;

    // get all stored random numbers
    const std::vector<double>& getRandomNumbers() const;
};