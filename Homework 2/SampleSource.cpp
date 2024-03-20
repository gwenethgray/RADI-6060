#include "MultiplicativeCongruentialGenerator.hh"
#include <cmath>
#include <vector>
#include <numeric>
#include <math.h>
#include <iostream>
#include <algorithm>


int main() {
    // parameters for the generator
    long long multiplier = pow(7, 5);
    long long modulus = pow(2, 31) - 1;
    long long initialSeed = modulus / 2;

    MultiplicativeCongruentialGenerator rng(multiplier, modulus, initialSeed);

    // physical parameters
    double xs = 0.1707; // cm^-1
    double radius = 10.0; // cm
    double trueMeanPosn = 3.0*radius/4.0;
    double trueVariance = 3.0*pow(radius, 2)/80.0;
    std::cout << "True mean radial position of 100 keV source particles in homogeneous water-equivalent sphere of radius 10 cm: "
              << trueMeanPosn << " ± " << pow(trueVariance, 0.5) << " cm (variance = " << trueVariance << " cm)" << std::endl;

    std::vector<double> meanPosns = {};
    std::vector<double> unbiasedVariances = {};
    std::vector<double> medianPosns = {};
    std::vector<int> histories = {100, 1000, 10000};

    // run simulations
    for (int N : histories) {
        // compute radial positions of randomly generated source particles
        std::vector<double> posns = {};
        for (int i = 0; i < N; i++) {
            double R = rng.getNextRandomNumber();
            double posn = radius*pow(R, 1.0/3.0);
            posns.push_back(posn);
        }

        // compute statistics
        auto posnSum = std::accumulate(posns.begin(), posns.end(), 0.0);
        double meanPosn = posnSum/(double)N;
        std::vector<double> residuals = {};
        for (double posn : posns) {
            double residual = pow(posn - meanPosn, 2);
            residuals.push_back(residual);
        }
        double variance = std::accumulate(residuals.begin(), residuals.end(), 0.0)/(double)N;
        meanPosns.push_back(meanPosn);
        unbiasedVariances.push_back(variance);

        std::sort(posns.begin(), posns.end());
        double medianPosn = posns[N/2];
        medianPosns.push_back(medianPosn);

        std::cout << "Monte Carlo estimate of mean radial position for N = " << N << " histories: "
                  << meanPosn << " ± " << pow(variance, 0.5) << " cm (variance = " << variance << " cm)" << std::endl;
        std::cout << "Monte Carlo estimate of median radial position for N = " << N << " histories: " << medianPosn << " cm" << std::endl;
    }

    return 0;
}
