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
    double trueMedianPosn = radius/pow(2.0, 1.0/3.0);
    double eta = trueMedianPosn/radius;
    double medVarCoeff = (3.0/5.0) - (3.0*eta/2.0) + pow(eta, 2);
    double trueMedianVariance = medVarCoeff*pow(radius, 2);
    std::cout << "Simulating radial positions of 100 keV source particles generated in a homogeneous water-equivalent sphere of radius 10 cm..." << std::endl
              << "True mean radial position: " << trueMeanPosn << " ± " << pow(trueVariance, 0.5) << " cm (variance = " << trueVariance << " cm^2)" << std::endl
              << "True median radial position: " << trueMedianPosn << " ± " << pow(trueMedianVariance, 0.5) << " cm (variance = " << trueMedianVariance << " cm^2)" << std::endl;

    std::vector<double> meanPosns = {};
    std::vector<double> unbiasedVariances = {};
    std::vector<double> medianPosns = {};
    std::vector<double> unbiasedMedianVariances = {};
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
        double posnSum = std::accumulate(posns.begin(), posns.end(), 0.0);
        double meanPosn = posnSum/(double)N;
        double meanPctDiff = 100*(trueMeanPosn - meanPosn)/meanPosn;
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
        double medianPctDiff = 100*(trueMedianPosn - medianPosn)/medianPosn;
        medianPosns.push_back(medianPosn);
        residuals = {};
        for (double posn : posns) {
            double residual = pow(posn - medianPosn, 2);
            residuals.push_back(residual);
        }
        double medianVariance = std::accumulate(residuals.begin(), residuals.end(), 0.0)/(double)N;
        unbiasedMedianVariances.push_back(medianVariance);

        std::cout << "Monte Carlo estimates for N = " << N << " histories:" << std::endl;
        std::cout << "\tMean radial position: " << meanPosn << " ± " << pow(variance, 0.5) << " cm (variance = " << variance << " cm^2)" << std::endl;
        std::cout << "\t\t" << meanPctDiff << "\% error from true mean" << std::endl;
        std::cout << "\tMedian radial position: " << medianPosn << " ± " << pow(medianVariance, 0.5) << " cm (variance = " << medianVariance << " cm^2)" << std::endl;
        std::cout << "\t\t" << medianPctDiff << "\% error from true median" << std::endl;
    }

    return 0;
}
