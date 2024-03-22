#include "MultiplicativeCongruentialGenerator.hh"
#include <cmath>
#include <vector>
#include <numeric>
#include <math.h>
#include <iostream>


int main() {
    // parameters for the generator
    long long multiplier = pow(7, 5);
    long long modulus = pow(2, 31) - 1;
    long long initialSeed = modulus / 2;

    MultiplicativeCongruentialGenerator rng(multiplier, modulus, initialSeed);

    // physical parameters
    double xs = 0.1707; // cm^-1
    double trueMfp = 1/xs;
    double trueVariance = pow(xs, -2);
    std::cout << "Simulating paths of 100 keV x-rays incident on a semi-infinite water-equivalent parallelepiped..." << std::endl
              << "True mean free path: " << trueMfp << " ± " << pow(trueVariance, 0.5) << " cm (variance = " << trueVariance << " cm^2)" << std::endl;

    std::vector<double> mfps = {};
    std::vector<double> unbiasedVariances = {};
    std::vector<int> histories = {100, 1000, 10000};

    // assume slab is pure absorber and 10 cm thick; estimate uncollided leakage flux
    double thickness = 10.0;
    double trueUncollidedFlux = exp(-xs*thickness);
    double trueUncollidedFluxVariance = trueUncollidedFlux;
    int uncollidedParticles = 0; // number of particles whose first collision occurs at a depth beyond 10 cm

    // run simulations
    for (int N : histories) {
        // compute path lengths travelled before first collision
        std::vector<double> paths = {};
        for (int i = 0; i < N; i++) {
            double R = rng.getNextRandomNumber();
            double distance = -log(R)/xs;
            paths.push_back(distance);

            if (N == 10000) {
                // count particle for uncollided leakage flux if its first collision occurs at a depth beyond 10 cm
                if (distance > 10.0) {
                    uncollidedParticles++;
                }
            }
        }

        // compute statistics
        auto pathSum = std::accumulate(paths.begin(), paths.end(), 0.0);
        double mfp = pathSum/(double)N;
        std::vector<double> residuals = {};
        for (double path : paths) {
            double residual = pow(path - mfp, 2);
            residuals.push_back(residual);
        }
        double variance = std::accumulate(residuals.begin(), residuals.end(), 0.0)/(double)N;
        mfps.push_back(mfp);
        unbiasedVariances.push_back(variance);

        std::cout << "Monte Carlo estimate of mean free path for N = " << N << " histories: "
                  << mfp << " ± " << pow(variance, 0.5) << " cm (variance = " << variance << " cm^2)" << std::endl
                  << "\t" << 100.0*(trueMfp - mfp)/mfp << "\% error from true mean" << std::endl
                  << "\t" << 100.0*(trueVariance - variance)/variance << "\% error from true variance" << std::endl;
        
        if (N == 10000) {
            double estUncollidedFlux = (double)uncollidedParticles / (double)N;
            double estUncollidedFluxVariance = ( pow(0.0 - estUncollidedFlux, 2)*(1.0 - estUncollidedFlux) + pow(1.0 - estUncollidedFlux, 2)*(estUncollidedFlux) );
            std::cout << "Assuming the slab is a 10 cm thick pure absorber, the true uncollided leakage flux is: "
                      << trueUncollidedFlux << " ± " << pow(trueUncollidedFluxVariance, 0.5) << std::endl;
            std::cout << "Monte Carlo estimate of mean uncollided leakage flux for N = 10000 histories: "
                      << estUncollidedFlux << " ± " << pow(estUncollidedFluxVariance, 0.5) << std::endl
                      << "\t" << 100.0*(trueUncollidedFlux - estUncollidedFlux)/estUncollidedFlux << "\% error from true mean" << std::endl
                      << "\t" << 100.0*(trueUncollidedFluxVariance - estUncollidedFluxVariance)/estUncollidedFluxVariance << "\% error from true variance" << std::endl;
        }
    }

    return 0;
}
