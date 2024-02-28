#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <random>

int main() {
    int Nvals[3] = { 100, 1000, 10000 };
    int numberOfBins = 10;
    int numberOfMoments = 5;

    // evaluate the generator's randomness using 10^2, 10^3, and 10^4 random numbers
    for (int N : Nvals) {

        std::cout << "Evaluating C++ intrinsic RNG algorithm at N = " << N << " seeds..." << std::endl;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        std::vector<double> randomNumbers = {};

        // generate N random numbers
        for (int i = 0; i < N; i++) {
            randomNumbers.push_back(dis(gen));
        }

        // show bin counts
        int nbins = 10;
        std::cout << "Binning random numbers..." << std::endl << "BIN:\t";
        for (int i = 1; i < 11; i++) { std::cout << i << "\t"; }
        std::cout << std::endl << "COUNT:\t";
        std::vector<int> counts(numberOfBins, 0);
        for (double R : randomNumbers) {
            int binIndex = static_cast<int>(R * numberOfBins);
            // handle the edge case when R = 1.0
            if (binIndex == numberOfBins) {
                binIndex--;
            }
            counts[binIndex]++;
        }
        for (int binCount : counts) {
            std::cout << binCount << "\t";
        }
        std::cout << std::endl;

        // perform frequency test
        std::cout << "Performing frequency test..." << std::endl;
        double expectedCount = static_cast<double>(randomNumbers.size()) / numberOfBins;
        double chiSquared = 0;
        // sum the squares of the residuals of the bin counts
        for (int count : counts) {
            double residual = count - expectedCount;
            chiSquared += pow(residual, 2) / expectedCount;
        }
        double criticalValue = 16.919; // N - 1 = 9 degrees of freedom, 0.05 significance level
        std::string cond;
        if (chiSquared < criticalValue) {
            cond = "passes the frequency test.";
        } else {
            cond = "fails the frequency test.";
        }
        std::cout << "Chi-squared = " << chiSquared << "\t\t->\tThis RNG " << cond << std::endl;

        // perform moments test
        const double threshold = 0.01; // upper limit of acceptable difference between computed and expected moment
        std::cout << "Performing moments test (threshold <" << threshold << ")..." << std::endl;
        std::cout << "Moment:\t\tExpected:\tComputed:\tDifference:\tResult:" << std::endl;
        bool allMomentsPass = true;
        for (int k = 1; k < numberOfMoments + 1; k++) {
            // compute k-th moment
            double moment = 0;
            for (double R : randomNumbers) {
                moment += pow(R, k);
            }
            moment /= randomNumbers.size();
            double expected = 1.0 / (k + 1);
            double difference = fabs(expected - moment);
            std::string result;
            if (difference < threshold) {
                result = "PASS";
            } else {
                result = "FAIL";
                allMomentsPass = false;
            }
            std::cout << " k=" << k << ":\t\t"
                      << std::setprecision(4) << std::fixed << expected << "\t\t"
                      << std::setprecision(4) << std::fixed << moment << "\t\t"
                      << std::setprecision(4) << std::fixed << difference << "\t\t"
                      << result << std::endl;
        }
        if (allMomentsPass) {
            std::cout << "This RNG passes the moments test." << std::endl;
        } else {
            std::cout << "This RNG fails the moments test." << std::endl;
        }

        std::cout << std::endl;
    }

    return 0;
}