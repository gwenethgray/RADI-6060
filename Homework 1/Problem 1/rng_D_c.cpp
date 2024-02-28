#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>

int main() {
    int Nvals[3] = { 100, 1000, 10000 };
    int numberOfBins = 10;
    int K = 10; // number of subdivisions to bin

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

        // perform the two-level frequency test
        std::cout << "Performing two-level frequency test with K=" << K << "..." << std::endl;
        std::vector<double> chiSquaredValues = {};
        double criticalValue = 16.919; // 9 degrees of freedom, 0.05 significance level

        // segment the set of random numbers into K groups of equal sizes
        int remainder = randomNumbers.size() % K;
        if (remainder > 0) {
            std::cout << "Current set of random numbers (size " << randomNumbers.size() << ") cannot be evenly split into "
                      << K << " groups. Neglecting the last " << remainder << " numbers generated." << std::endl;
        }
        int groupSize = randomNumbers.size() / K;

        for (int step = 0; step < K; step++) {
            // define subset of random numbers
            int start = step*groupSize;
            int end = (step + 1)*groupSize;
 
            // bin this subset
            std::vector<int> counts(numberOfBins, 0);
            for (int i = start; i < end; i++) {
                int binIndex = static_cast<int>(randomNumbers[i] * numberOfBins);
                // handle the edge case when R = 1.0
                if (binIndex == numberOfBins) {
                    binIndex--;
                }
                counts[binIndex]++;
            }

            double expectedCount = static_cast<double>(groupSize) / numberOfBins;
            double chiSquared = 0;

            // sum the squares of the residuals of the bin counts
            for (int count : counts) {
                double residual = count - expectedCount;
                chiSquared += pow(residual, 2) / expectedCount;
            }

            chiSquaredValues.push_back(chiSquared);

            std::string cond;
            if (chiSquared < criticalValue) {
                cond = "PASS";
            } else {
                cond = "FAIL";
            }

            std::cout << "Chi-squared value for group " << step << ":\t" << chiSquared << "\t\t" << cond << std::endl;
        }

        // calculate chi squared value of the set of chi squared values
        // maximum chi-squared value
        double maxChiSquared = *max_element(chiSquaredValues.begin(), chiSquaredValues.end());

        // bin the set of chi-squared values from 0 to maximum
        std::vector<int> counts(numberOfBins, 0);
        for (int i = 0; i < chiSquaredValues.size(); i++) {
            // scale this chi-squared value by the maximum to map it to [0,1) and multiply by the number of bins to turn it into an integer 
            int binIndex = static_cast<int>(chiSquaredValues[i] * numberOfBins / maxChiSquared);
            // handle the edge case when the mapping = 1.0
            if (binIndex == numberOfBins) {
                binIndex--;
            }
            counts[binIndex]++;
        }

        double expectedCount = 1.0;
        double level2chiSquared = 0;

        // sum the squares of the residuals of the bin counts
        for (int count : counts) {
            double residual = count - expectedCount;
            level2chiSquared += pow(residual, 2) / expectedCount;
        }

        std::string cond;
        if (level2chiSquared < criticalValue) {
            cond = "passes the two-level frequency test.";
        } else {
            cond = "fails the two-level frequency test.";
        }

        std::cout << "Second level chi-squared = " << level2chiSquared << "\t\t->\tThis RNG " << cond << std::endl;
        std::cout << std::endl;
    }

    return 0;
}