#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

class MultiplicativeCongruentialGenerator {
private:
    long long seed; // current random number sequences becomes seed for next
    long long I1; // multiplier
    long long M; // modulus
    double randomNumber; // current random number between 0 and 1
    std::vector<long long> seeds; // cumulative set of random number sequences between 0 and M
    std::vector<double> randomNumbers; // cumulative set of random numbers between 0 and 1

public:
    MultiplicativeCongruentialGenerator(long long multiplier, long long modulus, long long initialSeed)
    : I1(multiplier), M(modulus), seed(initialSeed % modulus), randomNumber(static_cast<double>(initialSeed % modulus) / modulus) {
        seeds = { seed };
        randomNumbers = { randomNumber };
    }

    // compute the next random number in the sequence and store it
    void computeNext() {
        seed = (I1 * seed) % M;
        seeds.push_back(seed);
        randomNumber = static_cast<double>(seed) / M; // next random number from 0 to 1
        randomNumbers.push_back(randomNumber);
    }

    // get current seed
    long long getSeed() const {
        return seed;
    }

    // get current random number
    double getRandomNumber() const {
        return randomNumber;
    }

    // get all stored seeds
    const std::vector<long long>& getSeeds() const {
        return seeds;
    }

    // get all stored random numbers
    const std::vector<double>& getRandomNumbers() const {
        return randomNumbers;
    }

    // divide stored random numbers into bins of equal width and return the bin counts
    std::vector<int> binRandomNumbers(int numberOfBins) const {
        std::vector<int> counts(numberOfBins, 0);
        
        for (double R : randomNumbers) {
            int binIndex = static_cast<int>(R * numberOfBins);
            // handle the edge case when R = 1.0
            if (binIndex == numberOfBins) {
                binIndex--;
            }
            counts[binIndex]++;
        }

        return counts;
    }

    // perform the moments test to evaluate the effectiveness of this random number generator
    void runMomentsTest(int numberOfMoments) {
        std::cout << "Performing moments test..." << std::endl;
        std::cout << "Moment:\t\tExpected:\tComputed:\tDifference:\tResult:" << std::endl;
        const double threshold = 0.01; // upper limit of acceptable difference between computed and expected moment
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
    }

    // perform the frequency test to evaluate the effectiveness of this random number generator
    void runFrequencyTest(int numberOfBins = 10) {
        std::cout << "Performing frequency test..." << std::endl;
        std::vector<int> counts = binRandomNumbers(numberOfBins);
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
    }

    // perform the two-level frequency test
    void runTwoLevelFrequencyTest(int K = 10, int numberOfBins = 10) {
        std::cout << "Performing two-level frequency test with K=" << K << "..." << std::endl;
        std::vector<double> chiSquaredValues = {};
        double criticalValue = 16.919; // N - 1 = 9 degrees of freedom, 0.05 significance level

        // segment the set of random numbers into K groups of equal sizes N
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
        double level2chiSquared = 0;
        for (double chiSquared : chiSquaredValues) {
            double residual = chiSquared - criticalValue;
            level2chiSquared += pow(residual, 2) / criticalValue;
        }

        std::string cond;
        if (level2chiSquared < criticalValue) {
            cond = "passes the two-level frequency test.";
        } else {
            cond = "fails the two-level frequency test.";
        }

        std::cout << "Second level chi-squared = " << level2chiSquared << "\t\t->\tThis RNG " << cond << std::endl;
    }
};

int main() {
    // parameters for the generator
    long long multiplier = pow(7, 5);
    long long modulus = pow(2, 31) - 1;
    long long initialSeed = modulus / 2;

    int Nvals[3] = { 100, 1000, 10000 };

    // evaluate the generator's randomness using 10^2, 10^3, and 10^4 random numbers
    for (int N : Nvals) {
        MultiplicativeCongruentialGenerator rng(multiplier, modulus, initialSeed);

        std::cout << "Evaluating RNG algorithm at N = " << N << " seeds..." << std::endl;

        // generate N random numbers
        for (int i = 1; i < N; i++) {
            rng.computeNext();
        }

        // perform two-level frequency test
        int K = 10;
        int nbins = 10;
        rng.runTwoLevelFrequencyTest(K, nbins);
        std::cout << std::endl;
    }

    return 0;
}