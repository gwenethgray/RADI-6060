#include <iostream>
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
};

int main() {
    // parameters for the generator
    long long multiplier = pow(7, 5);
    long long modulus = pow(2, 31) - 1;
    long long initialSeed = modulus / 2;

    MultiplicativeCongruentialGenerator rng(multiplier, modulus, initialSeed);

    // print parameters I1 (multiplier), M (modulus), and K(0) (initialSeed)
    std::cout << "Parameters:\n\tI1 = " << multiplier
              << "\n\tM = " << modulus << "\n\tK(0) = " << initialSeed << std::endl;

    // print table header
    std::cout << "Tabulated Random Numbers" << std::endl;
    std::cout << "n\tKn\t\tRn" << std::endl;

    // generate and print 10 random numbers
    for (int i = 0; i < 50; i++) {
        std::cout << i << "\t" << rng.getSeed() << "\t" << rng.getRandomNumber() << std::endl;
        rng.computeNext();
    }

    return 0;
}