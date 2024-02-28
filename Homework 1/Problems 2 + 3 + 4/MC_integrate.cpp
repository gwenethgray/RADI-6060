#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
using namespace std;

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
    int numberOfHistories = 100;
    double xl = 0;
    double xh = 2;


    // MEAN METHOD

    std::cout << "Computing the integral of f(x) = x^3 between xl = 0 and xh = 2 using the mean method..." << std::endl;

    MultiplicativeCongruentialGenerator rng(multiplier, modulus, initialSeed);

    // store mappings of random numbers to interval from xl to xh
    std::vector<double> x = {xl + rng.getRandomNumber()*(xh - xl)};
    // store values of x^3 at stored x values
    std::vector<double> y = {pow(x[0], 3)};
    // accumulate sum of function evaluations to compute mean
    double sum = y[0];

    // generate 99 more random numbers
    for (int i = 1; i < numberOfHistories; i++) {
        rng.computeNext();
        x.push_back(xl + rng.getRandomNumber()*(xh - xl));
        y.push_back(pow(x[i], 3));
        sum += y[i];
    }

    // compute mean of function evaluations
    double mean = sum / y.size();

    // compute integral
    double integral = (xh - xl) * mean;

    std::cout << "Computed value:\t\t" << integral << std::endl
              << "Expected value:\t\t" << pow(2, 4)/static_cast<double>(4) << std::endl;



    // MEAN METHOD

    std::cout << "Computing the integral of f(x) = x^3 between xl = 0 and xh = 2 using the rejection method..." << std::endl;

    MultiplicativeCongruentialGenerator rng(multiplier, modulus, initialSeed);

    

    return 0;
}