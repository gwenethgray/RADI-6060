#include <iostream>
#include <vector>

#ifndef MULTIPLICATIVE_CONGRUENTIAL_GENERATOR
#define MULTIPLICATIVE_CONGRUENTIAL_GENERATOR

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
    double getNextRandomNumber();

    // get current seed
    long long getSeed() const;

    // get current random number
    double getRandomNumber() const;

    // get all stored seeds
    const std::vector<long long>& getSeeds() const;

    // get all stored random numbers
    const std::vector<double>& getRandomNumbers() const;
};

#endif