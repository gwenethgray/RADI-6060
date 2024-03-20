#include "MultiplicativeCongruentialGenerator.hh"
#include <cmath>
#include <vector>

// MCG constructor
MultiplicativeCongruentialGenerator::MultiplicativeCongruentialGenerator(long long multiplier, long long modulus, long long initialSeed)
: I1(multiplier), M(modulus), seed(initialSeed % modulus), randomNumber(static_cast<double>(initialSeed % modulus) / modulus) {
    seeds = { seed };
    randomNumbers = { randomNumber };
}

// compute the next random number in the sequence and store it
void MultiplicativeCongruentialGenerator::computeNext() {
    seed = (I1 * seed) % M;
    seeds.push_back(seed);
    randomNumber = static_cast<double>(seed) / M; // next random number from 0 to 1
    randomNumbers.push_back(randomNumber);
}

// compute and return the next random number
double MultiplicativeCongruentialGenerator::getNextRandomNumber() {
    computeNext();
    return randomNumber;
}

// get current seed
long long MultiplicativeCongruentialGenerator::getSeed() const {
    return seed;
}

// get current random number
double MultiplicativeCongruentialGenerator::getRandomNumber() const {
    return randomNumber;
}

// get all stored seeds
const std::vector<long long>& MultiplicativeCongruentialGenerator::getSeeds() const {
    return seeds;
}

// get all stored random numbers
const std::vector<double>& MultiplicativeCongruentialGenerator::getRandomNumbers() const {
    return randomNumbers;
}