#include "MultiplicativeCongruentialGenerator.hh"
#include <cmath>
#include <vector>
#include <iostream>

int N_histories = 10;



int main() {
    // parameters for the generator
    long long multiplier = pow(7, 5);
    long long modulus = pow(2, 31) - 1;
    long long initialSeed = modulus / 2;

    MultiplicativeCongruentialGenerator rng(multiplier, modulus, initialSeed);

    // initiate histories
    for (int i = 0; i < N_histories; i++) {
        // compute starting coordinates
        rng.getNextRandomNumber();
    }

    return 0;
}
