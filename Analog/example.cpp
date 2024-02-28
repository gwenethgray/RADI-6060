#include "MultiplicativeCongruentialGenerator.hh"
#include <iostream>
#include <cmath>
#include <vector>

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