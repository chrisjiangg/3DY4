/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <chrono>

#define PI 3.14159265358979323846

// Macro for conditional compilation
#ifdef DOUBLE
	typedef double real;
#else
	typedef float real;
#endif

// Epsilon for unit test comparison
const real EPSILON = 1e-3;

// Input x is qualified as const in order not to be changed within the function
// Although C++ is case sensitive, it is "safer" to use Xf instead of X to avoid confusion
void DFT(const std::vector<real> &x, std::vector<std::complex<real>> &Xf) {
	// Remove all the elements in the vector container before resizing
	Xf.clear();
	// Allocate space for the vector holding the frequency bins
	// Initialize all the values in the frequency vector to complex 0
	// Only new elements are initialized (hence it is safe to clear before resizing)
	Xf.resize(x.size(), std::complex<real>(0));
	for (int m = 0; m < (int)Xf.size(); m++) {
		for (int k = 0; k < (int)x.size(); k++) {
			// Below is the declaration of a complex variable
			// that is initialized through its constructor (real, imag)
			std::complex<real> expval(0, -2 * PI * (k * m) / x.size());
			// Accumulate partial products to frequency bin "m"
			Xf[m] += x[k] * std::exp(expval);
		}
	}
}

void IDFT(const std::vector<std::complex<real>> &Xf, std::vector<real> &x) {
    int N = Xf.size();
    x.clear();
    x.resize(N, 0.0);
    
    for (int n = 0; n < N; n++) {
        std::complex<real> sum(0, 0);
        for (int m = 0; m < N; m++) {
            std::complex<real> expval(0, 2 * PI * (n * m) / N);
            sum += Xf[m] * std::exp(expval);
        }
        x[n] = std::real(sum) / N; 
    }
}

// Function to generate N random real values
// x is the input/output vector
// N is the number of samples
// Random values are between -max and max
// Precision should be capped to 3 (in the context of this type of experiments)
void generateRandomSamples(std::vector<real> &x, unsigned int N, unsigned short int max, unsigned char precision) {
	// Allocate space for the vector with random values
	x.clear();
	x.resize(N);
	int int_random_max = 2 * (max * int(pow(10, precision)));
	for (int i = 0; i < (int)x.size(); i++) {
		x[i] = int(std::rand() % int_random_max);
		x[i] = (x[i] - (int_random_max / 2)) / pow(10, precision);
	}
}

// Function to print a real vector
void printRealVector(const std::vector<real> &x) {
	std::cout << "Printing real vector of size " << x.size() << "\n";
	for (int i = 0; i < (int)x.size(); i++) {
		std::cout << x[i] << " ";
	}
	std::cout << "\n";
}

// Function to print a complex vector
void printComplexVector(const std::vector<std::complex<real>> &X) {
	std::cout << "Printing complex vector of size " << X.size() << "\n";
	for (int i = 0; i < (int)X.size(); i++) {
		std::cout << X[i] << " ";
	}
	std::cout << "\n";
}

// Function to check if two real vectors are approximately equal
bool areVectorsClose(const std::vector<real> &v1, const std::vector<real> &v2, const real &tolerance) {
	if (v1.size() != v2.size()) return false;
	for (int i = 0; i < (int)v1.size(); i++) {
		if (std::abs(v1[i] - v2[i]) > tolerance) return false;
	}
	return true;
}

// Function to compute the maximum absolute difference between two vectors
real maxDifference(const std::vector<real> &v1, const std::vector<real> &v2) {
	real max_diff = 0.0;
	for (int i = 0; i < (int)v1.size(); i++) {
		real diff = std::abs(v1[i] - v2[i]);
		max_diff = std::max(max_diff, diff);
	}
	return max_diff;
}


// Unit test for Fourier and Inverse Fourier transform
void fourierUnitTest(unsigned int N = 32, unsigned short int max = 10, unsigned char precision = 3) {
	// Generate random sequence
	std::vector<real> x;
	generateRandomSamples(x, N, max, precision);

	// Perform DFT and IDFT
	std::vector<std::complex<real>> Xf;
	DFT(x, Xf);
	std::vector<real> rx(N, 0.0);
	// NOTE: you will have implement IDFT to generate the reconstructed the real signal rx;

    auto start_time = std::chrono::high_resolution_clock::now();
    IDFT(Xf, rx);
    auto stop_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> IDFT_run_time = stop_time - start_time;

	// Check if the original and reconstructed signals match
	if (!areVectorsClose(x, rx, EPSILON)) {
		std::cout << "Comparison between original signal and the inverse Fourier transform of its frequency domain representation fails\n";
		std::cout << "Original signal: ";
		for (int i = 0; i < (int)x.size(); i++) std::cout << x[i] << " ";
		std::cout << "\nReconstructed signal: ";
		for (int i = 0; i < (int)rx.size(); i++) std::cout << rx[i] << " ";
		std::cout << "\nMaximum difference: " << maxDifference(x, rx) << "\n";
	} else {
		std::cout << "Unit test for Fourier/Inverse Fourier transform passed.\n";
	}
}

void measureExecutionTimes() {
    std::vector<unsigned int> sample_sizes = {1024, 2048, 4096, 8192, 16384}; //samples for testing
    for (unsigned int N : sample_sizes) {
        std::vector<real> x;
        generateRandomSamples(x, N, 10, 3);

        std::vector<std::complex<real>> Xf;
        auto start_time_dft = std::chrono::high_resolution_clock::now();
        DFT(x, Xf);
        auto stop_time_dft = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> DFT_run_time = stop_time_dft - start_time_dft;

        std::vector<real> rx(N, 0.0);
        auto start_time_idft = std::chrono::high_resolution_clock::now();
        IDFT(Xf, rx);
        auto stop_time_idft = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> IDFT_run_time = stop_time_idft - start_time_idft;

        std::cout << "N = " << N << ": DFT ran for " << DFT_run_time.count() << " ms, "
                  << "IDFT ran for " << IDFT_run_time.count() << " ms\n";
    }
}

int main() {
	// Initialize the seed for the random number generator using current time
	unsigned int seed = std::time(0x0);
	std::srand(seed);
	std::cout << "Starting from seed " << std::hex << seed << std::dec << "\n";

	// By default, we use floats; if you wish to use double precision, compile as shown below:
	// g++ -DDOUBLE -o program_name source_file.cpp
	// -DDOUBLE will define the DOUBLE macro for conditional compilation
	std::cout << "Working with reals on " << sizeof(real) << " bytes" << std::endl;

	// Declare a vector of real values; no memory is allocated at this time
	std::vector<real> x;
	// Generate 32 samples between -10 and 10; for extra flexibility
	// the last argument gives precision in terms of fraction digits
	// Note: Memory for x is allocated within the function called below
	generateRandomSamples(x, 32, 10, 2);

	// Declare a vector of complex values; no memory is allocated for it
	std::vector<std::complex<real>> Xf;
	// Perform DFT of x to produce Xf
	// Measure the execution time using the "chrono" class
	auto start_time = std::chrono::high_resolution_clock::now();
	DFT(x, Xf);
	// The "auto" keyword is used for type deduction (or inference) in C++
	// on some platforms you might need to provide -std=c++11 as an argument to g++
	auto stop_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> DFT_run_time = stop_time - start_time;
	std::cout << "DFT ran for " << DFT_run_time.count() << " milliseconds" << "\n";

	// Run the unit test after implementing IDFT
	// fourierUnitTest(x.size());
	measureExecutionTimes();
	// Finished
	return 0;
}
