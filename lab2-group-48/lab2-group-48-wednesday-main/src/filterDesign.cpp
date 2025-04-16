/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <complex>
#include <cmath>

#define PI 3.14159265358979323846

// Macro for conditional compilation
#ifdef DOUBLE
	typedef double real;
#else
	typedef float real;
#endif

// Function for DFT (reused from previous experiment)
void DFT(const std::vector<real> &x, std::vector<std::complex<real>> &Xf) {
	Xf.clear();
	Xf.resize(x.size(), std::complex<real>(0));
	for (int m = 0; m < (int)Xf.size(); m++) {
		for (int k = 0; k < (int)x.size(); k++) {
			std::complex<real> expval(0, -2 * PI * (k * m) / x.size());
			Xf[m] += x[k] * std::exp(expval);
		}
	}
}

// Function to generate a sine with Fs samples per second (scaled by interval duration)
void generateSin(std::vector<real> &t, std::vector<real> &x, real Fs, real interval, real frequency = 7.0, real amplitude = 5.0, real phase = 0.0)
{
	// We do NOT allocate memory space explicitly
	// for the time (t) vector and sample (x) vector
	t.clear();
	x.clear();
	real dt = 1 / Fs;
	for (real i = 0.0; i < interval; i += dt) {
		// Vector size increases when pushing new elements into it
		t.push_back(i);
		x.push_back(amplitude * std::sin(2 * PI * frequency * i + phase));
	}
}

real generateRandomSin(std::vector<real> &t, std::vector<real> &x, real Fs, real interval) {
	//random frequency between 25 and 75 Hz
	real frequency = 25 + static_cast<real>(std::rand()) / (RAND_MAX / (75 - 25));
	//random amplitude between 3 and 6
	real amplitude = 3 + static_cast<real>(std::rand()) / (RAND_MAX / (6 - 3));
	
	t.clear();
	x.clear();
	real dt = 1 / Fs;
	for (real i = 0.0; i < interval; i += dt) {
		t.push_back(i);
		x.push_back(amplitude * std::sin(2 * PI * frequency * i));
	}
	return frequency;
}


// Function to add an array of sines
void addSin(const std::vector<std::vector<real>> &sv, std::vector<real> &added)
{
	// Assumes at least one sine passed to this function AND all input sines are of the same size
	for (int i = 0; i < (int)sv[0].size(); i++) {
		real addval = 0.0;
		// Note: sv.size() returns the number of sines (or rows in 2D repr)
		// sv[0].size() returns the number of samples in a sine (or cols in 2D repr)
		for (int k = 0; k < (int)sv.size(); k++) {
			addval += sv[k][i];
		}
		added.push_back(addval);
	}
}

void multiplySines(const std::vector<real> &x1, const std::vector<real> &x2, std::vector<real> &result) {
	result.clear();
	for (size_t i = 0; i < x1.size(); i++) {
		result.push_back(x1[i] * x2[i]);
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

// Function to record data in a format to be read by GNUplot
void plotaddedSinesSpectrum(const std::vector<real> &t, const std::vector<std::vector<real>> &sv, const std::vector<real> &x, const std::vector<real> &y, const std::vector<std::complex<real>> &Xf, const std::vector<std::complex<real>> &Yf)
{
	// Write data in text format to be parsed by GNUplot
	const std::string filename = "../data/example.dat";
	std::fstream fd;  // file descriptor
	fd.open(filename, std::ios::out);
	fd << "#\tindex\tsine(0)\tsine(1)\tsine(2)\tdata in\tdata out\tspectrum in\tspectrum out\n";
	for (int i = 0; i < (int)t.size(); i++) {
		fd << "\t " << i << "\t";
		for (int k = 0; k < (int)sv.size(); k++) {
			fd << std::fixed << std::setprecision(3) << sv[k][i] << "\t ";
		}
		fd << x[i] << "\t " << y[i] << "\t\t ";
		fd << std::abs(Xf[i]) / Xf.size() << "\t\t " << std::abs(Yf[i]) / Yf.size() << "\n";
	}
	std::cout << "Generated " << filename << " to be used by GNUplot\n";
	fd.close();
}

// Function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(real Fs, real Fc, unsigned short int num_taps, std::vector<real> &h)
{
	// Allocate memory for the impulse response
	h.clear();
	h.resize(num_taps, 0.0);

	// The rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	real norm_freq = Fc / (Fs / 2); //cutoff freq
    int mid = num_taps / 2;

    for (int i = 0; i < num_taps; i++) {
        if (i == mid)
            h[i] = norm_freq;
        else {
            real x = PI * (i - mid);
            h[i] = norm_freq * (std::sin(norm_freq * x) / x);
        }
        h[i] *= (0.5 - 0.5 * std::cos((2 * PI * i) / (num_taps - 1)));
    }
}

// Function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
void convolveFIR(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h)
{
	// Allocate memory for the output (filtered) data
	y.clear();
	y.resize(x.size() + h.size() - 1, 0.0);

	// The rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab

	int N = x.size(); //input signal length
    int M = h.size(); //impulse response length

	for (int n = 0; n < N + M - 1; n++) {
        for (int m = 0; m < M; m++) {
            if (n - m >= 0 && n - m < N) {
                y[n] += h[m] * x[n - m];
            }
        }
    }
}

int main()
{
	// By default, we use floats; if you wish to use double precision, compile as shown below:
	// g++ -DDOUBLE -o program_name source_file.cpp
	// -DDOUBLE will define the DOUBLE macro for conditional compilation
	std::srand(std::time(0));
	std::cout << "Working with reals on " << sizeof(real) << " bytes" << std::endl;

	real Fs = 512.0;                     // samples per second
	real interval = 1.0;                 // number of seconds
	unsigned short int num_taps = 101;   // number of filter taps

	std::vector<real> t, sine1, sine2;
	real freq1 = generateRandomSin(t, sine1, Fs, interval);
	real freq2 = generateRandomSin(t, sine2, Fs, interval);

	std::vector<real> multiplied;
	multiplySines(sine1, sine2, multiplied);

	real Fc = std::max(freq1, freq2);

	// Declare a vector of vectors for multiple sines
	//std::vector<std::vector<real>> sv;
	// Declare time and sine vectors
	//std::vector<real> t, sine;

	// Generate and store the first tone
	// (check the function to understand the order of arguments)
	//generateSin(t, sine, Fs, interval, 10.0, 5.0, 0.0);
	//sv.push_back(sine);
	// Generate and store the second tone
	//generateSin(t, sine, Fs, interval, 60.0, 2.0, 0.0);
	//sv.push_back(sine);
	// Generate and store the third tone
	//generateSin(t, sine, Fs, interval, 80.0, 3.0, 0.0);
	//sv.push_back(sine);

	// Declare the added sine vector and add the three tones
	//std::vector<real> x;
	//addSin(sv, x);

	// Declare a vector of complex values for DFT; no memory is allocated for it
	//std::vector<std::complex<real>> Xf;
	//DFT(x, Xf);

	// Generate the impulse response h and convolve it with the input data x
	// in order to produce the output data y
	std::vector<real> h;
	impulseResponseLPF(Fs, Fc, num_taps, h);

	std::vector<real> filtered;
	convolveFIR(filtered, multiplied, h);

	// Compute DFT of the filtered data
	std::vector<std::complex<real>> Xf, Yf;
	DFT(multiplied, Xf);
	DFT(filtered, Yf);

	// Prepare the data for GNUplot
	std::vector<std::vector<real>> sv = {sine1, sine2};
	plotaddedSinesSpectrum(t, sv, multiplied, filtered, Xf, Yf);

	// There is an example.gnuplot script in the ../data sub-folder that can be used GNUplot
	// using the data that is logged by this program (stored in ../data/example.dat by default)
	// naturally, you can comment the line below once you are comfortable to run GNUplot
	std::cout << "Run: gnuplot -e 'set terminal png size 1024,768' ../data/example.gnuplot > ../data/example.png\n";

	return 0;
}
