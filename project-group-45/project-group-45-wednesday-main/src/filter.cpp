/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <complex>

// Function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(real Fs, real Fc, unsigned short int num_taps, std::vector<real> &h)
{
	// Bring your own functionality
	// Allocate memory for the impulse response
	h.clear();
	h.resize(num_taps, 0.0);

	// The rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	double normCutoff = Fc / (Fs/2);
	int middleIndex = (num_taps-1)/2;  

	for (int i = 0; i<num_taps; i++){
		if(i == middleIndex){
			h[i] = normCutoff;
		}else{
			h[i] = normCutoff*(sin(PI*normCutoff*(i - middleIndex)) / (PI*normCutoff*(i-middleIndex)));
		}
		h[i] = h[i] * (0.5 - cos(2*PI*i / (num_taps - 1))/2);
	}
}

// Function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF_amp(real Fs, real Fc, unsigned short int num_taps, std::vector<real> &h, int amp)
{
	// Bring your own functionality
	// Allocate memory for the impulse response
	h.clear();
	h.resize(num_taps, 0.0);

	// The rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	double normCutoff = Fc / (Fs/2);
	int middleIndex = (num_taps-1)/2;  

	for (int i = 0; i<num_taps; i++){
		if(i == middleIndex){
			h[i] = normCutoff;
		}else{
			h[i] = normCutoff*(sin(PI*normCutoff*(i - middleIndex)) / (PI*normCutoff*(i-middleIndex)));
		}
		h[i] = amp *h[i] * (0.5 - cos(2*PI*i / (num_taps - 1))/2);
	}
}

// Function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
void convolveFIR(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h)
{
	// Bring your own functionality
	// Allocate memory for the output (filtered) data
	y.clear();
	y.resize(x.size() + h.size() - 1, 0.0);

	// The rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	for (size_t n = 0; n < y.size(); ++n) {
        for (size_t k = 0; k < h.size(); ++k) {
            if (n >= k && n - k < x.size()) {
                y[n] += h[k] * x[n - k];
            }
        }
    }
}

//////////////////////////////////////////////////////////////
// New code as part of benchmarking/testing and the project
//////////////////////////////////////////////////////////////

void convolveFIR_inefficient(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h) {
    y.clear();
    y.resize(int(x.size() + h.size()-1), 0.0);
    for (auto n = 0; n < (int)y.size(); n++) {
        for (auto k = 0; k < (int)x.size(); k++) {
            if ((n - k >= 0) && (n - k) < (int)h.size())
                y[n] += x[k] * h[n - k];
        }
    }
}

void convolveFIR_reference(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h) {
    y.clear();
    y.resize(int(x.size() + h.size()-1), 0.0);

    for (auto n = 0; n < (int)y.size(); n++) {
        for (auto k = 0; k < (int)h.size(); k++) {
            if ((n - k >= 0) && (n - k) < (int)x.size())
                y[n] += h[k] * x[n - k];
        }
    }
}
