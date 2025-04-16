/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>

#define PI 3.14159265358979323846

// Function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(real Fs, real Fc, unsigned short int num_taps, std::vector<real> &h)
{
	// Bring your own functionality
	h.clear();
    	h.resize(num_taps, 0.0);

    	real cutoff = Fc/(Fs/2);
    	for(int i=0; i<num_taps; i++){
        	if(i == (num_taps-1)/2){
            		h[i] = cutoff;
        	}
        	else{
            		h[i] = cutoff * (sin(PI*cutoff*(i-(num_taps-1)/2)))/(PI*cutoff*(i-(num_taps-1)/2));
        	}
        	h[i] = h[i]*pow(sin((i*PI)/num_taps), 2);
    	}
}

// Function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
void convolveFIR(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h)
{
	// Bring your own functionality
	
	y.clear();
    	y.resize(x.size() + h.size() - 1, 0.0);

    	for(real i=0; i < y.size(); i++){
        	for(real j=0; j < h.size(); j++){
            		if(((i-j)>=0) && ((i-j)< x.size())){
                		y[i] += h[j] * x[i-j]; 
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

// Algorithm 7: Single-pass convolution (iterate over signal x)
void convolveFIR_signal(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h) {
    y.clear();
    y.resize(int(x.size() + h.size() - 1), 0.0);
    for (int n = 0; n < (int)y.size(); n++) {
        for (int k = 0; k < (int)x.size(); k++) {
            if ((n - k >= 0) && (n - k) < (int)h.size())
                y[n] += x[k] * h[n - k];
        }
    }
}

// Algorithm 8: Single-pass convolution (iterate over kernel h)
void convolveFIR_kernel(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h) {
    y.clear();
    y.resize(int(x.size() + h.size() - 1), 0.0);
    for (int n = 0; n < (int)y.size(); n++) {
        for (int k = 0; k < (int)h.size(); k++) {
            if ((n - k >= 0) && (n - k) < (int)x.size())
                y[n] += h[k] * x[n - k];
        }
    }
}

// Algorithm 9: Convolution with precomputed valid range
void convolveFIR_precomputed(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h) {
    y.clear();
    y.resize(int(x.size() + h.size() - 1), 0.0);
    for (int n = 0; n < (int)y.size(); n++) {
        int kmin = std::max(0, n - (int)x.size() + 1);
        int kmax = std::min((int)h.size(), n + 1);
        for (int k = kmin; k < kmax; k++) {
            y[n] += h[k] * x[n - k];
        }
    }
}

// Algorithm 10: Convolution with precomputed range and loop unrolling
void convolveFIR_unrolled(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h) {
    y.clear();
    y.resize(int(x.size() + h.size() - 1), 0.0);
    for (int n = 0; n < (int)y.size(); n++) {
        int kmin = std::max(0, n - (int)x.size() + 1);
        int kmax = std::min((int)h.size(), n + 1);
        int unrolled_limit = kmax - ((kmax - kmin) % 4);
        for (int k = kmin; k < unrolled_limit; k += 4) {
            y[n] += h[k] * x[n - k] + h[k + 1] * x[n - (k + 1)] +
                     h[k + 2] * x[n - (k + 2)] + h[k + 3] * x[n - (k + 3)];
        }
        for (int k = unrolled_limit; k < kmax; k++) {
            y[n] += h[k] * x[n - k];
        }
    }
}

// Algorithm 8 (Updated): Single-pass convolution with kernel h in outer loop
void convolveFIR_kernel_updated(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h) {
    y.clear();
    y.resize(int(x.size() + h.size() - 1), 0.0);
    for (int k = 0; k < (int)h.size(); k++) { // Outer loop over kernel h
        for (int n = 0; n < (int)y.size(); n++) {
            if ((n - k >= 0) && (n - k) < (int)x.size())
                y[n] += h[k] * x[n - k];
        }
    }
}

// Algorithm 9 (Updated): Precomputed valid range with loop-reordered convolution
void convolveFIR_precomputed_updated(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h) {
    y.clear();
    y.resize(int(x.size() + h.size() - 1), 0.0);
    for (int k = 0; k < (int)h.size(); k++) { // Outer loop over kernel h
        for (int n = std::max(0, k); n < std::min((int)y.size(), (int)x.size() + k); n++) {
            y[n] += h[k] * x[n - k];
        }
    }
}

// Algorithm 10 (Updated): Loop unrolling applied to Algorithm 9
void convolveFIR_unrolled_updated(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h) {
    y.clear();
    y.resize(int(x.size() + h.size() - 1), 0.0);
    for (int k = 0; k < (int)h.size(); k++) { // Outer loop over kernel h
        int nmin = std::max(0, k);
        int nmax = std::min((int)y.size(), (int)x.size() + k);
        int unrolled_limit = nmax - ((nmax - nmin) % 4);
        for (int n = nmin; n < unrolled_limit; n += 4) {
            y[n] += h[k] * x[n - k];
            y[n+1] += h[k] * x[n+1 - k];
            y[n+2] += h[k] * x[n+2 - k];
            y[n+3] += h[k] * x[n+3 - k];
        }
        for (int n = unrolled_limit; n < nmax; n++) {
            y[n] += h[k] * x[n - k];
        }
    }
}


