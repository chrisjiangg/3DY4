/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

// Source code for Fourier-family of functions
#include "dy4.h"
#include "fourier.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>

#define PI 3.14159265358979323846

// Just DFT function (no FFT)
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

// Function to compute the magnitude values in a complex vector
void computeVectorMagnitude(const std::vector<std::complex<real>> &Xf, std::vector<real> &Xmag)
{
	Xmag.clear();
	Xmag.resize(Xf.size(), real(0));
	for (int i = 0; i < (int)Xf.size(); i++) {
		Xmag[i] = std::abs(Xf[i]) / Xf.size();
	}
}

// Add your own code to estimate the PSD
void computePSD(const std::vector<real>& samples, std::vector<real>& PSD) {
    int N = samples.size();
    std::vector<std::complex<real>> freqDomain(N);
    
    // Compute DFT
    DFT(samples, freqDomain);

    // Compute Power Spectral Density
    PSD.resize(N / 2);
    for (int i = 0; i < N / 2; i++) {
        PSD[i] = std::norm(freqDomain[i]) / (N * N);
    }
}
//MAKE SURE TO ADD THIS IS PART OF TAKEHOME
//////////////////////////////////////////////////////////////
// New code as part of benchmarking/testing and the project
//////////////////////////////////////////////////////////////
//copy this template for new functions do not touch any given functions only add our own
void DFT_reference(const std::vector<real> &x, std::vector<std::complex<real>> &Xf) {

	Xf.clear();
	Xf.resize(x.size(), std::complex<real>(0));
	for (int m = 0; m < (int)Xf.size(); m++) {
		for (int k = 0; k < (int)x.size(); k++) {
			std::complex<real> expval(0, -2 * M_PI * (k * m) / x.size());
			Xf[m] +=  + x[k] * std::exp(expval);
		}
	}
}

void DFT_init_bins(const std::vector<real> &x, std::vector<std::complex<real>> &Xf) {

	int N = (int)x.size();
	std::fill(Xf.begin(), Xf.end(), std::complex<real>(0., 0.));
	for (int m = 0; m < N; m++) {
		for (int k = 0; k < N; k++) {
			std::complex<real> expval(0, -2 * M_PI * (k * m) / N);
			Xf[m] += x[k] * std::exp(expval);
		}
	}
}

void generate_DFT_twiddles(const int& N, std::vector<std::complex<real>> &Twiddle1D) {

	Twiddle1D.resize(N);
	for (int k = 0; k < N; k++) {
		std::complex<real> expval(0, -2 * M_PI * k / N);
		Twiddle1D[k] = std::exp(expval);
	}
}

void generate_DFT_matrix(const int& N, std::vector<std::vector<std::complex<real>>> &Twiddle2D) {

	Twiddle2D.resize(N, std::vector<std::complex<real>>(N));
    std::vector<std::complex<real>> Twiddle1D;
	generate_DFT_twiddles(N, Twiddle1D);

	for (int m = 0; m < N; m++) {
		for (int k = 0; k < N; k++) {
			Twiddle2D[m][k] = Twiddle1D[(k * m) % N];
		}
	}
}

//added functions
//algo 2

void DFT_code_motion(const std::vector<real> &x, std::vector<std::complex<real>> &Xf) {
    Xf.clear();
    Xf.resize(x.size(), std::complex<real>(0));
    
    // Precompute angular frequency step
    std::complex<real> ang_freq_step(0, -2 * M_PI / x.size());
    
    for (int m = 0; m < (int)Xf.size(); m++) {
        std::complex<real> intermediate_exp = std::exp(ang_freq_step * (real)m);
        for (int k = 0; k < (int)x.size(); k++) {
            std::complex<real> expval = std::pow(intermediate_exp, k); // Using multiplication instead of exp() calls
            Xf[m] += x[k] * expval;
        }
    }
}

//algo 3
// Optimized DFT using precomputed 1D Twiddle factors
void DFT_twiddle_vector(const std::vector<real> &x, std::vector<std::complex<real>> &Xf) {
    Xf.clear();
    Xf.resize(x.size(), std::complex<real>(0));
    
    // Precompute 1D Twiddle factors
    std::vector<std::complex<real>> Twiddle1D(x.size());
    for (int k = 0; k < (int)x.size(); k++) {
        Twiddle1D[k] = std::exp(std::complex<real>(0, -2 * M_PI * k / x.size()));
    }
    
    for (int m = 0; m < (int)Xf.size(); m++) {
        for (int k = 0; k < (int)x.size(); k++) {
            Xf[m] += x[k] * Twiddle1D[(k * m) % x.size()];
        }
    }
}

//algo 4
// DFT Matrix with Frequency Bins Traversed in the Outer Loop
void DFT_matrix_outer(const std::vector<real> &x, std::vector<std::complex<real>> &Xf) {
    Xf.clear();
    Xf.resize(x.size(), std::complex<real>(0));
    
    int N = (int)x.size();
    std::vector<std::vector<std::complex<real>>> Twiddle2D(N, std::vector<std::complex<real>>(N));
    
    // Precompute full 2D Twiddle Factor Matrix
    for (int m = 0; m < N; m++) {
        for (int k = 0; k < N; k++) {
            Twiddle2D[m][k] = std::exp(std::complex<real>(0, -2 * M_PI * k * m / N));
        }
    }
    
    // Compute the DFT using the precomputed matrix
    for (int m = 0; m < N; m++) {
        for (int k = 0; k < N; k++) {
            Xf[m] += x[k] * Twiddle2D[m][k];
        }
    }
}

//algo 5
// DFT Matrix with Frequency Bins Traversed in the Inner Loop
void DFT_matrix_inner(const std::vector<real> &x, std::vector<std::complex<real>> &Xf) {
    Xf.clear();
    Xf.resize(x.size(), std::complex<real>(0));
    
    int N = (int)x.size();
    std::vector<std::vector<std::complex<real>>> Twiddle2D(N, std::vector<std::complex<real>>(N));
    
    // Precompute full 2D Twiddle Factor Matrix
    for (int m = 0; m < N; m++) {
        for (int k = 0; k < N; k++) {
            Twiddle2D[k][m] = std::exp(std::complex<real>(0, -2 * M_PI * k * m / N));
        }
    }
    
    // Compute the DFT using the precomputed matrix with frequency bins in the inner loop
    for (int k = 0; k < N; k++) {
        for (int m = 0; m < N; m++) {
            Xf[m] += x[k] * Twiddle2D[k][m];
        }
    }
}

//algo 6
// DFT Matrix with Loop Unrolling
void DFT_loop_unrolling(const std::vector<real> &x, std::vector<std::complex<real>> &Xf) {
    Xf.clear();
    Xf.resize(x.size(), std::complex<real>(0));
    
    int N = (int)x.size();
    std::vector<std::vector<std::complex<real>>> Twiddle2D(N, std::vector<std::complex<real>>(N));
    
    // Precompute full 2D Twiddle Factor Matrix
    for (int m = 0; m < N; m++) {
        for (int k = 0; k < N; k++) {
            Twiddle2D[m][k] = std::exp(std::complex<real>(0, -2 * M_PI * k * m / N));
        }
    }
    
    // Compute the DFT using loop unrolling for better performance
    for (int m = 0; m < N; m++) {
        int k = 0;
        for (; k + 3 < N; k += 4) {  // Unrolling loop by a factor of 4
            Xf[m] += x[k] * Twiddle2D[m][k] +
                     x[k+1] * Twiddle2D[m][k+1] +
                     x[k+2] * Twiddle2D[m][k+2] +
                     x[k+3] * Twiddle2D[m][k+3];
        }
        
        // Process remaining elements
        for (; k < N; k++) {
            Xf[m] += x[k] * Twiddle2D[m][k];
        }
    }
}

void estimatePSD(const std::vector<double>& samples, double Fs, std::vector<double>& freq, std::vector<double>& psd_est) {
    int freq_bins = NFFT;
    double df = Fs / freq_bins;
    
    freq.clear();
    for (double num_freq = 0; num_freq < Fs / 2; num_freq += df) {
        freq.push_back(num_freq);
    }
    
    std::vector<double> hann(freq_bins);
    for (int i = 0; i < freq_bins; i++) {
        hann[i] = 0.5 * (1 - std::cos(2 * M_PI * i / (freq_bins - 1)));
    }
    
    int no_segments = static_cast<int>(std::floor(samples.size() / static_cast<double>(freq_bins)));
    std::vector<double> psd_list(freq_bins / 2, 0.0);
    
    for (int k = 0; k < no_segments; k++) {
        std::vector<std::complex<double>> Xf(freq_bins);
        for (int m = 0; m < freq_bins; m++) {
            Xf[m] = std::complex<double>(samples[k * freq_bins + m] * hann[m], 0);
        }
        
        std::vector<std::complex<double>> posSamples(Xf.begin(), Xf.begin() + freq_bins / 2);
        for (int w = 0; w < freq_bins / 2; w++) {
            psd_list[w] += 2.0 * (1 / (Fs * freq_bins / 2)) * std::norm(posSamples[w]);
        }
    }
    
    std::vector<double> psd_seg(freq_bins / 2, 0.0);
    for (int k = 0; k < freq_bins / 2; k++) {
        for (int l = 0; l < no_segments; l++) {
            psd_seg[k] += psd_list[k + l * (freq_bins / 2)];
        }
        psd_seg[k] /= no_segments;
    }
    
    psd_est.clear();
    psd_est.resize(freq_bins / 2, 0.0);
    for (int k = 0; k < freq_bins / 2; k++) {
        psd_est[k] = 10 * std::log10(psd_seg[k]);
    }
}

void matrixPSD(const std::vector<double>& samples, double Fs, std::vector<double>& freq, std::vector<double>& psd_est) {
    int freq_bins = NFFT;
    double df = Fs / freq_bins;
    
    freq.clear();
    for (double num_freq = 0; num_freq < Fs / 2; num_freq += df) {
        freq.push_back(num_freq);
    }
    
    int no_segments = static_cast<int>(std::floor(samples.size() / static_cast<double>(freq_bins)));
    
    std::vector<std::vector<std::complex<double>>> dft_matrix(freq_bins, std::vector<std::complex<double>>(freq_bins));
    for (int m = 0; m < freq_bins; m++) {
        for (int k = 0; k < freq_bins; k++) {
            dft_matrix[m][k] = std::exp(std::complex<double>(0, -2 * M_PI * k * m / freq_bins));
        }
    }
    
    std::vector<double> hann_window(freq_bins);
    for (int i = 0; i < freq_bins; i++) {
        hann_window[i] = 0.5 * (1 - std::cos(2 * M_PI * i / (freq_bins - 1)));
    }
    
    std::vector<std::vector<std::complex<double>>> Xf(no_segments, std::vector<std::complex<double>>(freq_bins, {0, 0}));
    for (int seg = 0; seg < no_segments; seg++) {
        for (int m = 0; m < freq_bins; m++) {
            for (int k = 0; k < freq_bins; k++) {
                Xf[seg][m] += samples[seg * freq_bins + k] * hann_window[k] * dft_matrix[m][k];
            }
        }
    }
    
    psd_est.clear();
    psd_est.resize(freq_bins / 2, 0.0);
    for (int m = 0; m < freq_bins / 2; m++) {
        double sum_power = 0.0;
        for (int seg = 0; seg < no_segments; seg++) {
            sum_power += (1 / ((Fs / 2) * (freq_bins / 2))) * std::norm(Xf[seg][m]);
        }
        psd_est[m] = 10 * std::log10(sum_power / no_segments);
    }
}




