/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

// Source code for Fourier-family of functions
#include "dy4.h"
#include "fourier.h"

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