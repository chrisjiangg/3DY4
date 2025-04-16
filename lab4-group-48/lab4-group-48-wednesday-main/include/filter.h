/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// Add headers as needed
#include <iostream>
#include <vector>

// Declaration of function prototypes
void impulseResponseLPF(real, real, unsigned short int, std::vector<real> &);
void convolveFIR(std::vector<real> &, const std::vector<real> &, const std::vector<real> &);

//////////////////////////////////////////////////////////////
// New code as part of benchmarking/testing and the project
//////////////////////////////////////////////////////////////

void convolveFIR_inefficient(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h);
void convolveFIR_reference(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h);

void convolveFIR_signal(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h);
void convolveFIR_kernel(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h);
void convolveFIR_precomputed(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h);
void convolveFIR_unrolled(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h);

void convolveFIR_kernel_updated(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h);
void convolveFIR_precomputed_updated(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h);
void convolveFIR_unrolled_updated(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h);

//void matrixPSD(const std::vector<real> &samples, int nfft, int Fs, std::vector<real> &psd);
//void estimatePSD(const std::vector<real> &samples, int nfft, int Fs, std::vector<real> &psd);

#endif // DY4_FILTER_H

