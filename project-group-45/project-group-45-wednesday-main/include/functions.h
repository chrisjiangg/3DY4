#ifndef DY4_FUNCTIONS_H
#define DY4_FUNCTIONS_H

#include "dy4.h"

// Add headers as needed
#include <iostream>
#include <vector>

// Declaration of function prototypes
void readBlockData(unsigned int num_samples, unsigned int block_id, std::vector<real> &block_data);
void split_IQ(const std::vector<real> &bin_data, std::vector<real> &i_data, std::vector<real> &q_data);
void my_own_conv_state(std::vector<real> &state, std::vector<real> &data, const std::vector<real> &coeff, std::vector<real> &y, int down_sample, int up_sample);
void my_own_fmDemodArctan(const std::vector<real>& I, const std::vector<real>& Q, std::vector<real>& fm_demod, real& prev_I_state, real& prev_Q_state);
void estimatePSD(const std::vector<real> &samples, unsigned int nfft, real Fs, std::vector<real> &freq, std::vector<real> &psd_est);
void upsample(std::vector<real> &data,int upsample);
void downsample(std::vector<real> &data,int downsample);
void conv_down(std::vector<real> &state, std::vector<real> data, const std::vector<real> coeff, std::vector<real> &y, int down_sample);
int findGCD(int a, int b);
void convolution(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h, int numBlocks);
void conv_ds_slow(std::vector<real> &y, const std::vector<real> x, const std::vector<real> h, int ds);
void impulseResponseBPF(real Fs, real Fc1, real Fc2, int num_taps, std::vector<real> &h);
void DelayMono(const std::vector<real> fm_demod_block, std::vector<real> &delay_state, std::vector<real> &mono_delayed);

void fmPll(std::vector<real> &pllIn, real freq, real fs, real ncoscale, real normBandwidth,  std::vector<real> &ncoOut, real &integrator, real &phaseEst, real &feedbackI, real &feedbackQ, real &ncoOut1, real &trigOffset);
//void conv_resampler(std::vector<real> &state, std::vector<real> &data, const std::vector<real> &coeff, std::vector<real> &y,int upsample, int downsample);

//////////////////////////////////////////////////////////////
// New code as part of benchmarking/testing and the project
//////////////////////////////////////////////////////////////

#endif // DY4_FUNCTIONS_H
