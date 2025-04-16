#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <complex>
#include "functions.h"


void readBlockData(unsigned int num_samples, unsigned int block_id, std::vector<real> &block_data){
	std::vector<char> raw_data(num_samples);
	std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));
	for(int k =0; k <(int)num_samples; k++){
		block_data[k] = real(((unsigned char)raw_data[k]-128)/128.0);
	}
}

int findGCD(int a, int b){
    while(b!=0){
        int out =b;
        b = a%b;
        a=out;
    }
    return a;

}
//split
void split_IQ(const std::vector<real> &iq_data, std::vector<real> &i_data, std::vector<real> &q_data){
    for (size_t i = 0; i < iq_data.size(); i+=2) {
            i_data.push_back(iq_data[i]);
            q_data.push_back(iq_data[i+1]);
    }
}
void upsample(std::vector<real> &data,int upsample){
    std::vector<real> data_out(data.size()*upsample,0.0);
    for(int i=0;i<data_out.size();i++){
        //std::cerr<<"i: "<<i<<" upsample: "<<upsample<<" : "<<i%upsample<<" : "<< data[i/2]<<std::endl;
        if(i%upsample==0){
            data_out[i] = data[i/upsample];
        }else{
            data_out[i] = 0;
        }
    }
    data.resize(data_out.size(),0.0);
    data = data_out;
}

void downsample(std::vector<real> &data,int downsample){
    std::vector<real> data_out(data.size()/downsample,0.0);
    for(int i=0;i<data_out.size();i+=1){
        data_out[i] = data[i*downsample];
    }
    data.resize(data_out.size(),0.0);
    data = data_out;
}
//downsample
void conv_down(std::vector<real> &state, std::vector<real> data, const std::vector<real> coeff, std::vector<real> &y, int down_sample)
{
    y.clear(); 
    y.resize((data.size()*1)/down_sample,0.0);

    int state_s = state.size();
    int data_s = data.size();
    int coeff_s = coeff.size();

    std::vector<real> x_ext(state_s + data_s,0.0);

    for(int i=0; i<state_s; i++){
        x_ext[i] = state[i];
    }

    for(int i=0; i<data_s; i++){
        x_ext[i + state_s] = data[i];
    }

    for (size_t n = 0; n < y.size(); ++n) {
        for (size_t k = 0; k < coeff.size(); ++k) {
            if (n >= k && n - k < x_ext.size()) {
                y[n] += coeff[k] * x_ext[n - k];
            }
        }
    }
    
    for(int k =0; k<(int)state_s;k++){
        state[k] = x_ext[x_ext.size() - state.size() + k];
    }

    //downsample(y, down_sample);
}

void conv_ds_slow(std::vector<real> &y, const std::vector<real> x, const std::vector<real> h, int ds){
    convolveFIR(y,x,h);
    for (int i=0; i< int(y.size()/ds); i++){
        y[i] = y[i*ds];
    }
    y.resize(int(y.size()/ds));
    y.shrink_to_fit();
}

void my_own_conv_state(std::vector<real> &state, std::vector<real> &data, const std::vector<real> &coeff, std::vector<real> &y, int down_sample, int up_sample)
{
    y.clear(); 
    
    y.resize(data.size()*up_sample/down_sample,0.0);

    for(int i = 0;  i < (int)y.size(); i++){
        int phase = (int)((i * down_sample) % up_sample);
        y[i] = 0;
        for(int j=phase; j<(int)coeff.size(); j+=up_sample){
            int k = (int)((i*down_sample-j)/(float)up_sample);
            if(k >= 0){
                y[i] += coeff[j] * data[k];
            }
            else{
                y[i]+= coeff[j]*state[state.size() + k];
            }
        }
    }
    for(int k =0; k<(int)state.size();k++){
        state[k]=data[data.size()-state.size()+k];
    }
}

void convolution(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h, int numBlocks) {
	// Allocate memory for the output (filtered) data
	y.clear();
	y.resize(0, 0.0);

	// The rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	std::vector<real> initial_state(h.size()-1);
	std::vector<real> currentBlock, extendedBlock, blockOutput;
	size_t blockOutputSize, extendedBlockSize, currentBlockSize, blockSize = ceil(x.size() / (float)numBlocks);
	double sum;
	int  i = 0;
	bool end = 0;

	while(i<x.size() || !end){
        currentBlockSize = std::min(blockSize, x.size() - i);
		//currentBlockSize = (currentBlockSize == 0 ) ? 1 : currentBlockSize;
		end = currentBlockSize ==0;
		currentBlock.clear();
		currentBlock.resize(currentBlockSize);
		for (int j=i; j<i+currentBlockSize; j++) currentBlock[j-i] = x[j];

		extendedBlock.clear();
		extendedBlock.resize(currentBlockSize + initial_state.size());
		for (int j=0; j<extendedBlock.size(); j++) ( j < initial_state.size()) ? extendedBlock[j] = initial_state[j] : extendedBlock[j] = currentBlock[j-initial_state.size()];
		extendedBlockSize = currentBlockSize + initial_state.size();
		
		//Output before
		/*
		std::cout<<"\nInitial State: "<<std::endl;
		for(int j=0; j<initial_state.size();j++) std::cout << initial_state[j]<<", ";
		std::cout<<"\nCurrent Block: "<<std::endl;
		for(int j=0; j<currentBlock.size();j++) std::cout << currentBlock[j]<<", ";
		std::cout<<"\nExtended Block: "<<std::endl;
		for(int j=0; j<extendedBlock.size();j++) std::cout << extendedBlock[j]<<", ";
		std::cout<<std::endl;
		*/

		//Convolution
		blockOutputSize = (i+currentBlockSize >=x.size()) ? extendedBlockSize: currentBlockSize;
		blockOutput.resize(blockOutputSize);
		//std::cout<<blockOutputSize<<std::endl;
		for (int j=0; j<blockOutputSize; j++){
			sum = 0.0;
			for (int k = 0; k<h.size(); k++){
				if(j + k < extendedBlockSize){
					sum += h[h.size()-1-k]*extendedBlock[j+k];
					//std::cout<<h[h.size()-1-k]<<" * " << extendedBlock[j+k] << " = " << h[h.size()-1-k] * extendedBlock[j+k] << std::endl;
				}
			}
			//std::cout<<"sum: "<<sum<<std::endl;
			blockOutput[j] = sum;
		}
		
		for(int j = 0; j<h.size()-1; j++){
			initial_state[j] = extendedBlock[extendedBlockSize - h.size() + 1 + j];
		}
		if(i>=x.size()){
			//std::cout <<i<<" : "<< x.size() << std::endl;
			y.insert(y.end(), blockOutput.begin(),blockOutput.begin());
		}else{
			//std::cout <<i<<" : "<< x.size() << std::endl;
			y.insert(y.end(), blockOutput.begin(), blockOutput.end());
		}
		i+=currentBlockSize;
		
	}
}

// Function to perform FM demodulation using arctangent method
void my_own_fmDemodArctan(const std::vector<real>& I, const std::vector<real>& Q, std::vector<real>& fm_demod, real& prev_I_state, real& prev_Q_state) {
    
    real dQ, dI;

    for (int k = 0; k < (int)(I.size()); k++) {
        if (k == 0) {
            dQ = Q[k] - prev_Q_state;
            dI = I[k] - prev_I_state;
        } else {
            dQ = Q[k] - Q[k - 1];
            dI = I[k] - I[k - 1];
        }
        if ((pow(I[k], 2) + pow(Q[k], 2)) == 0) {
            fm_demod[k] = 0;
        } else {
            fm_demod[k] = (1 / (pow(I[k], 2) + pow(Q[k], 2))) * ((I[k] * dQ) - (Q[k] * dI));
        }
    }
    prev_I_state = I[I.size() - 1];
    prev_Q_state = Q[Q.size() - 1];
}

void estimatePSD(const std::vector<real> &samples, unsigned int nfft, real Fs, std::vector<real> &freq, std::vector<real> &psd_est)
{
    unsigned int freq_bins = nfft;
    real df = 1.0*Fs / freq_bins;
    freq.clear();

    for(real i = 0; i < Fs/2; i+=df){
        freq.push_back(i);
    }

    std::vector<double> hann(freq_bins);
    for (unsigned int i = 0; i < freq_bins; ++i) {
        hann[i] = 0.5 * (1 - std::cos(2 * PI * i / (freq_bins - 1)));
    }

    unsigned int num_segments = std::floor(samples.size() / freq_bins);
    std::vector<double> psd_list;
    
    for (unsigned int k = 0; k < num_segments; k++) {
        std::vector<real> segment(samples.begin() + k * freq_bins, samples.begin() + (k + 1) * freq_bins);
        for (unsigned int i = 0; i < freq_bins; i++) {
            segment[i] *= hann[i];
        }
        
        std::vector<std::complex<real>> Xf;
        DFT(segment, Xf);
        Xf.resize(freq_bins/2);

        std::vector<double> psd_seg(freq_bins / 2, 0.0);
        double scale_factor = 1.0 / (Fs * (freq_bins / 2));

        for (unsigned int i = 0; i < freq_bins/2; ++i) {
            psd_seg[i] = scale_factor * std::norm(Xf[i]); 
        }
        for (size_t i = 0; i < psd_seg.size(); i++) {
            psd_seg[i] *= 2;
        }
        psd_list.insert(psd_list.end(), psd_seg.begin(), psd_seg.end());
    }
        
    unsigned int half_bins = freq_bins / 2;
    std::vector<double> psd_seg(half_bins, 0.0);

    for (unsigned int k = 0; k < half_bins; k++) {
        for (unsigned int l = 0; l < num_segments; l++) {
            psd_seg[k] += psd_list[k + l * half_bins];
        }
        psd_seg[k] /= num_segments; 
    }

    psd_est.resize(half_bins);
    for (unsigned int k = 0; k < half_bins; k++) {
        psd_est[k] = 10 * std::log10(psd_seg[k]); 
    }

}


void impulseResponseBPF(real Fs, real Fc1, real Fc2, int num_taps, std::vector<real> &h) {
    // Clear and resize output vector
    h.clear();
    h.resize(num_taps, 0.0);

    // Normalize frequencies to Nyquist frequency (Fs / 2)
    real normFc1 = Fc1 / (Fs / 2.0);
    real normFc2 = Fc2 / (Fs / 2.0);
    int middle = (num_taps - 1) / 2;
    real fmid = (Fc1 + Fc2)/2;
    real scale = 0;
    int n;

    for (int i = 0; i < num_taps; ++i) {
        n = static_cast<int>(i - middle);
        if (n==0) {
            // Special case to avoid division by 0
            h[i] = normFc2 - normFc1;
        } else {
            real sinc1 = sin(PI * normFc1 * n) / (PI * n);
            real sinc2 = sin(PI * normFc2 * n) / (PI * n);
            h[i] = sinc2 - sinc1;
        }

        // Apply Hann window
        h[i] *= 0.5 * (1.0 - cos(2.0 * PI * i / (num_taps - 1)));
        scale = scale + h[i] * cos(2.0 * PI * n * fmid / Fs );
    }
    for(int i=0; i<num_taps;i++){
        h[i] = h[i]/scale;
    }
}

void DelayMono(const std::vector<real> fm_demod_block, std::vector<real> &delay_state, std::vector<real> &mono_delayed) {
    int delay_size = delay_state.size();
    int total_size = fm_demod_block.size();

    mono_delayed.resize(total_size);

    // Step 1: Use old samples from delay_state for the start of mono_delayed
    for (int i = 0; i < delay_size; ++i) {
        mono_delayed[i] = delay_state[i];
    }

    // Step 2: Fill the rest of mono_delayed with current fm_demod samples
    for (int i = delay_size; i < total_size; ++i) {
        mono_delayed[i] = fm_demod_block[i - delay_size];
    }

    // Step 3: Update delay_state with the last 'delay_size' samples from current fm_demod
    for (int i = 0; i < delay_size; ++i) {
        delay_state[i] = fm_demod_block[total_size - delay_size + i];
    }
}


void fmPll(std::vector<real> &pllIn, real freq, real fs, real ncoscale, real normBandwidth, std::vector<real> &ncoOut, real &integrator, real &phaseEst, real &feedbackI, real &feedbackQ, real &ncoOut1, real &trigOffset){
    real Cp = 2.666;
    real Ci = 3.555;

    // gain for the porportional term
    real Kp = normBandwidth * Cp;
    // gain for the intergrator term
    real Ki = (normBandwidth * normBandwidth) * Ci;

    // output array for the NCO
    ncoOut.resize(pllIn.size()+1,0);

    // Initialize internal state
    real integrator1 = integrator;
    real phaseEst1 = phaseEst;
    real feedbackI1 = feedbackI;
    real feedbackQ1 = feedbackQ;
    ncoOut[0] = ncoOut1;
    real trigOffset1 = trigOffset;

    for (int k = 0; k < (int)pllIn.size(); k++) {
        // Phase detector
        real errorI = pllIn[k] * (+feedbackI1);
        real errorQ = pllIn[k] * (-feedbackQ1);

        // Phase error detection
        real errorD = std::atan2(errorQ, errorI);

        // Loop filter
        integrator1 = integrator1 + Ki * errorD;
        // Update phaseEst
        phaseEst1 = phaseEst1 + Kp * errorD + integrator1;

        // NCO (internal oscillator)
        trigOffset1 += 1;
        real trigArg = 2 * PI * (freq/fs) * trigOffset1 + phaseEst1;
        feedbackI1 = std::cos(trigArg);
        feedbackQ1 = std::sin(trigArg);
        ncoOut[k+1] = std::cos(trigArg * ncoscale);
    }
    integrator = integrator1;
    phaseEst = phaseEst1;
    feedbackI = feedbackI1;
    feedbackQ = feedbackQ1;
    ncoOut1 = ncoOut[ncoOut.size()-1];
    trigOffset = trigOffset1;
}
