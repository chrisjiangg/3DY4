/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

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
#include "queue.cpp"
#include <complex>
#include "functions.h"
#include <thread>

//structure for the mode settings
struct ModeSettings{
	real rf_Fs; //RF sample rate
	real if_Fs; //IF sample rate
	real audio_Fs; //audio sample rate
	int sps; //SPS for RDS path mode 0 and 2
};

//map to store the mode settings for our group
std::map<int, ModeSettings> modeSettings = {
	{0, {2400, 240, 48, 14}}, //mode 0
	{1, {1440, 360, 36, 0}}, //mode 1 (no RDS)
	{2, {2400, 240, 44.1, 23}}, //mode 2
	{3, {2304, 192, 44.1, 0}} //mode 3 (no RDS)
};


void RF_frontEnd(std::vector<real> iq_data, real rf_Fs, real if_Fs, real &prev_I_state,real &prev_Q_state,int Ntaps,std::vector<real> &i_sample,std::vector<real> &q_sample, std::vector<real> rf_coeff, std::vector<real> &fm_demod, std::vector<real> &fm_demod_block){
    std::vector<real> i_filter, q_filter;
    std::vector<real> i_data, q_data;
    split_IQ(iq_data, i_data, q_data);

    int gcd = findGCD(rf_Fs, if_Fs);
    int upsample_factor = if_Fs / gcd;
    int downsample_factor = rf_Fs / gcd;

    my_own_conv_state(i_sample, i_data, rf_coeff, i_filter, downsample_factor, upsample_factor);
    my_own_conv_state(q_sample, q_data, rf_coeff, q_filter, downsample_factor, upsample_factor);


    fm_demod_block.resize(i_filter.size(), 0.0f);
    my_own_fmDemodArctan(i_filter, q_filter, fm_demod_block, prev_I_state, prev_Q_state);
}

void mono(std::vector<real> &fm_demod, real if_Fs, real audio_Fs, std::vector<real> &mono_sample, std::vector<short int> &fmono_audio,int upsample_factor,int downsample_factor,int Ntaps_mono, std::vector<real> &mono_coeff,std::vector<real> &mono_audio) {
    
    mono_audio.clear();
	real s=0;
	
    my_own_conv_state(mono_sample, fm_demod, mono_coeff, mono_audio, downsample_factor, upsample_factor);

    fmono_audio.resize(mono_audio.size());
    for (long unsigned int i = 0; i < mono_audio.size(); i++) {
        if (std::isnan(mono_audio[i])) {
            fmono_audio[i] = 0; // Handle NaN values
        } else {
			fmono_audio[i] = static_cast<short int>(mono_audio[i]*16384);
        }
    }
	fwrite(&fmono_audio[0], sizeof(short int), fmono_audio.size(), stdout);
}

void stereo (std::vector<real> fm_demod_block,
	real if_Fs, 
	real audio_Fs,
	int downsample,
	int upsample,        
	std::vector<real> &delay_state,
	std::vector<real> & mono_state,
	std::vector<real> & mono_coeff,
	std::vector<real> & mono_delayed,
	std::vector<real> &pilot_coeff,
	std::vector<real> &pilot_state,
	std::vector<real> &stereo_coeff,
	std::vector<real> &stereo_state,
	std::vector<real> &mixed_coeff,
	std::vector<real> &mixed_state,
	std::vector<real> &mixed_out,
	std::vector<short int> &fstereo_audio,
	real &integrator, real &phaseEst, real &feedbackI, real &feedbackQ, real &ncoOut1, real &trigOffset
	){
	// Delay Mono Path
	mono_delayed.resize(fm_demod_block.size(), 0.0); //could this be 100?
	DelayMono(fm_demod_block, delay_state, mono_delayed);

	// Resample Mono Path (fm_demod delayed -> mono_audio)
	std::vector<real> mono_audio;
	my_own_conv_state(mono_state, mono_delayed, mono_coeff, mono_audio, downsample, upsample);

	// // Extract Pilot Tone (19 kHz bandpass filter)
	std::vector<real> pilot;
	my_own_conv_state(pilot_state, fm_demod_block, pilot_coeff, pilot, 1, 1);// -> pilot

	// PLL to recover stereo carrier (38 kHz from 19 kHz pilot)
	std::vector<real> stereo_carrier;
	real pll_freq = 19000;  // Pilot frequency
	fmPll(pilot, pll_freq, if_Fs, 2.0, 0.01, stereo_carrier, integrator, phaseEst, feedbackI, feedbackQ, ncoOut1, trigOffset);//-> stereo_carrier
	

	// Extract stereo frequency(22k - 54k)
	std::vector<real> stereo;
	my_own_conv_state(stereo_state, fm_demod_block, stereo_coeff, stereo, 1, 1);// -> stereo

	//MIXER
	std::vector<real> mixed_s_c;
	for (int i =0; i < stereo.size(); i++){
		mixed_s_c.push_back(stereo[i]*stereo_carrier[i]*2);
	}

	std::vector<real> mixed_filtered;
	my_own_conv_state(mixed_state, mixed_s_c, mixed_coeff, mixed_filtered, downsample, upsample);// -> mixed_filtered

	//stereo combiner
	mixed_out.clear();
	for(int i=0;i<mixed_filtered.size();i++){
		mixed_out.push_back((mono_audio[i] + mixed_filtered[i])/2); 
		mixed_out.push_back((mono_audio[i] - mixed_filtered[i])/2);
	}

	//scaling to PCM
	fstereo_audio.resize(mixed_out.size());
    for (long unsigned int i = 0; i < mixed_out.size(); i++) {
        if (std::isnan(mixed_out[i])) {
            fstereo_audio[i] = 0; // Handle NaN values
        } else {
			fstereo_audio[i] = static_cast<short int>(mixed_out[i]*16384);
        }
    }
	fwrite(&fstereo_audio[0], sizeof(short int), fstereo_audio.size(), stdout);
}


void produce(Queue &my_queue,int block_size, real rf_Fs, real if_Fs, int Ntaps){
	
	std::vector<real> rf_coeff;
	std::vector<real> fm_demod;
	std::vector<real> fm_demod_block;
	real prev_I_state = 0.0f;
	real prev_Q_state = 0.0f;
	std::vector<real> i_sample(Ntaps - 1, 0.0);
	std::vector<real> q_sample(Ntaps - 1, 0.0);
	
	//RF_FRONTEND
	int gcd = findGCD(rf_Fs, if_Fs);
	int upsample_factor = if_Fs / gcd;
	int downsample_factor = rf_Fs / gcd;

	impulseResponseLPF_amp(rf_Fs, 100e3, Ntaps, rf_coeff, upsample_factor);
	for(unsigned int block_id=0; ; block_id++){
		std::vector<real> block_data(block_size);
		//std::cerr << "Before readblockdata: block_data size = " << block_data.size() << std::endl;
		readBlockData(block_size, block_id, block_data);
		//std::cerr << "After readblockdata: block_data size = " << block_data.size() << std::endl;
		if ((std::cin.rdstate()) != 0) {
			std::cerr << "End of input stream reached" << std::endl;
			exit(1);
		}
		else {
			std::cerr << "Read block " << block_id << std::endl;
			//std::cerr << "[DEBUG] block_data[0]: " << block_data[0] << std::endl;
			RF_frontEnd(block_data, rf_Fs, if_Fs, prev_I_state, prev_Q_state, 101, i_sample, q_sample, rf_coeff, fm_demod, fm_demod_block);
			my_queue.enqueue(fm_demod_block);
		}
	}
}

void consumer(Queue &my_queue, std::string processing_path, real if_Fs, real audio_Fs, int Ntaps){

	//MONO
    std::vector<real> mono_coeff;
	std::vector<short int> fmono_audio;
	std::vector<real> mono_sample(Ntaps - 1, 0.0);
    std::vector<real> mono_audio;
	int gcd_m = findGCD(if_Fs, audio_Fs);
    int upsample_factor_m = audio_Fs / gcd_m;
    int downsample_factor_m = if_Fs / gcd_m;
    int Ntaps_mono = Ntaps * upsample_factor_m;

	impulseResponseLPF_amp(if_Fs * upsample_factor_m, 16e3, Ntaps_mono, mono_coeff, upsample_factor_m);

	//STEREO
	real carrier_UB = 19500;
	real carrier_LB = 18500;
	real stereo_UB = 54000;
	real stereo_LB = 22000;
	real integrator = 0.0;
	real phaseEst = 0.0;
	real feedbackI = 1.0;
	real feedbackQ = 0.0;
	real ncoOut1 = 1.0;
	real trigOffset = 0.0;

	std::vector<real> carrier_coeff;
	std::vector<real> carrier_state(Ntaps-1,0.0);
	std::vector<real> stereo_coeff;
	std::vector<real> stereo_state(Ntaps-1,0.0);
	std::vector<real> mixed_coeff;
	std::vector<real> mixed_state(Ntaps-1,0.0);
	std::vector<short int> fstereo_audio;
	std::vector<real> mixed_out;
	std::vector<real> delay_state(51,0.0);
	std::vector<real> mono_delayed;
	std::vector<real> fm_demod_block;

	impulseResponseBPF(if_Fs, carrier_LB, carrier_UB, Ntaps, carrier_coeff);
	impulseResponseBPF(if_Fs, stereo_LB, stereo_UB, Ntaps, stereo_coeff);
	impulseResponseLPF_amp(if_Fs * upsample_factor_m, 16e3, Ntaps_mono, mixed_coeff, 1); // not too sure about amplifying

	while(true){
		my_queue.dequeue(fm_demod_block);
		if(processing_path == "m"){
			mono(fm_demod_block, if_Fs, audio_Fs, mono_sample, fmono_audio, upsample_factor_m, downsample_factor_m, Ntaps_mono, mono_coeff, mono_audio);
		}else if(processing_path == "s"){
			stereo (fm_demod_block,if_Fs, audio_Fs, downsample_factor_m, upsample_factor_m, delay_state,mono_sample,mono_coeff,mono_delayed, carrier_coeff, carrier_state,stereo_coeff, stereo_state, mixed_coeff,mixed_state, mixed_out,fstereo_audio,integrator, phaseEst, feedbackI, feedbackQ, ncoOut1,trigOffset);
		}
		if (std::cin.rdstate()!=0){
			std::cerr<<"End of stereo/mono processing"<<std::endl;
			exit(1);
		}
	}
}

int main(int argc, char *argv[])
{
	//check the commandline arguments
	if(argc <2){
		std::cerr << "Usage: " << argv[0] << " <mode> [processing_path] [block_size]" << std::endl;
		std::cerr << "Valid modes: 0, 1, 2, 3" << std::endl;
		std::cerr << "Valid processing paths: m, s, r" << std::endl;
		return 1;
	}

	//get the commandline arguments
	int mode = std::stoi(argv[1]);
	int block_size = 0;
	//default mono mode
	std::string processing_path = "m"; 
	setvbuf(stdout, NULL, _IONBF, 0);

	if(mode <0 || mode >3){
		std::cerr << "Error: Invalid mode. Valid modes are 0, 1, 2, and 3." << std::endl;
		return 1;
	}
	
	if(argc >= 3){
		processing_path = argv[2];
		if(processing_path != "m" && processing_path != "s" && processing_path != "r"){
			std::cerr << "Error: Invalid processing path. Valid processing paths are m, s, and r." << std::endl;
			return 1;
		}
	}
	
	if(argc >= 4){
	 	block_size = std::stoi(argv[3]);
	}
	
	std::cerr << "Running Mode: " << mode << ", Processing Path: " << processing_path << std::endl;

	ModeSettings settings = modeSettings[mode];
	real rf_Fs = settings.rf_Fs * 1e3; //convert to samples/second
	real if_Fs = settings.if_Fs * 1e3; //convert to samples/second
	real audio_Fs = settings.audio_Fs * 1e3; //convert to samples/second
	int sps = settings.sps * 2375; // samples/symbol
	int Ntaps = 101;

	block_size = rf_Fs * 0.04 * 2; //the 2 is sus

	std::cerr << "Block Size: " << block_size <<std::endl;

	std::cerr << "Mode: " << mode << "\nRF Fs: " << rf_Fs << "\nIF FS: " << if_Fs << "\nAudio Fs: " << audio_Fs << "\nSPS: " << sps << std::endl;

	std::cerr << "Working with reals on " << sizeof(real) << " bytes" << std::endl;

	//Queues
	Queue my_queue;

	std::thread produce_RF(produce, std::ref(my_queue),block_size, rf_Fs, if_Fs, Ntaps);

	std::thread consume_mono_stereo(consumer, std::ref(my_queue), processing_path, if_Fs, audio_Fs, Ntaps);

	produce_RF.join();
	consume_mono_stereo.join();

	//RDS path modes 0 and 2
	if(mode == 0 || mode == 2){
		int rds_rate = sps;
		std::cerr << "RDS symbol rate: " << rds_rate << " samples/second" << std::endl;
	}
	else {
		std::cerr << "RDS path is not supported in mode " << mode << std::endl;
	}

	return 0;
}
