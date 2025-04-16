/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>

#define PI 3.14159265358979323846

// For this program we use ONLY 32-bit floats due to how the
// wavio.py Python script processes .wav files (see model sub-folder)
typedef float real;

// Function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(real Fs, real Fc, unsigned short int num_taps, std::vector<real> &h) {
	// Allocate memory for the impulse response
	h.clear();
	h.resize(num_taps, 0.0);

	// The rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
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
void convolveFIR(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h) {
	// Allocate memory for the output (filtered) data
	y.clear();
	y.resize(x.size() + h.size() - 1, 0.0);

	// The rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	for(int i=0; i < y.size(); i++){
		for(int j=0; j < h.size(); j++){
			if(((i-j)>=0) && ((i-j)< x.size())){
				y[i] += h[j] * x[i-j]; 
			}
		}
	}
}

//these functions are used in the unit testing
//block-based convolution
void blockConvolveFIR(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h, int blockSize, std::vector<real> &state) { 
    y.clear();
    y.resize(x.size(), 0.0);

    for (int i = 0; i < x.size(); i += blockSize) {
        int blockEnd = std::min(i + blockSize, (int)x.size());
        for (int j = i; j < blockEnd; j++) {
            y[j] = 0.0;
            for (int k = 0; k < h.size(); k++) {
                if (j - k >= 0) {
                    y[j] += h[k] * x[j - k];
                } else {
                    y[j] += h[k] * state[state.size() + (j - k)];
                }
            }
        }
        for (int j = 0; j < h.size() - 1; j++) {
            state[j] = x[blockEnd - 1 - j];
        }
    }
}

//provides random samples
void generateRandomSamples(std::vector<real> &x, int length, real maxValue) {
    x.resize(length);
    for (int i = 0; i < length; i++) {
        x[i] = static_cast<real>(rand()) / RAND_MAX * maxValue;
    }
}

bool areVectorsClose(const std::vector<real> &a, const std::vector<real> &b, real precision) {
    if (a.size() != b.size()) return false;
    for (int i = 0; i < a.size(); i++) {
        if (std::abs(a[i] - b[i]) > precision) return false;
    }
    return true;
}

real maxDifference(const std::vector<real> &a, const std::vector<real> &b) {
    real maxDiff = 0.0;
    for (int i = 0; i < a.size(); i++) {
        maxDiff = std::max(maxDiff, std::abs(a[i] - b[i]));
    }
    return maxDiff;
}

void unitTestBlockConvolution(int xLength, int hLength, int numBlocks, real maxValue, real precision) {
    std::vector<real> x, h;
    generateRandomSamples(x, xLength, maxValue);
    generateRandomSamples(h, hLength, maxValue);

    std::vector<real> yFull;
    convolveFIR(yFull, x, h);

    // Perform block convolution with the original number of blocks
    std::vector<real> yBlock, state(h.size() - 1, 0.0);
    blockConvolveFIR(yBlock, x, h, xLength / numBlocks, state);

    // Perform block convolution with double the number of blocks
    std::vector<real> yBlockDouble, stateDouble(h.size() - 1, 0.0);
    blockConvolveFIR(yBlockDouble, x, h, xLength / (numBlocks * 2), stateDouble);

    // Perform block convolution with half the number of blocks
    std::vector<real> yBlockHalf, stateHalf(h.size() - 1, 0.0);
    blockConvolveFIR(yBlockHalf, x, h, xLength / (numBlocks / 2), stateHalf);

    // Compare results
    if (!areVectorsClose(yFull, yBlock, precision) || !areVectorsClose(yFull, yBlockDouble, precision) || !areVectorsClose(yFull, yBlockHalf, precision)) {
        std::cout << "Unit test failed: Block convolution result does not match full convolution.\n";
        std::cout << "Max difference (original blocks): " << maxDifference(yFull, yBlock) << "\n";
        std::cout << "Max difference (double blocks): " << maxDifference(yFull, yBlockDouble) << "\n";
        std::cout << "Max difference (half blocks): " << maxDifference(yFull, yBlockHalf) << "\n";
    } else {
        std::cout << "Unit test passed.\n";
    }
}


// Function to read audio data from a binary file that contains raw samples
// represented as 32-bit reals; we also assume two audio channels
// Note: Check the Python script that can prepare this type of files
void read_audio_data(const std::string in_fname, std::vector<real> &audio_data) {
	// File descriptor for the input to be read
	std::ifstream fdin(in_fname, std::ios::binary);
	if (!fdin) {
		std::cout << "File " << in_fname << " not found ... exiting\n";
		exit(1);
	} else {
		std::cout << "Reading raw audio from \"" << in_fname << "\"\n";
	}
	// Search for end of file to count the number of samples to be read
	fdin.seekg(0, std::ios::end);
	// We assume the Python script has written data in 32-bit reals
	const unsigned int num_samples = fdin.tellg() / sizeof(real);

	// Allocate memory space to store all the samples
	audio_data.clear();
	audio_data.resize(num_samples);
	// Back to the beginning of the file to read all samples at once
	fdin.seekg(0, std::ios::beg);
	// Do a single read for audio data from the input file stream
	fdin.read(reinterpret_cast<char *>(&audio_data[0]), num_samples * sizeof(real));
	// Close the input file
	fdin.close();
}

// Function to split an audio data where the left channel is in even samples
// and the right channel is in odd samples
void split_audio_into_channels(const std::vector<real> &audio_data, std::vector<real> &audio_left, std::vector<real> &audio_right) {
	for (int i = 0; i < (int)audio_data.size(); i++) {
		if (i % 2 == 0) {
			audio_left.push_back(audio_data[i]);
		} else {
			audio_right.push_back(audio_data[i]);
		}
	}
}

// Function to write audio data to a binary file that contains raw samples
// represented as 32-bit reals; we also assume two audio channels
// Note: Check the Python script that can read this type of files
// and then reformat them to .wav files to be run on third-party players
void write_audio_data(const std::string out_fname, const std::vector<real> &audio_left, const std::vector<real> &audio_right) {
	// File descriptor for the output to be written
	if (audio_left.size() != audio_right.size()) {
		std::cout << "Something got messed up with audio channels\n";
		std::cout << "They must have the same size ... exiting\n";
		exit(1);
	} else {
		std::cout << "Writing raw audio to \"" << out_fname << "\"\n";
	}
	std::ofstream fdout(out_fname, std::ios::binary);
	for (int i = 0; i < (int)audio_left.size(); i++) {
		// We assume we have handled a stereo audio file
		// Hence, we must interleave the two channels
		// (Change as needed if testing with mono files)
		fdout.write(reinterpret_cast<const char *>(&audio_left[i]), sizeof(audio_left[i]));
		fdout.write(reinterpret_cast<const char *>(&audio_right[i]), sizeof(audio_right[i]));
	}
	fdout.close();
}

int main() {
    // By default, we use only floats in this program
    std::cout << "Working with reals on " << sizeof(real) << " bytes" << std::endl;

    // Assume the wavio.py script was run beforehand to produce a binary file
    const std::string in_fname = "../data/float32samples.bin";
    // Declare vector where the audio data will be stored
    std::vector<real> audio_data;
    // Note: We allocate memory for audio_data from within this read function
    read_audio_data(in_fname, audio_data);

    // Set up the filtering flow
    real Fs = 44100.0;  // sample rate for our "assumed" audio (change as needed for 48 ksamples/sec audio files)
    real Fc = 10000.0;  // cutoff frequency (explore ... but up-to Nyquist only!)
    // Number of FIR filter taps (feel free to explore ...)
    unsigned short int num_taps = 51;

    // Impulse response (reuse code from the previous experiment)
    std::vector<real> h;
    impulseResponseLPF(Fs, Fc, num_taps, h);
    // Note: Memory for the impulse response vector and output data vectors
    // should be allocated from within the corresponding functions
    // (as for the previous experiment, from where you should reuse your code)

    // There is one more point before filtering is done:
    // Recall we assume there are two channels in the audio data, hence
    // the channels must be handled separately by your DSP functions; therefore
    // split the audio_data into two channels (audio_left and audio_right)

    // Declare vectors where the audio left/right channels will be stored
    std::vector<real> audio_left, audio_right;
    // Note: We allocate the memory for the left/right channels
    // from within the split function that is called in the code below
    split_audio_into_channels(audio_data, audio_left, audio_right);

    // Convolution code for filtering (reuse from the previous experiment)
    std::vector<real> single_pass_left, single_pass_right;
    convolveFIR(single_pass_left, audio_left, h);
    convolveFIR(single_pass_right, audio_right, h);
    // Note: By default the above convolution produces zero on the output stream
    // YOU will need to update the convolveFIR and impulseResponseLPF functions

    // Create a binary file to be read by wavio.py script to produce a .wav file
    // Note: Small adjustments will need to be made to wavio.py, i.e., you should
    // match the filenames, no need for self-checks as default Python code, ...
    const std::string out_fname = "../data/float32filtered.bin";
    write_audio_data(out_fname, single_pass_left, single_pass_right);

    // Initialize state vectors for block processing
    std::vector<real> state_left(h.size() - 1, 0.0);
    std::vector<real> state_right(h.size() - 1, 0.0);

    // Perform block convolution on left and right channels
    std::vector<real> block_processed_left, block_processed_right;
    int blockSize = 1024; // Example block size, adjust as needed
    blockConvolveFIR(block_processed_left, audio_left, h, blockSize, state_left);
    blockConvolveFIR(block_processed_right, audio_right, h, blockSize, state_right);

    // Write the block-processed audio data to a new binary file
    const std::string block_out_fname = "../data/float32blockfiltered.bin";
    write_audio_data(block_out_fname, block_processed_left, block_processed_right);

    // Perform unit test for block convolution
    int xLength = 1000;
    int hLength = 51;
    int numBlocks = 10;
    real maxValue = 1.0;
    real precision = 1e-5;
    unitTestBlockConvolution(xLength, hLength, numBlocks, maxValue, precision);

    return 0;
}
