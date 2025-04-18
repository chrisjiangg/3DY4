#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import numpy as np
from scipy.io import wavfile

#
# Since manipulating .wav files is not the objective of the SDR project and
# we are using them solely for "assessing" the outcome of the DSP tasks
# while troubleshooting, we will avoid processing any .wav files in C++,
# mainly because of the error-prone nature of handling .wav file headers.
#
# For the reason above, the Python script below can be used to parse/format
# .wav files to/from binary files where the sample representation is known
# (or better said agreed on) by both the Python script and the C++ program.
#
# .wav files should be opened only in this Python script, and samples written
# in binary (e.g., assuming 32-bit floating point for this example) should be
# read by the C++ program in binary format (raw data, no headers). Subsequently,
# the C++ program should output the processed data also in binary format,
# which can be read back by this Python script to be formatted properly with a
# header into a .wav file that can then be used on a third-party audio player.
#

if __name__ == "__main__":

	# Parse an audio file
	audio_Fs, audio_data = wavfile.read("../data/audio_test.wav")
	print(' Audio sample rate = {0:f} \
		\n Number of channels = {1:d} \
		\n Number of samples = {2:d}' \
		.format(audio_Fs, audio_data.ndim, len(audio_data)))

	# Output binary file name (where samples are written from Python)
	out_fname = "../data/float32samples.bin"
	# Dump audio data in a binary file where each sample is a 32-bit float
	audio_data.astype('float32').tofile(out_fname)
	print(" Written binary data to \"" + out_fname + "\" in float32 format")

	# Input binary file name (from where samples are read into Python)
	# The default is just a self-check; of course, change filenames as needed
	in_fname = "../data/float32samples.bin"
	# in_fname = "../data/float32filtered.bin"
	# Read data from a binary file (assuming 32-bit floats)
	float_data = np.fromfile(in_fname, dtype='float32')
	print(" Read binary data from \"" + in_fname + "\" in float32 format")

	# We assume below there are two audio channels where data is
	# interleaved, i.e., left channel sample, right channel sample, ...
	# for mono .wav files the reshaping below is unnecessary
	reshaped_data = np.reshape(float_data, (-1, 2))

	# Self-check if the read and write are working correctly
	# not needed while working with data generated from C++
	print(" Are the two data sets identical ? " +
			str(np.array_equal(audio_data, reshaped_data)))

	wavfile.write("../data/audio_processed.wav", \
				audio_Fs, \
				reshaped_data.astype(np.int16))

	# Note: We can also dump audio data in other formats, if needed
	# audio_data.astype('int16').tofile('int16samples.bin')
