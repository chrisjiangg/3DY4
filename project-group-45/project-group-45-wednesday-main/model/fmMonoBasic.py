#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

"""
Ten prerecorded I/Q sample files can be downloaded from the link sent via email
and posted on Avenue.

There is NO need to copy all the sample files at once into your data sub-folder.
Instead, it is recommended to download one file at a time and overwrite the
iq_samples.raw file in the data sub-folder, which is the default path used by the script.

For those of you with the RF dongle and Raspberry Pi kit at home, the command-line
instructions for recording RF data are provided below. After installing the necessary
drivers to work with the RF dongle, the 8-bit unsigned values for the I/Q pairs can be
recorded using the following command:

rtl_sdr -f 99.9M -s 2.4M - > iq_samples.raw

The above command assumes that we are tuned to the FM station at 99.9 MHz, using an
RF sample rate of 2.4 Msamples/sec, and saving the data to a file named iq_samples.raw
(you can change the filename as needed).

For the above command, data acquisition runs indefinitely and must be stopped manually
by pressing Ctrl+C. If you want to stop recording after a predefined number of samples,
e.g., 12 million I/Q pairs (equivalent to 5 seconds at 2.4 Msamples/sec), you can
include an additional argument:

rtl_sdr -f 99.9M -s 2.4M -n 12000000 - > iq_samples.raw

To verify that the raw I/Q data has been recorded properly, place the file in your
project repository's data sub-folder and run the Python script from the model sub-folder.
The script should generate both .png image files (showing PSD estimates) and a .wav file.

In the source code below (at the beginning of the main function), you can observe where
the raw_data is read and where the normalization of the 8-bit unsigned I/Q samples to
32-bit float samples (in the range -1 to +1) is performed. While using 32-bit floats and
normalizing to the range -1 to +1 are optional choices (commonly adopted by many
third-party SDR software implementations), it is up to each project group to decide how
to handle the 8-bit unsigned I/Q samples in their Python model and C++ implementation.
For example, one can choose how to reformat the data in the range -1 to +1 in 64-bit double
format (as shown in the lines commented below the normalization to 32-bit floats). Note,
however, it is recommended to use the same data format in both the Python model and C++.

A final but critical note: the .gitignore file should NOT be modified to allow pushing
.raw files to GitHub. Including .raw files in the repository will result in very large
repositories that take excessive time to clone, pull, or push. As outlined in the
reference .gitignore file, only source files should be kept in your repositories.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import wavfile
from scipy import signal
import math
from fmSupportLib import fmDemodArctan, fmPlotPSD

# Custom FIR filter and filtering functions
def custom_firwin(num_taps, cutoff):
    n = np.arange(num_taps)
    h = np.sinc(2 * cutoff * (n - (num_taps - 1) / 2))
    h *= np.hanning(num_taps)
    h /= np.sum(h)
    return h

def custom_lfilter(coeff, x, state=None):
    if state is None:
        state = np.zeros(len(coeff) - 1)
    y = np.convolve(coeff, x)
    if len(state) > 0:
        y[:len(state)] += state
        state = y[-len(state):]
    return y, state

# RF parameters
rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 101
rf_decim = 10

# Audio parameters
audio_Fs = 48e3
audio_Fc = 16e3  # Cutoff frequency for mono audio (16 kHz)
audio_taps = 101  # Number of filter taps
audio_decim = 5  # Decimation factor to reach 48 Ksamples/sec

if __name__ == "__main__":

    # Read the raw IQ data
    in_fname = "../data/samples0.raw"
    raw_data = np.fromfile(in_fname, dtype='uint8')
    print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")

    # Normalize raw IQ data to 32-bit float format (-1 to +1)
    iq_data = (np.float32(raw_data) - 128.0) / 128.0
    print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

    # Coefficients for front-end low-pass filter
    rf_coeff = custom_firwin(rf_taps, rf_Fc / (rf_Fs / 2))

    # Filter to extract the FM channel
    i_filt= np.convolve(rf_coeff, iq_data[0::2]) #I/Q is interleaved together this seperates them as well.
    q_filt = np.convolve(rf_coeff, iq_data[1::2]) 

    # Downsample the FM channel
    i_ds = i_filt[::rf_decim] #DownSample to 240KSamples/sec
    q_ds = q_filt[::rf_decim]

    

    # Set up subfigures for plotting
    subfig_height = np.array([0.8, 2, 1.6])
    plt.rc('figure', figsize=(7.5, 7.5))
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
    fig.subplots_adjust(hspace=.6)

    # FM demodulator
    fm_demod, _ = fmDemodArctan(i_ds, q_ds)

    # PSD after FM demodulation
    fmPlotPSD(ax0, fm_demod, (rf_Fs / rf_decim) / 1e3, subfig_height[0], 'Demodulated FM (full recording)')

    # Coefficients for the audio low-pass filter
    audio_coeff = custom_firwin(audio_taps, audio_Fc / (rf_Fs / rf_decim / 2))

    # Extract mono audio data
    audio_filt = np.convolve(audio_coeff, fm_demod)

    # PSD for mono audio
    fmPlotPSD(ax1, audio_filt, (rf_Fs / rf_decim) / 1e3, subfig_height[1], 'Extracted Mono')

    # Downsample mono audio
    audio_data = audio_filt[::audio_decim]

    # PSD for downsampled audio
    fmPlotPSD(ax2, audio_data, audio_Fs / 1e3, subfig_height[2], 'Downsampled Mono Audio')

    # Save plots
    fig.savefig("../data/fmMonoBasic.png")
    plt.show()

    # Write audio data to file
    out_fname = "../data/fmMonoBasic.wav"
    wavfile.write(out_fname, int(audio_Fs), np.int16((audio_data / 2) * 32767))
    print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")