#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math
from fmSupportLib import *

# Custom FIR filter and filtering functions
def design_fir_filter(cutoff_freq, sample_rate, num_taps):
    normalized_cutoff = cutoff_freq / (sample_rate / 2)
    coeffs = np.zeros(num_taps)
    center = (num_taps - 1) // 2

    for idx in range(num_taps):
        if idx == center:
            coeffs[idx] = normalized_cutoff
        else:
            coeffs[idx] = (
                normalized_cutoff
                * np.sin(np.pi * normalized_cutoff * (idx - center))
                / (np.pi * normalized_cutoff * (idx - center))
            )
        coeffs[idx] *= (0.5 - 0.5 * np.cos(2 * np.pi * idx / (num_taps - 1)))

    return coeffs

def apply_fir_filter(coeffs, signal, prev_state=None, chunk_size=None):
    output = np.zeros(len(signal))

    if prev_state is None:
        prev_state = np.zeros(len(coeffs) - 1)

    if chunk_size is None:
        for i in range(len(signal)):
            output[i] = 0
            for j in range(len(coeffs)):
                if i - j >= 0:
                    output[i] += signal[i - j] * coeffs[j]
                else:
                    output[i] += prev_state[i - j] * coeffs[j]
        prev_state = signal[-(len(coeffs) - 1):]
    else:
        for start in range(0, len(signal), chunk_size):
            stop = min(start + chunk_size, len(signal))
            segment = signal[start:stop]
            for i in range(len(segment)):
                output[start + i] = 0
                for j in range(len(coeffs)):
                    if i - j >= 0:
                        output[start + i] += segment[i - j] * coeffs[j]
                    else:
                        output[start + i] += prev_state[i - j] * coeffs[j]
            prev_state = segment[-(len(coeffs) - 1):]

    return output, prev_state

# RF parameters
rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 101
rf_decim = 10

# Stereo parameters
pilot_lower = 18500
pilot_upper = 19500
stereo_lower= 22000
stereo_upper= 54000

# RDS parameters
rds_lower = 54000
rds_upper = 60000

# Audio parameters
audio_Fs = 48e3
audio_Fc = 16e3
audio_taps = 101
audio_decim = 5
audio_input = 240e3

# Mono and Stereo
mono = False
stereo = True  
rds = False

if __name__ == "__main__":
    in_fname = "../data/samples0.raw"
    raw_data = np.fromfile(in_fname, dtype='uint8')
    print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")

    iq_data = (np.float32(raw_data) - 128.0) / 128.0
    print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

    rf_coeff = signal.firwin(rf_taps, rf_Fc / (rf_Fs / 2), window=('hann'))
    audio_coeff = design_fir_filter(audio_Fc, audio_input, audio_taps)

    subfig_height = np.array([0.8, 2, 1.6])
    plt.rc('figure', figsize=(7.5, 7.5))
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
    fig.subplots_adjust(hspace=.6)

    i_filt = signal.lfilter(rf_coeff, 1.0, iq_data[0::2])
    q_filt = signal.lfilter(rf_coeff, 1.0, iq_data[1::2])
    i_ds = i_filt[::rf_decim]
    q_ds = q_filt[::rf_decim]

    fm_demod, _ = fmDemodArctan(i_ds, q_ds)
    fmPlotPSD(ax0, fm_demod, (rf_Fs / rf_decim) / 1e3, subfig_height[0], 'Demodulated FM (Full Recording)')

    if mono:
        audio_filt = signal.lfilter(audio_coeff, 1.0, fm_demod)
        audio_data = audio_filt[::audio_decim]
        
        # Create a new figure for mono plots
        fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
        fig.subplots_adjust(hspace=.6)
        
        # Plot 1: Demodulated FM
        fmPlotPSD(ax0, fm_demod, (rf_Fs / rf_decim) / 1e3, subfig_height[0], 'Demodulated FM (Full Recording)')
        
        # Plot 2: Extracted Mono
        fmPlotPSD(ax1, audio_filt, (rf_Fs / rf_decim) / 1e3, subfig_height[1], 'Extracted Mono (Full Recording)')
        
        # Plot 3: Downsampled Mono Audio
        fmPlotPSD(ax2, audio_data, audio_Fs / 1e3, subfig_height[2], 'Downsampled Mono Audio')
        
        # Save the mono plots to a PNG file
        plt.savefig('../data/fmMonoPlots.png', dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        wavfile.write("../data/fmMonoBlock.wav", int(audio_Fs), np.int16((audio_data / 2) * 32767))
        print("Written mono audio samples to \"../data/fmMonoBlock.wav\" in signed 16-bit format")
        print("Saved mono PSD plots to \"../data/fmMonoPlots.png\"")

    if stereo:
        pilot_coeff = signal.firwin(audio_taps, [pilot_lower/(rf_Fs/rf_decim/2), pilot_upper/(rf_Fs/rf_decim/2)], pass_zero=False, window='hann')
        pilot_signal = signal.lfilter(pilot_coeff, 1.0, fm_demod)

        stereo_coeff = signal.firwin(audio_taps, [stereo_lower/(rf_Fs/rf_decim/2), stereo_upper/(rf_Fs/rf_decim/2)], pass_zero=False, window='hann')
        stereo_sub = signal.lfilter(stereo_coeff, 1.0, fm_demod)

        mixed = stereo_sub * pilot_signal * 2

        diff_lp = signal.firwin(audio_taps, audio_Fc/(rf_Fs/rf_decim/2), window='hann')
        l_minus_r = signal.lfilter(diff_lp, 1.0, mixed)

        l_plus_r = signal.lfilter(audio_coeff, 1.0, fm_demod)

        left = 0.5 * (l_plus_r + l_minus_r)
        right = 0.5 * (l_plus_r - l_minus_r)

        left_ds = left[::audio_decim]
        right_ds = right[::audio_decim]

        # Create a new figure for stereo plots
        fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
        fig.subplots_adjust(hspace=.6)
        
        # Plot 1: Demodulated FM
        fmPlotPSD(ax0, fm_demod, (rf_Fs / rf_decim) / 1e3, subfig_height[0], 'Demodulated FM (Full Recording)')
        # Plot 2: Stereo average
        fmPlotPSD(ax1, l_plus_r, (rf_Fs / rf_decim) / 1e3, subfig_height[1], 'Stereo Average (L+R) Channel')
        # Plot 3: Left and Right channels
        fig2, (ax3, ax4) = plt.subplots(nrows=2, gridspec_kw={'height_ratios': [1, 1]})
        fig2.subplots_adjust(hspace=.6)
        fmPlotPSD(ax3, left_ds, audio_Fs / 1e3, 1, 'Stereo Left Channel')
        fmPlotPSD(ax4, right_ds, audio_Fs / 1e3, 1, 'Stereo Right Channel')

        # Save the stereo plots to PNG files
        plt.figure(fig.number)  # Switch back to the first figure
        plt.savefig('../data/fmStereoPlots_1.png', dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        plt.figure(fig2.number)  # Switch to the second figure (left/right channels)
        plt.savefig('../data/fmStereoPlots_2.png', dpi=300, bbox_inches='tight')
        plt.close(fig2)

        stereo_out = np.stack((left_ds, right_ds), axis=1)
        wavfile.write("../data/fmStereoBlock.wav", int(audio_Fs), np.int16(stereo_out * 32767))
        print("Wrote stereo .wav file => ../data/fmStereoBlock.wav")
        print("Saved stereo PSD plots to \"../data/fmStereoPlots_1.png\" and \"../data/fmStereoPlots_2.png\"")
    #if rds:
	    
