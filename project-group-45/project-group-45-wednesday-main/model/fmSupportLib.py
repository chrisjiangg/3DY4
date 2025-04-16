#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import numpy as np
import math, cmath

#
# you should use the demodulator based on arctan given below as a reference
#
# in order to implement your OWN FM demodulator without the arctan function,
# a very good and to-the-point description is given by Richard Lyons at:
#
# https://www.embedded.com/dsp-tricks-frequency-demodulation-algorithms/
#
# the demodulator boils down to implementing equation (13-117) from above, where
# the derivatives are nothing else but differences between consecutive samples
#
# needless to say, you should not jump directly to equation (13-117)
# rather try first to understand the entire thought process based on calculus
# identities, like derivative of the arctan function or derivatives of ratios
#

# use the four quadrant arctan function for phase detect between a pair of
# IQ samples; then unwrap the phase and take its derivative to FM demodulate
def fmDemodArctan(I, Q, prev_phase = 0.0):
    #
    # the default prev_phase phase is assumed to be zero, however
    # take note in block processing it must be explicitly controlled

    # empty vector to store the demodulated samples
    fm_demod = np.empty(len(I))

    # iterate through each of the I and Q pairs
    for k in range(len(I)):

        # use the atan2 function (four quadrant version) to detect angle between
        # the imaginary part (quadrature Q) and the real part (in-phase I)
        current_phase = math.atan2(Q[k], I[k])

        # we need to unwrap the angle obtained in radians through arctan2
        # to deal with the case when the change between consecutive angles
        # is greater than Pi radians (unwrap brings it back between -Pi to Pi)
        [prev_phase, current_phase] = np.unwrap([prev_phase, current_phase])

        # take the derivative of the phase
        fm_demod[k] = current_phase - prev_phase

        # save the state of the current phase
        # to compute the next derivative
        prev_phase = current_phase

    # return both the demodulated samples as well as the last phase
    # (the last phase is needed to enable continuity for block processing)
    return fm_demod, prev_phase

#
# Custom FM demodulator based on phase difference
#
def fmDemodCustom(I, Q, prev_phase=0.0):
    """
    Custom FM demodulator based on phase difference.
    :param I: In-phase component.
    :param Q: Quadrature component.
    :param prev_phase: Previous phase (for block processing).
    :return: Demodulated signal and updated phase.
    """
    # Initialize the demodulated signal
    # Initialize demodulated signal array
    demodulated_signal = np.zeros(len(I))

    for index in range(len(I)):
        delta_Q = Q[index] - Q[index - 1] if index > 0 else 0
        delta_I = I[index] - I[index - 1] if index > 0 else 0

        numerator = (I[index] * delta_Q) - (Q[index] * delta_I)
        denominator = (I[index] ** 2 + Q[index] ** 2 + 1e-10)  # Prevent division by zero

        phase_difference = numerator / denominator
        current_phase = prev_phase + phase_difference

        # Unwrap phase to maintain continuity
        prev_phase, current_phase = np.unwrap([prev_phase, current_phase])

        # Store phase derivative as the demodulated signal
        demodulated_signal[index] = phase_difference

        # Update previous phase
        prev_phase = current_phase

    return demodulated_signal, prev_phase

# custom function for DFT that can be used by the PSD estimate
def DFT(x):

    # number of samples
    N = len(x)

    # frequency bins
    Xf = np.zeros(N, dtype='complex')

    # iterate through all frequency bins/samples
    for m in range(N):
        for k in range(N):
            Xf[m] += x[k] * cmath.exp(1j * 2 * math.pi * ((-k) * m) / N)

    # return the vector that holds the frequency bins
    return Xf

# custom function to estimate PSD based on the Bartlett method
# this is less accurate than the Welch method used in some packages
# however, as the visual inspections confirm, the estimate gives
# the user a "reasonably good" view of the power spectrum
def estimatePSD(samples, NFFT, Fs):

    # rename the NFFT argument (notation consistent with matplotlib.psd)
    # to freq_bins (i.e., frequency bins for which we compute the spectrum)
    freq_bins = NFFT
    # frequency increment (or resolution of the frequency bins)
    df = Fs / freq_bins

    # create the frequency vector to be used on the X axis
    # for plotting the PSD on the Y axis (only positive freq)
    freq = np.arange(0, Fs / 2, df)

    # design the Hann window used to smoothen the discrete data in order
    # to reduce the spectral leakage after the Fourier transform
    hann = np.empty(freq_bins)
    for i in range(len(hann)):
        hann[i] = 0.5 * (1 - math.cos(2 * math.pi * i / (freq_bins - 1)))

    # create an empty list where the PSD for each segment is computed
    psd_list = []

    # samples should be a multiple of frequency bins, so
    # the number of segments used for estimation is an integer
    # note: for this to work you must provide an argument for the
    # number of frequency bins not greater than the number of samples!
    no_segments = int(math.floor(len(samples) / float(freq_bins)))

    # iterate through all the segments
    for k in range(no_segments):

        # apply the hann window (using pointwise multiplication)
        # before computing the Fourier transform on a segment
        windowed_samples = samples[k * freq_bins:(k + 1) * freq_bins] * hann

        # compute the Fourier transform using the built-in FFT from numpy
        Xf = np.fft.fft(windowed_samples, freq_bins)

        # note, you can check how MUCH slower is DFT vs FFT by replacing the
        # above function call with the one that is commented below
        #
        # Xf = DFT(windowed_samples)
        #
        # note: the slow implementation of the Fourier transform is not as
        # critical when computing a static power spectra when troubleshooting
        #
        # note also: time permitting a custom FFT can be implemented

        # since input is real, we keep only the positive half of the spectrum
        # however, we will also add the signal energy of negative frequencies
        # to have a better and more accurate PSD estimate when plotting
        Xf = Xf[0:int(freq_bins / 2)] # keep only positive freq bins
        psd_seg = (1 / (Fs * freq_bins / 2)) * (abs(Xf)**2) # compute signal power
        psd_seg = 2 * psd_seg # add the energy from the negative freq bins

        # append to the list where PSD for each segment is stored
        # in sequential order (first segment, followed by the second one, ...)
        psd_list.extend(psd_seg)

    # iterate through all the frequency bins (positive freq only)
    # from all segments and average them (one bin at a time ...)
    psd_seg = np.zeros(int(freq_bins / 2))
    for k in range(int(freq_bins / 2)):
        # iterate through all the segments
        for l in range(no_segments):
            psd_seg[k] += psd_list[k + l * int(freq_bins / 2)]
        # compute the estimate for each bin
        psd_seg[k] = psd_seg[k] / no_segments

    # translate to the decibel (dB) scale
    psd_est = np.zeros(int(freq_bins / 2))
    for k in range(int(freq_bins / 2)):
        psd_est[k] = 10 * math.log10(psd_seg[k])

    # the frequency vector and PSD estimate
    return freq, psd_est

# custom function to format the plotting of the PSD
def fmPlotPSD(ax, samples, Fs, height, title):

    x_major_interval = (Fs / 12)        # adjust grid lines as needed
    x_minor_interval = (Fs / 12) / 4
    y_major_interval = 20
    x_epsilon = 1e-3
    x_max = x_epsilon + Fs / 2        # adjust x/y range as needed
    x_min = 0
    y_max = 10
    y_min = y_max - 100 * height
    ax.psd(samples, NFFT=512, Fs=Fs)
    #
    # below is the custom PSD estimate, which is based on the Bartlett method
    # it is less accurate than the PSD from matplotlib, however it is sufficient
    # to help us visualize the power spectra on the acquired/filtered data
    #
    # freq, my_psd = estimatePSD(samples, NFFT=512, Fs=Fs)
    # ax.plot(freq, my_psd)
    #
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])
    ax.set_xticks(np.arange(x_min, x_max, x_major_interval))
    ax.set_xticks(np.arange(x_min, x_max, x_minor_interval), minor=True)
    ax.set_yticks(np.arange(y_min, y_max, y_major_interval))
    ax.grid(which='major', alpha=0.75)
    ax.grid(which='minor', alpha=0.25)
    ax.set_xlabel('Frequency (kHz)')
    ax.set_ylabel('PSD (db/Hz)')
    ax.set_title(title)

def fmPll(pllIn, freq, Fs, state, ncoScale=1.0, phaseAdjust=0.0, normBandwidth=0.01):
    Cp = 2.666
    Ci = 3.555

    Kp = normBandwidth * Cp
    Ki = (normBandwidth ** 2) * Ci

    ncoOut = np.empty(len(pllIn))

    # Load saved state
    integrator = state['integrator']
    phaseEst = state['phaseEst']
    feedbackI = state['feedbackI']
    feedbackQ = state['feedbackQ']
    trigOffset = state['trigOffset']

    for k in range(len(pllIn)):
        errorI = pllIn[k] * (+feedbackI)
        errorQ = pllIn[k] * (-feedbackQ)
        errorD = math.atan2(errorQ, errorI)

        integrator += Ki * errorD
        phaseEst += Kp * errorD + integrator

        trigOffset += 1
        trigArg = 2 * math.pi * (freq / Fs) * trigOffset + phaseEst
        feedbackI = math.cos(trigArg)
        feedbackQ = math.sin(trigArg)

        ncoOut[k] = math.cos(trigArg * ncoScale + phaseAdjust)

    # Save updated state
    state['integrator'] = integrator
    state['phaseEst'] = phaseEst
    state['feedbackI'] = feedbackI
    state['feedbackQ'] = feedbackQ
    state['trigOffset'] = trigOffset

    return ncoOut, state

    
def delayBlock(input_block, state_block):
    output_block = np.concatenate((state_block, input_block[:-len(state_block)]))
    state_block = input_block[-len(state_block):]
    return output_block, state_block

def design_custom_bandpass_filter(f_low, f_high, fs, N_taps):
    center = (N_taps - 1) / 2
    f_mid = (f_low + f_high) / 2
    h = np.zeros(N_taps)
    scale_factor = 0

    for k in range(N_taps):
        n = k - center
        if n == 0:
            h[k] = (2 * f_high / fs) - (2 * f_low / fs)
        else:
            h[k] = (math.sin(2 * math.pi * f_high * n / fs) / (math.pi * n)) - \
                   (math.sin(2 * math.pi * f_low * n / fs) / (math.pi * n))

        # Apply Hann window
        hann = 0.5 - 0.5 * math.cos(2 * math.pi * k / (N_taps - 1))
        h[k] *= hann

        # Update scale factor
        scale_factor += h[k] * math.cos(2 * math.pi * n * f_mid / fs)

    # Normalize
    h /= scale_factor

    return h

#RDS FUNCTIONS
def fmPll_RDS(pllIn, freq, Fs, state, ncoScale=1.0, phaseAdjust=0.0, normBandwidth=0.01):
    Cp = 2.666
    Ci = 3.555

    Kp = normBandwidth * Cp
    Ki = (normBandwidth ** 2) * Ci

    ncoOutI = np.empty(len(pllIn))
    ncoOutQ = np.empty(len(pllIn))

    # Load saved state
    integrator = state['integrator']
    phaseEst = state['phaseEst']
    feedbackI = state['feedbackI']
    feedbackQ = state['feedbackQ']
    trigOffset = state['trigOffset']

    for k in range(len(pllIn)):
        # Phase detector (input * conjugate of NCO)
        errorI = pllIn[k] * (+feedbackI)
        errorQ = pllIn[k] * (-feedbackQ)
        errorD = math.atan2(errorQ, errorI)

        # Loop filter
        integrator += Ki * errorD
        phaseEst += Kp * errorD + integrator

        # NCO update
        trigOffset += 1
        trigArg = 2 * math.pi * (freq / Fs) * trigOffset + phaseEst
        feedbackI = math.cos(trigArg)
        feedbackQ = math.sin(trigArg)

        # Store current output
        ncoOutI[k] = math.cos(trigArg * ncoScale + phaseAdjust)
        ncoOutQ[k] = math.sin(trigArg * ncoScale + phaseAdjust)

    # Save updated state
    state['integrator'] = integrator
    state['phaseEst'] = phaseEst
    state['feedbackI'] = feedbackI
    state['feedbackQ'] = feedbackQ
    state['trigOffset'] = trigOffset

    return ncoOutI, ncoOutQ, state

def square_and_pll(signal, Fs, state):
    """
    Squares the RDS signal and uses PLL to lock onto the 114kHz tone (2x57kHz).
    """
    squared = signal ** 2
    return fmPll_RDS(squared, freq=114000, Fs=Fs, state=state)

def bpsk_demod(signal, carrier):
    return signal * carrier

def threshold_and_sample(signal):
    return np.where(signal > 0, 1, 0)

def matched_filter_rds(symbol_rate, Fs, taps=101):
    symbol_period = int(Fs / symbol_rate)
    return np.ones(symbol_period) / symbol_period

def root_raised_cosine(num_taps, sps, alpha):
    # sps: samples per symbol
    t = np.arange(-num_taps//2, num_taps//2 + 1)
    t = np.where(t == 0, 1e-8, t)  # prevent div by 0
    numerator = np.sin(np.pi * t * (1 - alpha) / sps) + \
                4 * alpha * t / sps * np.cos(np.pi * t * (1 + alpha) / sps)
    denominator = np.pi * t * (1 - (4 * alpha * t / sps) ** 2) / sps
    rrc = numerator / denominator
    rrc /= np.sqrt(np.sum(rrc**2))  # normalize power
    return rrc
def impulseResponseRootRaisedCosine(Fs, N_taps):

    """
    Root raised cosine (RRC) filter

    Fs          sampling rate at the output of the resampler in the RDS path
                sampling rate must be an integer multipler of 2375
                this integer multiple is the number of samples per symbol

    N_taps      number of filter taps

    """

    # duration for each symbol - should NOT be changed for RDS!
    T_symbol = 1/2375.0

    # roll-off factor (must be greater than 0 and smaller than 1)
    # for RDS a value in the range of 0.9 is a good trade-off between
    # the excess bandwidth and the size/duration of ripples in the time-domain
    beta = 0.90

    # the RRC inpulse response that will be computed in this function
    impulseResponseRRC = np.empty(N_taps)

    for k in range(N_taps):
        t = float((k-N_taps/2))/Fs
        # we ignore the 1/T_symbol scale factor
        if t == 0.0: impulseResponseRRC[k] = 1.0 + beta*((4/math.pi)-1)
        elif t == -T_symbol/(4*beta) or t == T_symbol/(4*beta):
            impulseResponseRRC[k] = (beta/np.sqrt(2))*(((1+2/math.pi)* \
                    (math.sin(math.pi/(4*beta)))) + ((1-2/math.pi)*(math.cos(math.pi/(4*beta)))))
        else: impulseResponseRRC[k] = (math.sin(math.pi*t*(1-beta)/T_symbol) +  \
                    4*beta*(t/T_symbol)*math.cos(math.pi*t*(1+beta)/T_symbol))/ \
                    (math.pi*t*(1-(4*beta*t/T_symbol)*(4*beta*t/T_symbol))/T_symbol)

    # returns the RRC impulse response to be used by convolution
    return impulseResponseRRC
if __name__ == "__main__":

    # do nothing when this module is launched on its own
    pass
