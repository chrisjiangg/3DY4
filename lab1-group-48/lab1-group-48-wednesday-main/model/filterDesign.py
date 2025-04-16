#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
import sys

# use generateSin/plotTime from the fourierTransform module
from fourierTransform import generateSin, plotTime, plotSpectrum

def freqzPlot(coeff, Fs, msg):

	# find the frequency response using freqz from SciPy:
	# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.freqz.html
	w, h = signal.freqz(coeff)

	# Reminder: np.pi rad/sample is actually the Nyquist frequency
	w = w * Fs / (2 * np.pi)  # needed to draw the frequency on the X axis

	# plots the magnitude response where the x axis is normalized in rad/sample
	fig, ax1 = plt.subplots()
	ax1.set_title('Digital filter frequency response (' + msg + ')')
	ax1.plot(w, 20 * np.log10(abs(h)), 'b')
	ax1.set_ylabel('Amplitude [dB]', color = 'b')
	ax1.set_xlabel('Frequency [Hz]')

	# uncomment the lines below if you wish to inspect the phase response
	# Note: as important as the phase response, it is not critical within our context
	#
	# we expect a linear phase in the passband, i.e., no phase distortion because
	# all frequency components of the input delayed by the same amount

	# ax2 = ax1.twinx()
	# angles = np.unwrap(np.angle(h))
	# ax2.plot(w, angles, 'g')
	# ax2.set_ylabel('Angle (radians)', color = 'g')

def filterSin(Fs, Fc, coeff):

	# we can control the frequency relative to the filter cutoff
	time, x = generateSin(Fs, interval = 1.0, frequency = Fc * 0.4)
	plotTime(x, time)

	# use lfilter from SciPy for FIR filtering:
	# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.lfilter.html
	fx = signal.lfilter(coeff, 1.0, x)

	# you should clearly observe the effects (attenuation, delay) introduced by the filter
	plotTime(fx, time)

def taps(fc,fs,N_taps1):
	Norm = fc / (fs / 2)
	h = np.zeros(N_taps) 
	middle = (N_taps - 1) / 2
	for i in range(N_taps):
		if (i == middle):
			h[i] = Norm
		else:
			h[i] = Norm * (np.sin(np.pi * Norm * (i - middle)) / (np.pi * Norm * (i - middle)))
		h[i] = h[i] * ((1/2) - (1/2) * np.cos((2 * np.pi * i) / (N_taps - 1)))
	return h

def multitone(Fs, interval, freq, amp, phase):
	dt = 1.0/Fs
	time = np.arange(0, interval, dt)
	x = np.zeros(len(time))

	for i in range(len(freq)):
		x += amp[i] * np.sin(2 * np.pi * freq[i] * time + phase[i])

	return time, x

def cli_error_msg():

	# error message to provide the correct command line interface (CLI) arguments
	print('Valid arguments:')
	print('\trc:  reference code')
	print('\til1: in-lab 1')
	print('\til2: in-lab 2')
	print('\tth:  take-home')
	sys.exit()

if __name__ == "__main__":

	if len(sys.argv[0:]) != 2:
		cli_error_msg()

	Fs = 100.0           # sampling rate
	Fc = 15.0            # cutoff frequency
	N_taps = 41          # number of taps for the FIR

	if (sys.argv[1] == 'rc'):  # runs the reference code (rc)

		print('Reference code for the digital filter design')

		# derive filter coefficients using firwin from SciPy:
		# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.firwin.html
		# second argument is the normalized cutoff frequency, i.e., the
		# cutoff frequency divided by Nyquist frequency (half of sampling rate)
		firwin_coeff = signal.firwin(N_taps, Fc / (Fs / 2), window = ('hann'))

		# plot the frequency response obtained through freqz
		freqzPlot(firwin_coeff, Fs, 'firwin for ' + str(int(Fc)) + ' Hz cutoff with ' + str(N_taps) + ' taps')

	elif (sys.argv[1] == 'il1'):

		print('In-lab experiment 1 for the digital filter design')

		# implement your own method for finding the coefficients for a low pass filter
		# my_own_coeff = ... provide the following arguments: Fc, Fs and N_taps
		# compare through visual inspection the frequency response against firwin
		# freqzPlot(my_own_coeff, Fs, 'my own FIR design with ' + str(N_taps) + ' taps')
		
		my_own_coeff = taps(Fc, Fs, N_taps)
		freqzPlot(my_own_coeff, Fs, str(int(Fc)) + str(N_taps))


	elif (sys.argv[1] == 'il2'):

		print('In-lab experiment 2 for the digital filter design')

		# you can confirm that a single tone has been filtered
		# filterSin(Fs, Fc, my_own_coeff)
		my_own_coeff = taps(Fc, Fs, N_taps)
		frequencies = [5, 10, 35]  
		amplitudes = [1.0, 2.0, 0.5]  
		phases = [0, np.pi / 4, np.pi / 2] 

		#generates multitone signal
		time, multitone_signal = multitone(Fs, interval=1.0, freq=frequencies, amp=amplitudes, phase=phases)

		plotTime(multitone_signal, time)
		plotSpectrum(multitone_signal, Fs, type='FFT')

		filtered_signal = signal.lfilter(my_own_coeff, 1.0, multitone_signal)

		plotTime(filtered_signal, time)
		plotSpectrum(multitone_signal, Fs, type='FFT')

	elif (sys.argv[1] == 'th'):

		print('Take-home exercise for the digital filter design')

		# for specific details check the lab document

		interval = 1.0
		frequencies = [5, 16, 25, 35]
		amplitudes = [1.0, 2.0,3.0, 4.0]
		phases = [0, np.pi/4, np.pi/2, 0]
		low_cutoff = 15.0  
		high_cutoff = 25.0  

		time, multitone_signal = multitone(Fs, interval, frequencies, amplitudes, phases)

		#original signals
		plotTime(multitone_signal, time)
		plotSpectrum(multitone_signal, Fs, type='FFT')

		#applies the bandpass filter to multitone
		bandpassfilter = signal.firwin(N_taps, [low_cutoff, high_cutoff], pass_zero=False, fs=Fs)

		#freq response for bandpass filter
		freqzPlot(bandpassfilter, Fs, 'Bandpass filter (15-25)')
		filtered_signal = signal.lfilter(bandpassfilter, 1.0, multitone_signal)

		#filtered signals
		plotTime(filtered_signal, time)
		plotSpectrum(filtered_signal, Fs, type='FFT')

	else:

		cli_error_msg()

	plt.show()
