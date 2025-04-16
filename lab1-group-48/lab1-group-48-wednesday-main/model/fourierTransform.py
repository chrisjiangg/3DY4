#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

from scipy.signal import square
import matplotlib.pyplot as plt
import numpy as np
import cmath, math
import sys

def plotSpectrum(x, Fs, type = 'FFT'):

	n = len(x)             # length of the signal
	df = Fs / n            # frequency increment (width of freq bin)

	# compute Fourier transform, its magnitude and normalize it before plotting
	if type == 'FFT':
		Xfreq = np.fft.fft(x)
	elif type == 'dft':
		Xfreq = dft(x)

	XMag = abs(Xfreq) / n

	# Note: because x is real, we keep only the positive half of the spectrum
	# Note also: half of the energy is in the negative half (not plotted)
	XMag = XMag[0 : int(n / 2)]

	# freq vector up to Nyquist freq (half of the sample rate)
	freq = np.arange(0, Fs / 2, df)

	fig, ax = plt.subplots()
	ax.plot(freq, XMag)
	ax.set(xlabel = 'Frequency (Hz)', ylabel = 'Magnitude',
		   title = 'Frequency domain plot')
	# fig.savefig("freq.png")
	plt.show()

def dft(x):
    N = len(x) #length of signal
    X = []  #list to store freq domain representation
    for m in range(N):  
        sum_val = 0 + 0j  #initialize sum as complex number
        for k in range(N):  
            sum_val += x[k] * cmath.exp(-2j * cmath.pi * m * k / N)  
        X.append(sum_val) #appends sum_val to freq domain 
    return np.array(X)  #numpy array


def idft(X):
    N = len(X)  
    x = []  
    for k in range(N):  
        sum_val = 0 + 0j  
        for m in range(N):  
            sum_val += X[m] * cmath.exp(2j * cmath.pi * m * k / N)  
        x.append(sum_val.real / N)  
    return np.array(x)  

def plotTime(x, time):

	fig, ax = plt.subplots()
	ax.plot(time, x)
	ax.set(xlabel = 'Time (sec)', ylabel = 'Amplitude',
		   title = 'Time domain plot')
	# fig.savefig("time.png")
	plt.show()

def generateSin(Fs, interval, frequency = 7.0, amplitude = 5.0, phase = 0.0):

	dt = 1.0 / Fs                          # sampling period (increment in time)
	time = np.arange(0, interval, dt)      # time vector over interval

	# generate the sin signal
	x = amplitude * np.sin(2 * math.pi * frequency * time + phase)

	return time, x

def generatesquare(Fs, interval, frequency = 7.0, amplitude = 5.0, phase = 0.0):
	dt = 1.0 / Fs                          # sampling period (increment in time)
	time = np.arange(0, interval, dt)      # time vector over interval

	# generate the square signal
	x = amplitude * square(2 * math.pi * frequency * time, duty_cycle)

	return time, x

def fourierUnitTest(min = -1, max = 1, N = 32):

	# use a NumPy function to generate a random sequence of length "size"
	# of uniformly distributed numbers between min and max
	x = np.random.uniform(low = min, high = max, size = N)

	# the reconstructed signal (rx) must be the inverse Fourier transform of the frequency domain representation of x
	# below is a placeholder that helps us set up the unit test - eventually rx must be the real part of IDFT(DFT(x))
	rx = np.zeros(len(x), dtype = x.dtype)
	
	Xfreq = dft(x)
	rx = idft(Xfreq)

	# check if all the values are close to each within the given tolerance
	if not np.allclose(x, rx, atol = 1e-4):
		print(f"Comparison between original signal and the inverse Fourier transform of its frequency domain representation fails")
		print(f"Original signal:", x)
		print(f"Reconstructed signal:", rx)
		print("Maximum difference:", np.max(np.abs(x - rx)))
	else:
		print(f"Unit test for Fourier/Inverse Fourier transform passed.")

def cli_error_msg():

	# error message to provide the correct command line interface (CLI) arguments
	print('Valid arguments:')
	print('\trc:  reference code')
	print('\til1: in-lab 1')
	print('\til2: in-lab 2')
	print('\til3: in-lab 3')
	print('\tth:  take-home')
	sys.exit()

if __name__ == "__main__":

	if len(sys.argv[0:]) != 2:
		cli_error_msg()

	Fs = 100.0          # sampling rate
	interval = 1.0      # set up to one full second

	if (sys.argv[1] == 'rc'): # runs the reference code (rc)

		print('Reference code for the Fourier transform')

		# generate the user-defined sin function
		time, x = generateSin(Fs, interval)
		# plot the signal in time domain
		plotTime(x, time)
		# plot the signal in frequency domain
		plotSpectrum(x, Fs, type = 'FFT')

	elif (sys.argv[1] == 'il1'):

		print('In-lab experiment 1 for the Fourier transform')

		# compute the spectrum with your own DFT
		# you can use cmath.exp() for complex exponentials
		# plotSpectrum(x, Fs, type = 'your DFT name')
		time, x = generateSin(Fs, interval)
		Xfreq = dft(x)
		reconstructed_x = idft(Xfreq)
		plotTime(reconstructed_x, time)
		plotSpectrum(x, Fs, type='dft')
		fourierUnitTest()
		# confirm DFT/IDFT correctness by checking if x == IDFT(DFT(x))
		# the self check should be done within the fourierUnitTest function

		# for further details, if any, check the lab document

	elif (sys.argv[1] == 'il2'):

		print('In-lab experiment 2 for the Fourier transform')

		# you can use np.random.randn() to generate a random parameter
		# we can overwrite the default values
		# frequency = 8.0                      # frequency of the signal
		# amplitude = 3.0                      # amplitude of the signal
		# phase = 1.0                          # phase of the signal
		# time, x = generateSin(Fs, interval, frequency, amplitude, phase)

		# You should also numerically check if the signal energy
		# in time and frequency domains is identical
		N = 1000
		x = np.random.uniform(-10,10,size = N)

		E_time = np.sum(np.abs(x)**2)

		Xfreq = dft(x)
		E_freq = np.sum(np.abs(Xfreq)**2) / N


		# create your own Unit Test to validate Parseval's theorem
		# (based on the same principle as for the Fourier transform)
		#checker
		print(f"Energy in Time Domain: {E_time}")
		print(f"Energy in Frequency Domain: {E_freq}")
		if np.isclose(E_time, E_freq, atol=1e-4):
			print("match")
		else:
			print("no match")


		Fs = 100.0  
		time = np.linspace(0, (N - 1) / Fs, N) 
		plotTime(x, time)  
		plotSpectrum(x, Fs, type='dft')  
		# for further details, if any, check the lab document

	elif (sys.argv[1] == 'il3'):

		print('In-lab experiment 3 for the Fourier transform')

		# generate randomized multi-tone signals
		# plot them in both time and frequency domain
		Fs = 1000

		time, x1 = generateSin(Fs, interval, 50, 1.0, 0)
		time, x2 = generateSin(Fs, interval, 120, 0.5, np.pi / 4)
		time, x3 = generateSin(Fs, interval, 200, 1.0, np.pi / 2)
		TotalTone = x1 + x2 + x3


		#plotTime(x1, time)
		#plotTime(x2, time)
		#plotTime(x3, time)
		plotTime(TotalTone, time)
		plotSpectrum(TotalTone, Fs, type='FFT')

		#checker
		E_time = np.sum(np.abs(TotalTone)**2)
		Xfreq = dft(TotalTone) 
		E_freq = np.sum(np.abs(Xfreq)**2) / len(TotalTone)
		if np.isclose(E_time, E_freq, atol=1e-4):
			print("match")
		else:
			print("no match")


		# for further details, if any, check the lab document

	elif (sys.argv[1] == 'th'):

		print('Take-home exercise for the Fourier transform')

		# for specific details check the lab document
		Fs = 100.0 
		interval = 1.0  
		frequency = 7.0 
		amplitude = 5.0  
		duty_cycle = 0.25  

		time, square_wave = generatesquare(Fs, interval, frequency, amplitude, duty_cycle)
		plotTime(square_wave, time)
		plotSpectrum(square_wave, Fs, type='FFT')

	else:

		cli_error_msg()

	plt.show()
