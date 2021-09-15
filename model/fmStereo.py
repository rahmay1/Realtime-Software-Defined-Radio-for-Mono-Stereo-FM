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

# use "custom" fmDemodArctan
from fmSupportLib import fmDemodArctan
from fmSupportLib import customfmDemodArctan
from fmSupportLib import manualLPIR
from fmSupportLib import manualBPIR
#from fmSupportLib import manualLFilter
from fmSupportLib import downsampler
from fmSupportLib import manualBlockLFilter
from fmPll import fmPll

rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 151
rf_decim = 10

audio_Fs = 240e3
audio_Fc = 16e3
audio_taps = 151
audio_decim = 5
# add other settings for audio, like filter taps, ...

if __name__ == "__main__":

	# read the raw IQ data from the recorded file
	# IQ data is assumed to be in 8-bits unsigned (and interleaved)
	in_fname = "../data/my_samples_u8.raw"
	raw_data = np.fromfile(in_fname, dtype='uint8')
	# IQ data is normalized between -1 and +1
	iq_data = (raw_data - 128.0)/128.0
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")

	# coefficients for the front-end low-pass filter
	rf_coeff = signal.firwin(rf_taps, \
							rf_Fc/(rf_Fs/2), \
							window=('hann'))

	# coefficients for the filter to extract the mono audio
	audio_coeff = signal.firwin(audio_taps, \
		audio_Fc/((rf_Fs/rf_decim)/2), window=('hann'))

	# coefficients for the filter to extract stereo channel audio
	audio_coeffS = manualBPIR(audio_Fs, 22e3, 54e3, audio_taps)

	# coefficients for the filter to extract pilot tone audio
	audio_coeffP = manualBPIR(audio_Fs, 18.5e3, 19.5e3, audio_taps)

	# set up drawing
	fig, (ax0, ax1, ax2, ax3, ax4) = plt.subplots(nrows=5,figsize=(30,15))
	fig.subplots_adjust(hspace = 1.0)

	# select a block_size that is in KB and
	# a multiple of decimation factors
	block_size = 1024 * rf_decim * audio_decim * 2
	block_count = 0

	# states needed for continuity in block processing
	state_i_lpf_100k = np.zeros(rf_taps-1)
	state_q_lpf_100k = np.zeros(rf_taps-1)
	state_lpf_16k = np.zeros(audio_taps-1)
	state_lpf_16k2 = np.zeros(audio_taps-1)
	state_I = 0
	state_Q = 0
	# add state as needed for the mono channel filter
	state_monoS = np.zeros(audio_taps-1)
	state_monoP = np.zeros(audio_taps-1)

	# audio buffer that stores all the audio blocks
	left_audio = np.array([])
	right_audio = np.array([])

	# Stores the stereo data for each block
	stereo_data_block = np.array([])

	# State Initialize for PLL

	integrator = 0.0
	phaseEst = 0.0
	feedbackI = 1.0
	feedbackQ = 0.0
	trigOffset = 0

	# if the number of samples in the last block is less than the block size
	# it is fine to ignore the last few samples from the raw IQ file
	while (block_count+1)*block_size < len(iq_data):
		# if you wish to have shorter runtimes while troubleshooting
		# you can control the above loop exit condition as you see fit

		print('Processing block ' + str(block_count))

		# filter to extract the FM channel (I samples are even, Q samples are odd)
		i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
				iq_data[(block_count)*block_size:(block_count+1)*block_size:2],
				zi=state_i_lpf_100k)
		q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
				iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2],
				zi=state_q_lpf_100k)

		# downsample the FM channel
		i_ds = i_filt[::rf_decim]
		q_ds = q_filt[::rf_decim]

		# FM demodulator
		fm_demod, state_I, state_Q = customfmDemodArctan(i_ds, q_ds, state_I, state_Q)

		# extract the stereo channel with a bandpass filter
		#audio_stereo, state_monoS = manualBlockLFilter(fm_demod, audio_coeffS, state_monoS)

		audio_stereo, state_monoS = signal.lfilter(audio_coeffS, 1.0, \
				fm_demod, zi=state_monoS)

		# extract the stereo channel with a bandpass filter
		#audio_pilot, state_monoP = manualBlockLFilter(fm_demod, audio_coeffP, state_monoP)

		audio_pilot, state_monoP = signal.lfilter(audio_coeffP, 1.0, \
				fm_demod, zi=state_monoP)

		# Defining the ncoOut state variable
		if (block_count == 0):
			ncoOut = np.empty(len(audio_pilot)+1)
			ncoOut[0] = 1.0
		else:
			ncoOut[0] = ncoOut[len(ncoOut)-1];

		# Calling the PLL
		ncoOut, integrator, phaseEst, feedbackI, feedbackQ, trigOffset = fmPll(audio_pilot, 19e3, audio_Fs, 2, 0.0, 0.01, ncoOut, integrator, phaseEst, feedbackI, feedbackQ, trigOffset)
		trigOffset+=len(audio_pilot)
		ncoOut[0] = ncoOut[len(ncoOut)-1];

		#Stereo Mixer
		stereo_data_block = np.multiply(ncoOut[0:len(ncoOut)-1],audio_stereo)

		# Digital Filtering for stereo signal (assuming mono 0)
		audio_filt_stereo, state_lpf_16k2 = signal.lfilter(audio_coeff, 1.0, \
				stereo_data_block, zi=state_lpf_16k2)
		audio_block_stereo = audio_filt_stereo[::audio_decim]

		# downsample audio data (mono 0 operation)
		audio_filt_mono, state_lpf_16k = signal.lfilter(audio_coeff, 1.0, \
				fm_demod, zi=state_lpf_16k)
		audio_block_mono = audio_filt_mono[::audio_decim]

		left_audio_block = np.add(audio_block_stereo, audio_block_mono)
		right_audio_block = np.subtract(audio_block_mono, audio_block_stereo)
		left_audio_block = left_audio_block/2
		right_audio_block = right_audio_block/2

		left_audio = np.concatenate((left_audio, left_audio_block))
		right_audio = np.concatenate((right_audio, right_audio_block))

		# concatanete most recently processed audio_block
		# to the previous blocks stored in audio_data
		# audio_data = np.concatenate((audio_data, audio_filt))

		# to save runtime select the range of blocks to log iq_data
		# this includes both saving binary files as well plotting PSD
		# below we assume we want to plot for graphs for blocks 10 and 11
		if (block_count >= 10 and block_count < 12) or block_count == 0 or block_count == 1:

			# PSD after FM demodulation
			ax0.clear()
			ax0.psd(fm_demod, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
			ax0.set_ylabel('PSD (dB/Hz)')
			ax0.set_xlabel('Freq (kHz)')
			ax0.set_title('Demodulated FM (block ' + str(block_count) + ')')

			# output binary file name (where samples are written from Python)
			fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
			# create binary file where each sample is a 32-bit float
			fm_demod.astype('float32').tofile(fm_demod_fname)

			# PSD after extracting stereo audio
			ax1.clear()
			ax1.psd(audio_stereo, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
			ax1.set_ylabel('PSD (dB/Hz)')
			ax1.set_xlabel('Freq (kHz)')
			ax1.set_title('Stereo Channel Extraction')

			# PSD after extracting pilot tone audio
			ax2.clear()
			ax2.psd(audio_pilot, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
			ax2.set_ylabel('PSD (dB/Hz)')
			ax2.set_xlabel('Freq (kHz)')
			ax2.set_title('Pilot Tone Extraction')

			# PSD for ncoOut Signal for the irst few blocks
			# ax3.clear()
			# #ax2.psd(audio_pilot, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
			# ax3.plot(ncoOut[0:40])
			# ax3.set_ylabel('NCO Output values')
			# ax3.set_xlabel('Samples')
			# ax3.set_title('ncoOut Signal(Beginning of Block)')

			# PSD for left stereo output
			ax3.clear()
			ax3.psd(left_audio_block, NFFT=512, Fs=48e3/1e3)
			ax3.set_ylabel('PSD (dB/Hz)')
			ax3.set_xlabel('Freq (kHz)')
			ax3.set_title('Stereo Left Data')

			# PSD for ncoOut Signal
			# ax4.clear()
			# #ax0.psd(ncoOut, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
			# ax4.plot(ncoOut[len(ncoOut) - 40:len(ncoOut)])
			# ax4.set_ylabel('NCO Output values')
			# ax4.set_xlabel('Samples')
			# ax4.set_title('ncoOut Signal(End of Block)')

			# PSD for right stereo output
			ax4.clear()
			ax4.psd(right_audio_block, NFFT=512, Fs=48e3/1e3)
			ax4.set_ylabel('PSD (dB/Hz)')
			ax4.set_xlabel('Freq (kHz)')
			ax4.set_title('Stereo Right Data')

			# save figure to file
			fig.savefig("../data/fmStereoBlock" + str(block_count) + ".png")

		block_count += 1

	print('Finished processing the raw I/Q samples')

	audio_output = np.vstack((left_audio, right_audio)).T

	# write audio data to a .wav file (assumes audio_data samples are -1 to +1) Save as a 2 dimensional numpy array (Search up)
	#wavfile.write("../data/fmStereoLeft.wav", int(48e3), np.int16((left_audio/2)*32767))
	#wavfile.write("../data/fmStereoRight.wav", int(48e3), np.int16((right_audio/2)*32767))
	wavfile.write("../data/fmStereo.wav", int(48e3), np.int16((audio_output/2)*32767))

	# uncomment assuming you wish to show some plots
	plt.show()
