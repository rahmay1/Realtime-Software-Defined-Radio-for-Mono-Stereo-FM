#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.animation as animation
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math
import sys

# use "custom" estimatePSD and fmDemodArctan
from fmSupportLib import estimatePSD, fmDemodArctan

rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 151
rf_decim = 10

audio_Fs = 48e3
audio_decim = 5
# add other settings for audio, like filter taps, ...

# we need a dummy init for animate to avoid calling
# the main update function (animate_update) twice at the start
def animate_init():
	pass

# this is the main animation function called by "animation.FuncAnimation"
# the first argument is mandatory and it keeps track of how many times
# this function has been called
# the subsequent arguments are custom to this particular application
def animate_update(block_count, iq_data, block_size, rf_coeff, audio_coeff):

	global audio_data, state_i_lpf_100k, state_q_lpf_100k, state_phase

	if (block_count+1)*block_size > len(iq_data):
		print('Finished processing the raw I/Q samples')
		# write audio data to file (assumes audio_data samples are -1 to +1)
		wavfile.write("../data/fmMonoAnim.wav", int(audio_Fs), np.int16((audio_data/2)*32767))
		sys.exit()

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
	fm_demod, state_phase = fmDemodArctan(i_ds, q_ds, state_phase)

	# PSD after FM demodulation
	ax0.clear()
	ax0.psd(fm_demod, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
	# freq, my_psd = estimatePSD(fm_demod, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
	# ax0.plot(freq, my_psd)
	ax0.set_ylabel('PSD (dB/Hz)')
	ax0.set_xlabel('Freq (kHz)')
	ax0.set_title('Demodulated FM (block ' + str(block_count) + ')')

	# extract the mono audtio data through filtering
	# audio_filt = ... change as needed

	# PSD after extracting mono audio
	# ax1.clear()
	# ax1.psd(audio_filt, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
	# ax1.set_ylabel('PSD (dB/Hz)')
	# ax1.set_xlabel('Freq (kHz)')
	# ax1.set_title('Extracted Mono')

	# downsample audio data
	# audio_block =  ... change as needed

	# PSD after decimating mono audio
	# ax2.clear()
	# ax2.psd(audio_block, NFFT=512, Fs=audio_Fs/1e3)
	# ax2.set_ylabel('PSD (dB/Hz)')
	# ax2.set_xlabel('Freq (kHz)')
	# ax2.set_title('Mono Audio')

	# concatanete most recently processed audio_block
	# to the previous blocks stored in audio_data
	# audio_data = np.concatenate((audio_data, audio_block))

if __name__ == "__main__":

	# read the raw IQ data from the recorded file
	# IQ data is normalized between -1 and +1 and interleaved
	in_fname = "../data/iq_samples.raw"
	iq_data = np.fromfile(in_fname, dtype='float32')
	print("Read raw RF data from \"" + in_fname + "\" in float32 format")

	# coefficients for the front-end low-pass filter
	rf_coeff = signal.firwin(rf_taps, \
							rf_Fc/(rf_Fs/2), \
							window=('hann'))

	# coefficients for the filter to extract mono audio
	audio_coeff = np.array([]) # to be changed by you

	# set up drawing
	fig, (ax0, ax1, ax2) = plt.subplots(nrows=3)
	fig.subplots_adjust(hspace = 1.0)

	# select a block_size that is in KB and
	# a multiple of decimation factors
	block_size = 1024 * rf_decim * audio_decim * 2
	block_count = 0

	# states needed for continuity in block processing
	state_i_lpf_100k = np.zeros(rf_taps-1)
	state_q_lpf_100k = np.zeros(rf_taps-1)
	state_phase = 0
	# add state as needed for the mono channel filter

	# audio buffer that stores all the audio blocks
	audio_data = np.array([])

	try:
		# calls the animation function (animate_update) repeatedly
		# check matplotlib documentation for further details
		ani = animation.FuncAnimation(fig, animate_update, \
						interval=50, init_func=animate_init, \
						fargs=(iq_data, block_size, rf_coeff, audio_coeff, ))
		plt.show()

	except KeyboardInterrupt:
		pass
