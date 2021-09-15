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

#
# use the four quadrant arctan function for phase detect between a pair of
# IQ samples; then unwrap the phase and take its derivative to demodulate
#
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

def customfmDemodArctan(I, Q, prev_I = 0.0, prev_Q = 0.0):
#
# the default prev_phase phase is assumed to be zero, however
# take note in block processing it must be explicitly controlled

    # empty vector to store the demodulated samples
    fm_demod = np.empty(len(I))

    # iterate through each of the I and Q pairs
    for k in range(len(I)):

        # use the atan2 function (four quadrant version) to detect angle between
        # the imaginary part (quadrature Q) and the real part (in-phase I)
        current_I = I[k]
        current_Q = Q[k]

        # take the derivative of the phase
        if (current_I != 0 and current_Q != 0):
            fm_demod[k] = (current_I * (current_Q - prev_Q) - current_Q * (current_I - prev_I))/ (current_I**2 + current_Q**2)
        else:
            fm_demod[k] = 0
        # save the state of the current phas
        # to compute the next derivative
        prev_I = current_I
        prev_Q = current_Q

    # return both the demodulated samples as well as the last phase
    # (the last phase is needed to enable continuity for block processing)
    return fm_demod, prev_I, prev_Q

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
# this is less accurate than the Welch method from matplotlib
# however, as the visual inspections confirm, the estimate gives
# the user a "reasonably good" view of the power spectrum
def estimatePSD(samples, NFFT, Fs):

    # rename the NFFT argument (notation consistent with matplotlib.psd)
    # to freq_bins (i.e., frequency bins for which we compute the spectrum)
    freq_bins = NFFT
    # frequency increment (or resolution)
    df = Fs/freq_bins

    # create the frequency vector to be used on the X axis
    # for plotting the PSD on the Y axis (only positive freq)
    freq = np.arange(0, Fs/2, df)

    # design the Hann window used to smoothen the discrete data in order
    # to reduce the spectral leakage after the Fourier transform
    hann = np.empty(freq_bins)
    for i in range(len(hann)):
        hann[i] = pow(math.sin(i*math.pi/freq_bins),2)

    # create an empty list where the PSD for each segment is computed
    psd_list = []

    # samples should be a multiple of frequency bins, so
    # the number of segments used for estimation is an integer
    # note: for this to work you must provide an argument for the
    # number of frequency bins not greater than the number of samples!
    no_segments = int(math.floor(len(samples)/float(freq_bins)))

    # iterate through all the segments
    for k in range(no_segments):

        # apply the hann window (using pointwise multiplication)
        # before computing the Fourier transform on a segment
        windowed_samples = samples[k*freq_bins:(k+1)*freq_bins] * hann

        # compute the Fourier transform using the built-in FFT from numpy
        Xf = np.fft.fft(windowed_samples, freq_bins)

        # note, you can check how MUCH slower is DFT vs FFT by replacing the
        # above function call with the one that is commented below
        #
        # Xf = DFT(windowed_samples)
        #
        # note, the slow impelementation of the Fourier transform is not as
        # critical when computing a static power spectra when troubleshooting
        # note: time permitting a custom FFT can be implemented

        # since input is real, we keep only the positive half of the spectrum
        # however, we will also add the signal energy of negative frequencies
        # to have a better a more accurate PSD estimate when plotting
        Xf = Xf[0:int(freq_bins/2)] # keep only positive freq bins
        psd_seg = 1/(Fs*freq_bins/2) * abs(Xf)**2 # compute signal power
        psd_seg = 2*psd_seg # add the energy from the negative freq bins

        # translate to the decibel (dB) scale
        for i in range(len(psd_seg)):
            psd_seg[i] = 10*math.log10(psd_seg[i])

        # append to the list where PSD for each segment is stored
        # in sequential order (first segment, followed by the second one, ...)
        psd_list.extend(psd_seg)

    # compute the estimate to be returned by the function through averaging
    psd_est = np.zeros(int(freq_bins/2))

    # iterate through all the frequency bins (positive freq only)
    # from all segments and average them (one bin at a time ...)
    for k in range(int(freq_bins/2)):
        # iterate through all the segments
        for l in range(no_segments):
            psd_est[k] += psd_list[k + l*int(freq_bins/2)]
        # compute the estimate for each bin
        psd_est[k] = psd_est[k] / no_segments

    # the frequency vector and PSD estimate
    return freq, psd_est

# A manual low pass function
def manualLPIR(Fs, Fc, Ntaps):

    Normcutoff = Fc/(Fs/2)
    h = np.zeros(Ntaps, dtype=np.double)
    for i in range(Ntaps):
        if (i == (Ntaps-1)/2):
            h[i] = Normcutoff
        else:
            h[i] = Normcutoff*((np.sin(np.pi*Normcutoff*(i-(Ntaps-1)/2)))/(np.pi*Normcutoff*(i-(Ntaps-1)/2)))

        h[i] = h[i]*((np.sin((i*np.pi)/Ntaps))**2)

    return h

def manualBPIR(Fs, Fb, Fe, Ntaps):
    Normcenter = ((Fe+Fb)/2)/(Fs/2)
    Normpass = (Fe-Fb)/(Fs/2)
    h = np.zeros(Ntaps, dtype=np.double)
    for i in range(Ntaps):
        if (i == (Ntaps-1)/2):
            h[i] = Normpass
        else:
            h[i] = Normpass*((np.sin(np.pi*(Normpass/2)*(i-(Ntaps-1)/2)))/(np.pi*(Normpass/2)*(i-(Ntaps-1)/2)))
        h[i] = h[i]*(np.cos(i*np.pi*Normcenter))
        h[i] = h[i]*((np.sin((i*np.pi)/Ntaps))**2)

    return h

# A manual linear filtering function (convolution)
def manualLFilter(x, h):
    y = np.zeros(len(x) + len(h) - 1, dtype=np.double)

    for i in range(len(y)):
        for j in range(len(h)):
            if (i-j >= 0 and i-j < len(x)):
                y[i] = y[i] + h[j]*x[i-j]
    return y

# A manual blocking Linear filtering function (convolution)
def manualBlockLFilter(x, h, zi):
    y = np.zeros(len(x), dtype=np.double)

    for i in range(len(y)):
        for j in range(len(h)):
            if (i >= j):
                y[i] = y[i] + h[j]*x[i-j]
            else:
                y[i] = y[i] + h[j]*zi[i-j]

    zf = x[-(len(h)-1):]

    return y, zf

# def manualBlockLFilter(x, h, zi):
#
#     conv = np.zeros(len(x))
#     for i in range(len(x)):  # calculate y[i]
#         conv[i] = 0
#         zi_ctrl = 0
#         for j in range(len(h)):
#             if(j>i):
#                 conv[i] += zi[len(zi) - zi_ctrl - 1] * h[j]
#             else:
#                 conv[i] += x[i - j] * h[j]
#             zi_ctrl+=1
#     zi = x[-len(zi): ] #Create new zi that stores the values of x in reverse (thus -len(zi) term)
#     return conv, zi

#def manualBlockLFilter(x, h, zi):
#    yi = np.zeros(len(x) + len(h) - 1, dtype=np.double)
#    y = np.zeros(len(x), dtype=np.double)
#    zf = np.zeros(len(h) - 1, dtype=np.double)
#    xzf = np.concatenate((zi, x), axis = None)

#    for i in range(len(xzf)):
#        for j in range(len(h)):
#            if (i >= j):
#                yi[i] = yi[i] + h[j]*xzf[i-j]
#            else:
#                yi[i] = 0

#    for i in range(len(x)):
#        y[i] = yi[i]

#    for i in range(len(x), len(yi)):
#        zf[i - len(x)] = yi[i]

#    return y, zf

def downsampler(N, Arr):
    ArrNew = np.zeros(math.ceil(len(Arr)/N), dtype=np.double)
    count = 0
    for i in range(len(Arr)):
        if (i%N == 0):
            ArrNew[count] = Arr[i]
            count += 1

    return ArrNew

if __name__ == "__main__":

    # do nothing when this module is launched on its own
    pass
