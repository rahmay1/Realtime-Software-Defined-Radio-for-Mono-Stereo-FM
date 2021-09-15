#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import numpy as np
import math

def fmPll(pllIn, freq, Fs, \
		ncoScale, phaseAdjust, normBandwidth, ncoOut, integrator, phaseEst, feedbackI, feedbackQ, trigOffset):

	"""
	pllIn 	 		array of floats
					input signal to the PLL (assume known frequency)

	freq 			float
					reference frequency to which the PLL locks (19kHz)

	Fs  			float
					sampling rate for the input/output signals

	ncoScale		float
					frequency scale factor for the NCO output (2 for stereo, or 0.5 for RDS)

	phaseAdjust		float
					phase adjust to be added to the NCO only (For RDS only)

	normBandwidth	float
					normalized bandwidth for the loop filter (keep as default I think)
					(relative to the sampling rate)

	state 			to be added

	"""

	# print("During PLL: ")
	# print("Integrator: " +str(integrator))
	# print("phaseEst: " +str(phaseEst))
	# print("feedbackI: " +str(feedbackI))
	# print("feedbackQ: " +str(feedbackQ))
	# print("trigOffset: " +str(trigOffset))

	# scale factors for proportional/integrator terms
	# these scale factors were derived assuming the following:
	# damping factor of 0.707 (1 over square root of 2)
	# there is no oscillator gain and no phase detector gain
	Cp = 2.666
	Ci = 3.555

	# gain for the proportional term
	Kp = (normBandwidth)*Cp
	# gain for the integrator term
	Ki = (normBandwidth*normBandwidth)*Ci

	# output array for the NCO
	#ncoOut = np.empty(len(pllIn)+1)

	# initialize internal state
	# integrator = 0.0
	# phaseEst = 0.0
	# feedbackI = 1.0
	# feedbackQ = 0.0
	# ncoOut[0] = 1.0
	# trigOffset = 0
	# note: state saving will be needed for block processing

	# if (len(pllIn) != len(ncoOut)):
	# 	print("LOSING VALUES! PLLIN IS NOT THE SAME AS NCOOUT")

	for k in range(len(pllIn)):

		# phase detector
		errorI = pllIn[k] * (+feedbackI)  # complex conjugate of the
		errorQ = pllIn[k] * (-feedbackQ)  # feedback complex exponential

		# four-quadrant arctangent discriminator for phase error detection
		errorD = math.atan2(errorQ, errorI)

		# loop filter
		integrator = integrator + Ki*errorD

		# update phase estimate
		phaseEst = phaseEst + Kp*errorD + integrator

		# internal oscillator
		trigArg = 2*math.pi*(freq/Fs)*(trigOffset+k+1) + phaseEst
		feedbackI = math.cos(trigArg)
		feedbackQ = math.sin(trigArg)
		ncoOut[k+1] = math.cos(trigArg*ncoScale + phaseAdjust)


	# for stereo only the in-phase NCO component should be returned
	# for block processing you should also return the state
	return ncoOut, integrator, phaseEst, feedbackI, feedbackQ, trigOffset
	# for RDS add also the quadrature NCO component to the output

if __name__ == "__main__":

	pass
