/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include "iofunc.h"

// Function that computes the impulse response using the hanning window
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h, int upsample)
{
	auto Normcutoff = Fc/(Fs/2);
	for (auto i = 0; i < num_taps; i++) {
		if (i == (num_taps-1)/(2)){
			h[i] = Normcutoff;
		}else{
			h[i] = Normcutoff*((sin(PI*Normcutoff*(i-(num_taps-1)/(2))))/(PI*Normcutoff*(i-(num_taps-1)/(2))));
		}
		h[i] = h[i]*pow(sin(i*PI/num_taps),2);

	}
}

// Function that computes the impulse response using the hanning window
void impulseResponseBPF(float Fs, float Fb, float Fe, unsigned short int num_taps, std::vector<float> &h)
{
	auto Normcenter = ((Fe+Fb)/2)/(Fs/2);
	auto Normpass = (Fe-Fb)/(Fs/2);
	for (auto i = 0; i < num_taps; i++) {
		if (i == (num_taps-1)/2){
			h[i] = Normpass;
		}else{
			h[i] = Normpass*((sin(PI*(Normpass/2)*(i-(num_taps-1)/2)))/(PI*(Normpass/2)*(i-(num_taps-1)/2)));
		}
		h[i] = h[i]*(cos(i*PI*Normcenter));
		h[i] = h[i]*pow(sin(i*PI/num_taps),2);

	}
}

// Function that computes the convolution of x and h into y with a state variable zi (also downsamples by a constant amount)
void convolveBlockMode0(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, int dec, std::vector<float> &zi) {
	unsigned int j=0;
	unsigned int i;
	unsigned int k;
	for ( i = 0; i < x.size(); i+=dec) {
		y[j] = 0;
		for (k = 0; k < h.size(); k++) {
			if ((i-k >= 0) && ((i-k) < (x.size()))){
				y[j] += h[k] * x[i-k];
			}else{
				y[j] += zi[k-i-1] * h[k];
	    }
		}

	j++;
	}
	for (unsigned int m = 0; m<=zi.size(); m++){
		zi[m] = x[x.size()-1-m];
	}
}

// Function that computes the convolution of x and h into y with a state variable zi (also upsamples and downsamples by a constant amount)
void convolveBlockMode1(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, int expand, int dec, std::vector<float> &zi, int &prevIndex, float &yEnd) {
	unsigned int j = 0;
	unsigned int i;
	y[j] = yEnd;
	j = 1;
	for (i = (prevIndex)%(x.size()*expand); i < x.size()*expand; i+=dec) {
		y[j] = 0;
		for (unsigned int k = 0; k < h.size(); k+=expand) {
			if ((i/24)-(k/24) >= 0){
				y[j] += expand*h[(i%24)+k] * x[(i/24)-(k/24)];
			}else{
				y[j] += expand*zi[(k/24)-(i/24)] * h[(i%24)+k];
	    }
		}
		y[j] = y[j];
		j++;
  }
	yEnd = y[y.size()-1];
	prevIndex = i;

	for (unsigned int m = 0; m <= zi.size(); m++){
		zi[m] = x[x.size()-1-m];
	}
}

// Function that computes the PLL of a input signal with state saving for each block
void fmPLL(std::vector<float> pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth, std::vector<float> &ncoOut, float &integrator, float &phaseEst, float &feedbackI, float &feedbackQ, float &trigOffset){
	float Cp = 2.66;
	float Ci = 3.555;
	float Kp = normBandwidth*Cp;
	float Ki = normBandwidth*normBandwidth*Ci;

	for (int k = 0; k < pllIn.size(); k++){
		float errorI = pllIn[k] * (+feedbackI);
		float errorQ = pllIn[k] * (-feedbackQ);

		float errorD = atan2(errorQ, errorI);

		integrator = integrator + Ki*errorD;
		phaseEst = phaseEst + Kp*errorD + integrator;

		float trigArg = 2*PI*(freq/Fs)*(trigOffset+k+1) + phaseEst;
		feedbackI = cos(trigArg);
		feedbackQ = sin(trigArg);
		ncoOut[k+1] = cos(trigArg*ncoScale + phaseAdjust);
	}

}
