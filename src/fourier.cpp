/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

// source code for Fourier-family of functions
#include "dy4.h"
#include "fourier.h"

// Function that computes the Discrete Fourier Transform
void DFT(const std::vector<float> &x, std::vector<std::complex<float>> &Xf) {
	Xf.resize(x.size(), static_cast<std::complex<float>>(0, 0));
	for (auto m = 0; m < Xf.size(); m++) {
		for (auto k = 0; k < x.size(); k++) {
				std::complex<float> expval(0, -2*PI*(k*m) / x.size());
				Xf[m] += x[k] * std::exp(expval);
		}
	}
}

// Function that computes the magnitude values in a complex vector
void computeVectorMagnitude(const std::vector<std::complex<float>> &Xf, std::vector<float> &Xmag)
{
	// only the positive frequencies
	Xmag.resize(Xf.size(), static_cast<float>(0));
  for (auto i = 0; i < Xf.size(); i++) {
    Xmag[i] = std::abs(Xf[i])/Xf.size();
  }
}

// Function that computes the estimation of Power Spectral Density of an input vector
void estimatePSD(const std::vector<float> &samples, float Fs, std::vector<float> &psd_est, std::vector<float> &freq)
{

	int freq_bins = NFFT;

	float df = Fs/freq_bins;

	freq.resize((Fs/2)/df, static_cast<float>(0));

	for (int i = 1; i < (Fs/2)/df; i++){
		freq[i] = freq[i-1] + df;
	}

	std::vector<float> hann;
	hann.resize(freq_bins,static_cast<float>(0));

	for (int i = 0; i < hann.size(); i++){
		hann[i] = pow(sin(i*PI/freq_bins),2);
	}

	std::vector<float> psd_list;
	int no_segments = (int)std::floor(samples.size()/(float)freq_bins);

	for (int k = 0; k < no_segments; k++){

		std::vector<float> windowed_samples;
		windowed_samples.resize((k+1)*freq_bins - (k*freq_bins), static_cast<float>(0));
		std::vector<float> newSamples;
		newSamples.resize((k+1)*freq_bins - (k*freq_bins), static_cast<float>(0));
		newSamples = std::vector<float>(std::begin(samples) + (k*freq_bins), std::begin(samples) + (k+1)*freq_bins);


		for (int i = 0; i < windowed_samples.size(); i++){
			windowed_samples[i] = newSamples[i] * hann[i];
		}

		std::vector<std::complex<float>> Xf;
		DFT(windowed_samples, Xf);

		Xf = std::vector<std::complex<float>>(std::begin(Xf), std::begin(Xf) + (int)freq_bins/2);

		std::vector<float> psd_seg;
		psd_seg.resize(Xf.size(), static_cast<float>(0));

		for (int i = 0; i < psd_seg.size(); i++){
			psd_seg[i] = pow(abs(Xf[i]),2) * 1/(Fs*freq_bins/2);
		}

		for (int i = 0; i < psd_seg.size(); i++){
			psd_seg[i] = 2*psd_seg[i];
		}


		for (int i = 0; i < psd_seg.size(); i++){
			psd_seg[i] = 10*log10(psd_seg[i]);
		}

		psd_list.insert(std::end(psd_list), std::begin(psd_seg), std::end(psd_seg));

	}

	psd_est.resize((int)freq_bins/2, static_cast<float>(0));

	for (int k = 0; k < (int)freq_bins/2; k++){
		for (int l = 0; l < no_segments; l++){
			if (k + l*((int)(freq_bins/2)) < psd_list.size()){
				psd_est[k] += psd_list[k + l*((int)(freq_bins/2))];
			}
		}
		psd_est[k] = psd_est[k] / no_segments;
	}

}
