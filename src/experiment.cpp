/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "stdlib.h"
#include <iostream>
#include <fstream>
#include <chrono>

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Handles reading in information form fdin data stream file.
void readStdinBlockData(unsigned int num_samples, unsigned int block_id, std::vector<float> &i_data, std::vector<float> &q_data){
	std::vector<char> raw_data(num_samples);
	std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));
	for(unsigned int k = 0; k < num_samples/2; k++){
		i_data[k] = float(((unsigned char)raw_data[k*2] - 128)/128.0);
		q_data[k] = float(((unsigned char)raw_data[(k*2)+1] - 128)/128.0);

	}
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Handles the FmDemod Stuff:
void fmDemod(std::vector<float> &I, std::vector<float> &Q, float& prev_I, float& prev_Q, std::vector<float> &fm_demod){
	for (int k=0; k < I.size(); k++){
		if (I[k] != 0 and Q[k] != 0){
			fm_demod[k] = (I[k] * (Q[k] - prev_Q) - Q[k] * (I[k] - prev_I))/ (pow(I[k],2) + pow(Q[k],2));
		}else{
			fm_demod[k] = 0;
		}
		prev_I = I[k];
		prev_Q = Q[k];
	}
}

int main(int argc, char *argv[])
{
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Create Mode control variable, Default to mode.
	int mode = 0;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Mode Selection process
	if(argc < 2){
		std::cerr << "Operating in default mode 0" << std::endl;
	}else if(argc == 2){
		mode = atoi(argv[1]);
		if (mode == 1){
			std::cerr << "Operating in mode 1" << std::endl;
		}else if (mode == 0){
			std::cerr << "Operating in mode 0" << std::endl;
		}
		else{
			std::cerr << "Wrong mode " << mode << std::endl;
			exit(1);
		}
	}else{
		std::cerr << "Usage: " << argv[0] << std::endl;
		std::cerr << "or " << std::endl;
		std::cerr << "Usage: " << argv[0] << " 1" << std::endl;
	}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Initialization of variables and vectors
// Set up rf parameters for filtering 100kHz carrier frequency. Signal frequency is 2.4MHz in Mode 0 and 2.5 MHz in Mode 1
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Front end Path Parameters
	int rf_Fc = 100e3;
	int rf_Fs = 2.4e6;
	if (mode == 1){
		int rf_Fs = 2.5e6;
	}
	int rf_taps = 151;
	int rf_decim = 10;
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Mono Path Parameters
	int mono_Fc = 16e3;
	int mono_Fs;
	int mono_taps = 151;
	int mono_decim;
	int mono_upscale;
	int mono_prevIndex = 0;
	int stereo_prevIndex = 0;
	float mono_yEnd = 0.0;
	float stereo_yEnd = 0.0;
	if (mode == 1){
		mono_decim = 125;
		mono_upscale = 24;
		mono_Fs = 250e3;
	}else{
		mono_decim = 5;
		mono_upscale = 1;
		mono_Fs = 240e3;
	}

	// Block specifications aswell as states for the I and Q demodulator
	int block_index = 0;
	int block_size = 1024*rf_decim*(mono_decim/mono_upscale)*2;
	float state_I = 0.0;
	float state_Q = 0.0;
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // Initialize vectors for storing audio data and signal I/Q samples.
	std::vector<float> i_data(block_size/2,0.0); // Float Vector for In-Phase Data
	std::vector<float> q_data(block_size/2,0.0); // Float Vector for Quadrant
	std::vector<float> zi_i(rf_taps,0.0);
	std::vector<float> zi_q(rf_taps,0.0);
	std::vector<float> zi_mono0(mono_taps, 0.0);
	std::vector<float> zi_mono1(mono_taps, 0.0);
	std::vector<float> i_temp(block_size/20, 0.0);
	std::vector<float> q_temp(block_size/20, 0.0);
	std::vector<float> mono_Output(block_size/(2*rf_decim*(mono_decim/mono_upscale)), 0.0);
	std::vector<float> fm_demodblock(block_size/20, 0.0);

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	std::cerr << "Block_size is " << block_size << std::endl;

	// Coefficient vector rf_coeff created to store values for the front-end LPF: impulseResponseLPF
	std::vector<float> rf_coeff(rf_taps, 0.0);
	impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff, mono_upscale);
	std::cerr << "Size of rf coeff vector at 100 taps is " << rf_coeff.size() << std::endl;

	// Coefficient vector mono_coeff created to store values for the Mono path LPF: impulseResponseLPF
	std::vector<float> mono_coeffMode0(mono_taps, 0.0);
	impulseResponseLPF(mono_Fs, mono_Fc, mono_taps, mono_coeffMode0, mono_upscale);
	std::cerr << "Size of Mono0 coeff vector at 100 taps is " << mono_coeffMode0.size() << std::endl;
	std::vector<float> mono_coeffMode1(mono_taps*mono_upscale, 0.0);
	impulseResponseLPF(mono_Fs*mono_upscale, mono_Fc, mono_taps*mono_upscale, mono_coeffMode1, mono_upscale);
	std::cerr << "Size of Mono1 coeff vector at 2400 taps is " << mono_coeffMode1.size() << std::endl;

	// Coeffcient vector created to store values for the stereo paths
	std::vector<float> audio_coeffS(mono_taps, 0.0);
	impulseResponseBPF(mono_Fs, 22e3, 54e3, mono_taps, audio_coeffS);

	std::vector<float> audio_coeffP(mono_taps, 0.0);
	impulseResponseBPF(mono_Fs, 18.5e3, 19.5e3, mono_taps, audio_coeffP);

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Stereo variables for Mode 0 and 1
	std::vector<float> zi_stereoMode0(mono_taps, 0.0);
	std::vector<float> zi_stereoMode1(mono_taps, 0.0);
	std::vector<float> zi_stereoP(mono_taps-1, 0.0);
	std::vector<float> zi_stereoS(mono_taps-1, 0.0);
	std::vector<float> stereo_data_block(fm_demodblock.size(), 0.0);
	std::vector<float> LeftandRightStereoBlocks(block_size/(rf_decim*(mono_decim/mono_upscale)), 0.0);
	std::vector<float> audio_stereo(fm_demodblock.size(), 0.0);
	std::vector<float> audio_pilot(fm_demodblock.size(), 0.0);
	std::vector<float> ncoOut(fm_demodblock.size(), 0.0);
	std::vector<float> audio_block_stereo(block_size/(2*rf_decim*(mono_decim/mono_upscale)), 0.0);

	// Output vectors for the fWrite Command
	std::vector<short int> audio_Output_Stereo(block_size/(rf_decim*(mono_decim/mono_upscale)), 0.0);
	std::vector<short int> audio_data_mono(block_size/(rf_decim*(mono_decim/mono_upscale)*2));

	// State saving variables for PLL
	float integrator = 0.0;
	float phaseEst = 0.0;
	float feedbackI = 1.0;
	float feedbackQ = 0.0;
	float trigOffset = 0;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Beginning of Block Processing
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	for(unsigned int block_id = 0; ; block_id++){
		// Initialize vectors to hold normalized I/Q values
		readStdinBlockData(block_size, block_id, i_data, q_data);

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Exit Condition
		if((std::cin.rdstate()) != 0){
			std::cerr << "End of input stream reached" << std::endl;
			break;
		}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Runs between 3.5ms to 6.2ms
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Filtering for I data in blocks
		convolveBlockMode0(i_temp, i_data, rf_coeff, rf_decim, zi_i);
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Runs between 3.5ms to 6.2ms
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Filtering for Q data in blocks
		convolveBlockMode0(q_temp, q_data, rf_coeff, rf_decim, zi_q);
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Runs between 0.056ms to 0.11ms
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Demodulating the I and Q data each block
		fmDemod(i_temp, q_temp, state_I, state_Q, fm_demodblock);
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Runs between 0.056ms to 0.11ms
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		// Stereo Processing before Mono 0 and 1
		convolveBlockMode0(audio_stereo, fm_demodblock, audio_coeffS, 1, zi_stereoS);
		convolveBlockMode0(audio_pilot, fm_demodblock, audio_coeffP, 1, zi_stereoP);

		// Initializing ncoOut for state saving
		if (block_id == 0){
			ncoOut[0] = 1.0; // Predefined value at the beginning
		}else{
			ncoOut[0] = ncoOut[ncoOut.size()-1]; // Grabbing the previous value of the ncoOut as the first (state saving)
		}

		// Calling the PLL
		fmPLL(audio_pilot, 19e3, mono_Fs, 2, 0.0, 0.01, ncoOut, integrator, phaseEst, feedbackI, feedbackQ, trigOffset);

		// Incrementing trigOffset for state saving
		trigOffset += audio_pilot.size();

		// Stereo Mixer
		for (int i = 0; i < audio_stereo.size(); i++){
			stereo_data_block[i] = ncoOut[i]*audio_stereo[i]*2;
		}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		// Mode 0 Processing
		if (mode == 0){

			// Digital Filtering sample rate conversion for stereo (assuming mono 0)
			convolveBlockMode0(audio_block_stereo, stereo_data_block, mono_coeffMode0, mono_decim, zi_stereoMode0);

			// Convolving and Decimation for mono 0
			convolveBlockMode0(mono_Output, fm_demodblock, mono_coeffMode0, mono_decim, zi_mono0);

			// Stereo Combiner
			for (int i = 0; i < audio_block_stereo.size(); i++){
				LeftandRightStereoBlocks[i*2] = (audio_block_stereo[i] + mono_Output[i]) / 2;
				LeftandRightStereoBlocks[(i*2)+1] = (mono_Output[i] - audio_block_stereo[i]) / 2;
			}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		// Mode 1 Processing
		}else if (mode == 1){

			// Resizing the output vectors when in mode 1 for the convolve-upsampler-decimator function
			mono_Output.resize(984, 0.0);
			audio_block_stereo.resize(984, 0.0);

			// Digital Filtering sample rate conversion for stereo (assuming mono 1)
			convolveBlockMode1(audio_block_stereo, stereo_data_block, mono_coeffMode1, mono_upscale, mono_decim, zi_stereoMode1, stereo_prevIndex, stereo_yEnd);

			// Upsampling, Convolving and Decimation for mono 1
			convolveBlockMode1(mono_Output, fm_demodblock, mono_coeffMode1, mono_upscale, mono_decim, zi_mono1, mono_prevIndex, mono_yEnd);


			// Stereo Combiner
			for (int i = 0; i < audio_block_stereo.size(); i++){
				LeftandRightStereoBlocks[i*2] = (audio_block_stereo[i] + mono_Output[i]) / 2;
				LeftandRightStereoBlocks[(i*2)+1] = (mono_Output[i] - audio_block_stereo[i]) / 2;
			}

		}

		// Graphing and gerenting bin files for various vectors for debugging (in certain blocks)
		/*
		if (block_id == 11 || block_id == 10){

			printRealVector(mono_Output);

			std::vector<float> psdfmdemod_est;
			std::vector<float> freq_fmdemod;
			std::vector<float> psdmono_est;
			std::vector<float> freq_mono;
			std::vector<float> psdstereo_est;
			std::vector<float> freq_stereo;
			std::vector<float> psdpilot_est;
			std::vector<float> freq_pilot;
			// std::cerr << "Printing fm_demodblock:" << std::endl;
			// printRealVector(fm_demodblock);
			estimatePSD(fm_demodblock, mono_Fs / 1000, psdfmdemod_est, freq_fmdemod);
			logVector("fmdemod_psd", freq_fmdemod, psdfmdemod_est);

			estimatePSD(audio_stereo, mono_Fs / 1000, psdstereo_est, freq_stereo);
			logVector("fmstereo_psd", freq_stereo, psdstereo_est);

			estimatePSD(audio_pilot, mono_Fs / 1000, psdpilot_est, freq_pilot);
			logVector("fmpilot_psd", freq_pilot, psdpilot_est);

			estimatePSD(mono_Output, 48e3 / 1000, psdmono_est, freq_mono);
			logVector("demodmono_psd", freq_mono, psdmono_est);

			std::vector<float> psdstereoleft_est;
			std::vector<float> psdstereoright_est;
			std::vector<float> freq_stereoleft;
			std::vector<float> freq_stereoright;

			std::vector<float> indexVector;

			if (block_id == 11){
				std::vector<float> ncoOutBeg;
				ncoOutBeg = std::vector<float>(std::begin(ncoOut), std::begin(ncoOut) + 40);
				genIndexVector(indexVector, ncoOutBeg.size());
				logVector("ncoOut11_psd", indexVector, ncoOutBeg);
			}else{
				std::vector<float> ncoOutEnd;
				ncoOutEnd = std::vector<float>(std::begin(ncoOut) + ncoOut.size() - 40, std::begin(ncoOut) + ncoOut.size());
				genIndexVector(indexVector, ncoOutEnd.size());
				logVector("ncoOut10_psd", indexVector, ncoOutEnd);
			}

			std::cerr << "After outputting block 10" << std::endl;
		}
		*/

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// abou 0.004ms to 0.1ms
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Outputting the stereo audio

	// Value casting for audio output (for stereo)
	for (unsigned int k=0; k<LeftandRightStereoBlocks.size(); k++){
			if (std::isnan(LeftandRightStereoBlocks[k])){
				audio_Output_Stereo[k]=0;
				std::cerr << "NANs running " << std::endl;
			}else{
				audio_Output_Stereo[k]=(static_cast<short int>(LeftandRightStereoBlocks[k]*16384));
			}
	}

	// Audio Output (for stereo)
	fwrite(&audio_Output_Stereo[0], sizeof(short int), audio_Output_Stereo.size()-1, stdout);

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Outputting the mono audio instead of the stereo audio

	// // Value casting for audio output (for mono)
	// for (unsigned int k=0; k<mono_Output.size(); k++){
	// 		if (std::isnan(mono_Output[k])){
	// 			audio_data_mono[k]=0;
	// 			std::cerr << "NANs running " << std::endl;
	// 		}else{
	// 			audio_data_mono[k]=(static_cast<short int>(mono_Output[k]*16384));
	// 		}
	// }
	//
	// // Audio Output (for mono)
	// fwrite(&audio_data_mono[0], sizeof(short int), audio_data_mono.size()-1, stdout);

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	}

	return 0;
}
