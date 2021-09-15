/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <vector>

// declaration of a function prototypes
void impulseResponseLPF(float, float, unsigned short int, std::vector<float> &, int);
void convolveBlockMode0(std::vector<float> &, const std::vector<float> &, const std::vector<float> &, int, std::vector<float> &);
void convolveBlockMode1(std::vector<float> &, const std::vector<float> &, const std::vector<float> &, int, int, std::vector<float> &, int &, float &);
void impulseResponseBPF(float , float , float , unsigned short int, std::vector<float> &);
void fmPLL(std::vector<float> , float , float , float , float , float , std::vector<float> &, float &, float &, float &, float &, float &);

#endif // DY4_FILTER_H
