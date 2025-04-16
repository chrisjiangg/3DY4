/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>

#define PI 3.14159265358979323846
typedef float real;

// Function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(real Fs, real Fc, unsigned short int num_taps, std::vector<real> &h)
{
    // Bring your own functionality
    h.clear();
    h.resize(num_taps, 0.0);

    real cutoff = Fc/(Fs/2);
    for(int i=0; i<num_taps; i++){
        if(i == (num_taps-1)/2){
            h[i] = cutoff;
        }
        else{
            h[i] = cutoff * (sin(PI*cutoff*(i-(num_taps-1)/2)))/(PI*cutoff*(i-(num_taps-1)/2));
        }
        h[i] = h[i]*pow(sin((i*PI)/num_taps), 2);
    }
}

// Function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
void convolveFIR(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h)
{
    // Bring your own functionality
    y.clear();
    y.resize(x.size() + h.size() - 1, 0.0);

    for(int i=0; i < y.size(); i++){
        for(int j=0; j < h.size(); j++){
            if(((i-j)>=0) && ((i-j)< x.size())){
                y[i] += h[j] * x[i-j]; 
            }
        }
    }
}