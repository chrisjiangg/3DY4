/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// Add headers as needed
#include <iostream>
#include <vector>

// Declaration of function prototypes
void impulseResponseLPF(real, real, unsigned short int, std::vector<real> &);
void convolveFIR(std::vector<real> &, const std::vector<real> &, const std::vector<real> &);

#endif // DY4_FILTER_H
