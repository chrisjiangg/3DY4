/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_IOFUNC_H
#define DY4_IOFUNC_H

// Add headers as needed
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>

// Declaration of function prototypes
void printRealVector(const std::vector<real> &);

void printComplexVector(const std::vector<std::complex<real>> &);

void readBinData(const std::string, std::vector<real> &);

void writeBinData(const std::string, const std::vector<real> &);

#endif // DY4_IOFUNC_H
