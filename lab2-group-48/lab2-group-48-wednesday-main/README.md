## DSP Primitives in C++

The purpose of this 3DY4 lab is to provide a link between the modelling of (DSP) primitives used in Software Deﬁned Radios (SDRs) in an interpreted and scripting language like Python and their implementation in C++. The details are in the lab [document](doc/3dy4-lab2.pdf).

The start-up C++ code is in the [src](src/) sub-folder. To compile C++ source files from the terminal, use:

`g++ my_src.cpp -o my_sw`

If you wish to use a particular version of C++, e.g. C++17, you can employ the `-std` command line argument:

`g++ -std=c++17 my_src.cpp -o my_sw`

A few good pointers to get you started with GDB are available [here](doc/gdb-getting-started.md).

Both the in-lab experiments and the take-home exercises are due 16 hours before your next lab is scheduled.
