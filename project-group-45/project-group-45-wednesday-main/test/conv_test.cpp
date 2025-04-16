/*
   Comp Eng 3DY4 (Computer Systems Integration Project)

   Department of Electrical and Computer Engineering
   McMaster University
   Ontario, Canada
*/

// This file shows how to write convolution unit tests, using Google C++ test framework.
// (it is based on https://github.com/google/googletest/blob/main/docs/index.md)

#include <limits.h>
#include "dy4.h"
#include "iofunc.h"
#include "logfunc.h"
#include "filter.h"
#include "fourier.h"
#include "functions.h"
#include "gtest/gtest.h"
#include "algorithm"

namespace {

	class Convolution_Fixture: public ::testing::Test {

		public:

			const int N = 1024;	// signal size
			const int M = 101;	// kernel size
			const int lower_bound = -1;
			const int upper_bound = 1;
			const real EPSILON = 1e-4;

			std::vector<real> x, h, y_reference, y_test, state;

			Convolution_Fixture() {
				x.resize(N);
				h.resize(M);
				y_reference.resize(N + M - 1);
				y_test.resize(N + M - 1);
				state.resize(M-1, 0.0);
			}

			void SetUp() {
				generate_random_values(x, lower_bound, upper_bound);
				generate_random_values(h, lower_bound, upper_bound);
				convolveFIR_reference(y_reference, x, h);
			}

			void TearDown() {
			}

			~Convolution_Fixture() {
			}
	};

	TEST_F(Convolution_Fixture, convolveFIR_inefficient_NEAR) {

		convolveFIR_inefficient(y_test, x, h);

		ASSERT_EQ(y_reference.size(), y_test.size()) << "Output vector sizes for convolveFIR_reference and convolveFIR_inefficient are unequal";

		for (int i = 0; i < (int)y_reference.size(); ++i) {
			EXPECT_NEAR(y_reference[i], y_test[i], EPSILON) << "Original/convolveFIR_inefficient vectors differ at index " << i;
		}
	}
	TEST_F(Convolution_Fixture, my_own_conv_state) {

		my_own_conv_state(state, x, h, y_test, 1 ,1);
		
		//ASSERT_EQ(y_reference.size(), y_test.size()) << "Output vector sizes for convolveFIR_reference and convolveFIR_inefficient are unequal";

		for (int i = 0; i < (int)y_test.size(); ++i) {
			EXPECT_NEAR(y_reference[i], y_test[i], EPSILON) << "Original/my_own_conv_state vectors differ at index " << i;
		}
	}

	TEST_F(Convolution_Fixture, conv_down) {

		state.clear();
		conv_down(state, x, h, y_test,1);
		// for (int i=0;i<y_test.size();i++){
		// 	if(y_reference[i]+EPSILON >y_test[i] && y_reference[i]-EPSILON <y_test[i]){
		// 		std::cerr<< "Index: "<< i <<" ref : "<<y_reference[i]<<" test : "<<y_test[i]<<std::endl;
		// 	}else
		// 		std::cerr<< "Index: "<< i <<" ref : "<<y_reference[i]<<" test : "<<y_test[i]<<std::endl;
		// }
		//ASSERT_EQ(y_reference.size(), y_test.size()) << "Output vector sizes for convolveFIR_reference and convolveFIR_inefficient are unequal";

		for (int i = 0; i < (int)y_test.size(); ++i) {
			EXPECT_NEAR(y_reference[i], y_test[i], EPSILON) << "Original/conv_down vectors differ at index " << i;
		}
	}
	TEST_F(Convolution_Fixture, upsample) {

		std::vector<real> x_test;
		x_test.push_back(1);
		x_test.push_back(2);
		x_test.push_back(3);

		x.clear();
		x.push_back(1);
		x.push_back(0);
		x.push_back(2);
		x.push_back(0);
		x.push_back(3);
		x.push_back(0);
		upsample(x_test, 2);
		
		// for (int i=0;i<x_test.size();i++){
		// 	std::cerr<< "Index: "<< i <<" ref : "<<x[i]<<" test : "<<x_test[i]<<std::endl;
		// }
		
		ASSERT_EQ(x.size(), x_test.size()) << "Output vector sizes for x and x_upsample are unequal";

		for (int i = 0; i < (int)x_test.size(); ++i) {
			EXPECT_NEAR(x[i], x_test[i], EPSILON) << "Original/x_upsample vectors differ at index " << i;
		}
	}
	TEST_F(Convolution_Fixture, downsample) {

		std::vector<real> x_test;
		x_test.push_back(1);
		x_test.push_back(2);
		x_test.push_back(3);

		x.clear();
		x.push_back(1);
		x.push_back(0);
		x.push_back(2);
		x.push_back(0);
		x.push_back(3);
		x.push_back(0);
		downsample(x, 2);
		
		// for (int i=0;i<x_test.size();i++){
		// 	std::cerr<< "Index: "<< i <<" ref : "<<x[i]<<" test : "<<x_test[i]<<std::endl;
		// }
		
		ASSERT_EQ(x.size(), x_test.size()) << "Output vector sizes for x and x_downsample are unequal";

		for (int i = 0; i < (int)x_test.size(); ++i) {
			EXPECT_NEAR(x[i], x_test[i], EPSILON) << "Original/x_downsample vectors differ at index " << i;
		}
	}

	TEST_F(Convolution_Fixture, downsample_convolution) {

		x.resize(N);
		h.resize(M);
		y_reference.resize(N + M - 1);
		y_test.resize(N + M - 1);
		state.resize(M-1, 0.0);
		generate_random_values(x, lower_bound, upper_bound);
		generate_random_values(h, lower_bound, upper_bound);


		my_own_conv_state(state, x, h, y_test, 2 ,1);

		//upsample(x,2);
		convolveFIR_reference(y_reference, x, h);
		downsample(y_reference,2);
		
		
		// for (int i=0;i<y_test.size();i++){
		// 	if(y_reference[i]+EPSILON >y_test[i] && y_reference[i]-EPSILON <y_test[i]){
		// 		//std::cerr<< "Index: "<< i <<" ref : "<<y_reference[i]<<" test : "<<y_test[i]<<std::endl;
		// 	}else
		// 		std::cerr<< "Index: "<< i <<" ref : "<<y_reference[i]<<" test : "<<y_test[i]<<std::endl;
		// }
		
		//ASSERT_EQ(y_reference.size(), y_test.size()) << "Output vector sizes for x and x_downsample are unequal";

		for (int i = 0; i < (int)y_test.size(); ++i) {
			EXPECT_NEAR(y_reference[i], y_test[i], EPSILON) << "Original/x_downsample vectors differ at index " << i;
		}
	}
	TEST_F(Convolution_Fixture, upsample_convolution) {
		//int n = 500;

		x.resize(N);
		h.resize(M);
		state.clear();
		state.resize(M-1, 0.0);
		y_reference.clear();
		y_reference.resize(N + M - 1);
		y_test.clear();

		generate_random_values(x, lower_bound, upper_bound);
		generate_random_values(h, lower_bound, upper_bound);

		my_own_conv_state(state, x, h, y_test, 1,10);

		upsample(x,10);
		convolveFIR_reference(y_reference, x, h);

		//down_sample(y_reference,2);
		
		
		for (long unsigned int i=0;i<y_test.size();i++){
			if(y_reference[i]+EPSILON >y_test[i] && y_reference[i]-EPSILON <y_test[i]){
				//std::cerr<< "Index: "<< i <<" ref : "<<y_reference[i]<<" test : "<<y_test[i]<<std::endl;
			}else
				std::cerr<< "Index: "<< i <<" ref : "<<y_reference[i]<<" test : "<<y_test[i]<<std::endl;
		}
		
		//ASSERT_EQ(y_reference.size(), y_test.size()) << "Output vector sizes for x and x_downsample are unequal";

		for (int i = 0; i < (int)y_test.size(); ++i) {
			EXPECT_NEAR(y_reference[i], y_test[i], EPSILON) << "Original/x_downsample vectors differ at index " << i;
		}
	}
	TEST_F(Convolution_Fixture, PLL) {

		// Initialize test parameters
        real freq = 19000.0; // Example: RDS carrier frequency (19 kHz)
        real fs = 114000.0;   // Example: Sampling frequency (114 kHz)
        real ncoscale = 2.0;
        real normBandwidth = 0.01;
        
        // Initial state variables
        real integrator = 0.0;
        real phaseEst = 0.0;
        real feedbackI = 1.0;
        real feedbackQ = 0.0;
        real ncoOut1 = 1.0;
        real trigOffset = 0.0;

		// Generate a test input signal (simulated FM signal with noise)
		int numSamples = 500;
		std::vector<real> pllIn(numSamples);
		std::vector<real> ncoOut;
		
		for (int i = 0; i < numSamples; ++i) {
			pllIn[i] = 0.5*std::sin(2 * M_PI * freq * i / fs) + 0.1 * ((rand() % 100) / 100.0); // Add small noise
		}
		real max_val = *std::max_element(pllIn.begin(), pllIn.end());
		if (max_val > 0) {
			for (auto &p : pllIn) {
				p /= max_val;  // Normalize to 1.0
			}
		}

		fmPll(pllIn, freq, fs, ncoscale, normBandwidth, ncoOut, integrator, phaseEst, feedbackI, feedbackQ, ncoOut1, trigOffset);
		
		// for (int i=0;i<x_test.size();i++){
		// 	std::cerr<< "Index: "<< i <<" ref : "<<x[i]<<" test : "<<x_test[i]<<std::endl;
		// }
		
		ASSERT_EQ(ncoOut.size(), pllIn.size() + 1) << "NCO output size mismatch";

		// Check that the phase estimates and feedback values remain bounded
		EXPECT_LT(std::abs(phaseEst), 2 * M_PI) << "Phase estimation out of expected range";
		EXPECT_LT(std::abs(feedbackI), 1.1) << "FeedbackI out of range";
		EXPECT_LT(std::abs(feedbackQ), 1.1) << "FeedbackQ out of range";

		//plot output
		std::vector<real> vector_index;
		genIndexVector(vector_index, pllIn.size());
		logVector("PLL_Input", vector_index, pllIn);
		genIndexVector(vector_index, ncoOut.size());
		logVector("PLL_Output", vector_index, ncoOut);
		
	}

} // end of namespace
