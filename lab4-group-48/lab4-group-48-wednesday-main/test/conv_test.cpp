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
#include "filter.h"
#include "gtest/gtest.h"

namespace {

	class Convolution_Fixture: public ::testing::Test {

		public:

			const int N = 1024;	// signal size
			const int M = 101;	// kernel size
			const int lower_bound = -1;
			const int upper_bound = 1;
			const real EPSILON = 1e-4;

			std::vector<real> x, h, y_reference, y_test;

			Convolution_Fixture() {
				x.resize(N);
				h.resize(M);
				y_reference.resize(N + M - 1);
				y_test.resize(N + M - 1);
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
//copy this function for all additional algorithms
	TEST_F(Convolution_Fixture, convolveFIR_inefficient_NEAR) { //change the second input which is the name of the test
		//change the line below to the function we want to test depending on the function may need to change/add inputs
		convolveFIR_inefficient(y_test, x, h);

		ASSERT_EQ(y_reference.size(), y_test.size()) << "Output vector sizes for convolveFIR_reference and convolveFIR_inefficient are unequal";

		for (int i = 0; i < (int)y_reference.size(); ++i) {
			EXPECT_NEAR(y_reference[i], y_test[i], EPSILON) << "Original/convolveFIR_inefficient vectors differ at index " << i;
		}
	}
	//convolveFIR_signal test
	TEST_F(Convolution_Fixture, convolveFIR_signal_NEAR) { //change the second input which is the name of the test
		//change the line below to the function we want to test depending on the function may need to change/add inputs
		convolveFIR_signal(y_test, x, h);

		ASSERT_EQ(y_reference.size(), y_test.size()) << "Output vector sizes for convolveFIR_reference and convolveFIR_signal are unequal";

		for (int i = 0; i < (int)y_reference.size(); ++i) {
			EXPECT_NEAR(y_reference[i], y_test[i], EPSILON) << "Original/convolveFIR_signal vectors differ at index " << i;
		}
	}
	//convolveFIR_kernel test
	TEST_F(Convolution_Fixture, convolveFIR_kernel_NEAR) { //change the second input which is the name of the test
		//change the line below to the function we want to test depending on the function may need to change/add inputs
		convolveFIR_kernel(y_test, x, h);

		ASSERT_EQ(y_reference.size(), y_test.size()) << "Output vector sizes for convolveFIR_reference and convolveFIR_kernel are unequal";

		for (int i = 0; i < (int)y_reference.size(); ++i) {
			EXPECT_NEAR(y_reference[i], y_test[i], EPSILON) << "Original/convolveFIR_kernel vectors differ at index " << i;
		}
	}
	//convolveFIR_precomputed test
	TEST_F(Convolution_Fixture, convolveFIR_precomputed_NEAR) { //change the second input which is the name of the test
		//change the line below to the function we want to test depending on the function may need to change/add inputs
		convolveFIR_precomputed(y_test, x, h);

		ASSERT_EQ(y_reference.size(), y_test.size()) << "Output vector sizes for convolveFIR_reference and convolveFIR_precomputed are unequal";

		for (int i = 0; i < (int)y_reference.size(); ++i) {
			EXPECT_NEAR(y_reference[i], y_test[i], EPSILON) << "Original/convolveFIR_precomputed vectors differ at index " << i;
		}
	}
	//convolveFIR_unrolled test
	TEST_F(Convolution_Fixture, convolveFIR_unrolled_NEAR) { //change the second input which is the name of the test
		//change the line below to the function we want to test depending on the function may need to change/add inputs
		convolveFIR_unrolled(y_test, x, h);

		ASSERT_EQ(y_reference.size(), y_test.size()) << "Output vector sizes for convolveFIR_reference and convolveFIR_unrolled are unequal";

		for (int i = 0; i < (int)y_reference.size(); ++i) {
			EXPECT_NEAR(y_reference[i], y_test[i], EPSILON) << "Original/convolveFIR_unrolled vectors differ at index " << i;
		}
	}
	//convolveFIR_kernel_updated test
	TEST_F(Convolution_Fixture, convolveFIR_kernel_updated_NEAR) { //change the second input which is the name of the test
		//change the line below to the function we want to test depending on the function may need to change/add inputs
		convolveFIR_kernel_updated(y_test, x, h);

		ASSERT_EQ(y_reference.size(), y_test.size()) << "Output vector sizes for convolveFIR_reference and convolveFIR_kernel_updated are unequal";

		for (int i = 0; i < (int)y_reference.size(); ++i) {
			EXPECT_NEAR(y_reference[i], y_test[i], EPSILON) << "Original/convolveFIR_kernel_updated vectors differ at index " << i;
		}
	}
	//convolveFIR_precomputed_updated test
	TEST_F(Convolution_Fixture, convolveFIR_precomputed_updated_NEAR) { //change the second input which is the name of the test
		//change the line below to the function we want to test depending on the function may need to change/add inputs
		convolveFIR_precomputed_updated(y_test, x, h);

		ASSERT_EQ(y_reference.size(), y_test.size()) << "Output vector sizes for convolveFIR_reference and convolveFIR_precomputed_updated are unequal";

		for (int i = 0; i < (int)y_reference.size(); ++i) {
			EXPECT_NEAR(y_reference[i], y_test[i], EPSILON) << "Original/convolveFIR_precomputed_updated vectors differ at index " << i;
		}
	}
	//convolveFIR_unrolled_updated test
	TEST_F(Convolution_Fixture, cconvolveFIR_unrolled_updated_NEAR) { //change the second input which is the name of the test
		//change the line below to the function we want to test depending on the function may need to change/add inputs
		convolveFIR_precomputed_updated(y_test, x, h);

		ASSERT_EQ(y_reference.size(), y_test.size()) << "Output vector sizes for convolveFIR_reference and convolveFIR_unrolled_updated are unequal";

		for (int i = 0; i < (int)y_reference.size(); ++i) {
			EXPECT_NEAR(y_reference[i], y_test[i], EPSILON) << "Original/convolveFIR_unrolled_updated vectors differ at index " << i;
		}
	}
} // end of namespace
