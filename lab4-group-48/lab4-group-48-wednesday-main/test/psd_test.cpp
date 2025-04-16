#include "dy4.h"
#include "gtest/gtest.h"
#include "fourier.h"
#include <vector>
#include <cmath>
#include "iofunc.h"

namespace {

    class psd_Fixture: public ::testing::Test {
    public:
        const int N = 1024; // Number of samples
        const int lower_bound = -1; // Lower bound for random values
        const int upper_bound = 1; // Upper bound for random values
        const real EPSILON = 1e-4; // Tolerance for floating-point comparison
        const int Fs = 1000; // Sampling frequency

        std::vector<real> x; // Input samples
        std::vector<real> freq_ref, freq_mat; // Frequency bins
        std::vector<real> psd_est_ref; // PSD from estimatePSD
        std::vector<real> psd_mat; // PSD from matrixPSD

        psd_Fixture() {
            x.resize(N);
        }

        void SetUp() {
            // Generate random input samples
            generate_random_values(x, lower_bound, upper_bound);
            estimatePSD(x, Fs, freq_ref, psd_est_ref);
            matrixPSD(x, Fs, freq_mat, psd_mat);
        }

        void TearDown() {}
    };

    TEST_F(psd_Fixture, PSD_matrixPSD_NEAR) {
        // Ensure the output sizes match
        ASSERT_EQ(psd_est_ref.size(), psd_mat.size())
            << "Output vector sizes for estimatePSD and matrixPSD are unequal";

        // Compare the PSD values
        for (int i = 0; i < static_cast<int>(psd_est_ref.size()); i++) {
            EXPECT_NEAR(psd_est_ref[i], psd_mat[i], EPSILON)
                << "Original/matrixPSD vectors differ at index " << i;
        }
    }
}
