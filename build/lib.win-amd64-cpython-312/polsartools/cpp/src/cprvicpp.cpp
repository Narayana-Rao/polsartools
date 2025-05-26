#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include <stdexcept>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <limits> // For handling NaN and infinity
namespace py = pybind11;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

py::array_t<double> process_chunk_cprvicpp(
    const std::vector<py::array_t<double>>& chunks,
    int window_size,
    // const std::vector<std::string>& input_filepaths,
    double chi_in,
    double psi_in
) {
    if (chunks.size() < 4) {
        throw std::runtime_error("At least 4 chunks required");
    }

    py::buffer_info buf_c11 = chunks[0].request();
    py::buffer_info buf_c12r = chunks[1].request();
    py::buffer_info buf_c12i = chunks[2].request();
    py::buffer_info buf_c22 = chunks[3].request();

    int nrows = buf_c11.shape[1];
    int ncols = buf_c11.shape[0];

    const double* c11_T1 = static_cast<double*>(buf_c11.ptr);
    const double* c12r_T1 = static_cast<double*>(buf_c12r.ptr);
    const double* c12i_T1 = static_cast<double*>(buf_c12i.ptr);
    const double* c22_T1 = static_cast<double*>(buf_c22.ptr);

    auto fp22_out = py::array_t<double>({ncols, nrows});
    auto lambda_out = py::array_t<double>({ncols, nrows});
    double* fp22 = static_cast<double*>(fp22_out.request().ptr);
    double* l_lambda = static_cast<double*>(lambda_out.request().ptr);

    int inci = window_size / 2;
    int incj = window_size / 2;
    int starti = inci;
    int startj = incj;
    int stopi = nrows - inci - 1;
    int stopj = ncols - incj - 1;

    for (int ii = startj; ii <= stopj; ++ii) {
        for (int jj = starti; jj <= stopi; ++jj) {
            int count = 0;
            std::complex<double> C11c(0.0), C12c(0.0), C21c(0.0), C22c(0.0);

            for (int wi = -inci; wi <= inci; ++wi) {
                for (int wj = -incj; wj <= incj; ++wj) {
                    int x = ii + wi;
                    int y = jj + wj;
                    if (x >= 0 && x < ncols && y >= 0 && y < nrows) {
                        int idx = x * nrows + y;
                        C11c += c11_T1[idx];
                        C12c += std::complex<double>(c12r_T1[idx], c12i_T1[idx]);
                        C21c += std::conj(std::complex<double>(c12r_T1[idx], c12i_T1[idx]));
                        C22c += c22_T1[idx];
                        count++;
                    }
                }
            }

            if (count > 0) {
                C11c /= count;
                C12c /= count;
                C21c /= count;
                C22c /= count;
            }

            if (!std::isfinite(C11c.real()) || !std::isfinite(C12c.real()) || !std::isfinite(C21c.real()) || !std::isfinite(C22c.real())) {
                C11c = C12c = C21c = C22c = {0.0, 0.0};
            }

            std::complex<double> s0 = C11c + C22c;
            std::complex<double> s1 = C11c - C22c;
            std::complex<double> s2 = C12c + C21c;
            std::complex<double> s3 = chi_in >= 0 ? std::complex<double>(0, 1) * (C12c - C21c)
                                                  : -std::complex<double>(0, 1) * (C12c - C21c);

            std::complex<double> K[4][4] = {
                {s0, 0.0, s2, 0.0},
                {0.0, 0.0, 0.0, s1},
                {s2, 0.0, 0.0, 0.0},
                {0.0, s1, 0.0, s3}
            };

            std::complex<double> SC = (s0 - s3) * 0.5;
            std::complex<double> OC = (s0 + s3) * 0.5;

            double min_val = std::min(SC.real(), OC.real());
            double max_val = std::max(SC.real(), OC.real());

            std::complex<double> K_depol[4][4] = {};
            K_depol[0][0] = 1.0;

            std::complex<double> num = 0.0, den1 = 0.0, den2 = 0.0;

            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b) {
                    num += std::conj(K[b][a]) * K_depol[a][b];
                    den1 += std::conj(K[b][a]) * K[b][a];
                    den2 += std::conj(K_depol[b][a]) * K_depol[b][a];
                }

            double depol_ang = std::real(2 * std::acos(num.real() / std::sqrt(std::abs(den1.real() * den2.real()))) * 180.0 / M_PI);
            double GD_t1_depol = depol_ang / 180.0;
            l_lambda[ii * nrows + jj] = (3.0 / 2.0) * GD_t1_depol;

            fp22[ii * nrows + jj] = min_val / max_val;
        }
    }

    auto vi_c = py::array_t<double>({ncols, nrows});
    double* vi_c_ptr = static_cast<double*>(vi_c.request().ptr);

    for (int i = 0; i < ncols * nrows; ++i) {
        vi_c_ptr[i] = (1.0 - l_lambda[i]) * std::pow(fp22[i], 2 * l_lambda[i]);
    }

    return vi_c;
}

// Pybind11 module definition
PYBIND11_MODULE(cprvicpp, m) {
    m.doc() = "CpRVI Processing Module";
    m.def("process_chunk_cprvicpp", &process_chunk_cprvicpp, 
          "Process a C2 chunk to generate CPRVI"
        );
}