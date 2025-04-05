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

using Matrix = std::vector<std::vector<double>>;

// Function to calculate mean ignoring NaNs and infinity
double nanmean(const std::vector<std::vector<double>>& mat, int start_row, int end_row, int start_col, int end_col) {
    double sum = 0.0;
    int count = 0;

    for (int i = start_row; i <= end_row; ++i) {
        for (int j = start_col; j <= end_col; ++j) {
            if (!std::isnan(mat[i][j]) && !std::isinf(mat[i][j])) {
                sum += mat[i][j];
                count++;
            }
        }
    }
    return count > 0 ? sum / count : std::numeric_limits<double>::quiet_NaN();
}

// Function to calculate trace of a square matrix
double trace(const Matrix& mat) {
    double sum = 0.0;
    for (size_t i = 0; i < mat.size(); ++i) {
        sum += mat[i][i];
    }
    return sum;
}

// Matrix multiplication function
Matrix multiply(const Matrix& mat1, const Matrix& mat2) {
    size_t nrows = mat1.size();
    size_t ncols = mat2[0].size();
    size_t inner_dim = mat2.size();
    Matrix result(nrows, std::vector<double>(ncols, 0.0));

    for (size_t i = 0; i < nrows; ++i) {
        for (size_t j = 0; j < ncols; ++j) {
            for (size_t k = 0; k < inner_dim; ++k) {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
    return result;
}

// Transpose function for a square matrix
Matrix transpose(const Matrix& mat) {
    size_t nrows = mat.size();
    size_t ncols = mat[0].size();
    Matrix transposed(ncols, std::vector<double>(nrows, 0.0));

    for (size_t i = 0; i < nrows; ++i) {
        for (size_t j = 0; j < ncols; ++j) {
            transposed[j][i] = mat[i][j];
        }
    }
    return transposed;
}

// Function to process chunks with proper GIL handling
py::array_t<double> process_chunk_cprvicpp(
    const std::vector<py::array_t<double>>& chunks,
    int window_size,
    const std::vector<std::string>& input_filepaths,
    double chi_in,
    double psi_in) {

    if (chunks.size() < 4) {
        throw std::runtime_error("Insufficient chunks provided.");
    }

    auto buf0 = chunks[0].request();
    int nrows = buf0.shape[0];
    int ncols = buf0.shape[1];

    std::vector<const double*> chunk_ptrs(chunks.size());
    for (size_t i = 0; i < chunks.size(); ++i) {
        auto buf = chunks[i].request();
        if (buf.ndim != 2 || buf.shape[0] != nrows || buf.shape[1] != ncols)
            throw std::runtime_error("All chunks must be 2D arrays of the same shape.");
        chunk_ptrs[i] = static_cast<const double*>(buf.ptr);
    }

    auto fp22 = std::vector<std::vector<double>>(nrows, std::vector<double>(ncols, 0.0));
    auto l_lambda = fp22;

    int inci = window_size / 2;
    int incj = window_size / 2;
    int starti = inci;
    int stopi = nrows - inci - 1;
    int startj = incj;
    int stopj = ncols - incj - 1;

    {
        py::gil_scoped_release release;

        for (int ii = startj; ii <= stopj; ++ii) {
            for (int jj = starti; jj <= stopi; ++jj) {
                try {
                    auto nanmean_from_ptr = [&](const double* data) {
                        double sum = 0.0;
                        int count = 0;
                        for (int i = ii - inci; i <= ii + inci; ++i) {
                            for (int j = jj - incj; j <= jj + incj; ++j) {
                                double val = data[i * ncols + j];
                                if (!std::isnan(val) && !std::isinf(val)) {
                                    sum += val;
                                    ++count;
                                }
                            }
                        }
                        return count > 0 ? sum / count : std::numeric_limits<double>::quiet_NaN();
                    };

                    double C11c = nanmean_from_ptr(chunk_ptrs[0]);
                    std::complex<double> C12c(nanmean_from_ptr(chunk_ptrs[1]), nanmean_from_ptr(chunk_ptrs[2]));
                    std::complex<double> C21c = std::conj(C12c);
                    double C22c = nanmean_from_ptr(chunk_ptrs[3]);

                    if (std::isnan(C11c) || std::isnan(C22c) || std::isnan(C12c.real()) || std::isnan(C12c.imag())) {
                        throw std::runtime_error("NaN detected in input values.");
                    }

                    // Polarimetric decomposition calculations (unchanged)
                    double s0 = C11c + C22c;
                    double s1 = C11c - C22c;
                    double s2 = std::real(C12c + C21c);
                    double s3 = (chi_in >= 0) ? -std::imag(C12c - C21c) : std::imag(C12c - C21c);

                    double SC = (s0 - s3) / 2.0;
                    double OC = (s0 + s3) / 2.0;
                    double min_sc_oc = std::min(SC, OC);
                    double max_sc_oc = std::max(SC, OC);

                    // Build matrices and compute trace-based depolarization
                    Matrix K_T = {
                        {0.5 * s0, 0, 0.5 * s2, 0},
                        {0, 0, 0, 0.5 * s1},
                        {0.5 * s2, 0, 0, 0},
                        {0, 0.5 * s1, 0, 0.5 * s3}
                    };
                    Matrix K_depol = {
                        {1.0, 0, 0, 0},
                        {0, 0, 0, 0},
                        {0, 0, 0, 0},
                        {0, 0, 0, 0}
                    };

                    double num = trace(multiply(transpose(K_T), K_depol));
                    double den1 = std::sqrt(std::abs(trace(multiply(transpose(K_T), K_T))));
                    double den2 = std::sqrt(std::abs(trace(multiply(transpose(K_depol), K_depol))));
                    double den = den1 * den2;

                    double GD_t1_depol = (den == 0 || std::isnan(num) || std::isnan(den)) ? 0.0 :
                        2.0 * std::acos(std::clamp(num / den, -1.0, 1.0)) / M_PI;

                    l_lambda[ii][jj] = 1.5 * GD_t1_depol;
                    fp22[ii][jj] = (max_sc_oc > 0.0 && !std::isnan(max_sc_oc)) ? (min_sc_oc / max_sc_oc) : 0.0;

                } catch (const std::exception& e) {
                    l_lambda[ii][jj] = 0.0;
                    fp22[ii][jj] = 0.0;
                }
            }
        }
    }

    // Final CPRVI matrix computation
    py::array_t<double> result({nrows, ncols});
    auto r = result.mutable_unchecked<2>();
    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < ncols; ++j) {
            r(i, j) = (1.0 - l_lambda[i][j]) * std::pow(fp22[i][j], 2.0 * l_lambda[i][j]);
        }
    }

    return result;
}
// Pybind11 module definition

PYBIND11_MODULE(cprvicpp, m) {
    m.doc() = "CPRVI Processing Module";
    m.def("process_chunk_cprvicpp", &process_chunk_cprvicpp, 
          "Process a chunk using CPRVI algorithm",
          py::arg("chunks"), 
          py::arg("window_size"), 
          py::arg("input_filepaths"), 
          py::arg("chi_in"), 
          py::arg("psi_in")
        //   py::call_guard<py::gil_scoped_release>() // GIL is released here
        );
}