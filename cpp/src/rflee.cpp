#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include<pybind11/complex.h>
#include<pybind11/functional.h>
#include<pybind11/chrono.h>
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <iostream>
#include <map>        // Add this line to include the map container
#include <stdexcept>  // For std::invalid_argument




namespace py = pybind11;
using ComplexMatrix = std::vector<std::vector<std::complex<double>>>;
using Matrix = std::vector<std::vector<float>>;

Matrix compute_span(const std::string& PolTypeOut, const std::vector<ComplexMatrix>& M_in) {
    // Get dimensions of M_in (assuming all matrices in M_in are of the same size)
    int Nlig_padded = M_in[0].size();   // Number of rows
    int Ncol_padded = M_in[0][0].size(); // Number of columns
    
    // Create a 2D vector for span (matrix) initialized to 0
    Matrix span(Nlig_padded, std::vector<float>(Ncol_padded, 0.0f));

    // Fill span matrix based on PolTypeOut value
    if (PolTypeOut == "C2" || PolTypeOut == "C2pp1" || PolTypeOut == "C2pp2" || PolTypeOut == "C2pp3" ||
        PolTypeOut == "T2" || PolTypeOut == "T2pp1" || PolTypeOut == "T2pp2" || PolTypeOut == "T2pp3") {
        // Compute span for "C2" or "T2" - Add M_in[0] and M_in[3]
        for (int i = 0; i < Nlig_padded; ++i) {
            for (int j = 0; j < Ncol_padded; ++j) {
                span[i][j] = std::real(M_in[0][i][j]) + std::real(M_in[3][i][j]);
            }
        }
    } else if (PolTypeOut == "C3" || PolTypeOut == "T3") {
        // Compute span for "C3" or "T3" - Add M_in[0], M_in[5], and M_in[8]
        for (int i = 0; i < Nlig_padded; ++i) {
            for (int j = 0; j < Ncol_padded; ++j) {
                span[i][j] = std::real(M_in[0][i][j]) + std::real(M_in[5][i][j]) + std::real(M_in[8][i][j]);
            }
        }
    } else if (PolTypeOut == "C4" || PolTypeOut == "T4") {
        // Compute span for "C4" or "T4" - Add M_in[0], M_in[7], M_in[12], and M_in[15]
        for (int i = 0; i < Nlig_padded; ++i) {
            for (int j = 0; j < Ncol_padded; ++j) {
                span[i][j] = std::real(M_in[0][i][j]) + std::real(M_in[7][i][j]) + std::real(M_in[12][i][j]) + std::real(M_in[15][i][j]);
            }
        }
    } else {
        throw std::invalid_argument("Unsupported PolTypeOut");
    }

    return span;
}
// Get the window parameters from the window_size
std::pair<int, int> get_window_params(int window_size) {
    // Define window_params map as in the Python code
    std::map<int, std::pair<int, int>> window_params = {
        {3, {1, 1}},
        {5, {3, 1}},
        {7, {3, 2}},
        {9, {5, 2}},
        {11, {5, 3}},
        {13, {5, 4}},
        {15, {7, 4}},
        {17, {7, 5}},
        {19, {7, 6}},
        {21, {9, 6}},
        {23, {9, 7}},
        {25, {9, 8}},
        {27, {11, 8}},
        {29, {11, 9}},
        {31, {11, 10}}
    };

    // Check if the window_size exists in the map, if not throw an exception
    if (window_params.find(window_size) != window_params.end()) {
        return window_params[window_size];
    } else {
        throw std::invalid_argument("The window width window_size must be set to 3 to 31");
    }
}


std::vector<std::vector<std::vector<std::complex<double>>>> compose_C3_chunks(
    const std::vector<std::vector<std::vector<std::complex<double>>>>& chunks) {
    int nrows = chunks[0].size();
    int ncols = chunks[0][0].size();

    // Initialize M_in with 9 blocks of nrows x ncols complex numbers
    std::vector<std::vector<std::vector<std::complex<double>>>> M_in(9, 
        std::vector<std::vector<std::complex<double>>>(nrows, 
        std::vector<std::complex<double>>(ncols, std::complex<double>(0.0, 0.0))));

    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < ncols; ++j) {
            M_in[0][i][j] = chunks[0][i][j]; // t11_T1
            M_in[1][i][j] = chunks[1][i][j]; // t12_T1 real part
            M_in[2][i][j] = chunks[2][i][j]; // t12_T1 imaginary part
            
            // Extract real and imaginary parts to create a new complex number
            double real_part1 = chunks[1][i][j].real();
            double imag_part1 = chunks[2][i][j].real();
            std::complex<double> complex_num1(real_part1, imag_part1);
            M_in[3][i][j] = std::conj(complex_num1); // t21_T1
            
            M_in[4][i][j] = chunks[3][i][j]; // t22_T1
            M_in[5][i][j] = chunks[4][i][j]; // t23_T1 real part
            M_in[6][i][j] = chunks[5][i][j]; // t23_T1 imaginary part
            
            // Extract real and imaginary parts to create a new complex number
            double real_part2 = chunks[4][i][j].real();
            double imag_part2 = chunks[5][i][j].real();
            std::complex<double> complex_num2(real_part2, imag_part2);
            M_in[7][i][j] = std::conj(complex_num2); // t32_T1
            
            M_in[8][i][j] = chunks[6][i][j]; // t33_T1
        }
    }

    return M_in;
}

std::vector<std::vector<std::vector<std::complex<double>>>> compose_C2_chunks(
    const std::vector<std::vector<std::vector<std::complex<double>>>>& chunks) {
    int nrows = chunks[0].size();
    int ncols = chunks[0][0].size();
    std::vector<std::vector<std::vector<std::complex<double>>>> M_in(4, 
        std::vector<std::vector<std::complex<double>>>(nrows, 
        std::vector<std::complex<double>>(ncols, {0.0, 0.0})));

    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < ncols; ++j) {
            M_in[0][i][j] = chunks[0][i][j]; // t11_T1
            M_in[1][i][j] = chunks[1][i][j]; // t12_T1 real part
            M_in[2][i][j] = chunks[2][i][j]; // t12_T1 imaginary part
            M_in[3][i][j] = chunks[3][i][j]; // t22_T1
        }
    }

    return M_in;
}

void transpose_M_in(std::vector<std::vector<std::vector<std::complex<double>>>>& M_in) {
    int depth = M_in.size();
    int rows = M_in[0].size();
    int cols = M_in[0][0].size();

    std::vector<std::vector<std::vector<std::complex<double>>>> transposed(depth, 
        std::vector<std::vector<std::complex<double>>>(cols, 
        std::vector<std::complex<double>>(rows, {0.0, 0.0})));

    for (int d = 0; d < depth; ++d) {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                transposed[d][j][i] = M_in[d][i][j];
            }
        }
    }

    M_in = std::move(transposed);
}

std::vector<std::vector<std::vector<double>>> extract_C3_filtered_chunks(
    const std::vector<std::vector<std::vector<std::complex<double>>>>& M_out) {
    // Implement extraction for C3 chunks here.
    // Placeholder: just return an empty 3D vector of double.
    return std::vector<std::vector<std::vector<double>>>(M_out.size(), 
        std::vector<std::vector<double>>(M_out[0].size(), 
        std::vector<double>(M_out[0][0].size(), 0.0)));
}

std::vector<std::vector<std::vector<double>>> extract_C2_filtered_chunks(
    const std::vector<std::vector<std::vector<std::complex<double>>>>& M_out) {
    // Implement extraction for C2 chunks here.
    // Placeholder: just return an empty 3D vector of double.
    return std::vector<std::vector<std::vector<double>>>(M_out.size(), 
        std::vector<std::vector<double>>(M_out[0].size(), 
        std::vector<double>(M_out[0][0].size(), 0.0)));
}


std::vector<std::vector<std::vector<float>>> make_Mask(int window_size) {
    std::vector<std::vector<std::vector<float>>> Mask(8, 
        std::vector<std::vector<float>>(window_size, 
        std::vector<float>(window_size, 0.0f)));

    int Nmax = 0;
    for (int k = 0; k < window_size; ++k) {
        for (int l = (window_size - 1) / 2; l < window_size; ++l) {
            Mask[Nmax][k][l] = 1.0f;
        }
    }

    Nmax = 4;
    for (int k = 0; k < window_size; ++k) {
        for (int l = 0; l <= (window_size - 1) / 2; ++l) {
            Mask[Nmax][k][l] = 1.0f;
        }
    }

    Nmax = 1;
    for (int k = 0; k < window_size; ++k) {
        for (int l = k; l < window_size; ++l) {
            Mask[Nmax][k][l] = 1.0f;
        }
    }

    Nmax = 5;
    for (int k = 0; k < window_size; ++k) {
        for (int l = 0; l <= k; ++l) {
            Mask[Nmax][k][l] = 1.0f;
        }
    }

    Nmax = 2;
    for (int k = 0; k <= (window_size - 1) / 2; ++k) {
        for (int l = 0; l < window_size; ++l) {
            Mask[Nmax][k][l] = 1.0f;
        }
    }

    Nmax = 6;
    for (int k = (window_size - 1) / 2; k < window_size; ++k) {
        for (int l = 0; l < window_size; ++l) {
            Mask[Nmax][k][l] = 1.0f;
        }
    }

    Nmax = 3;
    for (int k = 0; k < window_size; ++k) {
        for (int l = 0; l < window_size - k; ++l) {
            Mask[Nmax][k][l] = 1.0f;
        }
    }

    Nmax = 7;
    for (int k = 0; k < window_size; ++k) {
        for (int l = window_size - 1 - k; l < window_size; ++l) {
            Mask[Nmax][k][l] = 1.0f;
        }
    }

    return Mask;
}

std::pair<std::vector<std::vector<std::complex<float>>>, std::vector<std::vector<int>>> make_Coeff(
    float sigma2,
    int Deplct,
    int Nwindow_size,
    int window_sizeM1S2,
    int Sub_Nlig,
    int Sub_Ncol,
    const std::vector<std::vector<float>>& span,
    const std::vector<std::vector<std::vector<float>>>& Mask
) {
    std::vector<std::vector<std::complex<float>>> coeff(Sub_Nlig, 
        std::vector<std::complex<float>>(Sub_Ncol, {0.0f, 0.0f}));
    std::vector<std::vector<int>> Nmax(Sub_Nlig, std::vector<int>(Sub_Ncol, 0));

    for (int lig = 0; lig < Sub_Nlig; ++lig) {
        for (int col = 0; col < Sub_Ncol; ++col) {
            // 3x3 average SPAN sub-window calculation
            std::vector<std::vector<std::complex<float>>> subwin(3, std::vector<std::complex<float>>(3, {0.0f, 0.0f}));

            for (int k = 0; k < 3; ++k) {
                for (int l = 0; l < 3; ++l) {
                    std::complex<float> sum_subwin = {0.0f, 0.0f};
                    for (int kk = 0; kk < Nwindow_size; ++kk) {
                        for (int ll = 0; ll < Nwindow_size; ++ll) {
                            int idx_k = k * Deplct + kk + lig;
                            int idx_l = l * Deplct + ll + col;
                            sum_subwin += span[idx_k][idx_l] / float(Nwindow_size * Nwindow_size);
                        }
                    }
                    subwin[k][l] = sum_subwin;
                }
            }

            // Directional gradient computation
            std::vector<std::complex<float>> Dist(4, {0.0f, 0.0f});
            Dist[0] = -subwin[0][0] + subwin[0][2] - subwin[1][0] + subwin[1][2] - subwin[2][0] + subwin[2][2];
            Dist[1] =  subwin[0][1] + subwin[0][2] - subwin[1][0] + subwin[1][2] - subwin[2][0] - subwin[2][1];
            Dist[2] =  subwin[0][0] + subwin[0][1] + subwin[0][2] - subwin[2][0] - subwin[2][1] - subwin[2][2];
            Dist[3] =  subwin[0][0] + subwin[0][1] + subwin[1][0] - subwin[1][2] - subwin[2][1] - subwin[2][2];

            // Choice of directional mask based on gradient
            float MaxDist = -std::numeric_limits<float>::infinity();
            int Nmax_lig_col = 0;

            for (int k = 0; k < 4; ++k) {
                if (MaxDist < std::abs(Dist[k])) {
                    MaxDist = std::abs(Dist[k]);
                    Nmax_lig_col = k;
                }
            }

            if (Dist[Nmax_lig_col].real() > 0.0f) {
                Nmax_lig_col += 4;
            }
            Nmax[lig][col] = Nmax_lig_col;

            // Calculate span statistics
            float m_span = 0.0f, m_span2 = 0.0f, Npoints = 0.0f;

            for (int k = -window_sizeM1S2; k <= window_sizeM1S2; ++k) {
                for (int l = -window_sizeM1S2; l <= window_sizeM1S2; ++l) {
                    float mask_value = Mask[Nmax_lig_col][window_sizeM1S2 + k][window_sizeM1S2 + l];
                    if (mask_value == 1.0f) {
                        int idx_k = window_sizeM1S2 + k + lig;
                        int idx_l = window_sizeM1S2 + l + col;
                        float s = span[idx_k][idx_l];
                        m_span += s;
                        m_span2 += s * s;
                        Npoints += 1.0f;
                    }
                }
            }

            m_span /= Npoints;
            m_span2 /= Npoints;

            // Variance and coefficient computation
            float v_span = m_span2 - m_span * m_span;
            float cv_span = std::sqrt(std::abs(v_span)) / (1e-8f + m_span);
            float coeff_lig_col = (cv_span * cv_span - sigma2) / (cv_span * cv_span * (1.0f + sigma2) + 1e-8f);

            if (coeff_lig_col < 0.0f) {
                coeff_lig_col = 0.0f;
            }
            coeff[lig][col] = coeff_lig_col;
        }
    }

    return {coeff, Nmax};
}

std::vector<std::vector<std::vector<double>>> process_chunk_rfleecpp(
    const std::vector<std::vector<std::vector<std::complex<double>>>>& chunks,
    int window_size
) {
    int Nlook = 1;
    // Ensure at least one look
    if (Nlook <= 0) Nlook = 1;

    // Adjust window size if even
    if (window_size % 2 == 0) window_size += 1;

    int pad_top_left = window_size / 2;
    int pad_bottom_right = pad_top_left + 1;

    // Check the number of chunks
    int num_chunks = chunks.size();

    std::string PolTypeOut;
    int NpolarOut;
    std::vector<std::vector<std::vector<std::complex<double>>>> M_in;

    // Prepare input chunks and polarimetric type
    if (num_chunks == 9) {
        // Compose M_in for "C3"
        PolTypeOut = "C3";
        NpolarOut = 9;
        M_in = compose_C3_chunks(chunks);
    } else if (num_chunks == 4) {
        // Compose M_in for "C2"
        PolTypeOut = "C2";
        NpolarOut = 4;
        M_in = compose_C2_chunks(chunks);
    } else {
        throw std::invalid_argument("Unsupported number of chunks");
    }

    // Transpose the input matrix (axes: [2,0,1])
    // transpose_M_in(M_in);

    // Initialize output matrix
    std::vector<std::vector<std::vector<std::complex<double>>>> M_out(NpolarOut,
        std::vector<std::vector<std::complex<double>>>(M_in[0].size(), 
        std::vector<std::complex<double>>(M_in[0][0].size(), {0.0, 0.0})));

    // Get dimensions of the padded matrix
    int Nlig_padded = M_in[0].size();
    int Ncol_padded = M_in[0][0].size();
    int Nlig = Nlig_padded - window_size;
    int Ncol = Ncol_padded - window_size;

    int window_sizeM1S2 = (window_size - 1) / 2;

    // Initialize Valid and span arrays
    std::vector<std::vector<float>> Valid(Nlig_padded, std::vector<float>(Ncol_padded, 1.0f));
    // std::vector<std::vector<float>> span(Nlig_padded, std::vector<float>(Ncol_padded, 0.0f));

    // Speckle variance
    double sigma2 = 1.0 / Nlook;
    auto [Nwindow_size, Deplct] = get_window_params(window_size);
    // Create mask and coefficients
    auto Mask = make_Mask(window_size);

    Matrix span = compute_span(PolTypeOut, M_in);

    auto [coeff, Nmax] = make_Coeff(sigma2, Deplct, Nwindow_size, window_sizeM1S2, Nlig, Ncol, span, Mask);

    // Initialize output matrix
    // std::vector<std::vector<std::vector<std::complex<double>>>> M_out(NpolarOut,
    //     std::vector<std::vector<std::complex<double>>>(Nlig, std::vector<std::complex<double>>(Ncol, {0.0, 0.0})));

    // Process each pixel
    for (int lig = 0; lig < Nlig; ++lig) {
        for (int col = 0; col < Ncol; ++col) {
            if (Valid[window_sizeM1S2 + lig][window_sizeM1S2 + col] == 1.0f) {
                for (int Np = 0; Np < NpolarOut; ++Np) {
                    std::complex<double> mean = 0.0;
                    double Npoints = 0.0;

                    for (int k = -window_sizeM1S2; k <= window_sizeM1S2; ++k) {
                        for (int l = -window_sizeM1S2; l <= window_sizeM1S2; ++l) {
                            int idx_k = window_sizeM1S2 + lig + k;
                            int idx_l = window_sizeM1S2 + col + l;

                            if (Mask[Nmax[lig][col]][window_sizeM1S2 + k][window_sizeM1S2 + l] == 1.0f) {
                                mean += M_in[Np][idx_k][idx_l];
                                Npoints += 1.0;
                            }
                        }
                    }

                    if (Npoints > 0.0) mean /= Npoints;

                    std::complex<double> center_pixel = M_in[Np][window_sizeM1S2 + lig][window_sizeM1S2 + col];
                    // M_out[Np][lig][col] = mean + coeff[lig][col] * (center_pixel - mean);
                    M_out[Np][lig][col] = mean + static_cast<std::complex<double>>(coeff[lig][col]) * (center_pixel - mean);

                }
            }
        }
    }

    std::vector<std::vector<std::vector<double>>> filtered_chunks;

    // Extract the filtered chunks based on the polarimetric type
    if (PolTypeOut == "C3") {
        filtered_chunks = extract_C3_filtered_chunks(M_out);
    } else if (PolTypeOut == "C2") {
        filtered_chunks = extract_C2_filtered_chunks(M_out);
    }

    return filtered_chunks;
}

PYBIND11_MODULE(rflee, m) {
    m.def("process_chunk_rfleecpp", &process_chunk_rfleecpp, "Apply refined Lee filter to complex chunk");
}
