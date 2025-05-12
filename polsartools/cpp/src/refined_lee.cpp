#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <vector>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace py = pybind11;

// Function to create directional masks
std::vector<std::vector<std::vector<float>>> make_Mask(int window_size) {
    std::vector<std::vector<std::vector<float>>> Mask(8, std::vector<std::vector<float>>(window_size, std::vector<float>(window_size, 0.0f)));
    
    int Nmax = 0;
    for (int k = 0; k < window_size; k++) {
        for (int l = (window_size - 1) / 2; l < window_size; l++) {
            Mask[Nmax][k][l] = 1.0f;
        }
    }
    
    Nmax = 4;
    for (int k = 0; k < window_size; k++) {
        for (int l = 0; l <= (window_size - 1) / 2; l++) {
            Mask[Nmax][k][l] = 1.0f;
        }
    }
    
    Nmax = 1;
    for (int k = 0; k < window_size; k++) {
        for (int l = k; l < window_size; l++) {
            Mask[Nmax][k][l] = 1.0f;
        }
    }
    
    Nmax = 5;
    for (int k = 0; k < window_size; k++) {
        for (int l = 0; l <= k; l++) {
            Mask[Nmax][k][l] = 1.0f;
        }
    }
    
    Nmax = 2;
    for (int k = 0; k <= (window_size - 1) / 2; k++) {
        for (int l = 0; l < window_size; l++) {
            Mask[Nmax][k][l] = 1.0f;
        }
    }
    
    Nmax = 6;
    for (int k = (window_size - 1) / 2; k < window_size; k++) {
        for (int l = 0; l < window_size; l++) {
            Mask[Nmax][k][l] = 1.0f;
        }
    }
    
    Nmax = 3;
    for (int k = 0; k < window_size; k++) {
        for (int l = 0; l < window_size - k; l++) {
            Mask[Nmax][k][l] = 1.0f;
        }
    }
    
    Nmax = 7;
    for (int k = 0; k < window_size; k++) {
        for (int l = window_size - 1 - k; l < window_size; l++) {
            Mask[Nmax][k][l] = 1.0f;
        }
    }
    
    return Mask;
}

// Function to compute coefficients for the refined Lee filter
std::pair<std::vector<std::vector<float>>, std::vector<std::vector<int>>> 
make_Coeff(float sigma2, int Deplct, int Nwindow_size, int window_sizeM1S2,
           int Sub_Nlig, int Sub_Ncol,
           const std::vector<std::vector<float>>& span,
           const std::vector<std::vector<std::vector<float>>>& Mask) {

    std::vector<std::vector<float>> coeff(Sub_Nlig, std::vector<float>(Sub_Ncol, 0.0f));
    std::vector<std::vector<int>> Nmax(Sub_Nlig, std::vector<int>(Sub_Ncol, 0));
    
    for (int lig = 0; lig < Sub_Nlig; lig++) {
        for (int col = 0; col < Sub_Ncol; col++) {
            std::vector<std::vector<float>> subwin(3, std::vector<float>(3, 0.0f));
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    float sum_subwin = 0.0f;
                    for (int kk = 0; kk < Nwindow_size; kk++) {
                        for (int ll = 0; ll < Nwindow_size; ll++) {
                            int idx_k = k * Deplct + kk + lig;
                            int idx_l = l * Deplct + ll + col;
                            sum_subwin += span[idx_k][idx_l] / (Nwindow_size * Nwindow_size);
                        }
                    }
                    subwin[k][l] = sum_subwin;
                }
            }

            std::vector<float> Dist(4, 0.0f);
            Dist[0] = -subwin[0][0] + subwin[0][2] - subwin[1][0] + subwin[1][2] - subwin[2][0] + subwin[2][2];
            Dist[1] = subwin[0][1] + subwin[0][2] - subwin[1][0] + subwin[1][2] - subwin[2][0] - subwin[2][1];
            Dist[2] = subwin[0][0] + subwin[0][1] + subwin[0][2] - subwin[2][0] - subwin[2][1] - subwin[2][2];
            Dist[3] = subwin[0][0] + subwin[0][1] + subwin[1][0] - subwin[1][2] - subwin[2][1] - subwin[2][2];

            float MaxDist = -std::numeric_limits<float>::infinity();
            int Nmax_lig_col = 0;
            for (int k = 0; k < 4; k++) {
                if (std::abs(Dist[k]) > MaxDist) {
                    MaxDist = std::abs(Dist[k]);
                    Nmax_lig_col = k;
                }
            }
            if (Dist[Nmax_lig_col] > 0.0f) {
                Nmax_lig_col += 4;
            }
            Nmax[lig][col] = Nmax_lig_col;

            float m_span = 0.0f, m_span2 = 0.0f, Npoints = 0.0f;
            for (int k = -window_sizeM1S2; k <= window_sizeM1S2; k++) {
                for (int l = -window_sizeM1S2; l <= window_sizeM1S2; l++) {
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
            float v_span = m_span2 - m_span * m_span;
            float cv_span = std::sqrt(std::abs(v_span)) / (1e-8f + m_span);
            // float coeff_lig_col = (cv_span * cv_span - sigma2) / (cv_span * cv_span * (1 + sigma2) + 1e-8f);
            float coeff_lig_col = (cv_span * cv_span - sigma2) / (std::max(cv_span * cv_span * (1 + sigma2), 1e-6f));
            coeff_lig_col = std::max(0.0f, coeff_lig_col);
            coeff[lig][col] = coeff_lig_col;
        }
    }

    return {coeff, Nmax};
}

// Main function to refine Lee filtering
py::list refined_lee(py::list chunks, int window_size) {
    // Parse the input chunks into C++ vectors
    size_t num_chunks = chunks.size();
    std::vector<std::vector<std::vector<std::complex<float>>>> M_in;
    
    if (num_chunks == 9 || num_chunks == 4) {
        for (size_t i = 0; i < num_chunks; i++) {
            py::array_t<std::complex<float>> chunk = chunks[i].cast<py::array_t<std::complex<float>>>();
            py::buffer_info buf = chunk.request();
            
            if (buf.ndim != 2) {
                throw std::runtime_error("Each chunk must be a 2D array.");
            }
            
            size_t rows = buf.shape[0];
            size_t cols = buf.shape[1];
            
            std::vector<std::vector<std::complex<float>>> chunk_vec(rows, std::vector<std::complex<float>>(cols, {0, 0}));
            std::complex<float>* ptr = static_cast<std::complex<float>*>(buf.ptr);
            
            for (size_t r = 0; r < rows; r++) {
                for (size_t c = 0; c < cols; c++) {
                    chunk_vec[r][c] = ptr[r * cols + c];
                }
            }
            
            M_in.push_back(chunk_vec);
        }
    } else {
        throw std::runtime_error("Number of chunks must be 4 or 9.");
    }
    
    // Transpose M_in for processing
    size_t NpolarOut = M_in.size();
    size_t Nlig_padded = M_in[0].size();
    size_t Ncol_padded = M_in[0][0].size();
    size_t Nlig = Nlig_padded - window_size;
    size_t Ncol = Ncol_padded - window_size;
    size_t window_sizeM1S2 = (window_size - 1) / 2;
    
    // Initialize Valid mask
    std::vector<std::vector<float>> Valid(Nlig_padded, std::vector<float>(Ncol_padded, 1.0f));
    
    // Speckle variance given by the number of looks
    float sigma2 = 1.0f; // Example value for Nlook = 1
    
    // Gradient window calculation parameters
    int Nwindow_size = 0, Deplct = 0;
    std::map<int, std::pair<int, int>> window_params = {
        {3, {1, 1}}, {5, {3, 1}}, {7, {3, 2}}, {9, {5, 2}},
        {11, {5, 3}}, {13, {5, 4}}, {15, {7, 4}}, {17, {7, 5}},
        {19, {7, 6}}, {21, {9, 6}}, {23, {9, 7}}, {25, {9, 8}},
        {27, {11, 8}}, {29, {11, 9}}, {31, {11, 10}}
    };
    
    if (window_params.find(window_size) != window_params.end()) {
        Nwindow_size = window_params[window_size].first;
        Deplct = window_params[window_size].second;
    } else {
        throw std::runtime_error("Invalid window size. Must be between 3 and 31.");
    }
    
    // Create Mask
    auto Mask = make_Mask(window_size);
    
    // Compute span
    std::vector<std::vector<float>> span(Nlig_padded, std::vector<float>(Ncol_padded, 0.0f));
    if (num_chunks == 4) {  // PolTypeOut == "C2"
        for (size_t i = 0; i < Nlig_padded; i++) {
            for (size_t j = 0; j < Ncol_padded; j++) {
                span[i][j] = std::real(M_in[0][i][j]) + std::real(M_in[3][i][j]);
                if (std::isnan(span[i][j])) {
                    std::cerr << "NaN detected in span at (" << i << ", " << j << ")" << std::endl;
                }
            }
        }
    } else if (num_chunks == 9) {  // PolTypeOut == "C3"
        for (size_t i = 0; i < Nlig_padded; i++) {
            for (size_t j = 0; j < Ncol_padded; j++) {
                span[i][j] = std::real(M_in[0][i][j]) + std::real(M_in[5][i][j]) + std::real(M_in[8][i][j]);
            }
        }
    } else {
        throw std::runtime_error("Unsupported PolTypeOut.");
    }
    
    // Compute coefficients
    auto [coeff, Nmax] = make_Coeff(sigma2, Deplct, Nwindow_size, window_sizeM1S2, Nlig, Ncol, span, Mask);
    
    // Initialize output
    std::vector<std::vector<std::vector<std::complex<float>>>> M_out(NpolarOut, 
        std::vector<std::vector<std::complex<float>>>(Nlig, std::vector<std::complex<float>>(Ncol, {0, 0})));
    
    // Filtering element by element
    for (size_t lig = 0; lig < Nlig; lig++) {
        for (size_t col = 0; col < Ncol; col++) {
            if (Valid[window_sizeM1S2 + lig][window_sizeM1S2 + col] == 1.0f) {
                for (size_t Np = 0; Np < NpolarOut; Np++) {
                    std::complex<float> mean = 0.0f;
                    float Npoints = 0.0f;
                    for (int k = -window_sizeM1S2; k <= window_sizeM1S2; k++) {
                        for (int l = -window_sizeM1S2; l <= window_sizeM1S2; l++) {
                            float mask_value = Mask[Nmax[lig][col]][window_sizeM1S2 + k][window_sizeM1S2 + l];
                            if (mask_value == 1.0f) {
                                int idx_k = window_sizeM1S2 + lig + k;
                                int idx_l = window_sizeM1S2 + col + l;
                                mean += M_in[Np][idx_k][idx_l];
                                Npoints += 1.0f;
                            }
                        }
                    }
                    mean /= Npoints;
                    auto center_pixel = M_in[Np][window_sizeM1S2 + lig][window_sizeM1S2 + col];
                    M_out[Np][lig][col] = mean + coeff[lig][col] * (center_pixel - mean);
                }
            }
        }
    }
    
    // Prepare output chunks
    py::list filtered_chunks;
    if (num_chunks == 9) {
        filtered_chunks.append(M_out[0]);  // Real part
        filtered_chunks.append(M_out[1]);  // Real part
        filtered_chunks.append(M_out[2]);  // Imaginary part
    }
    
    if (num_chunks == 4) {
        filtered_chunks.append(M_out[0]);
        filtered_chunks.append(M_out[1]);
        filtered_chunks.append(M_out[2]);
    }
    
    return filtered_chunks;
}


// Pybind11 module
PYBIND11_MODULE(refined_lee, m) {
    m.doc() = "Refined Lee Filtering module";
    m.def("refined_lee", &refined_lee, "Apply refined Lee filtering to input chunks", py::arg("chunks"), py::arg("window_size"));
}
