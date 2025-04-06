#include <vector>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <cstring> // for memset

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>


namespace py = pybind11;
using namespace std;

using cfloat = complex<float>;

typedef std::complex<float> cpx;

// Helper to access 2D array
inline int idx(int row, int col, int width) {
    return row * width + col;
}


vector<vector<vector<float>>> make_Mask(int window_size) {
    vector<vector<vector<float>>> Mask(8, vector<vector<float>>(window_size, vector<float>(window_size, 0.0f)));
    int half = (window_size - 1) / 2;

    for (int k = 0; k < window_size; ++k)
        for (int l = half; l < window_size; ++l)
            Mask[0][k][l] = 1.0f;

    for (int k = 0; k < window_size; ++k)
        for (int l = 0; l <= half; ++l)
            Mask[4][k][l] = 1.0f;

    for (int k = 0; k < window_size; ++k)
        for (int l = k; l < window_size; ++l)
            Mask[1][k][l] = 1.0f;

    for (int k = 0; k < window_size; ++k)
        for (int l = 0; l <= k; ++l)
            Mask[5][k][l] = 1.0f;

    for (int k = 0; k <= half; ++k)
        for (int l = 0; l < window_size; ++l)
            Mask[2][k][l] = 1.0f;

    for (int k = half; k < window_size; ++k)
        for (int l = 0; l < window_size; ++l)
            Mask[6][k][l] = 1.0f;

    for (int k = 0; k < window_size; ++k)
        for (int l = 0; l < window_size - k; ++l)
            Mask[3][k][l] = 1.0f;

    for (int k = 0; k < window_size; ++k)
        for (int l = window_size - 1 - k; l < window_size; ++l)
            Mask[7][k][l] = 1.0f;

    return Mask;
}

pair<py::array_t<cfloat>, py::array_t<int>> make_Coeff(
    float sigma2, int Deplct, int Nwindow_size, int window_sizeM1S2,
    int Sub_Nlig, int Sub_Ncol, py::array_t<float> span,
    const vector<vector<vector<float>>> &Mask
) {
    auto buf_span = span.unchecked<2>();
    py::array_t<cfloat> coeff({Sub_Nlig, Sub_Ncol});
    py::array_t<int> Nmax({Sub_Nlig, Sub_Ncol});
    auto buf_coeff = coeff.mutable_unchecked<2>();
    // auto buf_Nmax = Nmax.mutable_unchecked<2>();
    auto buf_Nmax = Nmax.unchecked<2>();
    int mask_idx = buf_Nmax(lig, col);

    for (int lig = 0; lig < Sub_Nlig; ++lig) {
        for (int col = 0; col < Sub_Ncol; ++col) {
            float subwin[3][3] = {};
            for (int k = 0; k < 3; ++k) {
                for (int l = 0; l < 3; ++l) {
                    float sum = 0.0f;
                    for (int kk = 0; kk < Nwindow_size; ++kk) {
                        for (int ll = 0; ll < Nwindow_size; ++ll) {
                            int idx_k = k * Deplct + kk + lig;
                            int idx_l = l * Deplct + ll + col;
                            sum += buf_span(idx_k, idx_l) / (Nwindow_size * Nwindow_size);
                        }
                    }
                    subwin[k][l] = sum;
                }
            }

            float Dist[4] = {};
            Dist[0] = -subwin[0][0] + subwin[0][2] - subwin[1][0] + subwin[1][2] - subwin[2][0] + subwin[2][2];
            Dist[1] = subwin[0][1] + subwin[0][2] - subwin[1][0] + subwin[1][2] - subwin[2][0] - subwin[2][1];
            Dist[2] = subwin[0][0] + subwin[0][1] + subwin[0][2] - subwin[2][0] - subwin[2][1] - subwin[2][2];
            Dist[3] = subwin[0][0] + subwin[0][1] + subwin[1][0] - subwin[1][2] - subwin[2][1] - subwin[2][2];

            float MaxDist = -1e10;
            int Nmax_lig_col = 0;
            for (int k = 0; k < 4; ++k) {
                if (fabs(Dist[k]) > MaxDist) {
                    MaxDist = fabs(Dist[k]);
                    Nmax_lig_col = k;
                }
            }
            if (Dist[Nmax_lig_col] > 0)
                Nmax_lig_col += 4;

            buf_Nmax(lig, col) = Nmax_lig_col;

            float m_span = 0.0f, m_span2 = 0.0f, Npoints = 0.0f;
            for (int k = -window_sizeM1S2; k <= window_sizeM1S2; ++k) {
                for (int l = -window_sizeM1S2; l <= window_sizeM1S2; ++l) {
                    if (Mask[Nmax_lig_col][window_sizeM1S2 + k][window_sizeM1S2 + l] == 1.0f) {
                        int idx_k = lig + window_sizeM1S2 + k;
                        int idx_l = col + window_sizeM1S2 + l;
                        float s = buf_span(idx_k, idx_l);
                        m_span += s;
                        m_span2 += s * s;
                        Npoints += 1.0f;
                    }
                }
            }

            m_span /= Npoints;
            m_span2 /= Npoints;
            float v_span = m_span2 - m_span * m_span;
            float cv_span = sqrtf(fabs(v_span)) / (1e-8f + m_span);
            float coeff_val = (cv_span * cv_span - sigma2) / (cv_span * cv_span * (1 + sigma2) + 1e-8f);
            buf_coeff(lig, col) = coeff_val < 0.0f ? 0.0f : coeff_val;
        }
    }

    return {coeff, Nmax};
}

std::vector<std::vector<float>> process_chunk_refined_lee(
    const std::vector<const float*>& chunks,
    int window_size,
    int height,
    int width
) {
    int Nlook = 1;
    std::string PolTypeOut;
    int NpolarOut;

    // int padded_height = height + window_size;
    // int padded_width = width + window_size;
    int padded_height = height + window_size/2;
    int padded_width = width + window_size/2;

    int window_half = (window_size - 1) / 2;

    std::vector<cpx*> M_in_vec;

    if (chunks.size() == 9) {
        // C3 case
        PolTypeOut = "C3";
        NpolarOut = 9;

        M_in_vec.resize(9);
        for (int i = 0; i < 9; ++i) {
            M_in_vec[i] = new cpx[padded_height * padded_width];
        }

        for (int r = 0; r < padded_height; ++r) {
            for (int c = 0; c < padded_width; ++c) {
                int id = idx(r, c, padded_width);
                M_in_vec[0][id] = chunks[0][id];
                M_in_vec[1][id] = cpx(chunks[1][id], chunks[2][id]);
                M_in_vec[2][id] = cpx(chunks[3][id], chunks[4][id]);
                M_in_vec[3][id] = std::conj(M_in_vec[1][id]);
                M_in_vec[4][id] = chunks[5][id];
                M_in_vec[5][id] = cpx(chunks[6][id], chunks[7][id]);
                M_in_vec[6][id] = std::conj(M_in_vec[2][id]);
                M_in_vec[7][id] = std::conj(M_in_vec[5][id]);
                M_in_vec[8][id] = chunks[8][id];
            }
        }

    } else if (chunks.size() == 4) {
        // C2 case
        PolTypeOut = "C2";
        NpolarOut = 4;

        M_in_vec.resize(4);
        for (int i = 0; i < 4; ++i) {
            M_in_vec[i] = new cpx[padded_height * padded_width];
        }

        for (int r = 0; r < padded_height; ++r) {
            for (int c = 0; c < padded_width; ++c) {
                int id = idx(r, c, padded_width);
                M_in_vec[0][id] = chunks[0][id];
                M_in_vec[1][id] = cpx(chunks[1][id], chunks[2][id]);
                M_in_vec[2][id] = std::conj(M_in_vec[1][id]);
                M_in_vec[3][id] = chunks[3][id];
            }
        }
    } else {
        throw std::runtime_error("Invalid number of chunks. Must be 4 or 9.");
    }

    float sigma2 = 1.0f / Nlook;
    auto Mask = make_Mask(window_size); // returns 8 x window x window
    float* span = new float[padded_height * padded_width];

    // Compute SPAN
    // for (int r = 0; r < padded_height; ++r) {
    //     for (int c = 0; c < padded_width; ++c) {
    //         int id = idx(r, c, padded_width);
    //         if (PolTypeOut == "C2") {
    //             span[id] = std::real(M_in_vec[0][id]) + std::real(M_in_vec[3][id]);
    //         } else {
    //             span[id] = std::real(M_in_vec[0][id]) + std::real(M_in_vec[4][id]) + std::real(M_in_vec[8][id]);
    //         }
    //     }
    // }

    py::array_t<float> span_arr({padded_height, padded_width});
    auto buf_span = span_arr.mutable_unchecked<2>();
    
    for (int r = 0; r < padded_height; ++r) {
        for (int c = 0; c < padded_width; ++c) {
            int id = idx(r, c, padded_width);
            buf_span(r, c) = (PolTypeOut == "C2")
                             ? std::real(M_in_vec[0][id]) + std::real(M_in_vec[3][id])
                             : std::real(M_in_vec[0][id]) + std::real(M_in_vec[4][id]) + std::real(M_in_vec[8][id]);
        }
    }

    // Coefficients
    int output_height = padded_height - window_size;
    int output_width = padded_width - window_size;
    // auto [coeff, Nmax] = make_Coeff(sigma2, window_size == 3 ? 1 : 2, 3, window_half,
    //                                 output_height, output_width, span, Mask);
    auto [coeff, Nmax] = make_Coeff(sigma2, window_size == 3 ? 1 : 2, 3, window_half,
        output_height, output_width, span_arr, Mask);

    // Filtered output
    std::vector<std::vector<float>> filtered_chunks(NpolarOut);
    for (int i = 0; i < NpolarOut; ++i) {
        filtered_chunks[i].resize(output_height * output_width);
    }

    // Main filtering loop
    for (int lig = 0; lig < output_height; ++lig) {
        for (int col = 0; col < output_width; ++col) {
            int center_r = lig + window_half;
            int center_c = col + window_half;
            int id_center = idx(center_r, center_c, padded_width);
            int id_out = idx(lig, col, output_width);
            int mask_idx = Nmax[lig * output_width + col];

            for (int p = 0; p < NpolarOut; ++p) {
                cpx mean = 0.0f;
                float Npoints = 0.0f;

                for (int k = -window_half; k <= window_half; ++k) {
                    for (int l = -window_half; l <= window_half; ++l) {
                        int mask_val = Mask[mask_idx][window_half + k][window_half + l];
                        if (mask_val == 1.0f) {
                            int r = center_r + k;
                            int c = center_c + l;
                            int id = idx(r, c, padded_width);
                            mean += M_in_vec[p][id];
                            Npoints += 1.0f;
                        }
                    }
                }

                mean /= Npoints;
                cpx center_val = M_in_vec[p][id_center];
                auto buf_coeff = coeff.unchecked<2>();
                cpx filtered = mean + buf_coeff(lig, col) * (center_val - mean);
                // cpx filtered = mean + coeff[lig * output_width + col] * (center_val - mean);

                filtered_chunks[p][id_out] = (p == 1 || p == 2 || p == 5 || p == 7) ? filtered.imag() : filtered.real();
            }
        }
    }

    // Cleanup
    for (auto ptr : M_in_vec) delete[] ptr;
    delete[] span;

    return filtered_chunks;
}

PYBIND11_MODULE(rflee, m) {
    m.doc() = "Refined Lee Filter module with pointer-based C++ implementation";
    m.def("process_chunk_refined_lee", &process_chunk_refined_lee, "Apply refined Lee filter to SAR image chunk");
}
