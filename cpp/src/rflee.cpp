#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

using namespace std;
namespace py = pybind11;

typedef complex<float> complexf;
typedef vector<vector<float>> MatrixF;
typedef vector<vector<complexf>> MatrixC;
typedef vector<vector<vector<complexf>>> CubeC;
typedef vector<vector<vector<float>>> CubeF;

MatrixC pad_matrix_c(const MatrixC& mat, int pad_top, int pad_bottom, int pad_left, int pad_right) {
    int rows = mat.size();
    int cols = mat[0].size();
    MatrixC padded(rows + pad_top + pad_bottom, vector<complexf>(cols + pad_left + pad_right, complexf(0.0f, 0.0f)));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            padded[i + pad_top][j + pad_left] = mat[i][j];
    return padded;
}

vector<vector<vector<float>>> make_Mask(int window_size) {
    vector<vector<vector<float>>> Mask(8, vector<vector<float>>(window_size, vector<float>(window_size, 0.0f)));
    int half = (window_size - 1) / 2;

    for (int i = 0; i < window_size; ++i) {
        for (int j = half; j < window_size; ++j)
            Mask[0][i][j] = 1.0f;
        for (int j = 0; j <= half; ++j)
            Mask[4][i][j] = 1.0f;
        for (int j = i; j < window_size; ++j)
            Mask[1][i][j] = 1.0f;
        for (int j = 0; j <= i; ++j)
            Mask[5][i][j] = 1.0f;
    }
    for (int i = 0; i <= half; ++i)
        for (int j = 0; j < window_size; ++j)
            Mask[2][i][j] = 1.0f;
    for (int i = half; i < window_size; ++i)
        for (int j = 0; j < window_size; ++j)
            Mask[6][i][j] = 1.0f;
    for (int i = 0; i < window_size; ++i)
        for (int j = 0; j < window_size - i; ++j)
            Mask[3][i][j] = 1.0f;
    for (int i = 0; i < window_size; ++i)
        for (int j = window_size - 1 - i; j < window_size; ++j)
            Mask[7][i][j] = 1.0f;

    return Mask;
}

pair<CubeF, vector<vector<int>>> make_Coeff(
    float sigma2, int Deplct, int Nwindow_size, int window_sizeM1S2,
    int Sub_Nlig, int Sub_Ncol, const MatrixF& span,
    const vector<vector<vector<float>>>& Mask) {

    CubeF coeff(Sub_Nlig, vector<vector<float>>(Sub_Ncol, vector<float>(1, 0.0f)));
    vector<vector<int>> Nmax(Sub_Nlig, vector<int>(Sub_Ncol, 0));

    for (int lig = 0; lig < Sub_Nlig; ++lig) {
        for (int col = 0; col < Sub_Ncol; ++col) {
            float subwin[3][3] = {0};
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l)
                    for (int kk = 0; kk < Nwindow_size; ++kk)
                        for (int ll = 0; ll < Nwindow_size; ++ll)
                            subwin[k][l] += span[lig + k * Deplct + kk][col + l * Deplct + ll] / (Nwindow_size * Nwindow_size);

            float Dist[4] = {
                -subwin[0][0] + subwin[0][2] - subwin[1][0] + subwin[1][2] - subwin[2][0] + subwin[2][2],
                 subwin[0][1] + subwin[0][2] - subwin[1][0] + subwin[1][2] - subwin[2][0] - subwin[2][1],
                 subwin[0][0] + subwin[0][1] + subwin[0][2] - subwin[2][0] - subwin[2][1] - subwin[2][2],
                 subwin[0][0] + subwin[0][1] + subwin[1][0] - subwin[1][2] - subwin[2][1] - subwin[2][2]
            };

            int Nmax_lig_col = distance(Dist, max_element(Dist, Dist + 4, [](float a, float b) { return abs(a) < abs(b); }));
            if (Dist[Nmax_lig_col] > 0.0f) Nmax_lig_col += 4;
            Nmax[lig][col] = Nmax_lig_col;

            float m_span = 0.0f, m_span2 = 0.0f;
            int Npoints = 0;
            for (int k = -window_sizeM1S2; k <= window_sizeM1S2; ++k) {
                for (int l = -window_sizeM1S2; l <= window_sizeM1S2; ++l) {
                    int i = window_sizeM1S2 + k;
                    int j = window_sizeM1S2 + l;
                    if (Mask[Nmax_lig_col][i][j] == 1.0f) {
                        float s = span[lig + i][col + j];
                        m_span += s;
                        m_span2 += s * s;
                        Npoints++;
                    }
                }
            }

            m_span /= Npoints;
            m_span2 /= Npoints;
            float v_span = m_span2 - m_span * m_span;
            float cv_span = sqrtf(fabsf(v_span)) / (1e-8f + m_span);
            float coeff_val = (cv_span * cv_span - sigma2) / (cv_span * cv_span * (1 + sigma2) + 1e-8f);
            coeff_val = max(coeff_val, 0.0f);
            coeff[lig][col][0] = coeff_val;
        }
    }

    return {coeff, Nmax};
}

std::vector<py::array_t<double>> process_chunk_rfleecpp(
    const std::vector<py::array_t<double>>& chunks,
    int window_size,
    const std::vector<std::string>& input_filepaths
    ) {
    int Deplct = 1;
    float sigma2 = 0.25
    if (window_size % 2 == 0) window_size += 1;
    int window_sizeM1S2 = (window_size - 1) / 2;
    int pad_top = window_sizeM1S2;
    int pad_bottom = window_sizeM1S2 + 1;
    int pad_left = window_sizeM1S2;
    int pad_right = window_sizeM1S2 + 1;

    CubeC chunk;
    for (const auto& arr : chunks) {
        auto buf = arr.unchecked<2>();
        MatrixC mat(buf.shape(0), vector<complexf>(buf.shape(1)));
        for (ssize_t i = 0; i < buf.shape(0); ++i)
            for (ssize_t j = 0; j < buf.shape(1); ++j)
                mat[i][j] = complexf(buf(i, j), 0.0f);
        chunk.push_back(pad_matrix_c(mat, pad_top, pad_bottom, pad_left, pad_right));
    }

    int Nwindow_size = 3;
    int nvec = chunk.size();
    int rows = chunk[0].size();
    int cols = chunk[0][0].size();

    MatrixF span(rows, vector<float>(cols, 0.0f));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            span[i][j] = real(chunk[0][i][j]) + real(chunk[2][i][j]);

    auto Mask = make_Mask(window_size);
    int Sub_Nlig = rows - window_size + 1;
    int Sub_Ncol = cols - window_size + 1;

    auto [coeff, Nmax] = make_Coeff(sigma2, Deplct, Nwindow_size, window_sizeM1S2, Sub_Nlig, Sub_Ncol, span, Mask);

    std::vector<py::array_t<double>> output_arrays;

    for (int v = 0; v < nvec; ++v) {
        py::array_t<double> out_array({Sub_Nlig, Sub_Ncol, 2});
        auto r = out_array.mutable_unchecked<3>();

        for (int lig = 0; lig < Sub_Nlig; ++lig) {
            for (int col = 0; col < Sub_Ncol; ++col) {
                int nmax = Nmax[lig][col];
                complexf moy(0.0f, 0.0f);
                int N = 0;
                for (int k = -window_sizeM1S2; k <= window_sizeM1S2; ++k) {
                    for (int l = -window_sizeM1S2; l <= window_sizeM1S2; ++l) {
                        int i = window_sizeM1S2 + k;
                        int j = window_sizeM1S2 + l;
                        if (Mask[nmax][i][j] == 1.0f) {
                            moy += chunk[v][lig + i][col + j];
                            N++;
                        }
                    }
                }
                moy /= static_cast<float>(N);
                complexf val = chunk[v][lig + window_sizeM1S2][col + window_sizeM1S2] * coeff[lig][col][0] + moy * (1.0f - coeff[lig][col][0]);
                r(lig, col, 0) = real(val);
                r(lig, col, 1) = imag(val);
            }
        }
        output_arrays.push_back(out_array);
    }

    return output_arrays;
}

PYBIND11_MODULE(rflee, m) {
    m.def("process_chunk_rfleecpp", &process_chunk_rfleecpp, "Apply refined Lee filter to complex chunk");
}
