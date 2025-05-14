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


    // Initialize output matrix
    std::vector<std::vector<std::vector<std::complex<double>>>> M_out(NpolarOut,
        std::vector<std::vector<std::complex<double>>>(M_in[0].size(), 
        std::vector<std::complex<double>>(M_in[0][0].size(), {0.0, 0.0})));

    std::vector<std::vector<std::vector<double>>> filtered_chunks;
    

    // Iterate over M_in to separate the real and imaginary parts into filtered_chunks
    for (size_t i = 0; i < M_out.size(); ++i) {
        // For each chunk, create 2 separate layers for real and imaginary parts
        std::vector<std::vector<double>> real_layer(num_rows, std::vector<double>(num_cols, 0.0));
        std::vector<std::vector<double>> imag_layer(num_rows, std::vector<double>(num_cols, 0.0));

        for (size_t j = 0; j < M_out[i].size(); ++j) {
            for (size_t k = 0; k < M_out[i][j].size(); ++k) {
                // Assign the real and imaginary parts to their respective layers
                real_layer[j][k] = M_out[i][j][k].real();
                imag_layer[j][k] = M_out[i][j][k].imag();
            }
        }

        // Add the real and imaginary layers to the filtered_chunks
        filtered_chunks.push_back(real_layer);
        filtered_chunks.push_back(imag_layer);
    }

    // std::cout << "Size of filtered_chunks: " << filtered_chunks.size() << " x "
    //           << filtered_chunks[0].size() << " x " << filtered_chunks[0][0].size() << "\n";

    return filtered_chunks;
}

PYBIND11_MODULE(rflee, m) {
    m.def("process_chunk_rfleecpp", &process_chunk_rfleecpp, "Apply refined Lee filter to complex chunk");
}
