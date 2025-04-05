#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>  // This is required for automatic STL conversions
#include <vector>
#include <stdexcept>
#include <cmath>  // For std::floor and std::ceil

namespace py = pybind11;

// Function to compute the sliding window mean of a 2D numpy array
py::array_t<double> sliding_window_mean(const py::array_t<double>& arr, int window_size) {
    py::buffer_info buf_info = arr.request();
    int rows = buf_info.shape[0];
    int cols = buf_info.shape[1];

    // Create a result array to store the sliding window means
    auto result = py::array_t<double>(buf_info.shape);
    py::buffer_info buf_info_result = result.request();
    double* ptr_result = static_cast<double*>(buf_info_result.ptr);

    // Iterate over the array, apply the sliding window, and compute the mean
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            // Get the window bounds
            int row_start = std::max(i - window_size / 2, 0);
            int row_end = std::min(i + window_size / 2 + 1, rows);
            int col_start = std::max(j - window_size / 2, 0);
            int col_end = std::min(j + window_size / 2 + 1, cols);

            // Sum the elements within the window and compute the mean
            double sum = 0.0;
            int count = 0;
            for (int r = row_start; r < row_end; ++r) {
                for (int c = col_start; c < col_end; ++c) {
                    sum += static_cast<double*>(buf_info.ptr)[r * cols + c];
                    ++count;
                }
            }
            ptr_result[i * cols + j] = sum / count;
        }
    }

    return result;
}

// Function to sum the sliding window means of a list of 2D numpy arrays
py::array_t<double> sum_filt(const std::vector<py::array_t<double>>& arrays, int window_size,  
    const std::vector<std::string>& input_filepaths,
    double chi_in,
    double psi_in) {
    if (arrays.empty()) {
        throw std::runtime_error("Input array list is empty");
    }

    py::buffer_info buf_info1 = arrays[0].request();
    int rows = buf_info1.shape[0];
    int cols = buf_info1.shape[1];

    // Ensure that all arrays have the same shape
    for (const auto& arr : arrays) {
        py::buffer_info buf_info = arr.request();
        if (buf_info.shape[0] != rows || buf_info.shape[1] != cols) {
            throw std::runtime_error("All arrays must have the same shape");
        }
    }

    // Create a result array of the same shape (2D array)
    auto result = py::array_t<double>(buf_info1.shape);
    py::buffer_info buf_info_result = result.request();
    double* ptr_result = static_cast<double*>(buf_info_result.ptr);

    // Initialize the result array to zero
    for (std::size_t i = 0; i < rows * cols; ++i) {
        ptr_result[i] = 0.0;
    }

    // Sum the sliding window means for each array
    for (const auto& arr : arrays) {
        // Compute the sliding window mean for the current array
        auto window_mean = sliding_window_mean(arr, window_size);
        py::buffer_info buf_info_window = window_mean.request();
        double* ptr_window = static_cast<double*>(buf_info_window.ptr);

        // Add the mean values to the result
        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                ptr_result[i * cols + j] += ptr_window[i * cols + j];
            }
        }
    }

    return result;
}

PYBIND11_MODULE(testcprvi, m) {
    m.def("sum_filt", &sum_filt, 
        "Sum a list of 2D numpy arrays element-wise after computing sliding window means");
}
