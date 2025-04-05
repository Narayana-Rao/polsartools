#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>  // This is required for automatic STL conversions
#include <vector>
#include <stdexcept>
#include <cmath>  // For std::floor and std::ceil

namespace py = pybind11;

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
        py::buffer_info buf_info_window = arr.request();
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
