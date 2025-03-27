#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>  // For handling Python lists

namespace py = pybind11;

// Function to process a list of 2D NumPy arrays
std::vector<py::array_t<double>> process_arrays(std::vector<py::array_t<double>> input_arrays) {
    std::vector<py::array_t<double>> output_arrays;

    for (auto& input_array : input_arrays) {
        // Request a buffer and get information about the array
        py::buffer_info buf = input_array.request();

        if (buf.ndim != 2) {
            throw std::runtime_error("Each input array must be 2D");
        }

        auto rows = buf.shape[0];
        auto cols = buf.shape[1];

        // Create an output array with the same shape
        auto output_array = py::array_t<double>(buf.shape);
        py::buffer_info out_buf = output_array.request();

        // Pointers to array data
        double* input_ptr = static_cast<double*>(buf.ptr);
        double* output_ptr = static_cast<double*>(out_buf.ptr);

        // Process the array (example: multiply each element by 2)
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                output_ptr[i * cols + j] = input_ptr[i * cols + j] * 2;
            }
        }

        output_arrays.push_back(output_array);  // Add the processed array to the output list
    }

    return output_arrays;  // Return the list of processed arrays
}

// Module definition
PYBIND11_MODULE(process_arrays, m) {
    m.doc() = "Efficiently process a list of 2D arrays using NumPy";  // Module docstring
    m.def("process_arrays", &process_arrays, "Process a list of 2D NumPy arrays");
}
