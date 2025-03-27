#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>

namespace py = pybind11;

// Function to sum a list of 2D arrays
py::list sum_2d_arrays(py::list input_arrays) {
    if (input_arrays.size() == 0) {
        throw std::runtime_error("Input list is empty.");
    }

    py::array_t<float> first_array = input_arrays[0].cast<py::array_t<float>>();
    py::buffer_info buf = first_array.request();

    if (buf.ndim != 2) {
        throw std::runtime_error("Each input must be a 2D array.");
    }

    int rows = buf.shape[0];
    int cols = buf.shape[1];
    std::vector<std::vector<float>> sum_matrix(rows, std::vector<float>(cols, 0.0f));

    // Iterate through input arrays and sum element-wise
    for (auto array : input_arrays) {
        py::array_t<float> numpy_array = array.cast<py::array_t<float>>();
        py::buffer_info buf = numpy_array.request();

        if (buf.ndim != 2 || buf.shape[0] != rows || buf.shape[1] != cols) {
            throw std::runtime_error("All arrays must have the same dimensions.");
        }

        float* ptr = static_cast<float*>(buf.ptr);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                sum_matrix[i][j] += ptr[i * cols + j];
            }
        }
    }

    // Convert sum_matrix to a NumPy array
    py::array_t<float> output_array({rows, cols});
    py::buffer_info out_buf = output_array.request();
    float* out_ptr = static_cast<float*>(out_buf.ptr);
    
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            out_ptr[i * cols + j] = sum_matrix[i][j];
        }
    }

    py::list output;
    output.append(output_array);
    return output;
}

// Pybind11 module
PYBIND11_MODULE(sum_arrays, m) {
    m.doc() = "Sum of 2D Arrays Module";
    m.def("sum_2d_arrays", &sum_2d_arrays, "Sum a list of 2D arrays and return as a list with one element", py::arg("input_arrays"));
}
