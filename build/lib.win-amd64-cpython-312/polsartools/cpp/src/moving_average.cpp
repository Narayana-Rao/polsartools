#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>

namespace py = pybind11;

// Function to apply a moving average filter
std::vector<std::vector<float>> moving_average(const std::vector<std::vector<float>>& input, int window_size) {
    int rows = input.size();
    int cols = input[0].size();
    int half_window = window_size / 2;
    std::vector<std::vector<float>> output(rows, std::vector<float>(cols, 0.0f));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            float sum = 0.0f;
            int count = 0;

            for (int k = -half_window; k <= half_window; k++) {
                for (int l = -half_window; l <= half_window; l++) {
                    int ni = i + k;
                    int nj = j + l;
                    if (ni >= 0 && ni < rows && nj >= 0 && nj < cols) {
                        sum += input[ni][nj];
                        count++;
                    }
                }
            }

            output[i][j] = sum / count;
        }
    }

    return output;
}

// Function to apply moving average filter to a list of 2D arrays
py::list apply_moving_average(py::list input_arrays, int window_size) {
    py::list output_arrays;
    for (auto array : input_arrays) {
        py::array_t<float> numpy_array = array.cast<py::array_t<float>>();
        py::buffer_info buf = numpy_array.request();

        if (buf.ndim != 2) {
            throw std::runtime_error("Each input must be a 2D array.");
        }

        int rows = buf.shape[0];
        int cols = buf.shape[1];
        std::vector<std::vector<float>> matrix(rows, std::vector<float>(cols));

        float* ptr = static_cast<float*>(buf.ptr);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                matrix[i][j] = ptr[i * cols + j];
            }
        }

        auto smoothed_matrix = moving_average(matrix, window_size);

        py::array_t<float> output_array({rows, cols});
        py::buffer_info out_buf = output_array.request();
        float* out_ptr = static_cast<float*>(out_buf.ptr);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                out_ptr[i * cols + j] = smoothed_matrix[i][j];
            }
        }

        output_arrays.append(output_array);
    }

    return output_arrays;
}

// Pybind11 module
PYBIND11_MODULE(moving_average, m) {
    m.doc() = "Moving Average Filter module";
    m.def("apply_moving_average", &apply_moving_average, "Apply moving average filtering", py::arg("input_arrays"), py::arg("window_size"));
}
