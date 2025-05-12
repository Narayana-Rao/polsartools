#include <pybind11/pybind11.h>
#include "processing.h"

namespace py = pybind11;

PYBIND11_MODULE(polsartools_python_bindings, m) {
    m.def("process", &process, "A sample processing function");
}
