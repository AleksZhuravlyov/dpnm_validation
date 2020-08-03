#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "convection.h"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(convection, m) {
      m.def("calculate", &calculate,
            "velocities"_a,"times"_a);
}