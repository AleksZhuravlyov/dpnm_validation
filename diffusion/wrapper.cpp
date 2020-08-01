#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "math/Props.h"
#include "Equation.h"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(diffusion, m) {

  py::class_<Props>(m, "Props")
      .def(py::init<const std::map<std::string, std::variant<int, double>> &>(),
           "params"_a)

      .def_readwrite("params", &Props::_params)
      .def("print_params", &Props::printParams);

  py::class_<Equation>(m, "Equation")
      .def(py::init<const std::map<std::string, std::variant<int, double>> &>(),
           "params"_a)

      .def("calculate", &Equation::calculate)

      .def_readwrite("concs", &Equation::concs)
      .def_readwrite("velocities", &Equation::velocities)
      .def_readwrite("time_steps", &Equation::timeSteps)
      .def_readwrite("times", &Equation::times);

}