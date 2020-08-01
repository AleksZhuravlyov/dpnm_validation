#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <Sgrid.h>

namespace py = pybind11;
using namespace pybind11::literals;


PYBIND11_MODULE(sgrid, m) {
  py::class_<Sgrid>(m, "Sgrid")
      .def(py::init<Eigen::Ref<Eigen::Vector3i>,
               Eigen::Ref<Eigen::Vector3d>,
               Eigen::Ref<Eigen::Vector3d>>(),
           "points_dims"_a,
           "points_origin"_a,
           "spacing"_a)
      .def(py::init<std::string>(), "file_name"_a)

      .def("save", &Sgrid::save, "file_name"_a)


      .def_property("points_dims",
                    &Sgrid::getPointsDims, &Sgrid::setPointsDims)
      .def_property("spacing",
                    &Sgrid::getSpacing, &Sgrid::setSpacing)
      .def_property("points_origin",
                    &Sgrid::getPointsOrigin, &Sgrid::setPointsOrigin)

      .def_readwrite("points_N", &Sgrid::_pointsN)

      .def_property("cells_dims",
                    &Sgrid::getCellsDims, &Sgrid::setCellsDims)
      .def_readwrite("cells_N", &Sgrid::_cellsN)

      .def_property("faces_dims",
                    &Sgrid::getFacesDims, &Sgrid::setFacesDims)
      .def_readwrite("faces_N", &Sgrid::_facesN)

      .def_readwrite("cell_V", &Sgrid::_cellV)
      .def_property("face_S",
                    &Sgrid::getFaceS, &Sgrid::setFaceS)

      .def_property("neighbors_faces",
                    &Sgrid::getNeighborsFaces, &Sgrid::setNeighborsFaces)
      .def_property("neighbors_cells",
                    &Sgrid::getNeighborsCells, &Sgrid::setNeighborsCells)

      .def_property("normals_neighbors_cells",
                    &Sgrid::getNormalsNeighborsCells,
                    &Sgrid::setNormalsNeighborsCells)
      .def_property("normals_neighbors_faces",
                    &Sgrid::getNormalsNeighborsFaces,
                    &Sgrid::setNormalsNeighborsFaces)

      .def_property("points_arrays",
                    &Sgrid::getPointsArrays, &Sgrid::setPointsArrays)
      .def_property("cells_arrays",
                    &Sgrid::getCellsArrays, &Sgrid::setCellsArrays)
      .def_property("faces_arrays",
                    &Sgrid::getFacesArrays, &Sgrid::setFacesArrays);

  m.def("save_files_collection_to_file", &saveFilesCollectionToFile,
        "file_name"_a,"files_names"_a,"files_descriptions"_a);
}

