#include "Props.h"
#include <vector>


Props::Props(
    const std::map<std::string, std::variant<int, double>> &params)

    : _params(params),
      timeEquil(std::get<double>(_params["time_equil"])),
      time(std::get<double>(_params["time"])),
      timeStep(std::get<double>(_params["time_step"])),
      YCoordIn(std::get<double>(_params["y_coord_in"])),
      YCoordOut(std::get<double>(_params["y_coord_out"])),
      lenX(std::get<double>(_params["len_x"])),
      lenZ(std::get<double>(_params["len_z"])),
      gridBlockN(std::get<int>(_params["grid_block_d_n"])),
      concLeft(std::get<double>(_params["conc_left"])),
      concInit(std::get<double>(_params["conc_init"])),
      diffusivity(std::get<double>(_params["diffusivity"])),
      density(std::get<double>(_params["density"])),
      iterativeAccuracy(std::get<double>(_params["it_accuracy"])) {}


void Props::printParams() {
  for (auto &ent : _params) {
    std::cout << ent.first << ": ";
    if (std::get_if<int>(&ent.second))
      std::cout << std::get<int>(ent.second) << std::endl;
    else if (std::get_if<double>(&ent.second))
      std::cout << std::get<double>(ent.second) << std::endl;
  }
}

