#include "Props.h"
#include <vector>


Props::Props(
        const std::map<std::string, std::variant<int, double>> &params)

        : _params(params),
          time(std::get<double>(_params["time"])),
          timeStep(std::get<double>(_params["time_step"])),
          XCoordIn(std::get<double>(_params["x_coord_in"])),
          XCoordOut(std::get<double>(_params["x_coord_out"])),
          lenY(std::get<double>(_params["len_y"])),
          lenZ(std::get<double>(_params["len_z"])),
          gridBlockN(std::get<int>(_params["grid_block_n"])),
          concLeft(std::get<double>(_params["conc_left"])),
          concRight(std::get<double>(_params["conc_right"])),
          diffusivity(std::get<double>(_params["diffusivity"])),
          iterativeAccuracy(std::get<double>(_params["it_accuracy"]))
          {}


void Props::printParams() {
    for (auto &ent : _params) {
        std::cout << ent.first << ": ";
        if (std::get_if<int>(&ent.second))
            std::cout << std::get<int>(ent.second) << std::endl;
        else if (std::get_if<double>(&ent.second))
            std::cout << std::get<double>(ent.second) << std::endl;
    }
}

