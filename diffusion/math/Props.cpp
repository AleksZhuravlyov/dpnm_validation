/* MIT License
 *
 * Copyright (c) 2020 Aleksandr Zhuravlyov and Zakhar Lanets
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


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

