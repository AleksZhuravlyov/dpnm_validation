#ifndef CONVECTION_CONVECTION_H
#define CONVECTION_CONVECTION_H

#include <vector>
#include <map>


std::map<std::string, std::vector<double>> calculate(const std::vector<double> &velocities,
               const std::vector<double> &times,
               const int &equilStepsN);

#endif // CONVECTION_CONVECTION_H
