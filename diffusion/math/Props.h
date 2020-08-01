#ifndef PROPSDIFFUSION_H
#define PROPSDIFFUSION_H

#include <iostream>
#include <map>
#include <variant>
#include <vector>

class Props {

public:

    explicit Props(
            const std::map<std::string, std::variant<int, double>> &params);

    virtual ~Props() {}

    std::map<std::string, std::variant<int, double>> _params;

    double time;
    double timeStep;
    double XCoordIn;
    double XCoordOut;
    double lenY;
    double lenZ;
    int gridBlockN;
    double concLeft;
    double concRight;
    double diffusivity;
    double iterativeAccuracy;

    void printParams();

private:

};

#endif


