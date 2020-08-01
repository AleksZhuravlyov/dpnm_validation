#ifndef MATH_CONVECTIVE_H
#define MATH_CONVECTIVE_H

#include "math/Props.h"
#include "math/Local.h"

class Convective {

public:

    explicit Convective(Props &props,
                        Local &local);

    virtual ~Convective() = default;

    void calcOmegaCartes();

    void calculateBeta(const double &radius,
                       const double &effRadius,
                       const double &length,
                       const double &diffusivity,
                       const int &gridBlockN,
                       const std::vector<double> &omega);

    std::vector<double> calc_diffusivityList(const int &gridBlockN,
                                             const double &diffusivity);

    Props &props;
    Local &local;
    std::vector<double> omegaCartesian;
    std::vector<double> beta;


};

#endif
