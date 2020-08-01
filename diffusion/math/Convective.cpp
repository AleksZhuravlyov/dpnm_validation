#include "Convective.h"
#include <cmath>

Convective::Convective(Props &props,
                       Local &local) :
    props(props),
    local(local),
    beta(props.gridBlockN + 1, 0),
    omegaCartesian(props.gridBlockN + 1, 0) {}

void Convective::calcOmegaCartes() {

    for (int i = 0; i < beta.size(); i++)
        omegaCartesian[i] = props.lenZ * props.lenY;
}

std::vector<double>
Convective::calc_diffusivityList(const int &gridBlockN,
                                 const double &diffusivity) {

    return std::vector<double>(gridBlockN, diffusivity);
}

void Convective::calculateBeta(const double &radius,
                               const double &effRadius,
                               const double &length,
                               const double &diffusivity,
                               const int &gridBlockN,
                               const std::vector<double> &omega) {

    auto dRadius = local.calcDelRadius(radius, effRadius);
    auto diffusivityList = calc_diffusivityList(beta.size(),
                                                diffusivity);

    for (int i = 0; i < beta.size(); i++)
        beta[i] = diffusivityList[i] * omega[i] / dRadius;
}


