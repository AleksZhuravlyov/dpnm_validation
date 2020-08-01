#include "math/Local.h"
#include <cmath>

Local::Local(Props &props) :
    props(props),
    dRadius(0),
    alpha(props.gridBlockN, 0),
    radiusCurr(props.gridBlockN + 1, 0),
    volCartes(props.gridBlockN, 0) {}

int Local::left(const int &index) {
    return index;
}

int Local::right(const int &index) {
    return index + 1;
}

double Local::calcDelRadius(const double &radius,
                            const double &effRadius) {

    return (effRadius - radius) / props.gridBlockN;
}

void Local::calcMatrCoordCurr(const double &radius,
                              const double &effRadius) {

    dRadius = calcDelRadius(radius, effRadius);

    for (int i = 0; i < props.gridBlockN + 1; i++)
        radiusCurr[i] = radius + i * dRadius;
}

void Local::calcVolCartesian() {

    for (int i = 0; i < alpha.size(); i++)
        volCartes[i] = props.lenY * props.lenZ *
                       (props.XCoordOut - props.XCoordIn) /
                       props.gridBlockN;
}

void Local::calculateAlpha(const double &dt,
                           const std::vector<double> &vol) {

    for (int i = 0; i < alpha.size(); i++)
        alpha[i] = vol[i] / dt;
}













