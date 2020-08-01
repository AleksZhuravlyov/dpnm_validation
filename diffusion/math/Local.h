#ifndef PNFLOW_LOCALDIFFUSION_H
#define PNFLOW_LOCALDIFFUSION_H

#include "math/Props.h"

class Local {

public:

    explicit Local(Props &props);

    virtual ~Local() = default;

    static int left(const int &index);

    static int right(const int &index);

    double calcDelRadius(const double &radius, const double &effRadius);

    void calcVol(const std::string &coordType);

    void calcVolCartesian();

    void calcMatrCoordCurr(const double &radius, const double &effRadius);

    void calculateAlpha(const double &dt,
                        const std::vector<double> &vol);

    Props &props;

    double dRadius;
    std::vector<double> radiusCurr;
    std::vector<double> volCartes;
    std::vector<double> alpha;

};

#endif
