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
        omegaCartesian[i] = props.lenX * props.lenZ;
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


