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
        volCartes[i] = props.lenX * props.lenZ *
                       (props.YCoordOut - props.YCoordIn) /
                       props.gridBlockN;
}

void Local::calculateAlpha(const double &dt,
                           const std::vector<double> &vol) {

    for (int i = 0; i < alpha.size(); i++)
        alpha[i] = vol[i] / dt;
}













