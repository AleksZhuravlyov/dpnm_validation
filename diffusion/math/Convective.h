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
