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


#ifndef EQUATION_H
#define EQUATION_H

//#include <pybind11/eigen.h>

#include <vector>
#include <Eigen/Sparse>

#include "math/Props.h"
#include "math/Local.h"
#include "math/Convective.h"


typedef Eigen::Triplet<double> Triplet;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> Matrix;
typedef Matrix::InnerIterator MatrixIterator;
typedef Eigen::VectorXd Vector;
typedef Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> BiCGSTAB;


class Equation {

public:

  explicit Equation(
      const std::map<std::string, std::variant<int, double>> &params);

  virtual ~Equation() = default;

  void calculate();

  void calculateMatrix();

  void calculateFreeVector(const double &conc_in);

  void calculateGuessVector();

  void calculateConc();

  void calcConcIni(const double &concIni);

  void calcTimeVector();

  void forceDirichletBound(const double &concIni);

  void cfdProcedureOneStep(const std::string &boundCond,
                           const double &concThrWall,
                           const double &radius,
                           const double &effRadius,
                           const double &thrLength,
                           const std::vector<double> &volumes,
                           const std::vector<double> &surfaces,
                           const double &dt);

  void cfdProcedure(const std::string &boundCond,
                    const std::vector<double> &volumes,
                    const std::vector<double> &surfaces);


  void calcVelocity();

  int &dim;

  std::vector<std::vector<double>> conc;
  std::vector<double> timeSteps;
  std::vector<double> times;
  std::vector<std::vector<double>> concs;

  Props props;
  Local local;
  Convective convective;

  double &time;
  int iCurr;
  int iPrev;

  double velocity;
  std::vector<double> velocities;

  Matrix matrix;

  Vector freeVector;
  Vector guessVector;
  Vector variable;

};


#endif
