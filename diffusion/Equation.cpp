#include "Equation.h"

Equation::Equation(
    const std::map<std::string, std::variant<int, double>> &params) :
    props(params),
    local(props),
    convective(props, local),
    dim(props.gridBlockN),
    conc(std::vector<std::vector<double>>()),
    time(props.time),
    iCurr(0), iPrev(1),
    matrix(dim, dim),
    freeVector(dim),
    guessVector(dim),
    variable(dim) {

  std::vector<Triplet> triplets;
  triplets.reserve(3 * dim - 4);
  triplets.emplace_back(0, 0);
  for (int i = 1; i < dim - 1; i++) {
    triplets.emplace_back(i, i - 1);
    triplets.emplace_back(i, i);
    triplets.emplace_back(i, i + 1);
  }
  triplets.emplace_back(dim - 1, dim - 2);
  triplets.emplace_back(dim - 1, dim - 1);
  matrix.setFromTriplets(triplets.begin(), triplets.end());

  for (int i = 0; i < dim; i++) {
    freeVector[i] = 0;
    guessVector[i] = 0;
    variable[i] = 0;
  }

  local.calcVolCartesian();
  convective.calcOmegaCartes();

}


void Equation::calculate() {

  cfdProcedure("newman",
               local.volCartes,
               convective.omegaCartesian);

}


void Equation::calculateMatrix() {

  MatrixIterator(matrix, 0).valueRef() = local.alpha[0];

  for (int i = 1; i < dim - 1; ++i) {
    MatrixIterator it(matrix, i);
    double &betaLeft = convective.beta[Local::left(i)];
    double &betaRight = convective.beta[Local::right(i)];
    double &alpha = local.alpha[i];
    it.valueRef() = -1 * betaLeft;
    ++it;
    it.valueRef() = alpha + betaLeft + betaRight;
    ++it;
    it.valueRef() = -1 * betaRight;
  }

  MatrixIterator it(matrix, dim - 1);
  double &betaLeft = convective.beta[Local::left(dim - 1)];
  double &alpha = local.alpha[dim - 1];
  it.valueRef() = -1 * betaLeft;
  ++it;
  it.valueRef() = alpha + betaLeft;
}

void Equation::calcConcIni(const double &concIni) {

  conc.emplace_back(std::vector<double>(dim, concIni));
  conc.emplace_back(std::vector<double>(dim, concIni));
}


void Equation::forceDirichletBound(const double &conc_in) {

  MatrixIterator it(matrix, dim - 1);
  it.valueRef() = 0;
  ++it;
  it.valueRef() = local.alpha[dim - 1];

  freeVector[dim - 1] = local.alpha[dim - 1] * 5 * conc_in;

}

void Equation::calculateFreeVector(const double &conc_in) {
  freeVector[0] = local.alpha[0] * conc_in;
  for (int i = 1; i < dim; i++)
    freeVector[i] = local.alpha[i] * conc[iPrev][i];
}

void Equation::calculateGuessVector() {
  for (int i = 0; i < dim; i++)
    guessVector[i] = conc[iPrev][i];
}

void Equation::calculateConc() {

  BiCGSTAB biCGSTAB;

  biCGSTAB.compute(matrix);

  variable = biCGSTAB.solveWithGuess(freeVector, guessVector);

  for (int i = 0; i < dim; i++)
    conc[iCurr][i] = variable[i];
}

void Equation::calcVelocity() {
  velocity = 0;
  for (int i = 0; i < conc[iCurr].size(); i++)
    velocity -= conc[iCurr][i] - conc[iPrev][i];
  velocity /= props.density * props.lenY * props.lenZ;
}

void Equation::calcTimeVector() {

  auto time = props.time;
  auto configTimeStep = props.timeStep;
  double division = time / configTimeStep;
  double fullStepsN;
  auto lastStep = std::modf(division, &fullStepsN);
  auto _timeSteps = std::vector<double>(fullStepsN, configTimeStep);
  if (lastStep > 0)
    _timeSteps.push_back(lastStep * configTimeStep);
  timeSteps = _timeSteps;

}

void Equation::cfdProcedureOneStep(const std::string &boundCond,
                                   const double &concThrWall,
                                   const double &radius,
                                   const double &effRadius,
                                   const double &thrLength,
                                   const std::vector<double> &volumes,
                                   const std::vector<double> &surfaces,
                                   const double &dt) {

  std::swap(iCurr, iPrev);

  local.calculateAlpha(props.timeStep, volumes);

  convective.calculateBeta(radius,
                           effRadius,
                           thrLength,
                           props.diffusivity,
                           props.gridBlockN,
                           surfaces);

  calculateGuessVector();
  calculateMatrix();
  calculateFreeVector(concThrWall);
  if (boundCond == "dirichlet")
    forceDirichletBound(concThrWall);
  calculateConc();
  calcVelocity();

}

void Equation::cfdProcedure(const std::string &boundCond,
                            const std::vector<double> &volumes,
                            const std::vector<double> &surfaces) {

  calcTimeVector();
  calcConcIni(props.concInit);

  times.push_back(0);
  concs.push_back(conc[iCurr]);
  velocities.push_back(0);


  double timeCurr;
  for (int i = 0; i <= timeSteps.size(); i++) {

    timeCurr += timeSteps[i];

    cfdProcedureOneStep(boundCond, props.concLeft,
                        props.XCoordIn,
                        props.XCoordOut,
                        props.lenY,
                        volumes, surfaces,
                        timeSteps[i]);

    times.push_back(timeCurr);
    concs.push_back(conc[iCurr]);
    velocities.push_back(velocity);

  }
}


