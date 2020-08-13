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


#include "convection.h"
#include <string>
#include <iostream>
#include <vector>
#include <map>

#include "fvCFD.H"
#include "pisoControl.H"

std::map<std::string, std::vector<double>>
calculate(const std::vector<double> &velocities,
          const std::vector<double> &times,
          const int &equilStepsN) {


  std::map<std::string, std::vector<double>> results;

  int argc = 3;

  char arg0[] = "program";
  char arg1[] = "-case";
  char arg2[] = "./inOut/convection";

  char **argv = new char *[3];
  argv[0] = arg0;
  argv[1] = arg1;
  argv[2] = arg2;

#include "setRootCaseLists.H"
#include "createTime.H"
#include "createMesh.H"

  pisoControl piso(mesh);

#include "createFields.h"
#include "initContinuityErrs.H"

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


  Info << "\nStarting time loop\n" << endl;

  const label &topId = mesh.boundaryMesh().findPatchID("top");
  const label &outletId = mesh.boundaryMesh().findPatchID("outlet");


  auto UTopSize = U.boundaryFieldRef()[topId].size();
  auto UOutletSize = U.boundaryFieldRef()[outletId].size();
  int timeInd = 0;

  double UCur = 0;


  while (runTime.loop()) {
    Info << "Time = " << runTime.timeName() << nl << endl;


    if (timeInd >= equilStepsN)
      for (int i = 0; i < UTopSize; i++)
        U.boundaryFieldRef()[topId][i][1] = -velocities[timeInd + 1 - equilStepsN];


#include "CourantNo.H"

    // Momentum predictor

    fvVectorMatrix UEqn
        (
            fvm::ddt(U)
            + fvm::div(phi, U)
            - fvm::laplacian(nu, U)
        );

    if (piso.momentumPredictor()) {
      solve(UEqn == -fvc::grad(p));
    }

    // --- PISO loop
    while (piso.correct()) {
      volScalarField rAU(1.0 / UEqn.A());
      volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p));
      surfaceScalarField phiHbyA
          (
              "phiHbyA",
              fvc::flux(HbyA)
              + fvc::interpolate(rAU) * fvc::ddtCorr(U, phi)
          );

      adjustPhi(phiHbyA, U, p);

      // Update the pressure BCs to ensure flux consistency
      constrainPressure(p, U, phiHbyA, rAU);

      // Non-orthogonal pressure corrector loop
      while (piso.correctNonOrthogonal()) {
        // Pressure corrector

        fvScalarMatrix pEqn
            (
                fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
            );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (piso.finalNonOrthogonalIter()) {
          phi = phiHbyA - pEqn.flux();
        }
      }

#include "continuityErrs.H"

      U = HbyA - rAU * fvc::grad(p);
      U.correctBoundaryConditions();
    }


    runTime.write();

    /*auto &vX = p.internalField();

    const label &ID = mesh.boundaryMesh().findPatchID("velocity-inlet-5");
    auto Ub = U.boundaryFieldRef()[ID];

    std::cout << Ub[0][0] << std::endl;
    std::cout << Ub[0][1] << std::endl;
    std::cout << Ub[0][2] << std::endl;
    std::cout << std::endl;


    std::cout << runTime.value() << std::endl;

    if (int(runTime.value()) == 10) {
      std::cout << "bang" << std::endl;
      for (int i = 0; i < U.boundaryFieldRef()[ID].size(); i++) {
        U.boundaryFieldRef()[ID][i][0] = 2;
      }
    }


    if (int(runTime.value()) == 13) {

      const label &out = mesh.boundaryMesh().findPatchID("pressure-outlet-7");
      auto Uout = U.boundaryFieldRef()[out];

      std::cout << "Uout.size() " << Uout.size() << std::endl;
      for (int i = 0; i < Uout.size(); i++)
        std::cout << Uout[i][1] << std::endl;

    }*/

    Info << "ᕦ(ò_óˇ)ᕤ ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime.elapsedClockTime() << " s"
         << nl << endl;

    if (timeInd >= (equilStepsN - 1)) {

      UCur = 0;
      for (int i = 0; i < UTopSize; i++)
        UCur -= U.boundaryFieldRef()[topId][i][1];
      results["U_release_av"].push_back(UCur / UTopSize);

      UCur = 0;
      for (int i = 0; i < UOutletSize; i++)
        UCur += U.boundaryFieldRef()[outletId][i][0];
      results["U_outlet_av"].push_back(UCur / UOutletSize);
    } else {
      UCur = 0;
      for (int i = 0; i < UOutletSize; i++)
        UCur += U.boundaryFieldRef()[outletId][i][0];
      results["U_outlet_equil_av"].push_back(UCur / UOutletSize);
    }

    timeInd++;

  }

  Info << "End\n" << endl;

  return results;
}
