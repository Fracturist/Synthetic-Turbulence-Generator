/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    spectralInverter

Description
    Extended low-divergence spectral method (inverter version).

Author
    Hao Guo

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"
#include "singlePhaseTransportModel.H"
#include "spectrumModel.H"
#include "auxiliaryFunctions.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Generate initial field for DNS/LES/DES using "
        "extended spectral method (inverter version)."
    );

    // Additional options
    {
        argList::addOption
        (
            "turb",
            "name",
            "Specify the name of the turbulent field to be processed, \n"
            "e.g. turbulent kinetic energy k or Reynolds stress tensor R "
            "(default is k)"
        );

        argList::addBoolOption
        (
            "readEpsilon",
            "Use turbulent dissipation field epsilon to calculate integral length instead of vice versa"
        );

        argList::addBoolOption
        (
            "robust",
            "Use primitive version of spectrum-based method extended "
            "by merely multiplying Chol decomposion of R."
        );
    }

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl
        << "Generating velocity field using spectrum-based method"
        << nl << endl;

    #include "determineScale.H"

    #include "createFields.H"

    const bool robust = args.found("robust");

    if (robust)
    {
        Info<< "Running primitve spectral method\n" << endl;
    }

    const scalar kappaGrowth = nModes > 1 ?
        Foam::pow( kappaMax/kappaMin, 1./(nModes-1.) ) : 1.0;

    const scalar deltaKappaConst = logKappa ?
        (kappaGrowth - 1./kappaGrowth)/2.0 :
        (
            nModes > 1 ?
            (kappaMax - kappaMin)/(nModes-1.0) : 0.0
        );

    for (label modei = 1; modei <= nModes; ++modei)
    {
        // Calculate wave mode and mode distance
        scalar kappaM = logKappa ?
            kappaMin*Foam::pow(kappaGrowth, modei-1) :
            kappaMin + deltaKappaConst*(modei-1);

        scalar deltaKappa = logKappa ?
            kappaM*deltaKappaConst :
            deltaKappaConst;

        Info<< "Processing mode:" << modei << " kappaM:" << kappaM << endl;

        // Normalized energy
        tmp<scalarField> tE;
        if (useSpectrumTable)
        {
            tE =
                tmp<scalarField>::New
                (
                    C.size(),
                    EkTab->value(kappaM)/Foam::sqr(uTur0)
                );
        }
        else
        {
            tE = EkMod->value(kappaM);
        }

        // Wave amplitude
        scalarField qmNorm = Foam::sqrt(tE()*deltaKappa);

        // Phase angle
        scalar psim = rndGen.globalPosition(-mathematical::pi, mathematical::pi);

        // Direction vector sigma
        const vector sigma(randomUnitVector(rndGen));

        // Intermediate random vector zeta
        const vector zeta(randomUnitVector(rndGen));

        // Wave vector direction
        vectorField kappaDir;
        if (!robust)
        {
            const tensorField& L = tL().primitiveField();

            kappaDir = (L & sigma) ^ zeta;

            forAll(kappaDir, i)
            {
                const scalar flag(mag(kappaDir[i]));
                kappaDir[i] =
                    flag > ROOTSMALL ? kappaDir[i]/flag : vector::zero;
            }
        }
        else
        {
            vector ksi = sigma ^ zeta;
            ksi /= mag(ksi);
            kappaDir.resize(C.size(), ksi);
        }

        // Add the velocity contribution per mode
        Uc += 2*qmNorm*cos(kappaM*(kappaDir & C) + psim)*sigma;
    }

    Uc = tL().primitiveField() & Uc;
    U.primitiveFieldRef() += Uc;
    U.correctBoundaryConditions();
    U.write();

    tL.clear();

    #include "writeDerivedFields.H"


    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
