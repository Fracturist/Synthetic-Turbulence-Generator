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
    whiteNoiseMaker

Description
    Create white noise disturbation added to mean velocity field.

Author
    Hao Guo

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "auxiliaryFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Create white noise disturbation added to mean velocity field"
        " for startup of DNS/LES/DES."
    );

    // Additional options
    {
        argList::addOption
        (
            "turb",
            "name",
            "Specify the name of the turbulent field to be processed,\ne.g. turbulent kinetic energy k or Reynolds stress tensor R"
        );
    }

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "Generating velocity field using white noise\n" << endl;

    // Random number engine: different seed to make randon number irrelevant
    // on different processors
    Random rndGen(1380649 + Pstream::myProcNo());

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    tmp<volScalarField> tk;
    const word tField = args.getOrDefault<word>("turb", "k");
    if (tField == "R")
    {
        Info<< "Reading Reynolds stress tensor field R," << nl
            << "generating turbulent kinetic energy field k" << nl;

        volSymmTensorField R
        (
            IOobject
            (
                "R",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        tk = 0.5*tr(R);
    }
    else
    {
        Info<< "Reading turbulent kinetic energy field k" << nl;

        tk =
            tmp<volScalarField>::New
            (
                IOobject
                (
                    "k",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );
    }
    Info<< endl;

    vectorField& Uc = U.primitiveFieldRef();

    forAll(Uc, celli)
    {
        Uc[celli] += Foam::sqrt( tk()[celli]/1.5 )*localGaussVector(rndGen);
    }

    tk.clear();

    U.correctBoundaryConditions();
    U.write();

    #include "writeDerivedFields.H"


    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
