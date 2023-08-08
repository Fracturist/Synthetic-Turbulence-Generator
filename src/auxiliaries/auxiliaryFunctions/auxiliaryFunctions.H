/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 AUTHOR,AFFILIATION
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

Description
    Some auxiliary functions for convenience.

SourceFiles
    auxiliaryFunctions.C

Author
    Hao Guo

\*---------------------------------------------------------------------------*/

#ifndef auxiliaryFunctions_H
#define auxiliaryFunctions_H

#include "fvCFD.H"
#include <functional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Random Functions

    //- Return an random unit vector (globally) generated
    //- based on standard Gauss/normal distribution (globalGaussNormal)
    vector randomUnitVector(Random& rndGen);


    //- Return an random vector (globally) with components subject to
    //- standard Gauss/normal distribution (globalGaussNormal) independently
    vector randomGaussVector(Random& rndGen);

    //- Return an random unit vector (locally) generated
    //- based on standard Gauss/normal distribution (GaussNormal)
    vector localUnitVector(Random& rndGen);


    //- Return an random vector (locally) with components subject to
    //- standard Gauss/normal distribution (GaussNormal) independently
    vector localGaussVector(Random& rndGen);


// Tensor Functions

    //- Return the result tensor of
    //- cholesky decomposition of a symmetric tensor
    autoPtr<tensor> Chol(const symmTensor& R);

    //- Return the result tensor field of
    //- cholesky decomposition of a symmetric tensor field
    tmp<tensorField> Chol(const symmTensorField& R);


// Non-linear Equation Iterative Solver

    //- Secant method
    scalar secantMethod
    (
        std::function<scalar(scalar)> const& fun,
        scalar x0 = 0.0,
        scalar x00 = -1.0,
        scalar tolerance = SMALL,
        label maxIter = 1000
    );

    //- Secant method with output message
    scalar secantMethodDebug
    (
        std::function<scalar(scalar)> const& fun,
        scalar x0 = 0.0,
        scalar x00 = -1.0,
        scalar tolerance = SMALL,
        label maxIter = 1000
    );

    //- Newton's method
    scalar NewtonMethod
    (
        std::function<scalar(scalar)> const& fun,
        std::function<scalar(scalar)> const& dfun,
        scalar x0 = 0.0,
        scalar tolerance = SMALL,
        label maxIter = 1000
    );

    //- Newton's method with output message
    scalar NewtonMethodDebug
    (
        std::function<scalar(scalar)> const& fun,
        std::function<scalar(scalar)> const& dfun,
        scalar x0 = 0.0,
        scalar tolerance = SMALL,
        label maxIter = 1000
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
