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

\*---------------------------------------------------------------------------*/

#include "auxiliaryFunctions.H"

// * * * * * * * * * * * * * * Random Functions  * * * * * * * * * * * * * * //

Foam::vector Foam::randomUnitVector(Random& rndGen)
{
    // Based on multivariate Gauss/Normal distribution
    scalar x = rndGen.globalGaussNormal<scalar>();
    scalar y = rndGen.globalGaussNormal<scalar>();
    scalar z = rndGen.globalGaussNormal<scalar>();
    scalar a = Foam::sqrt(x*x + y*y + z*z);

    return vector
    (
        x/a,
        y/a,
        z/a
    );
}


Foam::vector Foam::randomGaussVector(Random& rndGen)
{
    // Based on multivariate Gauss/Normal distribution
    scalar x = rndGen.globalGaussNormal<scalar>();
    scalar y = rndGen.globalGaussNormal<scalar>();
    scalar z = rndGen.globalGaussNormal<scalar>();

    return vector
    (
        x,
        y,
        z
    );
}


Foam::vector Foam::localUnitVector(Random& rndGen)
{
    // Based on multivariate Gauss/Normal distribution
    scalar x = rndGen.GaussNormal<scalar>();
    scalar y = rndGen.GaussNormal<scalar>();
    scalar z = rndGen.GaussNormal<scalar>();
    scalar a = Foam::sqrt(x*x + y*y + z*z);

    return vector
    (
        x/a,
        y/a,
        z/a
    );
}


Foam::vector Foam::localGaussVector(Random& rndGen)
{
    // Based on multivariate Gauss/Normal distribution
    scalar x = rndGen.GaussNormal<scalar>();
    scalar y = rndGen.GaussNormal<scalar>();
    scalar z = rndGen.GaussNormal<scalar>();

    return vector
    (
        x,
        y,
        z
    );
}


// * * * * * * * * * * * * * * Tensor Functions  * * * * * * * * * * * * * * //

Foam::autoPtr<tensor> Foam::Chol(const symmTensor& R)
{
    autoPtr<tensor> tL = autoPtr<tensor>::New(Zero);
    tensor& L = tL.ref();

    // If R is effectively diagonal tensor
    // if ((sqr(R.xy()) + sqr(R.xz()) + sqr(R.yz())) < SMALL)
    if ( sqrt(sqr(R.xy()) + sqr(R.xz()) + sqr(R.yz())) < SMALL )
    {
        // If R is effectively zero tensor
        if ( sqrt(sqr(R.xx()) + sqr(R.yy()) + sqr(R.zz())) < SMALL )
        {
            return tL;
        }

        L.xx() = pos(R.xx()) ? Foam::sqrt(R.xx()) : 0.0;
        L.yy() = pos(R.yy()) ? Foam::sqrt(R.yy()) : 0.0;
        L.zz() = pos(R.zz()) ? Foam::sqrt(R.zz()) : 0.0;
        return tL;
    }

    L.xx() = Foam::sqrt(R.xx());
    L.yx() = R.yx()/(L.xx() + ROOTSMALL);
    L.yy() = Foam::sqrt( R.yy() - Foam::sqr(L.yx()) );
    L.zx() = R.zx()/(L.xx() + ROOTSMALL);
    L.zy() = ( R.zy() - L.zx()*L.yx() )/(L.yy() + ROOTSMALL);
    L.zz() = Foam::sqrt( R.zz() - Foam::sqr(L.zx()) - Foam::sqr(L.zy()) );

    return tL;
}


Foam::tmp<tensorField> Foam::Chol(const symmTensorField& R)
{
    tmp<tensorField> tL = tmp<tensorField>::New(R.size(), Zero);
    tensorField& L = tL.ref();

    forAll(R, i)
    {
        L[i] = Chol(R[i]);
    }

    return tL;
}

// * * * * * * * * * Non-linear Equation Iterative Solver  * * * * * * * * * //

Foam::scalar Foam::secantMethod
(
    std::function<scalar(scalar)> const& fun,
    scalar x0,
    scalar x00,
    scalar tolerance,
    label maxIter
)
{
    scalar x = 0;

    for (label i = 0; i < maxIter; ++i)
    {
        x = x0 - (x0 - x00)/(fun(x0) - fun(x00))*fun(x0);

        if (std::abs(fun(x)) < tolerance)
        {
            return x;
        }

        x00 = x0;
        x0 = x;
    }

    FatalErrorInFunction
        << "Result is not converged when reaching max " << maxIter
        << " iterations. Final residual is " << std::abs(fun(x))
        << abort(FatalError);

    return -1;
}


Foam::scalar Foam::secantMethodDebug
(
    std::function<scalar(scalar)> const& fun,
    scalar x0,
    scalar x00,
    scalar tolerance,
    label maxIter
)
{
    scalar x = 0;

    Info<< nl;
    Info<< "Tol = " << tolerance << nl
        << "Max iter = " << maxIter << nl
        << "x00 = " << x00 << nl
        << "x0 = " << x0 << nl << endl;
    for (label i = 0; i < maxIter; ++i)
    {
        x = x0 - (x0 - x00)/(fun(x0) - fun(x00))*fun(x0);

        Info<< "Iter " << i+1 << ": x = " << x
            << ", fun(x) = " << fun(x) << endl;

        if (std::abs(fun(x)) < tolerance)
        {
            Info << endl;
            return x;
        }

        x00 = x0;
        x0 = x;
    }

    FatalErrorInFunction
        << "Result is not converged when reaching max " << maxIter
        << " iterations. Final residual is " << std::abs(fun(x))
        << abort(FatalError);

    return -1;
}


Foam::scalar Foam::NewtonMethod
(
    std::function<scalar(scalar)> const& fun,
    std::function<scalar(scalar)> const& dfun,
    scalar x0,
    scalar tolerance,
    label maxIter
)
{
    scalar x = 0;

    for (label i = 0; i < maxIter; ++i)
    {
        x = x0 - fun(x0)/dfun(x0);

        if (std::abs(fun(x)) < tolerance)
        {
            return x;
        }

        x0 = x;
    }

    FatalErrorInFunction
        << "Result is not converged when reaching max " << maxIter
        << " iterations. Final residual is " << std::abs(fun(x))
        << abort(FatalError);

    return -1;
}


Foam::scalar Foam::NewtonMethodDebug
(
    std::function<scalar(scalar)> const& fun,
    std::function<scalar(scalar)> const& dfun,
    scalar x0,
    scalar tolerance,
    label maxIter
)
{
    scalar x = 0;

    Info<< nl;
    Info<< "Tol = " << tolerance << nl
        << "Max iter = " << maxIter << nl
        << "x0 = " << x0 << nl << endl;
    for (label i = 0; i < maxIter; ++i)
    {
        x = x0 - fun(x0)/dfun(x0);

        Info<< "Iter " << i+1 << ": x = " << x
            << ", fun(x) = " << fun(x) << endl;

        if (std::abs(fun(x)) < tolerance)
        {
            Info << endl;
            return x;
        }

        x0 = x;
    }

    FatalErrorInFunction
        << "Result is not converged when reaching max " << maxIter
        << " iterations. Final residual is " << std::abs(fun(x))
        << abort(FatalError);

    return -1;
}


// ************************************************************************* //
