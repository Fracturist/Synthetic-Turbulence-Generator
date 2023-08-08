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

#include "spectrumModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::spectrumModel::spectrumModel
(
    scalar kappa0,
    scalar uTur0,
    label size,
    scalar C0,
    scalar p0,
    scalar q0,
    const Switch& unityModel,
    scalar beta,
    const Switch& simpleModel
)
:
    kappa0_(kappa0),
    uTur0_(uTur0),
    size_(size),
    C0_(C0),
    p0_(p0),
    q0_(q0),
    unityModel_(unityModel),
    beta_(beta),
    simpleModel_(simpleModel)
{
    if (!simpleModel_)
    {
        // General model

        if (!unityModel_)
        {
            beta_ = Foam::pow(2*C0_*q0_*Foam::tgamma(4.0/3.0*q0_), 3.0/4.0/q0_);
        }
        else
        {
            beta_ = 0.52;
            C0_ = 1.5;
            cEta_ = 0.401684789281759;
        }

        cL_ = Foam::pow3
        (
            Foam::tgamma(1.0/3.0)
        *Foam::tgamma(0.5 + p0_/2)
        /Foam::tgamma(5.0/6.0 + p0_/2)
        *C0_/3.0
        );
        kappa0TimesL_ = Foam::sqrt(0.6*p0_*cL_);

        scalar L0 = kappa0TimesL_/kappa0_;
        scalar epsilon0 = Foam::pow3(uTur0_)/L0;
        scalar eta0 = Foam::pow(epsilon0, -0.25)*Foam::pow(1e-5, 0.75);
        L_ = scalarField(size_, L0);
        eta_ = scalarField(size_, eta0);
    }
    else
    {
        // Simple (low Re_t) model

        C0_ = 16.0*Foam::sqrt(2.0/M_PI);

        beta_ = 2.0;
        q0_ = 0.5;
        p0_ = 4.0;
        kappa0TimesL_ = 1.0;

        L_ = scalarField(size_, kappa0TimesL_/kappa0);
        eta_ = L_;

        cL_ = 0;
        cEta_ = 0;
    }
}


Foam::spectrumModel::spectrumModel
(
    const scalarField& L,
    const scalarField& eta,
    scalar C0,
    scalar p0,
    scalar q0,
    const Switch& unityModel,
    scalar beta,
    const Switch& simpleModel
)
:
    C0_(C0),
    p0_(p0),
    q0_(q0),
    unityModel_(unityModel),
    beta_(beta),
    simpleModel_(simpleModel)
{
    if (L.size() != eta.size())
    {
        FatalErrorInFunction
            << "Size of L and eta don't match"
            << abort(FatalError);
    }
    size_ = L.size();

    if (!simpleModel_)
    {
        // General model

        L_ = L;
        eta_ = eta;

        if (!unityModel_)
        {
            beta_ = Foam::pow(2*C0_*q0_*Foam::tgamma(4.0/3.0*q0_), 3.0/4.0/q0_);
        }
        else
        {
            beta_ = 0.52;
            C0_ = 1.5;
            cEta_ = 0.401684789281759;
        }

        cL_ = Foam::pow3
        (
            Foam::tgamma(1.0/3.0)
        *Foam::tgamma(0.5 + p0_/2)
        /Foam::tgamma(5.0/6.0 + p0_/2)
        *C0_/3.0
        );
        kappa0TimesL_ = Foam::sqrt(0.6*p0_*cL_);
    }
    else
    {
        // Simple (low Re_t) model
        L_ = L;
        eta_ = L_;

        C0_ = 16.0*Foam::sqrt(2.0/M_PI);

        beta_ = 2.0;
        q0_ = 0.5;
        p0_ = 4.0;
        kappa0TimesL_ = 1.0;

        cL_ = 0;
        cEta_ = 0;
    }
}


Foam::spectrumModel::spectrumModel
(
    scalar kappa0,
    const scalarField& k,
    const scalarField& nu,
    scalar C0,
    scalar p0,
    scalar q0,
    const Switch& unityModel,
    scalar beta,
    const Switch& simpleModel
)
:
    kappa0_(kappa0),
    C0_(C0),
    p0_(p0),
    q0_(q0),
    unityModel_(unityModel),
    beta_(beta),
    simpleModel_(simpleModel)
{
    if (k.size() != nu.size())
    {
        FatalErrorInFunction
            << "Size of k and nu don't match"
            << abort(FatalError);
    }
    size_ = k.size();

    if (!simpleModel_)
    {
        // General model

        if (!unityModel_)
        {
            beta_ = Foam::pow(2*C0_*q0_*Foam::tgamma(4.0/3.0*q0_), 3.0/4.0/q0_);
        }
        else
        {
            beta_ = 0.52;
            C0_ = 1.5;
            cEta_ = 0.401684789281759;
        }

        cL_ = Foam::pow3
        (
            Foam::tgamma(1.0/3.0)
        *Foam::tgamma(0.5 + p0_/2)
        /Foam::tgamma(5.0/6.0 + p0_/2)
        *C0_/3.0
        );
        kappa0TimesL_ = Foam::sqrt(0.6*p0_*cL_);

        scalar L0 = kappa0TimesL_/kappa0_;
        scalarField epsilon = Foam::pow( 2.0/3.0*k, 1.5 )/L0;
        L_ = scalarField(size_, L0);
        eta_ = Foam::pow(epsilon, -0.25)*Foam::pow(nu, 0.75);
    }
    else
    {
        // Simple (low Re_t) model

        C0_ = 16.0*Foam::sqrt(2.0/M_PI);

        beta_ = 2.0;
        q0_ = 0.5;
        p0_ = 4.0;
        kappa0TimesL_ = 1.0;

        L_ = scalarField(size_, kappa0TimesL_/kappa0);
        eta_ = L_;

        cL_ = 0;
        cEta_ = 0;
    }
}


Foam::spectrumModel::spectrumModel
(
    const dictionary& dict,
    label size
)
:
    spectrumModel
    (
        dict.subDict("modeledSpectrum").getOrDefault<scalar>("kappa0", 50),
        dict.subDict("modeledSpectrum").getOrDefault<scalar>("uTur0", 0.22),
        size,
        dict.subDict("modeledSpectrum").getOrDefault<scalar>("C", 1.5),
        dict.subDict("modeledSpectrum").getOrDefault<scalar>("p0", 4.0),
        dict.subDict("modeledSpectrum").getOrDefault<scalar>("q0", 0.5),
        dict.subDict("modeledSpectrum").getOrDefault<Switch>("unityModel", true),
        dict.subDict("modeledSpectrum").getOrDefault<scalar>("beta", 5.2),
        dict.subDict("modeledSpectrum").getOrDefault<Switch>("simpleModel", false)
    )
{}


Foam::spectrumModel::spectrumModel
(
    const dictionary& dict,
    const scalarField& L,
    const scalarField& eta
)
:
    spectrumModel
    (
        L,
        eta,
        dict.subDict("modeledSpectrum").getOrDefault<scalar>("C", 1.5),
        dict.subDict("modeledSpectrum").getOrDefault<scalar>("p0", 4.0),
        dict.subDict("modeledSpectrum").getOrDefault<scalar>("q0", 0.5),
        dict.subDict("modeledSpectrum").getOrDefault<Switch>("unityModel", true),
        dict.subDict("modeledSpectrum").getOrDefault<scalar>("beta", 5.2),
        dict.subDict("modeledSpectrum").getOrDefault<Switch>("simpleModel", false)
    )
{}


Foam::spectrumModel::spectrumModel
(
    const scalarField& k,
    const scalarField& nv,
    const dictionary& dict
)
:
    spectrumModel
    (
        dict.subDict("modeledSpectrum").getOrDefault<scalar>("kappa0", 50),
        k,
        nv,
        dict.subDict("modeledSpectrum").getOrDefault<scalar>("C", 1.5),
        dict.subDict("modeledSpectrum").getOrDefault<scalar>("p0", 4.0),
        dict.subDict("modeledSpectrum").getOrDefault<scalar>("q0", 0.5),
        dict.subDict("modeledSpectrum").getOrDefault<Switch>("unityModel", true),
        dict.subDict("modeledSpectrum").getOrDefault<scalar>("beta", 5.2),
        dict.subDict("modeledSpectrum").getOrDefault<Switch>("simpleModel", false)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::spectrumModel::~spectrumModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<scalarField> Foam::spectrumModel::value(scalar kappa) const
{
    if (simpleModel_)
    {
        return tmp<scalarField>::New
        (
            C0_
           *Foam::pow4(kappa)
           *Foam::pow5(L_)
           *Foam::exp
            (
                -2.0
               *Foam::sqr(kappa*L_)
            )
        );
    }

    return tmp<scalarField>::New
    (
        inertialValue(kappa)*energyContFunc(kappa)*dissipationFunc(kappa)
    );
}


Foam::tmp<scalarField> Foam::spectrumModel::inertialValue(scalar kappa) const
{
    if (simpleModel_)
    {
        NotImplemented;
    }

    return tmp<scalarField>::New
    (
        C0_*Foam::pow(L_, -2.0/3)*Foam::pow(kappa, -5.0/3)
    );
}


Foam::tmp<scalarField> Foam::spectrumModel::energyContFunc(scalar kappa) const
{
    if (simpleModel_)
    {
        NotImplemented;
    }

    return tmp<scalarField>::New
    (
        Foam::pow
        (
            kappa*L_
           /Foam::sqrt(Foam::sqr(kappa*L_) + cL_),
            5.0/3 + p0_
        )
    );
}


Foam::tmp<scalarField> Foam::spectrumModel::dissipationFunc(scalar kappa) const
{
    if (simpleModel_)
    {
        NotImplemented;
    }

    if (!unityModel_)
    {
        return tmp<scalarField>::New
        (
            Foam::exp
            (
                -beta_
               *Foam::pow
                (
                    kappa*eta_,
                    1.0/q0_
                )
            )
        );
    }
    else
    {
        return tmp<scalarField>::New
        (
            Foam::exp
            (
                -beta_
               *(
                    Foam::pow
                    (
                        Foam::pow4(kappa*eta_) + Foam::pow4(cEta_),
                        0.25
                    )
                  - cEta_
                )
            )
        );
    }
}


// * * * * * * * * * * * * * * Istream Operator  * * * * * * * * * * * * * * //

Istream& Foam::operator>>(Istream& is, spectrumModel&)
{
    NotImplemented;
    return is;
}


// * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * * //

void Foam::spectrumModel::writeEntries(Ostream& os) const
{
    if (!simpleModel_)
    {
        // General model

        os.writeEntry("C0", C0_);
        os.writeEntry("p0", p0_);
        os.writeEntry("q0", q0_);

        os.writeEntry("cL", cL_);
        if (!unityModel_)
        {
            os.writeEntry("cEta", cEta_);
        }
        os.writeEntry("beta", beta_);
        os.writeEntry("k0TimesL", kappa0TimesL_);

        os.writeEntry("size", size_);
    }
    else
    {
        // Simple (low Re_t) model

        os.writeEntry("alpha", C0_);
        os.writeEntry("monomial kappa index", p0_);
        os.writeEntry("exp kappa index", 1./q0_);

        os.writeEntry("beta", beta_);
        os.writeEntry("k0TimesL", kappa0TimesL_);

        os.writeEntry("size", size_);
    }
}


Ostream& Foam::operator<<(Ostream& os, const spectrumModel& Ek)
{
    os.check(FUNCTION_NAME);

    os.beginBlock("Ek");
    Ek.writeEntries(os);
    os.endBlock();

    return os;
}


Ostream& Foam::operator<<(Ostream& os, const autoPtr<spectrumModel>& Ek)
{
    os.check(FUNCTION_NAME);

    os.beginBlock("Ek");
    Ek->writeEntries(os);
    os.endBlock();

    return os;
}


// ************************************************************************* //
