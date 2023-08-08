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

Class
    Foam::spectrumModel

Description
    Class for different spectrum models.

SourceFiles
    spectrumModel.C

Author
    Hao Guo

\*---------------------------------------------------------------------------*/

#ifndef spectrumModel_H
#define spectrumModel_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Istream;
class Ostream;
class spectrumModel;

Istream& operator>>(Istream&, spectrumModel&);
Ostream& operator<<(Ostream&, const spectrumModel&);
Ostream& operator<<(Ostream&, const autoPtr<spectrumModel>&);


/*---------------------------------------------------------------------------*\
                         Class spectrumModel Declaration
\*---------------------------------------------------------------------------*/

class spectrumModel
// :
//     public baseClassName
{
    // Private Data

        //- Peak energy kappa
        scalar kappa0_;

        //- Turbulent characteristic velocity
        scalar uTur0_;

        //- Field size of E, L and eta
        label size_;


        //- Inertial range model constant
        scalar C0_;

        //- Energe-containing range model constant
        scalar p0_;
        scalar cL_;

        //- Dissipation range range model constant
        scalar q0_;
        Switch unityModel_;
        scalar beta_;
        scalar cEta_;

        //- Low Re_t model switch
        Switch simpleModel_;

        //- kappa0*L constant
        scalar kappa0TimesL_;


        //- Integral length scale
        scalarField L_;

        //- Dissipation length scale
        scalarField eta_;


public:

    // Constructors

        //- Construct from components given kappa0, uTur0 and size
        spectrumModel
        (
            scalar kappa0 = 50,
            scalar uTur0 = 0.22,
            label size = 1,
            scalar C0 = 1.5,
            scalar p0 = 4,
            scalar q0 = 0.5,
            const Switch& unityModel = true,
            scalar beta = 5.2,
            const Switch& simpleModel = false
        );

        //- Construct from components given L and eta field
        spectrumModel
        (
            const scalarField& L,
            const scalarField& eta,
            scalar C0 = 1.5,
            scalar p0 = 4,
            scalar q0 = 0.5,
            const Switch& unityModel = true,
            scalar beta = 5.2,
            const Switch& simpleModel = false
        );

        //- Construct from components given TKE and nu field, and kappa0
        spectrumModel
        (
            scalar kappa0,
            const scalarField& k,
            const scalarField& nu,
            scalar C0 = 1.5,
            scalar p0 = 4,
            scalar q0 = 0.5,
            const Switch& unityModel = true,
            scalar beta = 5.2,
            const Switch& simpleModel = false
        );

        //- Construct from dictionary
        spectrumModel
        (
            const dictionary& dict,
            label size = 1
        );

        //- Construct from dictionary with L and eta field
        spectrumModel
        (
            const dictionary& dict,
            const scalarField& L,
            const scalarField& eta
        );

        //- Construct from dictionary with TKE and nu field
        spectrumModel
        (
            const scalarField& k,
            const scalarField& nu,
            const dictionary& dict
        );


    //- Destructor
    virtual ~spectrumModel();


    // Member Functions

        //- Return overall Ek value
        //- as a function of (scalar) independent variable
        virtual tmp<scalarField> value(scalar kappa) const;

        //- Return inertial range Ek value
        //- as a function of (scalar) independent variable
        virtual tmp<scalarField> inertialValue(scalar kappa) const;

        //- Return energy-containing range function value
        //- as a function of (scalar) independent variable
        virtual tmp<scalarField> energyContFunc(scalar kappa) const;

        //- Return dissipation range function value
        //- as a function of (scalar) independent variable
        virtual tmp<scalarField> dissipationFunc(scalar kappa) const;


    // Access

        //- Return constant C0
        inline scalar C0() const
        {
            return C0_;
        }

        //- Return power exponent p0 of energy-containing range function
        inline scalar p0() const
        {
            return p0_;
        }

        //- Return power exponent q0 of dissipation range function
        inline scalar q0() const
        {
            return q0_;
        }

        //- Return model constant c_L of energy-containing range function
        inline scalar cL() const
        {
            return cL_;
        }

        //- Return model constant c_eta of dissipation range function
        inline scalar cEta() const
        {
            return cEta_;
        }

        //- Return model constant beta of dissipation range function
        inline scalar beta() const
        {
            return beta_;
        }

        //- Return the constant product of extreme value kappa and L
        inline scalar kappa0TimesL() const
        {
            return kappa0TimesL_;
        }

        //- Use unity-consistent model or not
        inline const Switch& unityModel() const
        {
            return unityModel_;
        }

        //- Use simple low Re_t model or not
        inline const Switch& simpleModel() const
        {
            return simpleModel_;
        }

        //- Return size of the variable field
        inline label size() const
        {
            return size_;
        }

        //- Return integral length scale L
        inline const scalarField& L() const
        {
            return L_;
        }

        //- Return dissipation length scale eta
        inline const scalarField& eta() const
        {
            return eta_;
        }


    // Write

        //- Write in dictionary format
        void writeEntries(Ostream& os) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, spectrumModel&);
        friend Ostream& operator<<(Ostream&, const spectrumModel&);
        friend Ostream& operator<<(Ostream&, const autoPtr<spectrumModel>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
