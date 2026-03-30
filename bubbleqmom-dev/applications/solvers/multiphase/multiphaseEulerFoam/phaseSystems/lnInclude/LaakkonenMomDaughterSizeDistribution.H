/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::diameterModels::daughterSizeDistributionMomModels::
    LaakkonenMomDaughterSizeDistribution

Description
    Daughter size distribution model used by Laakkonen et al. (2006). Note that
    it must be multiplied by the number of daughter particles (2). Also, the
    coefficient C4 needs to be set to the value of 2 as to ensure number and
    mass conservation following Kumar and Ramkrishna (1996). Also the internal
    coordinate was changed to particle volume, which gives

    \f[
		\int_0^{d_i} d^k \beta(d|d_i) dd = \frac{3240 d_i^k}{(k+9)(k+12)(k+15)}
    \f]

    where

    \vartable
        v_i         |  Volume of daughter bubble i [m]
        v_j         |  Volume of mother bubble j [m]
    \endvartable

    References:
    \verbatim
        "Validation of bubble breakage, coalescence and mass transfer models for
         gas-liquid dispersion in agitated vessel"
        Laakkonen, M., Alopaeus, V., Aittamaa, J.
        Chemical Engineering Science, Vol. 61, 2006, pp. 218-228
        Eq. 3, p. 220
    \endverbatim

    \verbatim
        "On the solution of population balance equations by
         discretization - I. A fixed pivot technique."
        Kumar, S. & Ramkrishna, D.
        Chemical Engineering Science, Vol. 51(8) 1996, 1311-1332.
        Eq. 33, p. 1318
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        C1           | Coefficient C1          | no          | 6.0
        C2           | Coefficient C2          | no          | 0.04
        C3           | Coefficient C3          | no          | 0.01
    \endtable

SourceFiles
    LaakkonenAlopaeusAittamaaBinaryBreakup.C

\*---------------------------------------------------------------------------*/

#ifndef LaakkonenMomDaughterSizeDistribution_H
#define LaakkonenMomDaughterSizeDistribution_H

#include "daughterSizeDistributionMomModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace daughterSizeDistributionMomModels
{

/*---------------------------------------------------------------------------*\
                      Class LaakkonenMomDaughterSizeDistribution Declaration
\*---------------------------------------------------------------------------*/

class LaakkonenMomDaughterSizeDistribution
:
    public daughterSizeDistributionMomModel
{
public:

    //- Runtime type information
    TypeName("LaakkonenMom");


    // Constructor

        LaakkonenMomDaughterSizeDistribution
        (
            const breakupMomModel& breakup,
            const dictionary& dict
        );


    //- Destructor
    virtual ~LaakkonenMomDaughterSizeDistribution();


    // Member Functions

        virtual void calcIntegrald
        (
            volScalarField& daughterIntegral,
            const volScalarField& di,
            const label i,
            const label k
        ) const;

        virtual void calcIntegralm
        (
            volScalarField& daughterIntegral,
            const volScalarField& mi,
            const label i,
            const label k
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace daughterSizeDistributionMomModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
