/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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
    Foam::RASModels::kOmegaSSTSAS_ext

Description
    Scale-adaptive URAS model based on the k-omega-SST RAS model.

    References:
    \verbatim
        Egorov, Y., & Menter F.R. (2008).
        Development and Application of SST-SAS Model in the DESIDER Project.
        Advances in Hybrid RANS-LES Modelling,
        Notes on Num. Fluid Mech. And Multidisciplinary Design,
        Volume 97, 261-270.
    \endverbatim

    The model coefficients are
    \verbatim
        kOmegaSSTSAS_extCoeffs
        {
            // Default SST coefficients
            alphaK1     0.85;
            alphaK2     1.0;
            alphaOmega1 0.5;
            alphaOmega2 0.856;
            beta1       0.075;
            beta2       0.0828;
            betaStar    0.09;
            gamma1      5/9;
            gamma2      0.44;
            a1          0.31;
            b1          1.0;
            c1          10.0;
            F3          no;

            // Default SAS coefficients
            Cs          0.11;
            kappa       0.41;
            zeta2       3.51;
            sigmaPhi    2.0/3.0;
            C           2;

            // Delta must be specified for SAS e.g.
            delta cubeRootVol;

            cubeRootVolCoeffs
            {}
        }
    \endverbatim

SourceFiles
    kOmegaSSTSAS_ext.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSSTSAS_ext_H
#define kOmegaSSTSAS_ext_H

#include "kOmegaSST_ext.H"
#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                        Class kOmegaSSTSAS_ext Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class kOmegaSSTSAS_ext
:
    public kOmegaSST_ext<BasicMomentumTransportModel>
{
protected:

    // Protected data

        // Model constants

            dimensionedScalar nExponent_;
            dimensionedScalar Cs_;
            dimensionedScalar kappa_;
            dimensionedScalar zeta2_;
            dimensionedScalar sigmaPhi_;
            dimensionedScalar C_;


        // Fields

            //- Run-time selectable delta model
            autoPtr<Foam::LESdelta> delta_;

            //- Turbulent length scale
            volScalarField L_;
            
            //- von Karman length scale
            volScalarField Lvk_;
			
			// L / Lvk
			volScalarField L_by_Lvk_;

            //- SAS source term
            volScalarField Qsas_;


    // Protected Member Functions

        //- SAS omega source
        virtual tmp<fvScalarMatrix> Qsas
        (
            const volScalarField::Internal& S2,
            const volScalarField::Internal& gamma,
            const volScalarField::Internal& beta
        ) const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaSSTSAS_ext");


    // Constructors

        //- Construct from components
        kOmegaSSTSAS_ext
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& type = typeName
        );

        //- Disallow default bitwise copy construction
        kOmegaSSTSAS_ext(const kOmegaSSTSAS_ext&) = delete;


    //- Destructor
    virtual ~kOmegaSSTSAS_ext()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Access function to filter width
        inline const volScalarField& delta() const
        {
            return delta_();
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const kOmegaSSTSAS_ext&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kOmegaSSTSAS_ext.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
