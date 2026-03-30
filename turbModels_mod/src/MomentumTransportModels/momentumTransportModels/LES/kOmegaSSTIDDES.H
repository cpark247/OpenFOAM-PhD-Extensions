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

Description
    k-omega-SST IDDES turbulence model for
    incompressible and compressible flows.

    Reference:
    \verbatim
        Gritskevich, M. S., Garbaruk, A. V.,
        Schütze, J., & Menter, F. R. (2012).
        Development of DDES and IDDES formulations for
        the k-ω shear stress transport model.
        Flow, turbulence and combustion, 88(3), 431-449.
        DOI:10.1007/s10494-011-9378-4
    \endverbatim

SourceFiles
    kOmegaSSTIDDES.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSSTIDDES_H
#define kOmegaSSTIDDES_H

#include "kOmegaSSTDES.H"
#include "IDDESDelta.H"

// #include "LESModel.H"
// #include "LESeddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                          Class kOmegaSST Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class kOmegaSSTIDDES
:
    public kOmegaSSTDES<BasicMomentumTransportModel>    
{
    // Private data
    
        // Model constants
            dimensionedScalar Cdt1_;
            dimensionedScalar Cdt2_;
            dimensionedScalar Cl_;
            dimensionedScalar Ct_;
    
        // Fields
            const IDDESDelta& IDDESDelta_;
    
        // Private Member Functions
            // const IDDESDelta& setDelta() const;

            tmp<volScalarField> alpha() const;
            tmp<volScalarField> ft(const volScalarField& magGradU) const;
            tmp<volScalarField> fl(const volScalarField& magGradU) const;

            tmp<volScalarField> rd
            (
                const volScalarField& nur,
                const volScalarField& magGradU
            ) const;

            //- Delay function
            tmp<volScalarField> fdt(const volScalarField& magGradU) const;

            // Disallow default bitwise copy construct and assignment
            kOmegaSSTIDDES(const kOmegaSSTIDDES&);
            void operator=(const kOmegaSSTIDDES&);

protected:

    // Protected Member Functions

        //- Length scale
        virtual tmp<volScalarField> dTilda
        (
            const volScalarField& magGradU,
            const volScalarField& CDES
        ) const;
        
public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaSSTIDDES");


    // Constructors

        //- Construct from components
        kOmegaSSTIDDES
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            //const word& propertiesName = BasicMomentumTransportModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~kOmegaSSTIDDES()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "kOmegaSSTIDDES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
