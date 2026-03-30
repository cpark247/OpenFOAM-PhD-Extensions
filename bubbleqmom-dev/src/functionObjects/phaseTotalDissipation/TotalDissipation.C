/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "TotalDissipation.H"
#include "momentumTransportModel.H"
#include "dynamicMomentumTransportModel.H"
#include "phaseDynamicMomentumTransportModel.H"
//#include "turbulentFluidThermoModel.H"
#include "ThermoPhaseModel.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"

#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(TotalDissipation, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        TotalDissipation,
        dictionary
    );

    const Foam::word Foam::functionObjects::TotalDissipation::turbulenceModelName
    (
        Foam::momentumTransportModel::typeName
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::TotalDissipation::compressible()
{
    if (obr_.foundObject<phaseCompressible::momentumTransportModel>
			    (
				IOobject::groupName(momentumTransportModel::typeName, phaseName_)
			    )
			    )
    {
	return true;
    }
    else
    {
        FatalErrorInFunction
            << "Turbulence model not found in database, deactivating"
            << exit(FatalError);
    }

    return false;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::TotalDissipation::TotalDissipation
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldName_("epsilonTotal"),
    phaseName_(dict.lookupOrDefault<word>("phase", "water")),
    active_(true)
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "TotalDissipation::TotalDissipation"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating." << nl
            << endl;
    }

    read(dict);


    if (active_)
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        volScalarField* TotalDissipationPtr
        (
            new volScalarField
            (
                IOobject
                (
                    fieldName_+"."+phaseName_,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("0", sqr(dimLength)/pow(dimTime,3), scalar(0.0))
            )
        );
	
        mesh.objectRegistry::store(TotalDissipationPtr);
        Info << "Store " 
	     << fieldName_ 
	     << "."
	     << phaseName_ 
	     << " in object registry." 
	     << endl << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::TotalDissipation::~TotalDissipation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::TotalDissipation::read(const dictionary&)
{
    // Do nothing
    return true;
}


bool Foam::functionObjects::TotalDissipation::execute()
{
    bool comp = compressible();

    if (!active_)
    {
        return false;
    }

    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    const volVectorField& U = mesh.lookupObject<volVectorField>("U."+phaseName_);

    const volSymmTensorField Sres(symm(fvc::grad(U)));

    const volScalarField& kMod = mesh.lookupObject<volScalarField>("k."+phaseName_);

    const volScalarField& omegaMod = mesh.lookupObject<volScalarField>("omega."+phaseName_);

    volScalarField& epsilonTotal =
        const_cast<volScalarField&>
        (
	    mesh.lookupObject<volScalarField>(fieldName_+"."+phaseName_)
	);
	
    // Retrieve turbulence properties from model
    const phaseCompressible::momentumTransportModel& model =
        mesh.lookupObject<phaseCompressible::momentumTransportModel>
	(
	    IOobject::groupName(momentumTransportModel::typeName, phaseName_)
        );
 
    dimensionedScalar unitDynamicViscosity	
    (
   	dimensionSet(0,2,-1,0,0,0,0),
    	1.0
    );

   const scalar Cmu = 0.09;
   const volScalarField& epsilonMod = Cmu*kMod*omegaMod;

   const volScalarField& epsilonRes = 2*model.nu()*(Sres && Sres);

   epsilonTotal = epsilonRes + epsilonMod;
   
   return true;
    
}


bool Foam::functionObjects::TotalDissipation::end()
{
    if (active_)
    {
        execute();
    }

    return true;
}


bool Foam::functionObjects::TotalDissipation::write()
{
    if (active_)
    {
        regIOobject& obj =
            const_cast<regIOobject&>
            (
                obr_.lookupObject<regIOobject>(fieldName_+"."+phaseName_)
            );

        obj.write();
    }

    return true;
}


// ************************************************************************* //
