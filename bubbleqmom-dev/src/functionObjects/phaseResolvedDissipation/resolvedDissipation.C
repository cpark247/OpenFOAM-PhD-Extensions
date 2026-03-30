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

#include "resolvedDissipation.H"
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
    defineTypeNameAndDebug(resolvedDissipation, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        resolvedDissipation,
        dictionary
    );

    const Foam::word Foam::functionObjects::resolvedDissipation::turbulenceModelName
    (
        Foam::momentumTransportModel::typeName
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::resolvedDissipation::compressible()
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

Foam::functionObjects::resolvedDissipation::resolvedDissipation
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldName_("epsilonRes"),
    phaseName_(dict.lookupOrDefault<word>("phase", "water")),
    active_(true)
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "resolvedDissipation::resolvedDissipation"
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

        volScalarField* resolvedDissipationPtr
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
	
        mesh.objectRegistry::store(resolvedDissipationPtr);
        Info << "Store " 
	     << fieldName_ 
	     << "."
	     << phaseName_ 
	     << " in object registry." 
	     << endl << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::resolvedDissipation::~resolvedDissipation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::resolvedDissipation::read(const dictionary&)
{
    // Do nothing
    return true;
}


bool Foam::functionObjects::resolvedDissipation::execute()
{
    bool comp = compressible();

    if (!active_)
    {
        return false;
    }

    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    const volVectorField& U = mesh.lookupObject<volVectorField>("U."+phaseName_);
    
    const volSymmTensorField Sres(symm(fvc::grad(U)));
    volScalarField& epsilonRes =
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
  
   epsilonRes = 2*model.nu()*(Sres && Sres);
   
   return true;
    
}


bool Foam::functionObjects::resolvedDissipation::end()
{
    if (active_)
    {
        execute();
    }

    return true;
}


bool Foam::functionObjects::resolvedDissipation::write()
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
