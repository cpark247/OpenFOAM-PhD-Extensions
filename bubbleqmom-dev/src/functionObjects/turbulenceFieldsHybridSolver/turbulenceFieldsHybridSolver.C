/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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

#include "turbulenceFieldsHybridSolver.H"
//#include "kinematicMomentumTransportModel.H"
#include "momentumTransportModel.H"
#include "dynamicMomentumTransportModel.H"
#include "phaseDynamicMomentumTransportModel.H"
#include "thermophysicalTransportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseSystem.H"
#include "ThermoPhaseModel.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(turbulenceFieldsHybridSolver, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        turbulenceFieldsHybridSolver,
        dictionary
    );
}
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::turbulenceFieldsHybridSolver::compressibleField,
    8
>::names[] =
{
    "k",
    "epsilon",
    "omega",
    "mut",
    "muEff",
    "alphaEff",
    "R",
    "devTau"
};

const Foam::NamedEnum
<
    Foam::functionObjects::turbulenceFieldsHybridSolver::compressibleField,
    8
> Foam::functionObjects::turbulenceFieldsHybridSolver::compressibleFieldNames_;

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::turbulenceFieldsHybridSolver::incompressibleField,
    7
>::names[] =
{
    "k",
    "epsilon",
    "omega",
    "nut",
    "nuEff",
    "R",
    "devSigma"
};

const Foam::NamedEnum
<
    Foam::functionObjects::turbulenceFieldsHybridSolver::incompressibleField,
    7
> Foam::functionObjects::turbulenceFieldsHybridSolver::incompressibleFieldNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::turbulenceFieldsHybridSolver::turbulenceFieldsHybridSolver
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldSet_(),
    phaseName_(dict.lookupOrDefault<word>("phase", word::null))
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::turbulenceFieldsHybridSolver::~turbulenceFieldsHybridSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::turbulenceFieldsHybridSolver::read(const dictionary& dict)
{
    if (dict.found("field"))
    {
        fieldSet_.insert(word(dict.lookup("field")));
    }
    else
    {
        fieldSet_.insert(wordList(dict.lookup("fields")));
    }

    if (dict.lookupOrDefault<Switch>("prefix", false))
    {
        prefix_ = momentumTransportModel::typeName + ':';
    }
    else
    {
        prefix_ = word::null;
    }

    Info<< type() << " " << name() << ": ";
    if (fieldSet_.size())
    {
        Info<< "storing fields:" << nl;
        forAllConstIter(wordHashSet, fieldSet_, iter)
        {
            Info<< "    "
                << IOobject::groupName(prefix_ + iter.key(), phaseName_) << nl;
        }
        Info<< endl;
    }
    else
    {
        Info<< "no fields requested to be stored" << nl << endl;
    }

    return true;
}


bool Foam::functionObjects::turbulenceFieldsHybridSolver::execute()
{
    const word modelName
    (
        IOobject::groupName(momentumTransportModel::typeName, phaseName_)
    );
 
    
    if (obr_.foundObject<phaseCompressible::momentumTransportModel>(modelName))
    {
        const phaseCompressible::momentumTransportModel& model =
            obr_.lookupObject<phaseCompressible::momentumTransportModel>(modelName);

        forAllConstIter(wordHashSet, fieldSet_, iter)
        {
            const word& f = iter.key();

            switch (compressibleFieldNames_[f])
            {
                case compressibleField::k:
                {
		    processField<scalar>(f, model.k());
                    break;
                }
                case compressibleField::epsilon:
                {
                    processField<scalar>(f, model.epsilon());
                    break;
                }
                case compressibleField::omega:
                {
		    const volScalarField k = model.k();
		    const volScalarField epsilon = model.epsilon();
		    const volScalarField omegaCal = epsilon/(0.09*max(k,dimensionedScalar("mink", k.dimensions(), 1e-14)));
                    processField<scalar>(f, omegaCal);
                    break;
                }
                case compressibleField::mut:
                {
                    processField<scalar>(f, model.mut());
                    break;
                }
                case compressibleField::muEff:
                {
                    processField<scalar>(f, model.muEff());
                    break;
                }
                case compressibleField::R:
                {
                    processField<symmTensor>(f, model.sigma());
                    break;
                }
                case compressibleField::devTau:
                {
                    processField<symmTensor>(f, model.devTau());
                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Invalid field selection" << exit(FatalError);
                }
            }
        }
    }

    else
    {
        FatalErrorInFunction
            << "Turbulence model not found in database, deactivating"
            << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::turbulenceFieldsHybridSolver::write()
{
    forAllConstIter(wordHashSet, fieldSet_, iter)
    {
        const word fieldName
        (
            IOobject::groupName(prefix_ + iter.key(), phaseName_)
        );
        writeObject(fieldName);
    }

    return true;
}


// ************************************************************************* //
