/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "daughterSizeDistributionMomModel.H"
#include "breakupMomModel.H"

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(daughterSizeDistributionMomModel, 0);
    defineRunTimeSelectionTable(daughterSizeDistributionMomModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::diameterModels::daughterSizeDistributionMomModel>
Foam::diameterModels::daughterSizeDistributionMomModel::New
(
    const breakupMomModel& breakup,
    const dictionary& dict
)
{
    word daughterSizeDistributionMomModelType
    (
        dict.lookup("daughterSizeDistributionMomModel")
    );

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(daughterSizeDistributionMomModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown daughter size distribution model type "
            << daughterSizeDistributionMomModelType << endl << endl
            << "Valid daughter size distribution model types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(breakup, dict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::diameterModels::daughterSizeDistributionMomModel::
daughterSizeDistributionMomModel
(
    const breakupMomModel& breakup,
    const dictionary& dict
)
:
    breakup_(breakup)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::daughterSizeDistributionMomModel::
~daughterSizeDistributionMomModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::diameterModels::daughterSizeDistributionMomModel::precompute()
{
	// ...
}


// ************************************************************************* //
