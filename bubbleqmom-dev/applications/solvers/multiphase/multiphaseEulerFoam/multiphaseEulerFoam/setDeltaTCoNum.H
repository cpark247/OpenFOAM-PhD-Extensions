/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Global
    setDeltaT

Description
    Reset the timestep to maintain a constant maximum courant Number.
    Reduction of time-step is immediate, but increase is damped to avoid
    unstable oscillations.

\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{
    scalar tmpMaxCo = maxCo;

    HashTable<diameterModels::populationBalanceMomModel*> popBal = mesh.lookupClass<diameterModels::populationBalanceMomModel>();
    scalarList realizableCoList(popBal.size(), 1.0);

    forAllIter(HashTable<diameterModels::populationBalanceMomModel*>, popBal, iter)
    {
	realizableCoList.append(iter()->realizableCo());
    }
    scalar realizableCo = gMin(realizableCoList);

    if(realizableCo < maxCo )
    {
        tmpMaxCo = realizableCo;
        Info<<"correct Co number, new Co "<<tmpMaxCo<<"\n";
    }

    scalar maxDeltaTFact = tmpMaxCo/(CoNum + small);
    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    runTime.setDeltaT
    (
        min
        (
            deltaTFact*runTime.deltaTValue(),
            maxDeltaT
        )
    );
    Info<< "maxCo = " << maxCo << endl;
    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
}

// ************************************************************************* //
