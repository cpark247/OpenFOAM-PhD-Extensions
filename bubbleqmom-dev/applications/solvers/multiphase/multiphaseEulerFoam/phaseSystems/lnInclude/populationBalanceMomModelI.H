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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::phaseSystem&
Foam::diameterModels::populationBalanceMomModel::fluid() const
{
    return fluid_;
}


inline const Foam::fvMesh&
Foam::diameterModels::populationBalanceMomModel::mesh() const
{
    return mesh_;
}


inline const Foam::dictionary&
Foam::diameterModels::populationBalanceMomModel::dict() const
{
    return dict_;
}


inline const Foam::phaseModel&
Foam::diameterModels::populationBalanceMomModel::continuousPhase() const
{
    return continuousPhase_;
}

inline const Foam::diameterModels::populationBalanceMomModel::phasePairTable&
Foam::diameterModels::populationBalanceMomModel::phasePairs() const
{
    return phasePairs_;
}

inline const Foam::volScalarField&
Foam::diameterModels::populationBalanceMomModel::alphas() const
{
    return alphas_();
}

/*
inline const Foam::volVectorField&
Foam::diameterModels::populationBalanceMomModel::U() const
{
    return U_();
}
*/
inline Foam::monoKineticQuadratureApproximation&
Foam::diameterModels::populationBalanceMomModel::getQuadratureObj()
{
    return quadrature_();
}

inline const Foam::monoKineticQuadratureApproximation&
Foam::diameterModels::populationBalanceMomModel::getQuadratureObj() const
{
    return quadrature_();
}

inline label
Foam::diameterModels::populationBalanceMomModel::nNodes() const
{
    return quadrature_().nodes().size();
}

inline const volScalarField&
Foam::diameterModels::populationBalanceMomModel::nodeAlphas(const label nodei) const
{
	return nodeAlphas_[nodei];
}

inline volScalarField&
Foam::diameterModels::populationBalanceMomModel::nodeAlphas(const label nodei)
{
	return nodeAlphas_[nodei];
}

inline const volScalarField&
Foam::diameterModels::populationBalanceMomModel::ds(const label nodei) const
{
    return ds_[nodei];
}


inline const volScalarField&
Foam::diameterModels::populationBalanceMomModel::dsm() const
{
	return dsm_();
}



// ************************************************************************* //
