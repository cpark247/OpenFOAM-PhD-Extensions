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

\*---------------------------------------------------------------------------*/

#include "phaseSystem.H"
#include "phaseSystemHaensch.H"
#include "MomentumTransferPhaseSystem.H"
#include "MomentumTransferPhaseSystemExtend.H"
#include "OneResistanceHeatTransferPhaseSystem.H"
#include "TwoResistanceHeatTransferPhaseSystem.H"
#include "PhaseTransferPhaseSystem.H"
#include "InterfaceCompositionPhaseChangePhaseSystem.H"
#include "PopulationBalancePhaseSystem.H"
#include "PopulationBalanceMomPhaseSystem.H"
#include "ThermalPhaseChangePhaseSystem.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef
        PhaseTransferPhaseSystem
        <
            OneResistanceHeatTransferPhaseSystem
            <
                MomentumTransferPhaseSystem<phaseSystemHaensch>
            >
        >
        basicMultiphaseSystemHaensch;

    addNamedToRunTimeSelectionTable
    (
        phaseSystem,
        basicMultiphaseSystemHaensch,
        dictionary,
        basicMultiphaseSystemHaensch
    );

    typedef
        PopulationBalancePhaseSystem
        <
            PhaseTransferPhaseSystem
            <
                OneResistanceHeatTransferPhaseSystem
                <
                    MomentumTransferPhaseSystem<phaseSystemHaensch>
                >
            >
        >
        populationBalanceMultiphaseSystemHaensch;

    addNamedToRunTimeSelectionTable
    (
        phaseSystem,
        populationBalanceMultiphaseSystemHaensch,
        dictionary,
        populationBalanceMultiphaseSystemHaensch
    );

    typedef
        PhaseTransferPhaseSystem
        <
            OneResistanceHeatTransferPhaseSystem
            <
                MomentumTransferPhaseSystemExtend<phaseSystemHaensch>
            >
        >
        basicMultiphaseSystemHaenschExtend;

    addNamedToRunTimeSelectionTable
    (
        phaseSystem,
        basicMultiphaseSystemHaenschExtend,
        dictionary,
        basicMultiphaseSystemHaenschExtend
    );

    typedef
        PopulationBalancePhaseSystem
        <
            PhaseTransferPhaseSystem
            <
                OneResistanceHeatTransferPhaseSystem
                <
                    MomentumTransferPhaseSystemExtend<phaseSystemHaensch>
                >
            >
        >
        populationBalanceMultiphaseSystemHaenschExtend;

    addNamedToRunTimeSelectionTable
    (
        phaseSystem,
        populationBalanceMultiphaseSystemHaenschExtend,
        dictionary,
        populationBalanceMultiphaseSystemHaenschExtend
    );

    typedef
        PopulationBalanceMomPhaseSystem
        <
            PhaseTransferPhaseSystem
            <
                OneResistanceHeatTransferPhaseSystem
                <
                    MomentumTransferPhaseSystem<phaseSystem>
                >
            >
        >
        populationBalanceMomMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        phaseSystem,
        populationBalanceMomMultiphaseSystem,
        dictionary,
        populationBalanceMomMultiphaseSystem
    );

    typedef
        PopulationBalanceMomPhaseSystem
        <
            PhaseTransferPhaseSystem
            <
                OneResistanceHeatTransferPhaseSystem
                <
                    MomentumTransferPhaseSystem<phaseSystemHaensch>
                >
            >
        >
        populationBalanceMomMultiphaseSystemHaensch;

    addNamedToRunTimeSelectionTable
    (
        phaseSystem,
        populationBalanceMomMultiphaseSystemHaensch,
        dictionary,
        populationBalanceMomMultiphaseSystemHaensch
    );

    typedef
        PopulationBalanceMomPhaseSystem
        <
            PhaseTransferPhaseSystem
            <
                OneResistanceHeatTransferPhaseSystem
                <
                    MomentumTransferPhaseSystemExtend<phaseSystemHaensch>
                >
            >
        >
        populationBalanceMomMultiphaseSystemHaenschExtend;

    addNamedToRunTimeSelectionTable
    (
        phaseSystem,
        populationBalanceMomMultiphaseSystemHaenschExtend,
        dictionary,
        populationBalanceMomMultiphaseSystemHaenschExtend
    );
    



}

// ************************************************************************* //
