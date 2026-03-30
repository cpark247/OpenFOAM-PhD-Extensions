/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
-------------------------------------------------------------------------------
Description
    Updates time-averaged population balance rate fields.

    Uses size group names discovered by createPopBalRateFields.H.

    Include this in the time loop AFTER the population balance model has been
    solved (i.e., after fluid.solve() and fluid.correct()).

    Looks up instantaneous rate fields created by:
    - LuoSvendsenForSAS2: binaryBreakupRate_<sizeGroupName>, binaryBreakupRateTotal
    - PrinceBlanchForSAS2: coalescenceRate_<sizeGroupName>, coalescenceRateTotal

    Updates time-averaged fields created by createPopBalRateFields.H

    Requires from createPopBalRateFields.H:
    - sizeGroupNames (List<word>)
    - nSizeClasses (label)
    - binaryBreakupRate_Mean_ (PtrList<volScalarField>)
    - coalescenceRate_Mean_ (PtrList<volScalarField>)
    - binaryBreakupRateTotal_Mean_ (volScalarField)
    - coalescenceRateTotal_Mean_ (volScalarField)

\*---------------------------------------------------------------------------*/

// ========== POPULATION BALANCE RATE AVERAGING ==========

if (nSizeClasses > 0)
{
    // --- Per-class binary breakup rates ---
    forAll(sizeGroupNames, i)
    {
        // Instantaneous field name from LuoSvendsenForSAS2: binaryBreakupRate_f0.air, etc.
        word instFieldName = "binaryBreakupRate_" + sizeGroupNames[i];
        
        if (mesh.foundObject<volScalarField>(instFieldName))
        {
            if (!popBalRateFieldsFound)
            {
                Info<< "Found population balance rate fields - starting averaging" << endl;
                popBalRateFieldsFound = true;
            }
            
            const volScalarField& instField = 
                mesh.lookupObject<volScalarField>(instFieldName);
            
            averager.updateAverage(binaryBreakupRate_Mean_[i], instField, currentDeltaT);
        }
    }

    // --- Per-class coalescence rates ---
    forAll(sizeGroupNames, i)
    {
        // Instantaneous field name from PrinceBlanchForSAS2: coalescenceRate_f0.air, etc.
        word instFieldName = "coalescenceRate_" + sizeGroupNames[i];
        
        if (mesh.foundObject<volScalarField>(instFieldName))
        {
            const volScalarField& instField = 
                mesh.lookupObject<volScalarField>(instFieldName);
            
            averager.updateAverage(coalescenceRate_Mean_[i], instField, currentDeltaT);
        }
    }
}

// --- Total binary breakup rate (always check, even if nSizeClasses == 0) ---
if (mesh.foundObject<volScalarField>("binaryBreakupRateTotal"))
{
    const volScalarField& binaryBreakupRateTotal = 
        mesh.lookupObject<volScalarField>("binaryBreakupRateTotal");
    
    averager.updateAverage(binaryBreakupRateTotal_Mean_, binaryBreakupRateTotal, currentDeltaT);
}

// --- Total coalescence rate (always check, even if nSizeClasses == 0) ---
if (mesh.foundObject<volScalarField>("coalescenceRateTotal"))
{
    const volScalarField& coalescenceRateTotal = 
        mesh.lookupObject<volScalarField>("coalescenceRateTotal");
    
    averager.updateAverage(coalescenceRateTotal_Mean_, coalescenceRateTotal, currentDeltaT);
}

