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
    Updates population balance rate fields.

    Uses size group names discovered by createPopBalRateFields.H.

    Include this in the time loop AFTER the population balance model has been
    solved (i.e., after fluid.solve() and fluid.correct()).

    Looks up instantaneous rate fields created by:
    - LuoSvendsenForSAS2: binaryBreakupRate_<sizeGroupName>, binaryBreakupRateTotal
    - PrinceBlanchForSAS2: coalescenceRate_<sizeGroupName>, coalescenceRateTotal

    Per size class the INSTANTANEOUS rates are written out:
    - binaryBreakupRate_<sizeGroupName>_Inst is the instantaneous breakup rate
      scaled by the daughter size-group volume width deltaV_i [-> 1/s].
    - coalescenceRate_<sizeGroupName>_Inst is the instantaneous coalescence rate.
    The time-averaged TOTAL rates (_Mean_) are retained, and instantaneous totals
    are written out as well:
    - binaryBreakupRateTotal_Inst = sum_i (binaryBreakupRate_i * deltaV_i) [1/s].
    - coalescenceRateTotal_Inst   = instantaneous copy of coalescenceRateTotal.

    Requires from createPopBalRateFields.H:
    - sizeGroupNames (List<word>)
    - nSizeClasses (label)
    - sizeGroupVolumeWidths (List<scalar>)
    - binaryBreakupRate_Inst_ (PtrList<volScalarField>)
    - coalescenceRate_Inst_ (PtrList<volScalarField>)
    - binaryBreakupRateTotal_Mean_ (volScalarField)
    - coalescenceRateTotal_Mean_ (volScalarField)
    - binaryBreakupRateTotal_Inst (volScalarField)
    - coalescenceRateTotal_Inst (volScalarField)

\*---------------------------------------------------------------------------*/

// ========== POPULATION BALANCE RATE OUTPUT ==========

if (nSizeClasses > 0)
{
    // --- Per-class binary breakup rates (instantaneous, scaled by deltaV) ---
    // Also accumulate the deltaV-scaled total (sum over classes -> 1/s).
    binaryBreakupRateTotal_Inst ==
        dimensionedScalar("zero", binaryBreakupRateTotal_Inst.dimensions(), 0.0);

    forAll(sizeGroupNames, i)
    {
        // Instantaneous field name from LuoSvendsenForSAS2: binaryBreakupRate_f0.air, etc.
        word instFieldName = "binaryBreakupRate_" + sizeGroupNames[i];

        if (mesh.foundObject<volScalarField>(instFieldName))
        {
            if (!popBalRateFieldsFound)
            {
                Info<< "Found population balance rate fields - writing instantaneous rates" << endl;
                popBalRateFieldsFound = true;
            }

            const volScalarField& instField =
                mesh.lookupObject<volScalarField>(instFieldName);

            // Scale by the daughter size-group volume width deltaV_i so the
            // written quantity has dimensions 1/s. deltaV carries dimVolume so
            // the field dimensions resolve to inv(dimTime).
            const dimensionedScalar dV
            (
                "dV", dimVolume, sizeGroupVolumeWidths[i]
            );

            binaryBreakupRate_Inst_[i] = instField*dV;

            // Accumulate into the deltaV-scaled instantaneous total
            binaryBreakupRateTotal_Inst += binaryBreakupRate_Inst_[i];
        }
    }

    // --- Per-class coalescence rates (instantaneous) ---
    forAll(sizeGroupNames, i)
    {
        // Instantaneous field name from PrinceBlanchForSAS2: coalescenceRate_f0.air, etc.
        word instFieldName = "coalescenceRate_" + sizeGroupNames[i];

        if (mesh.foundObject<volScalarField>(instFieldName))
        {
            const volScalarField& instField =
                mesh.lookupObject<volScalarField>(instFieldName);

            coalescenceRate_Inst_[i] = instField;
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

    // Instantaneous copy of the registered total
    coalescenceRateTotal_Inst = coalescenceRateTotal;
}
