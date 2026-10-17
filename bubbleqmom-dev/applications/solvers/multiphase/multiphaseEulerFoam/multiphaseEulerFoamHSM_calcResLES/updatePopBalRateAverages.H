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
    solved (i.e., after fluid.solve() and fluid.correct()) and AFTER
    calculateD32.H, which provides the d_32 needed for the coalescence
    conversion below.

    Every field written here is a FREQUENCY [1/s], so the time scales of
    Buffo et al. (2016), doi:10.1016/j.ces.2016.07.032, are simply

        tau_breakup     = 1/binaryBreakupRateTotal_Inst
        tau_coalescence = 1/coalescenceRateTotal_Inst

    Breakup
        LuoSvendsenForSAS2 registers breakup FREQUENCIES already: it integrates
        its rate density Omega(v_i : v_j) [1/(m^3 s)] over the daughter volume
        bin inside the model, using the same volume widths that
        populationBalanceModel::calcDeltas() uses. The fields are therefore
        copied straight through.
            binaryBreakupRate_<sizeGroupName>  = g_j, frequency of parent j
            binaryBreakupRateTotal             = sum_j f_j*g_j

    Coalescence
        PrinceBlanchForSAS2 registers coalescence KERNELS [m^3/s], size-fraction
        weighted:
            coalescenceRate_<sizeGroupName>    = sum_j f_j*a_ij
            coalescenceRateTotal               = sum_i sum_j f_i*f_j*a_ij
        They become frequencies through the number density of the dispersed
        phase built from the Sauter mean diameter, exactly as in the reference:
            n = alpha_d/V_32,   V_32 = (pi/6)*d_32^3

    Requires from createPopBalRateFields.H:
    - sizeGroupNames (List<word>)
    - nSizeClasses (label)
    - dispersedPhaseName (word)
    - binaryBreakupRate_Inst_ (PtrList<volScalarField>)
    - coalescenceRate_Inst_ (PtrList<volScalarField>)
    - binaryBreakupRateTotal_Mean_ (volScalarField)
    - coalescenceRateTotal_Mean_ (volScalarField)
    - binaryBreakupRateTotal_Inst (volScalarField)
    - coalescenceRateTotal_Inst (volScalarField)

    Requires from createMultiphaseResolvedLESFields.H / calculateD32.H:
    - d_32 (volScalarField)

\*---------------------------------------------------------------------------*/

// ========== POPULATION BALANCE RATE OUTPUT ==========

if (nSizeClasses > 0)
{
    // Number density of the dispersed phase from the Sauter mean diameter,
    //     n = alpha_d/((pi/6)*d_32^3)
    const scalar piBy6Local = 3.14159265358979323846/6.0;

    const volScalarField& alphaDispersed =
        mesh.lookupObject<volScalarField>("alpha." + dispersedPhaseName);

    const volScalarField nDispersed
    (
        alphaDispersed
       /max
        (
            piBy6Local*Foam::pow3(d_32),
            dimensionedScalar("minV32", dimVolume, 1e-30)
        )
    );

    // --- Per-class breakup frequencies [1/s] ---
    forAll(sizeGroupNames, i)
    {
        // Field name from LuoSvendsenForSAS2: binaryBreakupRate_f0.air, etc.
        word instFieldName = "binaryBreakupRate_" + sizeGroupNames[i];

        if (mesh.foundObject<volScalarField>(instFieldName))
        {
            if (!popBalRateFieldsFound)
            {
                Info<< "Found population balance rate fields - writing instantaneous rates" << endl;
                popBalRateFieldsFound = true;
            }

            binaryBreakupRate_Inst_[i] =
                mesh.lookupObject<volScalarField>(instFieldName);
        }
    }

    // --- Per-class coalescence frequencies [1/s] ---
    forAll(sizeGroupNames, i)
    {
        // Field name from PrinceBlanchForSAS2: coalescenceRate_f0.air, etc.
        word instFieldName = "coalescenceRate_" + sizeGroupNames[i];

        if (mesh.foundObject<volScalarField>(instFieldName))
        {
            coalescenceRate_Inst_[i] =
                mesh.lookupObject<volScalarField>(instFieldName)*nDispersed;
        }
    }

    // --- Total coalescence frequency [1/s] ---
    if (mesh.foundObject<volScalarField>("coalescenceRateTotal"))
    {
        coalescenceRateTotal_Inst =
            mesh.lookupObject<volScalarField>("coalescenceRateTotal")
           *nDispersed;

        averager.updateAverage
        (
            coalescenceRateTotal_Mean_,
            coalescenceRateTotal_Inst,
            currentDeltaT
        );
    }
}

// --- Total breakup frequency [1/s] (no conversion needed) ---
if (mesh.foundObject<volScalarField>("binaryBreakupRateTotal"))
{
    binaryBreakupRateTotal_Inst =
        mesh.lookupObject<volScalarField>("binaryBreakupRateTotal");

    averager.updateAverage
    (
        binaryBreakupRateTotal_Mean_,
        binaryBreakupRateTotal_Inst,
        currentDeltaT
    );
}
