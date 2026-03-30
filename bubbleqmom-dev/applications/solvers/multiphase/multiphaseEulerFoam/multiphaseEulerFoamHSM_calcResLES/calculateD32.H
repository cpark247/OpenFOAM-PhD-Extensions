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
    Calculates Sauter mean diameter (d_32) from population balance size groups
    and updates time-averaged values for both d_32 and size fractions.

    Uses size group names and diameters discovered/read by createPopBalRateFields.H.

    d_32 = 1 / sum(f_i / d_i)

    Include this in the time loop AFTER the population balance has been solved.

    Requires from createPopBalRateFields.H:
    - sizeGroupNames (List<word>)
    - sizeGroupDiameters (List<scalar>)
    - nSizeClasses (label)
    - f_Mean_ (PtrList<volScalarField>)
    - fFieldDimensions (dimensionSet)

\*---------------------------------------------------------------------------*/

// ========== D32 AND SIZE FRACTION CALCULATIONS ==========

if (nSizeClasses > 0)
{
    // Compute the dimension for f/d
    dimensionSet denominatorDim = fFieldDimensions / dimLength;
    
    // Initialize denominator with correct dimensions
    volScalarField denominator
    (
        IOobject
        (
            "d32_denominator_tmp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false  // Don't register
        ),
        mesh,
        dimensionedScalar("zero", denominatorDim, 0.0)
    );

    // Loop through all size groups
    forAll(sizeGroupNames, i)
    {
        const word& fieldName = sizeGroupNames[i];
        const scalar diameter = sizeGroupDiameters[i];
        
        if (mesh.foundObject<volScalarField>(fieldName))
        {
            const volScalarField& fi = mesh.lookupObject<volScalarField>(fieldName);
            
            // Update time-averaged size fraction (dimensions match via fFieldDimensions)
            averager.updateAverage(f_Mean_[i], fi, currentDeltaT);
            
            // Add contribution to denominator: f_i / d_i
            denominator += fi / dimensionedScalar("di", dimLength, diameter);
        }
        else
        {
            WarningInFunction
                << "Size group field " << fieldName << " not found in mesh"
                << endl;
        }
    }

    // d_32 = 1 / sum(f_i / d_i)
    // Result dimension = dimless / denominatorDim
    //                  = dimless / (fFieldDimensions / dimLength)
    //                  = dimLength / fFieldDimensions
    // If fFieldDimensions = dimless: d_32 has dim [length] ✓
    // If fFieldDimensions = dimLength: d_32 has dim [dimless] - needs adjustment
    
    if (fFieldDimensions == dimless)
    {
        // Standard case: f is dimensionless (volume fraction)
        d_32 = dimensionedScalar("one", dimless, 1.0) / max
        (
            denominator,
            dimensionedScalar("smallDenom", denominatorDim, 1.0)
        );
    }
    else
    {
        // Non-standard case: need to adjust formula to get dimLength result
        // Compute sum(f_i) for numerator
        volScalarField sumF
        (
            IOobject
            (
                "sumF_tmp",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", fFieldDimensions, 0.0)
        );
        
        forAll(sizeGroupNames, i)
        {
            const word& fieldName = sizeGroupNames[i];
            if (mesh.foundObject<volScalarField>(fieldName))
            {
                const volScalarField& fi = mesh.lookupObject<volScalarField>(fieldName);
                sumF += fi;
            }
        }
        
        // d_32 = sum(f_i) / sum(f_i/d_i)
        // Dimension: fFieldDimensions / (fFieldDimensions/dimLength) = dimLength
        d_32 = sumF / max
        (
            denominator,
            dimensionedScalar("smallDenom", denominatorDim, 1e-10)
        );
    }

    // Time-averaged d_32
    averager.updateAverage(d_32_Mean_, d_32, currentDeltaT);
}
