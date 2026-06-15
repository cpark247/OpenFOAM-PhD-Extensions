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
    Creates time-averaged fields for population balance rates and size fractions.
    
    Discovers size groups from the mesh's registered fields (f*.air) and reads 
    diameters from phaseProperties dictionary.
    
    Filters out old time level fields (_0, _0_0 suffixes).
    
    Dynamically gets dimensions from actual mesh fields to avoid dimension
    mismatch errors.
    
    Include this AFTER createFields.H and createMultiphaseResolvedLESFields.H.

    Fields created:
    - binaryBreakupRate_<sizeGroupName>_Inst   (instantaneous, scaled by the
      daughter size-group volume width deltaV_i; e.g. binaryBreakupRate_f0.air_Inst)
    - coalescenceRate_<sizeGroupName>_Inst     (instantaneous)
    - binaryBreakupRateTotal_Mean_             (time-averaged)
    - coalescenceRateTotal_Mean_               (time-averaged)
    - binaryBreakupRateTotal_Inst              (instantaneous, sum of deltaV-scaled
      per-class breakup rates -> 1/s)
    - coalescenceRateTotal_Inst                (instantaneous copy of the total)
    - <sizeGroupName>_Mean_                     (time-averaged; e.g., f0.air_Mean_)

    Volume width per size group:
        deltaV_i = (pi/6)*(d_{i+1/2}^3 - d_{i-1/2}^3)
    with diameter boundaries taken as midpoints between representative
    diameters and the two edge bins extended symmetrically (lower bound
    clipped at zero).

\*---------------------------------------------------------------------------*/

// ========== STORAGE FOR SIZE GROUP INFORMATION ==========

label nSizeClasses = 0;
List<word> sizeGroupNames;
List<scalar> sizeGroupDiameters;

// Volume width deltaV_i of each size group (used to scale breakup rates on
// output); populated below once diameters are known.
List<scalar> sizeGroupVolumeWidths;

// Store the dimension of f fields (will be read from actual mesh fields)
dimensionSet fFieldDimensions(dimless);

// Flag to track if pop balance rate fields have been found
bool popBalRateFieldsFound = false;

// Declare PtrLists at outer scope (will be populated if size groups exist)
PtrList<volScalarField> binaryBreakupRate_Inst_;
PtrList<volScalarField> coalescenceRate_Inst_;
PtrList<volScalarField> f_Mean_;

// ========== DISCOVER SIZE GROUPS FROM MESH AND READ DIAMETERS ==========

{
    // Step 1: Find all f*.air fields registered in the mesh
    // These were created by createFields.H -> phaseSystem -> velocityGroup
    // Filter out old time level fields (ending with _0, _0_0, etc.)
    
    DynamicList<word> foundNames;
    
    // Search through mesh's object registry for volScalarFields
    HashTable<volScalarField*> scalarFields = 
        mesh.lookupClass<volScalarField>();
    
    forAllIter(HashTable<volScalarField*>, scalarFields, iter)
    {
        const word& fieldName = iter.key();
        
        // Check if field matches pattern f[0-9]+.air EXACTLY (no suffix)
        // Field must end with ".air" - filter out "_0", "_0_0" old time levels
        
        if (fieldName.size() >= 6)  // minimum: "f0.air" = 6 chars
        {
            // Check if name ENDS with ".air" (not just contains it)
            string suffix = ".air";
            if (fieldName.size() >= suffix.size())
            {
                bool endsWithAir = (fieldName.substr(fieldName.size() - suffix.size()) == suffix);
                
                if (endsWithAir && fieldName[0] == 'f')
                {
                    // Check that characters between 'f' and ".air" are all digits
                    size_t dotPos = fieldName.size() - suffix.size();  // position of '.'
                    
                    if (dotPos > 1)  // at least one digit between 'f' and '.'
                    {
                        bool allDigits = true;
                        for (size_t i = 1; i < dotPos; i++)
                        {
                            if (!std::isdigit(fieldName[i]))
                            {
                                allDigits = false;
                                break;
                            }
                        }
                        if (allDigits)
                        {
                            foundNames.append(fieldName);
                        }
                    }
                }
            }
        }
    }
    
    nSizeClasses = foundNames.size();
    
    if (nSizeClasses == 0)
    {
        Info << "No size group fields (f*.air) found in mesh - "
             << "population balance rate averaging disabled" << endl;
    }
    else
    {
        // Sort names to ensure consistent ordering (f0.air, f1.air, f2.air, ...)
        sizeGroupNames.setSize(nSizeClasses);
        forAll(foundNames, i)
        {
            sizeGroupNames[i] = foundNames[i];
        }
        
        // Sort by extracting the numeric part
        for (label i = 0; i < nSizeClasses - 1; i++)
        {
            for (label j = i + 1; j < nSizeClasses; j++)
            {
                // Extract numbers from names for comparison
                // Field format: f<number>.air
                label numI = std::stoi(sizeGroupNames[i].substr(1, sizeGroupNames[i].find(".") - 1));
                label numJ = std::stoi(sizeGroupNames[j].substr(1, sizeGroupNames[j].find(".") - 1));
                
                if (numJ < numI)
                {
                    word temp = sizeGroupNames[i];
                    sizeGroupNames[i] = sizeGroupNames[j];
                    sizeGroupNames[j] = temp;
                }
            }
        }
        
        Info << "Found " << nSizeClasses << " size group fields: " << sizeGroupNames << endl;
        
        // Get the dimensions of the first f field to use for f_Mean_ creation
        {
            const volScalarField& firstFField = 
                mesh.lookupObject<volScalarField>(sizeGroupNames[0]);
            fFieldDimensions = firstFField.dimensions();
            Info << "Size group field dimensions: " << fFieldDimensions << endl;
        }
        
        // Step 2: Read diameters from phaseProperties
        sizeGroupDiameters.setSize(nSizeClasses, 0.0);
        
        IOdictionary phaseProperties
        (
            IOobject
            (
                "phaseProperties",
                runTime.constant(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        
        if (phaseProperties.found("air"))
        {
            const dictionary& airDict = phaseProperties.subDict("air");
            
            if (airDict.found("velocityGroupCoeffs"))
            {
                const dictionary& velGroupDict = airDict.subDict("velocityGroupCoeffs");
                
                // Read the sizeGroups entry as a string and parse manually
                if (velGroupDict.found("sizeGroups"))
                {
                    // Get the entry and convert to string for parsing
                    OStringStream oss;
                    velGroupDict.lookupEntry("sizeGroups", false, false).write(oss);
                    string sizeGroupsStr = oss.str();
                    
                    Info << "Parsing size group diameters from phaseProperties..." << endl;
                    
                    // For each size group, find its dSph value
                    forAll(sizeGroupNames, i)
                    {
                        // Extract base name (f0, f1, etc.) from f0.air, f1.air, etc.
                        word baseName = sizeGroupNames[i].substr(0, sizeGroupNames[i].find("."));
                        
                        // Search for pattern: baseName{dSph VALUE; or baseName { dSph VALUE;
                        // Find the size group entry in the string
                        size_t pos = sizeGroupsStr.find(baseName);
                        
                        if (pos != std::string::npos)
                        {
                            // Find dSph after this position
                            size_t dSphPos = sizeGroupsStr.find("dSph", pos);
                            
                            if (dSphPos != std::string::npos)
                            {
                                // Find the value after dSph
                                // Skip "dSph" and any whitespace
                                size_t valueStart = dSphPos + 4;  // length of "dSph"
                                
                                // Skip whitespace
                                while (valueStart < sizeGroupsStr.size() && 
                                       (sizeGroupsStr[valueStart] == ' ' || 
                                        sizeGroupsStr[valueStart] == '\t'))
                                {
                                    valueStart++;
                                }
                                
                                // Find end of number (semicolon or whitespace)
                                size_t valueEnd = valueStart;
                                while (valueEnd < sizeGroupsStr.size() && 
                                       sizeGroupsStr[valueEnd] != ';' &&
                                       sizeGroupsStr[valueEnd] != ' ' &&
                                       sizeGroupsStr[valueEnd] != '\t' &&
                                       sizeGroupsStr[valueEnd] != '}')
                                {
                                    valueEnd++;
                                }
                                
                                if (valueEnd > valueStart)
                                {
                                    string valueStr = sizeGroupsStr.substr(valueStart, valueEnd - valueStart);
                                    sizeGroupDiameters[i] = std::stod(valueStr);
                                }
                            }
                        }
                        
                        Info << "    " << sizeGroupNames[i] 
                             << ": dSph = " << sizeGroupDiameters[i] << " m" << endl;
                    }
                }
            }
        }
        
        // Verify all diameters were found
        forAll(sizeGroupDiameters, i)
        {
            if (sizeGroupDiameters[i] <= 0)
            {
                WarningInFunction
                    << "Could not read diameter for " << sizeGroupNames[i]
                    << " - setting to 1e-3 m" << endl;
                sizeGroupDiameters[i] = 1e-3;
            }
        }

        // Step 3: Compute the volume width deltaV_i of each size group.
        // Diameter boundaries are midpoints between representative diameters;
        // the two edge bins are extended symmetrically and the lower bound is
        // clipped at zero. deltaV_i = (pi/6)*(dHi^3 - dLo^3).
        sizeGroupVolumeWidths.setSize(nSizeClasses, 0.0);

        const scalar piBy6 = 3.14159265358979323846/6.0;

        Info<< "Computing size-group volume widths (deltaV)..." << endl;

        forAll(sizeGroupNames, i)
        {
            const scalar di = sizeGroupDiameters[i];

            // Lower diameter boundary
            scalar dLo;
            if (i == 0)
            {
                const scalar dHi0 =
                    (nSizeClasses > 1)
                  ? 0.5*(di + sizeGroupDiameters[1])
                  : di;
                dLo = 2.0*di - dHi0;            // symmetric edge extension
            }
            else
            {
                dLo = 0.5*(sizeGroupDiameters[i-1] + di);
            }

            // Upper diameter boundary
            scalar dHi;
            if (i == nSizeClasses - 1)
            {
                const scalar dLoLast =
                    (nSizeClasses > 1)
                  ? 0.5*(sizeGroupDiameters[i-1] + di)
                  : di;
                dHi = 2.0*di - dLoLast;         // symmetric edge extension
            }
            else
            {
                dHi = 0.5*(di + sizeGroupDiameters[i+1]);
            }

            dLo = max(dLo, scalar(0));          // no negative volume

            sizeGroupVolumeWidths[i] =
                piBy6*(Foam::pow3(dHi) - Foam::pow3(dLo));

            Info<< "    " << sizeGroupNames[i]
                << ": dRange = [" << dLo << ", " << dHi << "] m"
                << ", deltaV = " << sizeGroupVolumeWidths[i] << " m^3" << endl;
        }
    }
}

// ========== CREATE TIME-AVERAGED FIELDS (only if size groups exist) ==========

if (nSizeClasses > 0)
{
    // ========== INSTANTANEOUS BREAKUP RATE FIELDS ==========
    // Written instantaneously and scaled by the daughter size-group volume
    // width deltaV_i, so dimensions are inv(dimTime) [1/s].

    binaryBreakupRate_Inst_.setSize(nSizeClasses);

    forAll(binaryBreakupRate_Inst_, i)
    {
        word fieldName = "binaryBreakupRate_" + sizeGroupNames[i] + "_Inst";

        binaryBreakupRate_Inst_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar(fieldName, inv(dimTime), 0.0)
            )
        );
    }

    // ========== INSTANTANEOUS COALESCENCE RATE FIELDS ==========

    coalescenceRate_Inst_.setSize(nSizeClasses);

    forAll(coalescenceRate_Inst_, i)
    {
        word fieldName = "coalescenceRate_" + sizeGroupNames[i] + "_Inst";

        coalescenceRate_Inst_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar(fieldName, dimVolume/dimTime, 0.0)
            )
        );
    }
    
    // ========== TIME-AVERAGED SIZE FRACTION FIELDS ==========
    // Use dimensions read from actual mesh fields (fFieldDimensions)
    
    f_Mean_.setSize(nSizeClasses);
    
    forAll(f_Mean_, i)
    {
        word fieldName = sizeGroupNames[i] + "_Mean_";
        
        f_Mean_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar(fieldName, fFieldDimensions, 0.0)
            )
        );
    }
    
    Info << "Created " << nSizeClasses << " instantaneous breakup rate fields (scaled by deltaV)" << endl;
    Info << "Created " << nSizeClasses << " instantaneous coalescence rate fields" << endl;
    Info << "Created " << nSizeClasses << " time-averaged size fraction fields" << endl;
}

// ========== TOTAL RATE FIELDS (always created) ==========

// Total breakup rate (time-averaged)
volScalarField binaryBreakupRateTotal_Mean_
(
    IOobject
    (
        "binaryBreakupRateTotal_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("binaryBreakupRateTotal_Mean_", inv(dimVolume*dimTime), 0.0)
);

// Total coalescence rate (time-averaged)
volScalarField coalescenceRateTotal_Mean_
(
    IOobject
    (
        "coalescenceRateTotal_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("coalescenceRateTotal_Mean_", dimVolume/dimTime, 0.0)
);

// Total breakup rate (instantaneous, sum of deltaV-scaled per-class -> 1/s)
volScalarField binaryBreakupRateTotal_Inst
(
    IOobject
    (
        "binaryBreakupRateTotal_Inst",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("binaryBreakupRateTotal_Inst", inv(dimTime), 0.0)
);

// Total coalescence rate (instantaneous)
volScalarField coalescenceRateTotal_Inst
(
    IOobject
    (
        "coalescenceRateTotal_Inst",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("coalescenceRateTotal_Inst", dimVolume/dimTime, 0.0)
);
