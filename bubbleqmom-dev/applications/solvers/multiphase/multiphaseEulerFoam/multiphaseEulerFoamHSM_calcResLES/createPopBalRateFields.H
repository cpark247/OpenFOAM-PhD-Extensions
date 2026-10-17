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

    All rate fields are FREQUENCIES [1/s], so that the time scales of
    Buffo et al. (2016), doi:10.1016/j.ces.2016.07.032, follow directly as the
    reciprocal:
        tau_breakup     = 1/binaryBreakupRateTotal_Inst
        tau_coalescence = 1/coalescenceRateTotal_Inst

    LuoSvendsenForSAS2 already emits breakup frequencies (its rate density is
    integrated over the daughter volume bin inside the model, using the same
    volume widths as populationBalanceModel::calcDeltas()).
    PrinceBlanchForSAS2 emits coalescence KERNELS [m^3/s]; updatePopBalRate-
    Averages.H converts them with the number density built from the Sauter
    mean diameter, n = alpha_d/((pi/6)*d_32^3), as given in the reference.

    Fields created:
    - binaryBreakupRate_<sizeGroupName>_Inst   (instantaneous, 1/s;
      breakup frequency of that PARENT class, e.g. binaryBreakupRate_f0.air_Inst)
    - coalescenceRate_<sizeGroupName>_Inst     (instantaneous, 1/s)
    - binaryBreakupRateTotal_Mean_             (time-averaged, 1/s)
    - coalescenceRateTotal_Mean_               (time-averaged, 1/s)
    - binaryBreakupRateTotal_Inst              (instantaneous, 1/s)
    - coalescenceRateTotal_Inst                (instantaneous, 1/s)
    - <sizeGroupName>_Mean_                     (time-averaged; e.g., f0.air_Mean_)

\*---------------------------------------------------------------------------*/

// ========== STORAGE FOR SIZE GROUP INFORMATION ==========

label nSizeClasses = 0;
List<word> sizeGroupNames;
List<scalar> sizeGroupDiameters;

// Name of the dispersed phase carrying the size groups (parsed from the size
// group field names, e.g. "f0.air" -> "air"). Used to look up alpha.<phase>
// for the Buffo et al. (2016) coalescence frequency conversion.
word dispersedPhaseName;

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

        // "f0.air" -> "air"
        dispersedPhaseName =
            sizeGroupNames[0].substr(sizeGroupNames[0].find(".") + 1);

        Info << "Dispersed phase carrying the size groups: "
             << dispersedPhaseName << endl;
        
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

    }
}

// ========== CREATE TIME-AVERAGED FIELDS (only if size groups exist) ==========

if (nSizeClasses > 0)
{
    // ========== INSTANTANEOUS BREAKUP RATE FIELDS ==========
    // Breakup frequency [1/s] of each PARENT size class, taken straight from
    // LuoSvendsenForSAS2 (no post-scaling needed).

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
    // Coalescence frequency [1/s]: the mean kernel of PrinceBlanchForSAS2
    // multiplied by the d_32-based number density (see
    // updatePopBalRateAverages.H).

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
                dimensionedScalar(fieldName, inv(dimTime), 0.0)
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
    
    Info << "Created " << nSizeClasses << " instantaneous breakup frequency fields [1/s]" << endl;
    Info << "Created " << nSizeClasses << " instantaneous coalescence frequency fields [1/s]" << endl;
    Info << "Created " << nSizeClasses << " time-averaged size fraction fields" << endl;
}

// ========== TOTAL RATE FIELDS (always created) ==========

// Total breakup frequency (time-averaged, 1/s)
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
    dimensionedScalar("binaryBreakupRateTotal_Mean_", inv(dimTime), 0.0)
);

// Total coalescence frequency (time-averaged, 1/s)
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
    dimensionedScalar("coalescenceRateTotal_Mean_", inv(dimTime), 0.0)
);

// Total breakup frequency (instantaneous, 1/s)
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

// Total coalescence frequency (instantaneous, 1/s)
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
    dimensionedScalar("coalescenceRateTotal_Inst", inv(dimTime), 0.0)
);
