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
    - binaryBreakupRate_<sizeGroupName>_Mean_  (e.g., binaryBreakupRate_f0.air_Mean_)
    - coalescenceRate_<sizeGroupName>_Mean_
    - binaryBreakupRateTotal_Mean_
    - coalescenceRateTotal_Mean_
    - <sizeGroupName>_Mean_                    (e.g., f0.air_Mean_)

\*---------------------------------------------------------------------------*/

// ========== STORAGE FOR SIZE GROUP INFORMATION ==========

label nSizeClasses = 0;
List<word> sizeGroupNames;
List<scalar> sizeGroupDiameters;

// Store the dimension of f fields (will be read from actual mesh fields)
dimensionSet fFieldDimensions(dimless);

// Flag to track if pop balance rate fields have been found
bool popBalRateFieldsFound = false;

// Declare PtrLists at outer scope (will be populated if size groups exist)
PtrList<volScalarField> binaryBreakupRate_Mean_;
PtrList<volScalarField> coalescenceRate_Mean_;
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
    }
}

// ========== CREATE TIME-AVERAGED FIELDS (only if size groups exist) ==========

if (nSizeClasses > 0)
{
    // ========== TIME-AVERAGED BREAKUP RATE FIELDS ==========
    
    binaryBreakupRate_Mean_.setSize(nSizeClasses);
    
    forAll(binaryBreakupRate_Mean_, i)
    {
        word fieldName = "binaryBreakupRate_" + sizeGroupNames[i] + "_Mean_";
        
        binaryBreakupRate_Mean_.set
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
                dimensionedScalar(fieldName, inv(dimVolume*dimTime), 0.0)
            )
        );
    }
    
    // ========== TIME-AVERAGED COALESCENCE RATE FIELDS ==========
    
    coalescenceRate_Mean_.setSize(nSizeClasses);
    
    forAll(coalescenceRate_Mean_, i)
    {
        word fieldName = "coalescenceRate_" + sizeGroupNames[i] + "_Mean_";
        
        coalescenceRate_Mean_.set
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
    
    Info << "Created " << nSizeClasses << " time-averaged breakup rate fields" << endl;
    Info << "Created " << nSizeClasses << " time-averaged coalescence rate fields" << endl;
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
