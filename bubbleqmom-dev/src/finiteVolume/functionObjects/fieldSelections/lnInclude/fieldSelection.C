/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2019 OpenCFD Ltd.
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

#include "fieldSelection.H"
#include "objectRegistry.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldSelection::fieldSelection
(
    const objectRegistry& obr,
    const bool includeComponents
)
:
    List<fieldInfo>(),
    obr_(obr),
    includeComponents_(includeComponents),
    selection_()
{}


bool Foam::functionObjects::fieldSelection::resetFieldFilters
(
    const HashSet<wordRe>& names
)
{
    static word cmptStr(".component(");
    static string::size_type len(cmptStr.size());

    DynamicList<fieldInfo> nameAndComponent(names.size());

    forAllConstIter (HashSet<wordRe>,names,nameIter)
    {
        string::size_type n = nameIter.key().find(cmptStr);
        if (n != string::npos)
        {
            // Field should be written <field>.component(i)

            if (!includeComponents_)
            {
                FatalErrorInFunction
                    << "Component specification not allowed for " << nameIter.key()
                    << exit(FatalError);
            }

            if (nameIter.key().isPattern())
            {
                FatalErrorInFunction
                    << "Cannot use \".component option\" in combination with "
                    << "wildcards for " << nameIter.key()
                    << exit(FatalError);
            }

            word baseName = nameIter.key().substr(0, n);

            // Extract the component - number between ()'s
            string::size_type closei = nameIter.key().find(')', n);

            if (closei == string::npos)
            {
                FatalErrorInFunction
                    << "Invalid field component specification for "
                    << nameIter.key() << nl
                    << ". Field should be expressed as <field>.component(i)"
                    << exit(FatalError);
            }

            string::size_type cmptWidth = closei - n - len;

            label component
            (
                readLabel(IStringStream(nameIter.key().substr(n+len, cmptWidth))())
            );

            nameAndComponent.append(fieldInfo(wordRe(baseName), component));
        }
        else
        {
            nameAndComponent.append(fieldInfo(nameIter.key()));
        }
    }

    this->transfer(nameAndComponent);

    return true;
}


bool Foam::functionObjects::fieldSelection::resetFieldFilters
(
    const wordRe& name
)
{
    HashSet<wordRe> tmp(0);
    tmp.insert(name);
    return resetFieldFilters(tmp);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldSelection::read(const dictionary& dict)
{
    HashSet<wordRe> fields(dict.lookup("fields"));

    return resetFieldFilters(fields);
}


bool Foam::functionObjects::fieldSelection::containsPattern() const
{
    forAll (*this,fi)
    {
        if ((*this)[fi].name().isPattern())
        {
            return true;
        }
    }

    return false;
}


void Foam::functionObjects::fieldSelection::clearSelection()
{
    selection_.clear();
}


bool Foam::functionObjects::fieldSelection::updateSelection()
{
    return false;
}


bool Foam::functionObjects::fieldSelection::checkSelection()
{
    bool ok = true;
    forAll (*this,fi)
    {
        if (!(*this)[fi].found())
        {
            WarningInFunction
                << "Field " << (*this)[fi].name() << " not found"
                << endl;

            ok = false;
        }
    }

    return ok;
}


// ************************************************************************* //
