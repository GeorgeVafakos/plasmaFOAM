/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) YEAR AUTHOR, AFFILIATION
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

Description
    Template for use with codeStream.

\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "Ostream.H"
#include "Pstream.H"
#include "pointField.H"
#include "tensor.H"
#include "unitConversion.H"

//{{{ begin codeInclude
#line 22 "/home/george/OpenFOAM/george-v2406/plasmaFOAM/actuator/AC-DBD/testGeom/test/0/rhoq/#codeStream"
#include "fvCFD.H"
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C" void codeStream_4346bbee0056a638bd3166d202091b38dd9fe5c2(Foam::Ostream& os, const Foam::dictionary& dict)
{
//{{{ begin code
    #line 39 "/home/george/OpenFOAM/george-v2406/plasmaFOAM/actuator/AC-DBD/testGeom/test/0/rhoq/#codeStream"
const IOdictionary& d = static_cast<const IOdictionary&>(dict);
        const fvMesh& mesh = refCast<const fvMesh>(d.db());
        // Access internal mesh information
        
        scalarField rhoq(mesh.nCells(), 0.);

        scalar rhoq0 = -0.00758;
        scalar PlasmaStart = 50e-6;
        scalar PlasmaWidth = 9.0e-3;
        scalar PlasmaThick = 2.4e-3;
        
        forAll(rhoq, i)
        {
            //Access cell centers coordinates
            const scalar x = mesh.C()[i][0];
            const scalar y = mesh.C()[i][1];
            // const scalar z = mesh.C()[i][2];
            
            // Triangular plasma region - linear distribution in the x-y plane
            if ( x>PlasmaStart && x<PlasmaWidth && y>0 && y<((PlasmaThick/(PlasmaWidth-PlasmaStart))*(PlasmaWidth-x)) )
            {
                rhoq[i] = rhoq0*(1.0 - (x-PlasmaStart)/PlasmaWidth - y/PlasmaThick);
            }
            /*
            // Triangular plasma region - linear distribution in the x direction
            if ( x>PlasmaStart && x<PlasmaWidth && y>0 && y<((PlasmaThick/(PlasmaWidth-PlasmaStart))*(PlasmaWidth-x)) )
            {
                rhoq[i] = rhoq0*((PlasmaWidth-(x-PlasmaStart))/PlasmaWidth);
            }
            
            // Elliptical plasma region - elliptical distribution in the x-y plane
            if ( x>PlasmaStart && x<PlasmaWidth && y>0 && y<sqrt(pow(PlasmaThick/(PlasmaWidth-PlasmaStart),2.0)*(pow(PlasmaWidth,2.0)-pow(x,2.0)) ))
            {
                rhoq[i] = rhoq0*(1.0 - (x-PlasmaStart)/PlasmaWidth - y/PlasmaThick);
            }
            
            // Rectangular plasma region - uniform distribution
            if ( x>PlasmaStart && x<PlasmaWidth && y>0 && y<PlasmaThick  )
            {
                rhoq[i] = rhoq0;
            }
            */
        }

        // writeEntry(os, "", rhoq);
        // rhoq.write();
        rhoq.writeEntry("", os);
//}}} end code
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

