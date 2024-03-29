/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

Class
    Foam::functionObjects::wallFlux

Group

Description
    This function object...

    Example of function object specification:
    \verbatim
    wallFlux1
    {
        type           wallFlux;
        libs ("libwallFluxFunctionObject.so");
        ...
        region      fluid;
        patches     (".*Wall");
    }
    \endverbatim

Usage
    \table
        Property | Description                | Required   | Default value
        type     | type name: wallHeatFlux    | yes        |
        patches  | list of patches to process | no         | all wall patches
        region   | region to be evaluated     | no         | default region
    \endtable

SourceFiles
    wallFlux.C

\*---------------------------------------------------------------------------*/

#ifndef migrationWallFlux_H
#define migrationWallFlux_H

#include "fvMeshFunctionObject.H"

#include "logFiles.H"
#include "writeLocalObjects.H"
#include "volFieldsFwd.H"
#include "HashSet.H"

//#include "volFields.H"
#include "fvOptionList.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class migrationWallFlux Declaration
\*---------------------------------------------------------------------------*/

class migrationWallFlux
:
    public fvMeshFunctionObject,
    public logFiles,
    public writeLocalObjects
{


protected:

    // Protected data

        //- Optional list of patches to process
        labelHashSet patchSet_;

    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);

        //- Calculate the flux
        void calcFlux
        (
            const volScalarField& kf_,
            const volScalarField& fi_, //const
	    const surfaceScalarField& SumGrad_C_,
            volScalarField& migrationWallFlux
        );

private:

//Private data

 //- Run-time selectable finite volume options, e.g. sources, constraints
        fv::optionList fvOptions_;

	volScalarField kf_;
	
        volScalarField fi_;

	//surfaceScalarField SumGrad_C_;

// Private Member Functions

        //- Disallow default bitwise copy construct
        migrationWallFlux(const migrationWallFlux&);

        //- Disallow default bitwise assignment
        void operator=(const migrationWallFlux&);


public:



    //- Runtime type information
    TypeName("migrationWallFlux");


    // Constructors

        //- Construct from Time and dictionary
        migrationWallFlux
        (
            const word& name,
            const Time& runTime,
	    
            const dictionary& dict
        );


    //- Destructor
    virtual ~migrationWallFlux();


    // Member Functions

        //- Read the migrationWallFlux data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Write the migrationWallFlux
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
