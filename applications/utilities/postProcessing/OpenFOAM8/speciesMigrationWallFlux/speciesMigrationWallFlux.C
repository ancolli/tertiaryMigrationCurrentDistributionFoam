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

\*---------------------------------------------------------------------------*/

#include "speciesMigrationWallFlux.H"

#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(speciesMigrationWallFlux, 0);
    addToRunTimeSelectionTable(functionObject, speciesMigrationWallFlux, dictionary);
}
}

void Foam::functionObjects::speciesMigrationWallFlux::writeFileHeader(const label i)
{
    // Add headers to output data
    writeHeader(file(), "Wall flux");
    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
    writeTabbed(file(), "min");
    writeTabbed(file(), "max");
    writeTabbed(file(), "integral");
    file() << endl;
}

void Foam::functionObjects::speciesMigrationWallFlux::calcFlux //void
(
    const volScalarField& fi_,
    const volScalarField& C_, 
    volScalarField& speciesMigrationWallFlux
)

{
    
    surfaceScalarField flux
    (
        nu_e_/nu_*F_*D_*fvc::snGrad(C_) + nu_e_/nu_*pow(F_,2)*u_*z_*fvc::interpolate(C_)*fvc::snGrad(fi_) // gradients are defined inwards
    );

    volScalarField::Boundary& speciesMigrationWallFluxBf =
        speciesMigrationWallFlux.boundaryFieldRef();

    const surfaceScalarField::Boundary& fluxBf =
        flux.boundaryField();

    forAll(speciesMigrationWallFluxBf, patchi)
    {
        speciesMigrationWallFluxBf[patchi] = fluxBf[patchi];
    }
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::functionObjects::speciesMigrationWallFlux::speciesMigrationWallFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),

    logFiles(obr_, name),
    writeLocalObjects(obr_, log),
    patchSet_(),
    fvOptions_(mesh_),
    
    fi_
    (
       IOobject
       (
           "fi",
           mesh_.time().timeName(),
           mesh_,
           IOobject::NO_READ,//MUST_READ,
           IOobject::NO_WRITE
       ),
       mesh_
       //dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), 0.0)
      
    ),

    C_
    (
       IOobject
       (
           type(),//"C_" + CName_,
           mesh_.time().timeName(),
           mesh_,
           IOobject::NO_READ,//MUST_READ,
           IOobject::NO_WRITE
       ),
       mesh_,
       dimensionedScalar("0", C_.dimensions(), 0.0)//dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), 0.0)//
    )
       
{
    volScalarField* speciesMigrationWallFluxPtr
    (
	new volScalarField
        (
            IOobject
            (
                type(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("0", C_.dimensions()/dimTime*dimLength, 0.0)
        )
    );


    
    mesh_.objectRegistry::store(speciesMigrationWallFluxPtr);
    
    read(dict);
    resetName(typeName);
    resetLocalObjectName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::speciesMigrationWallFlux::~speciesMigrationWallFlux()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::speciesMigrationWallFlux::read(const dictionary& dict)
{


    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);
     
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("patches", wordReList()))
        );

    Info<< type() << " " << name() << ":" << nl;

    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing wall patches: " << nl;
        labelHashSet filteredPatchSet;
        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                filteredPatchSet.insert(patchi);
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                WarningInFunction
                    << "Requested wall flux on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }


    if (dict.found("fvOptions"))
    {
        fvOptions_.reset(dict.subDict("fvOptions"));
    }

    
    dict.lookup("CName") >> CName_; 

        IOdictionary speciesProperties
    (
    	IOobject
    	(
        	"speciesProperties",
        	mesh_.time().constant(),//runTime.constant(),
        	mesh_,
        	IOobject::MUST_READ,//MUST_READ_IF_MODIFIED,
        	IOobject::NO_WRITE
    	)
    );


    dictionary speciesName
    (
    	speciesProperties.subDict(CName_)
    );

    speciesName.lookup("D") >> D_;
    speciesName.lookup("u") >> u_;
    speciesName.lookup("z") >> z_;

    IOdictionary transportProperties
    (
    	IOobject
    	(
        	"transportProperties",
		mesh_.time().constant(),//mesh_.time().timeName(),
        	mesh_,
        	IOobject::MUST_READ,
        	IOobject::NO_WRITE
   	)
   );

   transportProperties.lookup("F") >> F_;

   transportProperties.lookup("nu_e") >> nu_e_;

   transportProperties.lookup("nu_" + CName_) >> nu_;

   return true;
}


bool Foam::functionObjects::speciesMigrationWallFlux::execute()
{
    volScalarField& speciesMigrationWallFlux = lookupObjectRef<volScalarField>(type());

    volScalarField fi_
    (
       IOobject
       (
          "fi",
           mesh_.time().timeName(),
           mesh_,
           IOobject::MUST_READ,
           IOobject::NO_WRITE
       ),
       mesh_   
    );

    volScalarField C_
    (
       IOobject
       (
           "C_" + CName_,
           mesh_.time().timeName(),
           mesh_,
           IOobject::MUST_READ,
           IOobject::NO_WRITE
       ),
       mesh_ 
    );

    calcFlux(fi_, C_, speciesMigrationWallFlux);

    return true;
}


bool Foam::functionObjects::speciesMigrationWallFlux::end()
{

    return true;
}


bool Foam::functionObjects::speciesMigrationWallFlux::write()
{
    Log << type() << " " << name() << " write:" << nl;

    writeLocalObjects::write();

    logFiles::write();

    const volScalarField& speciesMigrationWallFlux =
        obr_.lookupObject<volScalarField>(type());

    const fvPatchList& patches = mesh_.boundary();

    const surfaceScalarField::Boundary& magSf =
        mesh_.magSf().boundaryField();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& pp = patches[patchi];

        const scalarField& hfp = speciesMigrationWallFlux.boundaryField()[patchi];

        const scalar minHfp = gMin(hfp);
        const scalar maxHfp = gMax(hfp);
        const scalar integralHfp = gSum(magSf[patchi]*hfp)/gSum(magSf[patchi]);
	const scalar current = gSum((magSf[patchi]*hfp)); //absolute value

        if (Pstream::master())
        {
            file()
                << mesh_.time().value()
                << tab << pp.name()
                << tab << minHfp
                << tab << maxHfp
                << tab << integralHfp
		<< tab << current
                << endl;
        }

        Log << "    min/max/average/current(" << pp.name() << ") = "
            << minHfp << ", " << maxHfp << ", " << integralHfp << ", " << current << endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
