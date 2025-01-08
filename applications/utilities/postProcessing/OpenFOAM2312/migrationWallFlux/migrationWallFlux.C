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

#include "migrationWallFlux.H"

#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(migrationWallFlux, 0);
    addToRunTimeSelectionTable(functionObject, migrationWallFlux, dictionary);
}
}

void Foam::functionObjects::migrationWallFlux::writeFileHeader(Ostream& os) const//void Foam::functionObjects::migrationWallFlux::writeFileHeader(const label i)
{
    // Add headers to output data
    writeHeader(os, "Wall flux");//writeHeader(file(), "Wall flux");
    writeCommented(os, "Time");//writeCommented(file(), "Time");
    writeTabbed(os, "patch");//writeTabbed(file(), "patch");
    writeTabbed(os, "min current density (A/m2)");//writeTabbed(file(), "min current density (A/m2)");
    writeTabbed(os, "max current density (A/m2)");//writeTabbed(file(), "max current density (A/m2)");
    writeTabbed(os, "Average current density (A/m2)");//writeTabbed(file(), "Average current density (A/m2)");
    writeTabbed(os, "Absolute current (A)");//writeTabbed(file(), "Absolute current (A)");
    os << endl;//file() << endl;
}

void Foam::functionObjects::migrationWallFlux::calcFlux //void
(
    const volScalarField& kf_,
    const volScalarField& fi_, //const
    const surfaceScalarField& SumGrad_C_,
    volScalarField& migrationWallFlux
)

{
    
    surfaceScalarField flux
    (
      fvc::interpolate(kf_)*fvc::snGrad(fi_) - SumGrad_C_ // gradients are defined inwards
    );

    volScalarField::Boundary& migrationWallFluxBf =
        migrationWallFlux.boundaryFieldRef();

    const surfaceScalarField::Boundary& fluxBf =
        flux.boundaryField();

    forAll(migrationWallFluxBf, patchi)
    {
        migrationWallFluxBf[patchi] = fluxBf[patchi];
    }
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::functionObjects::migrationWallFlux::migrationWallFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),

    //logFiles(obr_, name),
    writeFile(obr_, name, typeName, dict),//writeLocalObjects(obr_, log),
    patchSet_()//,
    //fvOptions_(mesh_),

    /*kf_
    (
       IOobject
       (
           "keff",
           mesh_.time().timeName(),
           mesh_,
           IOobject::MUST_READ,
           IOobject::NO_WRITE
       ),
       mesh_
    ),
 
    fi_
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
      
    )

    SumGrad_C_
    (
       IOobject
       (
           "SumGrad_C",
           mesh_.time().timeName(),
           mesh_,
           IOobject::MUST_READ,
           IOobject::NO_WRITE
       ),
       mesh_
    )*/
 
{
    volScalarField* migrationWallFluxPtr
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
            dimensionedScalar("0", dimCurrent/dimArea, 0.0)
        )
    );


    
    mesh_.objectRegistry::store(migrationWallFluxPtr);
    
    read(dict);
    //resetName(typeName);
    //resetLocalObjectName(typeName);
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::migrationWallFlux::~migrationWallFlux()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::migrationWallFlux::read(const dictionary& dict)
{


    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);//writeLocalObjects::read(dict);
     
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


    /*if (dict.found("fvOptions"))
    {
        fvOptions_.reset(dict.subDict("fvOptions"));
    }*/


    return true;
}


bool Foam::functionObjects::migrationWallFlux::execute()
{
    volScalarField& migrationWallFlux = lookupObjectRef<volScalarField>(type());
    
    volScalarField kf_
    (
       IOobject
       (
           "keff",
           mesh_.time().timeName(),
           mesh_,
           IOobject::MUST_READ,
           IOobject::NO_WRITE
       ),
       mesh_
    );
    
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
    
    surfaceScalarField SumGrad_C_
    (
       IOobject
       (
           "SumGrad_C",
           mesh_.time().timeName(),
           mesh_,
           IOobject::MUST_READ,
           IOobject::NO_WRITE
       ),
       mesh_
    );

    calcFlux(kf_, fi_, SumGrad_C_, migrationWallFlux);

    return true;
}


/*bool Foam::functionObjects::migrationWallFlux::end()
{

    return true;
}*/


bool Foam::functionObjects::migrationWallFlux::write()
{
    //Log << type() << " " << name() << " write:" << nl;

    //writeLocalObjects::write();

    //logFiles::write();

    //const volScalarField& migrationWallFlux =
    //    obr_.lookupObject<volScalarField>(type());
    
    const volScalarField& migrationWallFlux = lookupObject<volScalarField>(type());

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << migrationWallFlux.name() << endl;

    migrationWallFlux.write();

    const fvPatchList& patches = mesh_.boundary();

    const surfaceScalarField::Boundary& magSf =
        mesh_.magSf().boundaryField();

    for (const label patchi : patchSet_)//forAllConstIter(labelHashSet, patchSet_, iter)
    {
        //label patchi = iter.key();
        const fvPatch& pp = patches[patchi];

        const scalarField& hfp = migrationWallFlux.boundaryField()[patchi];

        const scalar minHfp = gMin(hfp);
        const scalar maxHfp = gMax(hfp);
        const scalar integralHfp = gSum(magSf[patchi]*hfp)/gSum(magSf[patchi]);
	const scalar currentAbs = gSum(mag(magSf[patchi]*hfp)); //absolute value

        if (Pstream::master())
        {
            writeCurrentTime(file());
        
            file()
                //<< mesh_.time().value()
                << tab << pp.name()
                << tab << minHfp
                << tab << maxHfp
                << tab << integralHfp
		<< tab << currentAbs
                << endl;
	}
        Log << "    min current density (A/m2)/max current density (A/m2)/Average current density (A/m2)/ Absolute current (A) (" << pp.name() << ") = "
            << minHfp << ", " << maxHfp << ", " << integralHfp << ", " << currentAbs << endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
