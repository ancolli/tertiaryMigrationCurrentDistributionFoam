// * * * * * * * * * * * * * Properties fields * * * * * * * * * * * * * * * //

Info<< "Reading speciesProperties\n" << endl;

IOdictionary speciesProperties
(
    IOobject
    (
        "speciesProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const wordList speciesNameList(speciesProperties.toc());

PtrList<species> speciesList(speciesNameList.size());

forAll(speciesNameList, i)
{
    const word& speciesName = speciesNameList[i];
    const dictionary& speciesDict = speciesProperties.subDict(speciesName);

    speciesList.set
    (
        i,
        new species
        (
            speciesDict,
            mesh
        )
    );

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);



Info<< "Reading turbulent Schmidt number\n" << endl;

dimensionedScalar Sct
(
    transportProperties.lookup("Sct")
);


Info<< "Reading electric conductivity kf\n" << endl;

dimensionedScalar kf
(
    transportProperties.lookup("kf")
);

Info<< "Reading Faraday constant F\n" << endl;

dimensionedScalar F
(
    transportProperties.lookup("F")
);


Info<< "Reading controlProperties\n" << endl;

IOdictionary controlProperties
(
    IOobject
    (
        "controlProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

IOdictionary Results
(
    IOobject
    (
        "Results",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::AUTO_WRITE
    )
);

// * * * * * * * * * * * * * Reading fields * * * * * * * * * * * * * * * * //

Info<< "Reading field U\n" << endl;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field nut\n" << endl;

volScalarField nut //nutPointer = new volScalarField//
(
    IOobject
    (
        "nut",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,//MUST_READ,//
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(0,2,-1,0,0,0,0), 0.0)
);

Info<< "Reading field fi\n" << endl;

    volScalarField fi
    (
        IOobject
        (
            "fi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField Flux
    (
        IOobject
        (
            "Flux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE//AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("0", dimensionSet(0,-2,0,0,0,1,0), 0.0)
    );

    surfaceScalarField SumGrad_C
    (
        IOobject
        (
            "SumGrad_C",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE//NO_WRITE//
        ),
        mesh,
	dimensionedScalar("0", dimensionSet(0,-2,0,0,0,1,0), 0.0)//(0,-2,-1,0,1,0,0)
    );

    volScalarField keff
    (
        IOobject
        (
            "keff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,//MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("keff", kf)
    );

    volScalarField electroneut
    (
        IOobject
        (
            "electroneut",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,//MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0,-3,0,0,1,0,0), 0.0)
    );

    surfaceScalarField electroCon
    (
        IOobject
        (
            "electroCon",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,//MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(-1,-5,4,0,0,2,0), 0.0)
    );

    volScalarField diffterm
    (
        IOobject
        (
            "diffterm",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,//MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0,-3,0,0,0,1,0), 0.0)
    );

#include "createPhi.H"

#include "createFvOptions.H"
