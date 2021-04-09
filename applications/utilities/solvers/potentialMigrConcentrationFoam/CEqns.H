{
        // concentration fields
	dictionary speciesSolveDict = mesh.solutionDict().subDict("solvers").subDict("C_species");
	scalar relax_init = readScalar(mesh.solutionDict().subDict("relaxationFactors").subDict("equations").lookup("relax_init"));
	scalar relax_end = readScalar(mesh.solutionDict().subDict("relaxationFactors").subDict("equations").lookup("relax_end"));
	scalar switch_iter = readScalar(mesh.solutionDict().subDict("relaxationFactors").subDict("equations").lookup("switch_iter"));
	//dictionary relaxation = mesh.solutionDict().subDict("relaxationFactors").subDict("equations");
	//dictionary relaxation_end = mesh.solutionDict().subDict("relaxationFactors").subDict("equations");
	dimensionedScalar MinSpecies("MinSpecies",dimensionSet(0,-3,0,0,1,0,0), SMALL);
	//dimensionedScalar MaxSpecies("MaxSpecies",dimensionSet(0,-3,0,0,1,0,0), 3000.0);

	
	//Results.set("iter", ++iter); 

	//turbulent diffusivity	
	volScalarField Dturb
	(
		"Dturb",
		nut/Sct
	);
		
	//solve for each species
        forAll(speciesNameList, i)
        {

            fvScalarMatrix MassTransport
            (
                  fvm::ddt(speciesList[i])
                + fvm::div(phi, speciesList[i], "div(phi,C_species)")
                - fvm::laplacian(speciesList[i].D() + Dturb, speciesList[i], "laplacian(D,C_species)") //     
		- F*speciesList[i].z()*fvm::SuSp(speciesList[i].u()*fvc::laplacian(fi),speciesList[i])    
		==
		  F*speciesList[i].z()*speciesList[i].u()*(fvc::grad(speciesList[i]) & fvc::grad(fi))
		+ fvOptions(speciesList[i])
	    );
		
	    if (runTime.value() < switch_iter)
	    {
		MassTransport.relax(relax_init);
	    }
	    else
	    {
		MassTransport.relax(relax_end);
	    }
	    
	    fvOptions.constrain(MassTransport);
	    MassTransport.solve(speciesSolveDict);
	    fvOptions.correct(speciesList[i]);
	    
	// impose electroneutrality
	    dimensionedScalar MINIMO("MINIMO",dimensionSet(0,-3,0,0,1,0,0), SMALL);
	    //dimensionedScalar MAXIMO("MAXIMO",dimensionSet(0,0,0,0,0,0,0), 0.99);
	    speciesList[i] == max(MINIMO, speciesList[i]*(1-F*F*speciesList[i].z()*speciesList[i].u()*electroneut/keff));
	    //bound(speciesList[i], MinSpecies);
	    speciesList[i].correctBoundaryConditions();
	    //bound(speciesList[i], MinSpecies);
	    speciesList[i] == max(MINIMO, speciesList[i]);

	    /*forAll(mesh.boundaryMesh(), p_ith)
    	    { 
		
			speciesList[i].boundaryFieldRef()[p_ith] == speciesList[i].boundaryFieldRef()[p_ith]*(1-F.value()*F.value()*speciesList[i].z().value()*speciesList[i].u().value()*electroneut.boundaryFieldRef()[p_ith]/keff.boundaryFieldRef()[p_ith]);//
            } */

	    //
	
	    Info<< "Min/max C_" << speciesNameList[i] << ":" << min(speciesList[i]).value() << ' '
            << max(speciesList[i]).value() << endl;

	    //Info<< "Min/max multiplicador_" << speciesNameList[i] << ":" << min((1-speciesList[i].z()*speciesList[i].D()*electroneut/electroconductivity)).value() << ' '
            //<< max((1-speciesList[i].z()*speciesList[i].D()*electroneut/electroconductivity)).value() << endl;
	  
	}
}

// ************************************************************************* //
