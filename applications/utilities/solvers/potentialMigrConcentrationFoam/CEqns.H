{
        // concentration fields
	dictionary speciesSolveDict = mesh.solutionDict().subDict("solvers").subDict("C_species");
	scalar relax_init = readScalar(mesh.solutionDict().subDict("relaxationFactors").subDict("equations").lookup("relax_init"));
	scalar relax_end = readScalar(mesh.solutionDict().subDict("relaxationFactors").subDict("equations").lookup("relax_end"));
	scalar switch_iter = readScalar(mesh.solutionDict().subDict("relaxationFactors").subDict("equations").lookup("switch_iter"));
	dimensionedScalar MinSpecies("MinSpecies",dimensionSet(0,-3,0,0,1,0,0), SMALL);
	
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
	    speciesList[i] == max(MINIMO, speciesList[i]*(1-F*F*speciesList[i].z()*speciesList[i].u()*electroneut/keff));
	    speciesList[i].correctBoundaryConditions();
	    speciesList[i] == max(MINIMO, speciesList[i]);
	
	    Info<< "Min/max C_" << speciesNameList[i] << ":" << min(speciesList[i]).value() << ' '
            << max(speciesList[i]).value() << endl;	  
	}
}

// ************************************************************************* //
