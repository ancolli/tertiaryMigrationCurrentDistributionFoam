{
//diffusion term
    volScalarField difterm(0*F*fvc::laplacian(speciesList[0].D(), speciesList[0], "laplacian(D,C_species)"));
    //volScalarField divterm(0*F*fvc::div(phi, speciesList[0].z()*speciesList[0], "div(phi,C_species)"));
    volScalarField kff(0*F*F*speciesList[0].z()*speciesList[0].z()*speciesList[0].u()*speciesList[0]);
    volScalarField carga(0*speciesList[0].z()*speciesList[0]);
    //volScalarField electroCon(0*speciesList[0].z()*speciesList[0].z()*speciesList[0].D()*speciesList[0]);
    
    forAll(speciesNameList, i)
    {
	carga += speciesList[i].z()*speciesList[i];
	//electroCon += speciesList[i].z()*speciesList[i].z()*speciesList[i].D()*speciesList[i];
	difterm += F*speciesList[i].z()*fvc::laplacian(speciesList[i].D(), speciesList[i], "laplacian(D,C_species)");
	//divterm += F*speciesList[i].z()*fvc::div(phi, speciesList[i], "div(phi,C_species)");
       	kff += F*F*speciesList[i].z()*speciesList[i].z()*speciesList[i].u()*speciesList[i];
	//electroCon += F*F*speciesList[i].z()*speciesList[i].z()*speciesList[i].u()*speciesList[i]/speciesList[i].D();
    };

//conductivity field  
    keff == kff;
    electroneut == carga;
    //electroconductivity == electroCon;
    diffterm == difterm;
    
    Info<< "Min/max keff:" << min(keff).value() << ' ' << max(keff).value() << endl;
    Info<< "Max electroneutrality:" << max(mag(electroneut)).value() << endl;//<< max(mag(electroneut.internalField())).value() << endl;
    //Info<< "Min/max difterm:" << min(difterm).value() << ' ' << max(difterm).value() << endl;
//potential field
    fi.storePrevIter();
    solve //fvScalarMatrix fi_eq//
    ( 
       fvm::laplacian(keff, fi) == - diffterm  
    );
    
    fi.relax();

    //fi_eq.relax();
    //fi_eq.solve();

scalar maxPot = readScalar(controlProperties.lookup("maxPot"));
scalar minPot = readScalar(controlProperties.lookup("minPot"));

dimensionedScalar Min_fi("Min_fi",dimensionSet(1,2,-3,0,0,-1,0), SMALL);

scalar Vcell = readScalar(Results.lookup("Vcell"));
dimensionedScalar Max_fi("Max_fi",dimensionSet(1,2,-3,0,0,-1,0), Vcell);

fi == max(min(fi,Max_fi),Min_fi);
fi.correctBoundaryConditions();
//fi == max(min(fi,Max_fi),Min_fi);

Info<< "Min/max fi:" << min(fi).value() << ' ' << max(fi).value() << endl;

//check if wall_anode exist
string patchNameNeeded = "wall_anode";
scalar patchOK = 0;

    forAll(mesh.boundaryMesh(), p_ith)
    { 
// checks if "wall_anode" is in the current patchName
      string patchName = mesh.boundaryMesh()[p_ith].name();	
      if (patchNameNeeded == patchName)  
      {
        patchOK = 1;
      }
    }  

if (patchOK==1)  // checks if "wall_anode" is in the current patchName
      {
// if wall_anode exists then: Calculates flux on the anode
        label patchi_A = mesh.boundaryMesh().findPatchID("wall_anode"); 
	label patchi_C = mesh.boundaryMesh().findPatchID("wall_cathode");  
        Flux.boundaryFieldRef()[patchi_A] = keff.boundaryFieldRef()[patchi_A]*fi.boundaryFieldRef()[patchi_A].snGrad() - SumGrad_C.boundaryFieldRef()[patchi_A]; //gradients are defined inwards
	Flux.boundaryFieldRef()[patchi_C] = keff.boundaryFieldRef()[patchi_C]*fi.boundaryFieldRef()[patchi_C].snGrad() - SumGrad_C.boundaryFieldRef()[patchi_C]; //gradients are defined inwards
// Calculates average Flux
	scalar sumField = 0.0;
	scalar sumField_C = 0.0;
        {
        	sumField = gSum
                (
                	mesh.magSf().boundaryField()[patchi_A]
                      * Flux.boundaryField()[patchi_A]
                );
		sumField_C = gSum
                (
                	mesh.magSf().boundaryField()[patchi_C]
                      * Flux.boundaryField()[patchi_C]
                );
        }

        Info<< "Total current anode = " << sumField << " A "<< endl;
	Info<< "Total current cathode = " << sumField_C << " A "<< endl;

// reading probeLocation (reference electrode)
        label cellNo = mesh.findCell(vector(controlProperties.lookup("refLocation")));
	scalar Vref = fi[cellNo];
// calculate cell potential
        scalar mult = 0.0;
        scalar Vcell_old = 0.0;
	scalar new_I = 0.0;	
        scalar var = 0.0;

        scalar sumField_old = readScalar(Results.lookup("I"));

        scalar control = readScalar(controlProperties.lookup("control"));  
        scalar I = readScalar(controlProperties.lookup("I"));
        
        scalar alpha = readScalar(controlProperties.lookup("alpha"));
        scalar tol_I = readScalar(controlProperties.lookup("tol_I"));
	scalar tol_jump = readScalar(controlProperties.lookup("tol_jump"));
// potentiostatic mode
        scalar maxvar = readScalar(controlProperties.lookup("maxvar"));
        scalar tol_ref = readScalar(controlProperties.lookup("tol_ref"));
	scalar ref = readScalar(controlProperties.lookup("ref"));

        var=(ref < 0.0) ? -ref-Vref : ref-(Vcell-Vref); // <0 cathodic overpotential, >0 anodic overpotential
        var=(var > 0.0) ? min(var,maxvar) : max(var,-maxvar);
        var=(mag(var/ref) > tol_ref) ? var : 0.0; // limito ref a x%
// galvanostatic
	new_I = (mag(sumField_old/sumField-1.0) < tol_jump) ? sumField : I; //damping variable
	mult = (mag(new_I/I-1.0) > tol_I) ? (I-sumField)*Vcell/sumField/3.0 : 0.0;
	mult = (mult > 0.0) ? min(mult,maxvar) : max(mult,-maxvar);
// overrelaxation 
        Vcell_old = Vcell;
	Vcell = (control==1) ? Vcell+var : ((control==2) ? maxPot : Vcell+mult);
	Vcell = Vcell_old+alpha*(Vcell-Vcell_old);
	Vcell = max(min(Vcell,maxPot), minPot); // limited potentials
    
// writing Vcell, I and Vref       
        Results.set("Vcell", Vcell);
        Results.set("I", sumField);
	Results.set("Vref", (ref < 0.0) ? -Vref : (Vcell-Vref));//

   }
   else
   {
	FatalErrorIn
        (
            " Foam::patchProbes::findElements(const fvMesh&)"
        )   << " No patch name: wall_anode found "
            << " the patch name: wall_anode should exist "
            << endl
            << exit(FatalError);     
   }



}

// ************************************************************************* //
