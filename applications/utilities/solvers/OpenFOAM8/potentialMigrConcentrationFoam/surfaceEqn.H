{

    forAll(mesh.boundaryMesh(), p_ith)
    { 
	// Perform calculation of SumGrad_C and electroCon   
	SumGrad_C.boundaryFieldRef()[p_ith] = 0.0*F.value()*speciesList[0].z().value()*speciesList[0].D().value()*speciesList[0].boundaryFieldRef()[p_ith].snGrad();
	electroCon.boundaryFieldRef()[p_ith] = 0.0*speciesList[0].boundaryFieldRef()[p_ith];
 	forAll(speciesNameList, i)
    	{
		SumGrad_C.boundaryFieldRef()[p_ith] -= F.value()*speciesList[i].z().value()*speciesList[i].D().value()*speciesList[i].boundaryFieldRef()[p_ith].snGrad(); //gradients are defined inwards
		electroCon.boundaryFieldRef()[p_ith] += F.value()*F.value()*speciesList[i].z().value()*speciesList[i].z().value()*speciesList[i].u().value()*speciesList[i].boundaryFieldRef()[p_ith]/speciesList[i].D().value();
	};
    } 

}

// ************************************************************************* //
