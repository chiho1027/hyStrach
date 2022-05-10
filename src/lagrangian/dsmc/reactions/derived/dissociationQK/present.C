	///////////////////////////start////////////////////////////////////////////////////		
		
	//- The collision energy is being subtracted the heat of reaction
        scalar collisionEnergy = translationalEnergy + EVibP_disso
	  // - heatOfReactionDissociationJoules_[nReac][vibModeDisso];
	  - heatOfReactionDissociationJoules_[nReac];
	
	const label iMax = collisionEnergy/ktheta;
	label  j = 0;	  	  
	scalar func  = 0.0;

	if(iMax > 0)
	{  
	  do // acceptance - rejection
	  {
	    j = cloud_.randomLabel(0, iMax);
	      
	    func =
	      pow(1.0 - j*ktheta/collisionEnergy ,1.5-omegaPQ );
		
	  }while(func < cloud_.rndGen().sample01<scalar>() );	    
	}

	collisionEnergy -= j*ktheta;
	
	const scalar postExcitVibE = j*ktheta;
	
	/* //not good the vibrational distributtion will not become boltzman
        //- Energy redistribution for particle Q
        cloud_.binaryCollision().redistribute
        (
            q, collisionEnergy, omegaPQ, true
        );       
	*/
	
	const scalar& dofP = cloud_.constProps(typeIdP).rotationalDegreesOfFreedom();
        const scalar mP = cloud_.constProps(typeIdP).mass();
        const scalar mQ = cloud_.constProps(typeIdQ).mass();
        const scalar mR = mP*mQ/(mP + mQ);
        scalar relVelNonDissoParticle = sqrt(2.0*collisionEnergy/mR);

	//temp
	//vector pU = p.U();
	//vector qU = q.U();
	
        cloud_.binaryCollision().postCollisionVelocities
        (
            typeIdP,
            typeIdQ,
            p.U(),//pU,
            q.U(),//qU,
            relVelNonDissoParticle
        );	
	
	//- Energy left for the 2 products
        const scalar ERotP = p.ERot();
	const scalar EVibP_tot =
	  cloud_.constProps(typeIdP).eVib_tot
	  (
	   p.vibLevel()
	   );
	
	//remain energy and dof
	const scalar TMacro = cloud_.fields().overallT(p.cell());
	scalar translationalEnergyLeft = ERotP + EVibP_tot - EVibP_disso;
	scalar remainDof  = dofP;
        forAll(p.vibLevel(), m)
	{
	  if(m != vibModeDisso)
	  {	    
	    const scalarList equiVibTheta(1,cloud_.constProps(typeIdP).thetaV()[m]);
	    remainDof    += 2.0*cloud_.temperatureFunction(TMacro, equiVibTheta);	    
	  }
	}
	  
	///////////////////////handle cleave partile internal energy///////////////	
        //- Post-reaction
        const label typeId1 = productIdsDissociation_[nReac][0];
        const label typeId2 = productIdsDissociation_[nReac][1];
   
        //- Mass of products 1 and 2
        const scalar mP1 = cloud_.constProps(typeId1).mass();
        const scalar mP2 = cloud_.constProps(typeId2).mass();
        const scalar mRproducts = mP1*mP2/(mP1 + mP2);
	
	const scalarList& thetaVProductMole   = cloud_.constProps(typeId1).thetaV();
	const scalarList& thetaVProductSecond = cloud_.constProps(typeId2).thetaV();
	const scalar& moleRDof   = cloud_.constProps(typeId1).rotationalDegreesOfFreedom();
	const scalar& secondRDof = cloud_.constProps(typeId2).rotationalDegreesOfFreedom();       
		
	labelList vibLevelMole(thetaVProductMole.size(), 0);
	labelList vibLevelSecond(thetaVProductSecond.size(), 0);
	const label numberVibMode = vibLevelMole.size() + vibLevelSecond.size();

	label mode = 0;
	scalarList theta    =  thetaVProductMole;//assume
	labelList* vibLevel = &vibLevelMole;//assume
	
	//forAll(totalThetaVProductList ,m)
	for(label m=0; m<numberVibMode; m++)
	{
	    mode = m;        
	    
	    if(mode < thetaVProductMole.size() )
	    {
	      theta    =  thetaVProductMole;
	      vibLevel =  &vibLevelMole;
	    }
	    else
	    {
	      mode    -=  thetaVProductMole.size();
	      theta    =  thetaVProductSecond;		
	      vibLevel =  &vibLevelSecond;
	    }

	    const scalarList thetaVMode(1, theta[mode]);
	    remainDof  -= 2.0*cloud_.temperatureFunction(TMacro, thetaVMode);
	    
	    cloud_.postReactionVibrationalRedistribution
	    (
	     mode,
	     remainDof,
	     theta,
	     vibLevel,
	     translationalEnergyLeft
	     );	    
	}        

	//temp
	/*
	if(vibLevelMole.size() > 0)
	{
	  Info << "12OHre = " << vibLevelMole[0]<< endl;
	}

	return;
	*/
	
	////////////////-first Trial L-B redistribution (rotation)
	scalar ERotProductMole   = 0.0;
	if( thetaVProductMole.size() > 0 )
	{
	  remainDof -= moleRDof; 
	  const scalar energyRatioMole = cloud_.postCollisionRotationalEnergy( moleRDof, remainDof/2.0 );	  
	  ERotProductMole   = energyRatioMole*translationalEnergyLeft;
	  translationalEnergyLeft -= ERotProductMole;//- Relative collisionEnergy energy after rotational energy redistribution	  
	}

	////////////////-second Trial L-B redistribution (rotation)
	scalar ERotProductSecond = 0.0;
	if( thetaVProductSecond.size() > 0 )
	{
	  remainDof -= secondRDof; 
	  const scalar energyRatioSecond = cloud_.postCollisionRotationalEnergy( secondRDof, remainDof/2.0 );
	  ERotProductSecond = energyRatioSecond*translationalEnergyLeft;
	  translationalEnergyLeft -= ERotProductSecond;//- Relative collisionEnergy energy after rotational energy redistribution
	}
	
	///////////////////////////////////////////
	translationalEnergyLeft += postExcitVibE;
	
	///////////////////////////////end/////////////////////////////////////////
