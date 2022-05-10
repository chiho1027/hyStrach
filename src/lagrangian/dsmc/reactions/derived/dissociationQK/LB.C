
	///////////////////////////start////////////////////////////////////////////////////		
	//- Energy left for the 2 products
        const scalar ERotP = p.ERot();
	const scalar EVibP_tot =
	  cloud_.constProps(typeIdP).eVib_tot
	  (
	   p.vibLevel()
	   );
	
	//remain energy
	scalar collisionEnergy =
	  ERotP + EVibP_tot + translationalEnergy
	  - heatOfReactionDissociationJoules_[nReac];
	
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
	
	const scalar reverseOmega =  0.5*(
					  cloud_.constProps(typeId1).omega()
					  + cloud_.constProps(typeId2).omega()
					 );
	
	const scalar moleSecondThirdomega =  0.5*(
						  reverseOmega
						  + cloud_.constProps(typeIdQ).omega()
						 );
	
	labelList vibLevelMole(thetaVProductMole.size(), 0);
	labelList vibLevelSecond(thetaVProductSecond.size(), 0);
	const label numberVibMode = vibLevelMole.size() + vibLevelSecond.size();

	//target
	scalar remainDOF  = moleRDof + secondRDof
	  + 2.0*(5.0 - reverseOmega - moleSecondThirdomega);

	label mode = 0;
	scalarList theta    =  thetaVProductMole;//assume
	labelList* vibLevel = &vibLevelMole;//assume
	
	//redistribute vibrational energy
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
	    
	    cloud_.postReactionVibrationalRedistribution
	    (
	     mode,
	     remainDOF,
	     theta,
	     vibLevel,
	     collisionEnergy
	    );	    
	}        
	
	////////////////-first Trial L-B redistribution (rotation)
	scalar ERotProductMole   = 0.0;
	if( thetaVProductMole.size() > 0 )
	{
	  //const scalar remainDofbyTwo  = (dofP-moleRDof)/2.0+(p.vibLevel().size()-numberVibMode);
	  remainDOF -= moleRDof; 
	  const scalar energyRatioMole = cloud_.postCollisionRotationalEnergy( moleRDof, remainDOF/2.0 );	  
	  ERotProductMole   = energyRatioMole*collisionEnergy;
	  collisionEnergy -= ERotProductMole;//- Relative collisionEnergy energy after rotational energy redistribution	  
	}

	////////////////-second Trial L-B redistribution (rotation)
	scalar ERotProductSecond = 0.0;
	if( thetaVProductSecond.size() > 0 )
	{
	  //const scalar remainDofbyTwo  = (dofP-moleRDof-secondRDof)/2.0+(p.vibLevel().size()-numberVibMode);
	  remainDOF -= secondRDof; 
	  const scalar energyRatioSecond = cloud_.postCollisionRotationalEnergy( secondRDof, remainDOF/2.0 );
	  ERotProductSecond = energyRatioSecond*collisionEnergy;
	  collisionEnergy  -= ERotProductSecond;//- Relative collisionEnergy energy after rotational energy redistribution
	}	

	//////////////redistribute left translational energy///////////////
	//const scalar remainDofbyTwo  = (dofP-moleRDof-secondRDof)/2.0+(p.vibLevel().size()-numberVibMode);
	remainDOF -= 5.0-2.0*reverseOmega; 
	const scalar energyRatioTranslation  = cloud_.postCollisionRotationalEnergy( 5.0-2.0*reverseOmega, remainDOF/2.0 );
	const scalar translationalEnergyLeft = energyRatioTranslation*collisionEnergy;
	collisionEnergy                     -= translationalEnergyLeft;

	/* //not good the vibrational distributtion will not become boltzman
        //- Energy redistribution for particle Q
        cloud_.binaryCollision().redistribute
        (
            q, collisionEnergy, omegaPQ, true
        );       
	*/

	
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
        
	///////////////////////////////end/////////////////////////////////////////
