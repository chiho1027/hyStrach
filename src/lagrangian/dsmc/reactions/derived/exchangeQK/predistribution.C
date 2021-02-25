	////////////////// pre-redistribution particle energy////////////////////////////////////
        const scalar reverseOmega =  0.5*(
				           cloud_.constProps(typeIdMole).omega()
					   + cloud_.constProps(typeIdSecond).omega()
					 );

	scalar temp = 0; //- record previbrational energy	
        scalar preTranslationalE = collisionEnergy*
	  ((5.0 - 2.0*reverseOmega)/( 5.0 - 2.0*reverseOmega + moleRDof + moleVDof + secondRDof + secondVDof) );

	//vibMole
	labelList vibLevelMole(cloud_.constProps(typeIdMole).nVibrationalModes(), 0);
	const scalarList& thetaVProductMole = cloud_.constProps(typeIdMole).thetaV();
	collisionEnergy -= preTranslationalE;
	temp = collisionEnergy;
	cloud_.postReactionVibrationalRedistribution
	  (
	    thetaVProductMole,
	    moleVDof/2.0,
	    moleRDof/2.0 + secondRDof/2.0 + secondVDof/2.0, //remain dof
	    vibLevelMole,
	    collisionEnergy
	  );
	scalar preVibEMole = temp - collisionEnergy;

	// vibSecond
	labelList vibLevelSecond(cloud_.constProps(typeIdSecond).nVibrationalModes(), 0);
	const scalarList& thetaVProductSecond = cloud_.constProps(typeIdSecond).thetaV();
	temp = collisionEnergy;
	cloud_.postReactionVibrationalRedistribution
	  (
	    thetaVProductSecond,
	    secondVDof/2.0,
	    moleRDof/2.0 + secondRDof/2.0, //remain dof
	    vibLevelSecond,
	    collisionEnergy
	  );
	scalar preVibESecond = temp - collisionEnergy;
	
	// rotMole
	scalar energyRatioMole = 0.0;
	scalar ERotProductMole = 0.0;
	scalar ERotProductSecond = 0.0;
	if( thetaVProductSecond.size() > 0 )
	{
	  energyRatioMole   = cloud_.postCollisionRotationalEnergy( moleRDof, secondRDof/2.0 );
	  ERotProductMole   = energyRatioMole*collisionEnergy;
	  ERotProductSecond = collisionEnergy - ERotProductMole ;
	}
	else
	{
	  ERotProductMole = collisionEnergy;
	  //ERotProductSecond = 0.0;
	}
	
	///////////////////////////////////////////////////////////////////////////////////
	/////////////////-first Trial L-B redistribution (vibration)
	collisionEnergy = preTranslationalE + preVibEMole;
	cloud_.postReactionVibrationalRedistribution
	  (
	    thetaVProductMole,
	    moleVDof/2.0,
	    2.5 - reverseOmega, //QK
	    vibLevelMole,
	    collisionEnergy
	  );
	preTranslationalE = collisionEnergy;

	/////////////////-second Trial L-B redistribution (vibration)
	collisionEnergy = preTranslationalE + preVibESecond;
	cloud_.postReactionVibrationalRedistribution
	  (
	    thetaVProductSecond,
	    secondVDof/2.0,
	    2.5 - reverseOmega, //QK
	    vibLevelSecond,
	    collisionEnergy
	  );
	preTranslationalE = collisionEnergy;
	
	////////////////-first Trial L-B redistribution (rotation)
	collisionEnergy = preTranslationalE + ERotProductMole;	
	energyRatioMole = cloud_.postCollisionRotationalEnergy( moleRDof, 2.5 - reverseOmega );	
	ERotProductMole = energyRatioMole*collisionEnergy;
	preTranslationalE = collisionEnergy - ERotProductMole;//- Relative collisionEnergy energy after rotational energy redistribution
	
	////////////////-second Trial L-B redistribution (rotation)
	if( thetaVProductSecond.size() > 0 )
	{
	  collisionEnergy = preTranslationalE + ERotProductSecond;	  
	  scalar energyRatioSecond = cloud_.postCollisionRotationalEnergy( secondRDof, 2.5 - reverseOmega );
	  ERotProductSecond = energyRatioSecond*collisionEnergy;
	  preTranslationalE = collisionEnergy - ERotProductSecond;//- Relative collisionEnergy energy after rotational energy redistribution
	}
