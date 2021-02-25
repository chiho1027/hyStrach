	////////////////////////////////////// start ///////////////////////////////////////////////
	////////////////// pre-redistribution particle energy////////////////////////////////////
	
	scalar temp = 0; //- record previbrational energy	
        scalar preTranslationalE = collisionEnergy*
	  ((5.0 - 2.0*reverseOmega)/( 5.0 - 2.0*reverseOmega + moleRDof + moleVDof + secondRDof + secondVDof) );

	// rotMole
	scalar energyRatio = 0.0;
	scalar ERotProductMole = 0.0;
	collisionEnergy -= preTranslationalE; //remain energy
	energyRatio      = cloud_.postCollisionRotationalEnergy( moleRDof, moleVDof/2.0 + secondRDof/2.0 + secondVDof/2.0 );
	ERotProductMole  = energyRatio*collisionEnergy;

	// rotSecond
	scalar ERotProductSecond = 0.0;
	collisionEnergy  -= ERotProductMole;
	if(secondRDof > 0.0)
	{
	  energyRatio       = cloud_.postCollisionRotationalEnergy( secondRDof, moleVDof/2.0 + secondVDof/2.0 );
	  ERotProductSecond = energyRatio*collisionEnergy;
	}
	
	//vibMole and vibSecond
	labelList vibLevelMole(cloud_.constProps(typeIdMole).nVibrationalModes(), 0);
	labelList vibLevelSecond(cloud_.constProps(typeIdSecond).nVibrationalModes(), 0);
	scalar preVibEMole   = 0.0;
	scalar preVibESecond = 0.0;
	collisionEnergy  -= ERotProductSecond;
	
	if( thetaVProductSecond.size() > 0 )
	{
	  temp = collisionEnergy;
	  cloud_.postReactionVibrationalRedistribution
	  (
	    thetaVProductMole,
	    moleVDof/2.0,
	    moleRDof/2.0 + secondVDof/2.0, //remain dof
	    vibLevelMole,
	    collisionEnergy
	  );
	  preVibEMole   = temp - collisionEnergy;
	  preVibESecond = collisionEnergy; 
	}
	else
	{
	  preVibEMole = collisionEnergy;
	}
	
	///////////////////////////////////////////////////////////////////////////////////
	////////////////-first Trial L-B redistribution (rotation)
	collisionEnergy   = preTranslationalE + ERotProductMole;	
	energyRatio       = cloud_.postCollisionRotationalEnergy( moleRDof, 2.5 - reverseOmega );	
	ERotProductMole   = energyRatio*collisionEnergy;
	preTranslationalE = collisionEnergy - ERotProductMole;//- Relative collisionEnergy energy after rotational energy redistribution

	////////////////-second Trial L-B redistribution (rotation)
	if( secondRDof > 0.0 )
	{
	  collisionEnergy   = preTranslationalE + ERotProductSecond;	  
	  energyRatio       = cloud_.postCollisionRotationalEnergy( secondRDof, 2.5 - reverseOmega );
	  ERotProductSecond = energyRatio*collisionEnergy;
	  preTranslationalE = collisionEnergy - ERotProductSecond;//- Relative collisionEnergy energy after rotational energy redistribution
	}
	
	/////////////////-first Trial L-B redistribution (vibration)
	collisionEnergy   = preTranslationalE + preVibEMole;
	cloud_.postReactionVibrationalRedistribution
	  (
	   thetaVProductMole,
	   preVibEMole/2.0,
	   2.5 - reverseOmega,
	   vibLevelMole,
	   collisionEnergy
	  );
	preTranslationalE = collisionEnergy; 

	/////////////////-second Trial L-B redistribution (vibration)
	collisionEnergy   = preTranslationalE + preVibESecond;
	cloud_.postReactionVibrationalRedistribution
	  (
	   thetaVProductSecond,
	   preVibESecond/2.0,
	   2.5 - reverseOmega,
	   vibLevelSecond,
	   collisionEnergy
	  );
	preTranslationalE = collisionEnergy; 
	
	//////////////////////////////////////  end ///////////////////////////////////////
