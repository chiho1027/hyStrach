	////////////////////////////////////// start ///////////////////////////////////////////////
	////////////////// pre-redistribution particle energy////////////////////////////////////		

	//vibMole
	labelList vibLevelMole(cloud_.constProps(typeIdMole).nVibrationalModes(), 0);
	cloud_.postReactionVibrationalRedistribution
	  (
	    thetaVProductMole,
	    moleVDof/2.0,
	    moleRDof/2.0 + secondRDof/2.0 + secondVDof/2.0 + 2.5 - reverseOmega, //remain dof
	    vibLevelMole,
	    collisionEnergy
	  );

	// vibSecond
	labelList vibLevelSecond(cloud_.constProps(typeIdSecond).nVibrationalModes(), 0);
	cloud_.postReactionVibrationalRedistribution
	  (
	    thetaVProductSecond,
	    secondVDof/2.0,
	    moleRDof/2.0 + secondRDof/2.0 +  2.5 - reverseOmega, //remain dof
	    vibLevelSecond,
	    collisionEnergy
	  );
	
	// rotMole
	scalar energyRatio = 0.0;
	scalar ERotProductMole = 0.0;
	scalar ERotProductSecond = 0.0;
	if( thetaVProductSecond.size() > 0 )
	{
	  energyRatio       = cloud_.postCollisionRotationalEnergy( moleRDof, secondRDof/2.0 + 2.5 - reverseOmega );
	  ERotProductMole   = energyRatio*collisionEnergy;
	  collisionEnergy  -= ERotProductMole;

	  energyRatio       = cloud_.postCollisionRotationalEnergy( secondRDof, 2.5 - reverseOmega );
	  ERotProductSecond = energyRatio*collisionEnergy;
	  collisionEnergy  -= ERotProductSecond;
	  
	}
	else
	{
	  energyRatio       = cloud_.postCollisionRotationalEnergy( moleRDof, 2.5 - reverseOmega );
	  ERotProductMole = energyRatio*collisionEnergy;
	  collisionEnergy  -= ERotProductMole;
	}

	// using collisionEnergy to translationalEnergy 
	//////////////////////////////////////  end ///////////////////////////////////////
