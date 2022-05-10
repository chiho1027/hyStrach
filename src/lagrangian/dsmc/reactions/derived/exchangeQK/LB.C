	////////////////////////////////////// start ///////////////////////////////////////////////
	const scalar& vibEP =
	  cloud_.constProps(p.typeId()).eVib_tot(p.vibLevel());	
	const scalar& vibEQ =
	  cloud_.constProps(q.typeId()).eVib_tot(q.vibLevel());
	
	scalar collisionEnergy = translationalEnergy + heat
	  + p.ERot() + q.ERot() + vibEP + vibEQ;
	
	scalar remainDof       = moleRDof + secondRDof + 2.0*(2.5 - reverseOmega);
	labelList vibLevelMole(thetaVProductMole.size(), 0);
	labelList vibLevelSecond(thetaVProductSecond.size(), 0);	
	forAll(vibLevelMole, m)
	{	  
	  cloud_.postReactionVibrationalRedistribution
	  (
	   m,
	   remainDof,
	   thetaVProductMole,
	   &vibLevelMole,
	   collisionEnergy
	  );
	}

	forAll(vibLevelSecond, m)
	{	  
	  cloud_.postReactionVibrationalRedistribution
	  (
	   m,
	   remainDof,
	   thetaVProductSecond,
	   &vibLevelSecond,
	   collisionEnergy
	  );
	}       	
	
	scalar ERotProductMole   = 0.0;
	scalar ERotProductSecond = 0.0;
	if(moleRDof>0.0 && secondRDof>0.0)
	{
	  remainDof                   -= moleRDof;
	  scalar energyRatioMole = cloud_.postCollisionRotationalEnergy( moleRDof, remainDof/2.0);
	  ERotProductMole              = energyRatioMole*collisionEnergy;
	  collisionEnergy             -= ERotProductMole;

	  remainDof                   -= secondRDof;
	  scalar energyRatioSecond     = cloud_.postCollisionRotationalEnergy( secondRDof, remainDof/2.0);	  
	  ERotProductSecond            = energyRatioSecond*collisionEnergy;
	  collisionEnergy             -= ERotProductSecond;
	}
	else if(moleRDof>0.0)
	{
	  remainDof                   -= moleRDof;
	  scalar energyRatioMole = cloud_.postCollisionRotationalEnergy( moleRDof, remainDof/2.0);
	  ERotProductMole              = energyRatioMole*collisionEnergy;
	  collisionEnergy             -= ERotProductMole;
	}
	
	//////////////////- asign energy and velocity ///////////////////////////////////////
        const scalar relVelExchMol = sqrt(2.0*collisionEnergy/mRExch);//collisionEnergy LB method  //reactedEnergy present
	
	//////////////////////////////////////  end ///////////////////////////////////////
