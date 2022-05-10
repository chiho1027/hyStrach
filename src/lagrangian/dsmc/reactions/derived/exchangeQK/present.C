	////////////////////////////////////// start ///////////////////////////////////////////////
	const scalar TMacro = cloud_.fields().overallT(p.cell());
	
	labelList vibLevelMole(thetaVProductMole.size(), 0);
	labelList vibLevelSecond(thetaVProductSecond.size(), 0);

	const labelList vibBolzMole =
	  cloud_.equipartitionVibrationalEnergyLevel(TMacro, vibLevelMole.size(), typeIdMole);
	const labelList vibBolzSecond =
	   cloud_.equipartitionVibrationalEnergyLevel(TMacro, vibLevelSecond.size(), typeIdSecond);
	
	// determin which vibration mode will participate reaction for product	
	DynamicList<scalar> productExciteP(0);      
	forAll(thetaVProductMole, u)
	{	  
	  scalar sum = 0.0;
	  const scalar equiTranE = cloud_.equipartitionRotationalEnergy(TMacro, 5.0-2.0*reverseOmega);
	  const scalar kThetaV = thetaVProductMole[u]*physicoChemical::k.value();
	  const scalar Ec      = equiTranE + vibBolzMole[u]*kThetaV;
	  
	  const label iMax     = Ec/kThetaV;
	
	  for(label i=0; i<=iMax; i++)
	  {
	    sum += 
	      pow
	      (
	       1.0 - (i*kThetaV)/Ec,
	       1.5 - reverseOmega
	      );
	  }

	  productExciteP.append(1.0/sum);
	  
	  /*
	  productExciteP.append( 1./(
				     0.5+ Ec/
				     (
				      (2.5-reverseOmega)*thetaVProductMole[u]*physicoChemical::k.value()
				     )
				    )
			       );
	  */	
	}
	
	forAll(thetaVProductSecond, u)
	{	  
	  scalar sum = 0.0;
	  const scalar equiTranE = cloud_.equipartitionRotationalEnergy(TMacro, 5.0-2.0*reverseOmega);
	  const scalar kThetaV = thetaVProductSecond[u]*physicoChemical::k.value();
	  const scalar Ec      = equiTranE + vibBolzSecond[u]*kThetaV;
	  
	  const label iMax     = Ec/kThetaV;
	  
	  for(label i=0; i<=iMax; i++)
	  {
	    sum += 
	      pow
	      (
	       1.0 - (i*kThetaV)/Ec,
	       1.5 - reverseOmega
	      );
	  }

	  productExciteP.append(1.0/sum);
	  
	  /*
	  productExciteP.append( 1./(
				     0.5+ Ec/
				     (
				      (2.5-reverseOmega)*thetaVProductSecond[u]*physicoChemical::k.value()
				     )
				    )
			       );
	  */   
	}

	const label postExciteMode = selectExciteMode(productExciteP);

	label      excite   =  0;
	scalarList theta    =  thetaVProductMole;//assume
	labelList* vibLevel = &vibLevelMole;//assume       
	
	if(postExciteMode < thetaVProductMole.size() )
	{
	  excite   =  postExciteMode;
	}
	else
	{
	  excite   =  postExciteMode-thetaVProductMole.size();
	  theta    =  thetaVProductSecond;	  
	  vibLevel =  &vibLevelSecond;
	}
	
	const scalar kBByThetaVP = physicoChemical::k.value()*theta[excite];
	const label iMax = reactedEnergy/kBByThetaVP;
	
	if(iMax == 0)
	{
	  (*vibLevel)[excite] = 0;
	}
	else
	{
	  scalar func  = 0.0;
	  label  j     =   0;
    
	  //select postProuct pre-reaction vibrational level
	  do // acceptance - rejection
	  {
	    j = cloud_.randomLabel(0, iMax);
	    
	    func = 
	      pow(1.0 - j*kBByThetaVP/reactedEnergy ,1.5 - reverseOmega);
            
	  }while(func < cloud_.rndGen().sample01<scalar>() );

	  (*vibLevel)[excite] = j;	  
	  reactedEnergy    -= j*kBByThetaVP;	  	    
	}

	//claculate post equilibrium energy
	
	// determine equilibrium energy
	const scalar& pRDof = cloud_.constProps(typeIdP).rotationalDegreesOfFreedom();
	const scalar& qRDof = cloud_.constProps(typeIdQ).rotationalDegreesOfFreedom();
	
	const scalar ERotP = p.ERot();
	const scalar ERotQ = q.ERot();

	scalar equilibriumE    = ERotP + ERotQ;
	scalar equilibriumRDof = pRDof + qRDof;

	
	const label index = productExciteP.size() - 2;
	
	scalar equilibriumVE   = 0.0;
	scalar equilibriumVDof = 0.0;
	
	//add equilibrium vibrational energy
	forAll(p.vibLevel(), m)
	{
	  if(m != selectMode)
	  {
	    const scalar& equiVibEP =
	      cloud_.constProps(typeIdP).eVib_m(m, p.vibLevel()[m]);	    
	    const scalarList equiVibTheta(1,cloud_.constProps(typeIdP).thetaV()[m]);
	    
	    equilibriumVDof += 2.0*cloud_.temperatureFunction(TMacro, equiVibTheta);

	    if(index == 0)
	    {
	      equilibriumVE += equiVibEP;
	    }
	    else
	    {
	      equilibriumE  += equiVibEP;
	    }
	  }
	}

       	const scalar& equiVibEQ =
	  cloud_.constProps(q.typeId()).eVib_tot(q.vibLevel());
	const scalarList& thetaVQ = cloud_.constProps(typeIdQ).thetaV();
	  
	equilibriumVDof  += 2.0*cloud_.temperatureFunction(TMacro, thetaVQ);

	if(index == 0)
	{
	  equilibriumVE += equiVibEQ;
	}
	else
	{
	  equilibriumE  += equiVibEQ;
	}        	
	
	forAll(productExciteP, u)
	{
	  label mode = u;
	  
	  if(mode != postExciteMode) 
	  {	    
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

	    if(index == 0)
	    {
	      equilibriumRDof -= 2.0;
	      const scalar energyRatio = cloud_.postCollisionRotationalEnergy( 2.0, equilibriumRDof/2.0);
	      const scalar kBByThetaV  = physicoChemical::k.value()*theta[mode];
	      (*vibLevel)[mode]        = int(energyRatio*equilibriumE/kBByThetaV);	      
	      equilibriumE            -= (*vibLevel)[mode]*kBByThetaV;
	      
	      const scalarList thetaVMode(1, theta[mode]);
	      equilibriumRDof += 2.0 - 2.0*cloud_.temperatureFunction(TMacro, thetaVMode);
	    }
	    else
	    {	      
	      const scalarList thetaVMode(1, theta[mode]);
	      equilibriumVDof  -= 2.0*cloud_.temperatureFunction(TMacro, thetaVMode);	     
	      const scalar remainDof = equilibriumVDof + equilibriumRDof;	      
	      
	      //method LB quantum level	      
	      cloud_.postReactionVibrationalRedistribution
	      (
	       mode,
	       remainDof,
	       theta,
	       vibLevel,
	       equilibriumE
	      );
	    }
	  }
	  else
	  {
	    continue;
	  }
	}//end for
	
	//add other non participate vibrational energy
	equilibriumE    += equilibriumVE;
	equilibriumRDof += equilibriumVDof;
	
	
	scalar ERotProductMole   = 0.0;
	scalar ERotProductSecond = 0.0;
	if(moleRDof>0.0 && secondRDof>0.0)
	{
	  equilibriumRDof -= moleRDof;
	  const scalar energyRatioMole = cloud_.postCollisionRotationalEnergy( moleRDof, equilibriumRDof/2.0);
	  ERotProductMole              = energyRatioMole*equilibriumE;
	  equilibriumE                -= ERotProductMole;
	  ERotProductSecond            = equilibriumE;
	}
	else if(moleRDof>0.0)
	{
	  ERotProductMole              = equilibriumE;
	}
	else
	{
	  ERotProductSecond            = equilibriumE;
	}


	//////////////////- asign energy and velocity ///////////////////////////////////////
        const scalar relVelExchMol = sqrt(2.0*reactedEnergy/mRExch);//collisionEnergy LB method  //reactedEnergy present
	
	//////////////////////////////////////  end ///////////////////////////////////////
