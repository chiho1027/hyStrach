	////////////////////////////////////// start ///////////////////////////////////////////////
	labelList vibLevelMole(thetaVProductMole.size(), 0);
	labelList vibLevelSecond(thetaVProductSecond.size(), 0);

	// determin which vibration mode will participate reaction for product	
	DynamicList<scalar> productExciteP(0);      
	forAll(thetaVProductMole, u)
	{
	  productExciteP.append( 1./(
				     0.5+ reactedEnergy/
				     (
				      (2.5-reverseOmega)*thetaVProductMole[u]*physicoChemical::k.value()
				     )
				    )
			       );				
	}
	forAll(thetaVProductSecond, u)
	{
	  productExciteP.append( 1./(
				     0.5+ reactedEnergy/
				     (
				      (2.5-reverseOmega)*thetaVProductSecond[u]*physicoChemical::k.value()
				     )
				    )
			       );
	}

	//Info << "p = " <<productExciteP << endl;
	
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
	const scalar TMacro = cloud_.fields().overallT(p.cell());
	
	// determine equilibrium energy
	const scalar& pRDof = cloud_.constProps(typeIdP).rotationalDegreesOfFreedom();
	const scalar& qRDof = cloud_.constProps(typeIdQ).rotationalDegreesOfFreedom();
	
	const scalar ERotP = p.ERot();
	const scalar ERotQ = q.ERot();

	scalar equilibriumE = ERotP + ERotQ;
	scalar remainDof    = pRDof + qRDof;

	DynamicList<scalar> equilibriumVE(0);
	DynamicList<scalar> equilibriumVDof(0);

	//add equilibrium vibrational energy
	const label index = productExciteP.size() - 2;
	forAll(p.vibLevel(), m)
	{
	  if(m != selectMode)
	  {
	    const scalar& equiVibEP =
	      cloud_.constProps(typeIdP).eVib_m(m, p.vibLevel()[m]);
	    equilibriumVE.append(equiVibEP);

	    const scalarList thetaVP(1,cloud_.constProps(typeIdP).thetaV()[m]);	    
	    equilibriumVDof.append(2.0*cloud_.temperatureFunction(TMacro, thetaVP));	    
	  }
	}

	forAll(q.vibLevel(), m)
	{
	  const scalar& equiVibEQ =
	    cloud_.constProps(q.typeId()).eVib_m(m, q.vibLevel()[m]);
	  equilibriumVE.append(equiVibEQ);

	  const scalarList thetaVQ(1,cloud_.constProps(typeIdQ).thetaV()[m]);	  
	  equilibriumVDof.append(2.0*cloud_.temperatureFunction(TMacro, thetaVQ));	  
	}	

	label count = 0;
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

	    // the last vibrational energy sample with continue energy other sample by quantum energy
	    if(index == 0)
	    {
	      Info << "a = " << remainDof<<  endl;
	      remainDof -= 2.0;
	      const scalar energyRatio = cloud_.postCollisionRotationalEnergy( 2.0, remainDof/2.0);
	      const scalar kBByThetaV  = physicoChemical::k.value()*theta[mode];
	      (*vibLevel)[mode]        = int(energyRatio*equilibriumE/kBByThetaV);	      
	      equilibriumE            -= (*vibLevel)[mode]*kBByThetaV;
	      //add dof ???

	      Info << "b = " << remainDof<< endl;
	    }
	    else
	    {
	      Info << "c = " << remainDof<< endl;
	      equilibriumE += equilibriumVE[count];
	      remainDof    += equilibriumVDof[count];
	      
	      const scalarList thetaVMode(1, theta[mode]);	      
	      remainDof -= 2.0*cloud_.temperatureFunction(TMacro, thetaVMode);	    	    
	 
	      //method LB quantum level
	      cloud_.postReactionVibrationalRedistribution
	      (
	       mode,
	       remainDof,
	       theta,
	       vibLevel,
	       equilibriumE
	      );

	      count += 1;
	    }	      
	  }	  
	  else
	  {
	    continue;
	  }
	}//end for		

	//add other non participate vibrational energy
	for(label i=index; i<equilibriumVE.size(); i++)
	{
	  equilibriumE += equilibriumVE[i];
	  remainDof    += equilibriumVDof[i];
	}
	Info << "d = " << remainDof<< endl;
	/*
	if(postExciteMode != 0)
	{
	  Info << "14O2re = " <<  vibLevelMole[0] << endl;
	}       	
	*/

	/*
	if(vibLevelMole.size() == 3)
	{
	  Info << "14HO2re = " << vibLevelMole[0] << " " << vibLevelMole[1] << " " << vibLevelMole[2] << endl;
	}
	else
	{
	  Info << "14O2re = " <<vibLevelMole[0] << endl;
	  Info << "14H2re = " <<vibLevelSecond[0] << endl;
	}
	*/
	
	/*
	if(postExciteMode != 1)
	{
	  Info << "13H2Ore1 = " << vibLevelMole[1] << endl;
	}

	if(postExciteMode != 2)
	{
	  Info << "13H2Ore2 = " << vibLevelMole[2] << endl;
	}
	*/

	//relax_ = true;
	//return; 
	
  
	scalar ERotProductMole   = 0.0;
	scalar ERotProductSecond = 0.0;
	if(moleRDof>0.0 && secondRDof>0.0)
	{
	  remainDof -= moleRDof;
	  const scalar energyRatioMole = cloud_.postCollisionRotationalEnergy( moleRDof, remainDof/2.0);
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
