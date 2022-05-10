    ////////////////////////start//////////////////////////////////////    
    //determine which one vibratinal mode will be activate    

    scalar TMacro = cloud_.fields().overallT(p.cell());
    
    label  exciteMode = 0;
    if(thetaVProduct.size() > 1)
    {
      DynamicList<scalar> productExciteP(0,0.0);
      
      forAll(thetaVProduct ,m)
      {	
	const label iMax = (heatOfReactionRecombinationJoules_)/(physicoChemical::k.value()*thetaVProduct[m]);
	
	scalar temp = 0.0;
	
	for(label i=0; i<=iMax; i++)
	{
	  temp += cloud_.normalizedIncomGamma(2.5-reverseOmega, (iMax-i+1)*thetaVProduct[m]/TMacro)
	    *exp(-i*thetaVProduct[m]/TMacro)*(1.0-exp(-thetaVProduct[m]/TMacro));
	}
	
	if(temp == 0.0)
	{
	  temp = VSMALL;
	}
	
	productExciteP.append(temp);	
      }

      exciteMode = selectExciteMode(productExciteP);
    }

    scalar preCollisionEnergy = productThirdBodyTranslationalEnergy + heatOfReactionRecombinationJoules_;
	
    //calculate excite mode vibrational energy using translational energy of product + thirdbody to redistribute 
    const scalar kBByThetaVP = physicoChemical::k.value()*thetaVProduct[exciteMode];
    const label iMax = preCollisionEnergy/kBByThetaVP;
    if(iMax == 0)
    {
      vibLevelProduct[exciteMode] = 0;
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
	    pow(1.0 - j*kBByThetaVP/preCollisionEnergy , 1.5 - reverseOmega);//reverseOmega
	  
      }while(func < cloud_.rndGen().sample01<scalar>() );

      vibLevelProduct[exciteMode] = j;
      preCollisionEnergy         -= j*kBByThetaVP;
    }

    
    //calculate post velocity using remain translational energy of product + thirdbody
    //vector thirdBodyU = thirdBody.U();//temp
    
    scalar relVelProuduct = sqrt(2.0*preCollisionEnergy/mR);
    cloud_.binaryCollision().postCollisionVelocities
    (
     typeIdThirdBody,
     typeIdRecombinedMole,
     thirdBody.U(),
     postCollisionU,
     relVelProuduct
    );      
    
    //redistribution nonExcite  internal collisional energy
    scalar collisionEnergy = translationalEnergy;    
    scalar remainDof       = 2.0;
    
    //first is vibrational energy , second is rotational energy
    //calculate pre-collid particles
    if( cloud_.constProps(typeIdP).type() >= 20)
    {
      const scalar ERotP = p.ERot();
      const scalar EVibP = cloud_.constProps(typeIdP).eVib_tot(p.vibLevel());                
      const scalar& pRDof = cloud_.constProps(typeIdP).rotationalDegreesOfFreedom();
      const scalarList& thetaVP = cloud_.constProps(typeIdP).thetaV();	 
      
      collisionEnergy += EVibP + ERotP;
      remainDof       += pRDof + 2.0*cloud_.temperatureFunction(TMacro, thetaVP);
    }
    
    if( cloud_.constProps(typeIdQ).type() >= 20)
    {	
      const scalar ERotQ = q.ERot();
      const scalar EVibQ = cloud_.constProps(typeIdQ).eVib_tot(q.vibLevel());
      const scalar& qRDof = cloud_.constProps(typeIdQ).rotationalDegreesOfFreedom();
      const scalarList& thetaVQ = cloud_.constProps(typeIdQ).thetaV();	  
      
      collisionEnergy += EVibQ + ERotQ;
      remainDof       += qRDof + 2.0*cloud_.temperatureFunction(TMacro, thetaVQ);
    }        

    labelList* vibLevel  = &vibLevelProduct;    
    forAll(thetaVProduct ,m)
    {
      if(m != exciteMode)
      {
	const scalarList thetaVMode(1, thetaVProduct[m]);
	remainDof -= 2.0*cloud_.temperatureFunction(TMacro, thetaVMode);
	
	//sample with discrete LB // present used
	cloud_.postReactionVibrationalRedistribution
	(
	 m,
	 remainDof,
	 thetaVProduct,
	 vibLevel,
	 collisionEnergy
	);			  
      }
      else
      {
	continue;
      }
    }

    const scalar ERotProduct = collisionEnergy;

    /*
    //temp    
    Info << "12H2Ore = "
	 << vibLevelProduct[0] << " "
	 << vibLevelProduct[1] << " "
	 << vibLevelProduct[2] << endl;
    

    if(exciteMode != 0)
    {
      Info <<"12H2Ore0 = " << vibLevelProduct[0] << endl;;
    }

    if(exciteMode != 1)
    {
      Info <<"12H2Ore1 = " << vibLevelProduct[1] << endl;;
    }
    if(exciteMode != 2)
    {
      Info <<"12H2Ore2 = " << vibLevelProduct[2] << endl;;
    }
    
    relax_ = true;
    return;
    */
      
    //Info <<"2OHre = " << vibLevelProduct[0] << endl;


    /////////////////////end////////////////////////////////
