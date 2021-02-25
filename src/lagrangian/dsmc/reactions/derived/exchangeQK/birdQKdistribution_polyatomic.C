assume all polyatomic poticle each mode of  vibrational temperature is same so select ramdomly one mode which only this mode "react"  and other mode will follow the L-B distribution.


////////////////////////////////////// start ///////////////////////////////////////////////
	
	/*
	if(aCoeff_[nExIndex][0] == 0.053)
	  {
	    pre = collisionEnergy/cloud_.constProps(typeIdMole).nVibrationalModes();
	      //((5.0 - 2.0*reverseOmega+ moleVDof)/( 5.0 - 2.0*reverseOmega + moleRDof + moleVDof + secondRDof + secondVDof) );
	  }
	*/
	//scalar pre = collisionEnergy/cloud_.constProps(typeIdMole).nVibrationalModes();

	//Info << "id = " << typeIdP << endl;
	//Info << "a  = " << ((moleVDof)/( 5.0 - 2.0*reverseOmega + moleRDof + moleVDof + secondRDof + secondVDof) ) << endl;

	
	collisionEnergy = translationalEnergy + heat + q.ERot() + p.ERot();

	//assume radom choose if vibrational temperature at each mode is same
	const scalarList& vibLevelpreP = cloud_.constProps(typeIdP).eVib(p.vibLevel());
	
        label excite = 0;
        scalar EcPlusVib = 0.0;
	const scalar tempEc = collisionEnergy; 
	do
	{
	  EcPlusVib = tempEc;
	  excite = cloud_.randomLabel(0, thetaVProductMole.size()-1);
	  EcPlusVib += vibLevelpreP[excite];
	  
	}while( EcPlusVib  < 0);

	collisionEnergy += vibLevelpreP[excite];
	
	//vibMole	
	labelList vibLevelMole(cloud_.constProps(typeIdMole).nVibrationalModes(), 0);
	postReactionVibrationalRedistribution
	(
	 excite,
	 thetaVProductMole,
	 TMacro,
	 activationEnergy,
	 reverseOmega,
	 vibLevelMole,
	 collisionEnergy
	);

	//temp is (remain dof of vib)/2 initial value is 2.0 //assumpe fully excite
	scalar temp = 2.0;
      
	forAll(vibLevelMole, u)
	{
	  if(u != excite)
	  {
	    collisionEnergy += vibLevelpreP[u];
	    
	    const scalar kBThetaVP = physicoChemical::k.value()*thetaVProductMole[u];	    
	    const scalar iMaxr     = collisionEnergy/kBThetaVP;

	    scalar w    = 1.5-reverseOmega + moleRDof/2.0 + secondRDof/2.0 + temp;
	    scalar func = 0.0;
	    label  j    = 0;
	    do // acceptance - rejection
	    {
	      //iDash = rndGen_.position<label>(0, iMax); OLD
	      j = cloud_.randomLabel(0, iMaxr);
	      
	      // - equation 5.61, Bird
	      func = pow(1.0 - j*kBThetaVP/collisionEnergy, w);
	      
	    }while(func < cloud_.rndGen().sample01<scalar>());
	    
	    vibLevelMole[u]  = j;
	    collisionEnergy -= j*kBThetaVP;

	    temp -= 1.0;
	  }
	  else
	  {
	    continue;
	  }
	}//end for

	
	// vibSecond
	labelList vibLevelSecond(cloud_.constProps(typeIdSecond).nVibrationalModes(), 0);
	postReactionVibrationalRedistribution
	(
	 excite,
	 thetaVProductSecond,
	 TMacro,
	 activationEnergy,
	 reverseOmega,
	 vibLevelSecond,
	 collisionEnergy
	);	
	
	////////////////-first Trial L-B redistribution (rotation)	
	scalar energyRatioMole   = cloud_.postCollisionRotationalEnergy( moleRDof, 2.5 - reverseOmega + secondRDof/2.0 );
        scalar ERotProductMole   = energyRatioMole*collisionEnergy;
	collisionEnergy -= ERotProductMole;//- Relative collisionEnergy energy after rotational energy redistribution

	////////////////-second Trial L-B redistribution (rotation)
	scalar ERotProductSecond = 0.0;
	if( thetaVProductSecond.size() > 0 )
	{
	  const scalar energyRatioSecond = cloud_.postCollisionRotationalEnergy( secondRDof, 2.5 - reverseOmega );
	  ERotProductSecond = energyRatioSecond*collisionEnergy;
	  collisionEnergy -= ERotProductSecond;//- Relative collisionEnergy energy after rotational energy redistribution
	}
	
	/*
	scalar temp = 0; //- record previbrational energy	
        scalar preTranslationalE = collisionEnergy*
	  ((5.0 - 2.0*reverseOmega)/( 5.0 - 2.0*reverseOmega + moleRDof + moleVDof + secondRDof + secondVDof) );

	//vibMole
	labelList vibLevelMole(cloud_.constProps(typeIdMole).nVibrationalModes(), 0);
	//const scalarList& thetaVProductMole = cloud_.constProps(typeIdMole).thetaV();
	collisionEnergy -= preTranslationalE;
	temp = collisionEnergy;
	cloud_.postReactionVibrationalRedistribution
	  (
	    thetaVProductMole,
	    TMacro,
	    moleVDof/2.0,
	    moleRDof/2.0 + secondRDof/2.0 + secondVDof/2.0, //remain dof
	    vibLevelMole,
	    collisionEnergy
	  );
	scalar preVibEMole = temp - collisionEnergy;

	// vibSecond
	labelList vibLevelSecond(cloud_.constProps(typeIdSecond).nVibrationalModes(), 0);
	//const scalarList& thetaVProductSecond = cloud_.constProps(typeIdSecond).thetaV();
	temp = collisionEnergy;
	cloud_.postReactionVibrationalRedistribution
	  (
	    thetaVProductSecond,
	    TMacro,
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
	    TMacro,
	    1.0,//moleVDof/2.0,
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
	    TMacro,
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
	*/
	//////////////////////////////////////  end ///////////////////////////////////////



postReactionVibraiontalRedistribution

void exchangeQK::postReactionVibrationalRedistribution
(
 const label excite,
 const scalarList thetaVProduct,
 const scalar TMacro,
 const scalar activationEnergy,
 const scalar reverseOmega,
 labelList& vibLevel,
 scalar& Ec
)
{
  if( thetaVProduct.size() > 0 )
  {
    //const scalar equiE = Ec/thetaVProduct.size();
    //scalar remain = 0.0;
    
    //or random choose
    //const scalar Et = (2.5-reverseOmega)*physicoChemical::k.value()*TMacro;
    //const scalar Evpar = (Ec-Et)/(thetaVProduct.size()-1.0);

    const label m = excite;    
    
    //forAll(vibLevel, m)
    //{    
      // const 
      const scalar thetaPrim   = thetaVProduct[m]/TMacro; 
      const scalar kBByThetaVP = physicoChemical::k.value()*thetaVProduct[m];
       
      const scalar w           = 1.5 - reverseOmega;
      const scalar EaPrim      = activationEnergy/kBByThetaVP;
      const scalar lnTwo       = log(2.0)/thetaPrim;
      const scalar three       = 3.0/thetaPrim;
      const label  iMax        = Ec/kBByThetaVP;

      //scalar kmax = 0.0;
      //scalar jmax = 0.0;
      //kmaxFunc(lnTwo, three, EaPrim , w, thetaPrim, kmax, jmax);

      //jmax at j = 0
      scalar kmax = (4.0*ki(lnTwo, 0, false, EaPrim, w) + ki(three, 0, false, EaPrim, w) );
 
      scalar func = 0.0;
      label j = 0;

      
      do // acceptance - rejection
      {
	j = cloud_.randomLabel(0, iMax);
	//j = (Ec/kBByThetaVP)*cloud_.rndGen().sample01<scalar>();
	
	if(j > EaPrim )
        {
	  func = (4.0*ki(lnTwo, j, true, EaPrim, w) + ki(three, j, true, EaPrim, w) )//*thetaVProduct.size()
	         *exp(-j*thetaPrim + EaPrim*thetaPrim)
	         /kmax;
	}
	else
	{
	  func = (4.0*ki(lnTwo, j, false, EaPrim, w) + ki(three, j, false, EaPrim, w) )//*thetaVProduct.size()
	    /kmax;
	}	
      } while(func < cloud_.rndGen().sample01<scalar>());

      
      vibLevel[m] = j;
      Ec -= j*kBByThetaVP;

      //}end for

      /*
      scalar temp = 1.0;
      
      forAll(vibLevel, u)
      {
	if(u != excite)
	{	  
	  const scalar kBThetaVP = physicoChemical::k.value()*thetaVProduct[u];
	  const scalar dof       = ;
	  
	  Ec += vibLevelpre[u]*kBThetaVP;
	  
	  const scalar iMaxr     = Ec/kBThetaVP;
	   
	  do // acceptance - rejection
	  {
            //iDash = rndGen_.position<label>(0, iMax); OLD
	    j = cloud_.randomLabel(0, iMaxr);

            // - equation 5.61, Bird
            func = pow(1.0 - j*kBThetaVP/Ec, w);

	  }while(func < cloud_.rndGen().sample01<scalar>());

	  vibLevel[u] = j;
	  Ec -= j*kBThetaVP;

	}
	else
	{
	  continue;
	}
      }//end for
      */
      

      /*
      if(aCoeff_[0] == 0.05)
      {
	j = iMax*0.9; //OH
      }
      else
      {
	j = iMax*0.9;
      }
      */

      /*
      scalar temp = 0;
      for(label i=0; i<=iMax; i++)
      {	
	if(i > EaPrim )
        {
	  func = (4.0*ki(lnTwo, i, true, EaPrim, w) + ki(three, i, true, EaPrim, w) )
	         *exp(EaPrim*thetaPrim)*(w+1.0)*pow(thetaPrim, w)
	         /(6.0*exp(lgamma(w+1.0))*(1.0-exp(-thetaPrim)))
		 *errorFactor;
	}
	else
	{
	  func = (4.0*ki(lnTwo, i, false, EaPrim, w) + ki(three, i, false, EaPrim, w) )
	         *exp(i*thetaPrim)*(w+1.0)*pow(thetaPrim, w)
	         /(6.0*exp(lgamma(w+1.0))*(1.0-exp(-thetaPrim)))
		 *errorFactor;
	}
	
	//Info <<"activation = " << activationEnergy/physicoChemical::k.value() << endl;
	//Info <<"ki/kf every = "<< func << endl;
	
	if( i == 0 || abs(func - 1.0) < temp )
	{
	  j = i;
	  temp = abs(func - 1.0);
	  //Info <<"ki/kf = "<< func << endl;
	}	  
      }
      */
      
      //Info <<"j real  == "<< j << endl;
      //Info <<"func    == "<< func << endl<< endl;
      /*
      if(aCoeff_[0] == 0.05)
      {
	j= iMax;
      }
      */
      
      //vibLevel[m] = j;
      //Ec -= j*kBByThetaVP;  

      // }//end for
    

    //Ec -=  vibLevel[i]*physicoChemical::k.value()*thetaVProduct[i];
  }//end if   
}
