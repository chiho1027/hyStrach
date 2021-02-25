/*********************************************************************
1. select react vibrational mode through QK exchange probability
2. sample vibrational level through QK exchange probablility for select mode
3. other mode  vibrational level sample is bolzman distribution
*************************************************************************/
// determin which vibration mode will participate reaction
	label      excite = -1;
	scalar     totalP = 0.0;
	scalarList exciteP(thetaVProductMole.size()+thetaVProductSecond.size(), 0.0);
	forAll(thetaVProductMole, u)
	{
	  exciteP[u] = 1./(1.+ 2.*(1.+1./(2.5-reverseOmega))/thetaVProductMole[u]*TMacro);
	  //exciteP[u] = 1./(1.+ 1000./thetaVProductMole[u]);
	  //exciteP[u] = thetaVProductMole[u];
	  totalP    += exciteP[u];
	}
	forAll(thetaVProductSecond, u)
	{
	  exciteP[thetaVProductMole.size()+u] = 1./(1.+ 2.*(1.+1./(2.5-reverseOmega))/thetaVProductMole[u]*TMacro);
	  //exciteP[u] = 1./(1.+ 1000./thetaVProductMole[u]);
	  //exciteP[thetaVProductMole.size()+u] = thetaVProductSecond[u];
	  totalP    += exciteP[u];
	}

	const scalarList normalisedP = exciteP/totalP;
	const labelList sortedNormalisedPIndices = decreasing_sort_indices(normalisedP);	
	const scalar rndm = cloud_.rndGen().sample01<scalar>();	
	scalar cumulativeP = 0.0;
	forAll(sortedNormalisedPIndices, idx)
        {
	  const label i = sortedNormalisedPIndices[idx];
	  cumulativeP += normalisedP[i];
            
	  if (cumulativeP > rndm)
	  {
	    excite = i;

	    const scalar preCollisionEnergy =
	      translationalEnergy
	      + p.vibLevel()[excite]
	      *thetaVProductMole[excite]
	      *physicoChemical::k.value()
	      + heatOfReactionExchangeJoules_[nExIndex];
	    
	    if(excite < thetaVProductMole.size() )
	    {
	      // vibMole
	      /*
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
	      */
	      if(thetaVProductMole.size() > 0)
		{
		  const scalar kBByThetaVP = physicoChemical::k.value()*thetaVProductMole[excite];
		  scalar func = 0.0;
		  label j = 0;
		  const label iaP = preCollisionEnergy/kBByThetaVP;

		  if(preCollisionEnergy > activationEnergy)
		    {
		      do // acceptance - rejection
			{
			  j = cloud_.randomLabel(0, iaP);
			  
			  scalar summation = 0.0;			  
			  if(activationEnergy < kBByThetaVP)
			    {
			      // this refers to the first sentence in Bird's QK paper after Eq.(12).
			      summation = 1.0; 
			    }
			  else
			    {
			      //const label iaP = preCollisionEnergy/kBByThetaVP;
			      
			      for(label i=0; i<=iaP; i++)
				{
				  summation += 
				    pow
				    (
				     1.0 - i*kBByThetaVP/preCollisionEnergy,
				     1.5 - reverseOmega
				     );
				}
			    }
			  
			  func =
			    pow
			    (
			     1.0 - activationEnergy/preCollisionEnergy,
			     1.5 - reverseOmega
			     )
			    /summation;
			  
			} while(func < cloud_.rndGen().sample01<scalar>());
		    }
		  
		  vibLevelMole[excite] = j;
		  
		  collisionEnergy -= j*kBByThetaVP;
		}
	    }
	    else
	    {
	      // vibSecond
	      /*
	      postReactionVibrationalRedistribution
	      (
	       excite-thetaVProductMole.size(),
	       thetaVProductSecond,
	       TMacro,
	       activationEnergy,
	       reverseOmega,
	       vibLevelSecond,
	       collisionEnergy
	      );
	      */
	    }
	    /*
	    if(collisionEnergy < 0.0)
	    {
	      label newQ =  candidateList[cloud_.randomLabel(0, candidateList.size()-1)];
	      
	      //again choose another q
	      //calculate new collision energy
	      //subtracte i*k*theta
	      nTotExchangeReactions_[nExIndex]--;
	      nExchangeReactionsPerTimeStep_[nExIndex]--;
	      recomputeList

		return;
	    }
	    */

	    break;
	  }
	}

	// sample other vibrational mode by boltzman distribution
	label j = 0;
        scalarList theta = thetaVProductMole;
	labelList* vibLevel =  &vibLevelMole;
	forAll(exciteP, u)
	{
	  label mode = u;
	  if(mode != excite) 
	  {
	    if(mode >= thetaVProductMole.size() )
	    {
	      theta    = thetaVProductSecond;
	      vibLevel = &vibLevelSecond;
	      mode    -= thetaVProductMole.size();
	    }
	    	    
	    const scalar kToTheta = physicoChemical::k.value()*theta[mode];
	    
	    j = -log(cloud_.rndGen().sample01<scalar>())*TMacro/theta[mode];

	    (*vibLevel)[mode]   = j;// vibLevelMole[mode]  = j;
	    
	    collisionEnergy -= j*kToTheta;
	    
	    if(collisionEnergy < 0.0)
	    {
	      collisionEnergy   += j*kToTheta;	      
	      (*vibLevel)[mode]  = int(collisionEnergy/kToTheta);
	      collisionEnergy   -= (*vibLevel)[mode]*kToTheta;
	      
	    }
	  }	  
	  else
	  {
	    continue;
	  }
	}//end for	
