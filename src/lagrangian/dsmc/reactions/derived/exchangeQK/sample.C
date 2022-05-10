	//test 8O2_2400K equilibrium total collisionenergy distribution
	/*
	if(heatOfReactionExchangeJoules_[0] < 0)
	{
	  scalar noOfKT = (
			   translationalEnergy +
			   p.vibLevel()[0]*physicoChemical::k.value()*cloud_.constProps(typeIdP).thetaV()[0]
			  )
	    /(0.05*physicoChemical::k.value()*2400.0);
			
	  scalar colliDifference = noOfKT-int(noOfKT);	        
	      
	  if(colliDifference < 0.1 )
	  {				    
	    Info << "eqTotalDis = "
		 << int(noOfKT)
		 << endl;
	  }
	  else if( colliDifference >  0.9 )
	  {	    
	    Info << "eqTotalDis = "
		 << int(noOfKT)+1
		 << endl;
	  }
	}
	*/


/*
	    if(heatOfReactionExchangeJoules_[0] < 0)
	    {
	      Info << "13H2Ofor = " << p.vibLevel()[0]<< " " << p.vibLevel()[1]<< " "  << p.vibLevel()[2] << endl;
	    }
	    else
	    {
	      Info << "13HO2for = " << p.vibLevel()[0]<< " " << p.vibLevel()[1]<< " "  << p.vibLevel()[2] << endl;
	    }
	    */
	    
	    //if(cloud_.constProps(typeIdP).thetaV_m(0) == 4950.0)
	    //{
	    // Info << "15HO2for = " << p.vibLevel()[0]<< " " << p.vibLevel()[1]<< " "  << p.vibLevel()[2] << endl;
	    //}

	    /*
	    if( heatOfReactionExchangeJoules_[0] < 0)
	    {
	      Info << "8O2for = " << p.vibLevel()[0] << endl;
	    }
	    else
	    {
	      Info << "8OHfor = " << p.vibLevel()[0] << endl;
	    }
	    */

	    //test 8O2_2400K forward total collisionenergy distribution
	    //if(heatOfReactionExchangeJoules_[0] < 0)
	    /*
	    if(cloud_.constProps(typeIdP).thetaV_m(0) == 2256.0)
	    {		     	      
	      scalar noOfKT = (
			       translationalEnergy +
			       p.vibLevel()[0]*physicoChemical::k.value()*cloud_.constProps(typeIdP).thetaV()[0]
			      )
		/(0.05*physicoChemical::k.value()*2400);
			
	      scalar colliDifference = noOfKT-int(noOfKT);	        
	      
	      if(colliDifference < 0.1 )
	      {				    
		Info << "forTotalDis = "
		     << int(noOfKT)
		     << endl;
	      }
	      else if( colliDifference >  0.9 )
	      {	    
		Info << "forTotalDis = "
		     << int(noOfKT)+1
		     << endl;
	      }
	    }
	    */




/*
	    if(heatOfReactionExchangeJoules_[0] < 0)
	    {
	      Info << "13H2Ofor = "
		   << q.vibLevel()[0]<< " "
		   << q.vibLevel()[1]<< " "
		   << q.vibLevel()[2]<< endl;
	    }
	    else
	    {
	      Info << "13HO2for = "
		   << q.vibLevel()[0]<< " "
		   << q.vibLevel()[1]<< " "
		   << q.vibLevel()[2]<< endl;
	    }
	    */
	    //if(cloud_.constProps(typeIdQ).thetaV_m(0) == 4950.0)
	    //{
	    // Info << "15HO2for = " << q.vibLevel()[0]<< " " << q.vibLevel()[1]<< " "  << q.vibLevel()[2] << endl;
	    //}

	    //if( heatOfReactionExchangeJoules_[0] < 0)
	    //{
	    //  Info << "8O2for = " << q.vibLevel()[0] << endl;
	    //}
	    //else
	    //{
	    //  Info << "8OHfor = " << q.vibLevel()[0] << endl;
	    //}	   
