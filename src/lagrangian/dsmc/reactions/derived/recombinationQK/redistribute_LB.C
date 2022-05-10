    ////////////////////////strat//////////////////////////////////////    
    //redistribution nonExcite  internal collisional energy
    scalar collisionEnergy = productThirdBodyTranslationalEnergy + translationalEnergy + heatOfReactionRecombinationJoules_;
    scalar remainDof       = 2.0*(2.5-reverseOmega) + rotDofProduct;

    //first is vibrational energy , second is rotational energy
    //calculate pre-collid particles
    if( cloud_.constProps(typeIdP).type() >= 20)
    {
      const scalar ERotP = p.ERot();
      const scalar EVibP = cloud_.constProps(typeIdP).eVib_tot(p.vibLevel());          

      collisionEnergy += (EVibP + ERotP);      
    }
    else if( cloud_.constProps(typeIdQ).type() >= 20)
    {	
      const scalar ERotQ = q.ERot();
      const scalar EVibQ = cloud_.constProps(typeIdQ).eVib_tot(q.vibLevel());

      collisionEnergy += (EVibQ + ERotQ);    
    }        

    labelList* vibLevel       = &vibLevelProduct;    
    forAll(thetaVProduct ,m)
    {      
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

    //redistribute pruduct of rotational energy
    remainDof                      -= rotDofProduct;
    const scalar energyRatioProduct = cloud_.postCollisionRotationalEnergy( rotDofProduct, remainDof/2.0 );
    const scalar ERotProduct        = energyRatioProduct*collisionEnergy;
    collisionEnergy                -= ERotProduct;
        
    //calculate post velocity using remain translational energy of product + thirdbody
    //vector thirdBodyU = thirdBody.U();//temp   
    scalar relVelProuduct = sqrt(2.0*collisionEnergy/mR);
    cloud_.binaryCollision().postCollisionVelocities
    (
     typeIdThirdBody,
     typeIdRecombinedMole,
     thirdBody.U(),
     postCollisionU,
     relVelProuduct
    );      
    /////////////////////end////////////////////////////////
