    ////////////////////////////////////////////////////////////////////////////////////////////////
    // initialized post recombinationproduct viblevel
    labelList vibLevelMole(cloud_.constProps(typeIdRecombinedMol).nVibrationalModes(), 0);
    
    //- Trial L-B redistribution (vibration)  //bu i ding
    forAll(thetaVProduct, m)
    {
      label iMaxProduct = collisionEnergy/(physicoChemical::k.value()*thetaVProduct[m]);
      
      label vibLevelProduct =
	cloud_.postCollisionVibrationalEnergyLevel
	(
	 true,
	 0,// vibrationlevel
	 iMaxProduct,
	 cloud_.constProps(typeIdRecombinedMol).thetaV()[m],
	 cloud_.constProps(typeIdRecombinedMol).thetaD()[m],
	 cloud_.constProps(typeIdRecombinedMol).TrefZv()[m],
	 omegaProduct,
	 cloud_.constProps(typeIdRecombinedMol).Zref()[m],
	 collisionEnergy
	 );
      
      vibLevelMole[m] =  vibLevelProduct;
      //- Relative translational energy after vibrational energy redistribution
      collisionEnergy -= (vibLevelProduct*cloud_.constProps(typeIdRecombinedMol).thetaV()[m]*physicoChemical::k.value()); 
    }
    
    //- Trial L-B redistribution (rotation)
    const scalar energyRatio =
      cloud_.postCollisionRotationalEnergy
      (
       cloud_.constProps(typeIdRecombinedMol).rotationalDegreesOfFreedom(),
       2.5 - omegaProduct
       );    
    const scalar ERotProduct = energyRatio*collisionEnergy;
    //- Relative translational energy after rotational energy redistribution
    collisionEnergy -= ERotProduct;

    ////////////////////////////////////////////////
