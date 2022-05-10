/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "LarsenBorgnakkeVariableHardSphere.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"
#include "exchangeQK.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LarsenBorgnakkeVariableHardSphere, 0);
    addToRunTimeSelectionTable
    (
        BinaryCollisionModel,
        LarsenBorgnakkeVariableHardSphere,
        dictionary
    );
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LarsenBorgnakkeVariableHardSphere::LarsenBorgnakkeVariableHardSphere
(
    const dictionary& dict,
    dsmcCloud& cloud
)
:
    VariableHardSphere(dict, cloud),
    coeffDictLB_
    (
        dict.isDict(typeName + "Coeffs")
        ? dict.subDict(typeName + "Coeffs")
        : dictionary()
    ),
    rotationalRelaxationCollisionNumber_
    (
        coeffDictLB_.lookupOrDefault<scalar>
        (
            "rotationalRelaxationCollisionNumber",
            5.0
        )
    ),
    vibrationalRelaxationCollisionNumber_
    (
        coeffDictLB_.lookupOrDefault<scalar>
        (
            "vibrationalRelaxationCollisionNumber", 
            0.0
        )
    ),
    invZvFormulation_(2),
    electronicRelaxationCollisionNumber_
    (
        coeffDictLB_.lookupOrDefault<scalar>
        (
            "electronicRelaxationCollisionNumber",
            0.0
        )
    )
{
    const word inverseZvFormulationVersion =
        coeffDictLB_.lookupOrDefault<word>
        (
            "inverseZvFormulation", 
            word::null
        );
        
    if (inverseZvFormulationVersion == "pre-2008")
    {
        invZvFormulation_ = 0;
    }
    else if (inverseZvFormulationVersion == "2008")
    {
        invZvFormulation_ = 1;
    }
    else if (inverseZvFormulationVersion == "1/5Z")
    {
       invZvFormulation_ = 2;
    }
    else if (inverseZvFormulationVersion == "correctFactorZc")
    {
      invZvFormulation_ = 3;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LarsenBorgnakkeVariableHardSphere::~LarsenBorgnakkeVariableHardSphere()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LarsenBorgnakkeVariableHardSphere::collide
(
    dsmcParcel& pP,
    dsmcParcel& pQ,
    const label cellI,
    scalar cR
)
{   
    const label typeIdP = pP.typeId();
    const label typeIdQ = pQ.typeId();
    
    vector& UP = pP.U();
    vector& UQ = pQ.U();
    
    const scalar mP = cloud_.constProps(typeIdP).mass();
    const scalar mQ = cloud_.constProps(typeIdQ).mass();
    const scalar mR = mP*mQ/(mP + mQ);
    
    const scalar cRsqr = magSqr(UP - UQ);

    //- Pre-collision relative translational energy
    scalar translationalEnergy = 0.5*mR*cRsqr;
    
    const scalar omegaPQ =
        0.5
       *(
            cloud_.constProps(typeIdP).omega()
          + cloud_.constProps(typeIdQ).omega()
        );

    //if Zr smaller than 5(assme), the overall probability could greater than 1
    if(rotationalRelaxationCollisionNumber_ < 5)
    {
      redistribute(pP, translationalEnergy, omegaPQ, typeIdQ);
      redistribute(pQ, translationalEnergy, omegaPQ, typeIdP);
    }
    else
    {
      redistributeOnlyOneMode(pP, pQ, translationalEnergy ,omegaPQ);
    }

    //- Rescale the translational energy
    cR = sqrt(2.0*translationalEnergy/mR);
    
    VariableHardSphere::scatter(pP, pQ, cellI, cR);
}


void Foam::LarsenBorgnakkeVariableHardSphere::redistribute
(
    dsmcParcel& p,
    scalar& translationalEnergy,
    const scalar omegaPQ,
    const label QId,
    const bool postReaction
)
{
    const label typeIdP = p.typeId();
    const dsmcParcel::constantProperties& cP = cloud_.constProps(typeIdP);
    
    if (cP.type() == 0)
    {
        //- The particle is an electron, no energy to redistribute
        return void();
    }
    
    if(electronicRelaxationCollisionNumber_ != 0)
    {
    
      const scalar inverseElectronicCollisionNumber =
	1.0/electronicRelaxationCollisionNumber_;   
      label& ELevelP = p.ELevel();
    
      //- Electronic energy mode for P
      if (inverseElectronicCollisionNumber > cloud_.rndGen().sample01<scalar>())
      { 
        const label jMaxP = cP.nElectronicLevels();    
        const scalarList& EElistP = cP.electronicEnergyList();    
        const labelList& gListP = cP.electronicDegeneracyList(); 
        const scalar preCollisionEEleP = EElistP[ELevelP];
	
        //- Collision energy of particle P: relative translational energy 
        //   + pre-collision electronic energy
        const scalar EcP = translationalEnergy + preCollisionEEleP;
        
        ELevelP = 
            cloud_.postCollisionElectronicEnergyLevel
            (
                EcP,
                jMaxP,
                omegaPQ,
                EElistP,
                gListP
            );
                        
        //- Relative translational energy after electronic energy exchange
        translationalEnergy = EcP - EElistP[ELevelP];
      }
    }

    //- Vibrational energy mode for P
    if (cP.nVibrationalModes() > 0)
    {
        const scalarList& thetaVP = cP.thetaV();  
        const scalarList& thetaDP = cP.thetaD();
        const List<scalarList>& ZrefP = cP.Zref();
        const List<scalarList>& refTempZvP = cP.TrefZv();
        const scalarList& preCollisionEVibP = cP.eVib(p.vibLevel());
	
	label partnerIndex = findIndex(cP.ZrefPartnerId(), QId);
	if(partnerIndex == -1)
	{
	  partnerIndex = cP.ZrefPartnerMIndex();
	}
	
	//calculate correction Z fector
	
	
	//if(1.0/vibrationalRelaxationCollisionNumber_ > cloud_.rndGen().sample01<scalar>())
	//{
	forAll(thetaVP, i)
	{
	  //label i = cloud_.randomLabel(0, thetaVP.size()-1);
	
            //- Collision energy of particle P: relative translational energy 
            //    + pre-collision vibrational energy
            const scalar EcP = translationalEnergy + preCollisionEVibP[i]; 

            //- Maximum possible quantum level (equation 3, Bird 2010)
            const label iMaxP = EcP/(physicoChemical::k.value()*thetaVP[i]);

	    const scalar correctZFactorP = cloud_.fields().correctZFactor(p.cell(), typeIdP, i);
	    
            if (iMaxP > 0)
            {	      	      
                p.vibLevel()[i] = 
                    cloud_.postCollisionVibrationalEnergyLevel
                    (
                        postReaction,
                        p.vibLevel()[i],
                        iMaxP,
                        thetaVP[i],
                        thetaDP[i],
                        refTempZvP[partnerIndex][i],
                        omegaPQ,
                        ZrefP[partnerIndex][i],
                        EcP,
			correctZFactorP,
                        vibrationalRelaxationCollisionNumber_,
                        invZvFormulation_,
                        p.cell()
                    );		
		
                translationalEnergy = EcP - cP.eVib_m(i, p.vibLevel()[i]);
            }	    
	}
	//}
    }    

    //- Rotational energy mode for P
    const scalar rotationalDofP = cP.rotationalDegreesOfFreedom();        
    // Larsen Borgnakke rotational energy redistribution part. Using the serial
    // application of the LB method, as per the INELRS subroutine in Bird's
    // DSMC0R.FOR
    if (rotationalDofP > 0)
    {
      /*
         scalar particleProbabilityP = 
             ((zeta_T + 2.0*rotationalDofP)/(2.0*rotationalDofP))
             *(
                 1.0 - sqrt(
                             1.0 - (rotationalDofP/zeta_T)
                             *((zeta_T+rotationalDofP)/(zeta_T+2.0*rotationalDofP))
                             *(4.0/rotationalRelaxationCollisionNumber_)
                           )
              );
            
         Info << "particleProbabilityP = " << particleProbabilityP << endl;
	 //if (particleProbabilityP > cloud_.rndGen().sample01<scalar>())*/

      scalar& ERotP = p.ERot();
      
      const scalar inverseRotationalCollisionNumber =
	1.0/(rotationalRelaxationCollisionNumber_);
      
      const scalar preCollisionERotP = ERotP;
      
      if (inverseRotationalCollisionNumber > cloud_.rndGen().sample01<scalar>())
      {
	  const scalar EcP = translationalEnergy + preCollisionERotP;
	  const scalar ChiB = 2.5 - omegaPQ;
          
	  const scalar energyRatio = 
	    cloud_.postCollisionRotationalEnergy(rotationalDofP, ChiB);
	  
	  ERotP = energyRatio*EcP;
	  
	  translationalEnergy = EcP - ERotP;
      }
      
    }
    
}

void Foam::LarsenBorgnakkeVariableHardSphere::redistributeOnlyOneMode
(
    dsmcParcel& p,
    dsmcParcel& q,
    scalar& translationalEnergy,
    const scalar omegaPQ
)
{
    const label typeIdP = p.typeId();
    const label typeIdQ = q.typeId();
    const dsmcParcel::constantProperties& cP = cloud_.constProps(typeIdP);
    const dsmcParcel::constantProperties& cQ = cloud_.constProps(typeIdQ);

    /*
    //calculate total P
    const scalar inverseRotationalCollisionNumber =
      1.0/(rotationalRelaxationCollisionNumber_);

    cost scalar inverseVibrationalCollisionNumberP =
      1.0/(vibrationalRelaxationCollisionNumber_);
    
    
    if(R < P)
    {
      
    }
    */

    //generate rndon variable R
    const scalar R = cloud_.rndGen().sample01<scalar>();

    scalar sumP    = 0.0;
    
    //- Rotational energy mode for P
    const scalar rotationalDofP = cP.rotationalDegreesOfFreedom();        
    // Larsen Borgnakke rotational energy redistribution part. Using the serial
    // application of the LB method, as per the INELRS subroutine in Bird's
    // DSMC0R.FOR
    if (rotationalDofP > 0)
    {

      scalar& ERotP = p.ERot();
      
      const scalar inverseRotationalCollisionNumber =
	1.0/(rotationalRelaxationCollisionNumber_);
      
      const scalar preCollisionERotP = ERotP;

      sumP += inverseRotationalCollisionNumber;
      //if (inverseRotationalCollisionNumber > cloud_.rndGen().sample01<scalar>())
      if (sumP > R)
      {
	  const scalar EcP = translationalEnergy + preCollisionERotP;
	  const scalar ChiB = 2.5 - omegaPQ;
          
	  const scalar energyRatio = 
	    cloud_.postCollisionRotationalEnergy(rotationalDofP, ChiB);
	  
	  ERotP = energyRatio*EcP;
	  
	  translationalEnergy = EcP - ERotP;

	  return;
      }  
    }

    //- Rotational energy mode for Q
    const scalar rotationalDofQ = cQ.rotationalDegreesOfFreedom();        
    // Larsen Borgnakke rotational energy redistribution part. Using the serial
    // application of the LB method, as per the INELRS subroutine in Bird's
    // DSMC0R.FOR
    if (rotationalDofQ > 0)
    {

      scalar& ERotQ = q.ERot();
      
      const scalar inverseRotationalCollisionNumber =
	1.0/(rotationalRelaxationCollisionNumber_);
      
      const scalar preCollisionERotQ = ERotQ;

      sumP += inverseRotationalCollisionNumber;
      //if (inverseRotationalCollisionNumber > cloud_.rndGen().sample01<scalar>())
      if (sumP > R)
      {
	  const scalar EcQ = translationalEnergy + preCollisionERotQ;
	  const scalar ChiB = 2.5 - omegaPQ;
          
	  const scalar energyRatio = 
	    cloud_.postCollisionRotationalEnergy(rotationalDofQ, ChiB);
	  
	  ERotQ = energyRatio*EcQ;
	  
	  translationalEnergy = EcQ - ERotQ;

	  return;
      }      
    }

    
    //- Vibrational energy mode for P
    if (cP.nVibrationalModes() > 0)
    {
        const scalarList& thetaVP = cP.thetaV();  
        const scalarList& thetaDP = cP.thetaD();
        const List<scalarList>& ZrefP = cP.Zref();
        const List<scalarList>& refTempZvP = cP.TrefZv();
        const scalarList& preCollisionEVibP = cP.eVib(p.vibLevel());
	label partnerIndex = findIndex(cP.ZrefPartnerId(), typeIdQ);
	if(partnerIndex == -1)
	{
	  partnerIndex = cP.ZrefPartnerMIndex();
	}

	forAll(thetaVP, i)
	{	
	  //- Collision energy of particle P: relative translational energy 
            //    + pre-collision vibrational energy
            const scalar EcP = translationalEnergy + preCollisionEVibP[i]; 

            //- Maximum possible quantum level (equation 3, Bird 2010)
            const label iMaxP = EcP/(physicoChemical::k.value()*thetaVP[i]);

	    const scalar correctZFactorP = cloud_.fields().correctZFactor(p.cell(), typeIdP, i);
	    
            if (iMaxP > 0)
            {
	      const scalar inverseVibrationalCollisionNumberP =
                    cloud_.postCollisionVibrationalEnergyLevelOneMode
                    (		       
                        iMaxP,
                        thetaVP[i],
                        thetaDP[i],
                        refTempZvP[partnerIndex][i],
                        omegaPQ,
                        ZrefP[partnerIndex][i],
                        0.0,//Ec not used
			correctZFactorP,
                        vibrationalRelaxationCollisionNumber_,
                        invZvFormulation_,
                        p.cell()		        
                    );	      	    
	      
	      sumP += inverseVibrationalCollisionNumberP;
	      if( sumP > R )
	      {
		// post-collision quantum number
		scalar func = 0.0;
		scalar EVib = 0.0;
		label iDash = 0;
		
		do // acceptance - rejection
		{
		  //iDash = rndGen_.position<label>(0, iMax); OLD
		  iDash = cloud_.randomLabel(0, iMaxP);
		  
		  EVib = iDash*physicoChemical::k.value()*thetaVP[i];
		  
		  // - equation 5.61, Bird
		  func = pow(1.0 - EVib/EcP, 1.5 - omegaPQ);
		  
		} while(func < cloud_.rndGen().sample01<scalar>());
		
		p.vibLevel()[i] = iDash;
	      
		translationalEnergy = EcP - cP.eVib_m(i, p.vibLevel()[i]);
		
		return;
	      }

	    }//end iMax
	    
	}//end for      
    }
    
    
    //- Vibrational energy mode for Q
    if (cQ.nVibrationalModes() > 0)
    {
      
        const scalarList& thetaVQ = cQ.thetaV();  
        const scalarList& thetaDQ = cQ.thetaD();
        const List<scalarList>& ZrefQ = cQ.Zref();
        const List<scalarList>& refTempZvQ = cQ.TrefZv();
        const scalarList& preCollisionEVibQ = cQ.eVib(q.vibLevel());
	label partnerIndex = findIndex(cQ.ZrefPartnerId(), typeIdP);
	if(partnerIndex == -1)
	{
	  partnerIndex = cQ.ZrefPartnerMIndex();
	}
		
	forAll(thetaVQ, i)
	{	
	  
	  //- Collision energy of particle Q: relative translational energy 
            //    + pre-collision vibrational energy
            const scalar EcQ  = translationalEnergy + preCollisionEVibQ[i]; 

            //- Maximum possible quantum level (equation 3, Bird 2010)
            const label iMaxQ = EcQ/(physicoChemical::k.value()*thetaVQ[i]);

	    const scalar correctZFactorQ = cloud_.fields().correctZFactor(q.cell(), typeIdQ, i);
	    
            if (iMaxQ > 0)
            {
	      const scalar inverseVibrationalCollisionNumberQ =
                    cloud_.postCollisionVibrationalEnergyLevelOneMode
                    (		       
                        iMaxQ,
                        thetaVQ[i],
                        thetaDQ[i],
                        refTempZvQ[partnerIndex][i],
                        omegaPQ,
                        ZrefQ[partnerIndex][i],
                        0.0,//Ec not used
			correctZFactorQ,
                        vibrationalRelaxationCollisionNumber_,
                        invZvFormulation_,
                        q.cell()		        
                    );

	      sumP += inverseVibrationalCollisionNumberQ;
	      
	      if( sumP > R)
	      {		
		// post-collision quantum number
		scalar func = 0.0;
		scalar EVib = 0.0;
		label iDash = 0;
		
		do // acceptance - rejection
		{
		  //iDash = rndGen_.position<label>(0, iMax); OLD
		  iDash = cloud_.randomLabel(0, iMaxQ);
	      
		  EVib = iDash*physicoChemical::k.value()*thetaVQ[i];
		  
		  // - equation 5.61, Bird
		  func = pow(1.0 - EVib/EcQ, 1.5 - omegaPQ);
		  
		} while(func < cloud_.rndGen().sample01<scalar>());
		
		q.vibLevel()[i] = iDash;
		
		translationalEnergy = EcQ - cQ.eVib_m(i, q.vibLevel()[i]);

		 return;
	      }	    
	    }	  
	} 

	/*
	forAll(thetaVQ, i)
	{
            //- Collision energy of particle Q: relative translational energy 
            //    + pre-collision vibrational energy
            const scalar EcQ = translationalEnergy + preCollisionEVibQ[i]; 

            //- Maximum possible quantum level (equation 3, Bird 2010)
            const label iMaxQ = EcQ/(physicoChemical::k.value()*thetaVQ[i]);

            if (iMaxQ > 0)
            {	      
	      const scalar inverseVibrationalCollisionNumberQ =
                    cloud_.postCollisionVibrationalEnergyLevelOneMode
                    (		       
                        iMaxQ,
                        thetaVQ[i],
                        thetaDQ[i],
                        refTempZvQ[partnerIndex][i],
                        omegaPQ,
                        ZrefQ[partnerIndex][i],
                        EcQ,			
                        vibrationalRelaxationCollisionNumber_,
                        invZvFormulation_,
                        q.cell()		        
                    );

	      label iDash = 0;
	      
	      sumP += inverseVibrationalCollisionNumberQ;		    
	      if(sumP > R)		      
	      {
		// post-collision quantum number
		scalar func = 0.0;
		scalar EVib = 0.0;
		      
		do // acceptance - rejection
		{
		  //iDash = rndGen_.position<label>(0, iMax); OLD
		  iDash = cloud_.randomLabel(0, iMaxQ);
		  
		  EVib = iDash*physicoChemical::k.value()*thetaVQ[i];
		  
		  // - equation 5.61, Bird
		  func = pow(1.0 - EVib/EcQ, 1.5 - omegaPQ);
		  
		} while(func < cloud_.rndGen().sample01<scalar>());

		q.vibLevel()[i] = iDash;

		translationalEnergy = EcQ - cQ.eVib_m(i, q.vibLevel()[i]);

		return;
	      }
            }
       }
	*/
	
    }

    
}

const Foam::dictionary&
Foam::LarsenBorgnakkeVariableHardSphere::coeffDict() const
{
    return coeffDictLB_;
}

// ************************************************************************* //
