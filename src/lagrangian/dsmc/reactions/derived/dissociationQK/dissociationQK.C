/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description

\*---------------------------------------------------------------------------*/

#include "dissociationQK.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dissociationQK, 0);

addToRunTimeSelectionTable(dsmcReaction, dissociationQK, dictionary);


// * * * * * * * * * * *  Protected Member functions * * * * * * * * * * * * //

void dissociationQK::setProperties()
{
    dsmcReaction::setProperties();
    
    if (reactantIds_.size() != 2)
    {
        //- There must be exactly 2 reactants
        FatalErrorIn("dissociationQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "There should be two reactants, instead of " 
            << reactantIds_.size() << nl 
            << exit(FatalError);
    }

    dissociationStr_.setSize(reactantIds_.size(), word::null);
    nTotDissociationReactions_.setSize(reactantIds_.size(), 0);
    nDissociationReactionsPerTimeStep_.setSize(reactantIds_.size(), 0);
    heatOfReactionDissociationJoules_.setSize(reactantIds_.size());
    
    bool moleculeFound = false;
    
    forAll(reactantIds_, r)
    {
        //- Check if this reactant is a neutal molecule
        if (reactantTypes_[r] == 20 or reactantTypes_[r] == 30)
        {
            moleculeFound = true;
            
            heatOfReactionDissociationJoules_[r] = physicoChemical::k.value()
                *cloud_.constProps(reactantIds_[r]).thetaD()[0];
        }
    }
    
    if (!moleculeFound)
    {
        FatalErrorIn("dissociationQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "None of the reactants is a molecule." << nl 
            << exit(FatalError);
    }
    
    //- Reading in dissociation products
    const List<wordList> productsDissociation(propsDict_.lookup("dissociationProducts"));

    if (productsDissociation.size() != 2)
    {
        FatalErrorIn("dissociationQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "There should be two lists of products, instead of " 
            << productsDissociation.size() << nl 
            << "NB: a list can be left empty" << nl
            << exit(FatalError);
    }
    
    productIdsDissociation_.setSize(productsDissociation.size());

    forAll(productIdsDissociation_, r)
    {
        if (productsDissociation[r].size() > 0)
        {
            //- Check that there is no dissociation products should the reactant
            //  not be a molecule
            if (reactantTypes_[r] != 20 and reactantTypes_[r] != 30)
            {
                FatalErrorIn("dissociationQK::setProperties()")
                    << "For reaction named " << reactionName_ << nl
                    << "Reactant " << cloud_.typeIdList()[reactantIds_[r]]
                    << " is not a molecule " << nl 
                    << "and therefore, there should be no dissociation products "
                    << "instead of " << productsDissociation[r]
                    << exit(FatalError);
            }
            
            //- Check that there are two products
            if (productsDissociation[r].size() != 2)
            {
                FatalErrorIn("dissociationQK::setProperties()")
                    << "For reaction named " << reactionName_ << nl
                    << "There should be 2 dissociation products for molecule " 
                    << cloud_.typeIdList()[reactantIds_[r]] << " instead of " 
                    << productsDissociation[r].size() << ", that is "
                    << productsDissociation[r]
                    << exit(FatalError);
            }
        }
        else
        {
            //- If this species is a molecule, it should have dissociation
            //  products
            if (reactantTypes_[r] == 20 or reactantTypes_[r] == 30)
            {
                FatalErrorIn("dissociationQK::setProperties()")
                    << "For reaction named " << reactionName_ << nl
                    << "Reactant " << cloud_.typeIdList()[reactantIds_[r]]
                    << " is a molecule " << nl 
                    << "and therefore, it should have dissociation products "
                    << "instead of " << productsDissociation[r]
                    << exit(FatalError);
            }
        }
            
        productIdsDissociation_[r].setSize(productsDissociation[r].size());
        
        forAll(productIdsDissociation_[r], p)
        {
            productIdsDissociation_[r][p] = 
                findIndex
                (
                    cloud_.typeIdList(), 
                    productsDissociation[r][p]
                );

            //- Check that products belong to the typeIdList as defined in 
            //  constant/dsmcProperties
            if (productIdsDissociation_[r][p] == -1)
            {
                FatalErrorIn("dissociationQK::setProperties()")
                    << "For reaction named " << reactionName_ << nl
                    << "Cannot find type id: " << productsDissociation[r][p] << nl 
                    << exit(FatalError);
            }

            //- Check that products of diatomic particle are atoms
            if (
                    reactantTypes_[r] == 20 
                 && cloud_.constProps(productIdsDissociation_[r][p]).type() != 10
               )
            {
                FatalErrorIn("dissociationQK::setProperties()")
                    << "For reaction named " << reactionName_ << nl
                    << "Dissociation product of a diatomic molecule must be an atom: " 
                    << productsDissociation[r][p] << nl 
                    << exit(FatalError);
            }

	    /*
            //- Check that products of polyatomic particle are 
            //  diatomic/polyatomic molecules
            if (
                    reactantTypes_[r] == 30 
                 && (cloud_.constProps(productIdsDissociation_[r][p]).type() != 20
                     or cloud_.constProps(productIdsDissociation_[r][p]).type() != 30)
               )
            {
                FatalErrorIn("dissociationQK::setProperties()")
                    << "For reaction named " << reactionName_ << nl
                    << "Dissociation product of a polyatomic molecule must be a "
                    << "diatomic/polyatomic molecule instead of " 
                    << productsDissociation[r][p] << nl 
                    << exit(FatalError);
            }
	    */
        }
    }
}


void dissociationQK::testDissociation
(
    const dsmcParcel& p,
    const scalar translationalEnergy,
    label& vibModeDisso,
    scalar& reactionProbability
)
{
    const label typeIdP = p.typeId();
    
    if( cloud_.constProps(typeIdP).type() >= 20 )
    {
        const labelList& vibLevelP = p.vibLevel();
	
	//random shuffle different vib mode
	//const scalarList equelProbabilities(vibLevelP.size(),0.0);
	//const labelList  randomShuffleIndices = decreasing_sort_indices(equelProbabilities);
	
        //forAll(vibLevelP, m)	  
	//{
	//const label  indice  = randomShuffleIndices[m];
	    const label   indice = cloud_.randomLabel(0, vibLevelP.size()-1);
            const scalar thetaVP = cloud_.constProps(typeIdP).thetaV_m(indice);
            const scalar EVibP_m = cloud_.constProps(typeIdP).eVib_m(indice, vibLevelP[indice]);
            const label  idP = cloud_.constProps(typeIdP).charDissQuantumLevel_m(indice);    
            const scalar collisionEnergy = translationalEnergy + EVibP_m;
            const label  imaxP = collisionEnergy/(physicoChemical::k.value()*thetaVP);
	    
            //- Condition for the dissociation of the molecule P, mode m
            if (imaxP > idP)
	    //if(collisionEnergy > heatOfReactionDissociationJoules_[1] )
            {	      
                //- Add reaction to the list of competing reactions with 
                //  probability reactionProbability
                reactionProbability = 1.0;
                //- The molecule dissociates in the vibrational energy mode
                vibModeDisso = indice;
		
		//break;
            }
	//}	
    }
}

void dissociationQK::dissociateParticleByPartner
(
    dsmcParcel& p,
    dsmcParcel& q,
    const label nR,
    const label vibModeDisso,
    const scalar translationalEnergy
)
{
    const label typeIdP = p.typeId();
    const label typeIdQ = q.typeId();
    
    label nReac = nR;
    if (typeIdP == typeIdQ)
    {
        nReac = 0;
    }
    
    //- Dissociation of parcel p !!!!!
    nTotDissociationReactions_[nReac]++;
    nDissociationReactionsPerTimeStep_[nReac]++;

    //return;
    
    if (allowSplitting_)
    {
      //if reaction happen relax is true
      relax_ = false;

      //temp
      //relax_ = true;
      
	const scalar EVibP_disso = cloud_.constProps(typeIdP).eVib_m
	(
	 vibModeDisso,
	 p.vibLevel()[vibModeDisso]
	);
	
	
        const scalar omegaPQ =
            0.5
            *(
                  cloud_.constProps(typeIdP).omega()
                + cloud_.constProps(typeIdQ).omega()
            );
	
	const scalar ktheta = physicoChemical::k.value()*cloud_.constProps(typeIdP).thetaV_m(vibModeDisso);
	
	//- The collision energy is being subtracted the heat of reaction
        scalar collisionEnergy = translationalEnergy + EVibP_disso
	  // - heatOfReactionDissociationJoules_[nReac][vibModeDisso];
	  - heatOfReactionDissociationJoules_[nReac];
	
	const label iMax = collisionEnergy/ktheta;
	label  j = 0;	  	  
	scalar func  = 0.0;

	if(iMax > 0)
	{  
	  do // acceptance - rejection
	  {
	    j = cloud_.randomLabel(0, iMax);
	      
	    func =
	      pow(1.0 - j*ktheta/collisionEnergy ,1.5-omegaPQ );
		
	  }while(func < cloud_.rndGen().sample01<scalar>() );	    
	}

	collisionEnergy -= j*ktheta;
	
	const scalar postExcitVibE = j*ktheta;
	
	/* //not good the vibrational distributtion will not become boltzman
        //- Energy redistribution for particle Q
        cloud_.binaryCollision().redistribute
        (
            q, collisionEnergy, omegaPQ, true
        );       
	*/
	
	//const scalar& dofP = cloud_.constProps(typeIdP).rotationalDegreesOfFreedom();
        const scalar mP = cloud_.constProps(typeIdP).mass();
        const scalar mQ = cloud_.constProps(typeIdQ).mass();
        const scalar mR = mP*mQ/(mP + mQ);
        scalar relVelNonDissoParticle = sqrt(2.0*collisionEnergy/mR);

        cloud_.binaryCollision().postCollisionVelocities
        (
            typeIdP,
            typeIdQ,
            p.U(),
            q.U(),
            relVelNonDissoParticle
        );	

	//- Energy left for the 2 products
        const scalar ERotP = p.ERot();
	const scalar EVibP_tot =
	  cloud_.constProps(typeIdP).eVib_tot
	  (
	   p.vibLevel()
	   );
	
	//remain energy
	scalar translationalEnergyLeft = ERotP + EVibP_tot - EVibP_disso + postExcitVibE;
	
	///////////////////////handle cleave partile internal energy///////////////	
        //- Post-reaction
        const label typeId1 = productIdsDissociation_[nReac][0];
        const label typeId2 = productIdsDissociation_[nReac][1];
   
        //- Mass of products 1 and 2
        const scalar mP1 = cloud_.constProps(typeId1).mass();
        const scalar mP2 = cloud_.constProps(typeId2).mass();
        const scalar mRproducts = mP1*mP2/(mP1 + mP2);

	
	const scalarList& thetaVProductMole   = cloud_.constProps(typeId1).thetaV();
	const scalarList& thetaVProductSecond = cloud_.constProps(typeId2).thetaV();
	const scalar& moleRDof   = cloud_.constProps(typeId1).rotationalDegreesOfFreedom();
	const scalar& secondRDof = cloud_.constProps(typeId2).rotationalDegreesOfFreedom();
	
	const scalar reverseOmega =  0.5*(
					  cloud_.constProps(typeId1).omega()
					  + cloud_.constProps(typeId2).omega()
					 );
	
	
	labelList vibLevelMole(thetaVProductMole.size(), 0);
	labelList vibLevelSecond(thetaVProductSecond.size(), 0);
	const label numberVibMode = vibLevelMole.size() + vibLevelSecond.size();

	/*
	DynamicList<label> totalVibLevelList(vibLevelMole);
	totalVibLevelList.append(vibLevelSecond);
	DynamicList<scalar> totalThetaVProductList(thetaVProductMole);
	totalThetaVProductList.append(thetaVProductSecond);
	
	/////////////////prepare nonExcite vibrational energy list from particle P//////////
	DynamicList<scalar> nonExciteVibE(0);   
	forAll(p.vibLevel(),m)
	{
	  if(m != vibModeDisso)
	  {
	    nonExciteVibE.append( cloud_.constProps(typeIdP).eVib_m
				  (
				   m,
				   p.vibLevel()[m]
				   )
				 );
	    
	    if(nonExciteVibE.size() == numberVibMode)	      
	      break;
	  }
	}
	*/

	//original
	//scalar     remainDOF  = dofP + 2.0*(p.vibLevel().size()-1);

	//target
	scalar     remainDOF  = moleRDof + secondRDof + 2.0*(numberVibMode + 2.5-reverseOmega);

	label mode = 0;
	scalarList theta    =  thetaVProductMole;//assume
	labelList* vibLevel = &vibLevelMole;//assume
	
	//forAll(totalThetaVProductList ,m)
	for(label m=0; m<numberVibMode; m++)
	{
	    mode = m;
	    remainDOF -= 2.0;
	    
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
	    
	    cloud_.postReactionVibrationalRedistribution
	    (
	     mode,
	     remainDOF,
	     theta,
	     vibLevel,
	     translationalEnergyLeft
	     );	    
	}        
        
	/*
	//combine rotE and other not excite vibE with L-B
	scalar ERotPAddVibE = ERotP;
	forAll(totalVibLevelList, c)
	{
	  ERotPAddVibE += nonExciteVibE[c];  

	  const scalar kBByThetaVP = physicoChemical::k.value()*totalThetaVProductList[c];
	  const label iMax = ERotPAddVibE/kBByThetaVP;
	  label  j = 0;	  	  
	  scalar func  = 0.0;

	  if(iMax > 0)
	  {
	    if(dofP == 2.0)
	    {
	      j = cloud_.randomLabel(0, iMax);
	    }
	    else
	    {
	      const scalar rotDofMinusOne = dofP/2.0-1.0;
	      //select postProuct pre-reaction vibrational level	    
	      do // acceptance - rejection
	      {
		j = cloud_.randomLabel(0, iMax);
		
		func =
		  pow(1.0 - j*kBByThetaVP/ERotPAddVibE ,rotDofMinusOne );
		
	      }while(func < cloud_.rndGen().sample01<scalar>() );	    
	    }
	  }

	  if(c < vibLevelMole.size())
	  {
	    vibLevelMole[c] = j;
	  }
	  else
	  {
	    vibLevelSecond[c-vibLevelMole.size()] = j;
	  }
	  
	  ERotPAddVibE -= j*kBByThetaVP;
	  translationalEnergyLeft -= j*kBByThetaVP;
	}
	*/

	/*
	if(vibLevelMole.size() > 0){
	  Info << "product = " << vibLevelMole[0]<< endl;
	}

	return;
	*/
	
	////////////////-first Trial L-B redistribution (rotation)
	scalar ERotProductMole   = 0.0;
	if( thetaVProductMole.size() > 0 )
	{
	  //const scalar remainDofbyTwo  = (dofP-moleRDof)/2.0+(p.vibLevel().size()-numberVibMode);
	  remainDOF -= moleRDof; 
	  const scalar energyRatioMole = cloud_.postCollisionRotationalEnergy( moleRDof, remainDOF/2.0 );	  
	  ERotProductMole   = energyRatioMole*translationalEnergyLeft;
	  translationalEnergyLeft -= ERotProductMole;//- Relative collisionEnergy energy after rotational energy redistribution	  
	}

	////////////////-second Trial L-B redistribution (rotation)
	scalar ERotProductSecond = 0.0;
	if( thetaVProductSecond.size() > 0 )
	{
	  //const scalar remainDofbyTwo  = (dofP-moleRDof-secondRDof)/2.0+(p.vibLevel().size()-numberVibMode);
	  remainDOF -= secondRDof; 
	  const scalar energyRatioSecond = cloud_.postCollisionRotationalEnergy( secondRDof, remainDOF/2.0 );
	  ERotProductSecond = energyRatioSecond*translationalEnergyLeft;
	  translationalEnergyLeft -= ERotProductSecond;//- Relative collisionEnergy energy after rotational energy redistribution
	}

        const scalar cRproducts = sqrt(2.0*translationalEnergyLeft/mRproducts);
	vector UP2 = vector::zero;
	cloud_.binaryCollision().postReactionVelocities
	(
	 typeId1,
	 typeId2,
	 p.U(),
	 UP2,
	 cRproducts
	 );	

	//Info << "Pbefore(dissociation) = " << p.typeId() << endl;
	//Info << "Qbefore(dissociation) = " << q.typeId() << endl;
	
        // Molecule P dissociates
        // The molecule is first transformed into the first product ...
        p.typeId()   = typeId1;
        p.vibLevel() = vibLevelMole;
        p.ERot()     = ERotProductMole;
        p.ELevel()   = 0;
		
        // ... and a second product is then inserted in the same cell
        const label  product2_ELevel = 0;
        cloud_.addNewParcel
        (
            p.position(),
            UP2,
            p.RWF(),
            ERotProductSecond,
            product2_ELevel,
            p.cell(),
            p.tetFace(),
            p.tetPt(),
            typeId2,
            -1,
            p.classification(),
            vibLevelSecond
        );

	//Info << "Pafter(dissociation) = " << p.typeId() << endl;
	//Info << "Qafter(dissociation) = " << q.typeId() << endl;
	
	/*
	scalar postE =
	   0.5*cloud_.constProps(typeId2).mass()*magSqr(UP2)
	  +0.5*cloud_.constProps(p.typeId()).mass()*magSqr(p.U())
	  +0.5*cloud_.constProps(q.typeId()).mass()*magSqr(q.U())
	  +p.ERot()+q.ERot()+ERotProductSecond
	  +cloud_.constProps(p.typeId()).eVib_tot(p.vibLevel())
	  +cloud_.constProps(q.typeId()).eVib_tot(q.vibLevel())
	  +cloud_.constProps(typeId2).eVib_tot(vibLevelSecond);
	Info << "postE = " << postE << endl;
	*/
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dissociationQK::dissociationQK
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    productIdsDissociation_(),
    dissociationStr_(),
    nTotDissociationReactions_(),
    nDissociationReactionsPerTimeStep_(),
    heatOfReactionDissociationJoules_(),
    volume_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dissociationQK::~dissociationQK()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dissociationQK::initialConfiguration()
{
    setProperties();
    
    const word& reactantA = cloud_.typeIdList()[reactantIds_[0]];
    const word& reactantB = cloud_.typeIdList()[reactantIds_[1]];
    
    if (productIdsDissociation_[0].size())
    {
        const word& productA = cloud_.typeIdList()[productIdsDissociation_[0][0]];
        const word& productB = cloud_.typeIdList()[productIdsDissociation_[0][1]];
        dissociationStr_[0] = "Dissociation reaction " + reactantA + " + " 
            + reactantB + " --> " + productA + " + " + productB + " + " 
            + reactantB;
    }
    if (productIdsDissociation_[1].size() and reactantA != reactantB)
    {
        const word& productC = cloud_.typeIdList()[productIdsDissociation_[1][0]];
        const word& productD = cloud_.typeIdList()[productIdsDissociation_[1][1]];
                
        dissociationStr_[1] = "Dissociation reaction " + reactantB + " + " 
            + reactantA + " --> " + productC + " + " + productD + " + " 
            + reactantA;    
    } 
}


bool dissociationQK::tryReactMolecules
(
    const label& typeIdP,
    const label& typeIdQ
)
{
  //- Function used when setting the pair addressing matrix
    const label reactantPId = findIndex(reactantIds_, typeIdP);
    const label reactantQId = findIndex(reactantIds_, typeIdQ);

    //- If both indices were found in the list of reactants, there Ids will be
    //  different from -1
    if ((reactantPId != -1) && (reactantQId != -1))
    {
        //- Case of similar species colliding
        if((reactantPId == reactantQId) and (reactantIds_[0] == reactantIds_[1]))
        {
            return true;
        }
        
        //- Case of dissimilar species colliding
        if((reactantPId != reactantQId) and (reactantIds_[0] != reactantIds_[1]))
        {
            return true;
        }
    }
        
    return false;
}


void dissociationQK::reaction
(
    dsmcParcel& p,
    dsmcParcel& q,
    dsmcParcel& thirdBody
    //const label candidateId
    //const DynamicList<label>& candidateList,
    //const List<DynamicList<label> >& candidateSubList,
    //const label& candidateP,
    //const List<label>& whichSubCell
)
{}


void dissociationQK::reaction(dsmcParcel& p, dsmcParcel& q)
{
    //- Reset the relax switch
    relax_ = true;
    
    const label typeIdP = p.typeId();
    const label typeIdQ = q.typeId();
 
    //- Dissociation reaction AB + M --> A + B + M 
    //  If P is the first reactant AB
    //  NB: Q is necessarily M otherwise this class would not have been selected
    if (typeIdP == reactantIds_[0]) 
    {
      /*
      scalar preE = 0.5*cloud_.constProps(p.typeId()).mass()*magSqr(p.U())
	+0.5*cloud_.constProps(q.typeId()).mass()*magSqr(q.U())
	+p.ERot()+q.ERot()
	+cloud_.constProps(p.typeId()).eVib_tot(p.vibLevel())
	+cloud_.constProps(q.typeId()).eVib_tot(q.vibLevel())
	- heatOfReactionDissociationJoules_[postheat position];
      Info << "preE = " << preE << endl;
      */
      
        const scalar mP = cloud_.constProps(typeIdP).mass();
        const scalar mQ = cloud_.constProps(typeIdQ).mass();
        const scalar mR = mP*mQ/(mP + mQ);
        
        const scalar cRsqr = magSqr(p.U() - q.U());
        const scalar translationalEnergy = 0.5*mR*cRsqr;

	label vibModeDissoP = -1;
        label vibModeDissoQ = -1;

	const label& pType = cloud_.constProps(typeIdP).type();
	const label& qType = cloud_.constProps(typeIdQ).type();
	
	label bothDiatomic = 0;
	if( pType >= 20 && qType >= 20 )
	{
	  bothDiatomic = 1;
	}
        
        //- Possible reactions:
        // 1. Dissociation of P
        // 2. Dissociation of Q
        
        //scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(2, 0.0);    
	
	//only one reaction will be consider
	const label rn = cloud_.randomLabel(0, bothDiatomic);//2dissociation
	
	if( rn == 0)
	{
	  testDissociation
	  (
	   p,
	   translationalEnergy,
	   vibModeDissoP,
	   reactionProbabilities[0]
	  );

	  if(reactionProbabilities[0] > cloud_.rndGen().sample01<scalar>())
	  {
	    dissociateParticleByPartner
	    (
	     p, q, 0, vibModeDissoP, translationalEnergy
	    );
	  }
	  
	  //totalReactionProbability += reactionProbabilities[rn];	
	}
	else if(rn == bothDiatomic)
	{
	  testDissociation
	  (
	   q,
	   translationalEnergy,
	   vibModeDissoQ,
	   reactionProbabilities[1]
	  );	        	    

	  if(reactionProbabilities[1] > cloud_.rndGen().sample01<scalar>())
	  {
	    dissociateParticleByPartner
	    (
	     q, p, 1, vibModeDissoQ, translationalEnergy
	    );
	  }
	  
	  //totalReactionProbability += reactionProbabilities[rn];
	}
	  
	/*
        //- Decide if a reaction is to occur
        if (totalReactionProbability > cloud_.rndGen().sample01<scalar>())
        {
            //- A chemical reaction is to occur, normalise probabilities
            const scalarList normalisedProbabilities =
	      reactionProbabilities/totalReactionProbability;

            //- Sort normalised probability indices in decreasing order
            //  for identical probabilities, random shuffle
            const labelList sortedNormalisedProbabilityIndices =
                decreasing_sort_indices(normalisedProbabilities);
            scalar cumulativeProbability = 0.0;
    
            forAll(sortedNormalisedProbabilityIndices, idx)
            {                
                const label i = sortedNormalisedProbabilityIndices[idx];

                //- If current reaction can't occur, end the search
                if (normalisedProbabilities[i] > SMALL)
		{
                    cumulativeProbability += normalisedProbabilities[i];
                    
                    if (cumulativeProbability > cloud_.rndGen().sample01<scalar>())
                    {
                        //- Current reaction is to occur
                        if (i == 0)
                        {
	*/
			  /*
			  scalar noOfKT =
			    (translationalEnergy+p.vibLevel()[0]*physicoChemical::k.value()*cloud_.constProps(typeIdP).thetaV()[0])
			    //translationalEnergy
			    /(0.05*physicoChemical::k.value()*4000.0);
			  
			scalar colliDifference = noOfKT-int(noOfKT);
			if(colliDifference < 0.1 )
			  {			    
			    Info << "EcPre = "
				 << int(noOfKT)
				 << endl;
			  }
			else if( colliDifference >  0.9 )
			  {	    
			    Info << "EcPre = "
				 << int(noOfKT)+1
				 << endl;
				 }
			  */		          

			  /*
			  Info << "heat = " 
			       << p.vibLevel()[0] << " "
			       << p.vibLevel()[1] << " "
			       << p.vibLevel()[2] << endl;
			  */
			  /*
			  Info << "thirdBodyy = " 
			       << q.vibLevel()[0] << endl;
			  */

	/*
			  //- Dissociation of P is to occur
			  dissociateParticleByPartner
			  (
			       p, q, i, vibModeDissoP, translationalEnergy
			  );
			  
			  //- There can't be another reaction: break
			  break;
                        }

			
                        if (i == 1)
                        {
	*/
			  /*
			  if(vibModeDissoQ == 0)
			  {
			  scalar noOfKT =
			    (translationalEnergy+q.vibLevel()[0]*physicoChemical::k.value()*cloud_.constProps(typeIdQ).thetaV()[0])
			    //translationalEnergy
			    /(0.05*physicoChemical::k.value()*4000.0);
			  
			   scalar colliDifference = noOfKT-int(noOfKT);
			   if(colliDifference < 0.1 )
			   {			    
			     Info << "EcPre = "
				  << int(noOfKT)
				  << endl;
			   }
			   else if( colliDifference >  0.9 )
			   {	    
			     Info << "EcPre = "
				  << int(noOfKT)+1
				  << endl;
			   }			  
			  } 
			  */

			  /*
			  Info << "cold = "
			       << q.vibLevel()[0] << " "
			       << q.vibLevel()[1] << " "
			       << q.vibLevel()[2]			  
			       << endl;	
			  */
			  
			  //Info << "thirdBodyy = " 
			  //     << p.vibLevel()[0] << endl;

	/*
			  //- Dissociation of Q is to occur
			  dissociateParticleByPartner
			  (
			     q, p, i, vibModeDissoQ, translationalEnergy
			  );
			  
			  //- There can't be another reaction: break
                          break;
                        }
			
                    }
		}
                else
                {
                    //- All the following possible reactions have a probability
                    //  of zero
                    break;
                }
            }
        }
	*/
    }
    else
    {
        //- Dissociation reaction AB + M --> A + B + M 
        //  If P is the second reactant M, then switch arguments in this
        //  function and P will be first
        dissociationQK::reaction(q, p);
    }
}

void dissociationQK::outputResults(const label& counterIndex)
{
    if (writeRatesToTerminal_)
    {
        //- measure density 
        const List<DynamicList<dsmcParcel*>>& cellOccupancy = cloud_.cellOccupancy();
            
        volume_ = 0.0;

        scalarList molsReactants(2, 0.0);

        forAll(cellOccupancy, c)
        {
            const List<dsmcParcel*>& parcelsInCell = cellOccupancy[c];

            forAll(parcelsInCell, pIC)
            {
                dsmcParcel* p = parcelsInCell[pIC];
                
                const label pos = findIndex(reactantIds_, p->typeId());

                if (pos != -1)
                {
                    molsReactants[pos]++;
                }
            }

            volume_ += mesh_.cellVolumes()[c];
        }
        
        scalar volume = volume_;
        if (Pstream::parRun())
        {
            reduce(volume, sumOp<scalar>());
        }
        
        scalarList numberDensities(2, cloud_.nParticle()/volume);
        numberDensities[0] *= molsReactants[0];
        numberDensities[1] *= molsReactants[1];
        
        labelList nTotDissociationReactions = nTotDissociationReactions_;
        labelList nDissociationReactionsPerTimeStep = nDissociationReactionsPerTimeStep_;

        const scalar deltaT = mesh_.time().deltaT().value();
        scalar factor = 0.0;
        
        if (reactantIds_[0] == reactantIds_[1] && numberDensities[0] > 0.0)
        {
            factor = cloud_.nParticle()/
                (
                    counterIndex*deltaT
                   *numberDensities[0]*numberDensities[0]
                   *volume*2.0
                );
        }
        else if (numberDensities[0] > 0.0 && numberDensities[1] > 0.0)
        {
            factor = cloud_.nParticle()/
                (
                    counterIndex*deltaT
                   *numberDensities[0]*numberDensities[1]
                   *volume
		 );
        }
        
        for (label k=0; k<2; k++)
        {   
            if (dissociationStr_[k].size())
            {
                if (Pstream::parRun())
                {
                    //- Parallel communication
                    reduce(molsReactants[k], sumOp<label>());
                    reduce(nTotDissociationReactions[k], sumOp<label>());
                    reduce(nDissociationReactionsPerTimeStep[k], sumOp<label>());
                }
                
                const scalar reactionRateDissociation = factor*nTotDissociationReactions[k];
                
                Info<< dissociationStr_[k] 
                    << ", reaction rate = " << reactionRateDissociation
                    << ", nReactions = " << nDissociationReactionsPerTimeStep[k]
                    << endl;
            }    
        }
    }
    /*
    else
    {
        labelList nTotDissociationReactions = nTotDissociationReactions_;   
        labelList nDissociationReactionsPerTimeStep = nDissociationReactionsPerTimeStep_;
        
        for (label k=0; k<2; k++)
        {   
            if (dissociationStr_[k].size())
            {
                if (Pstream::parRun())
                {
                    //- Parallel communication
                    reduce(nTotDissociationReactions[k], sumOp<label>());
                    reduce(nDissociationReactionsPerTimeStep[k], sumOp<label>());
                }
                
                if (nTotDissociationReactions[k] > 0)
                {
                    Info<< dissociationStr_[k] 
                        << " is active, nReactions this time step = " 
                        << nDissociationReactionsPerTimeStep[k] 
                        << endl;
                 }
             }
        }
    }
    */

    nDissociationReactionsPerTimeStep_ = 0;
}

}
// End namespace Foam

// ************************************************************************* //
