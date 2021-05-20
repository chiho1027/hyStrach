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

#include "dissociationExchangeQK.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dissociationExchangeQK, 0);

addToRunTimeSelectionTable
(
    dsmcReaction,
    dissociationExchangeQK,
    dictionary
);


// * * * * * * * * * * *  Protected Member functions * * * * * * * * * * * * //

void dissociationExchangeQK::setProperties()
{
    dissociationQK::setProperties();
    exchangeQK::setProperties();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dissociationExchangeQK::dissociationExchangeQK
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    dissociationQK(t, cloud, dict),
    exchangeQK(t, cloud, dict)
    //propsDict_(dict.subDict(typeName + "Properties")),
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dissociationExchangeQK::~dissociationExchangeQK()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dissociationExchangeQK::initialConfiguration()
{
    dissociationQK::initialConfiguration();
    exchangeQK::initialConfiguration();
}


bool dissociationExchangeQK::tryReactMolecules
(
    const label& typeIdP,
    const label& typeIdQ
)
{
    return dissociationQK::tryReactMolecules(typeIdP, typeIdQ);
}


void dissociationExchangeQK::reaction
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


void dissociationExchangeQK::reaction(dsmcParcel& p, dsmcParcel& q)
{
    //- Reset the relax switch
    relax_ = true;
    
    const label typeIdP = p.typeId();
    const label typeIdQ = q.typeId();
    
    if (typeIdP == reactantIds_[0]) 
    { 
        const scalar mP = cloud_.constProps(typeIdP).mass();
        const scalar mQ = cloud_.constProps(typeIdQ).mass();
        const scalar mR = mP*mQ/(mP + mQ);
        
        const scalar omegaPQ =
            0.5
            *(
                  cloud_.constProps(typeIdP).omega()
                + cloud_.constProps(typeIdQ).omega()
            );
        
        const scalar cRsqr = magSqr(p.U() - q.U());
        const scalar translationalEnergy = 0.5*mR*cRsqr;

	// dissociation initialized data
        label vibModeDissoP = -1;
        label vibModeDissoQ = -1;

	// exchange initialized data
	const label numberOfExchange = exchangeQK::numberOfExchange_;
	
	// record each reaction probability of different vibrational mode
	List<DynamicList<scalar> > pReactionPDiffVibMode(numberOfExchange);
	List<DynamicList<scalar> > qReactionPDiffVibMode(numberOfExchange);

	// dissociate initalized data
	const label& pType = cloud_.constProps(typeIdP).type();
	const label& qType = cloud_.constProps(typeIdQ).type();
	
	//just one reaction will consider
	label bothDiatomic = 0;
	if( pType >= 20 && qType >= 20 )
	{
	  bothDiatomic = 1;
	}
	
        //- Possible reactions:
        // 1. Dissociation of P
        // 2. Dissociation of Q
        // 3. Exchange
        
        //scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(2+numberOfExchange, 0.0);
	
	//numberofexchange + 2dissociation or numberofexchange + dissociation
	const label rn = cloud_.randomLabel(0, numberOfExchange + bothDiatomic );
	if (rn == 0)
	{
	  dissociationQK::testDissociation
	  (
	   p,
	   translationalEnergy,
	   vibModeDissoP,
	   reactionProbabilities[0]
	  );

	  if(reactionProbabilities[0] > cloud_.rndGen().sample01<scalar>())
	  {
	    dissociationQK::dissociateParticleByPartner
	    (
	     p, q, 0, vibModeDissoP, translationalEnergy
	    );
	  }	  
	}
	else if (rn == bothDiatomic)
	{
	  dissociationQK::testDissociation
	  (
	   q,
	   translationalEnergy,
	   vibModeDissoQ,
	   reactionProbabilities[1]
	  );

	  if(reactionProbabilities[1] > cloud_.rndGen().sample01<scalar>())
	  {
	    dissociationQK::dissociateParticleByPartner
	    (
	     q, p, 1, vibModeDissoQ, translationalEnergy
	    );
	  }	  	    
	}
	else
	{
	  if (exchangeQK::posAtomReactant_ == 1) //exchange exsis atom ABC + D in chemicalDict
	  {
	    exchangeQK::testExchange
	    (
	     p,
	     translationalEnergy,
	     omegaPQ,
	     rn-1-bothDiatomic,
	     reactionProbabilities[rn],
	     pReactionPDiffVibMode[rn-1-bothDiatomic]
	    );

	    if(reactionProbabilities[rn] > cloud_.rndGen().sample01<scalar>())
	    {
	      //*****generate probability List for 'all' vibrational mode****************
	      DynamicList<scalar> excitePList = pReactionPDiffVibMode[rn-1-bothDiatomic];

	      exchange(p, q, excitePList, rn-1-bothDiatomic, translationalEnergy);
	    }	    
	  }
	  else if (exchangeQK::posAtomReactant_ == 0) //exchange exsis atom D + ABC in chemicalDict
	  {	      
	    exchangeQK::testExchange
	    (
	     q,
	     translationalEnergy,
	     omegaPQ,
	     rn-1-bothDiatomic,
	     reactionProbabilities[rn],
	     qReactionPDiffVibMode[rn-1-bothDiatomic]
	    );

	    if(reactionProbabilities[rn] > cloud_.rndGen().sample01<scalar>())
	    {
	      //*****generate probability List for 'all' vibrational mode*************
	      DynamicList<scalar> excitePList = qReactionPDiffVibMode[rn-1-bothDiatomic];
	      
	      exchange(q, p, excitePList, rn-1-bothDiatomic, translationalEnergy);
	    }	    
	  }
	  else // both are moleculer
	  {
	    scalar probabilitiesP = 0.0;
	    exchangeQK::testExchange
	    (
	     p,
	     translationalEnergy,
	     omegaPQ,
	     rn-1-bothDiatomic,
	     probabilitiesP,
	     pReactionPDiffVibMode[rn-1-bothDiatomic]
	    );

	    scalar probabilitiesQ = 0.0;
	    exchangeQK::testExchange
	    (
	     q,
	     translationalEnergy,
	     omegaPQ,
	     rn-1-bothDiatomic,
	     probabilitiesQ,
	     qReactionPDiffVibMode[rn-1-bothDiatomic]
	    );
	      
	    reactionProbabilities[rn] = probabilitiesP + probabilitiesQ;

	    if(reactionProbabilities[rn] > cloud_.rndGen().sample01<scalar>())
	    {
	      //*****generate probability List for 'all' vibrational mode****************
	      DynamicList<scalar> excitePList = pReactionPDiffVibMode[rn-1-bothDiatomic];
	      forAll(qReactionPDiffVibMode[rn-1-bothDiatomic], m)
	      {
		excitePList.append(qReactionPDiffVibMode[rn-1-bothDiatomic][m]);
	      }
	      
	      exchange(p, q, excitePList, rn-1-bothDiatomic, translationalEnergy);
	    }	      	     
	  }
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
                            //- Dissociation of P is to occur
                            dissociationQK::dissociateParticleByPartner
                            (
                                p, q, i, vibModeDissoP, translationalEnergy
                            );
                            //- There can't be another reaction: break
                            break;
                        }
                        else if (i == 1)
                        {
                            //- Dissociation of Q is to occur
                            dissociationQK::dissociateParticleByPartner
                            (
                                q, p, i, vibModeDissoQ, translationalEnergy
                            );
                            //- There can't be another reaction: break
                            break;
                        }
			else
                        {
                            //- Exchange reaction
                            if (exchangeQK::posAtomReactant_ != 0)
                            {
                                exchangeQK::exchange
                                (
				 p, q, translationalEnergy, (i-2)
                                );
                            }
                            else
                            {
                                exchangeQK::exchange
                                (
				 q, p, translationalEnergy, (i-2)
                                );
                            }
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
        //  If P is the second reactant, then switch arguments in this
        //  function and P will be first
        dissociationExchangeQK::reaction(q, p);
    }
}


inline label dissociationExchangeQK::nReactionsPerTimeStep() const
{
    return dissociationQK::nReactionsPerTimeStep() 
        + exchangeQK::nReactionsPerTimeStep();
}


void dissociationExchangeQK::outputResults(const label& counterIndex)
{
  if (writeRatesToTerminal_)
  {
        dissociationQK::outputResults(counterIndex);
        exchangeQK::outputResults(counterIndex);
  }
}

}
// End namespace Foam

// ************************************************************************* //
