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

#include "dissociationExchangeRecombinationQK.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dissociationExchangeRecombinationQK, 0);

addToRunTimeSelectionTable
(
    dsmcReaction,
    dissociationExchangeRecombinationQK,
    dictionary
);


// * * * * * * * * * * *  Protected Member functions * * * * * * * * * * * * //

void dissociationExchangeRecombinationQK::setProperties()
{
    dissociationQK::setProperties();
    exchangeQK::setProperties();
    recombinationQK::setProperties();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dissociationExchangeRecombinationQK::dissociationExchangeRecombinationQK
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    dissociationQK(t, cloud, dict),
    exchangeQK(t, cloud, dict),
    recombinationQK(t, cloud, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dissociationExchangeRecombinationQK::~dissociationExchangeRecombinationQK()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dissociationExchangeRecombinationQK::initialConfiguration()
{
    dissociationQK::initialConfiguration();
    exchangeQK::initialConfiguration();
    recombinationQK::initialConfiguration();
}


bool dissociationExchangeRecombinationQK::tryReactMolecules
(
    const label& typeIdP,
    const label& typeIdQ
)
{
    return dissociationQK::tryReactMolecules(typeIdP, typeIdQ);
}


void dissociationExchangeRecombinationQK::reaction
( dsmcParcel& p, dsmcParcel& q )
{}

void dissociationExchangeRecombinationQK::reaction
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
        const scalar cRsqr = magSqr(p.U() - q.U());
        const scalar translationalEnergy = 0.5*mR*cRsqr;
	const scalar omegaPQ =
	  0.5
	  *(
	    cloud_.constProps(typeIdP).omega()
	    + cloud_.constProps(typeIdQ).omega()
           );
	
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
	// 3. Recombination
        // 4~?. Exchanges (many) 

        //scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(3+numberOfExchange, 0.0);

	//numberofexchange + 2dissociation + recomb OR numberofexchange + dissociation + recomb
	const label rn = cloud_.randomLabel(0, numberOfExchange + bothDiatomic + 1 );
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
	else if (rn == 1 + bothDiatomic)
	{
	  //find thirdBody Index 
	  const label  nR = findIndex( thirdBodyId_, thirdBody.typeId() );
	  
	  recombinationQK::testRecombination
	  (
	   p,
	   q,
	   translationalEnergy,
	   omegaPQ,
	   nR,
	   reactionProbabilities[2]
	  );

	  if (reactionProbabilities[2] > cloud_.rndGen().sample01<scalar>())
	  {
	    recombinationQK::recombination(p, q, thirdBody, nR, translationalEnergy);
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
	     rn-2-bothDiatomic,
	     reactionProbabilities[rn],
	     pReactionPDiffVibMode[rn-2-bothDiatomic]
	    );

	    if(reactionProbabilities[rn] > cloud_.rndGen().sample01<scalar>())
	    {
	      //*****generate probability List for 'all' vibrational mode****************
	      DynamicList<scalar> excitePList = pReactionPDiffVibMode[rn-2-bothDiatomic];

	      exchange(p, q, excitePList, rn-2-bothDiatomic, translationalEnergy);
	    }	    
	  }
	  else if (exchangeQK::posAtomReactant_ == 0) //exchange exsis atom D + ABC in chemicalDict
	  {	      
	    exchangeQK::testExchange
	    (
	     q,
	     translationalEnergy,
	     omegaPQ,
	     rn-2-bothDiatomic,
	     reactionProbabilities[rn],
	     qReactionPDiffVibMode[rn-2-bothDiatomic]
	    );

	    if(reactionProbabilities[rn] > cloud_.rndGen().sample01<scalar>())
	    {
	      //*****generate probability List for 'all' vibrational mode*************
	      DynamicList<scalar> excitePList = qReactionPDiffVibMode[rn-2-bothDiatomic];
	      
	      exchange(q, p, excitePList, rn-2-bothDiatomic, translationalEnergy);
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
	     rn-2-bothDiatomic,
	     probabilitiesP,
	     pReactionPDiffVibMode[rn-2-bothDiatomic]
	    );

	    scalar probabilitiesQ = 0.0;
	    exchangeQK::testExchange
	    (
	     q,
	     translationalEnergy,
	     omegaPQ,
	     rn-2-bothDiatomic,
	     probabilitiesQ,
	     qReactionPDiffVibMode[rn-2-bothDiatomic]
	    );
	      
	    reactionProbabilities[rn] = probabilitiesP + probabilitiesQ;

	    if(reactionProbabilities[rn] > cloud_.rndGen().sample01<scalar>())
	    {
	      //*****generate probability List for 'all' vibrational mode****************
	      DynamicList<scalar> excitePList = pReactionPDiffVibMode[rn-2-bothDiatomic];
	      forAll(qReactionPDiffVibMode[rn-2-bothDiatomic], m)
	      {
		excitePList.append(qReactionPDiffVibMode[rn-2-bothDiatomic][m]);
	      }
	      
	      exchange(p, q, excitePList, rn-2-bothDiatomic, translationalEnergy);
	    }	      	     
	  }
	}	

	
	/*
        dissociationQK::testDissociation
        (
            p,
            translationalEnergy,
            vibModeDissoP,
            reactionProbabilities[0]
        );
	totalReactionProbability += reactionProbabilities[0];
	
        dissociationQK::testDissociation
        (
            q,
            translationalEnergy,
            vibModeDissoQ,
            reactionProbabilities[1]
        );
	totalReactionProbability += reactionProbabilities[1];

	recombinationQK::testRecombination
        (
            p,
	    q,
            translationalEnergy,
	    omegaPQ,
	    nR,
            reactionProbabilities[2]
        );
	totalReactionProbability += reactionProbabilities[2];

	for(label i=0; i<numberOfExchange; i++)
	{
	    if (exchangeQK::posAtomReactant_ == 1)
	    {
	      reactionProbabilities[3+i] = 
		exchangeQK::testExchange
		(
		 p,
		 translationalEnergy,
		 omegaPQ,
		 i
		 );
	      totalReactionProbability += reactionProbabilities[3+i];
	    }
	    else if (exchangeQK::posAtomReactant_ == 0)
	    {
	      reactionProbabilities[3+i] = 
		exchangeQK::testExchange
		(
		 q,
		 translationalEnergy,
		 omegaPQ,
		 i
		 );
	      totalReactionProbability += reactionProbabilities[3+i];
	    }
	    else
	    {
	      const scalar probabilitiesP =
		exchangeQK::testExchange
		(
		 p,
		 translationalEnergy,
		 omegaPQ,
		 i
		 );

	      const scalar probabilitiesQ =
		exchangeQK::testExchange
		(
		 q,
		 translationalEnergy,
		 omegaPQ,
		 i
		 );
	      
	      reactionProbabilities[3+i] = probabilitiesP + probabilitiesQ;
	      totalReactionProbability += reactionProbabilities[3+i];
	    }
	}
	
        
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
			else if (i == 2)
                        {
                            //- Recombination of P is to occur
			    recombinationQK::recombination
                            (
			        p, q, nR, translationalEnergy 
                            );
                            //- There can't be another reaction: break
                            break;
                        }
                        else // i>2
                        {
                            //- Exchange reaction
                            if (exchangeQK::posAtomReactant_ != 0)
                            {
                                exchangeQK::exchange
                                (
				 p, q, translationalEnergy, (i-3)
                                );
                            }
                            else
                            {
                                exchangeQK::exchange
                                (
				 q, p, translationalEnergy, (i-3)
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
      dissociationExchangeRecombinationQK::reaction(q, p, thirdBody);
    }
}


inline label dissociationExchangeRecombinationQK::nReactionsPerTimeStep() const
{
  return dissociationQK::nReactionsPerTimeStep() 
    + exchangeQK::nReactionsPerTimeStep()
    + recombinationQK::nReactionsPerTimeStep();
}


void dissociationExchangeRecombinationQK::outputResults(const label& counterIndex)
{
  if (writeRatesToTerminal_)
  {
        dissociationQK::outputResults(counterIndex);
        exchangeQK::outputResults(counterIndex);
	recombinationQK::outputResults(counterIndex);
  }
}

}
// End namespace Foam

// ************************************************************************* //
