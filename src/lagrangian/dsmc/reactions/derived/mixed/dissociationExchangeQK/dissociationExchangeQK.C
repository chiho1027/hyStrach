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
		
	scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(2+numberOfExchange, 0.0); 

	const label rnDiss = cloud_.randomLabel(0, bothDiatomic);//2dissociation

	//dissociation probablility
	if( rnDiss == 0)
	{
	  dissociationQK::testDissociation
	  (
	   p,
	   translationalEnergy,
	   vibModeDissoP,
	   reactionProbabilities[0]
	  );
	  
	  totalReactionProbability += reactionProbabilities[0];	
	}
	else //if(rnDiss == bothDiatomic )//&& bothDiatomic == 1)
	{
	  dissociationQK::testDissociation
	  (
	   q,
	   translationalEnergy,
	   vibModeDissoQ,
	   reactionProbabilities[1]
	  );	        	    

	  totalReactionProbability += reactionProbabilities[1];
	}

	/*
	//exchange probablility
	// random choose one reaction	
	const label rnEx = cloud_.randomLabel(0, numberOfExchange-1);

	//method 1:random choose one vibrational mode probability
	const label numberOfVibMode = p.vibLevel().size() + q.vibLevel().size();
	const label randomSelectMode = cloud_.randomLabel(0, numberOfVibMode-1);

	label selectPMode = 0;
	label selectQMode = 0;
	bool isPreact = false;
	
	if(p.vibLevel().size() > 0 && randomSelectMode < p.vibLevel().size())
	{
	  isPreact = true;	  
	  selectPMode = randomSelectMode;
	  
	  exchangeQK::testExchange
	  (
	   p,
	   translationalEnergy,
	   omegaPQ,
	   rnEx,
	   selectPMode,
	   reactionProbabilities[2]
	  );

	  totalReactionProbability += reactionProbabilities[2];
	}
	else
	{	  	  
	  selectQMode = randomSelectMode-p.vibLevel().size();
	  
	  exchangeQK::testExchange
	  (
	   q,
	   translationalEnergy,
	   omegaPQ,
	   rnEx,
	   selectQMode,
	   reactionProbabilities[2]
	  );

	  totalReactionProbability += reactionProbabilities[2];
	}	
	*/
	
	//exchange method 2 	
	const label numberOfVibMode = p.vibLevel().size() + q.vibLevel().size();
	const label randomSelectMode = cloud_.randomLabel(0, numberOfVibMode-1);

	label selectPMode = 0;
	label selectQMode = 0;
	bool isPreact = false;
	      
	if(p.vibLevel().size() > 0 && randomSelectMode < p.vibLevel().size())
	{
	  isPreact = true;
	  selectPMode = randomSelectMode;
	  
	  for(label i=0; i<numberOfExchange; i++)
	  {	      	    	    
	    exchangeQK::testExchange
	    (
	     p,
	     translationalEnergy,
	     omegaPQ,
	     i,
	     selectPMode,
	     reactionProbabilities[i+2]
	    );
	      
	    totalReactionProbability += reactionProbabilities[i+2];
	  }
	}
	else
	{
	  selectQMode = randomSelectMode-p.vibLevel().size();
	  
	  for(label i=0; i<numberOfExchange; i++)
	  {	    	    
	    exchangeQK::testExchange
	    (
	     q,
	     translationalEnergy,
	     omegaPQ,
	     i,
	     selectQMode,
	     reactionProbabilities[i+2]
	    );

	    totalReactionProbability += reactionProbabilities[i+2];
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
			else
                        {			  
                            //- Exchange reaction
			  if( isPreact )
			  {			      
			      exchangeQK::exchange
			      (
			       p, q, selectPMode, i-2, translationalEnergy
			      );
			  }
			  else
                          {			      
			      exchangeQK::exchange
                              (
			       q, p, selectQMode, i-2, translationalEnergy
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
