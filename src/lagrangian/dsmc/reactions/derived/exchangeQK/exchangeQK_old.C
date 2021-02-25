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

#include "exchangeQK.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(exchangeQK, 0);

addToRunTimeSelectionTable(dsmcReaction, exchangeQK, dictionary);


// * * * * * * * * * * *  Protected Member functions * * * * * * * * * * * * //

void exchangeQK::setProperties()
{
    dsmcReaction::setProperties();
    
    if (reactantIds_.size() != 2)
    {
        //- There must be exactly 2 reactants
        FatalErrorIn("exchangeQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "There should be two reactants, instead of " 
            << reactantIds_.size() << nl 
            << exit(FatalError);
    }
    
    exchangeStr_.setSize(numberOfExchange_, word::null);
    nTotExchangeReactions_.setSize(numberOfExchange_, 0);
    nExchangeReactionsPerTimeStep_.setSize(numberOfExchange_, 0);

    bool moleculeFound = false;
    //bool atomFound = false;
    
    forAll(reactantIds_, r)
    {
        //- Check if this reactant is a molecule
        if (reactantTypes_[r] >= 20)
        {
            moleculeFound = true;
            posMolReactant_ = r;
        }
        //- Check if this reactant is an atom
        else if (reactantTypes_[r] == 10 or reactantTypes_[r] == 11)
        {
	  //atomFound = true;
        }
        else
        {
            FatalErrorIn("exchangeQK::setProperties()")
                << "For reaction named " << reactionName_ << nl
                << "Reactant " << cloud_.typeIdList()[reactantIds_[r]]
                << " is neither a molecule nor an atom" << nl 
                << exit(FatalError);
        }
    }
    
    if (!moleculeFound)
    {
        FatalErrorIn("exchangeQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "None of the reactants is a molecule." << nl 
            << exit(FatalError);
    }
    /*
    else if (!atomFound)
    {
        FatalErrorIn("exchangeQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "None of the reactants is an atom." << nl 
            << exit(FatalError);
    }
    */
    
    //- Reading in exchange products
    const List<wordList> productsExchange(propsDict_.lookup("exchangeProducts"));

    if (productsExchange.size() != numberOfExchange_)
    {
        FatalErrorIn("exchangeQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "There should be n products, instead of " 
            << productsExchange.size() << nl 
            << exit(FatalError);
    }

    productIdsExchange_.setSize(productsExchange.size());
    
    forAll(productIdsExchange_, r)
    {
      if (productsExchange[r].size() > 0)
      {     
	//- Check that there are two products
	if (productsExchange[r].size() != 2)
	{
	  FatalErrorIn("exchangeQK::setProperties()")
	    << "For reaction named " << reactionName_ << nl
	    << "There should be 2 exhchange products for molecule " 
	    << cloud_.typeIdList()[reactantIds_[r]] << " instead of " 
	    << productsExchange[r].size() << ", that is "
	    << productsExchange[r]
	    << exit(FatalError);
	}
      }
      else
      {
	//- it should have Exchange products
	  FatalErrorIn("exchagneQK::setProperties()")
	    << "For reaction named " << reactionName_ << nl
	    << "it should have Exchange products "
	    << "instead of " << productsExchange[r]
	    << exit(FatalError);
      }
      
      productIdsExchange_[r].setSize(productsExchange[r].size());
      
      //moleculeFound = false;
      //atomFound = false;
      
      forAll(productIdsExchange_[r], p)
      {
	const label productIndex =  
	  findIndex
	  (
	   cloud_.typeIdList(), 
	   productsExchange[r][p]
	  );
	
        //- Check that products belong to the typeIdList as defined in 
        //  constant/dsmcProperties
        if (productIndex == -1)
	{
            FatalErrorIn("exchangeQK::setProperties()")
                << "For reaction named " << reactionName_ << nl
                << "Cannot find type id: " << productsExchange[r][p] << nl 
                << exit(FatalError);
        }

        //- Check if this product is a molecule
        if (cloud_.constProps(productIndex).type() >= 20)
        {
            moleculeFound = true;
            //- The molecule is set to be the first product
            productIdsExchange_[r][0] = productIndex;
        }
        //- Check if this product is an atom
        else if
        (
            cloud_.constProps(productIndex).type() == 10
         or cloud_.constProps(productIndex).type() == 11
        )
        {
	  //atomFound = true;
            //- The atom is set to be the second product
            productIdsExchange_[r][1] = productIndex;
        }
        else
        {
            FatalErrorIn("exchangeQK::setProperties()")
                << "For reaction named " << reactionName_ << nl
                << "Product " << cloud_.typeIdList()[productIndex]
                << " is neither a molecule nor an atom" << nl 
                << exit(FatalError);
        }
    
	if (!moleculeFound)
	{
	  FatalErrorIn("exchangeQK::setProperties()")
	    << "For reaction named " << reactionName_ << nl
	    << "None of the products is a molecule." << nl 
	    << exit(FatalError);
	}
	/*
	else if (!atomFound)
	{
	  FatalErrorIn("exchangeQK::setProperties()")
	    << "For reaction named " << reactionName_ << nl
	    << "None of the products is an atom." << nl 
	    << exit(FatalError);
	}
	*/

      }// for interal forAll
    }// end forAll
}


void exchangeQK::testExchange
(
    const dsmcParcel& p,
    const scalar translationalEnergy,
    const scalar omegaPQ,
    const label nExIndex,
    scalar& collisionEnergy,
    scalar& totalReactionProbability,
    scalar& reactionProbability
    
)
{
    const label typeIdP = p.typeId();
    
    //- Collision temperature: Eq.(10) of Bird's QK paper.
    const scalar TColl = (translationalEnergy/physicoChemical::k.value())/(2.5 - omegaPQ); 
    
    const scalar aDash =
      //aCoeff_[nExIndex];
        aCoeff_[nExIndex]
       *(
            pow(2.5 - omegaPQ, bCoeff_[nExIndex])
           *exp(lgamma(2.5 - omegaPQ))
           /exp(lgamma(2.5 - omegaPQ + bCoeff_[nExIndex]))
        );
    
    scalar activationEnergy = 
        (
            aDash*pow(TColl/273.0, bCoeff_[nExIndex])
           *fabs(heatOfReactionExchangeJoules_[nExIndex])
        );
    
    if (heatOfReactionExchangeJoules_[nExIndex] < 0.0) 
    {
        //- forward (endothermic) exchange reaction
        activationEnergy -= heatOfReactionExchangeJoules_[nExIndex];
    }
    
    label m = 0;
    do
    {
        const label vibLevel_m = p.vibLevel()[m];
        const scalar kBByThetaVP = physicoChemical::k.value()*cloud_.constProps(typeIdP).thetaV_m(m);
        const scalar EVibP_m = cloud_.constProps(typeIdP).eVib_m(m, vibLevel_m);
        
        //- Total collision energy
        collisionEnergy = translationalEnergy + EVibP_m + p.ERot();
        
        //- Condition for the exchange reaction to possibly occur
        if(collisionEnergy > activationEnergy)
        {
            scalar summation = 0.0;

            if(activationEnergy < kBByThetaVP)
            {
                // this refers to the first sentence in Bird's QK paper after Eq.(12).
                summation = 1.0; 
            }
            else
            {
                const label iaP = collisionEnergy/kBByThetaVP;

                for(label i=0; i<=iaP; i++)
                {
                    summation += 
                        pow
                        (
                            1.0 - cloud_.constProps(typeIdP).eVib_m(m, i)/collisionEnergy,
                            1.5 - omegaPQ
                        );
                }
            }

            //- Based on modified activation energy
            reactionProbability =
                pow
                (
                    1.0 - activationEnergy/collisionEnergy,
                    1.5 - omegaPQ
                )
                /summation;
            
            m = p.vibLevel().size();
        }
        
        m += 1;
        
    } while (m < p.vibLevel().size());
    
    totalReactionProbability += reactionProbability;
}


void exchangeQK::exchange
(
    dsmcParcel& p,
    dsmcParcel& q,
    scalar collisionEnergy,
    const label nExIndex
)
{
    const label typeIdP = p.typeId();
    const label typeIdQ = q.typeId();
    
    nTotExchangeReactions_[nExIndex]++;
    nExchangeReactionsPerTimeStep_[nExIndex]++;
    
    if (allowSplitting_)
    {
        relax_ = false;
        
        vector UP = p.U();
        vector UQ = q.U();
        
        const scalar mP = cloud_.constProps(typeIdP).mass();
        const scalar mQ = cloud_.constProps(typeIdQ).mass();
        const scalar mR = mP*mQ/(mP + mQ);
        const scalar cRsqr = magSqr(UP - UQ);
        
        scalar translationalEnergy = 0.5*mR*cRsqr;
        
        //- Center of mass velocity (pre-exchange)
        const vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);
	
        const label typeIdMol = productIdsExchange_[nExIndex][0];
        const label typeIdAtom = productIdsExchange_[nExIndex][1];

        //- Change species properties
        const scalar mPExch = cloud_.constProps(typeIdAtom).mass();
        const scalar mQExch = cloud_.constProps(typeIdMol).mass();
        const scalar mRExch = mPExch*mQExch/(mPExch + mQExch);
        
        const scalar EVibP = cloud_.constProps(typeIdP).eVib_tot(p.vibLevel());
        //new const scalar EEleP = cloud_.constProps(typeIdP).electronicEnergyList()[p.ELevel()];
        //new const scalar EEleQ = cloud_.constProps(typeIdQ).electronicEnergyList()[q.ELevel()];

	const scalar omegaPQ =
	  0.5
	  *(
	    cloud_.constProps(typeIdP).omega()
	    + cloud_.constProps(typeIdQ).omega()
            );
        
        //  Assumption: no energy redistribution for both particles
        //  All the energy is stored in the translational mode
        translationalEnergy += p.ERot() + EVibP //+ EEleP + EEleQ 
            + heatOfReactionExchangeJoules_[nExIndex];

	//new
	
	// calculate redistribution post particle energy
	labelList vibLevel(cloud_.constProps(typeIdMol).nVibrationalModes(), 0);
	const scalarList& thetaVProduct = cloud_.constProps(typeIdMol).thetaV();
	
	//- Trial L-B redistribution (vibration)  //bu i ding
	forAll(thetaVProduct, m)
	  {
	    label iMaxProduct = (translationalEnergy/(physicoChemical::k.value()*thetaVProduct[m] ));
	    
	    label vibLevelProduct =
	      cloud_.postCollisionVibrationalEnergyLevel
	      (
	       true,
	       0,// vibrationlevel
	       iMaxProduct,
	       cloud_.constProps(typeIdMol).thetaV()[m],
	       cloud_.constProps(typeIdMol).thetaD()[m],
	       cloud_.constProps(typeIdMol).TrefZv()[m],
	       omegaPQ,
	       cloud_.constProps(typeIdMol).Zref()[m],
	       translationalEnergy
	       );
	    
	    vibLevel[m] =  vibLevelProduct;
	    //- Relative translational energy after vibrational energy redistribution
	    translationalEnergy -= vibLevelProduct*cloud_.constProps(typeIdMol).thetaV()[m]*physicoChemical::k.value(); 
	  }
	
	//- Trial L-B redistribution (rotation)
	const scalar energyRatio =
	  cloud_.postCollisionRotationalEnergy
	  (
	   cloud_.constProps(typeIdMol).rotationalDegreesOfFreedom(),
	   2.5 - omegaPQ
	   );
	
	scalar ERotProduct = energyRatio*translationalEnergy;
	
	//- Relative translational energy after rotational energy redistribution
	translationalEnergy -= ERotProduct;

	//Info << " translational Energy is " << translationalEnergy << endl;

	//end new
			
        const scalar relVelExchMol = sqrt(2.0*translationalEnergy/mRExch);
	
        //- Variable Hard Sphere collision part for collision of molecules
        const scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
        const scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
        const scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();
    
        const vector postCollisionRelU =
            relVelExchMol
           *vector
            (
                cosTheta,
                sinTheta*cos(phi),
                sinTheta*sin(phi)
            );
        
        UP = Ucm + postCollisionRelU*mQExch/(mPExch + mQExch);
        UQ = Ucm - postCollisionRelU*mPExch/(mPExch + mQExch);

        //- p is originally the molecule and becomes the atom
        p.typeId() = typeIdAtom;
        p.U() = UP;
        p.ERot() = 0.0; 
        p.vibLevel().setSize
        (
            cloud_.constProps
            (
                typeIdAtom
            ).nVibrationalModes(),
            0
        );
        //p.ELevel() = 0; //new

	// new 
	//- q is originally the atom and becomes the molecule
        q.typeId() = typeIdMol;
        q.U() = UQ;
        q.ERot() = ERotProduct;
        q.vibLevel() = vibLevel;

	/* old
        //- q is originally the atom and becomes the molecule
        q.typeId() = typeIdMol;
        q.U() = UQ;
        q.ERot() = 0.0;
        q.vibLevel().setSize
        (
            cloud_.constProps
            (
                typeIdMol
            ).nVibrationalModes(),
            0
        );
        q.ELevel() = 0;
	*/
	
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
exchangeQK::exchangeQK
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    posMolReactant_(-1),
    numberOfExchange_(readLabel(propsDict_.lookup("numberOfExchange"))),
    productIdsExchange_(),
    exchangeStr_(),
    nTotExchangeReactions_(),
    nExchangeReactionsPerTimeStep_(),
    heatOfReactionExchangeJoules_(),
    aCoeff_(propsDict_.lookup("aCoeff")),
    bCoeff_(propsDict_.lookup("bCoeff")),
    volume_(0.0)  
{
  const scalarList heatOfReactionExchangeJoules(propsDict_.lookup("heatOfReactionExchange"));
  heatOfReactionExchangeJoules_.setSize(numberOfExchange_, 0.0);
  
  forAll(heatOfReactionExchangeJoules_, r)
  {
    heatOfReactionExchangeJoules_[r] = heatOfReactionExchangeJoules[r]*physicoChemical::k.value();
  }
}
  
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

exchangeQK::~exchangeQK()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void exchangeQK::initialConfiguration()
{
    setProperties();
    
    const word& reactantA = cloud_.typeIdList()[reactantIds_[0]];
    const word& reactantB = cloud_.typeIdList()[reactantIds_[1]];

    forAll(productIdsExchange_, r)
    {
      const word& productA = cloud_.typeIdList()[productIdsExchange_[r][0]];
      const word& productB = cloud_.typeIdList()[productIdsExchange_[r][1]];
      
      exchangeStr_[r] = "Exchange reaction " + reactantA + " + " 
        + reactantB + " --> " + productA + " + " + productB;
    }
}


bool exchangeQK::tryReactMolecules
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
        //- Case of dissimilar species colliding (by definition)
        if(reactantPId != reactantQId)
        {
            return true;
        }
    }
        
    return false;
}

void exchangeQK::reaction
(
    dsmcParcel& p,
    dsmcParcel& q,
    const DynamicList<label>& candidateList,
    const List<DynamicList<label> >& candidateSubList,
    const label& candidateP,
    const List<label>& whichSubCell
)
{}

void exchangeQK::reaction(dsmcParcel& p, dsmcParcel& q)
{
    //- Reset the relax switch
    relax_ = true;
    
    const label typeIdP = p.typeId();
    const label typeIdQ = q.typeId();
    
    //- Exchange reaction AB + C --> A + BC 
    //  If P is the first reactant AB (i.e., not the atom)
    //  NB: Q is necessarily M otherwise this class would not have been selected
    if
    (
     typeIdP == reactantIds_[0]
     //  cloud_.constProps(typeIdP).type() != 10
     //&& cloud_.constProps(typeIdP).type() != 11
    ) 
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
        
        //- Possible reactions:
        // 1. Exchange reaction
        
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(numberOfExchange_, 0.0);
        scalarList collisionEnergies(numberOfExchange_, 0.0);

	for(label r=0; r<numberOfExchange_; r++)
	{
	  if (posMolReactant_ == 0)
	  {
	    testExchange
	    (
	     p,
	     translationalEnergy,
	     omegaPQ,
	     r,
	     collisionEnergies[r],
	     totalReactionProbability,
	     reactionProbabilities[r]
	    );
	  }
	  else
	  {
	    testExchange
	    (
	     q,
	     translationalEnergy,
	     omegaPQ,
	     r,
	     collisionEnergies[r],
	     totalReactionProbability,
	     reactionProbabilities[r]
	    );
	  }	  
	}// end for
        
        //- Decide if an exchange reaction is to occur
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
		for(label r=0; r<numberOfExchange_; r++)
		{
		  //- Current reaction is to occur
		  if (i == r)
		  {
		    //- Exchange reaction
		    if (posMolReactant_ == 0)
		    {
		      exchange(p, q, collisionEnergies[i], i);
		    }
		    else
		    {
		      exchange(q, p, collisionEnergies[i], i);
		    }
		    //- There can't be another reaction: break
		    break;
		  }
		}// end forAll
	      }// end cumulate if
	    }// end reaction occur
	    else
	    {
	      //- All the following possible reactions have a probability
	      //  of zero
	      break;
	    }
	  }// end forAll
        }//end decied
    }
    else
    {
        //- If P is the second reactant M, then switch arguments in this
        //  function and P will be first
        exchangeQK::reaction(q, p);
    }
}

void exchangeQK::outputResults(const label& counterIndex)
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
        
        labelList nTotExchangeReactions = nTotExchangeReactions_;
        labelList nExchangeReactionsPerTimeStep = nExchangeReactionsPerTimeStep_;

        const scalar deltaT = mesh_.time().deltaT().value();
        scalar factor = 0.0;
        
        if (reactantIds_[0] == reactantIds_[1] && numberDensities[0] > 0.0)
        {
            factor = cloud_.nParticle()/
                (
                    counterIndex*deltaT
                   *numberDensities[0]*numberDensities[0]
                   *volume
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

	for(label r=0; r<numberOfExchange_; r++)
	{  
	  if (exchangeStr_[r].size())
	  {
	    if (Pstream::parRun())
	      {
		//- Parallel communication
		reduce(molsReactants[0], sumOp<label>());
		reduce(molsReactants[1], sumOp<label>());
		reduce(nTotExchangeReactions[r], sumOp<label>());
		reduce(nExchangeReactionsPerTimeStep[r], sumOp<label>());
	      }
	    
	    const scalar reactionRateExchange = factor*nTotExchangeReactions[r];
	    
	    Info<< exchangeStr_[r]
		<< ", reaction rate = " << reactionRateExchange
		<< ", nReactions = " << nExchangeReactionsPerTimeStep[r]
		<< endl;
	  }
	}// end forAll
    }
    else
    {
        labelList nTotExchangeReactions = nTotExchangeReactions_;   
        labelList nExchangeReactionsPerTimeStep = nExchangeReactionsPerTimeStep_;

	for(label r=0; r<numberOfExchange_; r++)
	{  
	  if (exchangeStr_[r].size())
	  {
	    if (Pstream::parRun())
	    {
	      //- Parallel communication
	      reduce(nTotExchangeReactions[r], sumOp<label>());
	      reduce(nExchangeReactionsPerTimeStep[r], sumOp<label>());
	    }
        
	    if (nTotExchangeReactions[r] > 0)
	    {
	      Info<< exchangeStr_[r]
		  << " is active, nReactions this time step = " 
		  << nExchangeReactionsPerTimeStep[r]
		  << endl;
	    }
	  }
	}// end forAll
	
    }// end else
    
    nExchangeReactionsPerTimeStep_ = 0;
}
  
}
// End namespace Foam

// ************************************************************************* //
