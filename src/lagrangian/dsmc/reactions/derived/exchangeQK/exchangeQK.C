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
	  // posMolReactant_ = r;
        }
        //- Check if this reactant is an atom
        else if (reactantTypes_[r] == 10 or reactantTypes_[r] == 11)
        {
	  // atomFound = true;
	  posAtomReactant_ = r;
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
      
      moleculeFound = false;
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
            productIdsExchange_[r][p] = productIndex;
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
            productIdsExchange_[r][p] = productIndex;
        }
        else
        {
            FatalErrorIn("exchangeQK::setProperties()")
                << "For reaction named " << reactionName_ << nl
                << "Product " << cloud_.typeIdList()[productIndex]
                << " is neither a molecule nor an atom" << nl 
                << exit(FatalError);
        }

      }// for interal forAll

      if (!moleculeFound)
      {
	FatalErrorIn("exchangeQK::setProperties()")
	  << "For reaction named " << reactionName_ << nl
	  << "None of the products is a molecule." << nl 
	  << exit(FatalError);
      }

      
    }// end forAll
}

scalar exchangeQK::equilibriumDistribution
(
 const dsmcParcel& p,
 const scalar translationalEnergy,
 const scalar omegaPQ,
 const label nExIndex 
 )
{
  const label typeIdP = p.typeId();

    scalar TMacro = cloud_.fields().overallT(p.cell());
    
    //- Collision temperature: Eq.(10) of Bird's QK paper.
    /*
    const scalar TColl = (translationalEnergy/physicoChemical::k.value())/(2.5 - omegaPQ);    
    const scalar aDash = 
        aCoeff_[nExIndex][0]
       *(
            pow(2.5 - omegaPQ, bCoeff_[nExIndex][0])
           *exp(lgamma(2.5 - omegaPQ))
           /exp(lgamma(2.5 - omegaPQ + bCoeff_[nExIndex][0]))
        );
    */

    scalar activationEnergy = 
      (
       aCoeff_[nExIndex][0]*pow(TMacro/273, bCoeff_[nExIndex][0])// changed aDash
       *fabs(heatOfReactionExchangeJoules_[nExIndex])
      );    
    
    if (heatOfReactionExchangeJoules_[nExIndex] < 0.0) 
    {
        //- forward (endothermic) exchange reaction
        activationEnergy -= heatOfReactionExchangeJoules_[nExIndex];
    }

    const label vibLevel_m = p.vibLevel()[0];
    const scalar kBByThetaVP = physicoChemical::k.value()*cloud_.constProps(typeIdP).thetaV_m(0);
    const scalar EVibP_m = cloud_.constProps(typeIdP).eVib_m(0, vibLevel_m);
    const scalar collisionEnergy = translationalEnergy + EVibP_m;
    const label iMax = collisionEnergy/kBByThetaVP;

    label  j    = 0;
    scalar func = 0.0;
    do
    {
      func = 0.0;
      j = cloud_.randomLabel(0, iMax);

      if(j*kBByThetaVP < activationEnergy)
      {
	  func =
	    exp((-j*kBByThetaVP)/(physicoChemical::k.value()*TMacro));	
      }
      else
      {
	  func =
	  pow((collisionEnergy - j*kBByThetaVP)/(collisionEnergy-activationEnergy) ,1.5 - omegaPQ)
	  *exp(-activationEnergy/(physicoChemical::k.value()*TMacro));	
      }            
    }while(func < cloud_.rndGen().sample01<scalar>());
    
    //return modifiyTranslationalEnergy
    return collisionEnergy - j*kBByThetaVP;
}

void exchangeQK::postReactionVibrationalRedistribution
(
 const label excite,
 const scalarList thetaVProduct,
 const scalar TMacro,
 const scalar activationEnergy,
 const scalar reverseOmega,
 labelList& vibLevel,
 scalar& Ec
)
{
  if( thetaVProduct.size() > 0 )
  {
    const label m = excite;
    
      // const 
      const scalar thetaPrim   = thetaVProduct[m]/TMacro; 
      const scalar kBByThetaVP = physicoChemical::k.value()*thetaVProduct[m];
      
      const scalar w           = 1.5 - reverseOmega;
      const scalar EaPrim      = activationEnergy/kBByThetaVP;
      const scalar lnTwo       = log(2.0)/thetaPrim;
      const scalar three       = 3.0/thetaPrim;
      const label  iMax        = Ec/kBByThetaVP;
      
      //scalar kmax = 0.0;
      //scalar jmax = 0.0;
      //kmaxFunc(lnTwo, three, EaPrim , w, thetaPrim, kmax, jmax);

      //jmax at j = 0
      scalar kmax = (4.0*ki(lnTwo, 0, false, EaPrim, w) + ki(three, 0, false, EaPrim, w) );
      
      scalar func = 0.0;
      label j = 0;
      
      do // acceptance - rejection
      {
	j = pow(iMax, 1.3)*cloud_.rndGen().sample01<scalar>();
	//j = 25.0*cloud_.rndGen().sample01<scalar>();
	//j = iMax*cloud_.rndGen().sample01<scalar>();
	
	if(j > EaPrim )
        {
	  func = (4.0*ki(lnTwo, j, true, EaPrim, w) + ki(three, j, true, EaPrim, w) )// *thetaVProduct.size()
	    *exp(-j*thetaPrim + EaPrim*thetaPrim)
	    /kmax;
	}
	else
	{
	  func = (4.0*ki(lnTwo, j, false, EaPrim, w) + ki(three, j, false, EaPrim, w) )// *thetaVProduct.size()
	    /kmax;	  
	}	
      } while(func < cloud_.rndGen().sample01<scalar>());
      
      vibLevel[m] = j;
      Ec -= j*kBByThetaVP;
      
      if(Ec < 0.0)
      {	
	Ec          += j*kBByThetaVP;
	vibLevel[m]  = int(Ec/kBByThetaVP);
	Ec          -= vibLevel[m]*kBByThetaVP;
      }

  }//end if
}



scalar exchangeQK::ki
(
 scalar chai,
 scalar j,
 const bool   largeThanEa,
 const scalar Ea,
 const scalar w
)
{
  scalar ki = 0.0;
  
  if( largeThanEa )
  {
    ki = pow(chai-Ea/(1.0+j/chai), w)/(chai+j+(w+1.0)/2.0);   
    return ki;
  }
  else
  {
    ki = pow(chai-j/(1.0+Ea/chai), w)/(chai+Ea+(w+1.0)/2.0);
    return ki;
  }
}

void exchangeQK::kmaxFunc
(
 const scalar lnTwo,
 const scalar three,
 const scalar EaPrim,
 const scalar w,
 const scalar thetaPrim,
 scalar& kmax,
 scalar& jmax
)
{ 
  scalar jo   = 0.0;
  scalar j1   = EaPrim + 10.0;//upper j
  jmax = (jo+j1)/2.0;
  scalar tol  = 0.001;//tolerence
  
  while( fabs( 4.0*jmaxP(lnTwo, jmax, EaPrim, w, thetaPrim)
	          +jmaxP(three, jmax, EaPrim, w, thetaPrim) ) >= tol )
  {
    if(  ( 4.0*jmaxP(lnTwo, jmax, EaPrim, w, thetaPrim) + jmaxP(three, jmax, EaPrim, w, thetaPrim) )
	*( 4.0*jmaxP(lnTwo, jo,   EaPrim, w, thetaPrim) + jmaxP(three, jo,   EaPrim, w, thetaPrim) ) < 0.0 )
    {
      j1 = jmax;
    }
    else
    {
      jo = jmax; 
    }    
    jmax  = (jo+j1)/2.0;
    
    //if(jmax == EaPrim)
    //  break;
  }

  //scalar kmax = 0.0;
  if( jmax > EaPrim )
  {
    kmax = (4.0*ki(lnTwo, jmax, true, EaPrim, w) + ki(three, jmax, true, EaPrim, w) );    	        
  }
  else
  {
    kmax = (4.0*ki(lnTwo, jmax, false, EaPrim, w) + ki(three, jmax, false, EaPrim, w) )
      *exp(jmax*thetaPrim - EaPrim*thetaPrim);
  }      
  //return kmax;
}

scalar exchangeQK::jmaxP
(
 const scalar chai,
 const scalar j,
 const scalar EaPrim,
 const scalar w,
 const scalar thetaPrim
)
{
  scalar jmaxP = 0.0;
  
  if(j > EaPrim)
  {
    jmaxP =
      pow(chai - EaPrim/(1.0 + j/chai), w)
     /pow(chai + j + (w + 1.0)/2.0, 2.0)
     *( EaPrim*w*( 1.0 + (w + 1)/(2.0*(j + chai)) )/(chai + j - EaPrim) - 1.0 );
  }
  else
  {
    jmaxP =
      pow(chai - j/(1.0 + EaPrim/chai), w)
     /pow(chai + j + (w + 1.0)/2.0, 1.0)
     *( thetaPrim - w/(EaPrim + chai - j));
  }
  
  return jmaxP;
}

  /*
void exchangeQK::testExchange
(
    const dsmcParcel& p,
    const scalar translationalEnergy,
    const scalar omegaPQ,
    const label nExIndex,
    scalar& reactionProbability,
    DynamicList<scalar>& reactPDiffVibMode
)
{
    const label typeIdP = p.typeId();

    scalar TMacro = cloud_.fields().overallT(p.cell());
 
    
    //if(TMacro == 0.0)
    //{
    //  Pout << "T = 0 !!!!" << endl;
    //  TMacro = cloud_.fields().calculateTMacro(p.cell());
    //}
    
    
    //- Collision temperature: Eq.(10) of Bird's QK paper.
    
    //const scalar TColl = (translationalEnergy/physicoChemical::k.value())/(2.5 - omegaPQ);    
    //const scalar aDash = 
    //    aCoeff_[nExIndex][0]
    //   *(
    //        pow(2.5 - omegaPQ, bCoeff_[nExIndex][0])
    //       *exp(lgamma(2.5 - omegaPQ))
    //       /exp(lgamma(2.5 - omegaPQ + bCoeff_[nExIndex][0]))
    //   );    
    
    scalar activationEnergy = 
      (
       aCoeff_[nExIndex][0]*pow(TMacro/273.0, bCoeff_[nExIndex][0])// changed aDash
       *fabs(heatOfReactionExchangeJoules_[nExIndex])
      );    
    
    if (heatOfReactionExchangeJoules_[nExIndex] < 0.0) 
    {
        //- forward (endothermic) exchange reaction
        activationEnergy -= heatOfReactionExchangeJoules_[nExIndex];
    }
    
    // reset reactionProbability
    reactionProbability = 0.0;
 
    forAll(p.vibLevel(), m)//m=mode
    {
      scalar singleVibModeP = 0.0;
      
      const label vibLevel_m = p.vibLevel()[m];
      const scalar kBByThetaVP = physicoChemical::k.value()*cloud_.constProps(typeIdP).thetaV_m(m);
      const scalar EVibP_m = cloud_.constProps(typeIdP).eVib_m(m, vibLevel_m);
      //- Total collision energy// QK paper
      const scalar collisionEnergy = translationalEnergy + EVibP_m;
      
      //- Condition for the exchange reaction to possibly occur
      if(collisionEnergy > activationEnergy)
      {
	scalar summation = 0.0;
	
	//if(activationEnergy < kBByThetaVP)
	//{
	  // this refers to the first sentence in Bird's QK paper after Eq.(12).
	//  summation = 1.0;	  
	//}
	//else
	//{
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
	//}
	
	//- Based on modified activation energy
	singleVibModeP =
	  pow
	  (
	   1.0 - activationEnergy/collisionEnergy,
	   1.5 - omegaPQ
	  )
	  /summation;
	
	reactionProbability += singleVibModeP;
      
      }// end if

      reactPDiffVibMode.append(singleVibModeP);
      //reactPDiffVibMode[m] = singleVibModeP;
      
    }//end for
    
    //return reactionProbability;
}
  */

void exchangeQK::testExchange
(
    const dsmcParcel& p,
    const scalar translationalEnergy,
    const scalar omegaPQ,
    const label nExIndex,
    const label selectMode,
    scalar& reactionProbability
)
{  
    const label typeIdP = p.typeId();

    scalar activationEnergy = 0.0;
    if(bCoeff_[nExIndex][0] < -1.5)
    {
      scalar TMacro = cloud_.fields().overallT(p.cell());

      activationEnergy = 
      (
       aCoeff_[nExIndex][0]*pow(TMacro/273.0, bCoeff_[nExIndex][0])
       *fabs(heatOfReactionExchangeJoules_[nExIndex])
      );  
    }
    else
    {
      //- Collision temperature: Eq.(10) of Bird's QK paper.    
      const scalar TColl = (translationalEnergy/physicoChemical::k.value())/(2.5 - omegaPQ);    
      const scalar aDash = 
        aCoeff_[nExIndex][0]
	*(
	  pow(2.5 - omegaPQ, bCoeff_[nExIndex][0])
	  *exp(lgamma(2.5 - omegaPQ))
	  /exp(lgamma(2.5 - omegaPQ + bCoeff_[nExIndex][0]))
	 );    

      activationEnergy = 
      (
       aDash*pow(TColl/273.0, bCoeff_[nExIndex][0])// changee d aDash
       *fabs(heatOfReactionExchangeJoules_[nExIndex])
      );        
    }       
    
    //modify with probability in endothermic reaction(heat< 0)
    scalar probabilityModifyFactor = 1.0;      
    
    if (heatOfReactionExchangeJoules_[nExIndex] < 0.0) 
    {
        //- forward (endothermic) exchange reaction
        activationEnergy -= heatOfReactionExchangeJoules_[nExIndex];	
    }
    else
    {      
       probabilityModifyFactor = aCoeff_[nExIndex][1];
    }
 
    // reset reactionProbability
    reactionProbability = 0.0;
    
    const label m = selectMode;
    const label vibLevel_m = p.vibLevel()[m];
    const scalar kBByThetaVP = physicoChemical::k.value()*cloud_.constProps(typeIdP).thetaV_m(m);
    const scalar EVibP_m = cloud_.constProps(typeIdP).eVib_m(m, vibLevel_m);
    //- Total collision energy// QK paper
    const scalar collisionEnergy = translationalEnergy + EVibP_m;
	
	
    //- Condition for the exchange reaction to possibly occur
    if(collisionEnergy > activationEnergy)
    {
      scalar summation = 0.0;

      //if(activationEnergy < kBByThetaVP)
      //{
	// this refers to the first sentence in Bird's QK paper after Eq.(12).
      //	summation = 1.0;	  
      //}
      //else
      //{
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
      //}
	
      //- Based on modified activation energy
      reactionProbability =		
	probabilityModifyFactor*
	pow
	(
	 1.0 - activationEnergy/collisionEnergy,
	 1.5 - omegaPQ
	)
	/summation;
      
    }// end if
}

void exchangeQK::exchange
(
    dsmcParcel& p,
    dsmcParcel& q,
    const label selectMode,
    const label nExIndex,
    const scalar translationalEnergy
)
{
    const label typeIdP = p.typeId();
    const label typeIdQ = q.typeId();
    
    nTotExchangeReactions_[nExIndex]++;
    nExchangeReactionsPerTimeStep_[nExIndex]++;

    //temp
    //relax_ = true;
    //return;      
    
    if (allowSplitting_)
    {
      //if reaction happen relax is false
      relax_ = false;

      //relax_ = true;//temp
              
        const scalar& heat =
	  heatOfReactionExchangeJoules_[nExIndex];	

	scalar reactedEnergy =
	  translationalEnergy
	  + cloud_.constProps(typeIdP).eVib_m(selectMode, p.vibLevel()[selectMode])
	  + heat;      

        vector UP = p.U();
        vector UQ = q.U();
	
        const scalar mP = cloud_.constProps(typeIdP).mass();
        const scalar mQ = cloud_.constProps(typeIdQ).mass();
  
	
        //- Center of mass velocity (pre-exchange)
        const vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

        const label typeIdMole = productIdsExchange_[nExIndex][0];   //molecule
        const label typeIdSecond = productIdsExchange_[nExIndex][1]; //atom or molecule

        //- Change species properties
        const scalar mPExch = cloud_.constProps(typeIdSecond).mass();
        const scalar mQExch = cloud_.constProps(typeIdMole).mass();
        const scalar mRExch = mPExch*mQExch/(mPExch + mQExch);

	const scalarList& thetaVProductMole = cloud_.constProps(typeIdMole).thetaV();
	const scalarList& thetaVProductSecond = cloud_.constProps(typeIdSecond).thetaV();

        const scalar reverseOmega =  0.5*(
				    cloud_.constProps(typeIdMole).omega()
				    + cloud_.constProps(typeIdSecond).omega()
				    );
	
	//- Dof
	const scalar& moleRDof   = cloud_.constProps(typeIdMole).rotationalDegreesOfFreedom();
	const scalar& secondRDof = cloud_.constProps(typeIdSecond).rotationalDegreesOfFreedom();

	////////////////////////////////////// start ///////////////////////////////////////////////
	const scalar TMacro = cloud_.fields().overallT(p.cell());
	
	labelList vibLevelMole(thetaVProductMole.size(), 0);
	labelList vibLevelSecond(thetaVProductSecond.size(), 0);

	const labelList vibBolzMole =
	  cloud_.equipartitionVibrationalEnergyLevel(TMacro, vibLevelMole.size(), typeIdMole);
	const labelList vibBolzSecond =
	   cloud_.equipartitionVibrationalEnergyLevel(TMacro, vibLevelSecond.size(), typeIdSecond);
	
	// determin which vibration mode will participate reaction for product	
	DynamicList<scalar> productExciteP(0);      
	forAll(thetaVProductMole, u)
	{	  
	  scalar sum = 0.0;
	  const scalar equiTranE = cloud_.equipartitionRotationalEnergy(TMacro, 5.0-2.0*reverseOmega);
	  const scalar kThetaV = thetaVProductMole[u]*physicoChemical::k.value();
	  const scalar Ec      = equiTranE + vibBolzMole[u]*kThetaV;
	  
	  const label iMax     = Ec/kThetaV;
	
	  for(label i=0; i<=iMax; i++)
	  {
	    sum += 
	      pow
	      (
	       1.0 - (i*kThetaV)/Ec,
	       1.5 - reverseOmega
	      );
	  }

	  productExciteP.append(1.0/sum);
	  
	  /*
	  productExciteP.append( 1./(
				     0.5+ Ec/
				     (
				      (2.5-reverseOmega)*thetaVProductMole[u]*physicoChemical::k.value()
				     )
				    )
			       );
	  */	
	}
	
	forAll(thetaVProductSecond, u)
	{	  
	  scalar sum = 0.0;
	  const scalar equiTranE = cloud_.equipartitionRotationalEnergy(TMacro, 5.0-2.0*reverseOmega);
	  const scalar kThetaV = thetaVProductSecond[u]*physicoChemical::k.value();
	  const scalar Ec      = equiTranE + vibBolzSecond[u]*kThetaV;
	  
	  const label iMax     = Ec/kThetaV;
	  
	  for(label i=0; i<=iMax; i++)
	  {
	    sum += 
	      pow
	      (
	       1.0 - (i*kThetaV)/Ec,
	       1.5 - reverseOmega
	      );
	  }

	  productExciteP.append(1.0/sum);
	  
	  /*
	  productExciteP.append( 1./(
				     0.5+ Ec/
				     (
				      (2.5-reverseOmega)*thetaVProductSecond[u]*physicoChemical::k.value()
				     )
				    )
			       );
	  */   
	}

	const label postExciteMode = selectExciteMode(productExciteP);

	label      excite   =  0;
	scalarList theta    =  thetaVProductMole;//assume
	labelList* vibLevel = &vibLevelMole;//assume       
	
	if(postExciteMode < thetaVProductMole.size() )
	{
	  excite   =  postExciteMode;
	}
	else
	{
	  excite   =  postExciteMode-thetaVProductMole.size();
	  theta    =  thetaVProductSecond;	  
	  vibLevel =  &vibLevelSecond;
	}
	
	const scalar kBByThetaVP = physicoChemical::k.value()*theta[excite];
	const label iMax = reactedEnergy/kBByThetaVP;
	
	if(iMax == 0)
	{
	  (*vibLevel)[excite] = 0;
	}
	else
	{
	  scalar func  = 0.0;
	  label  j     =   0;
    
	  //select postProuct pre-reaction vibrational level
	  do // acceptance - rejection
	  {
	    j = cloud_.randomLabel(0, iMax);
	    
	    func = 
	      pow(1.0 - j*kBByThetaVP/reactedEnergy ,1.5 - reverseOmega);
            
	  }while(func < cloud_.rndGen().sample01<scalar>() );

	  (*vibLevel)[excite] = j;	  
	  reactedEnergy    -= j*kBByThetaVP;	  	    
	}

	//claculate post equilibrium energy
	
	// determine equilibrium energy
	const scalar& pRDof = cloud_.constProps(typeIdP).rotationalDegreesOfFreedom();
	const scalar& qRDof = cloud_.constProps(typeIdQ).rotationalDegreesOfFreedom();
	
	const scalar ERotP = p.ERot();
	const scalar ERotQ = q.ERot();

	scalar equilibriumE    = ERotP + ERotQ;
	scalar equilibriumRDof = pRDof + qRDof;

	
	const label index = productExciteP.size() - 2;
	
	scalar equilibriumVE   = 0.0;
	scalar equilibriumVDof = 0.0;
	
	//add equilibrium vibrational energy
	forAll(p.vibLevel(), m)
	{
	  if(m != selectMode)
	  {
	    const scalar& equiVibEP =
	      cloud_.constProps(typeIdP).eVib_m(m, p.vibLevel()[m]);	    
	    const scalarList equiVibTheta(1,cloud_.constProps(typeIdP).thetaV()[m]);
	    
	    equilibriumVDof += 2.0*cloud_.temperatureFunction(TMacro, equiVibTheta);

	    if(index == 0)
	    {
	      equilibriumVE += equiVibEP;
	    }
	    else
	    {
	      equilibriumE  += equiVibEP;
	    }
	  }
	}

       	const scalar& equiVibEQ =
	  cloud_.constProps(q.typeId()).eVib_tot(q.vibLevel());
	const scalarList& thetaVQ = cloud_.constProps(typeIdQ).thetaV();
	  
	equilibriumVDof  += 2.0*cloud_.temperatureFunction(TMacro, thetaVQ);

	if(index == 0)
	{
	  equilibriumVE += equiVibEQ;
	}
	else
	{
	  equilibriumE  += equiVibEQ;
	}        	
	
	forAll(productExciteP, u)
	{
	  label mode = u;
	  
	  if(mode != postExciteMode) 
	  {	    
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

	    if(index == 0)
	    {
	      equilibriumRDof -= 2.0;
	      const scalar energyRatio = cloud_.postCollisionRotationalEnergy( 2.0, equilibriumRDof/2.0);
	      const scalar kBByThetaV  = physicoChemical::k.value()*theta[mode];
	      (*vibLevel)[mode]        = int(energyRatio*equilibriumE/kBByThetaV);	      
	      equilibriumE            -= (*vibLevel)[mode]*kBByThetaV;
	      
	      const scalarList thetaVMode(1, theta[mode]);
	      equilibriumRDof += 2.0 - 2.0*cloud_.temperatureFunction(TMacro, thetaVMode);
	    }
	    else
	    {	      
	      const scalarList thetaVMode(1, theta[mode]);
	      equilibriumVDof  -= 2.0*cloud_.temperatureFunction(TMacro, thetaVMode);	     
	      const scalar remainDof = equilibriumVDof + equilibriumRDof;	      
	      
	      //method LB quantum level	      
	      cloud_.postReactionVibrationalRedistribution
	      (
	       mode,
	       remainDof,
	       theta,
	       vibLevel,
	       equilibriumE
	      );
	    }
	  }
	  else
	  {
	    continue;
	  }
	}//end for
	
	//add other non participate vibrational energy
	equilibriumE    += equilibriumVE;
	equilibriumRDof += equilibriumVDof;
	
	
	scalar ERotProductMole   = 0.0;
	scalar ERotProductSecond = 0.0;
	if(moleRDof>0.0 && secondRDof>0.0)
	{
	  equilibriumRDof -= moleRDof;
	  const scalar energyRatioMole = cloud_.postCollisionRotationalEnergy( moleRDof, equilibriumRDof/2.0);
	  ERotProductMole              = energyRatioMole*equilibriumE;
	  equilibriumE                -= ERotProductMole;
	  ERotProductSecond            = equilibriumE;
	}
	else if(moleRDof>0.0)
	{
	  ERotProductMole              = equilibriumE;
	}
	else
	{
	  ERotProductSecond            = equilibriumE;
	}


	//////////////////- asign energy and velocity ///////////////////////////////////////
        const scalar relVelExchMol = sqrt(2.0*reactedEnergy/mRExch);//collisionEnergy LB method  //reactedEnergy present
	
	
	//////////////////////////////////////  end ///////////////////////////////////////
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


	//Info << "Pbefore(exchange) = " << p.typeId() << endl;
	//Info << "Qbefore(exchange) = " << q.typeId() << endl;
	
        //- p is originally the molecule becomes the atom
        p.typeId() = typeIdSecond;
        p.U() = UP;
	if( cloud_.constProps(typeIdSecond).type() < 20 )
	{
	  p.ERot() = 0.0; 
	  p.vibLevel().setSize
	  (
	    cloud_.constProps(typeIdSecond).nVibrationalModes(),
            0
	  );
	}
	else
	{
	  //- p is molecule
	  p.ERot() = ERotProductSecond;
	  p.vibLevel() = vibLevelSecond;
	}

	//- q is originally the atom and becomes the molecule
        q.typeId() = typeIdMole;
        q.U() = UQ;
        q.ERot() = ERotProductMole;
        q.vibLevel() = vibLevelMole;

	
	//Info << "Pafter(exchange) = " << p.typeId() << endl;
	//Info << "Qafter(exchange) = " << q.typeId() << endl;
	
	/*
	scalar postE =
	   0.5*cloud_.constProps(p.typeId()).mass()*magSqr(p.U())
	  +0.5*cloud_.constProps(q.typeId()).mass()*magSqr(q.U())
	  +p.ERot()+q.ERot()
	  +cloud_.constProps(p.typeId()).eVib_tot(p.vibLevel())
	  +cloud_.constProps(q.typeId()).eVib_tot(q.vibLevel());
	Info << "postE = " << postE << endl;
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
    posAtomReactant_(-1),
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

void exchangeQK::reaction
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

void exchangeQK::reaction(dsmcParcel& p, dsmcParcel& q)
{
    //- Reset the relax switch
    relax_ = true;
    
    const label typeIdP = p.typeId();
    const label typeIdQ = q.typeId();
    
    //- Exchange reaction ABD + C --> productI + productII
    //  If Q is the second reactant C (i.e., the atom)
    //  Q is necessarily atom otherwise this class would not have been selected
    if
    (
     (     typeIdP == reactantIds_[0] ) // match order of chemicalDict
     //  or
     // (     (cloud_.constProps(typeIdP).type() >= 20)     //molecule
     //  and  (cloud_.constProps(typeIdQ).type() >= 20)  //molecule
     // )
    ) 
    {
      /*
      scalar preE =
	 0.5*cloud_.constProps(p.typeId()).mass()*magSqr(p.U())
	+0.5*cloud_.constProps(q.typeId()).mass()*magSqr(q.U())
	+p.ERot()+q.ERot()
	+cloud_.constProps(p.typeId()).eVib_tot(p.vibLevel())
	+cloud_.constProps(q.typeId()).eVib_tot(q.vibLevel())
	+heatOfReactionExchangeJoules_[0];
      Info << "preE = " << preE << endl;
      */
      
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
	
        //- Possible reactions:
        // 1. Exchange reaction	
	// random choose one reaction
	/*
        scalar reactionProbability = 0.0;
	const label rn = cloud_.randomLabel(0, numberOfExchange_-1); 

	//method 1:random choose one vibrational mode probability
	
	const label numberOfVibMode = p.vibLevel().size() + q.vibLevel().size();
	const label randomSelectMode = cloud_.randomLabel(0, numberOfVibMode-1);	
	
	if(p.vibLevel().size() > 0 && randomSelectMode < p.vibLevel().size())
	{
	  const label selectPMode = randomSelectMode;

	  testExchange
	  (
	   p,	   
	   translationalEnergy,
	   omegaPQ,
	   rn,
	   selectPMode,
	   reactionProbability
	  );

	  if(reactionProbability > cloud_.rndGen().sample01<scalar>())
	  {	    	    
	    exchange(p, q, selectPMode, rn, translationalEnergy);

	  }
	}
	else
	{
	  const label selectQMode = randomSelectMode-p.vibLevel().size();
	  
	  testExchange
	  (
	   q,	   
	   translationalEnergy,
	   omegaPQ,
	   rn,
	   selectQMode,
	   reactionProbability
	  );

	  if(reactionProbability > cloud_.rndGen().sample01<scalar>())
	  {    	    
	    exchange(q, p, selectQMode, rn, translationalEnergy);
	  }
	}
	*/

	scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(numberOfExchange_, 0.0); 

	const label numberOfVibMode = p.vibLevel().size() + q.vibLevel().size();
	const label randomSelectMode = cloud_.randomLabel(0, numberOfVibMode-1);

	label selectPMode = 0;
	label selectQMode = 0;
	bool isPreact = false;
	      
	if(p.vibLevel().size() > 0 && randomSelectMode < p.vibLevel().size())
	{
	  isPreact = true;
	  selectPMode = randomSelectMode;
	  
	  forAll(reactionProbabilities, i)
	  {	    
	    testExchange
	    (
	     p,
	     translationalEnergy,
	     omegaPQ,
	     i,
	     selectPMode,
	     reactionProbabilities[i]
	    );
	    
	    totalReactionProbability += reactionProbabilities[i];
	  }
	}
	else
	{
	  selectQMode = randomSelectMode-p.vibLevel().size();
	  
	  forAll(reactionProbabilities, i)
	  {	    
	    testExchange
	    (
	     q,
	     translationalEnergy,
	     omegaPQ,
	     i,
	     selectQMode,
	     reactionProbabilities[i]
	    );

	    totalReactionProbability += reactionProbabilities[i];
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
		//- Exchange reaction

		if( isPreact )
		{        
		  exchangeQK::exchange
		  (
		   p, q, selectPMode, i, translationalEnergy
		   ); 	  
		}
		else
		{        
		  exchangeQK::exchange
		  (
		   q, p, selectQMode, i, translationalEnergy
		   );	  
		}
		//- There can't be another reaction: break
		break;
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
        // match order of chemicalDict
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
    /*
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
    */

    nExchangeReactionsPerTimeStep_ = 0;
}
  
}
// End namespace Foam

// ************************************************************************* //
