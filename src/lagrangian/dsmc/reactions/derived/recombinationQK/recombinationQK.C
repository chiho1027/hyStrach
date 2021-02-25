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

#include "recombinationQK.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(recombinationQK, 0);

addToRunTimeSelectionTable(dsmcReaction, recombinationQK, dictionary);
  

// * * * * * * * * * * *  Protected Member functions * * * * * * * * * * * * //
  
void recombinationQK::setProperties()
{
  dsmcReaction::setProperties();

  // reactant could be atom or molecule
  if(reactantIds_.size() != 2)
  {
    FatalErrorIn("recombinationQK::setProperties()")
      << "There should be two reactants, instead of " 
      << reactantIds_.size() << nl 
      << exit(FatalError);
  }

  // reading in recombination product
  
  const word productMolecule (propsDict_.lookup("recombinationProducts"));
  
  productIdRecombination_ = findIndex(cloud_.typeIdList(), productMolecule);

  //const label totalDissModes = cloud_.constProps(productIdRecombination_).nVibrationalModes();

  // chose randomly of heat of recombination
  heatOfReactionRecombinationJoules_ =
    physicoChemical::k.value()
    *cloud_.constProps(productIdRecombination_).thetaD()[0];

  // check that heat of recombination is positive
  if(heatOfReactionRecombinationJoules_ < 0.0)
  {
    FatalErrorIn("recombinationQK::setProperties()")
      << " heat of recombination must be positive  " << heatOfReactionRecombinationJoules_ << nl 
      << exit(FatalError);
  }

  // check that reactants belong to the typeIdList (constant/dsmcProperties)
  if(productIdRecombination_ == -1)
  {
    FatalErrorIn("recombinationQK::setProperties()")
      << "Cannot find type id: " << productMolecule << nl 
      << exit(FatalError);
  }

  // check that the product is a 'MOLECULE' (not an 'ATOM') 
  if(cloud_.constProps(productIdRecombination_).type() < 20)
  {
    FatalErrorIn("recombinationQK::setProperties()")
      << "Recombination product must be a molecule (not an atom): " << productMolecule 
      << nl 
      << exit(FatalError);
  }

  const wordList thirdBody(propsDict_.lookup("thirdBody"));
  thirdBodyId_.setSize(thirdBody.size(), -1);
  
  forAll(thirdBodyId_, n)
  {
    thirdBodyId_[n] = findIndex(cloud_.typeIdList(), thirdBody[n]);

    //- Check that reactants belong to the typeIdList as defined in 
    //  constant/dsmcProperties
    if (thirdBodyId_[n] == -1)
    {
      FatalErrorIn("recombinaionQK::setProperties()")
	<< "Cannot find type id: " << thirdBody[n] << nl 
	<< exit(FatalError);
    }
  }

  //initialized other variable
  rhoC_.setSize(thirdBodyId_.size(), 0.0);
  recombinationStr_.setSize(thirdBodyId_.size(), word::null);
  nTotRecombinationReactions_.setSize(thirdBodyId_.size(), 0);
  nRecombinationReactionsPerTimeStep_.setSize(thirdBodyId_.size(), 0);

}// end setProperty

  /*
label recombinationQK::selectThirdBody()
{
  //the number of thirdBody which react with
  label nR = 0;
  const List< DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();
  const DynamicList<dsmcParcel*>& cellParcels(cellOccupancy[p.cell()]);

      // If there are two or more particle in a subCell, choose
      // another from the same cell.  If the same candidate is
      // chosen, choose again.
      dsmcParcel* candidate = parcelsInCell[cloud_.randomLabel(0, parcelsInCell.size() -1)];//error

  p.cell()
  
  label candidateId = candidate->typeId();
  nR = findIndex( thirdBodyId_, candidateId );
  
  return nR;
}
  */

scalar recombinationQK::computeDensity()
{
  scalar numberDensity = 0.0;
  const scalar volume = volume_;
  const List< DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();

  label totParticle = 0;
  forAll(cellOccupancy, c)
  {    
    totParticle += cellOccupancy[c].size();
  }
  /*
  volume_ = 0.0;
  forAll(cellOccupancy, c)
  {
    volume_ += mesh_.cellVolumes()[c];
  }

  scalar volume = volume_;
  label  molsC  = 0;
  */
  /*
  forAll(cellOccupancy, c)
  {
    const List<dsmcParcel*>& parcelsInCell = cellOccupancy[c];
    
    forAll(parcelsInCell, pIC)
    {
      dsmcParcel* p = parcelsInCell[pIC];

      if(p->typeId() == thirdBodyId_[nR])
      {
	molsC++;
      }
    }
  }
  
  //- Parallel communication
  if(Pstream::parRun())
  {
    reduce(volume, sumOp<scalar>());
  }
  */

  numberDensity = (totParticle*cloud().nParticle())/volume;
  
  return numberDensity;
  //rhoC_[nR] = (molsC*cloud().nParticle())/volume;
}

void recombinationQK::postReactionVibrationalRedistribution
(
 const scalarList thetaVProduct,
 const scalar TMacro,
 const scalar recombOmega,
 labelList& vibLevel,
 scalar& Ec
)
{
  if( thetaVProduct.size() > 0 )
  { 
    forAll(vibLevel, m)
    {
      const scalar thetaPrim = thetaVProduct[m]/TMacro;      
      const scalar kBByThetaVP = physicoChemical::k.value()*thetaVProduct[m];
      const label  iMax = Ec/kBByThetaVP;

      //Info << "Ec     = " << Ec << endl;
      //Info << "theaV  = " << thetaVProduct[m] << endl;
      //Info << "iMax   = " << iMax << endl;
      //Info << "TMacro = " << TMacro << endl;
      
      label j = 0;
      scalar func = 0.0;
      do // acceptance - rejection
      {
	j = cloud_.randomLabel(0, iMax);
	const scalar funcMax = normalizedIncomGamma( 2.5-recombOmega, iMax*thetaPrim)*exp(-j*thetaPrim);
	
	func = normalizedIncomGamma( 2.5-recombOmega, (iMax-j)*thetaPrim )
	  *exp(-j*thetaPrim)
	  /funcMax;

	//Info <<"funcM = "<< funcMax << endl;
	//Info <<" j    = "<< j << endl;
	//Info <<" func = "<< func << endl;
	
      }	while(func < cloud_.rndGen().sample01<scalar>());

      vibLevel[m] = j;
      Ec -= j*kBByThetaVP;
      
    } // end for 
  } // end if 
  else
  {
    FatalErrorIn("recombinationQK::postReactionVibrationRedistribution()");
  }
}

scalar recombinationQK::normalizedIncomGamma
(
 const scalar a,
 const scalar x
)
{
  scalar result = 0.0;
  scalar c      = 1.0;
  scalar cSum   = 1.0;
  scalar r      = a;
  scalar error  = 1e-10;

  /*
  if(x >10.0)
  {
    result = exp(-x)*pow(x, a-1)
      *(1.0 + (a-1.0)/x + (a-1.0)*(a-2.0)/pow(x, 2.0));

    if(result < VSMALL)
      result = VSMALL;

    return result;
  }
  */
  
  do
  {
    r    += 1.0;
    c    *= x/r;
    cSum += c;

    if(c > 1e304)
      break;
  } while( c/cSum > error );
  
  result = 1.0 - exp(-x)*pow(x, a)/exp(lgamma(a+1.0))*cSum;
  
  return result;
}

void recombinationQK::testRecombination
(
     const dsmcParcel& p,
     const dsmcParcel& q,
     const scalar translationalEnergy,
     const scalar omegaPQ,
     const label  nR,
     scalar& reactionProbability
)
{  
  // compute overall number density
  const scalar overallNumberDensity = computeDensity();
  
  const label typeIdP = p.typeId();
  const label typeIdQ = q.typeId();

  const scalar dP = cloud_.constProps(typeIdP).d();
  const scalar dQ = cloud_.constProps(typeIdQ).d();
  const scalar dThirdBody = cloud_.constProps(thirdBodyId_[nR]).d();
  const scalar VRef = (pi/6.0)*(pow((dP+dQ+dThirdBody),3.0));

  scalar VColl = 0.0;
  if(bCoeff_[nR] < -1.5)
  {
    scalar TMacro = cloud_.fields().overallT(p.cell());
    
    if(TMacro == 0.0)
    {
      cloud_.fields().calculateFields();
      TMacro = cloud_.fields().overallT(p.cell());
    }
    
    VColl = aCoeff_[nR]*pow(TMacro/273.0, bCoeff_[nR])*VRef;
  }
  else
  {
    //- Collision temperature: Eq.(10) of Bird's QK paper.
    const scalar TColl =
      (translationalEnergy/physicoChemical::k.value())
      /(2.5 - omegaPQ);
    
    const scalar aPrime =
      aCoeff_[nR]
      *(
	pow((2.5 - omegaPQ), bCoeff_[nR])
	*exp(lgamma(2.5 - omegaPQ))
	/exp(lgamma(2.5 - omegaPQ + bCoeff_[nR]))
	);

    VColl = aPrime*pow(TColl/273.0, bCoeff_[nR])*VRef;
  }
  /*
  if(mesh_.time().value() == mesh_.time().deltaTValue())
  {
    cloud_.evolve_fields();
  }
  */


  /*
  Info << "VRef  = " <<VRef << endl;
  Info << "VColl = " <<VColl << endl;
  Info << "rho   = " << rhoC_[nR] << endl;
  Info << "p     = " << rhoC_[nR]*VColl << endl;
  */
  
  //reactionProbability = rhoC_[nR]*VColl;
  reactionProbability = overallNumberDensity*VColl;

  if(reactionProbability > 1.0)
    {
      reactionProbability = 1.0;
    }
}// end test

void recombinationQK::recombination
(
   dsmcParcel& p,
   dsmcParcel& q,
   const label nR,
   const scalar translationalEnergy
)
{
  const label typeIdP = p.typeId();
  const label typeIdQ = q.typeId();
  
  nTotRecombinationReactions_[nR]++;
  nRecombinationReactionsPerTimeStep_[nR]++;
  
  if (allowSplitting_)
  {    
    relax_ = false;
    
    // determine collisional energy
    scalar collisionEnergy = 0.0;
    
    //calculate pre-collid particles
    if( cloud_.constProps(typeIdP).type() >= 20)
    {
      const scalar ERotP = p.ERot();
      const scalar EVibP = cloud_.constProps(typeIdP).eVib_tot(p.vibLevel());
    
      if( cloud_.constProps(typeIdQ).type() >= 20)
      {
	// A, B = molecule
	const scalar ERotQ = q.ERot();
	const scalar EVibQ = cloud_.constProps(typeIdQ).eVib_tot(q.vibLevel());
	
	collisionEnergy = translationalEnergy + ERotP + ERotQ + EVibP + EVibQ + heatOfReactionRecombinationJoules_;
      }
      else
      {
	// A = molecule, B = atom
	collisionEnergy = translationalEnergy + ERotP + EVibP + heatOfReactionRecombinationJoules_;
      }
    }
    else
    {
      if( cloud_.constProps(typeIdQ).type() >= 20)
      {
	// A = atom, B = molecule
	const scalar ERotQ = q.ERot();
	const scalar EVibQ = cloud_.constProps(typeIdQ).eVib_tot(q.vibLevel());
	
	collisionEnergy = translationalEnergy + ERotQ + EVibQ + heatOfReactionRecombinationJoules_;
      }
      else
      {
	// A, B are both atom
	collisionEnergy = translationalEnergy + heatOfReactionRecombinationJoules_;
      }
    }
    
    // calculate redistribution post particle energy
    const label& typeIdRecombinedMol = productIdRecombination_;
    const scalar& omegaProduct = cloud_.constProps(typeIdRecombinedMol).omega();
    const scalarList& thetaVProduct = cloud_.constProps(typeIdRecombinedMol).thetaV();

    // calculate TMacro
    scalar TMacro = cloud_.fields().overallT(p.cell());

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // initialized post recombinationproduct viblevel
    labelList vibLevelMole(cloud_.constProps(typeIdRecombinedMol).nVibrationalModes(), 0);    
    postReactionVibrationalRedistribution
    (
     thetaVProduct,
     TMacro,
     omegaProduct,
     vibLevelMole,
     collisionEnergy
    );

    const scalar ERotProduct = collisionEnergy;
    /*
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
    */

    ////////////////////////////////////////////////
    // determine pos velocity
    //const scalar relVelExchMol = sqrt(2.0*collisionEnergy/cloud_.constProps(typeIdRecombinedMol).mass());

    // Variable Hard Sphere collision part for collision of molecules
    //const scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
    //const scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
    //const scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();    
    //const vector& postCollisionU =
    //  relVelExchMol
    // *vector
    // (
    //  cosTheta,
    //   sinTheta*cos(phi),
    //   sinTheta*sin(phi)
    //  );
    
    const scalar mP = cloud_.constProps(typeIdP).mass();
    const scalar mQ = cloud_.constProps(typeIdQ).mass();

    //postVelocisty =  Ucm
    const vector& postCollisionU = (mP*p.U() + mQ*q.U())/(mP + mQ);

    
    const vector& position = p.position();    
    label cell = -1;
    label tetFace = -1;
    label tetPt = -1;
    
    mesh_.findCellFacePt
      (
       position,
       cell,
       tetFace,
       tetPt
      );

    const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();
    const DynamicList<dsmcParcel*>& cellParcels(cellOccupancy[cell]);
   
    
    //label idP = -1;
    //label idQ = -1;
    
    scalar RWF = 1.0;

    forAll(cellParcels, i)
      {
	const dsmcParcel& parcel = *cellParcels[i];
	
	if(parcel.position() == p.position())
	  {
	    //idP = i;
	    RWF = p.RWF();
	  }
	
	if(parcel.position() == q.position())
	  {
	    //idQ = i;
	    RWF = q.RWF();
	  }
      }

    label classification = 0;
    
    if(cloud_.rndGen().sample01<scalar>() > 0.5)
      {
	classification = p.classification();
      }
    else
      {
	classification = q.classification();
      }

      /*
    p.RWF() = p.position().x();
    q.RWF() = q.position().x();
    
    cloud_.addNewParcel
    (
     position,
     postCollisionU,
     RWF,
     ERotProduct,
     0,// elevel
     cell,
     tetFace,
     tetPt,
     typeIdRecombinedMol,
     -1,
     classification,
     vibLevelMole
    );
    */
    
    //cloud_.removeParcelFromCellOccupancy( typeIdQ, q.cell() );

    // delete Particle q
    forAllIter(dsmcCloud, cloud_, c)
    {	  
      if( c().position() == q.position() )
      {
	cloud_.deleteParticle(q);
      }
    }

    //cloud_.deleteParticle(q);
    
    cloud_.reBuildCellOccupancy();

    //remark this particle to delete in no time contur
    q.typeId() = -1;
    
    //- p is originally the atom and becomes the molecule
    p.typeId() = typeIdRecombinedMol;
    p.U() = postCollisionU;
    p.ERot() = ERotProduct;
    p.vibLevel() = vibLevelMole;
    //p.ELevel() = eLevel;
    p.RWF() = RWF;
    p.classification() = classification;

    
    /*
    scalar postE = 0.5*cloud_.constProps(p.typeId()).mass()*magSqr(p.U())
      +p.ERot()
      +cloud_.constProps(p.typeId()).eVib_tot(p.vibLevel());
    Info << "postE = " << postE << endl;
    */
    
  }// end splitting
  
}// end recombination
  
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
recombinationQK::recombinationQK
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    posMolReactant_(-1),
    thirdBodyId_(),
    nTotRecombinationReactions_(),
    nRecombinationReactionsPerTimeStep_(),
    productIdRecombination_(),
    recombinationStr_(),
    aCoeff_(propsDict_.lookup("aCoeff")),
    bCoeff_(propsDict_.lookup("bCoeff")),
    heatOfReactionRecombinationJoules_(),
    volume_(0.0),
    rhoC_()
{
    forAll(mesh_.cells(), c)
    {
	volume_ += mesh_.cellVolumes()[c];
    }
    
    if(Pstream::parRun())
    {
        reduce(volume_, sumOp<scalar>());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

recombinationQK::~recombinationQK()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void recombinationQK::initialConfiguration()
{
    setProperties();
    //computeDensity();

    const word& reactantA = cloud_.typeIdList()[reactantIds_[0]];
    const word& reactantB = cloud_.typeIdList()[reactantIds_[1]];
    const word& product = cloud_.typeIdList()[productIdRecombination_];

    forAll(thirdBodyId_, n)
    {      
      const word& thirdBody = cloud_.typeIdList()[thirdBodyId_[n]];
      
      recombinationStr_[n] = "Recombination reaction " + reactantA + " + " 
        + reactantB + " + " + thirdBody + " --> " + product + " + " + thirdBody;
    }
}

bool recombinationQK::tryReactMolecules(const label& typeIdP, const label& typeIdQ)
{
    //- Function used when setting the pair addressing matrix
    const label reactantPId = findIndex(reactantIds_, typeIdP);
    const label reactantQId = findIndex(reactantIds_, typeIdQ);

    //- If both indices were found in the list of reactants, there Ids will be
    //  different from -1
    if((reactantPId != -1) && (reactantQId != -1))
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
  
void recombinationQK::reaction
(
    dsmcParcel& p,
    dsmcParcel& q,
    const label candidateId
    //const DynamicList<label>& candidateList,
    //const List<DynamicList<label> >& candidateSubList,
    //const label& candidateP,
    //const List<label>& whichSubCell
)
{
  //- Reset the relax switch
  relax_ = true;

  if(candidateId == thirdBodyId_[0])
  {
  
  const label typeIdP = p.typeId();
  const label typeIdQ = q.typeId();

  // Recombination reaction A + B + M --> AB + M
  //  If P is the first reactant A
  //  Q is necessarily B otherwise this class would not have been selected
  if(typeIdP == reactantIds_[0])
  { 
    //find thirdBody Index 
    const label  nR = findIndex( thirdBodyId_, candidateId );
    //const label nR = 0;
    
    /*
    scalar preE = 0.5*cloud_.constProps(p.typeId()).mass()*magSqr(p.U())
      +0.5*cloud_.constProps(q.typeId()).mass()*magSqr(q.U())
      +p.ERot()+q.ERot()
      +cloud_.constProps(p.typeId()).eVib_tot(p.vibLevel())
      +cloud_.constProps(q.typeId()).eVib_tot(q.vibLevel())
      +heatOfReactionRecombinationJoules_;
    Info << "preE = " << preE << endl;
    */

    const scalar mP = cloud_.constProps(typeIdP).mass();
    const scalar mQ = cloud_.constProps(typeIdQ).mass();
    const scalar mR = mP*mQ/(mP + mQ);
    const scalar cRsqr = magSqr(p.U() - q.U());
    const scalar translationalEnergy = 0.5*mR*cRsqr;
    const scalar omegaPQ =
      0.5*(
	   cloud_.constProps(typeIdP).omega()
	   + cloud_.constProps(typeIdQ).omega()
	   );
    
    //- Possible reactions:
    // 1. Recombination reaction
    scalar totalReactionProbability = 0.0;
    scalarList reactionProbabilities(1, 0.0);
    testRecombination
    (
       p,
       q,
       translationalEnergy,
       omegaPQ,
       nR,
       reactionProbabilities[0]
    );
    totalReactionProbability += reactionProbabilities[0];
    
    // Decide if an recombinaiton reaction is to occur
    if (totalReactionProbability > cloud_.rndGen().sample01<scalar>())
    {
      //recombination(p, q, translationalEnergy);
      recombination(p, q, nR, translationalEnergy);
    }
  }
  else
  {
    //- Recombination reaction  A + B + M --> AB + M
    //  If P is the second reactant B, then switch arguments in this
    //  function and P will be first
    recombinationQK::reaction(q, p, candidateId);
  }
  }
}

void recombinationQK::reaction
(
    dsmcParcel& p,
    dsmcParcel& q
)
{
}

void  recombinationQK::outputResults(const label& counterIndex)
{   
  if(writeRatesToTerminal_)
  {
    // measure density 
    const List< DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();

    volume_ = 0.0;
    forAll(cellOccupancy, c)
    {
      volume_ += mesh_.cellVolumes()[c];
    }

    scalar volume = volume_;
    scalarList molsReactants(2, 0.0);
    scalarList molsC(thirdBodyId_.size(), 0.0);
    labelList nTotRecombinationReactions = nTotRecombinationReactions_;   
    labelList nRecombinationReactionsPerTimeStep = nRecombinationReactionsPerTimeStep_;
    
    forAll(cellOccupancy, c)
    {
      const List<dsmcParcel*>& parcelsInCell = cellOccupancy[c];
      
      forAll(parcelsInCell, pIC)
      {
	dsmcParcel* p = parcelsInCell[pIC];
	
	const label pos = findIndex(reactantIds_, p->typeId());
        
	if(pos != -1)
	{
	  molsReactants[pos]++;
	}
	
	forAll(molsC, n)
	{
	  if(p->typeId() == thirdBodyId_[n])
	  {
	    molsC[n]++;
	  }
	}
      }
    }

    if (Pstream::parRun())
    {
      reduce(molsReactants[0], sumOp<label>());
      reduce(molsReactants[1], sumOp<label>());
      reduce(volume, sumOp<scalar>());
      forAll(thirdBodyId_, n)
      {
	reduce(molsC[n], sumOp<scalar>());
	reduce(nTotRecombinationReactions[n], sumOp<label>());
	reduce(nRecombinationReactionsPerTimeStep[n], sumOp<label>());
      }
    }
    
    scalarList numberDensities(2, cloud_.nParticle()/volume);
    numberDensities[0] *= molsReactants[0];
    numberDensities[1] *= molsReactants[1];

    rhoC_ = (molsC*cloud().nParticle())/volume;

    const scalar deltaT = mesh_.time().deltaT().value();
    scalarList factor(thirdBodyId_.size(), 0.0);

    if(recombinationStr_.size() > 0)
    {
      forAll(thirdBodyId_, n)
      {
	if (reactantIds_[0] == reactantIds_[1] && numberDensities[0] > 0.0 && rhoC_[n] > 0.0)
	{	  
	  factor[n] = cloud_.nParticle()/
	  (
	   counterIndex*deltaT
	   *numberDensities[0] *numberDensities[0]*rhoC_[n]
	   *volume
	  );
	}
	else if (numberDensities[0] > 0.0 && numberDensities[1] > 0.0 && rhoC_[n] > 0.0 )
	{	  
	  factor[n] = cloud_.nParticle()/
	  (
	   counterIndex*deltaT
	   *numberDensities[0]*numberDensities[1]*rhoC_[n]
	   *volume
	  );

	}

      const scalar reactionRateRecombination = factor[n]*nTotRecombinationReactions[n];
      
      Info<< recombinationStr_[n] 
	  << ", reaction rate = " << reactionRateRecombination 
	  << ", nReactions = " 
	  << nRecombinationReactionsPerTimeStep[n]
	  << endl;
      }
    }// END IF

  }
  /*
  else
  {
    labelList nTotRecombinationReactions = nTotRecombinationReactions_;   
    labelList nRecombinationReactionsPerTimeStep = nRecombinationReactionsPerTimeStep_;
 
    if(Pstream::parRun())
    {
      forAll(thirdBodyId_, n)
      {
	//- Parallel communication
	reduce(nTotRecombinationReactions[n], sumOp<label>());
	reduce(nRecombinationReactionsPerTimeStep[n], sumOp<label>());
      }
    }

    forAll(thirdBodyId_, n)
    {
      if(nTotRecombinationReactions[n] > 0)
	{
	  Info<< recombinationStr_[n]
	      << " is active, nReactions this time step = " 
	      << nRecombinationReactionsPerTimeStep[n]
	      << endl;
	}
    }
  }// end else
  */
  
  nRecombinationReactionsPerTimeStep_ = 0;
}
  
}
// End namespace Foam

