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

#include "dsmcNeumannOutlet2.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using namespace Foam::constant::mathematical;

namespace Foam
{

defineTypeNameAndDebug(dsmcNeumannOutlet2, 0);

addToRunTimeSelectionTable(dsmcGeneralBoundary, dsmcNeumannOutlet2, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcNeumannOutlet2::dsmcNeumannOutlet2
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcGeneralBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    typeIds_(),    
    cellVolume_(faces_.size(), scalar(0.0)),    
    accumulatedParcelsToInsert_(),
    moleFractions_(),
    outletVelocity_(),
    outletNumberDensity_(),
    outletTemperature_(),
    totalMoleFractions_(),    
    totalNumberDensity_(),
    totalTemperature_(),
    totalVelocity_(),
    totalRotationalEnergy_(faces_.size(), scalar(0.0)),
    totalRotationalDof_(faces_.size(), scalar(0.0)),
    totalVibrationalEnergy_(),
    vibT_(),
    vDof_(),
    sampleNumber_(),
    sampleMoleFractions_(),
    sampleVelocity_(),
    sampleTemperature_(),
    sampleNumberDensity_(),
    nTimeSteps_(0)
{
    writeInTimeDir_ = false;
    writeInCase_ = true;

    setProperties();
        
    forAll(cellVolume_, c)
    {
        cellVolume_[c] = mesh_.cellVolumes()[cells_[c]];  //get volume of each boundary cell
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcNeumannOutlet2::~dsmcNeumannOutlet2()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcNeumannOutlet2::initialConfiguration()
{}

void dsmcNeumannOutlet2::calculateProperties()
{
}

void dsmcNeumannOutlet2::controlParcelsBeforeMove()
{    
    //////////////////////////////////////////////////////////////////////
  /*
    if(nTimeSteps_ == 0)
    {
      controlParcelsAfterCollisions();
    }
  */
    //if(nTimeSteps_ > 0)
    //{
        Random& rndGen = cloud_.rndGen();
        
        //label nTotalParcelsAdded = 0;
        //label nTotalParcelsToBeAdded = 0;
 
        //forAll(moleFractions_, iD)
        //{
        //    forAll(moleFractions_[iD], c)
        //    {    
        //        if(nTotalParcels_[c] > VSMALL)
        //        {		  
        //            moleFractions_[iD][c] = nTotalParcelsSpecies_[iD][c]/nTotalParcels_[c];
        //        }
        //    }
        //}

	//scalarField averageMoleFractions(moleFractions_.size(), 0.0);
	
        labelField parcelsInserted(typeIds_.size(), 0);
	labelField parcelsToBeAdded(typeIds_.size(), 0);
	
        forAll(accumulatedParcelsToInsert_, iD) 
        {
            // loop over all faces of the patch
            forAll(accumulatedParcelsToInsert_[iD], f)
            {
                vector faceVelocity = outletVelocity_[f];
                scalar faceTemperature = outletTemperature_[f];
		
                const label& faceI = faces_[f];
                const label& cellI = cells_[f];
                const vector& fC = mesh_.faceCentres()[faceI];
                const vector& sF = mesh_.faceAreas()[faces_[f]];
		
                scalar fA = mag(sF);

                List<tetIndices> faceTets = polyMeshTetDecomposition::faceTetIndices
                (
                    mesh_,
                    faceI,
                    cellI
                );

                // Cumulative triangle area fractions
                List<scalar> cTriAFracs(faceTets.size(), 0.0);

                scalar previousCummulativeSum = 0.0;
        
                forAll(faceTets, triI)
                {
                    const tetIndices& faceTetIs = faceTets[triI];

                    cTriAFracs[triI] =
                        faceTetIs.faceTri(mesh_).mag()/fA
                        + previousCummulativeSum;

                    previousCummulativeSum = cTriAFracs[triI];
                }

                // Force the last area fraction value to 1.0 to avoid any
                // rounding/non-flat face errors giving a value < 1.0
                cTriAFracs.last() = 1.0;

                //Normal unit vector *negative* so normal is pointing into the
                //domain
                vector n = sF;
                n /= -mag(n);

                // Wall tangential unit vector. Use the direction between the
                // face centre and the first vertex in the list
                vector t1 = fC - mesh_.points()[mesh_.faces()[faceI][0]]; 
                t1 /= mag(t1);

                // Other tangential unit vector.  Rescaling in case face is not
                //flat and n and t1 aren't perfectly orthogonal
                vector t2 = n^t1;
                t2 /= mag(t2);

                label nParcelsToInsert = max(label(accumulatedParcelsToInsert_[iD][f]), 0);
                
                if ((accumulatedParcelsToInsert_[iD][f] - nParcelsToInsert) > rndGen.sample01<scalar>())
                {
                    nParcelsToInsert++;
                }
		
                accumulatedParcelsToInsert_[iD][f] -= nParcelsToInsert; //remainder has been set
		parcelsToBeAdded[iD] += nParcelsToInsert;

		Pout << "nParcelsToInsert = " << nParcelsToInsert << endl;
		Pout << "nP = " << cloud_.cellOccupancy()[cellI-4].size() << endl;
		  
		//nTotalParcelsToBeAdded += nParcelsToInsert;

                const label& typeId = typeIds_[iD];
                scalar mass = cloud_.constProps(typeId).mass();

		vector Usum = vector::zero; 
		
                for (label i = 0; i < nParcelsToInsert; i++)
                {
                    // Choose a triangle to insert on, based on their relative
                    // area

                    scalar triSelection = rndGen.sample01<scalar>();

                    // Selected triangle
                    label selectedTriI = -1;

                    forAll(cTriAFracs, triI)
                    {
                        selectedTriI = triI;

                        if (cTriAFracs[triI] >= triSelection)
                        {
                            break;
                        }
                    }

                    // Randomly distribute the points on the triangle.

                    const tetIndices& faceTetIs = faceTets[selectedTriI];

                    point p = faceTetIs.faceTri(mesh_).randomPoint(rndGen);
                    
                    // Velocity generation
                    scalar mostProbableSpeed
                    (
                        cloud_.maxwellianMostProbableSpeed
                        (
                            faceTemperature,
                            mass
                        )
                    );

		    
		    //new		    
		    const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();
		    const List<dsmcParcel*>& parcelsInCell = cellOccupancy[cellI-2];

		    DynamicList<dsmcParcel*> parcelIsTypeId(0);
		    
		    forAll(parcelsInCell, pIC)
		    {
		      dsmcParcel* p = parcelsInCell[pIC];
		      
		      if(p->typeId() == typeId)
		      {
			parcelIsTypeId.append(p);
		      }
		    }
		    		    
		    const label s = cloud_.randomLabel(0, parcelIsTypeId.size()-1);

		    dsmcParcel* selectedParcel = parcelIsTypeId[s];
		    
		    vector U = selectedParcel->U();		    
		    
		    scalar ERot =selectedParcel->ERot();
		    labelList vibLevel =selectedParcel->vibLevel();
		    label ELevel =selectedParcel->ELevel();		    
		    
		    Usum += U;
		    
		    
                    scalar sCosTheta = (faceVelocity & n)/mostProbableSpeed;

                    // Coefficients required for Bird eqn 12.5
                    scalar uNormProbCoeffA = sCosTheta + sqrt(sqr(sCosTheta) + 2.0);

                    scalar uNormProbCoeffB =
                        0.5*
                        (
                            1.0
                            + sCosTheta*(sCosTheta - sqrt(sqr(sCosTheta) + 2.0))
                        );

                    // Equivalent to the QA value in Bird's DSMC3.FOR
                    scalar randomScaling = 3.0;

                    if (sCosTheta < -3)
                    {
                        randomScaling = mag(sCosTheta) + 1;
                    }

                    scalar P = -1;

                    // Normalised candidates for the normal direction velocity
                    // component
                    scalar uNormal;
                    scalar uNormalThermal;
        
                    if(abs(faceVelocity & n) > VSMALL)
                    {
                        // Select a velocity using Bird eqn 12.5
                        do
                        {
                            uNormalThermal =
                                randomScaling*(2.0*rndGen.sample01<scalar>() - 1);

                            uNormal = uNormalThermal + sCosTheta;

                            if (uNormal < 0.0)
                            {
                                P = -1;
                            }
                            else
                            {
                                P = 2.0*uNormal/uNormProbCoeffA
                                    *exp(uNormProbCoeffB - sqr(uNormalThermal));
                            }

                        } while (P < rndGen.sample01<scalar>());
                    }
                    else
                    {
                        uNormal = sqrt(-log(rndGen.sample01<scalar>()));
                    }		    		    		    
		    
		    
		    /*
                    vector U =
                        sqrt(physicoChemical::k.value()*faceTemperature/mass)
                        *(
                            rndGen.GaussNormal<scalar>()*t1
                            + rndGen.GaussNormal<scalar>()*t2
                        )
                        + (t1 & faceVelocity)*t1
                        + (t2 & faceVelocity)*t2
                        + mostProbableSpeed*uNormal*n;

		    Usum += U;		    
		    
                    scalar ERot = cloud_.equipartitionRotationalEnergy
                    (
                        faceTemperature,
                        cloud_.constProps(typeId).rotationalDegreesOfFreedom()
                    );
                    
                    labelList vibLevel = cloud_.equipartitionVibrationalEnergyLevel
                    (
                        faceTemperature,
                        cloud_.constProps(typeId).nVibrationalModes(),
                        typeId
                    );
                    
                    label ELevel = cloud_.equipartitionElectronicLevel
                    (
                        faceTemperature,
                        cloud_.constProps(typeId).electronicDegeneracyList(),
                        cloud_.constProps(typeId).electronicEnergyList()
                    );
		    */

                    label newParcel = patchId();                    
		      
                    const scalar& RWF = cloud_.coordSystem().RWF(cellI);
		 
                    cloud_.addNewParcel
                    (
                        p,
                        U,
                        RWF,
                        ERot,
                        ELevel,
                        cellI,
                        faces_[f],
                        faceTetIs.tetPt(),
                        typeId,
                        newParcel,
                        0,
                        vibLevel
                    );

                    //nTotalParcelsAdded++;
                    parcelsInserted[iD]++;
                }

		if(parcelsInserted[iD] > 0)
		{
		  Usum /= parcelsInserted[iD];
		  Pout << "pracelsInserted[iD] = " << parcelsInserted[iD] << endl;
		  Pout << "Usum = " << Usum << endl;
		}
		
            }

	    /*
            if (Pstream::parRun())
            {
                reduce(parcelsInserted[iD], sumOp<scalar>());

                Info<< "dsmcNeumannOutlet specie: " << typeIds_[iD]
                    <<", inserted parcels: " << parcelsInserted[iD]
                    << endl;
            }
            else
            {                
                Info<< "dsmcNeumannOutlet specie: " << typeIds_[iD]
                    <<", inserted parcels: " << parcelsInserted[iD]
                    << endl;
            }
	    */
	  
        }

	/*
        forAll(moleFractions_, iD)
        {
	  if(moleFractions_[iD].size() > 0)
	  {
            forAll(moleFractions_[iD], f)
            {
                averageMoleFractions[iD] += moleFractions_[iD][f];
            }

	    //Pout << "iiiii = " << moleFractions_[iD].size() << endl;
	    
            averageMoleFractions[iD] /= moleFractions_[iD].size();
	  }
	}

        if (Pstream::parRun())
        {            
            Info << "moleFractions = " << averageMoleFractions << endl;
        }
        else
        {
            Info << "moleFractions = " << averageMoleFractions << endl;
        }
	*/

  //}

}

void dsmcNeumannOutlet2::controlParcelsBeforeCollisions()
{
    
}

void dsmcNeumannOutlet2::controlParcelsAfterCollisions()
{    
    nTimeSteps_ += 1.0;
    
    const scalar sqrtPi = sqrt(pi);
    
    //scalarField molecularMass(cells_.size(), 0.0);
    //scalarField molarcontantPressureSpecificHeat(cells_.size(), 0.0);
    //scalarField molarcontantVolumeSpecificHeat(cells_.size(), 0.0);
    //scalarField gasConstant(cells_.size(), 0.0);
    //scalarField gamma(cells_.size(), 0.0);

    //forAll(cells_, c)
    //{     
    //  forAll(moleFractions_, j)  
    //  {
    //	const label& typeId = typeIds_[j];
	
    //	molecularMass[c] += cloud_.constProps(typeId).mass()*moleFractions_[j][c];
	//molarcontantPressureSpecificHeat[c] += (5.0 + cloud_.constProps(typeId).rotationalDegreesOfFreedom())*moleFractions_[j][c];
	//molarcontantVolumeSpecificHeat[c] += (3.0 + cloud_.constProps(typeId).rotationalDegreesOfFreedom())*moleFractions_[j][c];
    //  } 
	
      // R = k/m
      //gasConstant[c] = physicoChemical::k.value()/molecularMass[c];
      
      //gamma[c] = molarcontantPressureSpecificHeat[c]/molarcontantVolumeSpecificHeat[c];
    //} 
    
    // calculate properties in cells attached to each boundary face    
    const label lastIndex = sampleNumber_-1;
    
    vectorField momentum(faces_.size(), vector::zero);
    vectorField UCollected(faces_.size(), vector::zero);  
    scalarField mass(faces_.size(), scalar(0.0));    
    scalarField mcc(faces_.size(), scalar(0.0));    
    scalarField newNParcels(faces_.size(), scalar(0.0)); 
    List <scalarField> newParcelsSpecies(typeIds_.size());
    forAll(newParcelsSpecies, iD)
    {
      newParcelsSpecies[iD].setSize(nFaces_, 0.0);
    }
    
    /*
      scalarField nParcelsInt(faces_.size(), scalar(0.0));
      scalarField rotationalEnergy(faces_.size(), scalar(0.0));
      scalarField rotationalDof(faces_.size(), scalar(0.0));    
      scalarField vDoF(faces_.size(), scalar(0.0));    
      scalarField massDensity(faces_.size(), scalar(0.0));
      scalarField pressure(faces_.size(), scalar(0.0));
      scalarField speedOfSound(faces_.size(), scalar(0.0));
      scalarField velocityCorrection(faces_.size(), scalar(0.0));
      scalarField massDensityCorrection(faces_.size(), scalar(0.0));
    */

    const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();

    forAll(cells_, c)
    {
        const label celli = cells_[c];
	
        const List<dsmcParcel*>& parcelsInCell = cellOccupancy[celli];
	
        forAll(parcelsInCell, pIC)
        {
            dsmcParcel* p = parcelsInCell[pIC];
            
            label iD = findIndex(typeIds_, p->typeId());
            
            if(iD != -1)
            {
                momentum[c] += cloud_.nParticles(celli)*cloud_.constProps(p->typeId()).mass()*p->U();
                mass[c] += cloud_.nParticles(celli)*cloud_.constProps(p->typeId()).mass();
                mcc[c] += cloud_.nParticles(celli)*cloud_.constProps(p->typeId()).mass()*mag(p->U())*mag(p->U());
		newNParcels[c] += 1.0;
                UCollected[c] += p->U();

		//rotationalEnergy[c] += p->ERot();
                //rotationalDof[c] += cloud_.constProps(p->typeId()).rotationalDegreesOfFreedom();
		
                //if(cloud_.constProps(p->typeId()).rotationalDegreesOfFreedom() > VSMALL)
                //{
                //    nParcelsInt[c] += 1.0;
                //}
                
//              totalVibrationalEnergy_[iD][c] += p->vibLevel()*physicoChemical::k.value()*cloud_.constProps(p->typeId()).thetaV();
                newParcelsSpecies[iD][c] += 1.0;
            }
        }

	// delete the first collume value in sample data and add newTime property into laste collume zero for default
	forAll(sampleMoleFractions_, iD)
	{  	    
	  totalMoleFractions_[iD][c] -= sampleMoleFractions_[iD][c][0];
	  totalTemperature_[c]       -= sampleTemperature_[c][0];
	  totalNumberDensity_[c]     -= sampleNumberDensity_[c][0];
	  totalVelocity_[c]          -= sampleVelocity_[c][0];
	    
	  forAll(sampleMoleFractions_[iD][c], s)
	  {
	    if(s == sampleMoleFractions_[iD][c].size()-1)
	    {
	      sampleMoleFractions_[iD][c][s] = 0.0;
	      sampleTemperature_[c][s]       = 0.0;
	      sampleNumberDensity_[c][s]     = 0.0;
	      sampleVelocity_[c][s]          = vector::zero;
	    }
	    else
	    {
	      sampleMoleFractions_[iD][c][s] = sampleMoleFractions_[iD][c][s+1];
	      sampleTemperature_[c][s]       = sampleTemperature_[c][s+1];
	      sampleNumberDensity_[c][s]     = sampleNumberDensity_[c][s+1];
	      sampleVelocity_[c][s]          = sampleVelocity_[c][s+1];
	    }
	  }
	}
	
        //nTotalParcelsInt_[c] += nParcelsInt[c];
	
	//totalRotationalEnergy_[c] += rotationalEnergy[c];
        
        //totalRotationalDof_[c] += rotationalDof[c];

	//mcc_[c] += mcc[c];
	
        //nTotalParcels_[c] += nParcels[c];        
        
        //totalMomentum_[c] += momentum[c];
        
        //totalMass_[c] += mass[c];	
                	
        if(newNParcels[c] > 1)
        {                
	    forAll(moleFractions_, iD)
	    {		
	      scalar newMoleFractions = newParcelsSpecies[iD][c]/newNParcels[c];

	      if(nTimeSteps_ <= 2)
	      {
		newMoleFractions = sampleMoleFractions_[iD][c][0];
	      }
	      
	      totalMoleFractions_[iD][c] += newMoleFractions;
	      
	      sampleMoleFractions_[iD][c][lastIndex] = newMoleFractions;
	      
	      moleFractions_[iD][c] =
		dsmcGeneralBoundary::leastSquaresGradPredictValue(sampleMoleFractions_[iD][c], totalMoleFractions_[iD][c]);	     
	    }	   

	    vector newVelocity = momentum[c]/mass[c];      
	    
	    if(nTimeSteps_ <= 2)
	    {
	      newVelocity = sampleVelocity_[c][0];
	    }
	    
	    totalVelocity_[c] += newVelocity;

	    sampleVelocity_[c][lastIndex] = newVelocity;

	    scalar newNumberDensity = (newNParcels[c]*cloud_.nParticles(celli)) / (cellVolume_[c]);

	    if(nTimeSteps_ <= 2)
	    {
	      newNumberDensity = sampleNumberDensity_[c][0];
	    }
	    
	    totalNumberDensity_[c] += newNumberDensity;

	    sampleNumberDensity_[c][lastIndex] = newNumberDensity;

            scalar newTranslationalTemperature = (1.0/(3.0*physicoChemical::k.value()))
                  *(
		     (mcc[c]/(newNParcels[c]*cloud_.nParticles(celli)))
                      - (
			 (mass[c]/(newNParcels[c]*cloud_.nParticles(celli))
			 )*mag(newVelocity)*mag(newVelocity)
			)
		    );

	    Pout << "parcels =" << newNParcels[c] << endl;
	    
	    if(nTimeSteps_ <= 2)
	    {
	      newTranslationalTemperature = sampleTemperature_[c][0];
	    }

	    /*
            if(newTranslationalTemperature < VSMALL)
            {
                newTranslationalTemperature = 300.00;
            }
	    */

	    totalTemperature_[c] += newTranslationalTemperature;

	    sampleTemperature_[c][lastIndex] = newTranslationalTemperature;
	    
            //pressure[c] = numberDensity[c]*physicoChemical::k.value()*translationalTemperature[c];
                                    
            //speedOfSound[c] = sqrt(gamma[c]*gasConstant[c]*translationalTemperature[c]); 
            
            // Liou and Fang, 2000, equation 26 STEP 1
            //outletNumberDensity_[c] = totalMass_[c] / (molecularMass[c]*cellVolume_[c]*nTimeSteps_);

	    outletNumberDensity_[c] = dsmcGeneralBoundary::leastSquaresGradPredictValue(sampleNumberDensity_[c], totalNumberDensity_[c]);
	    //outletNumberDensity_[c] = 4.58684e24;
	    
	    outletTemperature_[c] = dsmcGeneralBoundary::leastSquaresGradPredictValue(sampleTemperature_[c], totalTemperature_[c]);
	    //outletTemperature_[c] = 1600.0;
	    
            outletVelocity_[c] = dsmcGeneralBoundary::leastSquaresGradPredictValue(sampleVelocity_[c], totalVelocity_[c]);

	    /*
	    if(c == 1)
	    {
	      Pout << "nTimeSteps = " << nTimeSteps_ << endl;
	      Pout << "newVelocity = " << newVelocity << endl;
	      Pout << "sample = " << sampleVelocity_[1] << endl;
	      Pout << "total = " << totalVelocity_[1] << endl;
	      Pout << "dsmcNeumann velocity = " << outletVelocity_[1] << endl;   
	    }
	    */
        }
    }
       
    if(faces_.size() > VSMALL)
    {
      Pout << "dsmcNeumann velocity = " << outletVelocity_[(faces_.size()/2)] << endl;      
    }
    
    forAll(accumulatedParcelsToInsert_, iD)
    {
        const label& typeId = typeIds_[iD];
        
        forAll(accumulatedParcelsToInsert_[iD], f)
        {
            const label& faceI = faces_[f];
            const vector& sF = mesh_.faceAreas()[faceI];//positive
            const scalar fA = mag(sF);
            
            const scalar deltaT = cloud_.deltaTValue(mesh_.boundaryMesh()[patchId_].faceCells()[faceI]);

            scalar mass = cloud_.constProps(typeId).mass();              

            scalar mostProbableSpeed
            (
                cloud_.maxwellianMostProbableSpeed
                (
                    outletTemperature_[f],
                    mass
                )
            );

             // Dotting boundary velocity with the face unit normal
            // (which points out of the domain, so it must be
            // negated), dividing by the most probable speed to form
            // molecularSpeedRatio * cosTheta
            
            scalar sCosTheta = (outletVelocity_[f] & -sF/fA )/mostProbableSpeed;
            
            //const scalar& RWF = cloud_.coordSystem().pRWF(patchId_, f);
            
            // From Bird eqn 4.22	    
	    
            accumulatedParcelsToInsert_[iD][f] += 
                moleFractions_[iD][f]*
                (
                    fA*outletNumberDensity_[f]*deltaT*mostProbableSpeed
                    *
                    (
                        exp(-sqr(sCosTheta)) + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                    )
                )
                /(2.0*sqrtPi*cloud_.nParticles(patchId_, f));
	    
	    Pout << "acc = " << accumulatedParcelsToInsert_[iD][f] << endl;
        } 
    }

}

void dsmcNeumannOutlet2::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
}

void dsmcNeumannOutlet2::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    //setProperties();
}

void dsmcNeumannOutlet2::setProperties()
{
    sampleNumber_ = propsDict_.lookupOrDefault<label>("sampleNumber", 10);
    
    const List<word> molecules (propsDict_.lookup("typeIds"));

    if(molecules.size() == 0)
    {
        FatalErrorIn("dsmcNeumannOutlet2::setProperties()")
            << "Cannot have zero typeIds being inserd." << nl << "in: "
            << mesh_.time().system()/"boundariesDict"
            << exit(FatalError);
    }
 
    DynamicList<word> moleculesReduced(0);

    forAll(molecules, i)
    {
        const word& moleculeName(molecules[i]);

        if(findIndex(moleculesReduced, moleculeName) == -1)
        {
            moleculesReduced.append(moleculeName);
        }
    }

    moleculesReduced.shrink();

    //  set the type ids

    typeIds_.setSize(moleculesReduced.size(), -1);

    forAll(moleculesReduced, i)
    {
        const word& moleculeName(moleculesReduced[i]);

        label typeId(findIndex(cloud_.typeIdList(), moleculeName));

        if(typeId == -1)
        {
            FatalErrorIn("dsmcNeumannOutlet2::dsmcNeumannOutlet2()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << mesh_.time().system()/"boundariesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }

    sampleTemperature_.clear();
    outletTemperature_.clear();
    totalTemperature_.clear();
    
    sampleTemperature_.setSize(nFaces_);
    outletTemperature_.setSize(nFaces_);
    totalTemperature_.setSize(nFaces_);
    const scalar temperature = readScalar(propsDict_.lookup("initialTemperature"));
    forAll(sampleTemperature_, m)
    {      
      outletTemperature_[m] = temperature;
      sampleTemperature_[m].setSize(sampleNumber_, temperature);
      totalTemperature_[m] = sum(sampleTemperature_[m]); 
    }

    sampleNumberDensity_.clear();
    outletNumberDensity_.clear();
    totalNumberDensity_.clear();
    
    sampleNumberDensity_.setSize(nFaces_);
    outletNumberDensity_.setSize(nFaces_);
    totalNumberDensity_.setSize(nFaces_);
    const scalar numberDensity = readScalar(propsDict_.lookup("initialNumberDensities"));
    forAll(sampleNumberDensity_, m)
    {
      outletNumberDensity_[m] = numberDensity;
      sampleNumberDensity_[m].setSize(sampleNumber_, numberDensity);
      totalNumberDensity_[m] = sum(sampleNumberDensity_[m]);
    }
    
    sampleVelocity_.clear();
    outletVelocity_.clear();
    totalVelocity_.clear();
  
    sampleVelocity_.setSize(nFaces_);
    outletVelocity_.setSize(nFaces_);
    totalVelocity_.setSize(nFaces_);
    const vector velocity = vector(propsDict_.lookup("initialVelocity"));
    forAll(sampleVelocity_, m)
    {
      outletVelocity_[m] = velocity;
      sampleVelocity_[m].setSize(sampleNumber_, velocity);
      totalVelocity_[m] = sum(sampleVelocity_[m]);
    }
  
    // read in the mole fraction per specie

    const dictionary& moleFractionsDict
    (
        propsDict_.subDict("moleFractions")
    );

    moleFractions_.clear();
    sampleMoleFractions_.clear();
    totalMoleFractions_.clear();

    moleFractions_.setSize(typeIds_.size());
    sampleMoleFractions_.setSize(typeIds_.size());
    totalMoleFractions_.setSize(typeIds_.size());
    forAll(moleFractions_, iD)
    {
        const scalar moleFractions = readScalar
        (
            moleFractionsDict.lookup(moleculesReduced[iD])
        );
      
        moleFractions_[iD].setSize(nFaces_);
	totalMoleFractions_[iD].setSize(nFaces_);
	sampleMoleFractions_[iD].setSize(nFaces_);

	forAll(sampleMoleFractions_[iD], c)
	{
	  moleFractions_[iD][c] = moleFractions;
	  sampleMoleFractions_[iD][c].setSize(sampleNumber_, moleFractions);
	  totalMoleFractions_[iD][c] = sum(sampleMoleFractions_[iD][c]);
	}	  
    }
        
    // set the accumulator  

    accumulatedParcelsToInsert_.setSize(typeIds_.size());

    forAll(accumulatedParcelsToInsert_, m)
    {
        accumulatedParcelsToInsert_[m].setSize(nFaces_, 0.0);
    }
    
    vibT_.setSize(typeIds_.size());

    forAll(vibT_, m)
    {
        vibT_[m].setSize(nFaces_, 0.0);
    }
    
    vDof_.setSize(typeIds_.size());

    forAll(vDof_, m)
    {
        vDof_[m].setSize(nFaces_, 0.0);
    }
   
//     totalVibrationalEnergy_.setSize(typeIds_.size());
// 
//     forAll(totalVibrationalEnergy_, m)
//     {
//         totalVibrationalEnergy_[m].setSize(nFaces_, 0.0);
//     }
    
//    nTotalParcelsSpecies_.setSize(typeIds_.size());

//    forAll(nTotalParcelsSpecies_, m)
//    {
//        nTotalParcelsSpecies_[m].setSize(nFaces_, 0.0);
//    }
           
}

void dsmcNeumannOutlet2::setNewBoundaryFields()
{
  /*
    patchId_ = mesh_.boundaryMesh().findPatchID(patchName_);

    const polyPatch& patch = mesh_.boundaryMesh()[patchId_];

    //- initialise data members
    faces_.setSize(patch.size());
    cells_.setSize(patch.size());

    //- loop through all faces and set the boundary cells
    //- no conflict with parallelisation because the faces are unique
    
    nFaces_ = 0;
    patchSurfaceArea_ = 0.0;

    for(label i = 0; i < patch.size(); i++)
    {
        label globalFaceI = patch.start() + i;

        faces_[i] = globalFaceI;
        cells_[i] = patch.faceCells()[i];
        nFaces_++;
        patchSurfaceArea_ += mag(mesh_.faceAreas()[globalFaceI]);
    }

    if(Pstream::parRun())
    {
        reduce(patchSurfaceArea_, sumOp<scalar>());
    }
    
   forAll(accumulatedParcelsToInsert_, m)
    {
        accumulatedParcelsToInsert_[m].setSize(nFaces_, 0.0);
    }
    */
}



} // End namespace Foam

// ************************************************************************* //
