/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

Class
    noTimeCounter

Description

\*----------------------------------------------------------------------------*/

#include "noTimeCounter.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noTimeCounter, 0);

addToRunTimeSelectionTable(collisionPartnerSelection, noTimeCounter, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
noTimeCounter::noTimeCounter
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    collisionPartnerSelection(mesh, cloud, dict),
    infoCounter_(0)
//     propsDict_(dict.subDict(typeName + "Properties"))
{}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noTimeCounter::~noTimeCounter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void noTimeCounter::initialConfiguration()
{

}

void noTimeCounter::collide()
{
    if (!cloud_.binaryCollision().active())
    {
        return;
    }

    // Temporary storage for subCells
    List< DynamicList<label> > subCells(8);    

    label collisionCandidates = 0;

    label collisions = 0;

    const List<DynamicList<dsmcParcel*>>& cellOccupancy = cloud_.cellOccupancy();

    const polyMesh& mesh = cloud_.mesh();
   
    forAll(cellOccupancy, cellI)
    {
        const scalar deltaT = cloud_.deltaTValue(cellI);
        
        const DynamicList<dsmcParcel*>& cellParcels(cellOccupancy[cellI]);
        
        const scalar& cellVolume = mesh.cellVolumes()[cellI];

        const label nC(cellParcels.size());

        if (nC > 1)
        {
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // Assign particles to one of 8 Cartesian subCells

            // Clear temporary lists
            forAll(subCells, i)
            {
                subCells[i].clear();
            }
	    
            // Inverse addressing specifying which subCell a parcel is in
            List<label> whichSubCell(cellParcels.size());

            const point& cC = mesh.cellCentres()[cellI];

            forAll(cellParcels, i)
            {
                const dsmcParcel& p = *cellParcels[i];

                vector relPos = p.position() - cC;

                label subCell =
		  pos(relPos.x()) + 2*pos(relPos.y()) + 4*pos(relPos.z());// pos: positive =1 negetive = 0

		subCells[subCell].append(i);

                whichSubCell[i] = subCell;
            }
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            scalar sigmaTcRMax = cloud_.sigmaTcRMax()[cellI];
            
            //scalar selectedPairs = 0.0;
            
            scalar selectedPairs =
                cloud_.collisionSelectionRemainder()[cellI]
                + 0.5*nC*(nC - 1)*cloud_.nParticles(cellI, true)*sigmaTcRMax*deltaT
                /cellVolume;
               
            const label nCandidates(selectedPairs);

            cloud_.collisionSelectionRemainder()[cellI] = selectedPairs - nCandidates;

            collisionCandidates += nCandidates;

           // list of candidates in cell
            DynamicList<label> candidateList(0);

            for (label c = 0; c < nC; c++)
            {
               candidateList.append(c);
            }
            candidateList.shrink();
	    
	    // list of candidates subcells
            List<DynamicList<label> > candidateSubList(subCells);
	    
	    // used for delet parcel, the original size of list 
	    //DynamicList<label> removeParcelId(0);  
	    
            for (label c = 0; c < nCandidates; c++)
            {
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // subCell candidate selection procedure
	      const label numberOfC = candidateList.size();
	      
                // Select the first collision candidate
	      const label candidateListPIndex = cloud_.randomLabel(0,numberOfC-1);
	      const label candidateP = candidateList[candidateListPIndex];
	      
	      // Declare the second collision candidate
	      label candidateQ = -1;
	      
	      DynamicList<label>& subCellPs = candidateSubList[whichSubCell[candidateP]];
	      
	      label candidateListQIndex = -1;
	      label nSC = subCellPs.size();	      
	      if (nSC > 1)
              {
		// If there are two or more particle in a subCell, choose
		// another from the same cell.  If the same candidate is
		// chosen, choose again.
		do
		{
		  //candidateListQIndex = cloud_.randomLabel(0, nSC-1);
		  candidateQ = subCellPs[cloud_.randomLabel(0, nSC-1)];
		  
		} while (candidateP == candidateQ);
		
		for(label i=candidateQ; i>-1; i--)
		{
		  if(candidateList[i] == candidateQ)
		  {
		    candidateListQIndex = i;
		    break;
		  }
		}
		
		if(candidateListQIndex == -1)
		{
		  Info << "Not find Q particle !" << endl;
		}

	      }
	      else //only himself
              {
		// Select a possible second collision candidate from the
		// whole cell.  If the same candidate is chosen, choose
		// again.		
		do
		{
		  candidateListQIndex = cloud_.randomLabel(0, numberOfC-1);
		  candidateQ = candidateList[candidateListQIndex];
		  
		} while (candidateP == candidateQ);
	      }
	      
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // uniform candidate selection procedure

                // // Select the first collision candidate
                // label candidateP = cloud_.randomLabel(0, nC-1);

                // // Select a possible second collision candidate
                // label candidateQ = cloud_.randomLabel(0, nC-1);

                // // If the same candidate is chosen, choose again
                // while (candidateP == candidateQ)
                // {
                //     candidateQ = cloud_.randomLabel(0, nC-1);
                // }
	      
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	      
		dsmcParcel& parcelP = *cellParcels[candidateListPIndex];
                dsmcParcel& parcelQ = *cellParcels[candidateListQIndex];

                label chargeP = -2;
                label chargeQ = -2;

                chargeP = cloud_.constProps(parcelP.typeId()).charge();
                chargeQ = cloud_.constProps(parcelQ.typeId()).charge();
	
                //do not allow electron-electron collisions		
                if(!(chargeP == -1 && chargeQ == -1))
                {
                    scalar sigmaTcR = cloud_.binaryCollision().sigmaTcR
                    (
                        parcelP,
                        parcelQ
                    );
                    
                    // Update the maximum value of sigmaTcR stored, but use the
                    // initial value in the acceptance-rejection criteria because
                    // the number of collision candidates selected was based on this
                    if (sigmaTcR > cloud_.sigmaTcRMax()[cellI])
                    {
                        cloud_.sigmaTcRMax()[cellI] = sigmaTcR;
                    }
	
                    if ((sigmaTcR/sigmaTcRMax) > rndGen_.sample01<scalar>())
                    {
                        // chemical reactions
                        // find which reaction model parcel p and q should use
                        label rMId = cloud_.reactions().returnModelId(parcelP, parcelQ);

    //                             Info << " parcelP id: " <<  parcelP.typeId() 
    //                                 << " parcelQ id: " << parcelQ.typeId()
    //                                 << " reaction model: " << rMId
    //                                 << endl;
                        if(rMId != -1)
                        {
			  // try to react molecules
			  if(cloud_.reactions().reactions()[rMId]->reactWithLists())
			  {
			    ///////////////////////////////////////////////////////////////////////////
			    label candidateThird          = -1;
			    label candidateListThirdIndex = -1;
			    if (nSC > 2) // p, q, thirdBody all in subcell 
			    {
			      do
			      {
				//candidateListThirdIndex = cloud_.randomLabel(0, nSC-1);
				candidateThird = subCellPs[cloud_.randomLabel(0, nSC-1)];
				
			      } while (candidateP == candidateThird || candidateQ == candidateThird);			    
			      
			      for(label i=candidateThird; i>-1; i--)
			      {
				if(candidateList[i] == candidateThird)
				{
				  candidateListThirdIndex = i;
				  break;
				}
			      }			      
			      
			      if(candidateListThirdIndex == -1)
			      {
				Info << "Not find Third particle !" << endl;
			      }
			      
			    }			    
			    else if(numberOfC > 2) // thridBody not in subcell and numberOfC > 2
			    {
			      do
			      {
				candidateListThirdIndex = cloud_.randomLabel(0, numberOfC-1);
				candidateThird = candidateList[candidateListThirdIndex];
				
			      } while (candidateP == candidateThird || candidateQ == candidateThird);
			    }

			    
			    if( candidateListThirdIndex == -1 )
			    {
			      cloud_.reactions().reactions()[rMId]->reaction
			      (
			       parcelP,
			       parcelQ			  
			      );
			    }
			    else
			    {
			      dsmcParcel& parcelThirdBody = *cellParcels[candidateListThirdIndex];		    
	        
			      // so far for recombination only
			      cloud_.reactions().reactions()[rMId]->reaction
			      (
			       parcelP,
			       parcelQ,
			       parcelThirdBody
			      );        
			    }
			    
			    /*
			    Info << "after = " << endl;
			    Info << "P Type                = " << parcelP.typeId() << endl;
			    Info << "Q Type                = " << parcelQ.typeId() << endl;
			    Info << "third Type            = " << parcelThirdBody.typeId() << endl;
			    Info << "candidateP            = " << candidateP << endl;
			    Info << "candidateQ            = " << candidateQ << endl;
			    Info << "candidateListSize     = " << candidateList.size() << endl;
			    Info << "cellParcelsSize       = " << cellParcels.size() << endl;
			    */
			   
			    //recombinationa occure and delete particle in list
			    //if(!cloud_.reactions().reactions()[rMId]->relax())
			    if(parcelP.typeId() == -1 || parcelQ.typeId() == -1)
			    {
			      /*
			      Info << "P Type           = " << parcelP.typeId() << endl;
			      Info << "Q Type           = " << parcelQ.typeId() << endl;
			      Info << "candidateP       = " << candidateP << endl;
			      Info << "candidateQ       = " << candidateQ << endl;
			      Info << "candidateList    = " << candidateList.size() << endl;
			      Info << "candidateSubList = " << candidateSubList << endl;
			      Info << "whichSubCell     = " << whichSubCell << endl;
			      Info << "/////////////////////////" << endl;
			      */

			      label deleteCandidate = 0;
			      if(parcelP.typeId() == -1)
			      {
				deleteCandidate = candidateP;
			      }
			      else
			      {
				deleteCandidate = candidateQ;
			      }
			      
			      const label deleteCandidateIndex = findIndex(candidateList, deleteCandidate);
			     
			      DynamicList<label> newCandidateList(0);
			      forAll(candidateList, i)
			      {
				if(i != deleteCandidateIndex)
				{
				  newCandidateList.append(candidateList[i]);
				}
			      }
			      candidateList.transfer(newCandidateList);
			      candidateList.shrink();

			      DynamicList<label>& subCellDelete = candidateSubList[whichSubCell[deleteCandidate]]; 
			      DynamicList<label> newSubCellPs(0);
			      const label newIndex = findIndex(subCellDelete, deleteCandidate);
			      
			      forAll(subCellDelete, i)
			      {
				if(i != newIndex)
				{
				  newSubCellPs.append(subCellDelete[i]);
				}
			      }
			      subCellDelete.transfer(newSubCellPs);
			      subCellDelete.shrink();

			      /*
			      Info << "after after = " << endl;
			      Info << "P Type             = " << parcelP.typeId() << endl;
			      Info << "Q Type             = " << parcelQ.typeId() << endl;
			      Info << "third Type         = " << parcelThirdBody.typeId() << endl;
			      Info << "candidateP         = " << candidateP << endl;
			      Info << "candidateQ         = " << candidateQ << endl;
			      Info << "candidateListSize  = " << candidateList.size() << endl;
			      Info << "cellParcelsSize    = " << cellParcels.size() << endl;
			      */
			    }
			  }
			  else
			  {
			    cloud_.reactions().reactions()[rMId]->reaction
			    (
			     parcelP,
			     parcelQ			     
			    );

			    /*
			    //if post collision energy is negtive, reselect Q particle and
			    //recompute with samd post choosen vibrational level
			    DynamicList<label> recomputeList(0);  		      			    
			    cloud_.reactions().reactions()[rMId]->reaction
			    (
			     parcelP,
			     parcelQ,
			     recomputeList 
			    );
			    
			    while(recomputeList.size() != 0)
			    {
			      label newQ = 0;
			    if (nSC > 2) // p, q, newQ all in subcell 
			    {				
			      do
			      {
				newQ = subCellPs[cloud_.randomLabel(0, nSC-1)];
				
			      } while (candidateP == candidateThird || candidateQ == candidateThird);

			      
			      
			      
			    }			    
			    else// subcell only 1 or 2 particle
			    {
			      do
			      {
				candidateThird = candidateList[cloud_.randomLabel(0, numberOfC-1)];
				
			      } while (candidateP == candidateThird || candidateQ == candidateThird);
			    }
			    }
			    */
			  }
				//                         }
				// if reaction unsuccessful use conventional collision model
			  if(cloud_.reactions().reactions()[rMId]->relax())
                          {
			    cloud_.binaryCollision().collide
			      (
			       parcelP,
			       parcelQ,
			       cellI
			      );
			  }
                        }
                        else // if reaction model not found, use conventional collision model
                        {
                            cloud_.binaryCollision().collide
                            (
                                parcelP,
                                parcelQ,
                                cellI
                            );
                        }
			
                        collisions++;
                    }
                }	
            }

        }
    }

    //Pout << "collisions individual1 = " << collisions << endl;
    reduce(collisions, sumOp<label>());
    
    reduce(collisionCandidates, sumOp<label>());
  
    cloud_.sigmaTcRMax().correctBoundaryConditions();
  
    infoCounter_++;

    if(infoCounter_ >= cloud_.nTerminalOutputs())
    {
        if (collisionCandidates)
        {
            Info<< "    Collisions                      = "
                << collisions << nl
    //             << "    Acceptance rate                 = "
    //             << scalar(collisions)/scalar(collisionCandidates) << nl
                << endl;
                
            infoCounter_ = 0;
        }
        else
        {
            Info<< "    No collisions" << endl;
            
            infoCounter_ = 0;
        }
    }
    
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
