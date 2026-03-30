/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

Application
    pimpleFoam

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids,
    with optional mesh motion and mesh topology changes.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "kinematicMomentumTransportModel.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "CorrectPhi.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }
   
    //const scalar NtimeStep = runTime.timeStamp; 
    const scalar firstTime = runTime.value();
    //Info<< "firstTime = "<< firstTime << nl << endl;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    
    volScalarField V
    (
         IOobject
         (
         mesh.V().name(),
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE,
         false
         ),
         mesh,
         dimensionedScalar(mesh.V().dimensions(), Zero),
         calculatedFvPatchField<scalar>::typeName
    );

    V.ref() = mesh.V();    
   
    while (pimple.run(runTime))
    {
        #include "readDyMControls.H"
	
        
        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;
	
        Info<< "Time = " << runTime.timeName() << nl << endl;
	Info<< "deltaT = "<< runTime.deltaTValue() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
            {
                fvModels.preUpdateMesh();

                mesh.update();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        #include "correctPhi.H"
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            fvModels.correct();

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }
 
	const scalar CurrentTime = runTime.value() - firstTime; // When the averaging should start from a non-zero time - usually the case
        //const scalar CurrentTime = runTime.value();
	Info<< "CurrentTime = " << CurrentTime << nl << endl;

	const scalar CurrentDeltaT = runTime.deltaTValue();

	UMean_ = ((CurrentTime-CurrentDeltaT)/CurrentTime)*UMean_ + (CurrentDeltaT/CurrentTime)*U;
        

        //// get modelled turbulent fields from RANS model
        kMod_       = turbulence->k();
        //sigmaMod_   = turbulence->sigma();
        nuMod_      = turbulence->nut(); // named nuMod_ to not confuse with nut from the turbulence model
        epsilonMod_ = turbulence->epsilon();
	//omegaMod_   = turbulence->omega(); // omega does not exist in this class
	omegaMod_   = epsilonMod_/(max(kMod_,dimensionedScalar("minkMod", kMod_.dimensions(), SMALL))*0.09); // SMALL for double precision is 1e-15

        //// time averaged modelled/SGS fields
	kMod_Mean_       = ((CurrentTime-CurrentDeltaT)/CurrentTime)*kMod_Mean_       + (CurrentDeltaT/CurrentTime)*kMod_;
        epsilonMod_Mean_ = ((CurrentTime-CurrentDeltaT)/CurrentTime)*epsilonMod_Mean_ + (CurrentDeltaT/CurrentTime)*epsilonMod_;
	omegaMod_Mean_   = ((CurrentTime-CurrentDeltaT)/CurrentTime)*omegaMod_Mean_   + (CurrentDeltaT/CurrentTime)*omegaMod_;
        nut_Mean_        = ((CurrentTime-CurrentDeltaT)/CurrentTime)*nut_Mean_        + (CurrentDeltaT/CurrentTime)*nuMod_; // changed back to 'nut' 

	// Inst. velocity
	UPrime_ = U - UMean_;

	//// Resolved TKE
	kRes_ = 0.5*Foam::magSqr(UPrime_);

	//// Resolved epsilon calculation
	volTensorField gradUPrime_ = fvc::grad(UPrime_);
	volSymmTensorField Sres_(symm(gradUPrime_));
	epsilonRes_ = 2.0*nu*Foam::magSqr(Sres_);
	omegaRes_ = epsilonRes_/(max(kRes_,dimensionedScalar("minkRes", kRes_.dimensions(), SMALL))*0.09);
  
	//// time averaged resolved fields
        kRes_Mean_       = ((CurrentTime-CurrentDeltaT)/CurrentTime)*kRes_Mean_       + (CurrentDeltaT/CurrentTime)*kRes_;
        epsilonRes_Mean_ = ((CurrentTime-CurrentDeltaT)/CurrentTime)*epsilonRes_Mean_ + (CurrentDeltaT/CurrentTime)*epsilonRes_;
        omegaRes_Mean_   = ((CurrentTime-CurrentDeltaT)/CurrentTime)*omegaRes_Mean_   + (CurrentDeltaT/CurrentTime)*omegaRes_;
        

        //// For verification: Another way of resolved epsilon calculation
	volTensorField tauPrimeByRho_ = nu*(Foam::dev2(T(fvc::grad(UPrime_))) + gradUPrime_);
        epsilonResControl_ = tauPrimeByRho_ && gradUPrime_;

	//// Total values
	kTot_       = kRes_ + kMod_;
	epsilonTot_ = epsilonRes_ + epsilonMod_;
	omegaTot_   = omegaRes_ + omegaMod_;

	//// time averaged total fields
        kTot_Mean_       = ((CurrentTime-CurrentDeltaT)/CurrentTime)*kTot_Mean_       + (CurrentDeltaT/CurrentTime)*kTot_;
        epsilonTot_Mean_ = ((CurrentTime-CurrentDeltaT)/CurrentTime)*epsilonTot_Mean_ + (CurrentDeltaT/CurrentTime)*epsilonTot_;
        omegaTot_Mean_   = ((CurrentTime-CurrentDeltaT)/CurrentTime)*omegaTot_Mean_   + (CurrentDeltaT/CurrentTime)*omegaTot_;


        //// Energy
	u2_ = Foam::sqr(U.component(vector::X));
        v2_ = Foam::sqr(U.component(vector::Y));
        w2_ = Foam::sqr(U.component(vector::Z));
        uPrime2_ = Foam::sqr(UPrime_.component(vector::X));
        vPrime2_ = Foam::sqr(UPrime_.component(vector::Y));
        wPrime2_ = Foam::sqr(UPrime_.component(vector::Z));

	//// Qsas calculation
	volScalarField y_ = wallDist(mesh).y();
	volScalarField delta_ = 1.0*pow(V, 1.0/3.0); //deltaCoeff=1.0
                
	volScalarField CDkOmegaPlus_ = max(((2.0*0.856)*(fvc::grad(kMod_) & fvc::grad(omegaMod_))/omegaMod_), dimensionedScalar(dimless/sqr(dimTime), 1.0e-10));
	
	volTensorField tgradU_  = fvc::grad(U);
	volScalarField S2_      = 2*magSqr(symm(tgradU_));
	volScalarField F1_      = tanh(pow4(min
                                              (
                                               min
                                                  (
                                                   max
                                                      (
                                                          (1.0/0.09)*sqrt(kMod_)/(omegaMod_*y_),
                                                          500.0*(nu)/(sqr(y_)*omegaMod_)
                                                      ),
                                                      (4.0*0.856)*kMod_/(CDkOmegaPlus_*sqr(y_))
                                                  ),
                                               scalar(10.0)
					      )
					    )
		                      );
			
	volScalarField beta_   = 0.075*F1_ + 0.0828*(1.0-F1_);
	volScalarField gamma_  = 5.0/9.0*F1_ + 0.44*(1.0-F1_);

	volScalarField L_   = sqrt(kMod_)/(pow025(0.09)*omegaMod_);
      	volScalarField Lvk_ = max
                              (
                                 0.41*Foam::sqrt(S2_)/(Foam::mag(fvc::laplacian(U)) + dimensionedScalar(dimensionSet(0, -1, -1, 0, 0, 0, 0), rootVSmall)),
                                 0.11*Foam::sqrt(0.41*3.51/(beta_/0.09 - gamma_))*delta_
                              );

        Qsas_ = min
        (
            max
            (
                3.51*0.41*S2_*Foam::sqr(L_/Lvk_)
              - (2.0*2.0/3.0*2.0)*kMod_
               *max
                (
                    Foam::magSqr(fvc::grad(omegaMod_))/Foam::sqr(omegaMod_),
                    Foam::magSqr(fvc::grad(kMod_))/Foam::sqr(kMod_)
                ),
                dimensionedScalar(dimensionSet(0, 0, -2, 0, 0, 0, 0), 0)
            ),
            // Limit SAS production of omega for numerical stability,
            // particularly during start-up
            omegaMod_/(0.1*CurrentDeltaT)*dimensionedScalar(dimensionSet(0, 0, -1, 0, 0, 0, 0), 1.0)
        );
         
	Qsas_Mean_ = ((CurrentTime-CurrentDeltaT)/CurrentTime)*Qsas_Mean_ + (CurrentDeltaT/CurrentTime)*Qsas_;
          
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
