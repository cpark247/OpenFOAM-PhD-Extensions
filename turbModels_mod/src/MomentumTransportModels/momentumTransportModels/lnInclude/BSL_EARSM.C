/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "BSL_EARSM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// Anmerkung: I == delta_ij = Einheitsmatrix == Kronecker-Delta

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
/* = spezielle Initialisierungs-Fkt., die aufgerufen wird um ein neues Objekt der Klasse zu erstellen; 
	kann auch zur Initialisierung von Attributen eines Objkts genutzt werden
*/

const dimensionedScalar smallOmega_ ("smallOmega",dimensionSet (0,0,-1,0,0,0,0), 1e-9); //Skalar mit Einheit 1/s
const dimensionedScalar smallDimless_ ("smallDimless",dimless, 1e-9); //dimensionsloses Skalar

//Vorlage, die alle Turbulenzmodelle umfasst
template<class BasicMomentumTransportModel> 
BSL_EARSM<BasicMomentumTransportModel>::BSL_EARSM
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    // const word& propertiesName,
    const word& type
)
:
    Foam::kOmegaSST // die folgenden Templates sind vom kOmegaSST-Template unmittelbar abhängig
    <
        eddyViscosity<RASModel<BasicMomentumTransportModel>>,
        BasicMomentumTransportModel
    >
	// erstelle die folgenden Variablen/Methoden aus den obigen Komponenten
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
        // propertiesName
    ),
    alphaK1_ // sigma_k1 und sigma_omega1
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_ // sigma_omega2
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856 //0.85616
        )
    ),
    // gamma1_ 
    // (
    //     dimensioned<scalar>::lookupOrAddToDict
    //     (
    //         "gamma1",
    //         this->coeffDict_,
    //         0.5532
    //     )
    // ),
    // gamma2_
    // (
    //     dimensioned<scalar>::lookupOrAddToDict
    //     (
    //         "gamma2",
    //         this->coeffDict_,
    //         0.4403
    //     )
    // ),
    // Prt_
    // (
    //     dimensioned<scalar>::lookupOrAddToDict
    //     (
    //         "Prt",
    //         this->coeffDict_,
    //         1.0
    //     )
    // ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.8
        )
    ),

    A1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A1",
            this->coeffDict_,
            1.245
        )
    ),

	// Zeitskala mit Kolmogorov-Limiter (Gl. 5) mit Cmu = betaStar
    tau_ ( max((1./(this->betaStar_*this->omega_)),6.*sqrt( this->nu() / (this->betaStar_*this->k_*this->omega_)) ) ),
	// // Trace-less vel-gradient
    gradU_( fvc::grad(this->U_)-(1./3.)*I*tr(fvc::grad(this->U_)) ), //  ---> why?? Ans: Isotropic part less important in incomp. flow
    // gradU_( fvc::grad(this->U_)),
	// dimensionsloser Spannungstensor der mittleren Strömung (Gl. 4_1)
    S_( tau_*symm(gradU_.T())),
    // S_( tau_*symm(gradU_)),
	// dimensionsloser Drehgeschwindigkeitstensor der mittleren Strömung (Gl. 4_2) 
    W_( tau_*skew(gradU_.T())),
    // W_( tau_*skew(gradU_)),
	// Tensor-Invariante (Bezug: S_, W_) (Gl. 6_1)
    IIs_( tr(S_&S_)), // why trace of a scalar or is S&S not a scalar?
    
	// Tensor-Invariante (Bezug: S_, W_) (Gl. 6_2)
    IIo_( tr(W_&W_)),
	// Tensor-Invariante (Bezug: S_, W_) (Gl. 6_3)
    IV_( tr((S_&W_)&W_)),
	// Konstante (Bezug: Gl. 9, S. 3 unten)
    C1prime_( (9./4.)*(C1_-1.)), 
	// Parameter aus (Gl. 11) definiert in (Gl. 12_1)
    P1_( C1prime_*((sqr(C1prime_)/27.)+(9.*IIs_/20.)-(2.*IIo_/3.))),
	// Parameter aus (Gl. 11) definiert in (Gl. 12_2)
    P2_( sqr(P1_)-pow(((sqr(C1prime_)/9.)+(9.*IIs_/10.)+(2.*IIo_/3.)),3)),

	// Objekt N wird definiert und durch die Gleichung in der letzten Zeile initialisiert 
	// (Gl. 9 mit Pk/epsilon = sqrt(2.*this->betaStar_*IIs_) aus Gl. 16 --> simplifiziertes BSL-EARSM-Modell)
    N_(
        IOobject
        (
            "N",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        C1prime_+(9./4.)*sqrt(2.*this->betaStar_*IIs_) // BSL-EARSM
    ),

//    N_(scalar(0.)*IIs_),

	// Parameter aus (Gl. 7) definiert in (Gl. 8_1)
    Q_( (sqr(N_)-2.*IIo_)/A1_),
	// Parameter aus (Gl. 7) definiert in (Gl. 8_2)
    Q1_( Q_*(2.*sqr(N_)-IIo_)/6.),
	// Koeffizient der Tensorbasis aus (Gl. 2) definiert in (Gl. 7_1)
    Beta1_( -1.*N_/(Q_+smallDimless_)),
	// Koeffizient der Tensorbasis aus (Gl. 2) definiert in (Gl. 7_2) mit beta2 = 0
    Beta2_( 0.*N_),         //0 work around
	// Koeffizient der Tensorbasis aus (Gl. 2) definiert in (Gl. 7_3)
    Beta3_( -2.*IV_/(N_*Q1_+smallDimless_)),
	// Koeffizient der Tensorbasis aus (Gl. 2) definiert in (Gl. 7_4)
    Beta4_( -1.*N_/(Q_+smallDimless_)),
	// Koeffizient der Tensorbasis aus (Gl. 2) definiert in (Gl. 7_5)
    Beta6_( -1.*N_/(Q1_+smallDimless_)),
	// Koeffizient der Tensorbasis aus (Gl. 2) definiert in (Gl. 7_1)
    Beta9_( 1./(Q1_+smallDimless_)),

	// Term T1 definiert in (Gl. 3_11)
    check0_(tr(S_)),        //T1 term
	// Term T2 definiert in (Gl. 3_12) 
    check1_(tr(((S_ & S_)-(1./3.)*IIs_*I))),        //T2 term
	// Term T3 definiert in (Gl. 3_13)
    check2_(tr(((W_ & W_)-(1./3.)*IIo_*I))),        //T3 term
	// Term T4 definiert in (Gl. 3_14)
    check3_(tr(((S_ & W_)-(W_ & S_)))),             //T4 term
	// Term T6 definiert in (Gl. 3_21)
    check4_(tr((((S_ & W_) & W_) + ((W_ & W_) & S_) - (2./3.)*IV_*I - IIo_*S_))),   //T6 term
	// Term T9 definiert in (Gl. 3_31)
    check5_(tr((((W_ & S_) & (W_ & W_)) - ((W_ & W_) & (S_ & W_)) + 0.5*IIo_*((S_ & W_) - (W_ & S_))))),  //T9 term

	// Objekt des Anisotropietensor wird definiert und durch (Gl. 2) initialisiert
    aij_(
        IOobject
        (
            "aij",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
            symm(
            Beta1_*S_									// beta1 * T1
											// beta2=0 * T2 == 0
            + Beta3_*((W_ & W_)-(1./3.)*IIo_*I)						// beta3 * T3
            + Beta4_*((S_ & W_)-(W_ & S_))						// beta4 * T4
            + Beta6_*(((S_ & W_) & W_) + ((W_ & W_) & S_) - (2./3.)*IV_*I - IIo_*S_)	// beta6 * T6
											// beta9=0 * T9 == 0 -> no visible alteration [Menter 2012]
            )
        ),

	// Reynolds-Spannungs-Tensor (Gl. 1)
    tauij_(
        "tauij",
        symm(
            this->k_*(aij_+(2./3.)*I)
            )
        ),
	// Non-linear part of Reynolds Stresses    
    Nij_(
        "Nij",
        symm(

            (this->k_*aij_+2.0*(this->nut_)*symm(gradU_))*this->rho_

            )
        )

{
    if (type == typeName)
    {
        this->printCoeffs(type);
        Info << "BSL-EARSM Turbulence Model." << endl;
    }
}

// Reynoldsspannungen in symmetrischem Tensor, oberes rechtes Dreieck besetzt, 6 statt 9 Einträge
template<class BasicMomentumTransportModel>
tmp<volSymmTensorField> BSL_EARSM<BasicMomentumTransportModel>::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*this->k_ - (this->nut_)*dev(twoSymm(fvc::grad(this->U_)))+Nij_/this->rho_,
            this->k_.boundaryField().types()
        )
    );
}

// Methode, die durch Template "CViscousStress" bereitgestellt wird
template<class BasicMomentumTransportModel>
tmp<volSymmTensorField> BSL_EARSM<BasicMomentumTransportModel>::devRhoReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),

           dev(Nij_) - this->rho_*this->nuEff()*dev(twoSymm(fvc::grad(this->U_)))
        )
    );
}

template<class BasicMomentumTransportModel>
tmp<fvVectorMatrix> BSL_EARSM<BasicMomentumTransportModel>::divDevRhoReff(volVectorField& U) const
{
    //Info << "Using divDevRhoReff from BSL-EARSM" << endl;
    return
    (
        fvc::div(this->alpha_*Nij_)
 	  - fvc::div((this->alpha_*this->rho_*this->nuEff())*dev2(T(fvc::grad(U))))
 	  - fvm::laplacian(this->alpha_*this->rho_*this->nuEff(), U)
    );
}

template<class BasicMomentumTransportModel>
void BSL_EARSM<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        this->nut_ = this->k_/(this->omega_+smallOmega_);
        this->nut_.correctBoundaryConditions();
        this->correctNut();
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    // const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    BasicMomentumTransportModel::correct();

	// skalares Feld mit Werten im Zellmittelpunkt
    volScalarField::Internal divU(fvc::div(fvc::absolute(this->phi(), U)));

	// Tensorfeld tgradU wird definiert und mit den Werten des Gradienten des Geschwindigkeitsfeldes U initialisiert
    tmp<volTensorField> tgradU = fvc::grad(U);

	// Zeitskala mit Kolmogorov-Limiter (Gl. 5) 
    tau_ = max((1./(this->betaStar_*this->omega_)),6.*sqrt( this->nu()  /(this->betaStar_*this->k_*this->omega_)));
	// 
    // gradU_ = fvc::grad(this->U_)-(1./3.)*I*tr(fvc::grad(this->U_));
	// dimensionsloser Spannungstensor der mittleren Strömung (Gl. 4_1)
    S_ = tau_*symm(gradU_.T());
	// dimensionsloser Drehgeschwindigkeitstensor der mittleren Strömung (Gl. 4_2)
    W_ = tau_*skew(gradU_.T());

	// Tensor-Invariante (Bezug: S_, W_) (Gl. 6_1)
    // Info << "S_&S_ = " << S_&S_ << endl;
    IIs_ = tr(S_ & S_);
	// Tensor-Invariante (Bezug: S_, W_) (Gl. 6_2)
    IIo_ = tr(W_ & W_);
	// Tensor-Invariante (Bezug: S_, W_) (Gl. 6_3)
    IV_  = tr((S_ & W_) & W_);
	// Konstante (Bezug: Gl. 9, S. 3 unten)
    C1prime_ = (9./4.)*(C1_-1.);           //constant
	// Parameter aus (Gl. 11) definiert in (Gl. 12_1)
    P1_ = C1prime_*((sqr(C1prime_)/27.)+(9.*IIs_/20.)-(2.*IIo_/3.));
	// Parameter aus (Gl. 11) definiert in (Gl. 12_2)
    P2_ = sqr(P1_)-pow(((sqr(C1prime_)/9.)+(9.*IIs_/10.)+(2.*IIo_/3.)),3);

	// Berechnung des Parameters N (Bezug Gl. 10 -> 9 -> 8 -> 7) definiert durch (Gl. 11)
    forAll(N_,cellI)
    {
	// Falls P2 >= 0, dann: (Gl. 11_1)
        if(P2_[cellI]  >= 0)
        {
        N_[cellI]= (C1prime_.value()/3.)+pow((P1_[cellI]+sqrt(P2_[cellI])),(1./3.))
                    +sign(P1_[cellI]-sqrt(P2_[cellI]))*pow(mag(P1_[cellI]-sqrt(P2_[cellI])),(1./3.));
        }
	// Falls P2 > 0, dann: (Gl. 11_2)
        else
        {
 
          N_[cellI]= (C1prime_.value()/3.)+2.*pow((sqr(P1_[cellI])-P2_[cellI]),(1./6.))
                    *cos((1./3.)*acos(P1_[cellI]/(sqrt(sqr(P1_[cellI])-P2_[cellI]))));
        }
    }

    //N_.correctBoundaryConditions(); 
	// Parameter aus (Gl. 7) definiert in (Gl. 8_1)
    Q_ = (sqr(N_)-2.*IIo_)/A1_;
	// Parameter aus (Gl. 7) definiert in (Gl. 8_2)
    Q1_ = Q_*(2.*sqr(N_)-IIo_)/6.;
	// Koeffizient der Tensorbasis aus (Gl. 2) definiert in (Gl. 7_1)
    Beta1_ = -1.*N_/(Q_+smallDimless_);
	// Koeffizient der Tensorbasis aus (Gl. 2) definiert in (Gl. 7_2) mit beta2 = 0
    Beta2_ = 0.*N_;     //0 work around
	// Koeffizient der Tensorbasis aus (Gl. 2) definiert in (Gl. 7_3)
    Beta3_ = -2.*IV_/(N_*Q1_+smallDimless_);
	// Koeffizient der Tensorbasis aus (Gl. 2) definiert in (Gl. 7_4)
    Beta4_ = -1.*N_/(Q_+smallDimless_);
	// Koeffizient der Tensorbasis aus (Gl. 2) definiert in (Gl. 7_5)
    Beta6_ = -1.*N_/(Q1_+smallDimless_);
	// Koeffizient der Tensorbasis aus (Gl. 2) definiert in (Gl. 7_1)
    Beta9_ = 1./(Q1_+smallDimless_);

	// Term T1 definiert in (Gl. 3_11)
    check0_=tr(S_);				//T1 term
	// Term T2 definiert in (Gl. 3_12)
    check1_=tr(((S_ & S_)-(1./3.)*IIs_*I));	//T2 term
	// Term T3 definiert in (Gl. 3_13)
    check2_=tr(((W_ & W_)-(1./3.)*IIo_*I));	//T3 term
	// Term T4 definiert in (Gl. 3_14)
    check3_=tr(((S_ & W_)-(W_ & S_)));		//T4 term
	// Term T6 definiert in (Gl. 3_21)
    check4_=tr((((S_ & W_) & W_) + ((W_ & W_) & S_) - (2./3.)*IV_*I - IIo_*S_));	//T6 term
	// Term T9 definiert in (Gl. 3_31)
    check5_=tr((((W_ & S_) & (W_ & W_)) - ((W_ & W_) & (S_ & W_)) + 0.5*IIo_*((S_ & W_) - (W_ & S_))));	//T9 term


    //Info<< "div(U_)  "<<" max:"<< max(fvc::div(this->U_)).value() <<
    //" min:"<< min(fvc::div(this->U_)).value() <<" avg:"<< average(fvc::div(this->U_)).value() << endl;
    Info<< "div(phi_)     "<<" max:"<< max(fvc::div(this->phi())).value() <<
    " min:"<< min(fvc::div(this->phi())).value() <<" avg:"<< average(fvc::div(this->phi())).value() << endl;
    Info<< "gradU_ TRACE  "<<" max:"<< max(tr(gradU_)).value() <<
    " min:"<< min(tr(gradU_)).value() <<" avg:"<< average(tr(gradU_)).value() << endl;
    Info<< "T1 TRACE      "<<" max:"<< max(check0_).value() <<
    " min:"<< min(check0_).value() <<" avg:"<< average(check0_).value() << endl;
    Info<< "T3 TRACE      "<<" max:"<< max(check2_).value() <<
    " min:"<< min(check2_).value() <<" avg:"<< average(check2_).value() << endl;
    Info<< "T4 TRACE      "<<" max:"<< max(check3_).value() <<
    " min:"<< min(check3_).value() <<" avg:"<< average(check3_).value() << endl;
    Info<< "T6 TRACE      "<<" max:"<< max(check4_).value() <<
    " min:"<< min(check4_).value() <<" avg:"<< average(check4_).value() << endl;
    Info<< "N             "<<" max:"<< max(N_).value() <<
    " min:"<< min(N_).value() <<" avg:"<< average(N_).value() << endl;

	// Objekt des Anisotropietensor wird definiert und durch (Gl. 2) initialisiert
    aij_ = symm(
            Beta1_*S_
            //+ Beta2_*((S_ & S_)-(1./3.)*IIs_*I)
            + Beta3_*((W_ & W_)-(1./3.)*IIo_*I)
            + Beta4_*((S_ & W_)-(W_ & S_))
            + Beta6_*(((S_ & W_) & W_) + ((W_ & W_) & S_) - (2./3.)*IV_*I - IIo_*S_)
            //+ Beta9_*(((W_ & S_) & (W_ & W_)) - ((W_ & W_) & (S_ & W_)) + 0.5*IIo_*((S_ & W_) - (W_ & S_)))
            );
    Info<< "aij_ TRACE    "<<" max:"<< max(tr(aij_)).value() <<
    " min:"<< min(tr(aij_)).value() <<" avg:"<< average(tr(aij_)).value() << endl;

	// Reynolds-Süannungs-Tensor (Gl. 1)
    tauij_ = symm(this->k_*(aij_+(2./3.)*I));
	//
    Nij_ = symm(this->k_*aij_+2.0*this->nut_*symm(gradU_))*this->rho_;       //Non-linear part of Reynolds Stresses    


	// skalares Feld mit Werten im Zellmittelpunkt 
		// GName = Hilfs-Fkt. zur Wiedergabe des Namens des turbulenten G-Feldes
    volScalarField::Internal G(this->GName(), -(tauij_&&gradU_));
    tgradU.clear();

    // Update omega and G at the wall
    this->omega_.boundaryFieldRef().updateCoeffs();

	// skalares Feld mit Werten im Zellmittelpunkt definiert durch (Gl. C-10)
		// alphaOmega2 = Sigma_(Omega 2) in Gl. 
    volScalarField CDkOmega
    (
        (2*this->alphaOmega2_)*(fvc::grad(this->k_) & fvc::grad(this->omega_))/this->omega_
    );

	// skalares Feld mit Werten im Zellmittelpunkt
		// F1 = Blending-Fkt. für sigma_k, definiert in kOmegaSSTBase.C ab Zeile 40
		// F1 als Blending-Fkt. mit F1 = 0 im Freistrom für k-epsilon und mit F1 = 1 in GS für k-omega
    volScalarField F1(this->F1(CDkOmega));
    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));

        // Turbulent frequency equation (Gl. 13_2)
        tmp<fvScalarMatrix> omegaEqn
        (
		// zeitliche Ableitung = D omega / Dt 			--> Zeitableitungs-Term
            fvm::ddt(alpha, rho, this->omega_)
		// Divergenz						
          + fvm::div(alphaRhoPhi, this->omega_)
		// Laplace, DomegaEff = effektive Diffusivität von omega --> Laplace-Term
          - fvm::laplacian(alpha*rho*this->DomegaEff(F1), this->omega_)
         ==
		// = (gamma * omaga / k) * P_k 				--> P_k-Term
			// mit P_k = min(G, ...) (Gl. 14) 
            alpha()*rho()*gamma
           *min
            (
                G,
                this->c1_*this->betaStar_*this->k_()*this->omega_()
            )*this->omega_()/this->k_()
		// = Quellterm = 2/3 * gamma * divU * omega		--> ????????????????
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, this->omega_)
		// = Quellterm = - beta * omega²			--> beta - omega²-Term
          - fvm::Sp(alpha()*rho()*beta*this->omega_(), this->omega_)
		// = Quellterm = - (F1 - 1) * CDkOmega / omega * omega	--> CDkOmega-Term
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/this->omega_(),
                this->omega_
            )
        );

        omegaEqn.ref().relax();
        fvConstraints.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(this->omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvConstraints.constrain(this->omega_);
        bound(this->omega_, this->omegaMin_);
    }

    // Turbulent kinetic energy equation (Gl. 13_1)
    tmp<fvScalarMatrix> kEqn
    (
	// zeitliche Ableitung = Dk / Dt				--> Zeitableitungs-Term
        fvm::ddt(alpha, rho, this->k_)
	// Divergenz							
      + fvm::div(alphaRhoPhi, this->k_)
	// Laplace, DkEff = effektive Diffusivität von k 		--> Laplace-Term
      - fvm::laplacian(alpha*rho*this->DkEff(F1), this->k_)
     ==
	// = P_k = min(G, ...) 						--> P_k-Term
        min
	(
		alpha()*rho()*G, 
		(this->c1_*this->betaStar_)*alpha()*rho()*this->k_()*this->omega_()
	)
	// = Quellterm = - 2/3 * divU * k				--> ???????????????
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, this->k_)
	// = Quellterm = - betaStar * omega * k				--> beta*-k-Omega-Term
      - fvm::Sp(alpha()*rho()*this->betaStar_*this->omega_(),this->k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(this->k_);
    bound(this->k_, this->kMin_);

    //this->correctNut(S2, F23);
    nut = this->k_/(this->omega_+smallOmega_);
    nut.correctBoundaryConditions();
    this->correctNut();
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

//
