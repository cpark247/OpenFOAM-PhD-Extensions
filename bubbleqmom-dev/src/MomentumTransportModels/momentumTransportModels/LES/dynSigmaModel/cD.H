tensor G;
scalar i1;
scalar i2;
scalar i3;

scalar alpha1;
scalar alpha2;
scalar alpha3;

scalarField sigma1(D.size());
scalarField sigma2(D.size());
scalarField sigma3(D.size());

label clipped_alpha1 = 0;
scalar min_alpha1 = 1e10;
label clipped_sigma1 = 0;
scalar min_sigma1 = 1e10;
// method B in appendix nicoud_2011
forAll(D,celli)
{
  G = D[celli].T()&D[celli];

  i1 = tr(G);
  i2 = 0.5*(sqr(tr(G)) - tr(G&G));
  i3 = det(G);

  alpha1 = sqr(i1)/9.0 - i2/3.0;
    if (alpha1<0.0)
    //if (alpha1<=0.0)
    {
        min_alpha1 = min(min_alpha1, alpha1);
        alpha1 = SMALL;
        ++clipped_alpha1;
        /*
        FatalErrorIn("sigmaModel") << "We have a problem with alpha1 it is smaller than zero: " << alpha1
          << " in cell: " << celli << endl
          << abort(FatalError);
          */
    }
  alpha2 = pow(i1,3.0)/27.0 - i1*i2/6.0 + i3/2.0;
  //alpha3 = (3.*acos(max(min(alpha2/(pow(alpha1, 1.5)+0.1),1),-1))+0.1);
  alpha3 = 1./3.*acos(max(min(alpha2/(pow(alpha1, 1.5)+SMALL),1),-1));

    sigma1[celli] = sqrt(max(i1/3.0 + 2.0*sqrt(max(alpha1,0))*cos(alpha3),0));
  sigma2[celli] = sqrt(max(i1/3.0 - 2.0*sqrt(max(alpha1,0))*cos(constant::mathematical::pi/3. + alpha3),0));
    sigma3[celli] = sqrt(max(i1/3.0 - 2.0*sqrt(max(alpha1,0))*cos(constant::mathematical::pi/3. - alpha3),0));
    if (sigma1[celli]<0.0)
    //if (sigma1[celli]<=0.0)
    {
        min_sigma1 = min(min_sigma1, sigma1[celli]);
        sigma1[celli] = SMALL;
        ++clipped_sigma1;
        /*
        FatalErrorIn("sigmaModel") << "We have a problem with sigma1 it is smaller than zero: " << sigma1[celli]
          << " in cell: " << celli << endl
          << abort(FatalError);
        */
    }
}
reduce(clipped_alpha1, sumOp<label>());
reduce(min_alpha1, minOp<scalar>());
if (clipped_alpha1 > 0)
{
    Info << "Bounding alpha1 in " << clipped_alpha1 << ": min " << min_alpha1 << endl;
}

reduce(clipped_sigma1, sumOp<label>());
reduce(min_sigma1, minOp<scalar>());
if (clipped_sigma1 > 0)
{
    Info << "Bounding sigma1 in " << clipped_sigma1 << ": min " << min_sigma1 << endl;
}

//tmp<volScalarField> tmpPtr(new volScalarField(this->DSigma_));
//tmpPtr().internalField() = sigma3*(sigma1-sigma2)*(sigma2-sigma3)/(sqr(sigma1)+SMALL);
volScalarField tmp=DSigma_;
//tmp = sigma3*(sigma1-sigma2)*(sigma2-sigma3)/(sqr(sigma1)+SMALL)
forAll(DSigma_,celli)
{
	tmp[celli] = sigma3[celli]*(sigma1[celli]-sigma2[celli])*(sigma2[celli]-sigma3[celli])/(sqr(sigma1[celli])+SMALL);
}
return tmp;
