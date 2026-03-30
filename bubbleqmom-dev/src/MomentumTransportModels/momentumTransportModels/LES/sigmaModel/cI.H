const volScalarField mm
(
    sqr(this->delta())*(4*sqr(mag(filter_(dev(symm(D))))) - filter_(sqr(mag(dev(symm(D))))))
);

dimensionedScalar mmmm = average(magSqr(mm));

if (mmmm.value() > VSMALL)
{
    tmp<volScalarField> KK =
        0.5*(filter_(magSqr(this->U_)) - magSqr(filter_(this->U_)));

    return average(KK*mm)/mmmm;
}
else
{
  tmp<volScalarField> KK =
      0.5*(filter_(magSqr(this->U_)) - magSqr(filter_(this->U_)));

  return 0.0;
}
