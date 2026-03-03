/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

Description
    Scalar specialisation of the BoundedOrnsteinUhlenbeck ItoEqn.

    Standard bounded random Ornstein-Uhlenbeck process:
        dω = -θ·ω·dt + sqrt(2θ)·dW
    clipped to [-wMax, +wMax] after each step.

    Reference:
        Sundaram & Klimenko, Proc. Combust. Inst. 36 (2017) 1937-1945.
        DOI: 10.1016/j.proci.2016.07.116

\*---------------------------------------------------------------------------*/

#include "scalarBoundedOrnsteinUhlenbeck.H"

// * * * * * * * * * * * * * Static Member Data * * * * * * * * * * * * * * //

#include "scalarBoundedOrnsteinUhlenbeckModels.C"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BoundedOrnsteinUhlenbeck<Foam::scalar>::BoundedOrnsteinUhlenbeck
(
    const dictionary& entry,
    const fvMesh& /*mesh*/
)
:
    scalarItoEqn(),

    theta_
    (
        readScalar(entry.lookup("theta"))
    ),

    wMax_
    (
        readScalar(entry.lookup("wMax"))
    )
{}


Foam::BoundedOrnsteinUhlenbeck<Foam::scalar>::BoundedOrnsteinUhlenbeck
(
    const BoundedOrnsteinUhlenbeck& bou
)
:
    scalarItoEqn(),
    theta_(bou.theta_),
    wMax_(bou.wMax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::BoundedOrnsteinUhlenbeck<Foam::scalar>::computeDriftCoeff
(
    const particle& /*p*/,
    const scalar& omega
) const
{
    // Drift: A = -θ·ω  (mean-reverting toward zero)
    this->A_ = -theta_ * omega;
}


void Foam::BoundedOrnsteinUhlenbeck<Foam::scalar>::computeDiffusionCoeff
(
    const particle& /*p*/
) const
{
    // Diffusion: D = sqrt(2θ)  (stationary variance = 1 for unit wMax)
    this->D_ = Foam::sqrt(2.0 * theta_);
}


Foam::scalar Foam::BoundedOrnsteinUhlenbeck<Foam::scalar>::integrate
(
    const particle& p,
    scalar dt,
    const scalar& omega
) const
{
    computeDriftCoeff(p, omega);
    computeDiffusionCoeff(p);

    scalar dw = this->rndGen_.Normal(0, 1) * Foam::sqrt(dt);

    // Euler-Maruyama step
    scalar omegaNew = omega + this->A_ * dt + this->D_ * dw;

    // Apply symmetric bound: clip to [-wMax_, +wMax_]
    omegaNew = Foam::max(-wMax_, Foam::min(wMax_, omegaNew));

    return omegaNew;
}


Foam::scalar Foam::BoundedOrnsteinUhlenbeck<Foam::scalar>::t0
(
    const particle& /*p*/,
    scalar /*dt*/,
    const scalar& /*omega*/
) const
{
    // OU process is initialised at its mean (zero)
    return 0.0;
}


// ************************************************************************* //
