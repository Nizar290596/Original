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

Typedef
    Foam::scalarBoundedOrnsteinUhlenbeck

Description
    Scalar specialisation of the BoundedOrnsteinUhlenbeck ItoEqn.

    Implements a bounded random Ornstein-Uhlenbeck (OU) process:

        dω = -θ·ω·dt + sqrt(2θ)·dW

    where θ is the relaxation rate and dW is a Wiener increment.
    After the Euler-Maruyama step the result is clipped to [-wMax, +wMax].

    This is a zero-mean stochastic process with no Eulerian relaxation target,
    making it distinct from the ShadowPosition SDE. It is used as a second
    conditioning variable in sparse MMC following:

        Sundaram & Klimenko, Proc. Combust. Inst. 36 (2017) 1937-1945.
        DOI: 10.1016/j.proci.2016.07.116

    The modified progress variable is then computed externally as:
        phi* = phi * exp(beta * omega_OU)

SourceFiles
    scalarBoundedOrnsteinUhlenbeck.C

\*---------------------------------------------------------------------------*/

#ifndef scalarBoundedOrnsteinUhlenbeck_H
#define scalarBoundedOrnsteinUhlenbeck_H

#include "scalar.H"
#include "ItoEqn.H"
#include "scalarItoEqn.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class T>
class BoundedOrnsteinUhlenbeck;

template<>
class BoundedOrnsteinUhlenbeck<scalar>;

typedef BoundedOrnsteinUhlenbeck<scalar> scalarBoundedOrnsteinUhlenbeck;

/*---------------------------------------------------------------------------*\
         Class BoundedOrnsteinUhlenbeck<scalar> Declaration
\*---------------------------------------------------------------------------*/

template<>
class BoundedOrnsteinUhlenbeck<scalar>
:
    public scalarItoEqn
{
    protected:

        // Protected data

            //- Relaxation rate θ (1/time)
            scalar theta_;

            //- Upper bound for ω_OU (symmetric: [-wMax_, +wMax_])
            scalar wMax_;

    public:

        //- Runtime type information
        TypeName("scalarBoundedOrnsteinUhlenbeck");

        // Constructors

            //- Construct from dictionary and mesh
            BoundedOrnsteinUhlenbeck
            (
                const dictionary& entry,
                const fvMesh& mesh
            );

            //- Construct copy
            BoundedOrnsteinUhlenbeck
            (
                const BoundedOrnsteinUhlenbeck& bou
            );

        //- Destructor
        virtual ~BoundedOrnsteinUhlenbeck()
        {};

        // Member functions

            //- Drift coefficient: A = -θ·ω
            virtual void computeDriftCoeff
            (
                const particle& p,
                const scalar& omega
            ) const;

            //- Diffusion coefficient: D = sqrt(2θ)
            virtual void computeDiffusionCoeff(const particle& p) const;

            //- Euler-Maruyama step with bounding
            virtual scalar integrate
            (
                const particle& p,
                scalar dt,
                const scalar& omega
            ) const;

            //- Initial condition: zero (OU process starts at its mean)
            virtual scalar t0
            (
                const particle& p,
                scalar dt,
                const scalar& omega
            ) const;

            //- Clone
            virtual autoPtr<ItoEqn<scalar>> clone() const
            {
                return autoPtr<ItoEqn<scalar>>
                (
                    new BoundedOrnsteinUhlenbeck<scalar>(*this)
                );
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
