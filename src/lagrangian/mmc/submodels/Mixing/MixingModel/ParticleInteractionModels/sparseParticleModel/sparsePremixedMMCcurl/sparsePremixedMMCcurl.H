/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenCFD Ltd.
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
    Foam::sparsePremixedMMCcurl

Description
    Sparse second-conditioning MMC mixing model for premixed flames.

    Implements the two-stage sequential conditioning method of
    Sundaram & Klimenko (PCI 36, 2017, 1937-1945) to reduce the required
    particle count per LES cell in sparse MMC.

    Reference variable space: [sPx, sPy, sPz, wOU]
        - sPx, sPy, sPz : shadow positions (primary conditioning)
        - wOU           : bounded Ornstein-Uhlenbeck variable

    Two sequential mixing stages per time step:

      Stage 1 — Shadow-position pairing
        Particles are paired in [position, sPx, sPy, sPz] space.
        Only the progress variable c is mixed between paired particles.

      Stage 2 — Modified progress variable pairing
        The modified progress variable is computed for each particle:
            phi* = c * exp(beta * wOU)
        Particles are re-paired in [position, phi*] space.
        All reactive scalars (all XiC except c) are mixed.

    The LPF (laminar premixed flame) mixing time scale is used in both stages.

    Coefficients (sparsePremixedMMCcurlCoeffs):
        CL              — mixing time scale coefficient (default 0.1)
        beta            — modulation coefficient for phi* (default 1.0)
        couplingVarName — name of the progress variable in mmcVariablesDefinitions
                          (default "c")
        wOUVarName      — name of the bounded OU reference variable
                          (default "wOU")
        numShadowPos    — number of shadow-position reference variables
                          (default 3)
        limitkdTree     — disable mixing when dXiR > 2*Xii (default true)
        meanTimeScale   — use harmonic mean (true) or min (false) of pair
                          time scales (default true)

    Reference:
        Sundaram & Klimenko, Proc. Combust. Inst. 36 (2017) 1937-1945.
        DOI: 10.1016/j.proci.2016.07.116

\*---------------------------------------------------------------------------*/

#ifndef sparsePremixedMMCcurl_H
#define sparsePremixedMMCcurl_H

#include "mixParticleModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class sparsePremixedMMCcurl Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class sparsePremixedMMCcurl
:
    public mixParticleModel<CloudType>
{
    // Typedef for particle type
    typedef typename CloudType::particleType particleType;
    typedef typename mixParticleModel<CloudType>::eulerianFieldData eulerianFieldData;

    // Private data

        //- CL coefficient for LPF mixing time scale
        const scalar CL_;

        //- Modulation coefficient β for phi* = c * exp(β * wOU)
        const scalar beta_;

        //- Name of the progress variable in XiC
        const word couplingVarName_;

        //- Name of the bounded OU reference variable in XiR
        const word wOUVarName_;

        //- Number of shadow-position reference variables (typically 3)
        const label numShadowPos_;

        //- Limit mixing of particles far in composition space
        const Switch limitkdTree_;

        //- Use harmonic mean (true) or minimum (false) of pair time scales
        const Switch meanTimeScale_;

        //- Index of c in particle XiC array (set at construction)
        label cIndex_;

        //- Index of wOU in particle XiR array (set at construction)
        label wOUIndex_;

    // Private Member Functions

        //- Apply Stage 1 mixing: mix only c for a list of pairs
        void smixStage1
        (
            const DynamicList<List<label>>& pairs,
            scalar deltaT
        );

        //- Apply Stage 2 mixing: mix all scalars except c for a list of pairs
        void smixStage2
        (
            const DynamicList<List<label>>& pairs,
            scalar deltaT
        );

        //- Compute mixing for a single pair (internal helper)
        //  stageFlag: 1 = mix only cIndex_, 2 = mix all except cIndex_
        void mixpairStaged
        (
            particleType& p,
            const eulerianFieldData& pEulFields,
            particleType& q,
            const eulerianFieldData& qEulFields,
            scalar deltaT,
            label stageFlag
        );

    public:

        //- Runtime type information
        TypeName("sparsePremixedMMCcurl");

    // Constructors

        //- Construct from dictionary
        sparsePremixedMMCcurl
        (
            const dictionary& dict,
            CloudType& owner,
            const mmcVarSet& Xi
        );

        //- Construct copy
        sparsePremixedMMCcurl(const sparsePremixedMMCcurl<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<CloudMixingModel<CloudType>> clone() const
        {
            return autoPtr<CloudMixingModel<CloudType>>
            (
                new sparsePremixedMMCcurl<CloudType>(*this)
            );
        }

    // Destructor
        virtual ~sparsePremixedMMCcurl() = default;


    // Member Functions

        //- Two-stage Smix: Stage 1 (shadow position) then Stage 2 (phi*)
        virtual void Smix();

        //- Required by base class (not used directly — Smix calls staged helpers)
        virtual void mixpair
        (
            particleType& p,
            const eulerianFieldData& pEulerianFields,
            particleType& q,
            const eulerianFieldData& qEulerianFields,
            scalar& deltaT
        );

        virtual const scalarField XiR0(label patch, label patchFace);
        virtual const scalarField XiR0(label celli);

    // Info

        void printInfo();
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "sparsePremixedMMCcurl.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
