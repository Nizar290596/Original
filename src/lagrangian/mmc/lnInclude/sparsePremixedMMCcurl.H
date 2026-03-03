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

    Progress variable
    -----------------
    c is computed every time step from each particle's temperature using

        c = (T_p - Tu) / (Tb - Tu)                              (1)

    clamped to [0, 1], where T_p is the particle temperature (p.T()),
    Tu is the unburned temperature, and Tb is the burned temperature.
    The result is stored in the particle's XiC array (XiC[cIndex_]).
    For remote particles in parallel runs, T is interpolated from the
    Eulerian "T" field at the particle's position.

    Two sequential mixing stages per time step
    ------------------------------------------
      Stage 1 — Shadow-position pairing  (full particle set)
        ALL particles are paired in [position, sPx, sPy, sPz] space.
        Only the progress variable c is mixed.

      Stage 2 — Modified progress variable pairing  (flame-zone subset)
        Only particles with c in [cMin, cMax] (default [0.05, 0.95])
        are included. These are the flame-zone particles; fully unburned
        or fully burned particles are excluded. For each selected particle:
            phi* = c * exp(beta * wOU)
        The subset is paired in [position, phi*] space and all reactive
        scalars (all XiC except c) are mixed.

        This flame-zone filter is the sparsity mechanism described in
        Sundaram & Klimenko 2017: the reactive-scalar mixing only needs
        the fraction of particles inside the flame, dramatically reducing
        the number required for Stage 2 (e.g. 2500 from 10000).

    Coefficients (sparsePremixedMMCcurlCoeffs)
    -------------------------------------------
        Tu              — unburned temperature [K]  (default 298 K)
        Tb              — burned temperature [K]    (default 2240 K)
        cMin            — lower c limit for Stage 2 subset (default 0.05)
        cMax            — upper c limit for Stage 2 subset (default 0.95)
        CL              — LPF mixing time scale coefficient (default 0.1)
        beta            — modulation coefficient for phi* (default 1.0)
        couplingVarName — name of the c entry in mmcVariablesDefinitions
                          (default "c")
        wOUVarName      — name of the bounded OU reference variable
                          (default "wOU")
        numShadowPos    — number of shadow-position reference variables
                          (default 3: sPx, sPy, sPz)
        limitkdTree     — disable mixing when dXiR > 2*Xii (default true)
        meanTimeScale   — harmonic mean (true) or min (false) of pair time
                          scales (default true)

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
    typedef typename CloudType::particleType particleType;
    typedef typename mixParticleModel<CloudType>::eulerianFieldData eulerianFieldData;

    // Private data

        //- Unburned temperature Tu [K]  (denominator lower bound in eq. 1)
        const scalar Tu_;

        //- Burned temperature Tb [K]  (denominator upper bound in eq. 1)
        const scalar Tb_;

        //- Lower c bound for Stage 2 subset (flame zone filter)
        const scalar cMin_;

        //- Upper c bound for Stage 2 subset (flame zone filter)
        const scalar cMax_;

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

        //- Use harmonic mean (true) or min (false) of pair time scales
        const Switch meanTimeScale_;

        //- Index of c in XiC array (set at construction)
        label cIndex_;

        //- Index of wOU in XiR array (set at construction)
        label wOUIndex_;

    // Private Member Functions

        //- Compute progress variable from temperature (eq. 1), clamped [0,1]
        inline scalar progressVar(scalar T) const
        {
            return Foam::max(0.0, Foam::min(1.0, (T - Tu_) / (Tb_ - Tu_)));
        }

        //- Stage 1: mix only c for a list of pairs.
        //  Pair indices resolve into eulerianFieldDataList_.
        void smixStage1
        (
            const DynamicList<List<label>>& pairs,
            scalar deltaT
        );

        //- Stage 2: mix all scalars except c for a list of pairs.
        //  Pair indices resolve into subsetList, NOT eulerianFieldDataList_.
        void smixStage2
        (
            const DynamicList<eulerianFieldData>& subsetList,
            const DynamicList<List<label>>& pairs,
            scalar deltaT
        );

        //- Perform staged Curl mixing for one pair.
        //  stageFlag 1 → mix only cIndex_; 2 → mix everything except cIndex_.
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

        TypeName("sparsePremixedMMCcurl");

    // Constructors

        sparsePremixedMMCcurl
        (
            const dictionary& dict,
            CloudType& owner,
            const mmcVarSet& Xi
        );

        sparsePremixedMMCcurl(const sparsePremixedMMCcurl<CloudType>& cm);

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

        //- Two-stage Smix.
        //  Computes c = (T_p - Tu)/(Tb - Tu) for all local particles, then
        //  runs Stage 1 (full set, shadow-position pairing, mix c) and
        //  Stage 2 (flame-zone subset, phi* pairing, mix reactive scalars).
        virtual void Smix();

        //- Required by base class interface (fallback: calls Stage 2 mixing)
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

        void printInfo();
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "sparsePremixedMMCcurl.C"
#endif

#endif

// ************************************************************************* //
