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
    Foam::secondCondMMCcurl

Description
    Second-conditioning mixing model for the MMC framework.

    Operates exclusively on the particle subset flagged with
    secondCondFlag == 1 and pairs particles in the 4-dimensional reference
    space

        (φ°, ξ_x, ξ_y, ξ_z)

    where
      - φ° = phiModified (modified progress variable, already computed from
                          the OU process)
      - ξ = particle position relative to a user-supplied shadowOrigin

    The k-d tree splitting is driven by normalisation factors read from the
    coefficient dictionary.  A tight normalisation for φ° (e.g. phiMod_m =
    0.01) causes the tree to split preferentially on the progress-variable
    axis; looser normalisations for the shadow-position coordinates (e.g.
    sPx_m = sPy_m = sPz_m = 0.05) provide spatial context with lower
    priority.  Physical coordinates are NOT included as independent splitting
    dimensions (consistent with the first-conditioning MMCcurl model).

    Note
      - The mixing rule applied to found pairs is deferred (mixpair is a
        stub) and will be added in a subsequent step.
      - Parallel pairing is not yet implemented; local pairing only.

    Usage in cloudProperties:
    \verbatim
    mixingModel    secondCondMMCcurl;

    secondCondMMCcurlCoeffs
    {
        // Required by base class (first-conditioning XiR normalisation)
        r_i     1e-3;
        Xim_i { fm_m 0.03; }

        pairingMethod   local;

        // Second-conditioning normalisation factors
        phiMod_m    0.01;   // tight  -> dominates k-d tree splitting
        sPx_m       0.05;
        sPy_m       0.05;
        sPz_m       0.05;

        // Names of the shadow-position XiR variables as defined in
        // cloudProperties referenceVariables (must match exactly)
        shadowPosNames  (sPx sPy sPz);
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef secondCondMMCcurl_H
#define secondCondMMCcurl_H

#include "mixParticleModel.H"
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class secondCondMMCcurl Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class secondCondMMCcurl
:
    public mixParticleModel<CloudType>
{
    // Typedefs (mirror MMCcurl conventions)
    typedef typename CloudType::particleType particleType;
    typedef typename mixParticleModel<CloudType>::eulerianFieldData
        eulerianFieldData;

    // Private data

        //- Normalisation factors for the 4D reference space
        //  XiR order: [0] phiModified, [1] ξ_x, [2] ξ_y, [3] ξ_z
        scalar phiMod_m_;
        scalar sPx_m_;
        scalar sPy_m_;
        scalar sPz_m_;

        //- Names of the three shadow-position XiR variables as declared in
        //  cloudProperties referenceVariables (e.g. sPx, sPy, sPz).
        //  Values are read from the particle's XiR field at each timestep
        //  using the same indexInXiR_ lookup used by the base solver.
        word sPxName_;
        word sPyName_;
        word sPzName_;


    // Private member functions

        //- Build the local particle list and find pairs.
        //  Only particles with secondCondFlag == 1 are included.
        //  Populates particleList_, eulerianFieldDataList_, particlePairs_.
        void buildSecondCondList();

        //- Pair-finding wrapper that drives secondCondKdTree.
        void findSecondCondPairs
        (
            const DynamicList<eulerianFieldData>& eulerianFieldList,
            DynamicList<List<label>>& pairs
        ) const;

        //- Recursive k-d tree search in the 4D space (φ°, ξ_x, ξ_y, ξ_z).
        //
        //  All four variables are stored in eulerianFieldData::XiR()[0..3].
        //  The dimension with the largest normalised spread is selected as
        //  the split axis.  The lessArg comparator from mixParticleModel is
        //  invoked with ncond = bestI + 3 so that ii_ >= 3 → XiR()[ii_-3]
        //  is used, bypassing the physical-position branch (ii_ < 3).
        //  Physical coordinates are therefore never used as split dimensions.
        void secondCondKdTree
        (
            const DynamicList<eulerianFieldData>& particleList,
            label l,
            label u,
            std::vector<label>& pInd,
            std::vector<label>& L,
            std::vector<label>& U
        ) const;


public:

    //- Runtime type information
    TypeName("secondCondMMCcurl");


    // Constructors

        //- Construct from dictionary
        secondCondMMCcurl
        (
            const dictionary& dict,
            CloudType& owner,
            const mmcVarSet& Xi
        );

        //- Construct copy
        secondCondMMCcurl(const secondCondMMCcurl<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<CloudMixingModel<CloudType>> clone() const
        {
            return autoPtr<CloudMixingModel<CloudType>>
            (
                new secondCondMMCcurl<CloudType>(*this)
            );
        }


    // Destructor
        virtual ~secondCondMMCcurl() = default;


    // Member Functions

        //- Override Smix: build flagged-particle list, pair in 4D, then mix
        virtual void Smix() override;

        //- Pair-wise mixing for second-conditioning pairs.
        //  Stub — mixing rule to be implemented in a subsequent step.
        virtual void mixpair
        (
            particleType& p,
            const eulerianFieldData& pEulerianFields,
            particleType& q,
            const eulerianFieldData& qEulerianFields,
            scalar& deltaT
        ) override;

        //- Boundary initialisation (returns empty — not applicable here)
        virtual const scalarField XiR0
        (
            label patch,
            label patchFace,
            particle& p
        ) override;

        virtual const scalarField XiR0(label celli, particle& p) override;

    // Info

        void printInfo();
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "secondCondMMCcurl.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
