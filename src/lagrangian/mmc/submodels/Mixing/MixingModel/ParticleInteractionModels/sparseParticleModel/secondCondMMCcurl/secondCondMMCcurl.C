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

\*---------------------------------------------------------------------------*/

#include "secondCondMMCcurl.H"
#include <numeric>
#include <algorithm>
#include <cmath>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::secondCondMMCcurl<CloudType>::secondCondMMCcurl
(
    const dictionary& dict,
    CloudType& owner,
    const mmcVarSet& Xi
)
:
    mixParticleModel<CloudType>(dict, owner, typeName, Xi),

    phiMod_m_
    (
        this->coeffDict().template lookupOrDefault<scalar>("phiMod_m", 0.01)
    ),
    sPx_m_
    (
        this->coeffDict().template lookupOrDefault<scalar>("sPx_m", 0.05)
    ),
    sPy_m_
    (
        this->coeffDict().template lookupOrDefault<scalar>("sPy_m", 0.05)
    ),
    sPz_m_
    (
        this->coeffDict().template lookupOrDefault<scalar>("sPz_m", 0.05)
    ),
    shadowOrigin_
    (
        this->coeffDict().template lookupOrDefault<vector>
        (
            "shadowOrigin",
            vector::zero
        )
    )
{
    printInfo();
}


template<class CloudType>
Foam::secondCondMMCcurl<CloudType>::secondCondMMCcurl
(
    const secondCondMMCcurl<CloudType>& cm
)
:
    mixParticleModel<CloudType>(cm),
    phiMod_m_(cm.phiMod_m_),
    sPx_m_   (cm.sPx_m_),
    sPy_m_   (cm.sPy_m_),
    sPz_m_   (cm.sPz_m_),
    shadowOrigin_(cm.shadowOrigin_)
{
    printInfo();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::secondCondMMCcurl<CloudType>::Smix()
{
    buildSecondCondList();
    this->SmixList();
}


template<class CloudType>
void Foam::secondCondMMCcurl<CloudType>::buildSecondCondList()
{
    this->particleList_.clear();
    this->eulerianFieldDataList_.clear();
    this->particlePairs_.clear();

    label particleInd = 0;

    forAllIters(this->owner(), iter)
    {
        // Only include particles flagged for second conditioning
        if (iter().secondCondFlag() != 1)
            continue;

        // Store pointer to the particle (no deallocation ownership)
        this->particleList_.append(iter.get());

        eulerianFieldData eulerianFields;

        const vector pos = iter().position();

        eulerianFields.particleIndex()  = particleInd++;
        eulerianFields.processorIndex() = Pstream::myProcNo();

        // Physical position stored for reference (not used as k-d tree split
        // axis — all splitting is done via XiR, see secondCondKdTree)
        eulerianFields.position() = pos;

        // Shadow coordinates: particle position relative to shadowOrigin.
        // When shadowOrigin = (0,0,0) these equal the physical coordinates.
        const vector xi = pos - shadowOrigin_;

        // 4D reference space for the k-d tree:
        //   XiR[0] = φ°  (phiModified, tight normalisation → dominant axis)
        //   XiR[1] = ξ_x (shadow x,    loose normalisation)
        //   XiR[2] = ξ_y (shadow y,    loose normalisation)
        //   XiR[3] = ξ_z (shadow z,    loose normalisation)
        eulerianFields.XiR().resize(4);
        eulerianFields.XiR()[0] = iter().phiModified();
        eulerianFields.XiR()[1] = xi.x();
        eulerianFields.XiR()[2] = xi.y();
        eulerianFields.XiR()[3] = xi.z();

        // magSqrRefVar is not used by secondCondMMCcurl (no time-scale
        // calculation based on Eulerian gradients); set to zero.
        eulerianFields.magSqrRefVar().resize(4, 0.0);

        // Eulerian transport properties — populated to satisfy the
        // eulerianFieldData interface; set to neutral values since mixpair
        // is a stub and does not use them yet.
        eulerianFields.DEff()   = 0.0;
        eulerianFields.D()      = 0.0;
        eulerianFields.Dt()     = 0.0;
        eulerianFields.DeltaE() = 1.0;
        eulerianFields.mu()     = 0.0;
        eulerianFields.vb()     = 1.0;
        eulerianFields.Rand()   = 0.0;

        this->eulerianFieldDataList_.append(std::move(eulerianFields));
    }

    // Build pairs from the flagged-particle list (local pairing only)
    findSecondCondPairs(this->eulerianFieldDataList_, this->particlePairs_);
}


template<class CloudType>
void Foam::secondCondMMCcurl<CloudType>::findSecondCondPairs
(
    const DynamicList<eulerianFieldData>& eulerianFieldList,
    DynamicList<List<label>>& pairs
) const
{
    pairs.clear();

    if (eulerianFieldList.size() < 2)
        return;

    std::vector<label> L, U;
    L.reserve(eulerianFieldList.size());
    U.reserve(eulerianFieldList.size());

    std::vector<label> pInd(eulerianFieldList.size());
    std::iota(pInd.begin(), pInd.end(), 0);

    secondCondKdTree
    (
        eulerianFieldList,
        1,
        eulerianFieldList.size(),
        pInd, L, U
    );

    pairs.reserve(static_cast<label>(std::ceil(0.5*eulerianFieldList.size())));

    for (std::size_t i = 0; i < L.size(); i++)
    {
        label p = L[i] - 1;
        label q = L[i];

        if (U[i] - L[i] < 2)
        {
            List<label> pair(2);
            pair[0] = pInd[p];
            pair[1] = pInd[q];
            pairs.append(std::move(pair));
        }
        else if (U[i] - L[i] == 2)
        {
            label r = L[i] + 1;
            List<label> pair(3);
            pair[0] = pInd[p];
            pair[1] = pInd[q];
            pair[2] = pInd[r];
            pairs.append(std::move(pair));
        }
    }
}


template<class CloudType>
void Foam::secondCondMMCcurl<CloudType>::secondCondKdTree
(
    const DynamicList<eulerianFieldData>& particleList,
    label l,
    label u,
    std::vector<label>& pInd,
    std::vector<label>& L,
    std::vector<label>& U
) const
{
    // Base case: group of 2 or 3 particles → record as a leaf pair/triple
    if (u - l <= 2)
    {
        L.push_back(l);
        U.push_back(u);
        return;
    }

    // Split point (adjusted so each sub-list has an even size)
    label m = (l + u) / 2;
    if ((u - m) % 2 != 0) m++;

    auto iterL = pInd.begin(); std::advance(iterL, l - 1);
    auto iterU = pInd.begin(); std::advance(iterU, u);

    // Normalisation factors in the same order as XiR:
    //   [0] phiMod_m_  (tight  → dominant)
    //   [1] sPx_m_
    //   [2] sPy_m_
    //   [3] sPz_m_
    const scalar Xii[4] = {phiMod_m_, sPx_m_, sPy_m_, sPz_m_};

    scalar maxXiR[4] = {-GREAT, -GREAT, -GREAT, -GREAT};
    scalar minXiR[4] = { GREAT,  GREAT,  GREAT,  GREAT};

    for (auto it = iterL; it != iterU; ++it)
    {
        for (label i = 0; i < 4; ++i)
        {
            const scalar val = particleList[*it].XiR()[i];
            maxXiR[i] = std::max(maxXiR[i], val);
            minXiR[i] = std::min(minXiR[i], val);
        }
    }

    // Select the dimension with the largest normalised range
    scalar disMax = 0.0;
    label  bestI  = 0;

    for (label i = 0; i < 4; ++i)
    {
        const scalar dis = mag(maxXiR[i] - minXiR[i]) / Xii[i];
        if (dis > disMax)
        {
            disMax = dis;
            bestI  = i;
        }
    }

    // Use lessArg with ncond = bestI + 3 so that the XiR() branch of the
    // comparator is always triggered (lessArg: ii_ < 3 → position, ii_ >= 3
    // → XiR()[ii_-3]).  Physical coordinates are therefore never used as
    // split dimensions, consistent with the first-conditioning MMCcurl.
    typename mixParticleModel<CloudType>::lessArg comp(bestI + 3);

    std::sort
    (
        iterL,
        iterU,
        [&](label A, label B) -> bool
        {
            return comp(particleList[A], particleList[B]);
        }
    );

    // Recursive calls on the two halves
    secondCondKdTree(particleList, l,   m,   pInd, L, U);
    secondCondKdTree(particleList, m+1, u,   pInd, L, U);
}


template<class CloudType>
void Foam::secondCondMMCcurl<CloudType>::mixpair
(
    particleType& p,
    const eulerianFieldData& /*pEulerianFields*/,
    particleType& q,
    const eulerianFieldData& /*qEulerianFields*/,
    scalar& /*deltaT*/
)
{
    // Second-conditioning mixing rule — to be implemented.
    // Particles p and q have been paired in (φ°, ξ_x, ξ_y, ξ_z) space;
    // the actual mixing operation on their scalars will be added here.
    (void)p;
    (void)q;
}


template<class CloudType>
const Foam::scalarField Foam::secondCondMMCcurl<CloudType>::XiR0
(
    label /*patch*/,
    label /*patchFace*/,
    particle& /*p*/
)
{
    // Not applicable: secondCondMMCcurl does not initialise boundary XiR values.
    return scalarField(0);
}


template<class CloudType>
const Foam::scalarField Foam::secondCondMMCcurl<CloudType>::XiR0
(
    label /*celli*/,
    particle& /*p*/
)
{
    // Not applicable: secondCondMMCcurl does not initialise cell XiR values.
    return scalarField(0);
}


template<class CloudType>
void Foam::secondCondMMCcurl<CloudType>::printInfo()
{
    Info<< "Mixing Model: " << this->modelType() << nl
        << token::TAB << "Reference space: (phiModified, xi_x, xi_y, xi_z)" << nl
        << token::TAB << "phiMod_m:      " << phiMod_m_    << nl
        << token::TAB << "sPx_m:         " << sPx_m_       << nl
        << token::TAB << "sPy_m:         " << sPy_m_       << nl
        << token::TAB << "sPz_m:         " << sPz_m_       << nl
        << token::TAB << "shadowOrigin:  " << shadowOrigin_ << nl
        << token::TAB << "---------------------------" << nl
        << token::TAB << "Particle filter: secondCondFlag == 1 only" << nl
        << token::TAB << "k-d tree: splits on XiR (no physical-coord axis)"
        << endl;
}


// ************************************************************************* //
