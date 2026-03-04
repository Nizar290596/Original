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

Description
    Two-stage sequential second-conditioning MMC mixing model for sparse
    premixed combustion simulations.

    Implementation of the method described in:
        Sundaram & Klimenko, Proc. Combust. Inst. 36 (2017) 1937-1945.
        DOI: 10.1016/j.proci.2016.07.116

    Stage 1 — full particle set, shadow-position pairing,
               standard MMCcurl mixing of all XiC except c.
               c is recomputed from T after Stage 1.
    Stage 2 — flame-zone subset (c in [cMin, cMax]),
               combined [sP, phi*] pairing, mix all XiC including c.

\*---------------------------------------------------------------------------*/

#include "sparsePremixedMMCcurl.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::sparsePremixedMMCcurl<CloudType>::sparsePremixedMMCcurl
(
    const dictionary& dict,
    CloudType& owner,
    const mmcVarSet& Xi
)
:
    mixParticleModel<CloudType>(dict, owner, typeName, Xi),

    Tu_(this->coeffDict().lookupOrDefault("Tu", 298.0)),
    Tb_(this->coeffDict().lookupOrDefault("Tb", 2240.0)),
    cMin_(this->coeffDict().lookupOrDefault("cMin", 0.05)),
    cMax_(this->coeffDict().lookupOrDefault("cMax", 0.95)),
    CL_(this->coeffDict().lookupOrDefault("CL", 0.1)),
    beta_(this->coeffDict().lookupOrDefault("beta", 1.0)),

    couplingVarName_
    (
        this->coeffDict().lookupOrDefault<word>("couplingVarName", "c")
    ),

    wOUVarName_
    (
        this->coeffDict().lookupOrDefault<word>("wOUVarName", "wOU")
    ),

    numShadowPos_
    (
        this->coeffDict().lookupOrDefault<label>("numShadowPos", 3)
    ),

    limitkdTree_
    (
        this->coeffDict().lookupOrDefault("limitkdTree", true)
    ),

    meanTimeScale_
    (
        this->coeffDict().lookup("meanTimeScale")
    ),

    cIndex_(-1),
    wOUIndex_(-1)
{
    // Resolve index of wOU in XiR
    const auto& rVarMap = this->XiR_.rVarInXiR();
    if (rVarMap.found(wOUVarName_))
    {
        wOUIndex_ = rVarMap[wOUVarName_];
    }
    else
    {
        FatalErrorIn("sparsePremixedMMCcurl::sparsePremixedMMCcurl")
            << "Reference variable '" << wOUVarName_
            << "' not found in mmcVariablesDefinitions." << nl
            << "Available reference variables: " << rVarMap.toc()
            << exit(FatalError);
    }

    // Resolve index of c in XiC (coupling variable set)
    const auto& cVarMap = owner.coupling().XiC().cVarInXiC();
    if (cVarMap.found(couplingVarName_))
    {
        cIndex_ = cVarMap[couplingVarName_];
    }
    else
    {
        FatalErrorIn("sparsePremixedMMCcurl::sparsePremixedMMCcurl")
            << "Coupling variable '" << couplingVarName_
            << "' not found in mmcVariablesDefinitions." << nl
            << "Available coupling variables: " << cVarMap.toc()
            << exit(FatalError);
    }

    printInfo();
}


template<class CloudType>
Foam::sparsePremixedMMCcurl<CloudType>::sparsePremixedMMCcurl
(
    const sparsePremixedMMCcurl<CloudType>& cm
)
:
    mixParticleModel<CloudType>(cm),
    Tu_(cm.Tu_),
    Tb_(cm.Tb_),
    cMin_(cm.cMin_),
    cMax_(cm.cMax_),
    CL_(cm.CL_),
    beta_(cm.beta_),
    couplingVarName_(cm.couplingVarName_),
    wOUVarName_(cm.wOUVarName_),
    numShadowPos_(cm.numShadowPos_),
    limitkdTree_(cm.limitkdTree_),
    meanTimeScale_(cm.meanTimeScale_),
    cIndex_(cm.cIndex_),
    wOUIndex_(cm.wOUIndex_)
{
    printInfo();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::sparsePremixedMMCcurl<CloudType>::Smix()
{
    // ------------------------------------------------------------------
    // Build the global particle list.
    // eulerianFieldDataList_ receives XiR = [sPx, sPy, sPz, wOU] and
    // XiC for every particle (local + gathered remote processors in parallel).
    // c is NOT pre-computed here; it is derived from temperature after
    // Stage 1 so that it reflects the post-Stage-1 mixed state.
    // ------------------------------------------------------------------
    this->buildParticleList();

    const label N = this->eulerianFieldDataList_.size();
    if (N < 2)
        return;

    const scalar deltaT = this->owner().mesh().time().deltaT().value();

    // ------------------------------------------------------------------
    // STAGE 1 — Full particle set, shadow-position pairing.
    //           Standard MMCcurl mixing of all scalars EXCEPT c.
    //
    // XiR is restricted to [sPx, sPy, sPz] (first numShadowPos_ entries)
    // so the pairing kdTree operates purely in shadow-position space.
    // c is deliberately excluded from mixing here; it will be re-derived
    // from the mixed temperature field immediately after this stage.
    // ------------------------------------------------------------------
    DynamicList<eulerianFieldData> stage1List(this->eulerianFieldDataList_);

    forAll(stage1List, i)
    {
        scalarField sp(numShadowPos_);
        for (label k = 0; k < numShadowPos_; k++)
            sp[k] = this->eulerianFieldDataList_[i].XiR()[k];
        stage1List[i].XiR() = sp;
    }

    DynamicList<List<label>> stage1Pairs;
    this->findPairs(stage1List, stage1Pairs);

    smixStage1(stage1Pairs, deltaT);

    // ------------------------------------------------------------------
    // Post-Stage-1 — Recompute c from the now-mixed particle temperatures.
    //
    //   c = (T_p - Tu) / (Tb - Tu)  clamped to [0, 1]      (eq. 1)
    //
    // Update local particle objects, then sync the new c values into
    // eulerianFieldDataList_ so Stage 2 uses post-Stage-1 c values when
    // computing phi* and filtering the flame-zone subset.
    // Remote particle entries in eulerianFieldDataList_ retain their
    // pre-Stage-1 c (their temperature was mixed on their home processor).
    // ------------------------------------------------------------------
    forAllIter(CloudType, this->owner(), pIter)
    {
        particleType& p = pIter();
        p.XiC()[cIndex_] = progressVar(p.T());
    }

    forAll(this->eulerianFieldDataList_, i)
    {
        if (this->eulerianFieldDataList_[i].local())
        {
            const label pIdx =
                this->eulerianFieldDataList_[i].particleIndex();
            this->eulerianFieldDataList_[i].XiC()[cIndex_] =
                this->particleList_[pIdx]->XiC()[cIndex_];
        }
    }

    // ------------------------------------------------------------------
    // STAGE 2 — Flame-zone subset, second-conditioning pairing.
    //           Pairing in combined [sPx, sPy, sPz, phi*] space (4D).
    //           Mixes all coupling scalars (including c).
    //
    //   Subset criterion : c in [cMin_, cMax_]
    //   Modified progress variable: phi* = c * exp(beta_ * wOU)
    //
    // Pairing uses BOTH shadow positions and phi* simultaneously so that
    // only particles close in both physical reference space and flame
    // structure are mixed together (second conditioning, S&K 2017).
    // The flame-zone filter [cMin, cMax] is the sparsity mechanism.
    // ------------------------------------------------------------------
    DynamicList<eulerianFieldData> stage2List;

    forAll(this->eulerianFieldDataList_, i)
    {
        const eulerianFieldData& ef = this->eulerianFieldDataList_[i];

        const scalar c_p   = ef.XiC()[cIndex_];
        const scalar wOU_p = ef.XiR()[wOUIndex_];

        if (c_p < cMin_ || c_p > cMax_)
            continue;

        eulerianFieldData efSubset = ef;  // copy (preserves particleIndex)

        const scalar phiStar = c_p * Foam::exp(beta_ * wOU_p);

        // Combined reference: [sPx, sPy, sPz, phi*]
        scalarField combinedRef(numShadowPos_ + 1);
        for (label k = 0; k < numShadowPos_; k++)
            combinedRef[k] = ef.XiR()[k];
        combinedRef[numShadowPos_] = phiStar;
        efSubset.XiR() = combinedRef;

        stage2List.append(efSubset);
    }

    if (stage2List.size() < 2)
        return;

    DynamicList<List<label>> stage2Pairs;
    this->findPairs(stage2List, stage2Pairs);

    // Pair indices in stage2Pairs resolve into stage2List, NOT
    // eulerianFieldDataList_. smixStage2 receives stage2List explicitly.
    smixStage2(stage2List, stage2Pairs, deltaT);
}


// -----------------------------------------------------------------------
// Stage 1 helper — pairs indexed into eulerianFieldDataList_
// -----------------------------------------------------------------------

template<class CloudType>
void Foam::sparsePremixedMMCcurl<CloudType>::smixStage1
(
    const DynamicList<List<label>>& pairs,
    scalar deltaT
)
{
    for (const auto& pair : pairs)
    {
        const label sz = pair.size();
        if (sz == 2)
        {
            const eulerianFieldData& e1 =
                this->eulerianFieldDataList_[pair[0]];
            const eulerianFieldData& e2 =
                this->eulerianFieldDataList_[pair[1]];

            if (e1.local() || e2.local())
                mixpairStaged
                (
                    *this->particleList_[e1.particleIndex()], e1,
                    *this->particleList_[e2.particleIndex()], e2,
                    deltaT, 1
                );
        }
        else if (sz == 3)
        {
            const eulerianFieldData& e1 =
                this->eulerianFieldDataList_[pair[0]];
            const eulerianFieldData& e2 =
                this->eulerianFieldDataList_[pair[1]];
            const eulerianFieldData& e3 =
                this->eulerianFieldDataList_[pair[2]];

            if (e1.local() || e2.local())
                mixpairStaged
                (
                    *this->particleList_[e1.particleIndex()], e1,
                    *this->particleList_[e2.particleIndex()], e2,
                    deltaT, 1
                );
            if (e2.local() || e3.local())
                mixpairStaged
                (
                    *this->particleList_[e2.particleIndex()], e2,
                    *this->particleList_[e3.particleIndex()], e3,
                    deltaT, 1
                );
        }
    }
}


// -----------------------------------------------------------------------
// Stage 2 helper — pairs indexed into subsetList (the flame-zone subset)
// -----------------------------------------------------------------------

template<class CloudType>
void Foam::sparsePremixedMMCcurl<CloudType>::smixStage2
(
    const DynamicList<eulerianFieldData>& subsetList,
    const DynamicList<List<label>>& pairs,
    scalar deltaT
)
{
    for (const auto& pair : pairs)
    {
        const label sz = pair.size();
        if (sz == 2)
        {
            // Pairs are indexed into subsetList, not eulerianFieldDataList_.
            // particleIndex() still refers into particleList_.
            const eulerianFieldData& e1 = subsetList[pair[0]];
            const eulerianFieldData& e2 = subsetList[pair[1]];

            if (e1.local() || e2.local())
                mixpairStaged
                (
                    *this->particleList_[e1.particleIndex()], e1,
                    *this->particleList_[e2.particleIndex()], e2,
                    deltaT, 2
                );
        }
        else if (sz == 3)
        {
            const eulerianFieldData& e1 = subsetList[pair[0]];
            const eulerianFieldData& e2 = subsetList[pair[1]];
            const eulerianFieldData& e3 = subsetList[pair[2]];

            if (e1.local() || e2.local())
                mixpairStaged
                (
                    *this->particleList_[e1.particleIndex()], e1,
                    *this->particleList_[e2.particleIndex()], e2,
                    deltaT, 2
                );
            if (e2.local() || e3.local())
                mixpairStaged
                (
                    *this->particleList_[e2.particleIndex()], e2,
                    *this->particleList_[e3.particleIndex()], e3,
                    deltaT, 2
                );
        }
    }
}


// -----------------------------------------------------------------------
// Low-level staged mixing for one particle pair
// -----------------------------------------------------------------------

template<class CloudType>
void Foam::sparsePremixedMMCcurl<CloudType>::mixpairStaged
(
    particleType& p,
    const eulerianFieldData& /*pEulFields*/,
    particleType& q,
    const eulerianFieldData& /*qEulFields*/,
    scalar deltaT,
    label stageFlag
)
{
    if (p.wt() + q.wt() <= 0)
        return;

    // LPF mixing time scale from particle's stored mixing time field
    scalar tauP = CL_ * p.mixTime();
    scalar tauQ = CL_ * q.mixTime();

    if (tauP >= 1e30 || tauQ >= 1e30)
        return;

    // Optionally skip pairs too far apart in reference space
    if (limitkdTree_)
    {
        forAll(p.dXiR(), i)
        {
            if (p.dXiR()[i] > 2.0 * this->Xii_[i])
                return;
        }
    }

    scalar tauMix = 0.0;
    if (meanTimeScale_)
        tauMix = 2.0 / (1.0/(tauP + VSMALL) + 1.0/(tauQ + VSMALL));
    else
        tauMix = Foam::min(tauP, tauQ);

    const scalar mixExtent = 1.0 - Foam::exp(-deltaT / (tauMix + VSMALL));
    const scalar wtSum     = p.wt() + q.wt();

    if (stageFlag == 1)
    {
        // Stage 1: standard MMCcurl step — mix all coupling scalars except c.
        // c is not touched here; it is recomputed from temperature after Stage 1.
        forAll(p.XiC(), nn)
        {
            if (nn == cIndex_)
                continue;

            scalar av =
                (p.wt() * p.XiC()[nn] + q.wt() * q.XiC()[nn]) / wtSum;

            p.XiC()[nn] += mixExtent * (av - p.XiC()[nn]);
            q.XiC()[nn] += mixExtent * (av - q.XiC()[nn]);
        }
    }
    else
    {
        // Stage 2: second-conditioning step — mix all coupling scalars,
        // including c (which was freshly derived from the post-Stage-1
        // temperature before this stage runs).
        forAll(p.XiC(), nn)
        {
            scalar av =
                (p.wt() * p.XiC()[nn] + q.wt() * q.XiC()[nn]) / wtSum;

            p.XiC()[nn] += mixExtent * (av - p.XiC()[nn]);
            q.XiC()[nn] += mixExtent * (av - q.XiC()[nn]);
        }
    }

    p.mixTime() = tauMix;
    q.mixTime() = tauMix;
    p.mixExt()  = mixExtent;
    q.mixExt()  = mixExtent;
}


template<class CloudType>
void Foam::sparsePremixedMMCcurl<CloudType>::mixpair
(
    particleType& p,
    const eulerianFieldData& pEulFields,
    particleType& q,
    const eulerianFieldData& qEulFields,
    scalar& deltaT
)
{
    // Fallback when called directly (not through the two-stage Smix).
    mixpairStaged(p, pEulFields, q, qEulFields, deltaT, 2);
}


template<class CloudType>
const Foam::scalarField
Foam::sparsePremixedMMCcurl<CloudType>::XiR0(label patch, label patchFace)
{
    const mmcVarSet& setOfXi    = this->XiR();
    const labelHashTable& XiIdx = setOfXi.rVarInXi();
    const labelHashTable& XiRIdx= setOfXi.rVarInXiR();

    scalarField XXi(this->numXiR(), 0.0);

    for (const word& varName : this->XiRNames_)
    {
        if (setOfXi.Vars(XiIdx[varName]).type() == "evolved")
            XXi[XiRIdx[varName]] = 0.0;
        else
            XXi[XiRIdx[varName]] =
                setOfXi.Vars(XiIdx[varName])
                    .field().boundaryField()[patch][patchFace];
    }
    return XXi;
}


template<class CloudType>
const Foam::scalarField
Foam::sparsePremixedMMCcurl<CloudType>::XiR0(label celli)
{
    const mmcVarSet& setOfXi    = this->XiR();
    const labelHashTable& XiIdx = setOfXi.rVarInXi();
    const labelHashTable& XiRIdx= setOfXi.rVarInXiR();

    scalarField XXi(this->numXiR(), 0.0);

    for (const word& varName : this->XiRNames_)
    {
        if (setOfXi.Vars(XiIdx[varName]).type() == "evolved")
            XXi[XiRIdx[varName]] = 0.0;
        else
            XXi[XiRIdx[varName]] =
                setOfXi.Vars(XiIdx[varName]).field()[celli];
    }
    return XXi;
}


template<class CloudType>
void Foam::sparsePremixedMMCcurl<CloudType>::printInfo()
{
    Info << "Mixing Model: sparsePremixedMMCcurl" << nl
         << token::TAB
         << "(two-stage second conditioning, Sundaram & Klimenko 2017)" << nl
         << token::TAB << "Stage 1 — full set, shadow-position pairing:" << nl
         << token::TAB << "  CL:                  " << CL_              << nl
         << token::TAB << "  mixes all XiC except: " << couplingVarName_ << nl
         << token::TAB << "  (c recomputed from T after Stage 1)"        << nl
         << token::TAB << "Stage 2 — flame-zone subset, [sP+phi*] pairing:" << nl
         << token::TAB << "  Tu (unburned) [K]:   " << Tu_              << nl
         << token::TAB << "  Tb (burned)   [K]:   " << Tb_              << nl
         << token::TAB << "  cMin (subset lower): " << cMin_            << nl
         << token::TAB << "  cMax (subset upper): " << cMax_            << nl
         << token::TAB << "  beta (phi* coeff):   " << beta_            << nl
         << token::TAB << "  progress variable:   " << couplingVarName_ << nl
         << token::TAB << "    -> XiC index:      " << cIndex_          << nl
         << token::TAB << "  OU variable:         " << wOUVarName_       << nl
         << token::TAB << "    -> XiR index:      " << wOUIndex_         << nl
         << token::TAB << "  shadow pos dims:     " << numShadowPos_     << nl
         << token::TAB << "  pairing space dims:  " << numShadowPos_ + 1 << nl
         << token::TAB << "ri:                    " << this->ri_         << nl
         << token::TAB << "Xii:                   " << this->Xii_        << nl
         << token::TAB << "---------------------------"                   << nl;
    if (meanTimeScale_)
        Info << token::TAB
             << "Time scale: harmonic mean of pair values" << nl;
    else
        Info << token::TAB
             << "Time scale: minimum of pair values" << nl;
    if (limitkdTree_)
        Info << token::TAB
             << "Particles with dXiR > 2*Xii are NOT mixed" << endl;
    else
        Info << token::TAB
             << "No limiting of mixing based on dXiR" << endl;
}


// ************************************************************************* //
