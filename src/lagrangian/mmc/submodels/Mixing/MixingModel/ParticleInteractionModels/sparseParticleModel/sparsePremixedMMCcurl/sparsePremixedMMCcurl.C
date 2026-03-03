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
    Two-stage sequential conditioning MMC mixing model for sparse premixed
    simulations.

    Implementation of the method described in:
        Sundaram & Klimenko, Proc. Combust. Inst. 36 (2017) 1937-1945.
        DOI: 10.1016/j.proci.2016.07.116

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
    // The coupling variable set is accessible through the cloud.
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
    // Build the global particle list (handles parallel communication and
    // sets eulerianFieldDataList_ with XiR = [sPx, sPy, sPz, wOU]).
    // buildParticleList() also calls findPairs() internally, but we
    // override the pairs below for each stage.
    // ------------------------------------------------------------------
    this->buildParticleList();

    const scalar deltaT = this->owner().mesh().time().deltaT().value();

    // ------------------------------------------------------------------
    // Stage 1: pair in [position + shadow-position] space, mix only c
    // ------------------------------------------------------------------

    // Build a modified copy of the eulerian list with XiR restricted to
    // the shadow-position components only (indices 0..numShadowPos_-1).
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
    // Stage 2: pair in [position + phi*] space, mix reactive scalars
    //
    //   phi* = c * exp(beta * wOU)
    //
    // c is the Eulerian progress variable interpolated at each particle's
    // position (consistent for both local and remote particles in
    // parallel runs).
    // ------------------------------------------------------------------

    // Interpolate the Eulerian c field at each particle's position.
    const volScalarField& cField =
        this->owner().mesh().objectRegistry::
            lookupObject<volScalarField>(couplingVarName_);

    interpolationCellPoint<scalar> cInterp(cField);

    DynamicList<eulerianFieldData> stage2List(this->eulerianFieldDataList_);

    forAll(stage2List, i)
    {
        const vector& pos  = stage2List[i].position();
        const label   cell = this->owner().mesh().findCell(pos);
        const label   face = -1;  // face not needed for cellPoint interpolation

        scalar c_p = cInterp.interpolate(pos, cell, face);
        c_p = Foam::max(c_p, VSMALL);   // ensure positivity

        scalar wOU_p = this->eulerianFieldDataList_[i].XiR()[wOUIndex_];

        scalar phiStar = c_p * Foam::exp(beta_ * wOU_p);

        scalarField ps(1, phiStar);
        stage2List[i].XiR() = ps;
    }

    DynamicList<List<label>> stage2Pairs;
    this->findPairs(stage2List, stage2Pairs);

    smixStage2(stage2Pairs, deltaT);
}


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
            const eulerianFieldData& e1 = this->eulerianFieldDataList_[pair[0]];
            const eulerianFieldData& e2 = this->eulerianFieldDataList_[pair[1]];
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
            const eulerianFieldData& e1 = this->eulerianFieldDataList_[pair[0]];
            const eulerianFieldData& e2 = this->eulerianFieldDataList_[pair[1]];
            const eulerianFieldData& e3 = this->eulerianFieldDataList_[pair[2]];
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


template<class CloudType>
void Foam::sparsePremixedMMCcurl<CloudType>::smixStage2
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
            const eulerianFieldData& e1 = this->eulerianFieldDataList_[pair[0]];
            const eulerianFieldData& e2 = this->eulerianFieldDataList_[pair[1]];
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
            const eulerianFieldData& e1 = this->eulerianFieldDataList_[pair[0]];
            const eulerianFieldData& e2 = this->eulerianFieldDataList_[pair[1]];
            const eulerianFieldData& e3 = this->eulerianFieldDataList_[pair[2]];
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

    // LPF mixing time scale (from particle's stored mixing time field)
    scalar tauP = CL_ * p.mixTime();
    scalar tauQ = CL_ * q.mixTime();

    if (tauP >= 1e30 || tauQ >= 1e30)
        return;

    // Optionally skip pairs that are too far apart in reference space
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

    scalar mixExtent = 1.0 - Foam::exp(-deltaT / (tauMix + VSMALL));

    const scalar wtSum = p.wt() + q.wt();

    if (stageFlag == 1)
    {
        // Stage 1: mix only the progress variable c
        scalar cAv =
            (p.wt() * p.XiC()[cIndex_] + q.wt() * q.XiC()[cIndex_]) / wtSum;

        p.XiC()[cIndex_] += mixExtent * (cAv - p.XiC()[cIndex_]);
        q.XiC()[cIndex_] += mixExtent * (cAv - q.XiC()[cIndex_]);
    }
    else
    {
        // Stage 2: mix all coupling scalars except c
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
    // Default single-stage full mixing (used when Smix() is not overridden
    // by derived classes). Calls Stage 2 (all scalars) as fallback.
    mixpairStaged(p, pEulFields, q, qEulFields, deltaT, 2);
}


template<class CloudType>
const Foam::scalarField
Foam::sparsePremixedMMCcurl<CloudType>::XiR0(label patch, label patchFace)
{
    const mmcVarSet& setOfXi = this->XiR();

    const labelHashTable& XiIndexes  = setOfXi.rVarInXi();
    const labelHashTable& XiRIndexes = setOfXi.rVarInXiR();

    scalarField XXi(this->numXiR(), 0.0);

    for (const word& varName : this->XiRNames_)
    {
        if (setOfXi.Vars(XiIndexes[varName]).type() == "evolved")
        {
            XXi[XiRIndexes[varName]] = 0.0;
        }
        else
        {
            XXi[XiRIndexes[varName]] =
                setOfXi.Vars(XiIndexes[varName])
                    .field().boundaryField()[patch][patchFace];
        }
    }
    return XXi;
}


template<class CloudType>
const Foam::scalarField
Foam::sparsePremixedMMCcurl<CloudType>::XiR0(label celli)
{
    const mmcVarSet& setOfXi = this->XiR();

    const labelHashTable& XiIndexes  = setOfXi.rVarInXi();
    const labelHashTable& XiRIndexes = setOfXi.rVarInXiR();

    scalarField XXi(this->numXiR(), 0.0);

    for (const word& varName : this->XiRNames_)
    {
        if (setOfXi.Vars(XiIndexes[varName]).type() == "evolved")
        {
            XXi[XiRIndexes[varName]] = 0.0;
        }
        else
        {
            XXi[XiRIndexes[varName]] =
                setOfXi.Vars(XiIndexes[varName]).field()[celli];
        }
    }
    return XXi;
}


template<class CloudType>
void Foam::sparsePremixedMMCcurl<CloudType>::printInfo()
{
    Info << "Mixing Model: sparsePremixedMMCcurl (two-stage second conditioning)" << nl
         << token::TAB << "CL:                  " << CL_              << nl
         << token::TAB << "beta (phi* mod):     " << beta_            << nl
         << token::TAB << "Progress variable:   " << couplingVarName_ << nl
         << token::TAB << "  -> XiC index:      " << cIndex_          << nl
         << token::TAB << "OU variable:         " << wOUVarName_       << nl
         << token::TAB << "  -> XiR index:      " << wOUIndex_         << nl
         << token::TAB << "Shadow pos dims:     " << numShadowPos_     << nl
         << token::TAB << "ri:                  " << this->ri_         << nl
         << token::TAB << "Xi:                  " << this->Xii_        << nl
         << token::TAB << "---------------------------"                 << nl;
    if (meanTimeScale_)
        Info << token::TAB
             << "Use harmonic mean of two particles' mixing time scales" << nl;
    else
        Info << token::TAB
             << "Use the min of two particles' mixing time scales" << nl;
    if (limitkdTree_)
        Info << token::TAB
             << "Particles with dXiR > 2*Xii are NOT mixed" << endl;
    else
        Info << token::TAB
             << "No limiting of mixing based on dXiR" << endl;
}


// ************************************************************************* //
