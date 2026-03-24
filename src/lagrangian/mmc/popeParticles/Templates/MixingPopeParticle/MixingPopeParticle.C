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

\*---------------------------------------------------------------------------*/

#include "MixingPopeParticle.H"
#include "physicoChemicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * Member Functions * * * * * * * * * * * * //

template<class ParticleType>
template<class TrackCloudType>
void Foam::MixingPopeParticle<ParticleType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const label cellI
)
{
    ParticleType::setCellValues(cloud, td, dt, cellI);
    
    tetIndices tetIs = this->currentTetIndices();
    
    if (cloud.mixing().XiRNames().size())
    {
        const mmcVarSet& setOfXi(cloud.mixing().XiR());
        
        const labelHashTable& indexInXi  = setOfXi.rVarInXi();
        const labelHashTable& indexInXiR = setOfXi.rVarInXiR();
        
        forAllConstIter(wordList, cloud.mixing().XiRNames(), iter)
        {

            //- This part requires review
            //- (For now it works fine just for interpolated)

            //if (setOfXi.Vars(indexInXi[*iter]).refType()=="interpolated")
            //    XiR()[indexInXiR[*iter]] = td.XiRInterp()[indexInXiR[*iter]].interpolate(this->coordinates(),tetIs);
            //else if (setOfXi.Vars(indexInXi[*iter]).refType()=="evolved")
            //{
                // Do nothing since it is evolved by MMC model
		scalar XiROld = XiR(*iter);
		scalar XiRNew = setOfXi.Vars(indexInXi[*iter]).evolveMethod().compute(*this,dt,XiROld);
		XiR(*iter) = XiRNew;
		//Info << "XiROld " << XiROld << endl;
		//Info << "XiRNew " << XiRNew << endl;
           // }
           // else
           //     XiR()[indexInXiR[*iter]] = 0.0;//this->position();
        }
    }

    // Advance Ornstein-Uhlenbeck state for second-conditioning particles
    applyOUProcessUpdate(*this, cloud, dt);

    // W(φ): chemical source term  dφ/dt = A·(1−φ)·exp[Z·(φ−1)]
    {
        const scalar A = cloud.secondCondAPhi();
        const scalar Z = cloud.secondCondZPhi();
        phi_ += dt * A * (1.0 - phi_) * Foam::exp(Z * (phi_ - 1.0));
        phi_ = max(0.0, min(1.0, phi_));
    }

    // Recompute φ° = φ·exp(β·ω_OU) after OU and W(φ) updates
    phiModified_ = phi_ * Foam::exp(cloud.secondCondBeta() * omegaOU_);
}


template<class ParticleType>
template<class TrackCloudType>
void Foam::MixingPopeParticle<ParticleType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const label cellI
)
{
    ParticleType::calc(cloud, td, dt, cellI);
}

// * * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * * //

// Helper: apply IEM phi mixing and clamp to [0,1]
namespace
{
    template<class P>
    void mixPhiIEM(P& p, P& q, scalar mixExtent)
    {
        const scalar phiMean =
            (p.wt() * p.phi() + q.wt() * q.phi()) / (p.wt() + q.wt());
        p.phi() += mixExtent * (phiMean - p.phi());
        q.phi() += mixExtent * (phiMean - q.phi());
        p.phi() = max(0.0, min(1.0, p.phi()));
        q.phi() = max(0.0, min(1.0, q.phi()));
    }
}


template<class ParticleType>
void Foam::MixingPopeParticle<ParticleType>::mixProperties
(
    MixingPopeParticle<ParticleType>& p,
    MixingPopeParticle<ParticleType>& q,
    scalar mixExtent
)
{
    ParticleType::mixProperties(p, q, mixExtent);

    // S(φ): first-conditioning IEM mixing of progress variable
    mixPhiIEM(p, q, mixExtent);
}


template<class ParticleType>
void Foam::MixingPopeParticle<ParticleType>::mixProperties
(
    MixingPopeParticle<ParticleType>& p,
    MixingPopeParticle<ParticleType>& q,
    const scalar& mixExtent,
    const scalar& mixExtentSoot
)
{
    ParticleType::mixProperties(p, q, mixExtent,mixExtentSoot);

    // S(φ): first-conditioning IEM mixing of progress variable
    mixPhiIEM(p, q, mixExtent);
}


template<class ParticleType>
void Foam::MixingPopeParticle<ParticleType>::mixProperties
(
    MixingPopeParticle<ParticleType>& p,
    MixingPopeParticle<ParticleType>& q,
    scalar mixExtent,
    scalarList ScaledExtent
)
{
    ParticleType::mixProperties(p, q, mixExtent,ScaledExtent);

    // S(φ): first-conditioning IEM mixing of progress variable
    mixPhiIEM(p, q, mixExtent);
}


template<class ParticleType>
template<class CloudType>
void Foam::MixingPopeParticle<ParticleType>::setStaticProperties(CloudType& c)
{
    ParticleType::setStaticProperties(c);
    
    if (c.mixing().numXiR())
    {
        indexInXiR_ = c.mixing().XiR().rVarInXiR();
        XiRNames_ = c.mixing().XiRNames();
    }
}


template<class ParticleType>
Foam::DynamicList<Foam::scalar> 
Foam::MixingPopeParticle<ParticleType>::getStatisticalData
(
    const wordList& vars
) const
{
    return ParticleType::getStatisticalData(vars);
}


template<class ParticleType>
Foam::DynamicList<Foam::word>
Foam::MixingPopeParticle<ParticleType>::getStatisticalDataNames
(
    const wordList& vars
) const
{
    return ParticleType::getStatisticalDataNames(vars);
}


template<class ParticleType>
void Foam::MixingPopeParticle<ParticleType>::initStatisticalSampling()
{
    ParticleType::initStatisticalSampling();

    // Add physical-space mixing distance
    this->nameVariableLookUpTable().addNamedVariable("dx", dx_);

    // Second conditioning scalars
    this->nameVariableLookUpTable().addNamedVariable("omegaOU",     omegaOU_);
    this->nameVariableLookUpTable().addNamedVariable("phi",         phi_);
    this->nameVariableLookUpTable().addNamedVariable("phiModified", phiModified_);

    // Note: XiRNames is initialized as a static variablebefore 
    // XiR_ is set. Therefore we need to check if they are set to avoid
    // accessing elements out of bounds

    // Add reference variables
    if (XiRNames_.size() == XiR_.size())
    {
        for (const word& name : XiRNames_)
        {
            this->nameVariableLookUpTable().addNamedVariable
            (
                name,
                XiR(name)
            );

            word dXirName = "d"+name;
            this->nameVariableLookUpTable().addNamedVariable
            (
                dXirName,
                dXiR(name)
            );
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::MixingPopeParticle<ParticleType>::MixingPopeParticle
(
    const MixingPopeParticle<ParticleType>& p
)
:
    ParticleType(p),
    XiR_(p.XiR_),
    dXiR_(p.dXiR_),
    dx_(p.dx_),
    secondCondFlag_(p.secondCondFlag_),
    omegaOU_(p.omegaOU_),
    phi_(p.phi_),
    phiModified_(p.phiModified_)
{
    initStatisticalSampling();
}


template<class ParticleType>
Foam::MixingPopeParticle<ParticleType>::MixingPopeParticle
(
    const MixingPopeParticle<ParticleType>& p,
    const polyMesh& mesh
)
:
    ParticleType(p, mesh),
    XiR_(p.XiR_),
    dXiR_(p.dXiR_),
    dx_(p.dx_),
    secondCondFlag_(p.secondCondFlag_),
    omegaOU_(p.omegaOU_),
    phi_(p.phi_),
    phiModified_(p.phiModified_)
{
    initStatisticalSampling();
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "MixingPopeParticleIO.C"


// ************************************************************************* //
