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

#include "MixingPopeCloud.H"
#include "CloudMixingModel.H"


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::MixingPopeCloud<CloudType>::setModels(const mmcVarSet& Xi)
{
    mixingModel_.reset
    (
        CloudMixingModel<MixingPopeCloud<CloudType> >::New
        (
            this->subModelProperties(),
            *this,
            Xi //Since reference variables are inside submodel!!!
        ).ptr()
    );
}


template<class CloudType>
void Foam::MixingPopeCloud<CloudType>::cloudReset(MixingPopeCloud<CloudType>& c)
{
    CloudType::cloudReset(c);
    
    mixingModel_.reset(c.mixingModel_.ptr());
    
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MixingPopeCloud<CloudType>::MixingPopeCloud
(
    const word& cloudName,
    const fvMesh& mesh,
    const volVectorField& U,
    const volScalarField& DEff,
    const volScalarField& rho,
    const volVectorField& gradRho,
    const mmcVarSet& Xi,
    const Switch initAtCnstr,
    bool readFields
)
:
    CloudType
    (
        cloudName,
        mesh,
        U,
        DEff,
        rho,
        gradRho,
        Xi,
        false,  // Only the top level cloud will call initAtCnstr
        readFields
    ),

    mixingPopeCloud(),

    cloudCopyPtr_(nullptr),

    mixingModel_(nullptr),

    secondCondR_(0.0),
    secondCondBeta_(1.0),
    secondCondTauOU_(1.0)

{
    // Read second conditioning parameters — done in the body so that
    // subOrEmptyDict() returns a concrete dictionary& (no dependent-type
    // template-disambiguation issue as would occur in the initializer list).
    {
        const dictionary& scDict =
            this->cloudProperties_.subOrEmptyDict("secondConditioning");
        secondCondR_     = scDict.lookupOrDefault("R",     scalar(0.0));
        secondCondBeta_  = scDict.lookupOrDefault("beta",  scalar(1.0));
        secondCondTauOU_ = scDict.lookupOrDefault("tauOU", scalar(1.0));
    }

    Info << "Creating mixing Pope Particle Cloud." << nl << endl;

    setModels(Xi); // passing mmcVarSet

    Info << nl << "Mixing model constructed." << endl;

    setEulerianStatistics();

    if(readFields)
    {
        if(this->size()>0)
        {
            Info << nl << "Reading Mixing Pope particle cloud data from file." << endl;

            particleType::mixingParticleIOType::readFields(*this, this->mixing());
        }

        else
        {
            if(initAtCnstr)
            {
                Info << "Initial realease of Pope particles into the finite volume field." << nl << endl;

                this->initReleaseParticles();
            }
        }
    }
}


template<class CloudType>
Foam::MixingPopeCloud<CloudType>::MixingPopeCloud
(
    const word& cloudName,
    const fvMesh& mesh,
    const mmcVarSet& Xi,
    const Switch initAtCnstr,
    bool readFields
)
:
    CloudType
    (
        cloudName,
        mesh,
        Xi,
        false,      // Only the top level cloud will call initAtCnstr
        readFields
    ),

    mixingPopeCloud(),

    cloudCopyPtr_(nullptr),

    mixingModel_(nullptr),

    secondCondR_(0.0),
    secondCondBeta_(1.0),
    secondCondTauOU_(1.0)

{
    // Read second conditioning parameters (see comment in first constructor)
    {
        const dictionary& scDict =
            this->cloudProperties_.subOrEmptyDict("secondConditioning");
        secondCondR_     = scDict.lookupOrDefault("R",     scalar(0.0));
        secondCondBeta_  = scDict.lookupOrDefault("beta",  scalar(1.0));
        secondCondTauOU_ = scDict.lookupOrDefault("tauOU", scalar(1.0));
    }

    Info << "Creating mixing Pope Particle Cloud." << nl << endl;

    setModels(Xi); // passing mmcVarSet

    Info << nl << "Mixing model constructed." << endl;

    setEulerianStatistics();

    if(readFields)
    {
        if(this->size()>0)
        {
            Info << nl << "Reading Mixing Pope particle cloud data from file." << endl;

            particleType::mixingParticleIOType::readFields(*this, this->mixing());
        }

        else
        {
            if(initAtCnstr)
            {
                Info << "Initial realease of Pope particles into the finite volume field." << nl << endl;

                this->initReleaseParticles();
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MixingPopeCloud<CloudType>::~MixingPopeCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::MixingPopeCloud<CloudType>::setParticleProperties
(
    particleType& particle,
    const scalar& mass,
    const scalar& wt,
    const scalar& patchI,
    const scalar& patchFace,
    const bool& iniRls
)
{
    CloudType::setParticleProperties(particle, mass, wt, patchI, patchFace,iniRls);    
    
    /*if 
    (
        particleType::indexInXiR_.empty()
     && particleType::XiRNames_.empty()
    )
    {
        particleType::indexInXiR_ = mixing().XiR().rVarInXiR();
        particleType::XiRNames_   = mixing().XiRNames();
    }*/


    particle.dx() = 0;

    // Assign second-conditioning subset flag based on fraction R
    particle.secondCondFlag() =
        (secondCondR_ > 0 && this->rndGen_.Random() < secondCondR_) ? 1 : 0;

    label numXiR = mixing().numXiR();

    //- Extension for a set of variables    
    particle.dXiR().setSize(numXiR,0.0);
    
    if (!iniRls)
    {
        particle.XiR() = mixing().XiR0(patchI,patchFace,particle);
    }
    else
    {
        particle.XiR() = mixing().XiR0(particle.cell(),particle);
    } 

    particle.initStatisticalSampling();
}


template<class CloudType>
void Foam::MixingPopeCloud<CloudType>::setEulerianStatistics()
{
    
    if (this->eulerianStatsDict().found("dx"))
    {
        const dimensionSet dim = dimLength;
        
        this->eulerianStats().newProperty("dx",dim);
    }
    
    //- statistics of reference variables (mixing distances)
    forAll(this->mixing().XiRNames(),XiRI)
    {
        word refVarName = "d" + mixing().XiRNames()[XiRI];
        
        if (this->eulerianStatsDict().found(refVarName))
        {
            const dimensionSet dim = dimless;// they will be considered dimless 
            
            this->eulerianStats().newProperty(refVarName,dim);   
        }
    }
}


template<class CloudType>
void Foam::MixingPopeCloud<CloudType>::updateEulerianStatistics()
{
    CloudType::updateEulerianStatistics();

    forAllIters(*this, iter)
    {
        this->eulerianStats().findCell(iter().position());

        //- statistics of distance for reference variables 
        forAll(mixing().XiRNames(),XiRI)
        {
            const word& refVarName = mixing().XiRNames()[XiRI];
            const word& dRefVarName = "d" + refVarName;
            
            if (this->eulerianStatsDict().found(dRefVarName))
                this->eulerianStats().calculate(dRefVarName,iter().wt(),iter().dXiR(refVarName));
        }
            
        if (this->eulerianStatsDict().found("dx"))
            this->eulerianStats().calculate("dx",iter().wt(),iter().dx());
    }
}


template<class CloudType>
void Foam::MixingPopeCloud<CloudType>::writeFields() const
{
    CloudType::writeFields();
    
    if (this->size())
    {
        particleType::mixingParticleIOType::writeFields(*this, this->mixing());
    }
}

// ************************************************************************* //
