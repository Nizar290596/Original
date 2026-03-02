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

#include "KernelEstimation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class CloudType>
void Foam::KernelEstimation<CloudType>::buildLESParticleList()
{
    const HashTable<label, word>& fieldIndexes = this->XiC().cVarInXi();
    const HashTable<label, word>& pfieldIndexes = this->XiC().cVarInXiC();

    LESList_.setSize(mesh_.cells().size());

    forAll(mesh_.cells(),celli)
    {
        densParticle p(numP_);

        p[nI_]  = celli;
        p[nx_]  = mesh_.C()[celli][0];
        p[ny_]  = mesh_.C()[celli][1];
        p[nz_]  = mesh_.C()[celli][2];
        p[nwt_] = 0.;
        p[nT_]  = 0.;

        //- Radial distance to origin
        p[nDist_] = mag(mesh_.C()[celli]);

        //- Species
        label n = nYs_;

        forAllIters(Yindexes_, Yi)
        {
            p[n] = 0.;//YEqvETarget[*Yi][celli];
            n++;
        }

        forAllConstIter(wordList,this->XiCNames(), nameI)
        {
            p[nXiCs_ +  pfieldIndexes[*nameI]] =
                this->XiC().Vars(fieldIndexes[*nameI]).field()[celli];
        }

        LESList_[celli] = p;
    }
}

template <class CloudType>
void Foam::KernelEstimation<CloudType>::buildParticleList()
{
    particleList_.clear();

    scalar nParticles = this->owner().size();
    scalar nLESCells  = mesh_.cells().size();
    scalar prob       = 1. - (0.5*nLESCells/nParticles);

    forAllIters(this->owner(), iter)
    {
        densParticle p(numP_);

        p[nI_] = 0.;

        p[nx_] = iter().position().x();

        p[ny_] = iter().position().y();

        p[nz_] = iter().position().z();

        p[nwt_] = iter().m();

        p[nT_] = iter().T();

        p[nDist_] = VGREAT;

        //- Species
        label n = 0;
        forAllIter(labelList,Yindexes_, Yi)
        {
            p[nYs_ + n] = iter().Y()[*Yi];

            //- Find min and max of particle for each YEqvE
            maxVal_[n] = max(maxVal_[n],iter().Y()[*Yi]);
            minVal_[n] = min(minVal_[n],iter().Y()[*Yi]);
            n++;
        }

        //- Coupling(state) Variables
        forAll(iter().XiC(), j)
        {
            p[nXiCs_ + j] = iter().XiC()[j];
            maxXiCVal_[j] = max(maxXiCVal_[j],iter().XiC()[j]);
            minXiCVal_[j] = min(minXiCVal_[j],iter().XiC()[j]);
        }

        //- Min and Max of temperature
        maxVal_[numYEqv_] = max(maxVal_[numYEqv_],iter().T());
        minVal_[numYEqv_] = min(minVal_[numYEqv_],iter().T());

        if(nParticles > nLESCells)
        {
            if (this->owner().rndGen().Random() >= prob)
                particleList_.append(p);
        }
        else
            particleList_.append(p);
    }
}

template <class CloudType>
void Foam::KernelEstimation<CloudType>::computeTargets
(
    const mmcVarSet& XiC,
    PtrList<volScalarField>& YEqvETarget,
    volScalarField& TEqvETarget
)
{
    label N2Index = 0;

    volScalarField Yt(0.0*YEqvETarget[0]);

    //- Indicator (used to make zero source terms)
    this->Indicator() = Yt;

    TEqvETarget = this->owner().T();

    scalarField dTEqvETarget(TEqvETarget.size(),0.0);

    List<scalarField> dYEqvETarget;

    forAll(YEqvETarget, specieI)
    {
        if (!this->solveEqvSpecie()[specieI])
            YEqvETarget[specieI].primitiveFieldRef() = 0.0;
        else
        {
            YEqvETarget[specieI] = this->owner().composition().carrier().Y()[specieI];
           dYEqvETarget.append(scalarField(TEqvETarget.size(),0.0));

            Yt += YEqvETarget[specieI];
        }

        if (YEqvETarget[specieI].name()==inertSpecie_)
            N2Index = specieI;
    }

    //- Index of coupling variable used
        const label CVIndexinXiC = nXiCs_ + XiC.cVarInXiC()[this->cVarName()];

    //- Min and Max values of the coupling variable
    const label CVIndexinXi = XiC.cVarInXi()[this->cVarName()];

    scalar minf = min(XiC.Vars(CVIndexinXi).field().primitiveField());
    scalar maxf = max(XiC.Vars(CVIndexinXi).field().primitiveField());

    scalar minz = minXiCVal_[XiC.cVarInXiC()[this->cVarName()]];
    scalar maxz = maxXiCVal_[XiC.cVarInXiC()[this->cVarName()]];

    if (debug_)
        Info << "max and min values of couplingVar in LES are: "
             << maxf << " and "<< minf
             << "\nmax and min values of couplingVar in particles are: "
             << maxz << " and "<< minz
             << "\nmax and min values of YEqv in particles are: \n"
             << maxVal_ << nl << minVal_ << endl;

    //- Pointer Lists
        //- LES Particles List

        List<densParticle*> LESPtrList(LESList_.size());

        forAll(LESList_,ii)
               LESPtrList[ii] = &LESList_[ii];

    //- Iterators for sorting
    List<densParticle*>::iterator iterL2 = LESPtrList.begin();
    List<densParticle*>::iterator iterU2 = LESPtrList.end();

    if (debug_)
        Info << "Iterators initialized" << endl;

    //- Sort particleList

        labelList dims(4);
        dims[0] = nx_;
        dims[1] = ny_;
        dims[2] = nz_;
        dims[3] = CVIndexinXiC; 

        scalarList wts(4);
        wts[0] = 1.;
        wts[1] = 1.;
        wts[2] = 1.;
        wts[3] = fm_;

        // Convert the particleList to a List<List<scalar>> with only the 
        // dimensions of dims for kdTree 
        // create the list
        List<List<scalar>> pList(particleList_.size());

        forAll(pList,i)
        {
            List<scalar> temp(4);
            for (int j=0; j < 4; ++j)
                temp[j] = particleList_[i][dims[j]];
            pList[i] =std::move(temp);
        }

        //- Construct kd-Tree of particle list
        kdTree<List<scalar>> particleTree(pList,wts); // 4 dimensions for kdtree x,y,z,XiC

        //Sort LES list based on distance from origin
        std::sort(iterL2,iterU2,lessArg(nDist_));

    if (debug_)
        Info << "Lists sorted" << endl;

    //- Field with cell centres
    vectorField cellCentres = mesh_.C().internalField();

    //- LES field of coupling variable
    const volScalarField& f = XiC.Vars(CVIndexinXi).field();//

    volVectorField gradf = fvc::grad(f);

    if (debug_)
        Info << "Start to compute cell averages" << endl;

    label computedLES = 0;// counter for kernel computed cells

    //- Index of last cell computed with kernel
    label lastCCell = (*LESPtrList.first())[nI_];

    scalar DELTAd = 0.0;
    scalar DELTAf = 0.0;

    //- Start computation of Cell target values
    forAllIter(DynamicList<densParticle*>,LESPtrList,LESi)
    {
        scalar fLES  = (**LESi)[CVIndexinXiC];

        label celli  = (**LESi)[nI_];

        scalar f0LES = f[lastCCell];
        scalar df    = fLES - f0LES;

        //- Do not compute cell value, source term is zero
        //- Used to avoid biased kernel at boundaries in coupling space
        //- It also helps to speed up computations
        if ((fLES <= fLow_) || (fLES >= fHigh_))
            continue;

        //- Cell centre
        vector cCentre  = cellCentres[celli];
        vector c0Centre = cellCentres[lastCCell];

        //- Physical distance
        scalar dd = mag(cCentre - c0Centre);

        //- compute using kernel if cell is far in coupling space (df>= DELTAf)
        //- or physical space (dd>= DELTAd)
        if (mag(df) >= DELTAf || dd >= DELTAd || computedLES == 0)
        {
            //- Find the k-nn particles to compute kernel !!!!!!!!
            label nn = 20;//50;
            scalarList qv(4); //coupling var + 3 dimensions

            qv[0] = cCentre[0];
            qv[1] = cCentre[1];
            qv[2] = cCentre[2];
            qv[3] = fLES;

            auto result = particleTree.nNearest
            (
                qv,             // query vector (x,y,z,XiC)
                nn              // number of nearest neighbours
            );

            //- Define kernel radius in physical space based on k-NN

            vector disVector;
            disVector[0] = result.begin()->disReal[0];
            disVector[1] = result.begin()->disReal[1];
            disVector[2] = result.begin()->disReal[2];

            scalar rMax = min(max(0.0,mag(disVector)),rMaxMax_);

            //- Start to compute target values
            scalar sumWt   = 0.0;
            scalar sumdWt  = 0.0;
            scalar sumWtT  = 0.0;
            scalar sumdWtT = 0.0;
            scalarField sumWtYi(Yindexes_.size(),0.0);
            scalarField sumdWtYi(Yindexes_.size(),0.0);

            //- Loop over k-NN Paricles

            scalar h  = 0.25*fm_;
            scalar h2 = 0.5*(rMax+SMALL);

            scalar alpha   = 1.0/6.0/h;
            scalar alpha2  = alpha/(0.5*h);

            for(label nni = 0; nni < nn; nni++)
            {
                const densParticle& p = particleList_(result[nni].idx);

                //- Distance in mix fraction
                scalar cpfDistance = -result[nni].disReal[3];

                scalar NfDistance = cpfDistance/h; //Normalized distance

                vector disXP;
                disXP[0] = result[nni].disReal[0];
                disXP[1] = result[nni].disReal[1];
                disXP[2] = result[nni].disReal[2];

                scalar NDistance  = mag(disXP)/h2; //////!!!!!!!!

                //- 1-D cubic smoothing kernel weight
                //- J J Monaghan Rep. Prog. Phys.68 (2005) 1703–1759
                //- doi:10.1088/0034-4885/68/8/R01

                scalar  IDWf   = 0.0;
                scalar dIDWf   = 0.0;

                scalar q    = mag(NfDistance);
                scalar sign = (NfDistance > 0.0) ? 1.0 : ((NfDistance < 0.0) ? -1.0 : 0.0);

                if (q >= 0.0 && q < 1.0)
                {
                    IDWf = alpha*(pow(2.0 - q, 3.0) - 4.0 * pow(1.0 - q, 3.0));
                   dIDWf = -sign*alpha2*(3.0*sqr(2.0 - q) - 12.0*sqr(1.0 - q));
                }
                else if (q >= 1.0 && q < 2.0)
                {
                    IDWf = alpha*(pow(2.0 - q, 3.0));
                   dIDWf = -sign*alpha2*(3.0*sqr(2.0 - q));
                }

                //- 3-D cubic smoothing kernel weight
                //- J J Monaghan Rep. Prog. Phys.68 (2005) 1703–1759
                //- doi:10.1088/0034-4885/68/8/R01

                scalar IDWd   = 0.0;
                scalar alphad = 1.0/(4.0*pi*pow(h2,3.0));

                if (NDistance < 1.0)
                    IDWd = alphad*(pow(2.0 - NDistance,3.0) - 4.0 * pow(1.0 - NDistance,3.0));
                else if ((NDistance >= 1.0) && (NDistance < 2.0))
                    IDWd = alphad*(pow(2.0 - NDistance,3.0));

                scalar mWti = IDWd * p[nwt_];
                scalar  Wti = mWti *  IDWf;
                scalar dWti = mWti * dIDWf;

                sumWt   = sumWt   +  Wti;
                sumdWt  = sumdWt  + dWti;

                sumWtT  = sumWtT  + (p[nT_]  *  Wti);
                sumdWtT = sumdWtT + (p[nT_]  * dWti);

                label ii = 0;
                for (label n = nYs_; n <= nYe_ ; n++)
                {
                    sumWtYi[ii]  = sumWtYi[ii]  + (p[n] *  Wti);
                    sumdWtYi[ii] = sumdWtYi[ii] + (p[n] * dWti);
                    ii++;
                }
            }

            if (sumWt < VSMALL)
                continue;

            Yt[celli]      = 0.0;
            label I        = 0;
            scalar maxGrad = 0.0;

            forAllIter(labelList, Yindexes_, specieI)
            {
                scalar Yi = sumWtYi[I]/sumWt;

                YEqvETarget[*specieI][celli] = Yi;
               dYEqvETarget[I][celli] = (sumdWtYi[I] - (Yi * sumdWt))/sumWt;

                //- Find maximum gradient to determine distance allowed to extrapolate
                scalar diffMag  = mag(maxVal_[I]-minVal_[I]);
                scalar Ngrad    = (diffMag > ROOTVSMALL) ? diffMag/h : 1e30;
                scalar maxGradi = mag(dYEqvETarget[I][celli]/Ngrad);

                if (maxGradi>maxGrad)
                    maxGrad = maxGradi;

                if (*specieI != N2Index)
                    Yt[celli] = Yt[celli] + YEqvETarget[*specieI][celli];

                I++;
            }

            TEqvETarget[celli]  = sumWtT/sumWt;
            dTEqvETarget[celli] = (sumdWtT - (TEqvETarget[celli] * sumdWt))/sumWt;

            //- Find maximum gradient to determine distance allowed to extrapolate
            scalar diffMag  = mag(maxVal_[I]-minVal_[I]);
            scalar Ngrad    = (diffMag > ROOTVSMALL) ? diffMag/h : 1e30;
            scalar maxGradi = mag(dTEqvETarget[celli]/Ngrad);

            if (maxGradi>maxGrad)
            maxGrad = maxGradi;

            DELTAf = C2_*h/(maxGrad+ROOTVSMALL);// define C2 !!! Needs checking

            if (maxGrad >= 10.0)
                DELTAf = 0.0;
            else
                DELTAf = min(dfMax_,DELTAf);// define DELTA_MAX

            scalar gradFi = cmptMax(cmptMag(gradf[celli]));

            DELTAd = DELTAf/(gradFi +ROOTVSMALL);

            this->Indicator()[celli] = 1.0;

            if (debug_ && (TEqvETarget[celli]>3000.0 || TEqvETarget[celli] < 290.0))
                Pout << "Target Temp is: " << TEqvETarget[celli] << endl;

            lastCCell = celli;

            computedLES++;
        }
        else
        {
        //- Compute target value using Linear Extrapolation from
        //- last cell computed
            Yt[celli] = 0.0;

            label I = 0;
            forAllIter(labelList,Yindexes_,specieI)
            {
                scalar Yi = YEqvETarget[*specieI][lastCCell];
                scalar deltaYi = dYEqvETarget[I][lastCCell]*df;

                YEqvETarget[*specieI][celli] = Yi + deltaYi;

                if (*specieI != N2Index)
                    Yt[celli] = Yt[celli] + YEqvETarget[*specieI][celli];

                I++;
            }

            scalar deltaT = dTEqvETarget[lastCCell]*df;

            TEqvETarget[celli] = TEqvETarget[lastCCell] + deltaT;

            this->Indicator()[celli] = 1.0;

            if (debug_ && (TEqvETarget[celli]>2200.0 || TEqvETarget[celli] < 290.0))
            {
                    Info << "fLES: " << fLES << endl;
                    Info << "Target Temp is: " << TEqvETarget[celli] << endl;
            }
        }
    }

    //- Make species add up to one
    YEqvETarget[N2Index] = scalar(1) - Yt;
    YEqvETarget[N2Index].max(0.0);

    if (debug_)
    {
        Pstream::gather(computedLES,sumOp<label>());
        Info << "Target Values Computed \nThe number of cells with a computed "
             << "target value using kernel estimator is: " << computedLES << endl;
    }
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::KernelEstimation<CloudType>::KernelEstimation
(
    const dictionary& dict,
    CloudType& owner,
    const mmcVarSet& Xi
)
:
    ThermoPhysicalCouplingModel<CloudType>(dict,owner,typeName,Xi),

    mesh_(owner.mesh()),

    Yindexes_(),

    numYEqv_(0),

    nI_(0),

    nx_(nI_ + 1),

    ny_(nx_ + 1),

    nz_(ny_ + 1),

    nwt_(nz_ + 1),

    nT_(nwt_ + 1),

    nDist_(nT_ + 1),

    nYs_(nDist_ + 1),

    nYe_(nYs_ + 1),

    nXiCs_(nYe_ + 1),

    nXiCe_(nXiCs_ + 1),

    numP_(0),

    fm_(readScalar(this->coeffDict().lookup("fm"))),

    rMaxMax_(this->coeffDict().lookupOrDefault("rMax", 1.0e9)),

    particleList_(),

    LESList_(),

    maxVal_(),

    minVal_(),

    maxXiCVal_(),

    minXiCVal_(),

    fLow_(readScalar(this->coeffDict().lookup("fLow"))),

    fHigh_(readScalar(this->coeffDict().lookup("fHigh"))),

    inertSpecie_(this->owner().thermo().lookup("inertSpecie")),

    dfMax_(readScalar(this->coeffDict().lookup("dfMax"))),

    C2_(readScalar(this->coeffDict().lookup("C2"))),

    debug_(this->coeffDict().lookupOrDefault("debug",false)) //Define from dictionary
{
    forAll(this->solveEqvSpecie(), i)
    {
        if (this->solveEqvSpecie()[i])
            Yindexes_.append(i);
    }

    numYEqv_ = Yindexes_.size();
    numXiC_  = this->XiCNames().size();

    nYe_    = nYs_ + numYEqv_ - 1;
    nXiCs_  = nYe_ + 1;
    nXiCe_  = nXiCs_ + numXiC_ - 1;
    numP_   = nXiCe_ + 1;

    maxVal_.setSize(numYEqv_+1,0.0);
    minVal_.setSize(numYEqv_+1,GREAT);

    maxXiCVal_.setSize(numXiC_,0.0);
    minXiCVal_.setSize(numXiC_,GREAT);

    if (debug_)
        Info << "The number of equivalent species is: " << numYEqv_
             << " and the number of properties is: " << numP_ << nl
             << "The indices of the equivalents species are: "
             << Yindexes_ << endl;
    Info << "maximum rMax is: " << rMaxMax_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::KernelEstimation<CloudType>::~KernelEstimation()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::KernelEstimation<CloudType>::EqvETargetValues
(
    const mmcVarSet& XiC,
    PtrList<volScalarField>& YEqvETarget,
    volScalarField& TEqvETarget
)
{
    //- Create lists
    buildParticleList();

    buildLESParticleList();

    //- Gather particles of neighbour processors to avoid bias at 
    //- processor boundaries.

    labelList pNeighbours;

    forAll(mesh_.boundaryMesh(),patchi)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& ppp = 
                refCast<const processorPolyPatch>(pp);

            pNeighbours.append(ppp.neighbProcNo());
        }
    }

    //- Use particles in processor plus neighbouring processors
    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        for (label procI : pNeighbours)
        {
            UOPstream toBuffer(procI, pBufs);
            toBuffer << particleList_;
        }
        
        pBufs.finishedSends();

        for (label procI : pNeighbours)
        {
            DynamicList<densParticle> tempList;
            UIPstream fromBuffer(procI, pBufs);
            fromBuffer >> tempList;
            // Add particles of processor to current list
            particleList_.append(tempList);
        }
    }

    computeTargets(XiC,YEqvETarget,TEqvETarget);

    LESList_.clearStorage();
    particleList_.clearStorage();
}

// ************************************************************************* //

