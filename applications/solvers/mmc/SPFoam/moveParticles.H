/*---------------------------------------------------------------------------*\
|                               Move Particles                                |
|-----------------------------------------------------------------------------|
| Move particles and solve the reaction.                                      |
| Sample particle statistics if activated in the cloud dictionary             |
|                                                                             |
| Timestep sequence                                                           |
|   1. inflowBoundary().inflow()          particle injection                  |
|   2. mixing().Smix()                    first conditioning (all particles)  |
|   3. updatePhiReaction(dt)              W(phi) source (all particles)       |  <- guarded
|   4. updateOUProcess(dt)                OU advance + phi0 update (subset)   |  <- guarded
|   5. secondCondMixing().Smix()          second conditioning (subset)        |  <- guarded
|   6. solve()                            Lagrangian transport + chemistry    |
|   7. Sreact()                           load-balanced chemistry             |
|                                                                             |
| Steps 3-5 are only executed when secondCondMixingEnabled() is true.         |
\*---------------------------------------------------------------------------*/

if(pdfMethodOn)
{
    Info << "Evolve Particles -- Mixing, Moving, and Reaction" << endl;

    forAll(pSets, pSi)
    {
        basicReactingPopeParticle::trackingData td(pSets[pSi]);

        // 1. Particle injection at inflow boundaries
        pSets[pSi].inflowBoundary().inflow();

        // 2. First-conditioning mixing — all particles
        Info << "\tFirst-conditioning mixing..." << endl;
        const scalar t0mix = runTime.elapsedCpuTime();

        pSets[pSi].mixing().Smix();

        const scalar t1mix = runTime.elapsedCpuTime();
        Info << "\tFinish first-conditioning mixing in "
             << t1mix - t0mix << "s" << endl;

        // 3-5. Second-conditioning block — only when the model is enabled
        if (pSets[pSi].secondCondMixingEnabled())
        {
            const scalar deltaT = runTime.deltaTValue();

            // 3. W(phi) reaction for all particles
            pSets[pSi].updatePhiReaction(deltaT);

            // 4. Advance omega_OU and recompute phi0 for flagged particles
            pSets[pSi].updateOUProcess(deltaT);

            // 5. Second-conditioning mixing
            Info << "\tSecond-conditioning mixing..." << endl;
            const scalar t0sc = runTime.elapsedCpuTime();

            pSets[pSi].secondCondMixing().Smix();

            const scalar t1sc = runTime.elapsedCpuTime();
            Info << "\tFinish second-conditioning mixing in "
                 << t1sc - t0sc << "s" << endl;
        }

        // 6. Lagrangian transport + per-particle chemistry
        Info << "\tMove particle and solve reaction..." << endl;
        const scalar t0solve = runTime.elapsedCpuTime();

        pSets[pSi].solve(pSets[pSi], td);

        // 7. Load-balanced chemistry
        if(pSets[pSi].balanceReactionLoad())
        {
            scalar t0react = runTime.elapsedCpuTime();

            pSets[pSi].reaction().resetCPUTimeAnalysis();

            pSets[pSi].reaction().Sreact();

            pSets[pSi].reaction().writeCPUTimeAnalysis();

            pSets[pSi].reaction().ISATCntIncrmnt();

            scalar t1react = runTime.elapsedCpuTime();
            Info << "\tFinish balanced reaction in "
                 << t1react - t0react << "s" << endl;
        }
        const scalar t1solve = runTime.elapsedCpuTime();
        Info << "\tFinish solve reaction and move in "
             << t1solve - t0solve << "s" << endl;

        pSets[pSi].composition().correctMassFractions();

        pStats.particleSampling(pSets[pSi]);

        Info << "DONE" << endl;
    }

}
