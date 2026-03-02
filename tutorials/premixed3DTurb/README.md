# Case: Tin = 700; eqRatio = 0.6; turbulent

**Solver**: mmcPremixedFoam (OpenFOAM v2212)

**Additional**: inlet boundary condition for U field needs to be downloaded from git (ITV/syntheticTurbulenceInletUMeanControl)

**Modelling**: MMC-LES (LES with ATF; the density and progress variable source term are provided by the flamelet table (chemistry/flameletTab.dat))

0. (do not change the name of the starting folder to "0" because then the turbulence forcing is deactivated)
1. Install the inlet boundary condition for U field (link to a git repository is created in mmcFoam/thirdParty/syntheticTurbulenceInletUMeanControl)
2. To create all meshes run ./createMesh.sh
3. To decompose all meshes run ./decomposeMeshes.sh
4. I recommend to start the simulation, give the particle flame ~5ms to evolve and after that start sampling by changing a switch in constant/cloudProperties:
     
        particleSampling
        {
            enabled               true;
        ...
        }
      
5. In terminal run "gnuplot plotRelaxation" to track the progress variable relaxation over time
