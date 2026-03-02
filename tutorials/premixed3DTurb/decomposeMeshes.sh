#!/bin/bash


# Check if folder for log files exists
if [ ! -d "logs" ]; then
    mkdir logs
fi

# Decompose domain
if [ -e system/decomposeParDict ]; then

    printf "Decomposing the meshes ==============================\n"
    starttime=$(date +%s)

    printf "    LESMesh ............."
    decomposePar > logs/decomposeParDNSMesh.log
    printf "... done\n"
    
    printf "    evalLaplace .........."  
    decomposePar -region evalLaplace > logs/decomposeParEvalLaplace.log
    printf "... done\n"

    endtime=$(date +%s)
    runtime=$((endtime-starttime))
    printf "Done in $runtime seconds\n"
    printf "=====================================================\n"

fi



