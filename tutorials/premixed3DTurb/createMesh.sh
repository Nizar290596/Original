#!/bin/bash


# Check if folder for log files exists
if [ ! -d "logs" ]; then
    mkdir logs
fi


# Create meshes
printf "Creating and checking the meshes ====================\n"
starttime=$(date +%s)

printf "    LESMesh ............."
blockMesh > logs/blockMeshLESMesh.log
checkMesh > logs/checkMeshLESMesh.log
printf "... done\n"

printf "    evalLaplaceMesh .........." 
blockMesh -region evalLaplace > logs/blockMeshEvalLaplaceMesh.log
checkMesh -region evalLaplace > logs/checkMeshEvalLaplaceMesh.log
printf "... done\n"

printf "    superMesh ..........." 
blockMesh -region super > logs/blockMeshSuperMesh.log
checkMesh -region super > logs/checkMeshSuperMesh.log
printf "... done\n"

endtime=$(date +%s)
runtime=$((endtime-starttime))
printf "Done in $runtime seconds\n"
printf "=====================================================\n"


printf "\n"



