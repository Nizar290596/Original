#!/bin/bash
cd "${0%/*}" || exit   # Run from this directory

# Run the flameD test cases

# source support functions 
source ../../etc/supportFunctions.sh


basePath=${PWD}

__bannerByFilePath "$basePath/sectionalMethods/sectionalHomogeneousNucleationCondensation/jet/"
cd sectionalMethods/sectionalHomogeneousNucleationCondensation/jet/ 
./Allrun > log.run 2>&1
if [[ $? == 0 ]]; then 
    echo "    Success"
else
    echo "    Failure"
fi

__bannerByFilePath "$basePath/sectionalMethods/sectionalHomogeneousNucleationCondensation/psr"
cd $basePath
cd sectionalMethods/sectionalHomogeneousNucleationCondensation/psr 
./Allrun > log.run 2>&1
if [[ $? == 0 ]]; then 
    echo "    Success"
else
    echo "    Failure"
fi

__bannerByFilePath "$basePath/sectionalMethods/sectionalParticleFlameSynthesis/jet"
cd $basePath
cd sectionalMethods/sectionalParticleFlameSynthesis/jet
./Allrun > log.run 2>&1
if [[ $? == 0 ]]; then 
    echo "    Success"
else
    echo "    Failure"
fi

__bannerByFilePath "$basePath/sectionalMethods/sectionalParticleFlameSynthesis/psr"
cd $basePath
cd sectionalMethods/sectionalParticleFlameSynthesis/psr
./Allrun > log.run 2>&1
if [[ $? == 0 ]]; then 
    echo "    Success"
else
    echo "    Failure"
fi

