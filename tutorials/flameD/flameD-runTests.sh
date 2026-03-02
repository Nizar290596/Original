#!/bin/bash
cd "${0%/*}" || exit   # Run from this directory
# Run the flameD test cases

# source support functions 
source ../../etc/supportFunctions.sh

basePath=${PWD}

__bannerByFilePath "$basePath/GRI"
cd GRI/ && ./Allrun > log.run
if [[ $? == 0 ]]; then 
    echo "    Success"
else
    echo "    Failure"
fi

__bannerByFilePath "$basePath/DRM19"
cd $basePath
cd DRM19/ && ./Allrun > log.run
if [[ $? == 0 ]]; then 
    echo "    Success"
else
    echo "    Failure"
fi

__bannerByFilePath "$basePath/DRM19-TDAC"
cd $basePath
cd DRM19-TDAC/ && ./Allrun > log.run
if [[ $? == 0 ]]; then 
    echo "    Success"
else
    echo "    Failure"
fi

__bannerByFilePath "$basePath/NewtonLinODE"
cd $basePath
cd NewtonLinODE/ && ./Allrun > log.run
if [[ $? == 0 ]]; then 
    echo "    Success"
else
    echo "    Failure"
fi

__bannerByFilePath "$basePath/stdCurlsMixingModel"
cd $basePath
cd stdCurlsMixingModel/ && ./Allrun > log.run
if [[ $? == 0 ]]; then 
    echo "    Success"
else
    echo "    Failure"
fi

__bannerByFilePath "$basePath/RANS_DRM19"
cd $basePath
cd RANS_DRM19/ && ./Allrun > log.run
if [[ $? == 0 ]]; then 
    echo "    Success"
else
    echo "    Failure"
fi

