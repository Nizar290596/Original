#!/bin/bash
cd "${0%/*}" || exit   # Run from this directory

# source support functions 
source ../../etc/supportFunctions.sh

basePath=${PWD}

__bannerByFilePath "$basePath/flameSheet"
cd flameSheet/ 
./Allrun > log.run
if [[ $? == 0 ]]; then 
    echo "    Success"
else
    echo "    Failure"
fi
