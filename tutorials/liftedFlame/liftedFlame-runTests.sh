#!/bin/bash
cd "${0%/*}" || exit   # Run from this directory

# Run the flameD test cases

# source support functions 
source ../../etc/supportFunctions.sh


basePath=${PWD}

__bannerByFilePath "$basePath/Mueller"
cd Mueller
./Allrun > log.run 2>&1
if [[ $? == 0 ]]; then 
    echo "    Success"
else
    echo "    Failure"
fi

