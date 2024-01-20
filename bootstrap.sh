#!/bin/bash

mkdir build
cd build
cmake .. -DCAPD_ENABLE_MULTIPRECISION=OFF

start=`date +%s%N`

        make -j 5

end=`date +%s%N`
let "build_time = ($end - $start) / 1000000"
echo Build time: $build_time ms.
