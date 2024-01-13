# Oscillatory collision approach in the Earth-Moon restricted three body problem

## Overview
This repository contains the source code for the computer assisted proofs.

## Building and launching
### Prerequisites
In order to build and launch the program it is necessary to clone this repository into the PC with Linux OS (other operating systems are not supported). It is necessary to have autotools, cmake, and c++ compiler that supports C++17 standard installed.

### Launching the program

In order to build and launch the program the following commands have to be executed:

    git clone https://github.com/AleksanderPasiut/pcr3bp_code
    cd pcr3bp_code
    cmake .. -DCAPD_ENABLE_MULTIPRECISION=OFF
    make
    ./pcr3bp_code --gtest_filter=Pcr3bp_proof.*
