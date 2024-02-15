# Oscillatory collision approach in the Earth-Moon restricted three body problem

## Overview
This repository contains the source code for the interval arithmetic validation used in the paper "Oscillatory collision approach in the Earth-Moon restricted three body problem" by M. J. Capiński and A. Pasiut. 

The following source files contain the proofs for the respective theorems and lemmas:

* periodic_orbit_parameters_test.cpp
  * Theorem 7
* covering_relations_test.cpp
  * Theorem 8
  * Lemma 8
  * Lemma 11
* covering_relations_test.parallelogram_covering_derivative_check.cpp
  * Lemma 10

In the source files we use the term "periodic" to refer to the "ejection/collision" orbit. Also, we use the term "homoclinic" to refer to the "outer" orbit.

## Building and launching
### Prerequisites
In order to build and launch the program it is necessary to clone this repository into the PC with Linux OS (other operating systems are not supported). It is necessary to have cmake, git, and c++ compiler that supports C++17 standard installed.

### Launching the program

In order to build and launch the program the following commands have to be executed:

    git clone https://github.com/AleksanderPasiut/pcr3bp_code
    cd pcr3bp_code
    bash build_and_run.sh
