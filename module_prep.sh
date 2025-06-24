#!/bin/bash

module purge
module load emacs
module load gcc-libs
module load gcc-libs/10.2.0
module load openblas/0.3.13-openmp/gnu-10.2.0

echo 'Modules successfully loaded:' 
module list 
