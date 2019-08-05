#!/bin/bash

export THORNADO_MACHINE=$1

if [[ $THORNADO_MACHINE == mac* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == titan* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

  source ${MODULESHOME}/init/bash

  module unload fftw cray-hdf5 cray-petsc silo subversion
  module unload pgi gcc cce pathscale
  module unload PrgEnv-pgi PrgEnv-gnu PrgEnv-cray PrgEnv-pathscale PrgEnv-intel

elif [[ $THORNADO_MACHINE == darter* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

  module unload fftw cray-hdf5 cray-petsc
  module unload pgi gcc cce intel
  module unload PrgEnv-pgi PrgEnv-gnu PrgEnv-cray PrgEnv-intel

elif [[ $THORNADO_MACHINE == beacon* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

  module unload pgi gcc cce intel
  module unload PE-gnu PE-intel

elif [[ $THORNADO_MACHINE == acf* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

  module unload PE-intel PE-gnu 

elif [[ $THORNADO_MACHINE == sn1987b* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == bbarker* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == juliana* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == kristopher* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == sjdunham* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == accre* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == ranchu* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == zelledge* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == ranchuair* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == mcarpe21* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == paullaiu_gnu* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

fi

if [[ $THORNADO_MACHINE == mac_gnu ]]; then

  echo
  export WEAKLIB_DIR=
  export POSEIDON_DIR=

elif [[ $THORNADO_MACHINE == titan_gnu ]]; then

  echo

  module load PrgEnv-gnu
  module load cray-hdf5
  module load subversion

elif [[ $THORNADO_MACHINE == titan_cray ]]; then

  echo

  module load PrgEnv-cray
  module load cray-hdf5
  module load subversion

elif [[ $THORNADO_MACHINE == darter_gnu ]]; then

  echo

  module load PrgEnv-gnu
  module load cray-hdf5

elif [[ $THORNADO_MACHINE == darter_cray ]]; then

  echo

  module load PrgEnv-cray
  module load cray-hdf5

elif [[ $THORNADO_MACHINE == beacon_intel ]]; then

  echo

  module load PE-intel
  module load hdf5/1.8.14

elif [[ $THORNADO_MACHINE == acf_gnu ]]; then

  echo 

  module load PE-gnu
  module swap intel-mpi/2017.2.174 openmpi/2.1.1-gcc4.8.2
  module load hdf5
  module load lapack 

elif [[ $THORNADO_MACHINE == sn1987b ]]; then

  echo

elif [[ $THORNADO_MACHINE == bbarker ]]; then

  echo

elif [[ $THORNADO_MACHINE == juliana ]]; then

  echo
  export WEAKLIB_DIR=~/Desktop/thornado_files/weaklib
  export THORNADO_DIR=~/Desktop/thornado_files/thornado_new

elif [[ $THORNADO_MACHINE == kristopher ]]; then

  echo
  export THORNADO_DIR=~/Desktop/Work/ORNL/thornado
  export POSEIDON_DIR=~/Desktop/Work/ORNL/poseidon

elif [[ $THORNADO_MACHINE == sjdunham ]]; then

  echo

elif [[ $THORNADO_MACHINE == accre ]]; then

  echo

elif [[ $THORNADO_MACHINE == ranchu ]]; then

  echo

elif [[ $THORNADO_MACHINE == zelledge ]]; then

  echo
  export WEAKLIB_DIR=/mnt/c/Users/Zack/weaklib/
  export THORNADO_DIR=/mnt/c/Users/Zack/thornado/

elif [[ $THORNADO_MACHINE == ranchuair ]]; then

  echo

elif [[ $THORNADO_MACHINE == mcarpe21 ]]; then

  echo
  export THORNADO_DIR=~/thornado
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/hdf5/lib/

elif [[ $THORNADO_MACHINE == paullaiu_gnu ]]; then

  echo
  export THORNADO_DIR=/home/plc/Documents/Git/thornado
  export WEAKLIB_DIR=/home/plc/Documents/Git/Neutrino/weaklib

else

  echo "  WARNING: Unknown machine " $THORNADO_MACHINE

fi
