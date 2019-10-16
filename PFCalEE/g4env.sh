#!/bin/bash
export USERBASE=`pwd`
ARCH=x86_64-slc6-gcc46-opt
source /cvmfs/sft.cern.ch/lcg/external/gcc/4.6.3/x86_64-slc6/setup.sh
export QTHOME=/cvmfs/sft.cern.ch/lcg/external/qt/4.8.4/${ARCH}/
export DAWNHOME=/afs/cern.ch/sw/lcg/external/dawn/3_88a/x86_64-slc5-gcc43-opt/
export G4DAWNFILE_DEST_DIR=${USERBASE}/DawnFiles/
export HEPMC_DIR=/cvmfs/sft.cern.ch/lcg/external/HepMC/2.06.08/${ARCH}/
export FASTJET_INSTALL=/cvmfs/sft.cern.ch/lcg/external/fastjet/3.0.3/${ARCH}/

export G4BASE=/data/users/eno
cd /data/users/eno/geant4.9.6.p04-install/share/Geant4-9.6.4/geant4make/
source geant4make.sh
cd ${USERBASE} &> /dev/null
#for boost latest version
export BOOSTSYS=/cvmfs/cms.cern.ch/slc6_amd64_gcc462/external/boost/1.47.0-cms
export XERCESCROOT=/cvmfs/sft.cern.ch/lcg/external/XercesC/3.1.1/${ARCH}/
export PYTHONDIR=/cvmfs/sft.cern.ch/lcg/external/Python/2.7.3/${ARCH}/
export PYTHONPATH=${PYTHONDIR}:${ROOTSYS}/lib

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$XERCESCROOT/lib:$HEPMC_DIR/lib:$USERBASE/userlib/lib:$USERBASE/analysis/lib:$FASTJET_INSTALL/lib:$BOOSTSYS/lib:$ROOTSYS/lib:$PYTHONDIR/lib
cd /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.18/${ARCH}/root/
source bin/thisroot.sh
cd - &> /dev/null
export PATH=$DAWNHOME/bin:$PATH:$FASTJET_INSTALL/bin:$BOOSTSYS:/cvmfs/sft.cern.ch/lcg/external/CMake/2.8.6/${ARCH}/bin:/cvmfs/sft.cern.ch/lcg/external/CMake/2.8.6/${ARCH}
