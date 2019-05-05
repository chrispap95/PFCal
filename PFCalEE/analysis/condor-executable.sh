#!/bin/bash

# requires 4 argument inputs:   
# 1: UNIQUE_ID - any unique string identifier  
# 2: CONDOR_PROCESS - condor process number  
# RUN_DIR - running directory (CMSSW_X_Y_Z/subdir)   
# mode.  should be 0 for background or 1 for signal

#
# header 
#

UNIQUE_ID=$1
CONDOR_PROCESS=$2
RUN_DIR=$3


START_TIME=`/bin/date`
echo "started at $START_TIME"


#
# setup CMSSW software environment at UMD
#
cd $RUN_DIR

cd ..
pwd
export USERBASE=`pwd`

#
# g4env.sh
#

export ARCH=x86_64-slc6-gcc46-opt
export ARCH2=x86_64-slc6-gcc46-opt
export ARCHd=x86_64-slc6-gcc46-dbg
source /cvmfs/sft.cern.ch/lcg/external/gcc/4.6.3/x86_64-slc6/setup.sh ""

locate gfortran

export QTHOME=/cvmfs/sft.cern.ch/lcg/external/qt/4.8.4/${ARCH}/
export DAWNHOME=/afs/cern.ch/sw/lcg/external/dawn/3_88a/x86_64-slc5-gcc43-opt/
export G4DAWNFILE_DEST_DIR=${USERBASE}/DawnFiles/
export HEPMC_DIR=/cvmfs/sft.cern.ch/lcg/external/HepMC/2.06.08/${ARCH}/
export FASTJET_INSTALL=/cvmfs/sft.cern.ch/lcg/external/fastjet/3.0.3/${ARCH}/
export BOOSTSYS=/cvmfs/cms.cern.ch/slc6_amd64_gcc462/external/boost/1.47.0-cms
export path=$path:/cvmfs/sft.cern.ch/lcg/external/CMake/2.8.6/${ARCH}/bin
export path=$path:/cvmfs/sft.cern.ch/lcg/external/expat/2.0.1/${ARCH}

export G4BASE=/data/users/data
cd /data/users/eno/geant4.9.6.p04-install/share/Geant4-9.6.4/geant4make/
source geant4make.sh
cd ${USERBASE}

cd /cvmfs/sft.cern.ch/lcg/external/ROOT/5.34.00/${ARCH}/root/bin
source thisroot.sh
cd ${USERBASE}

export XERCESCROOT=/cvmfs/sft.cern.ch/lcg/external/XercesC/3.1.1/${ARCH}/
export PYTHONDIR=/cvmfs/sft.cern.ch/lcg/external/Python/2.7.3/${ARCH}/
export PYTHONPATH=${PYTHONDIR}:${ROOTSYS}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$XERCESCROOT/lib:$HEPMC_DIR/lib:$USERBASE/userlib/lib:$USERBASE/analysis/lib:$FASTJET_INSTALL/lib:$BOOSTSYS/lib:$ROOTSYS/lib:$PYTHONDIR/lib
export path=$DAWNHOME/bin:$path:$FASTJET_INSTALL/bin:$BOOSTSYS

pwd
#
#
#


FINAL_PREFIX_NAME=`echo ${UNIQUE_ID}_${CONDOR_PROCESS}`
FINAL_LOG=`echo $FINAL_PREFIX_NAME.log`

#
# run c
#
cd $RUN_DIR
pwd
echo $LD_LIBRARY_PATH
locate gfortran
./bin/simpleBH -c/data/users/chpapage/standalone/PFCal/PFCalEE/analysis/scripts/simpleBHet$4_eta1.7_mlsample.cfg >> $FINAL_LOG 2>&1



#
# end run
#

echo ""
END_TIME=`/bin/date`
echo "finished at $END_TIME"
exit $exitcode
