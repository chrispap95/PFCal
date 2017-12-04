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
OPT_MODE1=$4

START_TIME=`/bin/date`
echo "started at $START_TIME"


#
# setup CMSSW software environment at UMD
#
#export VO_CMS_SW_DIR=/sharesoft/cmssw
#. $VO_CMS_SW_DIR/cmsset_default.sh
cd $RUN_DIR

cd ..
pwd
export USERBASE=`pwd`

ARCH=x86_64-slc6-gcc46-opt
ARCH2=x86_64-slc6-gcc46-opt
ARCHd=x86_64-slc6-gcc46-dbg
gcc_config_version=4.6.3
mpfr_config_version=2.4.2
gmp_config_version=4.3.2
LCGPLAT=x86_64-slc6
LCG_lib_name=lib64
LCG_arch=x86_64
LCG_contdir=/cvmfs/sft.cern.ch/lcg/external
LCG_gcc_home=${LCG_contdir}/gcc/${gcc_config_version}/${LCGPLAT}
LCG_mpfr_home=${LCG_contdir}/mpfr/${mpfr_config_version}/${LCGPLAT}
LCG_gmp_home=${LCG_contdir}/gmp/${gmp_config_version}/${LCGPLAT}
export PATH=${LCG_gcc_home}/bin:${PATH}
export COMPILER_PATH=${LCG_gcc_home}/lib/gcc/${LCG_arch}-unknown-linux-gnu/${gcc_config_version}
export LD_LIBRARY_PATH=${LCG_gcc_home}/${LCG_lib_name}:${LCG_mpfr_home}/lib:${LCG_gmp_home}/lib:${LD_LIBRARY_PATH}


echo $LCG_contdir
export QTHOME=/cvmfs/sft.cern.ch/lcg/external/qt/4.8.4/${ARCH}/
export DAWNHOME=/afs/cern.ch/sw/lcg/external/dawn/3_88a/x86_64-slc5-gcc43-opt/
export G4DAWNFILE_DEST_DIR=${USERBASE}/DawnFiles/
export HEPMC_DIR=/cvmfs/sft.cern.ch/lcg/external/HepMC/2.06.08/${ARCH}/
export FASTJET_INSTALL=/cvmfs/sft.cern.ch/lcg/external/fastjet/3.0.3/${ARCH}/
export BOOSTSYS=/cvmfs/cms.cern.ch/slc6_amd64_gcc462/external/boost/1.47.0-cms

export G4BASE=/cvmfs/geant4.cern.ch/geant4
cd $G4BASE/9.6.p02/${ARCHd}/share/Geant4-9.6.2/geant4make/
source geant4make.sh
cd ${USERBASE}

cd /cvmfs/sft.cern.ch/lcg/external/ROOT/5.34.00/${ARCH}/root/bin
source thisroot.sh
cd ${USERBASE}

export PYTHONDIR=/cvmfs/sft.cern.ch/lcg/external/Python/2.7.3/${ARCH}/
export  PYTHONPATH=${PYTHONDIR}:${ROOTSYS}/lib
ls $USERBASE/userlib/lib
ls $USERBASE/analysis/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$XERCESCROOT/lib:$HEPMC_DIR/lib:$USERBASE/userlib/lib:$USERBASE/analysis/lib:$FASTJET_INSTALL/lib:$BOOSTSYS/lib:$ROOTSYS/lib:$PYTHONDIR/lib
echo $LD_LIBRARY_PATH
export path=$DAWNHOME/bin:$path:$FASTJET_INSTALL/bin:$BOOSTSYS




FINAL_PREFIX_NAME=`echo ${UNIQUE_ID}_${CONDOR_PROCESS}`
FINAL_LOG=`echo $FINAL_PREFIX_NAME.log`

#
# run c
#
cd $RUN_DIR
pwd
./bin/simpleBH $OPT_MODE1 >> $FINAL_LOG 2>&1



#
# end run
#

echo ""
END_TIME=`/bin/date`
echo "finished at $END_TIME"
exit $exitcode
