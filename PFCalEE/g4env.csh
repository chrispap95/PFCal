#!/bin/tcsh -f
#set echo
#set verbose
setenv USERBASE `pwd`
setenv ARCH x86_64-slc6-gcc46-opt
setenv ARCH2 x86_64-slc6-gcc46-opt
setenv ARCHd x86_64-slc6-gcc46-dbg
source /cvmfs/sft.cern.ch/lcg/external/gcc/4.6.3/x86_64-slc6/setup.csh
setenv QTHOME /cvmfs/sft.cern.ch/lcg/external/qt/4.8.4/${ARCH}/
setenv DAWNHOME /afs/cern.ch/sw/lcg/external/dawn/3_88a/x86_64-slc5-gcc43-opt/
setenv G4DAWNFILE_DEST_DIR ${USERBASE}/DawnFiles/
setenv HEPMC_DIR /cvmfs/sft.cern.ch/lcg/external/HepMC/2.06.08/${ARCH}/
setenv FASTJET_INSTALL /cvmfs/sft.cern.ch/lcg/external/fastjet/3.0.3/${ARCH}/
setenv BOOSTSYS /cvmfs/cms.cern.ch/slc6_amd64_gcc462/external/boost/1.47.0-cms

setenv G4BASE /cvmfs/geant4.cern.ch/geant4
cd $G4BASE/9.6.p02/${ARCHd}/share/Geant4-9.6.2/geant4make/
source geant4make.csh
cd ${USERBASE}

cd /cvmfs/sft.cern.ch/lcg/external/ROOT/5.34.00/${ARCH}/root/bin
source thisroot.csh
cd ${USERBASE}

setenv PYTHONDIR /cvmfs/sft.cern.ch/lcg/external/Python/2.7.3/${ARCH}/
setenv PYTHONPATH ${PYTHONDIR}:${ROOTSYS}/lib
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:$XERCESCROOT/lib:$HEPMC_DIR/lib:$USERBASE/userlib/lib:$USERBASE/analysis/lib:$FASTJET_INSTALL/lib:$BOOSTSYS/lib:$ROOTSYS/lib:$PYTHONDIR/lib
set path = ($DAWNHOME/bin $path $FASTJET_INSTALL/bin $BOOSTSYS)
#echo done
