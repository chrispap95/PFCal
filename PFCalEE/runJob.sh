#!/bin/bash
source /data/users/eno/StandAlone/PFCal/PFCalEE/g4env.sh
eos cp /eos/cms/store/cmst3/group/hgcal/HGCalMinbias/Pythia8//data/example_MyPythia.dat data/example_MyPythia.dat
cp /data/users/eno/StandAlone/PFCal/PFCalEEgit_V08-01-00/version_63/model_2/pi-/BON/et_10/eta_1.700/g4steer.mac .
PFCalEE g4steer.mac 63 2 1.700000 1 1.75,1.75,1.75,1.75,1.75,2.8,2.8,2.8,2.8,2.8,4.2,4.2,4.2,4.2,4.2 1,1,1,1,1,2.1,2.1,2.1,2.1,2.1,4.4,4.4,4.4,4.4  | tee g4.log
mv PFcal.root HGcal__version63_model2_BON_et10_eta1.700.root
localdir=`pwd`
echo "--Local directory is " $localdir >> g4.log
ls -ltrh * >> g4.log
eos mkdir -p /eos/cms/data/users/eno/StandAlone/PFCal/PFCalEE/gitV08-01-00/pi-
eos cp HGcal__version63_model2_BON_et10_eta1.700.root /eos/cms/data/users/eno/StandAlone/PFCal/PFCalEE/gitV08-01-00/pi-/HGcal__version63_model2_BON_et10_eta1.700.root
if (( "$?" != "0" )); then
echo " --- Problem with copy of file PFcal.root to EOS. Keeping locally." >> g4.log
else
eossize=`eos ls -l /eos/cms/data/users/eno/StandAlone/PFCal/PFCalEE/gitV08-01-00/pi-/HGcal__version63_model2_BON_et10_eta1.700.root | awk '{print $5}'`
localsize=`ls -l HGcal__version63_model2_BON_et10_eta1.700.root | awk '{print $5}'`
if [ $eossize != $localsize ]; then
echo " --- Copy of sim file to eos failed. Localsize = $localsize, eossize = $eossize. Keeping locally..." >> g4.log
else
echo " --- Size check done: Localsize = $localsize, eossize = $eossize" >> g4.log
echo " --- File PFcal.root successfully copied to EOS: /eos/cms/data/users/eno/StandAlone/PFCal/PFCalEE/gitV08-01-00/pi-/HGcal__version63_model2_BON_et10_eta1.700.root" >> g4.log
rm HGcal__version63_model2_BON_et10_eta1.700.root
fi
fi
echo "--deleting core files and hepmc files: too heavy!!"
rm core.*
rm data/example_MyPythia.dat
cp * /data/users/eno/StandAlone/PFCal/PFCalEEgit_V08-01-00/version_63/model_2/pi-/BON/et_10/eta_1.700/
echo "All done"
