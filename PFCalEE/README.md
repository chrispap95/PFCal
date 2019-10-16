# PFCalEE

Geant4 simulation of a Si-base sampling calorimeter

Check https://twiki.cern.ch/twiki/bin/view/CMS/HGCalPerformanceStudiesWithG4

Geometry implementation is instantiated by detector versions in an enum - cf. src/DetectorConstruction.cc and src/SamplingSection.cc

A small ntuple is stored with the energy deposits - cf. src/EventAction.cc

When changing the TTree content, adding homemade classes rather than
simple objects, the classes should be added to userlib. The dictionary
for root to understand the classes also need to be remade. Use "make
dictionary" before make, inside of userlib/. Follow instructions in
https://twiki.cern.ch/twiki/bin/view/Sandbox/AnnemarieMagnanSandbox
for what needs to be put in the homemade classes and makefile for root to
understand them (see also example in class userlib/include/HGCSSSimHit.hh).

## Setup the environment
In lxplus, do:
```bash
source g4envLXPLUS.sh
```
In hepcms, do:
```bash
source g4env.csh
```
## Compile
```bash
mkdir -p userlib/{lib,obj,bin} && cd userlib && make dictionary && make -j 5 && cd - && make -j 5
```
Then, to run analysis code:
```bash
cd analysis && mkdir lib obj bin && make -j 5
```
Anytime you change a .cpp file, you have to redo:
```bash
make -j 5
```

## Creating samples
### Submit in parallel the runs submitProd.py using particle gun
The samples are made in two steps:
-Sim files use this command in PFCalEE:
```bash
for run in `seq 0 09`; do ./submitProdCondor.py -s 2nd -q 1nw -g -t V08-00-00 -r $run -v 63 -m 2 -a 1.7 -b 3.8 -d gamma -f "" -F "" -n 250 -o /eos/home-c/chpapage/gamma -e /eos/home-c/chpapage/gamma; done
```

-Digi files use this command in PFCalEE/userlib:
```bash
for run in `seq 0 09`; do ./submitDigiNoenergy.py -s 1nw -q 1nw -g -t testV8 -r $run -v 63 -m 2 -a 1.7 -b 3.8 -d gamma -n -1 -o /eos/home-c/chpapage/gamma -e /eos/home-c/chpapage/gamma -E /eos/home-c/chpapage/gamma; done
```

These will produce 10 sim and digi files. what you need to change are:
-q is the long queue. since you want to produce high stat samples you need to set long time
-t  use the git tag of the code you dowloaded. if you're using my code don't change it.
-a the eta of the projected sample
-p phi of sample
-n number of events per sample
-f and -F and -o change to your own directory

To set the energy of the sample you need to edit submitProd.py and submitDigi.py.
These scripts will produce a flat energy distribution in (0,570)GeV.
If you need discrete energies you need to uncomment the scripts and edit the PrimaryGeneratorAction.cc.
(Hint: see the parent repositories in this case!)
For loop is to generate several samples with same stat in parallel.

## Running simpleBH analysis
Go to analysis and do:
```bash
./simpleBH -c scripts/simpleBH.cfg
```
For the .cfg configuration files see the scripts/CfgGenerator.py script. It can generate multiple files in one run.
