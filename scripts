#!/usr/bin/python

import os, sys

sys.path.append(os.path.abspath(os.path.curdir))

if __name__ == '__main__':
	energies = [5,10,20,30,40,60,80,100,150,200]
	eta = [1.7]
	identifier = "_mlsample_normalized"
	
	for i in energies:
		for j in eta:
			file1 = open("simpleBHet"+str(i)+"_eta"+str(j)+identifier+".cfg","w")
			file1.write("pNevts=0\n")
			file1.write("filePath=/data/users/snabili/gammaparticle/gamma/\n")
			file1.write("digifilePath=/data/users/snabili/gammaparticle/gamma/\n")
			file1.write("simFileName=HGcal__version63_model2_BON_et"+str(i)+"_eta"+str(j)+"00\n")
			file1.write("recoFileName=DigiIC3__version63_model2_BON_et"+str(i)+"_eta"+str(j)+"00\n")
			file1.write("outFilePath=/dev/null\n") #/data/users/chpapage/standalone/PFCal/PFCalEE/analysis/test.root")
			file1.write("debug=0\n")
			file1.write("nRuns=9\n")
			file1.write("deadfrac=0.01\n")
			file1.write("adjacent=1\n")
			file1.write("MLsample=1\n")
			file1.write("MLFilePath=/data/users/chpapage/standalone/PFCal/PFCalEE/analysis/training_sample_et"+str(i)+"_eta"+str(j)+"_2.root\n")
			file1.close()
