echo $SHELL
cd /afs/cern.ch/work/a/amartell/Linearity/Linearity/
pwd 
cd ../CMSSW_5_3_5/src
eval `scramv1 runtime -sh`
cd -
pwd 
source scripts/setup.sh 
root-config --version 
echo $PATH
unbuffer studyLinearity_MZ.exe /afs/cern.ch/work/a/amartell/Linearity/Linearity//data/2012/Winter2013/MZ_stdCat_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15/jobs//params_cat2.cfg >& /afs/cern.ch/work/a/amartell/Linearity/Linearity//data/2012/Winter2013/MZ_stdCat_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15/jobs//out_cat2.txt
