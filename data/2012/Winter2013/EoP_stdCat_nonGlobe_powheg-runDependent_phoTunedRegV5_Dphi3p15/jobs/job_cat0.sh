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
unbuffer studyLinearity_EoP.exe /afs/cern.ch/work/a/amartell/Linearity/Linearity//data/2012/Winter2013/EoP_stdCat_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15/jobs//params_cat0.cfg >& /afs/cern.ch/work/a/amartell/Linearity/Linearity//data/2012/Winter2013/EoP_stdCat_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15/jobs//out_cat0.txt
