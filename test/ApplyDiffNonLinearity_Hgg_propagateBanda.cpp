//g++ -Wall -o ApplyDiffNonLinearity_Hgg_propagateBanda `root-config --cflags --glibs` ApplyDiffNonLinearity_Hgg_propagateBanda.cpp

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TMath.h"
#include "TPaveStats.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TTree.h"
#include "TChain.h"
#include "TVirtualFitter.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>


// Cholesky decomposition of a matrix
double* cholesky(double* A, int n);

int GetSingleCategory(const float& scEta1, const float& R91);

int GetMVAcat_7TeV(const float& dipho_mva);
int GetMVAcat_8TeV(const float& dipho_mva);

double BandaCat0(double* x, double* par)
{
  if(x[0] <= 35.)                    return (1.00012 + par[0]* 0.000581163);
  if(35. < x[0] && x[0] <= 40.)      return (1.00071 + par[0]* 0.000523539);
  if(40. < x[0] && x[0] <= 45.)      return (1.00027 + par[0]* 0.000211144);
  if(45. < x[0] && x[0] <= 50.)      return (0.999811 + par[0]* 0.000180878);
  if(50. < x[0] && x[0] <= 60.)      return (1.00029 + par[0]* 0.000414001);
  if(60. < x[0] && x[0] <= 75.)      return (1.0002 + par[0]* 0.000532232);
  if(75. < x[0] && x[0] <= 100.)     return (1.00049 + par[0]* 0.000897605);
  if(x[0] > 100.)                    return (1.00231 + par[0]* 0.00176987);
  //  if(60. < x[0] && x[0] <= 75.)      return(1.0002 + par[0]* 0.000564709);
  //   if(x[0] > 75.) return(1.00118 + par[0]* 0.00109502);
  return 1.;
}

double BandaCat1(double* x, double* par)
{
  if(x[0] <= 35.)                    return (1.00089 + par[0]* 0.000752811);
  if(35. < x[0] && x[0] <= 40.)      return (1.00065 + par[0]* 0.000524134);
  if(40. < x[0] && x[0] <= 45.)      return (0.999998 + par[0]* 6.68271e-05);
  if(45. < x[0] && x[0] <= 50.)      return (0.999691 + par[0]* 0.000273031);
  if(50. < x[0] && x[0] <= 60.)      return (0.999164 + par[0]* 0.000759375);
  if(60. < x[0] && x[0] <= 75.)      return (0.999977 + par[0]* 0.000430609);
  if(75. < x[0] && x[0] <= 100.)     return (0.99935 + par[0]* 0.00104642);
  if(x[0] > 100.)                    return (1.00156 + par[0]* 0.00199785);
  //   if(60. < x[0] && x[0] <= 75.) return(0.999977 + par[0]* 0.000427457);
  //   if(x[0] > 75.) return(0.999892 + par[0]* 0.000602402);
  return 1.;
}

double BandaCat2(double* x, double* par)
{
  if(x[0] <= 35.)                    return (0.991248 + par[0]* 0.00737642);
  if(35. < x[0] && x[0] <= 40.)      return (0.99888 + par[0]* 0.00156878);
  if(40. < x[0] && x[0] <= 45.)      return (0.999493 + par[0]* 0.000429543);
  if(45. < x[0] && x[0] <= 50.)      return (0.999855 + par[0]* 0.000274104);
  if(50. < x[0] && x[0] <= 60.)      return (0.997784 + par[0]* 0.00167762);
  if(60. < x[0] && x[0] <= 75.)      return (1.00304 + par[0]* 0.00260998);
  if(75. < x[0] && x[0] <= 100.)     return (0.997809 + par[0]* 0.00263123);
  if(x[0] > 100.)                    return (1.004 + par[0]* 0.00425648);
  //   if(60. < x[0] && x[0] <= 75.) return(1.00304 + par[0]* 0.0025946);
  //   if(x[0] > 75.) return(0.999712 + par[0]* 0.00147001);
  return 1.;
}

double BandaCat3(double* x, double* par)
{
  if(x[0] <= 35.)                    return (0.994715 + par[0]* 0.00455526);
  if(35. < x[0] && x[0] <= 40.)      return (1.00219 + par[0]* 0.00164928);
  if(40. < x[0] && x[0] <= 45.)      return (1.00036 + par[0]* 0.000328358);
  if(45. < x[0] && x[0] <= 50.)      return (1.00005 + par[0]* 0.000160501);
  if(50. < x[0] && x[0] <= 60.)      return (1.00028 + par[0]* 0.000636193);
  if(60. < x[0] && x[0] <= 75.)      return (1.00305 + par[0]* 0.00251856);
  if(75. < x[0] && x[0] <= 100.)     return (1.00111 + par[0]* 0.00220828);
  if(x[0] > 100.)                    return (1.00799 + par[0]* 0.00676785);
  //   if(60. < x[0] && x[0] <= 75.) return(1.00305 + par[0]* 0.0026342);
  //   if(x[0] > 75.) return(1.0027 + par[0]* 0.00288338);
  return 1.;
}

  

int main()
{
  int nSamples = 4;
  int nCats = 4;
  int nPar = 2;
  int nTrials = 1000;
  
  //  std::string methodName = "exp3par";  
  std::string methodName = "pol1";  
  //  std::string TF1_folder = "../TF1_"+methodName+"_MZ";
  std::string TF1_folder = "../TF1_"+methodName+"_EtScale";
  std::string TF1_folderSys = "../TF1_"+methodName+"_syst";

  //  std::string HggCatType = "MVA_7TeV";
  //  std::string HggCatType = "MVA_8TeV";
  //std::string HggCatType = "CiC_7TeV";
  std::string HggCatType = "CiC_8TeV";

  //std::string analysis  = "CiC"; 
  std::string analysis = "stdCat";
  
  bool checkSystematic = false;
  //////////////// da sotto qui ok 
  //////////////// cambia sopra in caso

  std::string Energy;

  if(HggCatType == "CiC_7TeV" || HggCatType == "MVA_7TeV") Energy = "7TeV";
  if(HggCatType == "CiC_8TeV" || HggCatType == "MVA_8TeV") Energy = "8TeV";
  
  std::cout << "analysis = " << analysis << std::endl;
  std::cout << "HggCatType = " << HggCatType << std::endl;

  double extHtBinEdges[8] = {60.,70.,80.,90.,100.,120.,150.,1000.};
  int nHT = 8;
  std::vector<float> HtBinEdges;
  for(int posVec = 0; posVec< nHT; ++posVec){
      HtBinEdges.push_back( extHtBinEdges[posVec]);
      std::cout << " HtBinEdges.at("<<posVec<<") = " << HtBinEdges.at(posVec) << std::endl;
  }
  std::cout << " OK done " << std::endl;

  
  //----------------
  // input functions
  TFile* funcFile;
  
  TF1** f_std = new TF1*[nCats];
  for(int iCat = 0; iCat < nCats; ++iCat)
  {
    if(iCat == 0)      f_std[iCat] = new TF1(Form("f_std_cat%d_std",iCat), BandaCat0,20., 1000., 1);
    if(iCat == 1)      f_std[iCat] = new TF1(Form("f_std_cat%d_std",iCat), BandaCat1,20., 1000., 1);
    if(iCat == 2)      f_std[iCat] = new TF1(Form("f_std_cat%d_std",iCat), BandaCat2,20., 1000., 1);
    if(iCat == 3)      f_std[iCat] = new TF1(Form("f_std_cat%d_std",iCat), BandaCat3,20., 1000., 1);
    double* parameter = new double[1];
    parameter[0] = 0.;
    f_std[iCat]->SetParameters(parameter);
  }
  TF1** f_stat = new TF1*[nCats];
  for(int iCat = 0; iCat < nCats; ++iCat)
  {
    funcFile = TFile::Open((TF1_folder+"/TF1_"+methodName+"_"+analysis+"_"+Energy+"_std_"+Form("cat%d.root",iCat)).c_str(),"READ");
    f_stat[iCat] = (TF1*)( funcFile->Get(Form("HT_cat%d",iCat)) );
    f_stat[iCat] -> SetName(Form("f_stat_cat%d_std",iCat));
  }
  std::cout << " Stat presa " << std::endl;  

  double** corMatrix = new double*[nCats];
  if(methodName == "exp"){
    nPar = 2;
  corMatrix[0] = new double[nPar*nPar];
  (corMatrix[0])[0] = 1;
  (corMatrix[0])[1] = -0.993354;
  (corMatrix[0])[2] = -0.993354;
  (corMatrix[0])[3] = 1;
  corMatrix[1] = new double[nPar*nPar];
  (corMatrix[1])[0] = 1.0000e+00;
  (corMatrix[1])[1] = -9.9800e-01;
  (corMatrix[1])[2] = -9.9800e-01;
  (corMatrix[1])[3] = 1.0000e+00;
  corMatrix[2] = new double[nPar*nPar];
  (corMatrix[2])[0] = 1.0000e+00;
  (corMatrix[2])[1] = -9.9801e-01;
  (corMatrix[2])[2] = -9.9801e-01;
  (corMatrix[2])[3] = 1.0000e+00;
  corMatrix[3] = new double[nPar*nPar];
  (corMatrix[3])[0] = 1.0000e+00;
  (corMatrix[3])[1] = -9.9800e-01;
  (corMatrix[3])[2] = -9.9800e-01;
  (corMatrix[3])[3] = 1.0000e+00;
  }
  if(methodName == "exp3par"){
    nPar = 3;
  corMatrix[0] = new double[nPar*nPar];
  (corMatrix[0])[0] = 1;
  (corMatrix[0])[1] = -0.994128;
  (corMatrix[0])[2] = 0.368711;
  (corMatrix[0])[3] = -0.994128;
  (corMatrix[0])[4] = 1;
  (corMatrix[0])[5] = -0.339976;
  (corMatrix[0])[6] = 0.368711;
  (corMatrix[0])[7] = -0.339976;
  (corMatrix[0])[8] = 1;
  corMatrix[1] = new double[nPar*nPar];
  (corMatrix[1])[0] = 1.0000e+00;
  (corMatrix[1])[1] = -9.8999e-01;
  (corMatrix[1])[2] = -4.5653e-01;
  (corMatrix[1])[3] = -9.8999e-01;
  (corMatrix[1])[4] = 1.0000e+00;
  (corMatrix[1])[5] = 4.8328e-01;
  (corMatrix[1])[6] = -4.5653e-01;
  (corMatrix[1])[7] = 4.8328e-01;
  (corMatrix[1])[8] = 1.0000e+00;
  corMatrix[2] = new double[nPar*nPar];
  (corMatrix[2])[0] = 1.0000e+00;
  (corMatrix[2])[1] = -9.9900e-01;
  (corMatrix[2])[2] = 1.3125e-09;
  (corMatrix[2])[3] = -9.9900e-01;
  (corMatrix[2])[4] = 1.0000e+00;
  (corMatrix[2])[5] = 1.2810e-09;
  (corMatrix[2])[6] = 1.3125e-09;
  (corMatrix[2])[7] = 1.2810e-09;
  (corMatrix[2])[8] = 1.0000e+00;
  corMatrix[3] = new double[nPar*nPar];
  (corMatrix[3])[0] = 1.0000e+00;
  (corMatrix[3])[1] = -9.9737e-01;
  (corMatrix[3])[2] = 1.1856e-02;
  (corMatrix[3])[3] = -9.9737e-01;
  (corMatrix[3])[4] = 1.0000e+00;
  (corMatrix[3])[5] = 1.6137e-02;
  (corMatrix[3])[6] = 1.1856e-02;
  (corMatrix[3])[7] = 1.6137e-02;
  (corMatrix[3])[8] = 1.0000e+00;
  }
  if(methodName == "pol1" && analysis == "CiC"){
    nPar = 2;
  corMatrix[0] = new double[nPar*nPar];
  (corMatrix[0])[0] = 1;
  (corMatrix[0])[1] = -0.0680281;
  (corMatrix[0])[2] = -0.0680281;
  (corMatrix[0])[3] = 1;
  corMatrix[1] = new double[nPar*nPar];
  (corMatrix[1])[0] = 1.0000e+00;
  (corMatrix[1])[1] = 0.1428;
  (corMatrix[1])[2] = 0.1428;
  (corMatrix[1])[3] = 1.0000e+00;
  corMatrix[2] = new double[nPar*nPar];
  (corMatrix[2])[0] = 1.0000e+00;
  (corMatrix[2])[1] = 0.2883;
  (corMatrix[2])[2] = 0.2883;
  (corMatrix[2])[3] = 1.0000e+00;
  corMatrix[3] = new double[nPar*nPar];
  (corMatrix[3])[0] = 1.0000e+00;
  (corMatrix[3])[1] = 0.3860;
  (corMatrix[3])[2] = 0.3860;
  (corMatrix[3])[3] = 1.0000e+00;
  }
  if(methodName == "pol1" && analysis == "stdCat" && (HggCatType == "MVA_8TeV" || HggCatType == "CiC_8TeV")){
    nPar = 2;
    corMatrix[0] = new double[nPar*nPar];
  (corMatrix[0])[0] = 1;
  (corMatrix[0])[1] = 0.195665;
  (corMatrix[0])[2] = 0.195665;
  (corMatrix[0])[3] = 1;
  corMatrix[1] = new double[nPar*nPar];
  (corMatrix[1])[0] = 1.0000;
  (corMatrix[1])[1] = 0.3196;
  (corMatrix[1])[2] = 0.3196;
  (corMatrix[1])[3] = 1.0000;
  corMatrix[2] = new double[nPar*nPar];
  (corMatrix[2])[0] = 1.0000;
  (corMatrix[2])[1] = -0.2715;
  (corMatrix[2])[2] = -0.2715;
  (corMatrix[2])[3] = 1.0000;
  corMatrix[3] = new double[nPar*nPar];
  (corMatrix[3])[0] = 1.0000;
  (corMatrix[3])[1] = -0.0063;
  (corMatrix[3])[2] = -0.0063;
  (corMatrix[3])[3] = 1.0000;
  }

  if(methodName == "pol1" && analysis == "stdCat" && (HggCatType == "CiC_7TeV" || HggCatType == "MVA_7TeV")){
    nPar = 2;
  corMatrix[0] = new double[nPar*nPar];
  (corMatrix[0])[0] = 1;
  (corMatrix[0])[1] = -0.190282;
  (corMatrix[0])[2] = -0.190282;
  (corMatrix[0])[3] = 1;
  corMatrix[1] = new double[nPar*nPar];
  (corMatrix[1])[0] = 1.0000;
  (corMatrix[1])[1] = 0.1473;
  (corMatrix[1])[2] = 0.1473;
  (corMatrix[1])[3] = 1.0000;
  corMatrix[2] = new double[nPar*nPar];
  (corMatrix[2])[0] = 1.0000;
  (corMatrix[2])[1] = -0.2385;
  (corMatrix[2])[2] = -0.2385;
  (corMatrix[2])[3] = 1.0000;
  corMatrix[3] = new double[nPar*nPar];
  (corMatrix[3])[0] = 1.0000;
  (corMatrix[3])[1] = -0.3937;
  (corMatrix[3])[2] = -0.3937;
  (corMatrix[3])[3] = 1.0000;
  }

    
  double** C = new double*[nCats];
  for(int iCat = 0; iCat < nCats; ++iCat)
    C[iCat] = cholesky(corMatrix[iCat],nPar);

  std::cout << " Now Covarianze fatte " << std::endl;
  
  
  //------------
  // input trees
  TFile** inFile = new TFile*[nSamples];
  TTree** trees = new TTree*[nSamples];

  if(HggCatType == "CiC_8TeV"){
    inFile[0] = new TFile("/afs/cern.ch/user/m/musella/public/higgs/forArabella/histograms_CMS-HGG.root", "read");
    inFile[1] = new TFile("/afs/cern.ch/user/m/musella/public/higgs/forArabella/histograms_CMS-HGG.root", "read");
    inFile[2] = new TFile("/afs/cern.ch/user/m/musella/public/higgs/forArabella/histograms_CMS-HGG.root", "read");
    inFile[3] = new TFile("/afs/cern.ch/user/m/musella/public/higgs/forArabella/histograms_CMS-HGG.root", "read");
    trees[0] = (TTree*)inFile[0]->Get("ggh_m125_8TeV");
    trees[1] = (TTree*)inFile[0]->Get("vbf_m125_8TeV");
    trees[2] = (TTree*)inFile[0]->Get("wzh_m125_8TeV");
    trees[3] = (TTree*)inFile[0]->Get("tth_m125_8TeV");
  }
  if(HggCatType == "MVA_8TeV"){
    inFile[0] = new TFile("/afs/cern.ch/user/m/musella/public/higgs/forArabella/histograms_CMS-HGG.root", "read");
    inFile[1] = new TFile("/afs/cern.ch/user/m/musella/public/higgs/forArabella/histograms_CMS-HGG.root", "read");
    inFile[2] = new TFile("/afs/cern.ch/user/m/musella/public/higgs/forArabella/histograms_CMS-HGG.root", "read");
    inFile[3] = new TFile("/afs/cern.ch/user/m/musella/public/higgs/forArabella/histograms_CMS-HGG.root", "read");
    trees[0] = (TTree*)inFile[0]->Get("full_mva_trees/ggh_m125_8TeV");
    trees[1] = (TTree*)inFile[0]->Get("full_mva_trees/vbf_m125_8TeV");
    trees[2] = (TTree*)inFile[0]->Get("full_mva_trees/wzh_m125_8TeV");
    trees[3] = (TTree*)inFile[0]->Get("full_mva_trees/tth_m125_8TeV");
  }

  if(HggCatType == "CiC_7TeV"){
    inFile[0] = new TFile("../data/Hgg_signalTrees/Hgg_Dump_CiC_7TeV.root", "read");
//     inFile[1] = new TFile("../data/histograms_CMS-HGG.root", "read");
//     inFile[2] = new TFile("../data/histograms_CMS-HGG.root", "read");
//     inFile[3] = new TFile("../data/histograms_CMS-HGG.root", "read");
    trees[0] = (TTree*)inFile[0]->Get("opttree");
//     trees[1] = (TTree*)inFile[1]->Get("vbf_m125_8TeV"); 
//     trees[2] = (TTree*)inFile[2]->Get("wzh_m125_8TeV"); 
//     trees[3] = (TTree*)inFile[3]->Get("tth_m125_8TeV"); 
  }
  if(HggCatType == "MVA_7TeV"){
    inFile[0] = new TFile("../data/Hgg_signalTrees/Hgg_Dump_MVA_7TeV.root", "read");
//     inFile[1] = new TFile("../data/copiaHgg_MVA_125100.root", "read");
//     inFile[2] = new TFile("../data/copiaHgg_MVA_125400.root", "read");
//     inFile[3] = new TFile("../data/copiaHgg_MVA_125500.root", "read");
    trees[0] = (TTree*)inFile[0]->Get("opttree");
//     trees[1] = (TTree*)inFile[1]->Get("opttree");
//     trees[2] = (TTree*)inFile[2]->Get("opttree");
//     trees[3] = (TTree*)inFile[3]->Get("opttree");
  }

  std::cout << " Input presi " << std::endl;

  float full_weight;
  float full_cat;
  float dipho_mva;
  float dipho_pt;
  float mass;
  float et1, et2, eta1, eta2, r91, r92;
  int recoflag1, recoflag2;

  for(int iSample = 0; iSample < nSamples; ++iSample)
  {
    trees[iSample]->SetBranchStatus("*",0);
    trees[iSample]->SetBranchStatus("full_weight", 1);    trees[iSample]->SetBranchAddress("full_weight", &full_weight);
    trees[iSample]->SetBranchStatus("full_cat", 1);       trees[iSample]->SetBranchAddress("full_cat", &full_cat);
    trees[iSample]->SetBranchStatus("dipho_mva", 1);       trees[iSample]->SetBranchAddress("dipho_mva", &dipho_mva);
    trees[iSample]->SetBranchStatus("dipho_pt", 1);       trees[iSample]->SetBranchAddress("dipho_pt", &dipho_pt);
    trees[iSample]->SetBranchStatus("mass", 1);           trees[iSample]->SetBranchAddress("mass", &mass);
    trees[iSample]->SetBranchStatus("et1", 1);             trees[iSample]->SetBranchAddress("et1", &et1);
    trees[iSample]->SetBranchStatus("et2", 1);             trees[iSample]->SetBranchAddress("et2", &et2);
    trees[iSample]->SetBranchStatus("eta1", 1);           trees[iSample]->SetBranchAddress("eta1", &eta1);
    trees[iSample]->SetBranchStatus("eta2", 1);           trees[iSample]->SetBranchAddress("eta2", &eta2);
    trees[iSample]->SetBranchStatus("r91", 1);           trees[iSample]->SetBranchAddress("r91", &r91);
    trees[iSample]->SetBranchStatus("r92", 1);           trees[iSample]->SetBranchAddress("r92", &r92);
    trees[iSample]->SetBranchStatus("recoflag1", 1);           trees[iSample]->SetBranchAddress("recoflag1", &recoflag1);
    trees[iSample]->SetBranchStatus("recoflag2", 1);           trees[iSample]->SetBranchAddress("recoflag2", &recoflag2);
  }
  
  if(HggCatType == "MVA_8TeV")   nCats = 5;

  nCats = 8;
  
  std::cout << " Branches letti now OUTPUT HISTOS" << std::endl;

  //-----------------------------
  //OutputHistogram - all samples
  TH1F** Hgg_original = new TH1F*[nCats];
  for(int catIt = 0; catIt < nCats; ++catIt)
  {
    Hgg_original[catIt] = new TH1F(Form("Hgg_original_cat%d",catIt), "", 2400, 0., 300.);
    Hgg_original[catIt]->Sumw2();
  }
  TH1F** Hgg_original_pTLow = new TH1F*[nCats];
  for(int catIt = 0; catIt < nCats; ++catIt)
  {
    Hgg_original_pTLow[catIt] = new TH1F(Form("Hgg_original_pTLow_cat%d",catIt), "", 2400, 0., 300.);
    Hgg_original_pTLow[catIt]->Sumw2();
  }
  TH1F** Hgg_original_pTHigh = new TH1F*[nCats];
  for(int catIt = 0; catIt < nCats; ++catIt)
  {
    Hgg_original_pTHigh[catIt] = new TH1F(Form("Hgg_original_pTHigh_cat%d",catIt), "", 2400, 0., 300.);
    Hgg_original_pTHigh[catIt]->Sumw2();
  }
  
  TH1F** Hgg_measuredDiff = new TH1F*[nCats];
  for(int catIt = 0; catIt < nCats; ++catIt)
  {
    Hgg_measuredDiff[catIt] = new TH1F(Form("Hgg_measuredDiff_cat%d",catIt), "", 2400, 0., 300.);
    Hgg_measuredDiff[catIt]->Sumw2();
  }
  TH1F** Hgg_measuredDiff_pTLow = new TH1F*[nCats];
  for(int catIt = 0; catIt < nCats; ++catIt)
  {
    Hgg_measuredDiff_pTLow[catIt] = new TH1F(Form("Hgg_measuredDiff_pTLow_cat%d",catIt), "", 2400, 0., 300.);
    Hgg_measuredDiff_pTLow[catIt]->Sumw2();
  }
  TH1F** Hgg_measuredDiff_pTHigh = new TH1F*[nCats];
  for(int catIt = 0; catIt < nCats; ++catIt)
  {
    Hgg_measuredDiff_pTHigh[catIt] = new TH1F(Form("Hgg_measuredDiff_pTHigh_cat%d",catIt), "", 2400, 0., 300.);
    Hgg_measuredDiff_pTHigh[catIt]->Sumw2();
  }

 std::cout << " Branches letti now OUTPUT HISTOS  stat" << std::endl;
  
  TH1F** Hgg_measuredDiff_stat = new TH1F*[nCats*nTrials];
  for(int iTrial = 0; iTrial < nTrials; ++iTrial)
  {
    for(int catIt = 0; catIt < nCats; ++catIt)
    {
      Hgg_measuredDiff_stat[iTrial*nCats+catIt] = new TH1F(Form("Hgg_measuredDiff_stat_cat%d_trial%d",catIt,iTrial), "", 2400, 0., 300.);
      Hgg_measuredDiff_stat[iTrial*nCats+catIt]->Sumw2();
    }
  }
  TH1F** Hgg_measuredDiff_stat_pTLow = new TH1F*[nCats*nTrials];
  for(int iTrial = 0; iTrial < nTrials; ++iTrial)
  {
    for(int catIt = 0; catIt < nCats; ++catIt)
    {
      Hgg_measuredDiff_stat_pTLow[iTrial*nCats+catIt] = new TH1F(Form("Hgg_measuredDiff_stat_pTLow_cat%d_trial%d",catIt,iTrial), "", 2400, 0., 300.);
      Hgg_measuredDiff_stat_pTLow[iTrial*nCats+catIt]->Sumw2();
    }
  }
  TH1F** Hgg_measuredDiff_stat_pTHigh = new TH1F*[nCats*nTrials];
  for(int iTrial = 0; iTrial < nTrials; ++iTrial)
  {
    for(int catIt = 0; catIt < nCats; ++catIt)
    {
      Hgg_measuredDiff_stat_pTHigh[iTrial*nCats+catIt] = new TH1F(Form("Hgg_measuredDiff_stat_pTHigh_cat%d_trial%d",catIt,iTrial), "", 2400, 0., 300.);
      Hgg_measuredDiff_stat_pTHigh[iTrial*nCats+catIt]->Sumw2();
    }
  }


  std::cout << " Branches letti altri HISTOS " << std::endl;
  TH1F** SingleCat = new TH1F*[nCats];
  for(int catIt = 0; catIt < nCats; ++catIt)
    {
      SingleCat[catIt] = new TH1F(Form("SingleCat%d",catIt), "", 5, -1., 4.);
      SingleCat[catIt]->Sumw2();
    }

  TH1F** MassCorrection = new TH1F*[nCats];
  for(int catIt = 0; catIt < nCats; ++catIt)
    {
      MassCorrection[catIt] = new TH1F(Form("MassCorrection%d",catIt), "", 10000, 0.99, 1.01);
      MassCorrection[catIt]->Sumw2();
    }

  TH1F** ETSinglePho = new TH1F*[nCats];
  for(int catIt = 0; catIt < nCats; ++catIt)
    {
      ETSinglePho[catIt] = new TH1F(Form("ETSinglePho%d",catIt), "", 1000, 0., 1000.);
      ETSinglePho[catIt]->Sumw2();
    }

  TH1F** HTDoublePho = new TH1F*[nCats];
  for(int catIt = 0; catIt < nCats; ++catIt)
    {
      HTDoublePho[catIt] = new TH1F(Form("HTDoublePho%d",catIt), "", 1000, 0., 1000.);
      HTDoublePho[catIt]->Sumw2();
    }
  int totPho_cat0;
  int totGS1_cat0 = 0;
  int totGS6_cat0 = 0;
  int totGS1_b150_cat0 = 0;
  int totGS6_b150_cat0 = 0;

  int totPho_cat1;
  int totGS1_cat1 = 0;
  int totGS6_cat1 = 0;
  int totGS1_b150_cat1 = 0;
  int totGS6_b150_cat1 = 0;


  std::cout << " OK done 2 " << std::endl;

  //-------------------
  // loop - all samples
  for(int iSample = 0; iSample < nSamples; ++iSample)
  {
    std::cout << " iSample =  " << iSample << std::endl;
    inFile[iSample]->cd();
    for(int entry = 0; entry < trees[iSample]->GetEntries(); ++entry)
    {
      if( entry%10000 == 0 ) std::cout << " >>> reading entry " << entry << " / " << trees[iSample]->GetEntries() << "\r" << std::endl;
      trees[iSample]->GetEntry(entry);
      
      float theta1 = 2*atan(exp(-eta1));
      float theta2 = 2*atan(exp(-eta2));
      float Rt1 = sin(theta1);
      float Rt2 = sin(theta2);
      
      int eventMVAcat = -1;
      if( HggCatType == "MVA_7TeV") eventMVAcat = GetMVAcat_7TeV(dipho_mva);
      if( HggCatType == "MVA_8TeV") eventMVAcat = GetMVAcat_8TeV(dipho_mva);
      if( HggCatType == "CiC_7TeV" || HggCatType == "CiC_8TeV") eventMVAcat = full_cat;
      
      for(int iCat = 0; iCat < nCats; ++iCat)
      {
        if(eventMVAcat == iCat)
        {
	  int pho1Cat = iCat;
	  if(analysis == "stdCat") pho1Cat = GetSingleCategory(eta1,r91);
	  int pho2Cat = iCat;
	  if(analysis == "stdCat") pho2Cat = GetSingleCategory(eta2,r92);

// 	  std::cout << " pho1Cat = " << pho1Cat << std::endl;
// 	  std::cout << " pho2Cat = " << pho2Cat << std::endl;
	  if( (fabs(eta1) >= 1.4442 && fabs(eta1) <= 1.5600) || (fabs(eta2) >= 1.4442 && fabs(eta2) <= 1.5600) ) continue;
	  if(pho1Cat == -1 || pho2Cat == -1) continue;
	  SingleCat[iCat]->Fill(pho1Cat);
	  SingleCat[iCat]->Fill(pho2Cat);


	  ///////////////
	  if(iCat == 0){
	    totPho_cat0 += 2;
	    if(recoflag1 == 17)	    totGS1_cat0 += 1;
	    if(recoflag1 == 16)	    totGS6_cat0 += 1;
	    if(recoflag1 == 17 && dipho_pt < 150.)     totGS1_b150_cat0 += 1;
	    if(recoflag1 == 16 && dipho_pt < 150.)     totGS6_b150_cat0 += 1;
	    if(recoflag2 == 17)	    totGS1_cat0 += 1;
	    if(recoflag2 == 16)	    totGS6_cat0 += 1;
	    if(recoflag2 == 17 && dipho_pt < 150.)     totGS1_b150_cat0 += 1;
	    if(recoflag2 == 16 && dipho_pt < 150.)     totGS6_b150_cat0 += 1;
	  }
	  if(iCat == 1){
	    totPho_cat1 += 2;
	    if(recoflag1 == 17)	    totGS1_cat1 += 1;
	    if(recoflag1 == 16)	    totGS6_cat1 += 1;
	    if(recoflag1 == 17 && dipho_pt < 150.)     totGS1_b150_cat1 += 1;
	    if(recoflag1 == 16 && dipho_pt < 150.)     totGS6_b150_cat1 += 1;
	    if(recoflag2 == 17)	    totGS1_cat1 += 1;
	    if(recoflag2 == 16)	    totGS6_cat1 += 1;
	    if(recoflag2 == 17 && dipho_pt < 150.)     totGS1_b150_cat1 += 1;
	    if(recoflag2 == 16 && dipho_pt < 150.)     totGS6_b150_cat1 += 1;
	  }
	  ///////////////

	  //	  if(iCat == 0 )std::cout << " dipho_pt = " << dipho_pt << std::endl;

	  //	  if( (dipho_pt/mass) < (40./125.) && iCat == 0) std::cout << " evento low pT cat 0 " << std::endl;

	  Hgg_original[iCat]->Fill(mass,full_weight);
	  if( (dipho_pt/mass) < (40./125.) ) Hgg_original_pTLow[iCat]->Fill(mass,full_weight);
	  else Hgg_original_pTHigh[iCat]->Fill(mass,full_weight);
  
//	    if(et1 > 150. || et2 > 150.) continue;

	  MassCorrection[iCat]->Fill(sqrt(f_std[pho1Cat]->Eval(et1) * f_std[pho2Cat]->Eval(et2)), full_weight);
	  ETSinglePho[iCat]->Fill(et1);
	  ETSinglePho[iCat]->Fill(et2);
	  HTDoublePho[iCat]->Fill(et1+et2);


	  Hgg_measuredDiff[iCat]->Fill(mass*sqrt(f_std[pho1Cat]->Eval(et1) * f_std[pho2Cat]->Eval(et2)), full_weight);
	  if( (dipho_pt/mass) < (40./125.) ) 
	    Hgg_measuredDiff_pTLow[iCat]->Fill(mass*sqrt(f_std[pho1Cat]->Eval(et1) * f_std[pho2Cat]->Eval(et2)), full_weight);
	  else Hgg_measuredDiff_pTHigh[iCat]->Fill(mass*sqrt(f_std[pho1Cat]->Eval(et1) * f_std[pho2Cat]->Eval(et2)), full_weight);

// 	  std::cout << " f_std[pho1Cat]->Eval(et1) = " << f_std[pho1Cat]->Eval(et1) << std::endl;
// 	  std::cout << " f_std[pho2Cat]->Eval(et2) = " << f_std[pho2Cat]->Eval(et2) << std::endl;
// 	  std::cout << " corr = " << sqrt(f_std[pho1Cat]->Eval(et1) * f_std[pho2Cat]->Eval(et2)) << std::endl;


	  //	  std::cout << " fino a qui ok iCat = " << iCat << std::endl;
        }
      }
    }
  
    std::cout << " ALL integral sample " << iSample << std::endl;
    for(int iCat = 0; iCat < nCats; ++iCat)
      std::cout << " cat" << iCat << ": " << Hgg_original[iCat]->Integral() << std::endl;

    std::cout << " ok passa a sample dopo " << std::endl;
  }

  /*
  std::cout << "STATISTICA GAIN SWITCH" << std::endl;
  std::cout << " >>>>> totPho_cat0 = " << totPho_cat0 << std::endl;
  std::cout << " totGS1_cat0 = " << totGS1_cat0 << std::endl;
  std::cout << " totGS6_cat0 = " << totGS6_cat0 << std::endl;
  std::cout << " totGS1_b150_cat0 = " << totGS1_b150_cat0 << std::endl;
  std::cout << " totGS6_b150_cat0 = " << totGS6_b150_cat0 << std::endl;
  std::cout << " >>>>> totPho_cat1 = " << totPho_cat1 << std::endl;
  std::cout << " totGS1_cat1 = " << totGS1_cat1 << std::endl;
  std::cout << " totGS6_cat1 = " << totGS6_cat1 << std::endl;
  std::cout << " totGS1_b150_cat1 = " << totGS1_b150_cat1 << std::endl;
  std::cout << " totGS6_b150_cat1 = " << totGS6_b150_cat1 << std::endl;
  */
  std::cout << " OK done 3 " << std::endl;

  std::cout << " ALL integral" << std::endl;
  for(int iCat = 0; iCat < nCats; ++iCat)
    std::cout << " cat" << iCat << ": " << Hgg_original[iCat]->Integral() << std::endl;
  
  TProfile* envelop_cat0 = new TProfile("envelop_cat0", "", 1000, 0., 500., "s");
  TProfile* envelop_cat1 = new TProfile("envelop_cat1", "", 1000, 0., 500., "s");
  TProfile* envelop_cat2 = new TProfile("envelop_cat2", "", 1000, 0., 500., "s");
  TProfile* envelop_cat3 = new TProfile("envelop_cat3", "", 1000, 0., 500., "s");
  
  nCats = 4;

  TF1** f_stat_Trial = new TF1*[nCats];
  for(int iCat = 0; iCat < nCats; ++iCat)
    {
      if(iCat == 0)      f_stat_Trial[iCat] = new TF1(Form("f_stat_Trial_cat%d_std",iCat), BandaCat0,20., 1000., 1);
      if(iCat == 1)      f_stat_Trial[iCat] = new TF1(Form("f_stat_Trial_cat%d_std",iCat), BandaCat1,20., 1000., 1);
      if(iCat == 2)      f_stat_Trial[iCat] = new TF1(Form("f_stat_Trial_cat%d_std",iCat), BandaCat2,20., 1000., 1);
      if(iCat == 3)      f_stat_Trial[iCat] = new TF1(Form("f_stat_Trial_cat%d_std",iCat), BandaCat3,20., 1000., 1);
    }
  TF1** f_stat_Trial_ET = new TF1*[nCats];
  for(int iCat = 0; iCat < nCats; ++iCat)
    {
      if(iCat == 0)      f_stat_Trial_ET[iCat] = new TF1(Form("f_stat_Trial_ET_cat%d_std",iCat), BandaCat0 ,20., 1000., 1);
      if(iCat == 1)      f_stat_Trial_ET[iCat] = new TF1(Form("f_stat_Trial_ET_cat%d_std",iCat), BandaCat1 ,20., 1000., 1);
      if(iCat == 2)      f_stat_Trial_ET[iCat] = new TF1(Form("f_stat_Trial_ET_cat%d_std",iCat), BandaCat2 ,20., 1000., 1);
      if(iCat == 3)      f_stat_Trial_ET[iCat] = new TF1(Form("f_stat_Trial_ET_cat%d_std",iCat), BandaCat3 ,20., 1000., 1);
    }


  //-----------------------------------
  // loop - propagate statistical error
  for(int iTrial = 0; iTrial < nTrials; ++iTrial)
  {
    std::cout << ">>> iTrial: " << iTrial << " / " << nTrials << "\r" << std::flush;
    
    for(int iCat = 0; iCat < nCats; ++iCat)
    {
      double* fiParRndm = new double[1];
      fiParRndm[0] = gRandom->Gaus(0.,1.);
      
      f_stat_Trial[iCat] -> SetParameters( fiParRndm);
      f_stat_Trial_ET[iCat] -> SetParameters( fiParRndm);

     
      delete fiParRndm;

      if(iCat == 0)
	for(int nBin = 1; nBin<envelop_cat0->GetNbinsX(); ++nBin) {
	    envelop_cat0->Fill(envelop_cat0->GetBinCenter(nBin), f_stat_Trial_ET[iCat]->Eval(envelop_cat0->GetBinCenter(nBin)));
	}
      if(iCat == 1)
	for(int nBin = 1; nBin<envelop_cat1->GetNbinsX(); ++nBin) {
	  envelop_cat1->Fill(envelop_cat1->GetBinCenter(nBin), f_stat_Trial_ET[iCat]->Eval(envelop_cat1->GetBinCenter(nBin)));
	}
      if(iCat == 2)
	for(int nBin = 1; nBin<envelop_cat2->GetNbinsX(); ++nBin) {
	    envelop_cat2->Fill(envelop_cat2->GetBinCenter(nBin), f_stat_Trial_ET[iCat]->Eval(envelop_cat2->GetBinCenter(nBin)));
	}
      if(iCat == 3)
	for(int nBin = 1; nBin<envelop_cat3->GetNbinsX(); ++nBin) {
	  envelop_cat3->Fill(envelop_cat3->GetBinCenter(nBin), f_stat_Trial_ET[iCat]->Eval(envelop_cat3->GetBinCenter(nBin)));
	}

    } 
    //       continue;

    for(int iSample = 0; iSample < nSamples; ++iSample)
    {
      for(int entry = 0; entry < trees[iSample]->GetEntries(); ++entry)
      {
        //if( entry%10000 == 0 ) std::cout << " >>> reading entry " << entry << " / " << trees[iSample]->GetEntries() << "\r" << std::endl;
        trees[iSample]->GetEntry(entry);
        
        float theta1 = 2*atan(exp(-eta1));
        float theta2 = 2*atan(exp(-eta2));
        float Rt1 = sin(theta1);
        float Rt2 = sin(theta2);
        
	int eventMVAcat = -1;
	if( HggCatType == "MVA_7TeV" ) eventMVAcat = GetMVAcat_7TeV(dipho_mva);
	if( HggCatType == "MVA_8TeV" ) eventMVAcat = GetMVAcat_8TeV(dipho_mva);
	if( HggCatType == "CiC_7TeV" || HggCatType == "CiC_8TeV") eventMVAcat = full_cat;

	if(HggCatType == "MVA_8TeV") nCats = 5;

	nCats = 8;
        for(int iCat = 0; iCat < nCats; ++iCat)
        {
          if(eventMVAcat == iCat)
          {
	    int pho1Cat = iCat;
	    if(analysis == "stdCat") pho1Cat = GetSingleCategory(eta1,r91);
	    int pho2Cat = iCat;
	    if(analysis == "stdCat") pho2Cat = GetSingleCategory(eta2,r92);
	    if( (fabs(eta1) >= 1.4442 && fabs(eta1) <= 1.5600) || (fabs(eta2) >= 1.4442 && fabs(eta2) <= 1.5600) ) continue;
	    if(pho1Cat == -1 || pho2Cat == -1) continue;


	      Hgg_measuredDiff_stat[iTrial*nCats+iCat]->Fill(mass*sqrt(f_stat_Trial_ET[pho1Cat]->Eval(et1) * 
								       f_stat_Trial_ET[pho2Cat]->Eval(et2)), full_weight);
	      if( (dipho_pt/mass) < (40./125.) )
		Hgg_measuredDiff_stat_pTLow[iTrial*nCats+iCat]->Fill(mass*sqrt(f_stat_Trial_ET[pho1Cat]->Eval(et1) * 
									       f_stat_Trial_ET[pho2Cat]->Eval(et2)), full_weight);
	      else 
		Hgg_measuredDiff_stat_pTHigh[iTrial*nCats+iCat]->Fill(mass*sqrt(f_stat_Trial_ET[pho1Cat]->Eval(et1) * 
										f_stat_Trial_ET[pho2Cat]->Eval(et2)), full_weight);

          }
        }
	nCats = 4;
      }
    }
  }
  std::cout << std::endl;

  std::cout << " Hgg_measuredDiff[0] = " << Hgg_measuredDiff[0]->GetMean() << std::endl;
  std::cout << " Hgg_measuredDiff_stat[0] = " << Hgg_measuredDiff_stat[0]->GetMean() << std::endl;
  
  TFile outputFile("rescaledHgg_histos.root","RECREATE");
  outputFile.cd();
  if(HggCatType == "MVA_8TeV") nCats = 5;

  nCats = 8;
  for(int iCat = 0; iCat < nCats; ++iCat)
  {
    SingleCat[iCat]->Write();
    MassCorrection[iCat]->Write();
    ETSinglePho[iCat]->Write();
    HTDoublePho[iCat]->Write();

    Hgg_original[iCat]->Write();
    Hgg_measuredDiff[iCat]->Write();

    Hgg_original_pTLow[iCat]->Write();
    Hgg_original_pTHigh[iCat]->Write();
    Hgg_measuredDiff_pTLow[iCat]->Write();
    Hgg_measuredDiff_pTHigh[iCat]->Write();

  }

  std::cout << " OK done 4 " << std::endl;

  envelop_cat0->Write();  
  envelop_cat1->Write();
  envelop_cat2->Write();
  envelop_cat3->Write();

  outputFile.mkdir("stat");
  outputFile.cd("stat");
  //  nCats = 4;
  for(int iTrial = 0; iTrial < nTrials; ++iTrial)
  {
    for(int iCat = 0; iCat < nCats; ++iCat)
    {
      Hgg_measuredDiff_stat[iTrial*nCats+iCat]->Write();  
      Hgg_measuredDiff_stat_pTLow[iTrial*nCats+iCat]->Write();  
      Hgg_measuredDiff_stat_pTHigh[iTrial*nCats+iCat]->Write();  
    }
  }
  outputFile.Close();

  return 0;
}



double* cholesky(double* A, int n)
{
  double *L = (double*)calloc(n * n, sizeof(double));
  if (L == NULL)
    exit(EXIT_FAILURE);
 
  for (int i = 0; i < n; i++)
    for (int j = 0; j < (i+1); j++) {
      double s = 0;
      for (int k = 0; k < j; k++)
	s += L[i * n + k] * L[j * n + k];
      L[i * n + j] = (i == j) ?
	sqrt(A[i * n + i] - s) :
	(1.0 / L[j * n + j] * (A[i * n + j] - s));
    }
  
  return L;
}




int GetSingleCategory(const float& scEta1, const float& R91)
{
//   std::cout << " in GetSingleCategory scEta1 = " << scEta1 << std::endl;
//   std::cout << " in GetSingleCategory R91 = " << R91 << std::endl;

  if( fabs(scEta1) < 1.4442 )
    {
      if( R91 > 0.94 ) return 0;
      if( R91 < 0.94 ) return 1;
    }

  if( fabs(scEta1) > 1.5600 )
    {
      if( R91 > 0.94 ) return 2;
      if( R91 < 0.94 ) return 3;
    }

  return -1;
}



int GetMVAcat_7TeV(const float& dipho_mva){
  // 0.93, 0.85, 0.70, 0.19
  if( dipho_mva > 0.93) return 0;

  if( dipho_mva > 0.85 && dipho_mva < 0.93 ) return 1;

  if( dipho_mva > 0.70 && dipho_mva < 0.85 ) return 2;

  if( dipho_mva > 0.19 && dipho_mva < 0.70 ) return 3;

  return -1;
}


int GetMVAcat_8TeV(const float& dipho_mva){
  // 0.76,0.36,0.00,-0.42,-0.78
  if( dipho_mva > 0.76) return 0;

  if( dipho_mva > 0.36 && dipho_mva < 0.76 ) return 1;

  if( dipho_mva > 0.00 && dipho_mva < 0.36 ) return 2;

  if( dipho_mva > -0.42 && dipho_mva < 0.00 ) return 3;

  if( dipho_mva > -0.78 && dipho_mva < -0.42 ) return 4;

  return -1;
}
