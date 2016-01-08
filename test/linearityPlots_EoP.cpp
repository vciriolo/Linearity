#include "setTDRStyle.h"
#include "ntpleUtils.h"

#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>

#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TVirtualFitter.h"
#include "TMatrixDSym.h"
#include "TFitResult.h"


bool fitTemplate = true;
bool MCClosure   = false;
bool drawFitFunc = true;
bool rescaleErrors = true;
bool addSyst = true;
bool drawScaleSys = false;
bool drawSmearSys = false;
bool writeFunctions = false;
//std::string analysis  = "CiC";
std::string analysis = "stdCat";

//std::string fitMethod = "exp";
//std::string fitMethod = "exp3par";
std::string fitMethod = "pol1";
//std::string fitMethod = "pol0";




int main(int argc, char** argv)
{
  //Check if all nedeed arguments to parse are there
  if( argc < 2 )
  {
    std::cerr << " >>> linearityPlots_Eop::usage: " << argv[0] << " evtsPerPoint_cat0 (evtsPerPoint_cat1 evtsPerPoint_cat2 ...)" << std::endl;
    return -1;
  }
  
  int nCat = argc - 1;
  int* evtsPerPoint = new int[nCat];
  for(int iCat = 0; iCat < nCat; ++iCat)
    evtsPerPoint[iCat] = atoi(argv[iCat+1]);
  
  
  std::vector<std::string> directories;
  //  directories.push_back( ("data/2012/Winter2013/EoP_"+analysis+"_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15_10m4_BinErrsOK").c_str() );
  //  directories.push_back( "data/2012/Winter2013/EoP_finalPlots" );
  directories.push_back( ("data/2012/Winter2013/EoP_"+analysis+"_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15").c_str() );
  if(drawSmearSys == true){
directories.push_back(("data/2012/Winter2013/EoP_"+analysis+"_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15_10m4_BinErrsOK_ExtraSmearingPlus1Err").c_str());
directories.push_back(("data/2012/Winter2013/EoP_"+analysis+"_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15_10m4_BinErrsOK_ExtraSmearingMinus1Err").c_str());
  }
  if(drawScaleSys == true){
directories.push_back(("data/2012/Winter2013/EoP_"+analysis+"_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15_10m4_BinErrsOK_ScalePlus1Err").c_str());
directories.push_back(("data/2012/Winter2013/EoP_"+analysis+"_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15_10m4_BinErrsOK_ScaleMinus1Err").c_str());
  }
  //directories.push_back( "data/2012/Winter2013/EoP_CiC_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15_10m4_BinErrsOK_P0" );
  //directories.push_back( "data/2012/Winter2013/EoP_CiC_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15_10m4_BinErrsOK_EXP" );
  unsigned int nDir = directories.size();
  
  float SysError[nCat][6];
  SysError[0][0] = -0.0101562;
  SysError[0][1] = -0.0193683;
  SysError[0][2] = -0.0229157;
  SysError[0][3] = 9.82711e-05;
  SysError[0][4] = -0.0209955;
  SysError[0][5] = 0.00840904;

  SysError[1][0] = -0.0309876;
  SysError[1][1] = -0.0224256;
  SysError[1][2] = -0.0424731;
  SysError[1][3] = -0.0186554;
  SysError[1][4] = -0.0585795;
  SysError[1][5] = 0.;

  SysError[2][0] = -0.0335435;
  SysError[2][1] = -0.0138287;
  SysError[2][2] = -0.0234936;
  SysError[2][3] = 0.027082;
  SysError[2][4] = -0.0402911;
  SysError[2][5] = 0.;

  SysError[3][0] = -0.0377126;
  SysError[3][1] = -0.070253;
  SysError[3][2] = -0.0345543;
  SysError[3][3] = -0.0577678;
  SysError[3][4] = 0.130369;
  SysError[3][5] = 0.;
  
  //-----------------------------
  // Decide which methods to draw
  
  std::vector<std::string> methods;
      methods.push_back( "fit" );
  //methods.push_back( "gausFit" );
  //  methods.push_back( "mean" );
  //  methods.push_back( "recursiveMean" );
  //  methods.push_back( "smallestInterval" );
  int nMethods = methods.size();
  
  std::vector<std::string> nameTrial;
  nameTrial.push_back("std");
  if(drawSmearSys == true){
    nameTrial.push_back("smearPlus1");
    nameTrial.push_back("smearMinus1");
  }
  if(drawScaleSys == true){
    nameTrial.push_back("scalePlus1");
    nameTrial.push_back("scaleMinus1");
  }

  //--------------------
  // Set fitting options
  
  TVirtualFitter::SetDefaultFitter("Minuit2");
  
  setTDRStyle();
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.17);
  gStyle->SetLabelSize(0.04,"XYZ");
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  std::vector<int> colors;
  std::vector<int> linestyles;
  std::vector<int> markerstyles;
  colors.push_back(kRed+1);
  colors.push_back(kOrange+1);
  colors.push_back(kGreen+1);
  colors.push_back(kBlue+1);
  linestyles.push_back(1);
  linestyles.push_back(2);
  linestyles.push_back(2);
  markerstyles.push_back(20);
  markerstyles.push_back(26);
  markerstyles.push_back(32);
  
  
  
  //----------------
  // Define canvases
  
  TCanvas* c_all = new TCanvas("c_all","c_all");
  c_all -> Divide(2,2);
  
  TCanvas** c = new TCanvas*[nCat];
  
  
  
  //----------------------------
  // Define infiles and canvases
  for(unsigned int iDir = 0; iDir < directories.size(); ++iDir)
  {
    std::string directory = directories.at(iDir);
    std::string* inFileNames = new std::string[nCat];
    
    std::string baseFileName = "studyLinearity_EoP_";
    if( MCClosure == true) directory += "_MCClosure";
    
    for(int iMeth = 0; iMeth < nMethods; ++iMeth)
    {
      std::cout << ">>> method: " << methods.at(iMeth) << std::endl;
      
      TLegend** legend = new TLegend*[nCat];
      for(int iCat = 0; iCat < nCat; ++iCat)
      {
        legend[iCat] = new TLegend(0.16, 0.77, 0.30, 0.92);
        legend[iCat] -> SetFillColor(kWhite);
        legend[iCat] -> SetFillStyle(1001);  
        legend[iCat] -> SetTextFont(42);  
        legend[iCat] -> SetTextSize(0.05);
      }
      
      TLatex* latex = new TLatex(0.14,0.96,Form("analysis: %s   scale estimator: %s",analysis.c_str(),(methods.at(iMeth)).c_str()));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      
      double* scale_EoP = new double[nCat];
      
      TGraphAsymmErrors** g = new TGraphAsymmErrors*[nCat];
      
      TF1** f_fit = new TF1*[nCat];
      TH1F** hint = new TH1F*[nCat];
      
      
      for(int iCat = 0; iCat < nCat; ++iCat)
      {
        std::cout << "\n***************** cat: " << iCat << " *****************" << std::endl;
        
        char EvtString[50];
        sprintf(EvtString,"cat%d_%devtsPerPoint",iCat,evtsPerPoint[iCat]);
        
        inFileNames[iCat] = directory + "/" + baseFileName + std::string(EvtString) + ".root";
        TFile* f = TFile::Open((inFileNames[iCat]).c_str(),"READ");
        std::cout << ">>> inFileName: " << inFileNames[iCat] << std::endl;
        
        char graphName[50];
        sprintf(graphName,"step1/scale_%s_DAOverMC",methods.at(iMeth).c_str());
        g[iCat] = (TGraphAsymmErrors*)( f->Get(graphName) );
        
	if(addSyst == true && fitTemplate == false){
	  for(int point = 0; point < g[iCat]->GetN(); ++point){
	    double ey = g[iCat] -> GetErrorY(point);
	    double x,y;
	    g[iCat] -> GetPoint(point, x, y);
// 	    g[iCat] -> SetPointEYhigh(point,ey/100.);
// 	    g[iCat] -> SetPointEYlow(point,ey/100.);
	    double statE2 = pow(ey,2);
//  	    double sysE2 = pow(SysCat[iCat], 2);                     
//  	    g[iCat] -> SetPointEYhigh(point, sqrt(statE2 + sysE2) ); 
//  	    g[iCat] -> SetPointEYlow (point, sqrt(statE2 + sysE2) ); 
	    double sysE2 = pow(0.01*SysError[iCat][point], 2);
	    g[iCat] -> SetPointEYhigh(point, sqrt(statE2 + sysE2) );
	    g[iCat] -> SetPointEYlow (point, sqrt(statE2 + sysE2) );
	  }
	}
	if(fitTemplate){
	  for(int point = 0; point < g[iCat]->GetN(); ++point){
	    double ey = g[iCat] -> GetErrorY(point);
 	    g[iCat] -> SetPointEYhigh(point, ey/100.);
 	    g[iCat] -> SetPointEYlow (point, ey/100.);
	  }
	}


        if( iDir == 0 ) c[iCat] = new TCanvas();
        c[iCat] -> cd();
        c[iCat] -> SetGridx();
        c[iCat] -> SetGridy();
        //c[iCat] -> SetLogx();
        
        c_all -> cd(iCat+1);
        gPad -> SetGridx();
        gPad -> SetGridy();
        
        g[iCat] -> SetPoint(g[iCat]->GetN(),1500.,1.);
        g[iCat] -> GetXaxis() -> SetRangeUser(20,500.);        
        g[iCat] -> GetXaxis() -> SetMoreLogLabels();
        g[iCat] -> GetXaxis() -> SetTitle("E_{T} [GeV]");
        g[iCat] -> GetYaxis() -> SetTitle("#LTEoP#GT^{data} / #LTEoP#GT^{MC}");
        g[iCat] -> GetXaxis() -> SetTitleSize(0.05);
        g[iCat] -> GetYaxis() -> SetTitleSize(0.05);
        g[iCat] -> GetXaxis() -> SetTitleOffset(1.10);
        g[iCat] -> GetYaxis() -> SetTitleOffset(1.28);
        g[iCat] -> GetXaxis() -> SetLabelSize(0.04);
        g[iCat] -> GetYaxis() -> SetLabelSize(0.04);
        g[iCat] -> GetYaxis() -> SetNdivisions(405);
        g[iCat] -> SetMarkerSize(1.5);
        g[iCat] -> SetMarkerColor(colors.at(iCat));
        g[iCat] -> SetMarkerStyle(markerstyles.at(iDir));
        g[iCat] -> SetLineColor(kBlack);
        g[iCat] -> SetLineWidth(1);
        
	g[iCat] -> SetMinimum(0.99);
	g[iCat] -> SetMaximum(1.01);
	if(iCat == 2 || iCat == 3){
	  g[iCat] -> SetMinimum(0.95);
	  g[iCat] -> SetMaximum(1.05);
	}
	if(MCClosure == true){
	  g[iCat] -> SetMinimum(0.99);
	  g[iCat] -> SetMaximum(1.01);
        }
        if( drawFitFunc == false )
        {
          if( iDir == 0 )
          {
            c[iCat] -> cd();
            g[iCat] -> Draw("AP");
            c_all -> cd(iCat+1);
            g[iCat] -> Draw("AP"); 
          }
          else
          {
            c[iCat] -> cd();
            g[iCat] -> Draw("P,same");
            c_all -> cd(iCat+1);
            g[iCat] -> Draw("P,same");
          }
          
          if(MCClosure)
          {
            TF1* f_scaleVsEt = new TF1("f_scaleVsEt", "1. + [0] * (1 - exp(-[1] * (x-45.)) )",0., 1000.);
            f_scaleVsEt -> SetParameters(7.50e-03,2.00e-02);
	    //TF1* f_scaleVsEt = new TF1("f_scaleVsEt", "1.+0.000",0., 1000.);
            f_scaleVsEt->SetLineColor(kBlack);
            c[iCat] -> cd();
            f_scaleVsEt->Draw("same");
            c_all -> cd(iCat+1);
            f_scaleVsEt->Draw("same");
          }
        }
        else
        {
          if( iDir == 0 )
          {
            c[iCat] -> cd();
            g[iCat] -> Draw("AP");
            c_all -> cd(iCat+1);
            g[iCat] -> Draw("AP");
	  }
          else
          {
            c[iCat] -> cd();
            g[iCat] -> Draw("P,same");
            c_all -> cd(iCat+1);
            g[iCat] -> Draw("P,same");
	  }
        }
        
        
        char funcName[150];
        sprintf(funcName,"f_prefit_%s_%d",methods.at(iMeth).c_str(),iCat);
        TF1* f_prefit;
        if( fitMethod == "pol0")
        {
          f_prefit = new TF1(funcName,"[0]",20., 230.);
          f_prefit -> SetParameter(0,1.);         
        }
        if( fitMethod == "pol1")
        {
          f_prefit = new TF1(funcName,"[0]+[1]*(x-45.)",20., 5000.);
          f_prefit -> SetParameters(1.,0.00002);         
        }
        if( fitMethod == "exp3par")
        {
          f_prefit = new TF1(funcName,"1.+[0]*(1.-exp(-1.*[1]*(x-45.)))+[2]",20., 5000.);
          f_prefit -> SetParameters(0.005,0.02,0.);
        }
        if( fitMethod == "exp")
        {
	  f_prefit = new TF1(funcName,"1.+[0]*(1.-exp(-1.*[1]*(x-45.)) )",20., 5000.);
          f_prefit -> SetParameters(0.005,0.02);
          f_prefit -> SetParLimits(0,0.,1.);
          f_prefit -> SetParLimits(1,0.,0.05);
        }
        
        TFitResultPtr fitResult1;
        int fitStatus1 = -1;
        g[iCat] -> Fit(funcName,"QRHNS");
        fitStatus1 = fitResult1;
        
        //std::cout << " fitStatus1 = " << fitStatus1 << std::endl;
        //std::cout << " f_prefit->GetParameter(0) = " << f_prefit->GetParameter(0) << std::endl;
        //std::cout << " f_prefit->GetParameter(1) = " << f_prefit->GetParameter(1) << std::endl;
	std::cout << " f_prefit->GetChisquare()/f_prefit->GetNDF() = " << f_prefit->GetChisquare()/f_prefit->GetNDF() << std::endl;
        
        if( drawFitFunc == true && rescaleErrors == true )
        {
          for(int point = 0; point < g[iCat]->GetN(); ++point)
          {
 	    double ey = g[iCat] -> GetErrorY(point);
// 	    if(addSys == false){
// 	      g[iCat] -> SetPointEYhigh(point,ey*sqrt(f_prefit->GetChisquare()/f_prefit->GetNDF()));
// 	      g[iCat] -> SetPointEYlow (point,ey*sqrt(f_prefit->GetChisquare()/f_prefit->GetNDF()));
// 	    }
// 	    if(SysError[iCat][point] != 0. && addSys == true){
// 	      double statE2 = pow(ey,2);
// 	      double sysE2 = pow(SysError[iCat][point], 2);
// 	      double errorN2 =  (statE2) * (f_prefit->GetChisquare()/f_prefit->GetNDF()) / sysE2 - statE2/sysE2;
// 	      g[iCat] -> SetPointEYhigh(point, sqrt( statE2 + errorN2*sysE2) );
// 	      g[iCat] -> SetPointEYlow (point, sqrt( statE2 + errorN2*sysE2) );
// 	    }
	  }
        }
        if( drawFitFunc == false && rescaleErrors == true )      
        {
          if( MCClosure == true )
          {
 	    for(int point = 0; point < g[iCat]->GetN(); ++point)
            {
              double ey = g[iCat] -> GetErrorY(point);
	      if(addSyst == false){
		g[iCat] -> SetPointEYhigh(point,ey*sqrt(f_prefit->GetChisquare()/f_prefit->GetNDF()));
		g[iCat] -> SetPointEYlow (point,ey*sqrt(f_prefit->GetChisquare()/f_prefit->GetNDF()));
	      }
	      if(SysError[iCat][point] != 0. && addSyst == true){
		double statE2 = pow(ey,2);
		double sysE2 = pow(SysError[iCat][point], 2);
		double errorN2 =  (statE2) * (f_prefit->GetChisquare()/f_prefit->GetNDF()) / sysE2 - statE2/sysE2;
		g[iCat] -> SetPointEYhigh(point, sqrt( statE2 + errorN2*sysE2) );
		g[iCat] -> SetPointEYlow (point, sqrt( statE2 + errorN2*sysE2) );
	      }
            }
          }
        }

	////////////////                                                                                                                         
	int fitStatus1b = -1;
	g[iCat] -> Fit(funcName,"QRHNS");
	fitStatus1b = fitResult1;
	std::cout << " f_prefit->GetChisquare()/f_prefit->GetNDF() = " << f_prefit->GetChisquare()/f_prefit->GetNDF() << std::endl;
	if( rescaleErrors == true )
	  {
	    for(int point = 0; point < g[iCat]->GetN(); ++point)
	      {
		double ey = g[iCat] -> GetErrorY(point);
		if(addSyst == false || fitTemplate == true){
		  g[iCat] -> SetPointEYhigh(point,ey*sqrt(f_prefit->GetChisquare()/f_prefit->GetNDF()));
		  g[iCat] -> SetPointEYlow (point,ey*sqrt(f_prefit->GetChisquare()/f_prefit->GetNDF()));
		}
	      }
	  }
	////////////                                                                                                                             

        
        sprintf(funcName,"f_fit_%s_%d",methods.at(iMeth).c_str(),iCat);
        if( fitMethod == "pol0")
        {
          f_fit[iCat] = new TF1(funcName,"[0]",20., 5000.);
          f_fit[iCat] -> SetParameter(0,1.);         
        }
        if( fitMethod == "pol1" )
        {
	  f_fit[iCat] = new TF1(funcName,"[0]+[1]*(x-45.)",20., 5000.);
          f_fit[iCat] -> SetParameters(1.,0.00001);
        }
        if( fitMethod == "exp3par")
        {
          f_fit[iCat] = new TF1(funcName,"1.+[0]*(1.-exp(-1.*[1]*(x-45.)))+[2]",20., 5000.);
          f_fit[iCat] -> SetParameters(0.005,0.02,0.);
        }
        if( fitMethod == "exp")
        {
          f_fit[iCat] = new TF1(funcName,"1.+[0]*(1.-exp(-1.*[1]*(x-45.)))",20., 5000.);
          f_fit[iCat] -> SetParameters(0.005,0.02);
        }
        f_fit[iCat] -> SetLineColor(kBlue+2);
        f_fit[iCat] -> SetLineWidth(3);
        f_fit[iCat] -> SetLineStyle(linestyles.at(iDir));
        
        TFitResultPtr fitResult;
        int fitStatus = -1;
        int nTrials = 0;
        while( (fitStatus != 0) && (nTrials < 10) )
        {
          fitResult = g[iCat] -> Fit(funcName,"QRHNS");
          fitStatus = fitResult;
          if( fitStatus == 1 ) break;
          ++nTrials;
        }
        std::cout << " >>>> fitStatus = " << fitStatus << std::endl;
        std::cout << " >>>> nTrials = " << nTrials << std::endl;
        if( fitStatus == 0 && MCClosure == false )
        {
          TMatrixDSym cov = fitResult->GetCovarianceMatrix();
          TMatrixDSym cor = fitResult->GetCorrelationMatrix();
          for(int aa=0; aa<2; ++aa){
	    for(int bb=0; bb<2; ++bb){
	    std::cout << " (corMatrix[" << aa << "])[" << bb << "] = " << cor[aa][bb] << "; " << std::endl;
            }
          }
        }
        
        std::cout << " FINALE >>> f_fit[iCat]->GetChisquare()/f_fit[iCat]->GetNDF() = " << f_fit[iCat]->GetChisquare()/f_fit[iCat]->GetNDF() << std::endl;
        
        if( drawFitFunc == true )
        {
          hint[iCat] = new TH1F("hint","",5000,20.,1020.);
          (TVirtualFitter::GetFitter()) -> GetConfidenceIntervals(hint[iCat],0.68);
          hint[iCat] -> SetMarkerSize(0);
          hint[iCat] -> SetFillColor(kAzure-9);
          hint[iCat] -> SetFillStyle(3001);
          c[iCat] -> cd();
          f_fit[iCat] -> Draw("same");
          if( iDir == 0 )
          {
            hint[iCat] -> Draw("same,E4");
            f_fit[iCat] -> Draw("same");
	  }
	  g[iCat]->Draw("P,same");
          c_all -> cd(iCat+1);
          f_fit[iCat] -> Draw("same");
          if( iDir == 0 )
          {
            hint[iCat] -> Draw("same,E4");
            f_fit[iCat] -> Draw("same");
	  }
	  g[iCat]->Draw("P,same");
        }
        
        if( MCClosure == false && writeFunctions == true)
	  {
	    TFile TF1_defaultDiffNonLin((("TF1_pol1EoP/TF1_"+fitMethod+"_"+analysis+"_"+nameTrial.at(iDir)+Form("_cat%d.root",iCat))).c_str(), "recreate");
          
	    TF1* defaultDiffNonLinET;
	    TF1* defaultDiffNonLinHT;
	    if(fitMethod == "exp")
	      {
		defaultDiffNonLinET = new TF1(Form("cat%d", iCat), "1.+[0]*(1.-exp(-1.*[1]*(x-45.)))", 20., 1000.);
		defaultDiffNonLinET->SetParameter(0, f_fit[iCat]->GetParameter(0));
		defaultDiffNonLinET->SetParError(0, f_fit[iCat]->GetParError(0));
		defaultDiffNonLinET->SetParameter(1, 2.*f_fit[iCat]->GetParameter(1));
		defaultDiffNonLinET->SetParError(1, 2.*f_fit[iCat]->GetParError(1));
	      }
	    if(fitMethod == "exp3par")
	      {
		defaultDiffNonLinET = new TF1(Form("cat%d", iCat), "1.+[0]*(1.-exp(-1.*[1]*(x-45.)))+[2]", 20., 1000.);
		defaultDiffNonLinET->SetParameter(0, f_fit[iCat]->GetParameter(0));
		defaultDiffNonLinET->SetParError(0, f_fit[iCat]->GetParError(0));
		defaultDiffNonLinET->SetParameter(1, 2.*f_fit[iCat]->GetParameter(1));
		defaultDiffNonLinET->SetParError(1, 2.*f_fit[iCat]->GetParError(1));
		defaultDiffNonLinET->SetParameter(2, f_fit[iCat]->GetParameter(2));
		defaultDiffNonLinET->SetParError(2, f_fit[iCat]->GetParError(2));
	      }
	    if( fitMethod == "pol1" )
	      {
		defaultDiffNonLinET = new TF1(Form("ET_cat%d", iCat), "[0]+ [1] * (x-45.)", 20., 1000.);
		defaultDiffNonLinET->SetParameter(0, f_fit[iCat]->GetParameter(0));
		defaultDiffNonLinET->SetParError(0, f_fit[iCat]->GetParError(0));
		defaultDiffNonLinET->SetParameter(1, f_fit[iCat]->GetParameter(1));
		defaultDiffNonLinET->SetParError(1, f_fit[iCat]->GetParError(1));
		
		defaultDiffNonLinHT = new TF1(Form("HT_cat%d", iCat), "[0]+ [1] * (x-90.)", 20., 1000.);
		defaultDiffNonLinHT->SetParameter(0, f_fit[iCat]->GetParameter(0));
		defaultDiffNonLinHT->SetParError(0, f_fit[iCat]->GetParError(0));
		defaultDiffNonLinHT->SetParameter(1, 0.5*f_fit[iCat]->GetParameter(1));
		defaultDiffNonLinHT->SetParError(1, 0.5*f_fit[iCat]->GetParError(1));
	      }
	    
	    defaultDiffNonLinET -> Write();
	    defaultDiffNonLinHT -> Write();
	    TF1_defaultDiffNonLin.Close();
	  }
        
        legend[iCat] -> AddEntry(g[iCat],Form("cat. %d",iCat),"PL");
        
        
        TH1F* h_Ht_MC = (TH1F*)( f->Get("h_Ht_MC") );
        double HT_Z = h_Ht_MC->GetMean();
        scale_EoP[iCat] = f_fit[iCat] -> Eval(HT_Z);
        
        std::cout << std::fixed << "rel. scale(EoP): "  << std::setprecision(4) << scale_EoP[iCat] << std::endl;
        
        
        c[iCat] -> cd();
        legend[iCat] -> Draw("same");
        latex -> Draw("same");
        c_all -> cd(iCat+1);
        legend[iCat] -> Draw("same");
        latex -> Draw("same");
             
        std::string MCClosureLabel;
        if( MCClosure == true ) MCClosureLabel = "MCClosure_";
        else                    MCClosureLabel = "_";
        
        if( iDir == nDir-1 )
        {
          char pdfName[250];
          sprintf(pdfName,"%s/scale_%s_%s_%s_cat%d_%sDAOverMC.pdf",
                  directories[0].c_str(),analysis.c_str(),(methods.at(iMeth)).c_str(),fitMethod.c_str(),iCat,MCClosureLabel.c_str());
	  if(MCClosure == true)           sprintf(pdfName,"%s_MCClosure/scale_%s_%s_%s_cat%d_%sDAOverMC.pdf",
                  directories[0].c_str(),analysis.c_str(),(methods.at(iMeth)).c_str(),fitMethod.c_str(),iCat,MCClosureLabel.c_str());
          std::cout << ">>> saving file " << pdfName << std::endl;
          c[iCat] -> Print(pdfName,"pdf");
          
          char pngName[250];
          sprintf(pngName,"%s/scale_%s_%s_%s_cat%d_%sDAOverMC.png",
                  directories[0].c_str(),analysis.c_str(),(methods.at(iMeth)).c_str(),fitMethod.c_str(),iCat,MCClosureLabel.c_str());
	  if(MCClosure == true)           sprintf(pngName,"%s_MCClosure/scale_%s_%s_%s_cat%d_%sDAOverMC.png",
                  directories[0].c_str(),analysis.c_str(),(methods.at(iMeth)).c_str(),fitMethod.c_str(),iCat,MCClosureLabel.c_str());
          std::cout << ">>> saving file " << pngName << std::endl;
          c[iCat] -> Print(pngName,"png");
          
          if( iCat == nCat-1 )
          {
            sprintf(pdfName,"%s/scale_%s_%s_%s_allCat_%sDAOverMC.pdf",
                    directories[0].c_str(),analysis.c_str(),(methods.at(iMeth)).c_str(),fitMethod.c_str(),MCClosureLabel.c_str());
	  if(MCClosure == true)           sprintf(pdfName,"%s_MCClosure/scale_%s_%s_%s_allCat_%sDAOverMC.pdf",
                  directories[0].c_str(),analysis.c_str(),(methods.at(iMeth)).c_str(),fitMethod.c_str(),MCClosureLabel.c_str());
            std::cout << ">>> saving file " << pdfName << std::endl;
            c_all -> Print(pdfName,"pdf");
            
            sprintf(pngName,"%s/scale_%s_%s_allCat_%sDAOverMC.png",
                    directories[0].c_str(),analysis.c_str(),(methods.at(iMeth)).c_str(),fitMethod.c_str(),MCClosureLabel.c_str());
	  if(MCClosure == true)           sprintf(pngName,"%s_MCClosure/scale_%s_%s_%s_allCat_%sDAOverMC.png",
                  directories[0].c_str(),analysis.c_str(),(methods.at(iMeth)).c_str(),fitMethod.c_str(),MCClosureLabel.c_str());
            std::cout << ">>> saving file " << pngName << std::endl;
            c_all -> Print(pngName,"png"); 
          }
        }
      }

      TFile perTommaso("perTommaso_EoP_Etdep.root","recreate");
      for(int i =0; i<4; ++i){
	g[i]->Write(Form("graph_cat%d",i));
      }
      perTommaso.Close();

    }
  }    
}
  
