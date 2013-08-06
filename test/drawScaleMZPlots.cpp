// g++ -Wall -o drawScaleMZPlots `root-config --cflags --glibs` setTDRStyle.cc drawScaleMZPlots.cpp

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


bool MCClosure   = false;
bool drawFitFunc = true;

std::string analysis  = "CiC";
std::string fitMethod = "exp";
std::string directory = "data/2012/Moriond2013/CiC_nonGlobe_eleTunedReg_Dphi3p15/";



int main(int argc, char** argv)
{
  //Check if all nedeed arguments to parse are there
  if(argc != 5)
  {
    std::cerr << ">>>drawScaleMZPlots need configFiles for every category." << std::endl;
    return -1;
  }



  //----------------------
  // Parse the config file

  int category[4], evtsPerPoint[4];
  for( int cat = 0; cat < 4; cat++ )
  { 
    parseConfigFile(argv[cat+1]);
    category[cat]        = gConfigParser -> readIntOption("Options::category");
    evtsPerPoint[cat]    = gConfigParser -> readIntOption("Options::evtsPerPoint");
  }


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
    
  std::vector<std::string> methods;
  //methods.push_back( "fit" );
  //methods.push_back( "gausFit" );
  //methods.push_back( "mean" );
  methods.push_back( "recursiveMean" );
  //methods.push_back( "smallestInterval" );
  int nMethods = methods.size();
  
  std::vector<int> colors;
  colors.push_back(kBlue+0);
  colors.push_back(kRed+0);
  colors.push_back(kBlue+1);
  colors.push_back(kRed+1);

  std::vector<int> linestyles;
  linestyles.push_back(1);
  linestyles.push_back(1);
  linestyles.push_back(1);
  linestyles.push_back(1);  
  
  
  std::string* inFileNames     = new std::string[4];
  
  std::string baseFileName = "studyLinearity_MZ_";
  if( MCClosure == true) baseFileName += "_MCClosure";

  TCanvas** c01 = new TCanvas*[2];
  TCanvas** c23 = new TCanvas*[2];
  
  
  for(int iMeth = 0; iMeth < nMethods; ++iMeth)
  {
    std::cout << ">>> method: " << methods.at(iMeth) << std::endl;
    
    
    TLegend* legend01 = new TLegend(0.16, 0.77, 0.30, 0.92);
    legend01 -> SetFillColor(kWhite);
    legend01 -> SetFillStyle(1001);  
    legend01 -> SetTextFont(42);  
    legend01 -> SetTextSize(0.05);
    
    TLegend* legend23 = new TLegend(0.16, 0.77, 0.30, 0.92);
    legend23 -> SetFillColor(kWhite);
    legend23 -> SetFillStyle(1001);  
    legend23 -> SetTextFont(42);  
    legend23 -> SetTextSize(0.05);    
    
    
    char latexText[50];
    //sprintf(latexText,"scale - %s",methods.at(iMeth).c_str());
    sprintf(latexText,"analysis: %s   scale estimator: %s",analysis.c_str(),(methods.at(iMeth)).c_str());
    TLatex* latex = new TLatex(0.14,0.96,latexText);                                                                                                                              
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.04);
    
    
    
    double* scale_MZ = new double[4];
    //double* syst     = new double[4];
    //double* systErr  = new double[4];
    TGraphAsymmErrors** g = new TGraphAsymmErrors*[4];
    
    TF1** f_fit = new TF1*[4];
    TH1F** hint = new TH1F*[4];
    
    for(int cat = 0; cat < 4; ++cat)
    {
      std::cout << "\n*****************" << std::endl;
      char EvtString[50];
      sprintf(EvtString,"%devtsPerPoint_cat%d",evtsPerPoint[cat],category[cat]);
      
      inFileNames[cat] = directory + "/" + baseFileName + std::string(EvtString) + ".root";
      TFile* f     = TFile::Open(     (inFileNames[cat]).c_str(),"READ" );
      //std::cout << ">>> inFileName: " << inFileNames[cat] << std::endl;
      
      char graphName[50];
      sprintf(graphName,"step1/scale_%s_DAOverMC",methods.at(iMeth).c_str());
      g[cat] = (TGraphAsymmErrors*)( f->Get(graphName) );
      
      
      if( cat == 0 )
      {
        c01[iMeth] = new TCanvas();
        c01[iMeth] -> cd();
        c01[iMeth] -> SetGridx();
        c01[iMeth] -> SetGridy();
        //c01[iMeth] -> SetLogx();
      }
      if( cat == 1 )
      {
        c01[iMeth] -> cd();
      }
      if( cat == 2 )
      {
        c23[iMeth] = new TCanvas();
        c23[iMeth] -> cd();
        c23[iMeth] -> SetGridx();
        c23[iMeth] -> SetGridy();
        //c23[iMeth] -> SetLogx();
      }
      if( cat == 3 )
      {
        c23[iMeth] -> cd();
      }
      
      g[cat] -> SetPoint(g[cat]->GetN(),1500.,1.);
      g[cat] -> GetXaxis() -> SetRangeUser(60,200.);        
      g[cat] -> GetXaxis() -> SetMoreLogLabels();
      g[cat] -> GetXaxis() -> SetTitle("H_{T} [GeV]");
      g[cat] -> GetYaxis() -> SetTitle("scale_{data} / scale_{MC}");
      g[cat] -> GetXaxis() -> SetTitleSize(0.05);
      g[cat] -> GetYaxis() -> SetTitleSize(0.05);
      g[cat] -> GetXaxis() -> SetTitleOffset(1.10);
      g[cat] -> GetYaxis() -> SetTitleOffset(1.10);
      g[cat] -> GetXaxis() -> SetLabelSize(0.04);
      g[cat] -> GetYaxis() -> SetLabelSize(0.04);
      g[cat] -> SetMarkerSize(0.7);
      g[cat] -> SetMarkerStyle(20+cat%2);
      g[cat] -> SetMarkerColor(colors.at(cat));
      g[cat] -> SetLineColor(colors.at(cat));
      g[cat] -> SetLineStyle(linestyles.at(cat));
      
      g[cat] -> SetMinimum(0.980);
      g[cat] -> SetMaximum(1.020);
      if( drawFitFunc == false )
      {
        if( cat == 0 || cat == 2 ) g[cat] -> Draw("APL");
        else                       { g[cat] -> Draw("PL,same"); g[cat-1] -> Draw("PL,same"); }
      }
      else
      {
        if( cat == 0 || cat == 2 ) g[cat] -> Draw("AP");
        else                       { g[cat] -> Draw("P,same"); g[cat-1] -> Draw("P,same"); }
      }
      
      char funcName[150];
      sprintf(funcName,"f_prefit_%s_%d",methods.at(iMeth).c_str(),cat);
      TF1* f_prefit;
      if( fitMethod == "pol1")
      {
        f_prefit = new TF1(funcName,"[0]+[1]*x",65.,210.);
        f_prefit -> SetParameters(1.,1.);         
      }
      if( fitMethod == "exp")
      {
        f_prefit = new TF1(funcName,"1.+[0]*(1.-exp(-1.*[1]*(x-90.)))+[2]",65.,210.);
        f_prefit -> SetParameters(0.005,0.02,0.);
        f_prefit -> SetParLimits(0,0.,1.);
        f_prefit -> SetParLimits(1,0.,0.05);
      }
      g[cat] -> Fit(funcName,"QERNS");
      
      
      if( drawFitFunc == true )
      {
        for(int point = 0; point < g[cat]->GetN(); ++point)
        {
          double ey = g[cat] -> GetErrorY(point);
          g[cat] -> SetPointEYhigh(point,ey*sqrt(f_prefit->GetChisquare()/f_prefit->GetNDF()));
          g[cat] -> SetPointEYlow (point,ey*sqrt(f_prefit->GetChisquare()/f_prefit->GetNDF()));
        }
      }
      
      
      sprintf(funcName,"f_fit_%s_%d",methods.at(iMeth).c_str(),cat);
      if( fitMethod == "pol1" )
      {
        f_fit[cat] = new TF1(funcName,"[0]+[1]*x",65.,210.);
        f_fit[cat] -> SetParameters(1.,1.);
      }
      if( fitMethod == "exp")
      {
        f_fit[cat] = new TF1(funcName,"1.+[0]*(1.-exp(-1.*[1]*(x-90.)))+[2]",65.,210.);
        f_fit[cat] -> SetParameters(0.005,0.02,0.);
        f_fit[cat] -> SetParLimits(0,-1.,1.);
        f_fit[cat] -> SetParLimits(1,0.,0.05);
      }
      f_fit[cat] -> SetLineColor(colors.at(cat));
      f_fit[cat] -> SetLineWidth(2);
      f_fit[cat] -> SetLineStyle(linestyles.at(cat));
      
      TFitResultPtr fitResult;
      int fitStatus = -1;
      int nTrials = 0;
      
      while( (fitStatus != 1) && (nTrials < 2) )
      {
        fitResult = g[cat] -> Fit(funcName,"QENRS+");
        fitStatus = fitResult;
        if( fitStatus == 1 ) break;
        ++nTrials;
      }
      
      

      if( drawFitFunc == true )
      {
        f_fit[cat] -> Draw("same");
        hint[cat] = new TH1F("hint","",1000,65.,1065.);
        (TVirtualFitter::GetFitter()) -> GetConfidenceIntervals(hint[cat],0.68);
        hint[cat] -> SetMarkerSize(0);
        hint[cat] -> SetFillColor(colors.at(cat));
        hint[cat] -> SetFillStyle(3001);
        hint[cat] -> Draw("same,E4");
      }  
      
      
            
      char legendText[50];
      sprintf(legendText,"cat. %d",cat);
      if( cat == 0 || cat == 1) legend01 -> AddEntry(g[cat],legendText,"PL");
      else                      legend23 -> AddEntry(g[cat],legendText,"PL");
      
      
      
      TH1F* h_Ht_MC = (TH1F*)( f->Get("h_Ht_MC") );
      double HT_Z = h_Ht_MC->GetMean();
      scale_MZ[cat] = f_fit[cat] -> Eval(HT_Z);
      
      std::cout << ">>>>>> cat: " << cat
                << std::fixed
                << "   rel. scale(MZ): "  << std::setprecision(4) << scale_MZ[cat]
                << std::endl;
    }
    
    c01[iMeth] -> cd();
    legend01 -> Draw("same");
    latex -> Draw("same");
    c23[iMeth] -> cd();
    legend23 -> Draw("same");
    latex -> Draw("same");
    
    
    char pdfName[250];
    char pngName[250];
    
    if( MCClosure == false )
    {
      sprintf(pdfName,"%s/scale_%s_%s_%s_cat01_DAOverMC.pdf",directory.c_str(),analysis.c_str(),(methods.at(iMeth)).c_str(),fitMethod.c_str());
      std::cout << ">>> saving file " << pdfName << std::endl;
      c01[iMeth] -> Print(pdfName,"pdf");
      sprintf(pngName,"%s/scale_%s_%s_%s_cat01_DAOverMC.png",directory.c_str(),analysis.c_str(),(methods.at(iMeth)).c_str(),fitMethod.c_str());
      std::cout << ">>> saving file " << pngName << std::endl;
      c01[iMeth] -> Print(pngName,"png");
      
      sprintf(pdfName,"%s/scale_%s_%s_%s_cat23_DAOverMC.pdf",directory.c_str(),analysis.c_str(),(methods.at(iMeth)).c_str(),fitMethod.c_str());
      std::cout << ">>> saving file " << pdfName << std::endl;
      c23[iMeth] -> Print(pdfName,"pdf");
      sprintf(pngName,"%s/scale_%s_%s_%s_cat23_DAOverMC.png",directory.c_str(),analysis.c_str(),(methods.at(iMeth)).c_str(),fitMethod.c_str());
      std::cout << ">>> saving file " << pngName << std::endl;
      c23[iMeth] -> Print(pngName,"png");
    }
    else
    {
      sprintf(pdfName,"%s/scale_%s_%s_%s_cat01_MCClosure_DAOverMC.pdf",directory.c_str(),analysis.c_str(),(methods.at(iMeth)).c_str(),fitMethod.c_str());
      std::cout << ">>> saving file " << pdfName << std::endl;
      c01[iMeth] -> Print(pdfName,"pdf");
      sprintf(pngName,"%s/scale_%s_%s_%s_cat01_MCClosure_DAOverMC.png",directory.c_str(),analysis.c_str(),(methods.at(iMeth)).c_str(),fitMethod.c_str());
      std::cout << ">>> saving file " << pngName << std::endl;
      c01[iMeth] -> Print(pngName,"png");

      sprintf(pdfName,"%s/scale_%s_%s_%s_cat23_MCClosure_DAOverMC.pdf",directory.c_str(),analysis.c_str(),(methods.at(iMeth)).c_str(),fitMethod.c_str());
      std::cout << ">>> saving file " << pdfName << std::endl;
      c23[iMeth] -> Print(pdfName,"pdf");
      sprintf(pngName,"%s/scale_%s_%s_%s_cat23_MCClosure_DAOverMC.png",directory.c_str(),analysis.c_str(),(methods.at(iMeth)).c_str(),fitMethod.c_str());
      std::cout << ">>> saving file " << pngName << std::endl;
      c23[iMeth] -> Print(pngName,"png");
    }
  }
}




     //if( fitMethod == "exp" )
      //{
      //  std::cout << "qui11" << std::endl;
      //  f_fit[cat] = new TF1(funcName,"[0]+[1]*x+[2]",65.,200.);
      //  //f_fit[cat] = new TF1(funcName,"[0]*( 1-exp(-(x-[1])/[2]))",65.,200.);
      //  std::cout << "qui12" << std::endl;
      //  f_fit[cat] -> SetParameters(1.01,10.,10.);
      //  std::cout << "qui13<" << std::endl;
      //}
 
