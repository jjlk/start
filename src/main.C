/*
  Main program
*/

// STL
#include <iostream>
#include <vector>

// START
#include "Config.hh"
#include "BandsFactory.hh"
#include "Band.hh"
#include "DataSummary.hh"
#include "HandleResolArea.hh"
#include "MinimizeFactory.hh"
#include "ResidualsFactory.hh"
#include "Hypothesis.hh"
#include "PlotFactory.hh"
#include "PowerLaw.hh"
#include "ExpoCutOffPowerLaw.hh"
#include "LogParabolic.hh"
#include "BrokenPowerLaw.hh"

// Utilities
#define DEBUG 0
#include "debugging.hh"

// ROOT
#include <TString.h>
#include <TFile.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TRint.h>
#include <TRandom3.h>

// TMP
#include <TGraph.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <map>
#include <utility>
#include <TH2.h>

/**
 * \brief Main Program
 *
 * Do minimization and plot spectrum and residuals (also contours if specified).
 * Hypothesis are automatically rebinned at 2 sigmas
 *
 * \param fileconfig config file
 * \param savefile name + path of a rootfile where canvas will be saved
 * \param if true contours for phi and gamma are computed
 */
int main(int argc, char **argv) {
	TRandom3 rg(0);
  
	for (int i=0;i<10;++i) {
		std::cout << "\033[3" << Int_t(rg.Uniform(0.,10.)) << ";1m";
		std::cout << "******************************************************************************************" << std::endl;
		std::cout << "********************* DEPRECATED use scripts/startfit.C instead !!!! *********************" << std::endl;
		std::cout << "******************************************************************************************" ;
		std::cout << "\033[0m" << std::endl;
	}
	return 0;
  
	INFO_OUT << "START is launched!" << std::endl;

	if(argc<2) {
		INFO_OUT << "How to use it :" << std::endl;
		INFO_OUT << "param1 is the path + the config file name (ex: \"./config/myconfig.config\" for example)" << std::endl;
		INFO_OUT << "param2 is the name (not the path) of the rootfile where canvas will be saved (ex: \"bob.root\")" << std::endl;
		INFO_OUT << "param3 is \"true\" if you want contours and something else if not" << std::endl;
		return 0;
	}

	TString fileconfig = argv[1];
	TString savefile = argv[2];
	TString docontours = argv[3];

	// fileconfig="conffile/Crab.config";
	// savefile="Crab.root";
	// Creation of an object Config. It will contain all infos for the minimization & Co

	START::Config myConfig(fileconfig);

	// Creation of a vector myData which will contain all the data in object Bands

	std::vector<START::Band> myData;

	// Get IrfOpt
	std::string StrIrfAnalysis(myConfig.GetUserIrfOpt());
	START::STARTUtils::IrfOpt IrfType;
	if(StrIrfAnalysis=="HapFr_Hess_Stereo")
		IrfType = START::STARTUtils::HapFr_Hess_Stereo;
	else if(StrIrfAnalysis=="HapFr_Hess_Hybrid")
		IrfType = START::STARTUtils::HapFr_Hess_Hybrid;
	else if(StrIrfAnalysis=="HapFr_Hess_Mono")
		IrfType = START::STARTUtils::HapFr_Hess_Mono;
	else {
		std:: cout << " - HapFr_Hess_Stereo" << std::endl;
		std:: cout << " - HapFr_Hess_Hybrid" << std::endl;
		std:: cout << " - HapFr_Hess_Mono" << std::endl;
		exit(0);
	}


	// Creation of an object BandsFactory which will copy the data in the object myData
	START::BandsFactory BandsFact(myConfig, IrfType);
	BandsFact.MakeBands(myData);
	BandsFact.PrintBands(myData);

	// Creation of an object DataSummary which contain general informations on the bands
	START::DataSummary DataSum;

	/// Creation of an object HandleResolArea wich will store the effective areas and 
	/// resolutions useful to determine the energy threshold of the bands
	TString conf_for_area = (myConfig.GetUserForceAreaAnalysisConfig().CompareTo("")==0) ? myConfig.GetUserAnalysisConfig() : myConfig.GetUserForceAreaAnalysisConfig();
	bool usetruedistribution = (conf_for_area.Contains("thsq64"))  ? true : false;
	START::HandleResolArea *HandleIRF = 0;
	HandleIRF = new START::HandleResolArea(myConfig,IrfType); // JLK change
	//Handle->SetESafeThresholdCondition(HandleResolArea::AreaMaxMethod,5,0.); //paris analysis method
	HandleIRF->SetESafeThresholdCondition(START::HandleResolArea::AreaAndBiaisResolMethod,myConfig.GetUserAreaMin(),2.);
	HandleIRF->SetEnergyBinThresholdCondition(START::HandleResolArea::Safe);

	//print infos
	DataSum.PrintDataInfo(myData);
  
	// Declare one or more hypothesis
	std::vector<START::Hypothesis*> HypothesisArray;
	HypothesisArray.push_back(new START::PowerLaw("PWL1",START::Hypothesis::Differential,1.e-11,2.0));
	//HypothesisArray.push_back(new ExpoCutOffPowerLaw("EXP1"));
	//HypothesisArray.push_back(new LogParabolic("LOG1"));
	//HypothesisArray[0]->FixParameter("phi0",3.5e-11,true);
	//HypothesisArray[0]->FixParameter("Gamma",3.,true);
	//HypothesisArray[0]->FixParameter("beta",1/16.5,false);

	// Minimization
	START::MinimizeFactory *Minimizer = new START::MinimizeFactory(HypothesisArray, myConfig);
	Minimizer->MakeMinimization(myData);

	if(docontours=="true") { // do contours for phi0 against gamma
		for (unsigned int icont = 0; icont<HypothesisArray.size(); icont++) {
			//Minimizer->MakeContours(myData,HypothesisArray[icont]->GetName(),0,1,30); // contours phi0/gamma
			Minimizer->MakeContours(myData,HypothesisArray[icont]->GetName(),0,1,5); // contours phi0/gamma
		}
	}
	//Minimizer->MakeAllContours(myData,20); // all contours for all pairs of params are computed but it's very long...

	//Minimizer->ScanLikelihoodAfterMinimization(*HypothesisArray[1],myData,"Gamma",40);

	// Residuals
	START::ResidualsFactory *Residus = new START::ResidualsFactory(HypothesisArray);
	Residus->ComputeResidualsRolke(); // Compute residuals and stock them in the hypothesis

	// Plot
	START::PlotFactory *Plotter = 0;
	Plotter = new START::PlotFactory(HypothesisArray);

	if(docontours=="true") Plotter->PlotContours();
	//Plotter->PlotSpectrumAndResiduals(PlotFactory::LinearButterfly);
	Plotter->PlotSpectrumAndResiduals(START::PlotFactory::LogarithmButterfly);
	Plotter->PlotLikelihoodScans();
	for (unsigned int ireb = 0; ireb<HypothesisArray.size(); ++ireb) { // rebin hypothesis
		Plotter->RebinHypothesis(*HypothesisArray[ireb],2.0);
	}

	TString pathrootfile = myConfig.GetUserOutputFolderName();
	pathrootfile+="/";
	pathrootfile+=savefile;
	Plotter->SaveObjectsInROOTFile(pathrootfile); // Save objects in ROOT file

	delete HandleIRF;
	delete Minimizer;
	delete Residus;
	delete Plotter;
  
	for (unsigned int i = 0 ; i<HypothesisArray.size(); ++i) { // clean
		delete HypothesisArray[i];
	}
  
	return 0;
}

/*
  myData[15].PrintAll();
  
  TCanvas *canvas = new TCanvas("canvas","canvas");
  canvas->SetFillColor(0);

  double emin(0.02), emax(130.);

  std::vector<TString> nameconf;
  nameconf.push_back(Handle->MakeKeyRun(60.,0.,0.,""));
  nameconf.push_back(Handle->MakeKeyRun(60.,0.,0.5,""));
  nameconf.push_back(Handle->MakeKeyRun(60.,18.,0.5,""));
  nameconf.push_back(Handle->MakeKeyRun(60.,18.,0.,""));
  nameconf.push_back(Handle->MakeKeyRun(70,0.,0.,""));
  nameconf.push_back(Handle->MakeKeyRun(70.,0.,0.5,""));
  nameconf.push_back(Handle->MakeKeyRun(70.,18.,0.5,""));
  nameconf.push_back(Handle->MakeKeyRun(70.,18.,0.,""));

  std::map<TString, std::pair<std::vector<double >, std::vector<double > > > map = Handle->GetMapVectorsArea();

  if(map.size()==0) {
  std::cout << "empty ==> exit" << std::endl;
  exit(0);
  }

  for(int i(0); i<nameconf.size(); i++) std::cout << nameconf[i] << std::endl;

  std::vector<TGraph*> graph;

  TMultiGraph *multi = new TMultiGraph("multi","multi");

  //TH2D *histo = new TH2D("histo","Effective area for a band",100,0.,125.,100,200.,1.e6);
  //histo->Draw("AXIS");
  for(int iconf(0); iconf<nameconf.size(); iconf++) {

  int count(0);

  std::vector<double> x = map[nameconf[iconf]].first;
  std::vector<double> y = map[nameconf[iconf]].second;

  std::cout << "Remplissage de" << nameconf[iconf] << std::endl; 

  TString graphname = "graph_";
  graphname+=iconf;
  //3 4 6 7
  graph.push_back(new TGraph(x.size()));
  graph[iconf]->SetName(graphname);
  graph[iconf]->SetTitle(nameconf[iconf]);
  if(nameconf[iconf].Contains("eff60")) {
  graph[iconf]->SetMarkerStyle(22);
  if(nameconf[iconf].Contains("zen0") && nameconf[iconf].Contains("off0.0"))
  graph[iconf]->SetMarkerColor(3);
  else if(nameconf[iconf].Contains("zen0") && nameconf[iconf].Contains("off0.5"))
  graph[iconf]->SetMarkerColor(4);
  else if(nameconf[iconf].Contains("zen18") && nameconf[iconf].Contains("off0.0"))
  graph[iconf]->SetMarkerColor(6);
  else graph[iconf]->SetMarkerColor(7);
  }
  else if(nameconf[iconf].Contains("eff70")){
  graph[iconf]->SetMarkerStyle(23);
  if(nameconf[iconf].Contains("zen0") && nameconf[iconf].Contains("off0.0"))
  graph[iconf]->SetMarkerColor(3);
  else if(nameconf[iconf].Contains("zen0") && nameconf[iconf].Contains("off0.5"))
  graph[iconf]->SetMarkerColor(4);
  else if(nameconf[iconf].Contains("zen18") && nameconf[iconf].Contains("off0.0"))
  graph[iconf]->SetMarkerColor(6);
  else graph[iconf]->SetMarkerColor(7);
  }

  graph[iconf]->SetMarkerSize(1);
  graph[iconf]->SetDrawOption("P");
  graph[iconf]->SetFillStyle(0);
  graph[iconf]->GetXaxis()->SetLimits(emin,emax);
  graph[iconf]->SetFillColor(0);

  for(int i(0); i<x.size(); i++) {
      
  if(x.size()!=y.size()) exit(0);
      
  if(y[i]==0) continue;

  graph[iconf]->SetPoint(count,TMath::Log10(x[i]),TMath::Log10(y[i]));
  //graph[iconf]->SetPoint(count,x[i],y[i]);
      
  count++;
      
  }
    
  x.clear(); y.clear(); 
  
  }

  //for(int iconf(0); iconf<nameconf.size(); iconf++) graph[iconf]->Draw("P");

  std::vector<double> xband, yband; 

  int n=100;

  for(double ien(TMath::Log10(emin)); ien<TMath::Log10(emax); ien+=(TMath::Log10(emax)-TMath::Log10(emin))/n) { 
  double area = myData[15].GetInterpolatedArea(TMath::Power(10,ien));
  if(area<10.) continue;
  xband.push_back(TMath::Power(10,ien));
  yband.push_back(area);
  }
  
  TGraph *gband = new TGraph();
  gband->SetName("gband");
  TString bandname = "AreaBand_eff";
  char eff[10];
  sprintf(eff,"%.2f",myData[15].GetEff());
  bandname+=eff;
  bandname+="_zen";
  char zen[10];
  sprintf(zen,"%.2f",myData[15].GetZenON());
  bandname+=zen;
  bandname+="_off";
  char off[10];
  sprintf(off,"%.2f",myData[15].GetOffset());
  bandname+=off;

  gband->SetTitle(bandname);
  gband->SetMarkerStyle(7);
  gband->SetMarkerColor(1);
  gband->SetLineColor(2);
  gband->SetFillStyle(0);
  gband->SetFillColor(0);
  gband->GetXaxis()->SetTitle("bob");
  gband->GetYaxis()->SetTitle("bob");

  for(int ip(0); ip<xband.size(); ip++) {
  if(yband[ip]<10) continue;
  gband->SetPoint(ip,TMath::Log10(xband[ip]),TMath::Log10(yband[ip]));
  //gband->SetPoint(ip,xband[ip],yband[ip]);
  }

  //multi->Add(gband);
  for(int iconf(0); iconf<nameconf.size(); iconf++) multi->Add(graph[iconf]);//,"P");
  gband->Draw("APL");
  multi->Draw("P");

  TFile *file = new TFile("/Users/julien/Desktop/These/Tools/hap-11-02/analysis/bin/canvas.root","RECREATE");
  canvas->BuildLegend();
  canvas->Draw();
  canvas->Update();
  canvas->Write();
  //multi->Write();
  file->Close();

*/
