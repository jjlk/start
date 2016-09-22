// STL
#include <iostream>
#include <algorithm>
#include <sstream>

// ROOT
#include <TObject.h>
#include <TString.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH2.h>
#include <TFile.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TGraphPainter.h>
#include <TMultiGraph.h>
#include <TGraphAsymmErrors.h>
#include <TPad.h>
#include <TArrow.h>
#include <TPaveText.h>
#include <TMath.h>
#include <TLatex.h>
#include <TDirectory.h>
#include <TROOT.h>
#include <TMarker.h>
#include <TDatime.h>

// START
#include "PlotFactory.hh"
#include "ResidualsFactory.hh"
#include "Band.hh"
#include "TimeBinVector.hh"
#include "Hypothesis.hh"
#include "STARTUtils.hh"
#include "Residuals.hh"

// Utilities
#define DEBUG 1
#include "debugging.hh"

#define INFO std::cout << INFOCOLOR << "PlotFactory> " << RESETCOLOR
#define WARNING std::cout << WARNINGCOLOR << "PlotFactory> " << RESETCOLOR

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::PlotFactory)
#endif


/**
 * \brief Constructor to do a MWL plot
 */
START::PlotFactory::PlotFactory(std::string sourcename, PlotStyle Style, ResidualsStyle ResStyle, bool verbose)
:fverbose(verbose),
	fsourcename(sourcename),
	fdrawfit(false),
	fminimalinfosonpad(true),
	fSaveRootFile(0),
	fscalearrowfactor(0.05),
	fLegendMwl(0)
{

	fHypothesisArray.clear();
	fRebinnedHypothesisArray.clear();
	fPlotStyle = Style;
	fResidualsStyle = ResStyle;

	fresiduals.clear();
	fresiduals_on.clear();
	fresiduals_off.clear();
	fresidualssigmaplus.clear();
	fresidualssigmaminus.clear();
	fresiduals3sigmaplus.clear();
	fresiduals3sigmaminus.clear();
	fflux.clear();
	ffluxsigmaplus.clear();
	ffluxsigmaminus.clear();
	fflux3sigmaplus.clear();
	fflux3sigmaminus.clear();
	fenergybin.clear();

	CleanCanvasMwl();
	CleanPadMwl();
	CleanTH2Mwl();
	CleanGraphErrorsMwl();
	CleanTF1Mwl();
	CleanArrowsMwl();
	CleanButterflyMwl();

}

/**
 * \brief Constructor for an hypothesis' vector
 *
 * \param Hypopthesis array
 */
START::PlotFactory::PlotFactory(std::vector<Hypothesis*> &hypo, std::string sourcename,
								PlotStyle Style, ResidualsStyle ResStyle, bool verbose)
	:fverbose(verbose),
	 fsourcename(sourcename),
	 fdrawfit(false),
	 fminimalinfosonpad(true),
	 fSaveRootFile(0),
	 fscalearrowfactor(0.05),
	 fLegendMwl(0)
{

	fHypothesisArray.clear();

	for(std::vector<Hypothesis*>::iterator it=hypo.begin(); it!=hypo.end(); ++it) {
		fHypothesisArray.push_back(&(**it));
	}
  
	fPlotStyle = Style;
	fResidualsStyle = ResStyle;
  
	fresiduals.clear();
	// JLK ADD FOR ADA
	fresiduals_on.clear();
	fresiduals_off.clear();
	fresidualssigmaplus.clear();
	fresidualssigmaminus.clear();
	fresiduals3sigmaplus.clear();
	fresiduals3sigmaminus.clear();
	fflux.clear();
	ffluxsigmaplus.clear();
	ffluxsigmaminus.clear();
	fflux3sigmaplus.clear();
	fflux3sigmaminus.clear();
	fenergybin.clear();
	fRebinnedHypothesisArray.clear();
}

/**
 * \brief Constructor for an hypothesis
 *
 * \param Hypopthesis
 */
START::PlotFactory::PlotFactory(Hypothesis &hypo, std::string sourcename,
								PlotStyle Style, ResidualsStyle ResStyle,
								bool verbose)
	:fverbose(verbose),
	 fsourcename(sourcename),
	 fdrawfit(false),
	 fminimalinfosonpad(true),
	 fSaveRootFile(0),
	 fscalearrowfactor(0.05)
{

	fHypothesisArray.clear();
	fHypothesisArray.push_back(&hypo);

	fPlotStyle = Style;
	fResidualsStyle = ResStyle;
  
	fresiduals.clear();
	fresiduals_on.clear();
	fresiduals_off.clear();
	fresidualssigmaplus.clear();
	fresidualssigmaminus.clear();
	fresiduals3sigmaplus.clear();
	fresiduals3sigmaminus.clear();
	fflux.clear();
	ffluxsigmaplus.clear();
	ffluxsigmaminus.clear();
	fflux3sigmaplus.clear();
	fflux3sigmaminus.clear();
	fenergybin.clear();

	fRebinnedHypothesisArray.clear();

}

/**
 * \brief Destructor
 */
START::PlotFactory::~PlotFactory()
{
  
	if(fSaveRootFile!=0) delete fSaveRootFile;
	fSaveRootFile = 0;

	// JLK: I don't delete the pointers because I don't know
	// how to keep th plots hanging after ROOT finish his jobs
	/*
	// spectrum
	CleanPolyLines();
	CleanGraphSpectrum();
	CleanTF1FitSpectrum();
	CleanTH2();
	CleanTF1Residuals();
	CleanPads();
	CleanCanvasSpectrum();
	// contours
	CleanMultiGraphs();
	CleanCanvasContours();
	//Likelihood scans
	CleanGraphScanLikelihood();
	CleanCanvasLikelihoodScan();
	CleanGraphArraySpectrum();
	// Light curves
	CleanGraphLightCurves();
	CleanCanvasLightCurves();
	CleanTH2LightCurves();
	CleanPadLightCurves();
	CleanTF1MeanFluxLightCurves();
	CleanTF1ZeroFluxLightCurves();
	*/
	for (std::vector<Hypothesis*>::iterator ithypo = fRebinnedHypothesisArray.begin();ithypo!=fRebinnedHypothesisArray.end();++ithypo) {
		delete *ithypo;
		*ithypo=0;
	}

}


/**
 * \brief Save objects in root file 'rootfile'
 *
 * \param rootfile file where ROOT objects are saved
 */
void START::PlotFactory::SaveObjectsInROOTFile(TString rootfile)
{

	fSaveRootFile = new TFile(rootfile,"RECREATE");

	SaveSpectrumAndResidualsROOTObjects();

	SaveLightCurvesROOTObjects();

	SaveContoursROOTObjects();

	SaveLikelihoodScansROOTObjects();

	fSaveRootFile->Close();

	INFO << "Saving objects in rootfile : "<< rootfile << "... ok" << std::endl;
  
	for(std::map<std::string,TCanvas*>::iterator canvas=fMapCanvasSpectrumAndResiduals.begin(); canvas!=fMapCanvasSpectrumAndResiduals.end(); canvas++) {
		DEBUG_OUT << "canvas " << canvas->first << " exists" << std::endl;
	}

}


/**
 * \brief Save spectrum and residuals canvas in rootfile
 */
void START::PlotFactory::SaveSpectrumAndResidualsROOTObjects() {

	fSaveRootFile->cd();

	for(std::map<std::string,TCanvas*>::iterator canvas=fMapCanvasSpectrumAndResiduals.begin(); canvas!=fMapCanvasSpectrumAndResiduals.end(); canvas++) {
		DEBUG_OUT << "canvas " << canvas->first << " exists" << std::endl;
		if(canvas->second!=0) canvas->second->Write();
	}

}

/**
 * \brief Save spectrum and residuals canvas in rootfile
 */
void START::PlotFactory::SaveContoursROOTObjects() {

	fSaveRootFile->cd();

	for(std::map<std::string,TCanvas*>::iterator canvas=fMapCanvasContours.begin(); canvas!=fMapCanvasContours.end(); canvas++) {
		DEBUG_OUT << "canvas " << canvas->first << " exists" << std::endl;
		if(canvas->second!=0) canvas->second->Write();
	}

}

/**
 * \brief Save spectrum and residuals canvas in rootfile
 */
void START::PlotFactory::SaveLikelihoodScansROOTObjects() {

	fSaveRootFile->cd();

	for(std::map<std::string,TCanvas*>::iterator canvas=fMapCanvasLikelihoodScan.begin(); canvas!=fMapCanvasLikelihoodScan.end(); canvas++) {
		DEBUG_OUT << "canvas " << canvas->first << " exists" << std::endl;
		if(canvas->second!=0) canvas->second->Write();
	}

}

/**
 * \brief Plot main canvas with spectrum and residuals
 */
void START::PlotFactory::PlotSpectrumAndResiduals(ButterflyOption butterfly, bool drawfit)
{

	switch(butterfly) { // butterfly
	case NoButterfly : 
		fButterfly = NoButterfly;
		break;
	case LinearButterfly : 
		fButterfly = LinearButterfly;
		break;
	case LogarithmButterfly : 
		fButterfly = LogarithmButterfly;
		break;
	case ContoursButterfly : 
		fButterfly = ContoursButterfly;
		break;
	case CausticButterfly : 
		fButterfly = CausticButterfly;
		break;
	default :
		fButterfly = NoButterfly;
	}

	fdrawfit = drawfit; // fit 

	/*
	  CleanPolyLines();
	  CleanGraphArraySpectrum();
	  CleanTF1FitSpectrum();
	  CleanGraphSpectrum();
	  CleanTH2();
	  CleanTF1Residuals();
	  CleanPads();
	  CleanCanvasSpectrum();
	*/

	for(std::vector<Hypothesis*>::iterator hypo=fHypothesisArray.begin(); hypo!=fHypothesisArray.end(); hypo++) {

		if((*hypo)->GetConvergence()) {

			DEBUG_OUT << "Plotting " <<  (*hypo)->GetName() << std::endl;

			PlotHypothesis(**hypo);

		}
		else {
			WARNING << "Skip hypothesis " << (*hypo)->GetName() << " because it didn't converge!" << std::endl;
			continue;
		}

	}

	switch(butterfly) { // butterfly
	case NoButterfly : 
		INFO << "Plotting NoButterfly Spectrums... ok" << std::endl;
		break;
	case LinearButterfly : 
		INFO << "Plotting LinearButterfly Spectrums... ok" << std::endl;
		break;
	case LogarithmButterfly : 
		INFO << "Plotting LogarithmButterfly Spectrums... ok" << std::endl;
		break;
	case ContoursButterfly : 
		INFO << "Plotting ContoursButterfly Spectrums... ok" << std::endl;
		break;
	case CausticButterfly : 
		INFO << "Plotting CausticButterfly Spectrums... ok" << std::endl;
		break;
	default :
		INFO << "Plotting Spectrum... ok" << std::endl;
	}

}

/**
 * \brief Plot light curves
 */
void START::PlotFactory::PlotLightCurves(LightCurvesTimeAxis TimeAxisOption,LightCurvesErrors ErrorsOption) {

	fLightCurvesTimeAxis=TimeAxisOption;
	fLightCurvesErrors=ErrorsOption;

	//CleanCanvasLightCurves();
	//CleanGraphLightCurves();

	std::string errorshandling;
	switch(fLightCurvesErrors) {
	case Rolke:
		errorshandling="Rolke";
		break;
	case Gaussian:
		errorshandling="Gaussian";
	}

	for(std::vector<Hypothesis*>::iterator hypo=fHypothesisArray.begin(); hypo!=fHypothesisArray.end(); hypo++) {

		if((*hypo)->GetConvergence()) {

			DEBUG_OUT << (*hypo)->GetName() << " gets " << (*hypo)->GetMapTimeBinVector().size() << " lightcurve(s)" << std::endl;

			std::map<std::string,TimeBinVector> MapLC = (*hypo)->GetMapTimeBinVector();

			for(std::map<std::string,TimeBinVector>::const_iterator MapTBin=MapLC.begin(); 
				MapTBin!=MapLC.end(); ++MapTBin) {
	
				std::string lcname(MapTBin->first);

				if(MapTBin->second.tbin.size()==0) {
					WARNING << "Time bin vector is empty! ==> Skip this one: " << lcname << "!" << std::endl;
					continue;
				}

				InitCanvasLightCurves(**hypo,lcname);
				InitTH2LightCurves(**hypo,lcname);
				InitPadLightCurves(**hypo,lcname);
				InitGraphLightCurves(**hypo,lcname);
				InitTF1MeanFluxLightCurves(**hypo,lcname);
				InitTF1ZeroFluxLightCurves(**hypo,lcname);
				DrawLightCurves(**hypo,lcname);

				INFO << "LightCurves " << lcname << " with " << errorshandling 
					 << " errors for hypothesis " << (*hypo)->GetName() << "... done" << std::endl;  

			}
      
		}
		else {
			INFO << "Skip hypothesis " << (*hypo)->GetName() << " because it didn't converge!" << std::endl;
			continue;
		}
	}

	INFO << "Plotting LightCurves... ok" << std::endl;

}

/**
 * \brief Init Canvas for light curves
 */
void START::PlotFactory::InitCanvasLightCurves(const Hypothesis &hypo,std::string lctype) {

	DEBUG_OUT << "START!" << std::endl;

	int canvaswidth(0), canvasheight(0); 
  
	switch(fPlotStyle) {
	case Default:
	case Paper:
	case UserFriendly:
		canvaswidth=800;
		canvasheight=360;
	}

	std::string canvasname = MakeKeyCanvasLightCurves(hypo.GetName(),lctype);
	fMapCanvasLightCurves[canvasname] = new TCanvas(canvasname.c_str(),canvasname.c_str(),canvaswidth,canvasheight);

	DEBUG_OUT << "END!" << std::endl;

}

/**
 * \brief Initialize LC pad
 */
void START::PlotFactory::InitPadLightCurves(Hypothesis &hypo,std::string lctype)
{

	DEBUG_OUT << "START!" << std::endl;

	// Definition of spectrum's and residuals' pad size

	switch(fPlotStyle) {
	case Default:
	case UserFriendly:
	case Paper:
		fx1padlc=0.01;
		fy1padlc=0.01;
		fx2padlc=0.99;
		fy2padlc=0.99;
	}

	std::string padlcname = MakeKeyPadLightCurves(hypo.GetName(),lctype);
	fMapPadLightCurves[padlcname] = new TPad(padlcname.c_str(),padlcname.c_str(),fx1padlc,fy1padlc,fx2padlc,fy2padlc,0);
	fMapPadLightCurves[padlcname]->SetGridx();
	fMapPadLightCurves[padlcname]->SetGridy();
	fMapPadLightCurves[padlcname]->SetFillColor(0);

	switch(fPlotStyle) {
	case Default:
	case Paper:
	case UserFriendly:
		fMapPadLightCurves[padlcname]->SetTopMargin(0.1);
		fMapPadLightCurves[padlcname]->SetBottomMargin(0.15);
		fMapPadLightCurves[padlcname]->SetRightMargin(0.04);
	}

	DEBUG_OUT << "END!" << std::endl;

}

/**
 * \brief Initialize LC pad
 */
void START::PlotFactory::InitTH2LightCurves(Hypothesis &hypo,std::string lctype)
{

	DEBUG_OUT << "START!" << std::endl;

	std::string th2lcname = MakeKeyTH2LightCurves(hypo.GetName(),lctype);
	std::pair<double,double> timerange = GetLightCurvesPlotTimeRange(hypo,lctype);
	std::pair<double,double> fluxrange = GetLightCurvesPlotFluxRange(hypo,lctype);

	fMapTH2LightCurves[th2lcname] = new TH2D(th2lcname.c_str(),th2lcname.c_str(),
											 100,timerange.first,timerange.second,
											 100,fluxrange.first,fluxrange.second);
	fMapTH2LightCurves[th2lcname]->SetStats(kFALSE);
	fMapTH2LightCurves[th2lcname]->SetTitle("");

	switch(fPlotStyle) {
	case Default:
	case Paper:
	case UserFriendly:
		std::ostringstream xaxisname;
		xaxisname << "Time "; 
		xaxisname << "(MJD)";
		fMapTH2LightCurves[th2lcname]->GetXaxis()->SetTitle(xaxisname.str().c_str());
		fMapTH2LightCurves[th2lcname]->GetXaxis()->SetTitleSize(0.06);
		fMapTH2LightCurves[th2lcname]->GetXaxis()->SetTitleOffset(1.10);
		fMapTH2LightCurves[th2lcname]->GetXaxis()->SetLabelSize(0.05);
		fMapTH2LightCurves[th2lcname]->GetXaxis()->SetLabelOffset(0.025);
    

		std::ostringstream yaxisname;
		double energymin = hypo.GetLightCurveIntegratedFluxEnergyRange().first;
		//yaxisname << "#Phi_{I}(E>" << yaxisname.precision(0) << energymin*1000.;
		yaxisname << "#Phi_{I}(E>" << energymin*1000.;
		yaxisname << " GeV)";
		yaxisname << " " << hypo.GetIntegratedFluxUnits();
      
		fMapTH2LightCurves[th2lcname]->GetYaxis()->SetTitle(yaxisname.str().c_str());
		fMapTH2LightCurves[th2lcname]->GetYaxis()->SetTitleSize(0.06);
		fMapTH2LightCurves[th2lcname]->GetYaxis()->SetTitleOffset(0.54);
		fMapTH2LightCurves[th2lcname]->GetYaxis()->SetLabelSize(0.05);
	}

	// Set time axis if needed
	TDatime *TimeHandler = 0;
	switch(fLightCurvesTimeAxis) {
	case MJD:
		break;
	case YearMonth:
		fMapTH2LightCurves[th2lcname]->GetXaxis()->SetTimeDisplay(1);
		fMapTH2LightCurves[th2lcname]->GetXaxis()->SetTimeFormat("#splitline{%b}{%Y}");
		TimeHandler = new TDatime(2000,1,1,0,0,0); // 1st january 2000
		fMapTH2LightCurves[th2lcname]->GetXaxis()->SetTimeOffset(TimeHandler->Convert());
		fMapTH2LightCurves[th2lcname]->GetXaxis()->SetTitle("");
		fMapTH2LightCurves[th2lcname]->GetXaxis()->SetLabelOffset(0.057);
		break;
	case YearMonthDay:
		fMapTH2LightCurves[th2lcname]->GetXaxis()->SetTimeDisplay(1);
		fMapTH2LightCurves[th2lcname]->GetXaxis()->SetTimeFormat("%d/%m/%y");
		TimeHandler = new TDatime(2000,1,1,0,0,0); // 1st january 2000
		fMapTH2LightCurves[th2lcname]->GetXaxis()->SetTimeOffset(TimeHandler->Convert());
		fMapTH2LightCurves[th2lcname]->GetXaxis()->SetTitle("");
		fMapTH2LightCurves[th2lcname]->GetXaxis()->SetLabelOffset(0.057);
		break;
	case HourMinuteSecond:
		fMapTH2LightCurves[th2lcname]->GetXaxis()->SetTimeDisplay(1);
		fMapTH2LightCurves[th2lcname]->GetXaxis()->SetTimeFormat("%H:%M:%S");
		TimeHandler = new TDatime(2000,1,1,0,0,0); // 1st january 2000
		fMapTH2LightCurves[th2lcname]->GetXaxis()->SetTimeOffset(TimeHandler->Convert());
		fMapTH2LightCurves[th2lcname]->GetXaxis()->SetTitle("");
		fMapTH2LightCurves[th2lcname]->GetXaxis()->SetLabelOffset(0.057);
	}

	DEBUG_OUT << "END!" << std::endl;

}

/**
 * \brief Set time range for TH2's
 */
std::pair<double,double> START::PlotFactory::GetLightCurvesPlotTimeRange(const Hypothesis &hypothesis,std::string lctype)
{

	DEBUG_OUT << "START!" << std::endl;

	std::vector<TimeBin> TBin = hypothesis.GetMapTimeBinVector()[lctype].tbin;
	std::vector<double> vtmin;
	std::vector<double> vtmax;

	for(std::vector<TimeBin>::const_iterator bin=TBin.begin(); bin!=TBin.end(); ++bin) {
		vtmin.push_back(bin->GetTimeMin());
		vtmax.push_back(bin->GetTimeMax());
	}
	double tmin(0.), tmax(0.);
	std::vector<double>::const_iterator itmin, itmax;
	if(vtmax.size()>0) {
		itmin = min_element(vtmin.begin(),vtmin.end());
		itmax = max_element(vtmax.begin(),vtmax.end());
		tmin=*itmin;
		tmax=*itmax;
	}
	else {
		WARN_OUT << "Null size!! ==> Bug expected..." << std::endl;
	}

	double mjdreference(STARTUtils::GetLightCurveMJDReferenceTime()); // 1st january 2000
	double dayinseconds(60.*60.*24.);

	switch(fLightCurvesTimeAxis) {
	case MJD:
		break;
	case YearMonth:
	case YearMonthDay:
	case HourMinuteSecond:
		tmin=tmin-mjdreference;
		tmin*=dayinseconds;
		tmax=tmax-mjdreference;
		tmax*=dayinseconds;
	}

	double shift(0.1*(tmax-tmin));

	DEBUG_OUT << "END!" << std::endl;

	return std::make_pair(tmin-shift,tmax+shift);
}

/**
 * \brief Set flux range for TH2's
 */ 
std::pair<double,double> START::PlotFactory::GetLightCurvesPlotFluxRange(const Hypothesis &hypothesis,std::string lctype)
{

	DEBUG_OUT << "START!" << std::endl;

	std::vector<TimeBin> TimeBinnedData = hypothesis.GetMapTimeBinVector()[lctype].tbin;
	std::vector<double> vfluxmin;
	std::vector<double> vfluxmax;
  
	std::vector<double> flux, gausserr, rolke1min, rolke1max, rolke3min, rolke3max;
	std::vector<double> isupperlimit;

	for(std::vector<TimeBin>::const_iterator bin=TimeBinnedData.begin(); bin!=TimeBinnedData.end(); ++bin) {
		flux.push_back(bin->GetIntegratedFlux());
		gausserr.push_back(bin->GetIntegratedFluxError());

		rolke1min.push_back(bin->GetIntegratedFluxError1SigmaMinus());
		rolke1max.push_back(bin->GetIntegratedFluxError1SigmaPlus());
		rolke3min.push_back(bin->GetIntegratedFluxError3SigmaMinus());
		rolke3max.push_back(bin->GetIntegratedFluxError3SigmaPlus());
    
		isupperlimit.push_back(bin->GetIsUpperLimit());
	}

	for(unsigned int ipoint(0); ipoint<flux.size(); ipoint++) {
		switch(fLightCurvesErrors) {
		case Rolke:
			if(isupperlimit[ipoint]==true) { // upperlimit
				vfluxmax.push_back(rolke3max[ipoint]);
				vfluxmin.push_back(rolke3max[ipoint]*fscalearrowfactor);
			}
			else {
				vfluxmax.push_back(rolke1max[ipoint]);
				vfluxmin.push_back(rolke1min[ipoint]);
			}
			break;
		case Gaussian:
			vfluxmax.push_back(flux[ipoint]+gausserr[ipoint]);
			vfluxmin.push_back(flux[ipoint]-gausserr[ipoint]);
		}
	}

	double fluxmin(0.), fluxmax(0.);
	std::vector<double>::const_iterator itmin, itmax;
	if(vfluxmax.size()>0) {
		itmin = min_element(vfluxmin.begin(),vfluxmin.end());
		itmax = max_element(vfluxmax.begin(),vfluxmax.end());
		fluxmin=*itmin;
		fluxmax=*itmax;
	}
	else {
		WARN_OUT << "Null size!! ==> Bug expected..." << std::endl;
	}

	double shift(0.1*(fluxmax-fluxmin));

	//if(fluxmin-shift>0.) fluxmin=shift; // in order to have zero on the plot

	DEBUG_OUT << "END!" << std::endl;

	return std::make_pair(fluxmin-shift,fluxmax+shift);

}

/**
 * \brief Init graphs for light curves
 */
void START::PlotFactory::InitGraphLightCurves(const Hypothesis &hypo,std::string lctype) {

	DEBUG_OUT << "START!" << std::endl;
    
	if(DEBUG) hypo.GetMapTimeBinVector()[lctype].Print();
	std::vector<double> time, timemin, timemax;
	std::vector<double> flux, gausserr, rolke1min, rolke1max, rolke3min, rolke3max;
	std::vector<double> isupperlimit;

	TimeBinVector TimeBinnedData = hypo.GetMapTimeBinVector()[lctype];

	double mjdreference(STARTUtils::GetLightCurveMJDReferenceTime());
	double dayinseconds(60.*60.*24.);
	for(std::vector<TimeBin>::const_iterator bin=TimeBinnedData.tbin.begin(); bin!=TimeBinnedData.tbin.end(); ++bin) {
		time.push_back(bin->GetTimeMean());
		timemin.push_back(bin->GetTimeMin());
		timemax.push_back(bin->GetTimeMax());
    
		flux.push_back(bin->GetIntegratedFlux());
		gausserr.push_back(bin->GetIntegratedFluxError());

		rolke1min.push_back(bin->GetIntegratedFluxError1SigmaMinus());
		rolke1max.push_back(bin->GetIntegratedFluxError1SigmaPlus());
		rolke3min.push_back(bin->GetIntegratedFluxError3SigmaMinus());
		rolke3max.push_back(bin->GetIntegratedFluxError3SigmaPlus());

		isupperlimit.push_back(bin->GetIsUpperLimit());

		switch(fLightCurvesTimeAxis) {
		case MJD:
			break;
		case YearMonth: // convert in second because of time axis
		case YearMonthDay:
		case HourMinuteSecond:
			time.back()-=mjdreference;
			time.back()*=dayinseconds;
			timemin.back()-=mjdreference;
			timemin.back()*=dayinseconds;
			timemax.back()-=mjdreference;
			timemax.back()*=dayinseconds;
		}

	}
	DEBUG_OUT << "After filling vectors" << std::endl;
	std::string graphname = MakeKeyGraphLightCurves(hypo.GetName(),lctype);
	fMapGraphLightCurves[graphname] = new TGraphAsymmErrors(time.size());
	fMapGraphLightCurves[graphname]->SetName(graphname.c_str());
	std::string graphtitle = "Light Curve";
	fMapGraphLightCurves[graphname]->SetTitle(graphtitle.c_str());
	fMapGraphLightCurves[graphname]->SetMarkerColor(kRed);
	fMapGraphLightCurves[graphname]->SetMarkerStyle(20);
	fMapGraphLightCurves[graphname]->SetMarkerSize(0.6);
  
	std::string arrowlcname = MakeKeyArrowsLightCurves(hypo.GetName(),lctype);

	DEBUG_OUT << "After TGraphErrors initialization" << std::endl;

	unsigned int graphpoints(0);
	switch(fLightCurvesErrors) {
	case Rolke:
		for(unsigned int ipoint(0); ipoint<time.size(); ipoint++) {
      
			if(isupperlimit[ipoint]==false) {
				fMapGraphLightCurves[graphname]->SetPoint(graphpoints,time[ipoint],flux[ipoint]);
				fMapGraphLightCurves[graphname]->SetPointError(graphpoints,time[ipoint]-timemin[ipoint],
															   timemax[ipoint]-time[ipoint],
															   flux[ipoint]-rolke1min[ipoint],
															   rolke1max[ipoint]-flux[ipoint]);
				graphpoints++;
				DEBUG_OUT << "point " << ipoint << " x=" << time[ipoint] << " y=" << flux[ipoint] << std::endl;
			}
			else { // upper limit

				fMapArrowsLightCurves[arrowlcname].push_back(new TArrow(time[ipoint],
																		rolke3max[ipoint],
																		time[ipoint],
																		//(rolke3max[ipoint])*0.8,//fscalearrowfactor,
																		0.,
																		0.01,
																		"|-|>"));

				DEBUG_OUT_L(2) << "LC : adding upperlimit en= " << time[ipoint] << " upperlimit = " << rolke3max[ipoint] << " bottom arrow = " << 
					(rolke3max[ipoint])*fscalearrowfactor << std::endl;
				fMapArrowsLightCurves[arrowlcname].back()->SetLineWidth(1);
				fMapArrowsLightCurves[arrowlcname].back()->SetLineColor(kBlack);
				fMapArrowsLightCurves[arrowlcname].back()->SetFillColor(kBlack);
			}

		}
		break;
	case Gaussian:
		for(unsigned int ipoint(0); ipoint<time.size(); ipoint++) {
			fMapGraphLightCurves[graphname]->SetPoint(graphpoints,time[ipoint],flux[ipoint]);
			fMapGraphLightCurves[graphname]->SetPointError(graphpoints,time[ipoint]-timemin[ipoint],
														   timemax[ipoint]-time[ipoint],
														   gausserr[ipoint],
														   gausserr[ipoint]);
			graphpoints++;
		}
    
	}
  
	DEBUG_OUT << "END!" << std::endl;

}


/**
 * \brief Draw light curves
 */
void START::PlotFactory::DrawLightCurves(const Hypothesis &hypo,std::string lctype) {

	DEBUG_OUT << "START!" << std::endl;

	std::string canvasname = MakeKeyCanvasLightCurves(hypo.GetName(),lctype);
	DEBUG_OUT << "canvas=" << fMapCanvasLightCurves[canvasname] << std::endl;
	fMapCanvasLightCurves[canvasname]->Draw();
	fMapCanvasLightCurves[canvasname]->cd();

	std::string padname = MakeKeyPadLightCurves(hypo.GetName(),lctype);
	DEBUG_OUT << "pad=" << fMapPadLightCurves[padname] << std::endl;
	fMapPadLightCurves[padname]->Draw();
	fMapPadLightCurves[padname]->cd();

	std::string th2name = MakeKeyTH2LightCurves(hypo.GetName(),lctype);
	DEBUG_OUT << "th2=" << fMapTH2LightCurves[th2name] << std::endl;
	fMapTH2LightCurves[th2name]->Draw("");

	std::string tf1meanfluxname = MakeKeyTF1MeanFluxLightCurves(hypo.GetName(),lctype);
	DEBUG_OUT << "tf1mean=" << fMapTF1MeanFluxLightCurves[tf1meanfluxname] << std::endl;
	gStyle->SetOptFit(1111);
	std::string graphname = MakeKeyGraphLightCurves(hypo.GetName(),lctype);
	DEBUG_OUT << "graph=" << fMapGraphLightCurves[graphname] << std::endl;
	fMapGraphLightCurves[graphname]->Draw("P same");
	fMapGraphLightCurves[graphname]->Fit(fMapTF1MeanFluxLightCurves[tf1meanfluxname],"Q0");

	fMapTF1MeanFluxLightCurves[tf1meanfluxname]->Draw("same");

	std::string tf1zerofluxname = MakeKeyTF1ZeroFluxLightCurves(hypo.GetName(),lctype);
	DEBUG_OUT << "tf1zero=" << fMapTF1ZeroFluxLightCurves[tf1zerofluxname] << std::endl;
	fMapTF1ZeroFluxLightCurves[tf1zerofluxname]->Draw("same");

	// draw upper limits
	std::string arrowslcname = MakeKeyArrowsLightCurves(hypo.GetName(),lctype);
	DEBUG_OUT << "arrows = " << arrowslcname << std::endl;
	for(unsigned int i(0); i<fMapArrowsLightCurves[arrowslcname].size(); i++) {
		if(fMapArrowsLightCurves[arrowslcname][i]!=0) fMapArrowsLightCurves[arrowslcname][i]->Draw(""); // draw spectrum's upperlimit
	}

	// Draw source name
	fMapPadLightCurves[padname]->cd();
	std::ostringstream plotname;
	plotname << fsourcename;
	TLatex *drawsourcename = new TLatex();
	drawsourcename->SetNDC();
	drawsourcename->SetTextSize(0.045);
	drawsourcename->SetTextColor(1);
	drawsourcename->SetTextFont(52);
	switch(fPlotStyle) {
	case Paper:
		break;
	case Default:
	case UserFriendly:
		drawsourcename->DrawLatex(0.16,0.9,plotname.str().c_str());
	}
	DEBUG_OUT << "END!" << std::endl;

}

/**
 * \brief Init LC TF1 Mean flux
 */
void START::PlotFactory::InitTF1MeanFluxLightCurves(Hypothesis &hypo,std::string lctype) {

	DEBUG_OUT << "START!" << std::endl;

	double integratedflux=hypo.GetFluxIntegralFitParams(hypo.GetLightCurveIntegratedFluxEnergyRange().first,
														hypo.GetLightCurveIntegratedFluxEnergyRange().second);
	std::pair<double,double> timerange = GetLightCurvesPlotTimeRange(hypo,lctype);
	std::string tf1name = MakeKeyTF1MeanFluxLightCurves(hypo.GetName(),lctype);
	fMapTF1MeanFluxLightCurves[tf1name] = new TF1(tf1name.c_str(),"[0]",timerange.first,timerange.second);
	fMapTF1MeanFluxLightCurves[tf1name]->SetParameter(0,integratedflux);
	fMapTF1MeanFluxLightCurves[tf1name]->SetLineColor(kGreen-3);
	fMapTF1MeanFluxLightCurves[tf1name]->SetLineStyle(2);
	fMapTF1MeanFluxLightCurves[tf1name]->SetLineWidth(1);

	DEBUG_OUT << "END!" << std::endl;
}

/**
 * \brief Init LC TF1 Zero flux
 */
void START::PlotFactory::InitTF1ZeroFluxLightCurves(Hypothesis &hypo,std::string lctype) {

	DEBUG_OUT << "START!" << std::endl;

	std::string tf1name = MakeKeyTF1ZeroFluxLightCurves(hypo.GetName(),lctype);
	std::pair<double,double> timerange = GetLightCurvesPlotTimeRange(hypo,lctype);
	fMapTF1ZeroFluxLightCurves[tf1name] = new TF1(tf1name.c_str(),"[0]",timerange.first,timerange.second);
	fMapTF1ZeroFluxLightCurves[tf1name]->SetParameter(0,0.);
	fMapTF1ZeroFluxLightCurves[tf1name]->SetLineColor(kBlack);
	fMapTF1ZeroFluxLightCurves[tf1name]->SetLineStyle(1);
	fMapTF1ZeroFluxLightCurves[tf1name]->SetLineWidth(1);

	DEBUG_OUT << "END!" << std::endl;

}

/**
 * \brief Clean LC canvas
 */
void START::PlotFactory::CleanCanvasLightCurves() {
	for(std::map<std::string,TCanvas*>::iterator canvas=fMapCanvasLightCurves.begin(); canvas!=fMapCanvasLightCurves.end(); canvas++) {
		if(canvas->second!=0) delete canvas->second;
		canvas->second=0;
	}
	fMapCanvasLightCurves.clear();
}

/**
 * \brief Clean memory
 */
void START::PlotFactory::CleanTH2LightCurves()
{

	for(std::map<std::string,TH2D*>::iterator th2=fMapTH2LightCurves.begin(); th2!=fMapTH2LightCurves.end(); th2++) {
		if(th2->second!=0) delete th2->second;
		th2->second=0;
	}
	fMapTH2LightCurves.clear();

}

/**
 * \brief Clean memory
 */
void START::PlotFactory::CleanPadLightCurves()
{
  
	for(std::map<std::string,TPad*>::iterator pad=fMapPadLightCurves.begin(); pad!=fMapPadLightCurves.end(); pad++) {
		if(pad->second!=0) delete pad->second;
		pad->second=0;
	}
	fMapPadLightCurves.clear();

}

/**
 * \brief Clean LC graph
 */
void START::PlotFactory::CleanGraphLightCurves() {



	for(std::map<std::string,TGraphAsymmErrors*>::iterator graph=fMapGraphLightCurves.begin(); graph!=fMapGraphLightCurves.end(); graph++) {
		if(graph->second!=0) delete graph->second;
		graph->second=0;
	}
	fMapGraphLightCurves.clear();

	for(std::map<std::string,std::vector<TArrow*> >::iterator arrow=fMapArrowsLightCurves.begin(); arrow!=fMapArrowsLightCurves.end(); arrow++) {
		for(unsigned int i(0); i<fMapArrowsLightCurves[arrow->first].size(); i++) {
			if(arrow->second[i]!=0) delete arrow->second[i];
			arrow->second[i]=0;
		}
		fMapArrowsLightCurves[arrow->first].clear();
	}

}

/**
 * \brief Clean LC TF1 mean flux
 */
void START::PlotFactory::CleanTF1MeanFluxLightCurves() {
	for(std::map<std::string,TF1*>::iterator tf1=fMapTF1MeanFluxLightCurves.begin(); tf1!=fMapTF1MeanFluxLightCurves.end(); tf1++) {
		if(tf1->second!=0) delete tf1->second;
		tf1->second=0;
	}
	fMapTF1MeanFluxLightCurves.clear();
}

/**
 * \brief Clean LC TF1 zero flux
 */
void START::PlotFactory::CleanTF1ZeroFluxLightCurves() {
	for(std::map<std::string,TF1*>::iterator tf1=fMapTF1ZeroFluxLightCurves.begin(); tf1!=fMapTF1ZeroFluxLightCurves.end(); tf1++) {
		if(tf1->second!=0) delete tf1->second;
		tf1->second=0;
	}
	fMapTF1ZeroFluxLightCurves.clear();
}

/**
 * \brief Save Light curve objetcs in root file
 */
void START::PlotFactory::SaveLightCurvesROOTObjects() {

	fSaveRootFile->cd();
  
	for(std::map<std::string,TCanvas*>::iterator canvas=fMapCanvasLightCurves.begin(); canvas!=fMapCanvasLightCurves.end(); canvas++) {
		DEBUG_OUT << "canvas " << canvas->first << " exists" << std::endl;
		if(canvas->second!=0) canvas->second->Write();
	}

}


/**
 * \brief Init TPaveText
 */
void START::PlotFactory::InitPaveTextSpectrum(Hypothesis &hypo)
{

	double x1(0.),y1(0.),x2(0.),y2(0.);

	switch(fPlotStyle) {
	case Default:
	case UserFriendly:
	case Paper:
		x1 = 0.6;
		y1 = 0.60;
		x2 = 0.93;
		y2 = 0.93;
	}

	std::string pavetextname = MakeKeyPaveTextSpectrum(hypo.GetName());
	fMapPaveTextSpectrum[pavetextname] = new TPaveText(x1,y1,x2,y2,"NDC");

	fMapPaveTextSpectrum[pavetextname]->SetShadowColor(0);
	fMapPaveTextSpectrum[pavetextname]->SetFillColor(0);
	fMapPaveTextSpectrum[pavetextname]->SetFillStyle(1001);
	fMapPaveTextSpectrum[pavetextname]->SetLineColor(0);
	fMapPaveTextSpectrum[pavetextname]->SetBorderSize(0);
	fMapPaveTextSpectrum[pavetextname]->SetTextAlign(11);
	fMapPaveTextSpectrum[pavetextname]->SetTextFont(62);
	fMapPaveTextSpectrum[pavetextname]->SetTextSize(0.028);
	fMapPaveTextSpectrum[pavetextname]->SetMargin(0.);

	//fMapPaveTextSpectrum[pavetextname]->SetLabel(fsourcename.c_str());

	// print first lines

	//std::ostringstream introhypo1;
	//introhypo1 << fsourcename;
	//fMapPaveTextSpectrum[pavetextname]->AddText(fsourcename.str().c_str());

	std::ostringstream introhypo2;
	introhypo2 << "Minimization summary for " << hypo.GetModelName() << " : "; 
	fMapPaveTextSpectrum[pavetextname]->AddText(introhypo2.str().c_str());

	// print fitted parameters
	for(unsigned int ipar(0); ipar<hypo.GetParametersNb(); ipar++) {

		std::vector<std::pair<unsigned int,double> > normalizedparam(hypo.GetVectorNormalizedParameters());
		bool isnormalized(false);

		for(std::vector<std::pair<unsigned int,double> >::const_iterator normpar=normalizedparam.begin(); normpar!=normalizedparam.end(); ++normpar) {
			if(normpar->first==ipar) {
				isnormalized=true;
				break;
			}
		}

		int exponent(STARTUtils::GetIntegerExposantFromDecimalNumber(hypo.GetFittedParameters()[ipar]));

		std::ostringstream fittedparam;
		fittedparam << hypo.GetLatexParametersNames()[ipar] << " = (";
		if(isnormalized) {
			fittedparam.precision(3);
			fittedparam.setf(std::ios::fixed);
			fittedparam.width(4);
			fittedparam << hypo.GetFittedParameters()[ipar]/TMath::Power(10.,exponent) << " +/- ";
			fittedparam.precision(3);
			fittedparam.setf(std::ios::fixed);
			fittedparam.width(4);
			fittedparam << hypo.GetFittedParametersErrors()[ipar]/TMath::Power(10.,exponent) << ") 10^{";
			fittedparam.precision(0);
			fittedparam << exponent << "} ";
		}
		else {
			fittedparam.precision(3);
			fittedparam.setf(std::ios::fixed);
			fittedparam.width(4);
			fittedparam << hypo.GetFittedParameters()[ipar] << " +/- ";
			fittedparam.precision(3);
			fittedparam.setf(std::ios::fixed);
			fittedparam.width(4);
			fittedparam << hypo.GetFittedParametersErrors()[ipar] << ") ";
		}
		fittedparam << hypo.GetParametersUnits()[ipar];
    
		if( std::string(hypo.GetName()).find("EXP")!=std::string::npos
			&& hypo.GetParametersNames()[ipar].find("beta")!=std::string::npos ) {
			float Ec = 1./hypo.GetFittedParameters()[ipar];
			float eEc = pow(Ec,2) * hypo.GetFittedParametersErrors()[ipar];
			fittedparam << "  Ec = (";
			fittedparam.precision(3);
			fittedparam.setf(std::ios::fixed);
			fittedparam.width(4);
			fittedparam << Ec << " +/-";
			fittedparam.precision(3);
			fittedparam.setf(std::ios::fixed);
			fittedparam.width(4);
			fittedparam << eEc << ") TeV";
		}
    
		fMapPaveTextSpectrum[pavetextname]->AddText(fittedparam.str().c_str());
	}

	/*
	// print fitted parameters
	for(unsigned int ipar(0); ipar<hypo.GetParametersNb(); ipar++) {
    std::ostringstream fittedparam;
    fittedparam << hypo.GetLatexParametersNames()[ipar] << " = (";
    fittedparam.precision(3);
    fittedparam << hypo.GetFittedParameters()[ipar] << " +/- ";
    fittedparam.precision(3);
    fittedparam << hypo.GetFittedParametersErrors()[ipar];
    fittedparam << " ) " << hypo.GetParametersUnits()[ipar];
    fMapPaveTextSpectrum[pavetextname]->AddText(fittedparam.str().c_str());
	}
	*/

	// print flux at 1 TeV
	std::ostringstream fluxref;
	double eref = 1.;
	int exponent(STARTUtils::GetIntegerExposantFromDecimalNumber(hypo.GetFluxFitParams(eref)));
	fluxref << "#Phi(1 TeV) = (";
	fluxref.precision(3);
	fluxref.setf(std::ios::fixed);
	fluxref.width(4);
	fluxref << hypo.GetFluxFitParams(eref)/TMath::Power(10.,exponent) << " +/- ";
	fluxref.precision(3);
	fluxref.setf(std::ios::fixed);
	fluxref.width(4);
	fluxref << hypo.GetSigmaFlux(eref)/TMath::Power(10.,exponent) << ") 10^{";
	fluxref.precision(0);
	fluxref<< exponent << "} " << hypo.GetFluxUnits();
	fMapPaveTextSpectrum[pavetextname]->AddText(fluxref.str().c_str());

	if(hypo.GetSpectralType()==Hypothesis::Differential) {

		// decorrelation energy
		double edecorrelation=-1;
		switch(fButterfly) {
		case NoButterfly:
		case LinearButterfly:
			edecorrelation = hypo.GetLinearDecorrelationEnergy();
			break;
		case LogarithmButterfly:
			edecorrelation = hypo.GetLogarithmDecorrelationEnergy();
			break;
		case ContoursButterfly:
		case CausticButterfly:
			edecorrelation = hypo.GetContoursDecorrelationEnergy();
		}

		// decorrelation energy
		std::ostringstream decorrenergy;
		decorrenergy << "E_{Dec} = ";
		decorrenergy.precision(3);
		decorrenergy.setf(std::ios::fixed);
		decorrenergy.width(4);
		decorrenergy << edecorrelation << " TeV";
		fMapPaveTextSpectrum[pavetextname]->AddText(decorrenergy.str().c_str());
    
		// flux at decorrelation energy
		std::ostringstream fluxdec;
		exponent = STARTUtils::GetIntegerExposantFromDecimalNumber(hypo.GetFluxFitParams(edecorrelation));
		fluxdec << "#Phi(";
		fluxdec.precision(3);
		fluxdec.setf(std::ios::fixed);
		fluxdec.width(4);
		fluxdec << edecorrelation << " TeV) = (";
		fluxdec.precision(3);
		fluxdec.setf(std::ios::fixed);
		fluxdec.width(4);
		fluxdec << hypo.GetFluxFitParams(edecorrelation)/TMath::Power(10.,exponent) << " +/- ";
		fluxdec.precision(3);
		fluxdec.setf(std::ios::fixed);
		fluxdec.width(4);
		fluxdec << hypo.GetSigmaFlux(edecorrelation)/TMath::Power(10.,exponent) << ") 10^{";
		fluxdec.precision(0);
		fluxdec << exponent << "} " << hypo.GetFluxUnits();
		fMapPaveTextSpectrum[pavetextname]->AddText(fluxdec.str().c_str());
    
		// print integrated flux and energy flux
		if(hypo.GetFittedIntegratedFlux().first!=-1 && hypo.GetFittedIntegratedFlux().second!=-1) {
			std::ostringstream integratedflux;
			exponent = STARTUtils::GetIntegerExposantFromDecimalNumber(hypo.GetFittedIntegratedFlux().first);
			integratedflux << "I_{F}(E>";
			integratedflux.precision(3);
			integratedflux.setf(std::ios::fixed);
			integratedflux.width(4);
			integratedflux << hypo.GetIntegratedFitEnergyRange().first << " TeV) = (";
			integratedflux.precision(3);
			integratedflux.setf(std::ios::fixed);
			integratedflux.width(4);
			integratedflux << hypo.GetFittedIntegratedFlux().first/TMath::Power(10.,exponent) << " +/- ";
			integratedflux.precision(3);
			integratedflux.setf(std::ios::fixed);
			integratedflux.width(4);
			integratedflux << hypo.GetFittedIntegratedFlux().second/TMath::Power(10.,exponent) << ") 10^{";
			integratedflux.precision(0);
			integratedflux << exponent << "} " << hypo.GetFluxUnits();
			fMapPaveTextSpectrum[pavetextname]->AddText(integratedflux.str().c_str());
		}
    
		if(hypo.GetFittedEnergyFlux().first!=-1 && hypo.GetFittedEnergyFlux().second!=-1) {
			std::ostringstream energyflux;
			exponent = STARTUtils::GetIntegerExposantFromDecimalNumber(hypo.GetFittedEnergyFlux().first);
			energyflux << "E_{F}(E>";
			energyflux.precision(3);
			energyflux.setf(std::ios::fixed);
			energyflux.width(4);
			energyflux << hypo.GetIntegratedFitEnergyRange().first << " TeV) = (";
			energyflux.precision(3);
			energyflux.setf(std::ios::fixed);
			energyflux.width(4);
			energyflux << hypo.GetFittedEnergyFlux().first/TMath::Power(10.,exponent) << " +/- ";
			energyflux.precision(3);
			energyflux.setf(std::ios::fixed);
			energyflux.width(4);
			energyflux << hypo.GetFittedEnergyFlux().second/TMath::Power(10.,exponent) << ") 10^{";
			energyflux.precision(0);
			energyflux << exponent << "} " << hypo.GetEnergyFluxUnits();
			fMapPaveTextSpectrum[pavetextname]->AddText(energyflux.str().c_str());
		}
    
	}

	if(!fminimalinfosonpad) {
		std::ostringstream mnz_max;
		mnz_max << "Log L_{MAX} = ";
		mnz_max.precision(3);
		mnz_max << hypo.GetMaximumLikelihood();
		fMapPaveTextSpectrum[pavetextname]->AddText(mnz_max.str().c_str());
    
		std::ostringstream mnz_edm;
		mnz_edm << "Edm = ";
		mnz_edm.precision(3);
		mnz_edm << hypo.GetEDM();
		fMapPaveTextSpectrum[pavetextname]->AddText(mnz_edm.str().c_str());
	}

}

/**
 * \brief Clean memory
 */
void START::PlotFactory::CleanPaveTextSpectrum()
{
  
	for(std::map<std::string,TPaveText*>::iterator pave=fMapPaveTextSpectrum.begin(); pave!=fMapPaveTextSpectrum.end(); pave++) {
		if(pave->second!=0) delete pave->second;
		pave->second=0;
	}
	fMapPaveTextSpectrum.clear();
}

/**
 * \brief Init Canvas with
 */
void START::PlotFactory::InitCanvasSpectrum(const Hypothesis &hypo)
{

	int canvaswidth(0), canvasheight(0); 

	switch(fPlotStyle) {
	case Default:
	case UserFriendly:
	case Paper:
		canvaswidth=800;
		canvasheight=600;
	}

	std::string canvasname = MakeKeyCanvasSpectrum(hypo.GetName());
	fMapCanvasSpectrumAndResiduals[canvasname] = new TCanvas(canvasname.c_str(),canvasname.c_str(),canvaswidth,canvasheight);

}


/**
 * \brief Clean memory
 */
void START::PlotFactory::CleanCanvasSpectrum()
{
  
	for(std::map<std::string,TCanvas*>::iterator canvas=fMapCanvasSpectrumAndResiduals.begin(); canvas!=fMapCanvasSpectrumAndResiduals.end(); canvas++) {
		if(canvas->second!=0) delete canvas->second;
		canvas->second=0;
	}
	fMapCanvasSpectrumAndResiduals.clear();
}

/**
 * \brief Initialize pads
 */
void START::PlotFactory::InitPads(const Hypothesis &hypo)
{

	// Definition of spectrum's and residuals' pad size

	switch(fPlotStyle) {
	case Default:
	case UserFriendly:
	case Paper:
		fx1padspec=0.01;
		fy1padspec=0.25;
		fx2padspec=0.99;
		fy2padspec=0.95;
		fx1padres=fx1padspec;
		fy1padres=0.01;
		fx2padres=fx2padspec;
		fy2padres=0.25;
	}
  
	// Spectrum

	std::string padspecname = MakeKeyPadSpectrum(hypo.GetName());
	fMapPadSpectrum[padspecname] = new TPad(padspecname.c_str(),padspecname.c_str(),fx1padspec,fy1padspec,fx2padspec,fy2padspec,0);
	fMapPadSpectrum[padspecname]->SetLogy();
	fMapPadSpectrum[padspecname]->SetLogx();
	fMapPadSpectrum[padspecname]->SetGridx();
	fMapPadSpectrum[padspecname]->SetGridy();
	//fMapPadSpectrum[padspecname]->SetFillColor(2);
	// Residuals
  
	std::string padresname = MakeKeyPadResiduals(hypo.GetName());
	fMapPadResiduals[padresname] = new TPad(padresname.c_str(),padresname.c_str(),fx1padres,fy1padres,fx2padres,fy2padres,0);
	fMapPadResiduals[padresname]->SetLogx();
	fMapPadResiduals[padresname]->SetGridx();
	fMapPadResiduals[padresname]->SetGridy();
	//fMapPadResiduals[padresname]->SetFillColor(3);

	switch(fPlotStyle) {
	case Default:
	case UserFriendly:
	case Paper:
		fMapPadSpectrum[padspecname]->SetTopMargin(0.04);
		fMapPadSpectrum[padspecname]->SetBottomMargin(0.04);
		fMapPadSpectrum[padspecname]->SetRightMargin(0.035);

		fMapPadResiduals[padresname]->SetTopMargin(0.08);
		fMapPadResiduals[padresname]->SetBottomMargin(0.30);
		fMapPadResiduals[padresname]->SetRightMargin(0.035);
	}

}

/**
 * \brief Clean memory
 */
void START::PlotFactory::CleanPads()
{
  
	for(std::map<std::string,TPad*>::iterator pad=fMapPadSpectrum.begin(); pad!=fMapPadSpectrum.end(); pad++) {
		if(pad->second!=0) delete pad->second;
		pad->second=0;
	}
	fMapPadSpectrum.clear();
	for(std::map<std::string,TPad*>::iterator pad=fMapPadResiduals.begin(); pad!=fMapPadResiduals.end(); pad++) {
		if(pad->second!=0) delete pad->second;
		pad->second=0;
	}
	fMapPadResiduals.clear();

}

/**
 * \brief Initialize TH2
 */
void START::PlotFactory::InitTH2(Hypothesis &hypo)
{
  
	// Get residuals points & Co
	// Wee take residuals & co to define th2's limits in energy anf flux.
  
	fresiduals.clear();
	fresidualssigmaplus.clear();
	fresidualssigmaminus.clear();
	fresiduals3sigmaplus.clear();
	fresiduals3sigmaminus.clear();
	fflux.clear();
	ffluxsigmaplus.clear();
	ffluxsigmaminus.clear();
	fflux3sigmaplus.clear();
	fflux3sigmaminus.clear();
	fenergybin.clear();
  
	fenergybin = hypo.GetMeanEnergy();
	fresiduals = hypo.GetResiduals();
	fresidualssigmaplus = hypo.GetResidualsSigmaPlus();
	fresidualssigmaminus = hypo.GetResidualsSigmaMinus();
	fresiduals3sigmaplus = hypo.GetResiduals3SigmaPlus();
	fresiduals3sigmaminus = hypo.GetResiduals3SigmaMinus();
	fflux = hypo.GetResidualsFlux();
	ffluxsigmaplus = hypo.GetResidualsFluxSigmaPlus();
	ffluxsigmaminus = hypo.GetResidualsFluxSigmaMinus();
	fflux3sigmaplus = hypo.GetResidualsFlux3SigmaPlus();
	fflux3sigmaminus = hypo.GetResidualsFlux3SigmaMinus();     

	// Determine spectrum's and residuals' energy range
	std::pair<double,double> erange = GetPlotEnergyRange(hypo);
	double emin(0.), emax(0.);
	emin = erange.first;
	emax = erange.second;
  
	// kill crazy points (points with an error to big or if residuals is zero). It tells which points are killed
	//BadPointsKiller(fenergybin, fresiduals, fresidualssigmaplus, fresidualssigmaminus,fresiduals3sigmaplus, 
	//		  fresiduals3sigmaminus, fflux, ffluxsigmaplus, ffluxsigmaminus, fflux3sigmaplus, fflux3sigmaminus,false);
  
	// Determination of spectrum's flux range
  
	std::pair<double,double> fluxrange = GetPlotFluxRange(hypo,fenergybin,fflux,ffluxsigmaminus,ffluxsigmaplus,fflux3sigmaplus);
	double fluxmin = fluxrange.first;
	double fluxmax = fluxrange.second;
  
	// Determine residuals' plot range
  
	double limitresidualsth2 = GetPlotResidualsMaximum(fresidualssigmaminus,fresidualssigmaplus);
  
	// Spectrum
  
	std::string th2specname = MakeKeyTH2Spectrum(hypo.GetName());
	fMapTH2Spectrum[th2specname] = new TH2D(th2specname.c_str(),th2specname.c_str(),fenergybin.size(),emin,emax,1000,fluxmin,fluxmax);
	fMapTH2Spectrum[th2specname]->SetStats(kFALSE);
	fMapTH2Spectrum[th2specname]->SetTitle("");

	switch(fPlotStyle) {
	case Default:
	case Paper:
	case UserFriendly:
		fMapTH2Spectrum[th2specname]->GetXaxis()->SetTickLength(0.05);
		fMapTH2Spectrum[th2specname]->GetXaxis()->SetLabelSize(0.00);
		fMapTH2Spectrum[th2specname]->GetXaxis()->SetTitleFont(42);
		fMapTH2Spectrum[th2specname]->GetXaxis()->SetLabelFont(42);

		fMapTH2Spectrum[th2specname]->GetYaxis()->SetTitle("dN/dE [cm^{-2}.s^{-1}.TeV^{-1}]");
		fMapTH2Spectrum[th2specname]->GetYaxis()->SetTitleSize(0.05);
		fMapTH2Spectrum[th2specname]->GetYaxis()->SetTitleOffset(0.81);
		fMapTH2Spectrum[th2specname]->GetYaxis()->SetTitleFont(42);
		fMapTH2Spectrum[th2specname]->GetYaxis()->SetLabelFont(42);
	}
  
	// Residuals
  
	std::string th2resname = MakeKeyTH2Residuals(hypo.GetName());
	fMapTH2Residuals[th2resname] = new TH2D(th2resname.c_str(),th2resname.c_str(),fenergybin.size(),emin,emax,1000,-limitresidualsth2,limitresidualsth2);
	fMapTH2Residuals[th2resname]->SetStats(kFALSE);
	fMapTH2Residuals[th2resname]->SetTitle("");
  
	switch(fPlotStyle) {
	case Default:
	case Paper:
	case UserFriendly:
		fMapTH2Residuals[th2resname]->GetXaxis()->SetTickLength(0.15);
		fMapTH2Residuals[th2resname]->GetXaxis()->SetLabelSize(0.13);
		fMapTH2Residuals[th2resname]->GetXaxis()->SetLabelOffset(0.037);
		fMapTH2Residuals[th2resname]->GetXaxis()->SetTitle("Energy [TeV]");
		fMapTH2Residuals[th2resname]->GetXaxis()->SetTitleSize(0.15);
		fMapTH2Residuals[th2resname]->GetXaxis()->SetTitleOffset(0.94);
		fMapTH2Residuals[th2resname]->GetXaxis()->SetNoExponent();
		fMapTH2Residuals[th2resname]->GetXaxis()->SetMoreLogLabels();
		fMapTH2Residuals[th2resname]->GetXaxis()->SetTitleFont(42);
		fMapTH2Residuals[th2resname]->GetXaxis()->SetLabelFont(42);
		fMapTH2Residuals[th2resname]->GetXaxis()->SetTickLength(0.15);

		fMapTH2Residuals[th2resname]->GetYaxis()->SetLabelSize(0.08);
		fMapTH2Residuals[th2resname]->GetYaxis()->SetNdivisions(10,3,0);
		fMapTH2Residuals[th2resname]->GetYaxis()->CenterTitle();
		fMapTH2Residuals[th2resname]->GetYaxis()->SetTitle("(n_{exp}-n_{th})/n_{th}");
		//fMapTH2Residuals[th2resname]->GetYaxis()->SetTitle("Residuals");
		fMapTH2Residuals[th2resname]->GetYaxis()->SetTitleSize(0.16);
		fMapTH2Residuals[th2resname]->GetYaxis()->SetTitleOffset(0.24);
		fMapTH2Residuals[th2resname]->GetYaxis()->SetTitleFont(42);
		fMapTH2Residuals[th2resname]->GetYaxis()->SetLabelFont(42);
	}

}

/**
 * \brief Clean memory
 */
void START::PlotFactory::CleanTH2()
{

	for(std::map<std::string,TH2D*>::iterator th2=fMapTH2Spectrum.begin(); th2!=fMapTH2Spectrum.end(); th2++) {
		if(th2->second!=0) delete th2->second;
		th2->second=0;
	}
	fMapTH2Spectrum.clear();

	for(std::map<std::string,TH2D*>::iterator th2=fMapTH2Residuals.begin(); th2!=fMapTH2Residuals.end(); th2++) {
		if(th2->second!=0) delete th2->second;
		th2->second=0;
	}
	fMapTH2Residuals.clear();

}


/**
 * \brief Initialize TGraphs and fill them 
 */
void START::PlotFactory::InitAndFillGraphs(const Hypothesis &hypo)
{

	// Get residuals points & Co
	// Wee take residuals & co to define th2's limits in energy and flux.

	std::cout << "salut "<< std::endl;
	
	fresiduals.clear();
	fresiduals_on.clear();
	fresiduals_off.clear();
	fresidualssigmaplus.clear();
	fresidualssigmaminus.clear();
	fresiduals3sigmaplus.clear();
	fresiduals3sigmaminus.clear();
	fflux.clear();
	ffluxsigmaplus.clear();
	ffluxsigmaminus.clear();
	fflux3sigmaplus.clear();
	fflux3sigmaminus.clear();
	fenergybin.clear();
  
	fenergybin = hypo.GetMeanEnergy();
	fresiduals = hypo.GetResiduals();
	// JLK ADD FOR ADA
	fresiduals_on = hypo.GetResidualsOn();
	fresiduals_off = hypo.GetResidualsOff();
	fresidualssigmaplus = hypo.GetResidualsSigmaPlus();
	fresidualssigmaminus = hypo.GetResidualsSigmaMinus();
	fresiduals3sigmaplus = hypo.GetResiduals3SigmaPlus();
	fresiduals3sigmaminus = hypo.GetResiduals3SigmaMinus();
	fflux = hypo.GetResidualsFlux();
	ffluxsigmaplus = hypo.GetResidualsFluxSigmaPlus();
	ffluxsigmaminus = hypo.GetResidualsFluxSigmaMinus();
	fflux3sigmaplus = hypo.GetResidualsFlux3SigmaPlus();
	fflux3sigmaminus = hypo.GetResidualsFlux3SigmaMinus();  
  
	// kill crazy points to fill graph points with an error to big or if residuals is zero). It doesn't tell which points are killed
	//BadPointsKiller(fenergybin, fresiduals, fresidualssigmaplus, fresidualssigmaminus,fresiduals3sigmaplus, 
	//		  fresiduals3sigmaminus, fflux, ffluxsigmaplus, ffluxsigmaminus, fflux3sigmaplus, fflux3sigmaminus,true);
  
	// take care of bin width for UserFriendly style
	Band InfoBand(hypo.GetBandArray()[0]);
	InfoBand.ClearBandInfo();
	InfoBand.AddInfoFromBands(hypo.GetBandArray(),true);

	// looking for the first selected energy bin
	unsigned int firstbin(0);
	for(std::vector<EnergyBin>::iterator bin=InfoBand.ebin.begin(); bin!=InfoBand.ebin.end(); ++bin) {
		if(bin->GetKeepBin()==1) break;
		firstbin++;
	}

	// Spectrum
  
	// init graph
	std::string graphspecname = MakeKeyGraphErrorsSpectrum(hypo.GetName());
	fMapGraphErrorsSpectrum[graphspecname] = new TGraphAsymmErrors();
	fMapGraphErrorsSpectrum[graphspecname]->SetName(graphspecname.c_str());
	fMapGraphErrorsSpectrum[graphspecname]->SetMarkerColor(kRed);
	fMapGraphErrorsSpectrum[graphspecname]->SetMarkerStyle(20);
	fMapGraphErrorsSpectrum[graphspecname]->SetMarkerSize(0.8);
	fMapGraphErrorsSpectrum[graphspecname]->SetLineWidth(1);
	fMapGraphErrorsSpectrum[graphspecname]->SetLineColor(kBlack);
  
	// arrows
	std::string arrowspecname = MakeKeyArrowsSpectrum(hypo.GetName());
  
	// Fill graph
	if(fenergybin.size()>0 && fflux.size()>0 && ffluxsigmaplus.size()>0 && ffluxsigmaminus.size()>0 && fflux3sigmaplus.size()>0
	   && fflux3sigmaminus.size()>0) {
		unsigned int graphpoint(0);
		for(unsigned int ipoint(0); ipoint<fenergybin.size(); ipoint++) {
			if(fresiduals[ipoint]>0. && ffluxsigmaminus[ipoint]!=0){ // JLK
				double binemin(InfoBand.ebin[ipoint+firstbin].GetEmin()); // taking into account first bin from bands
				double binemax(InfoBand.ebin[ipoint+firstbin].GetEmax()); // taking into account first bin from bands
				double emean = fenergybin[ipoint];
				fMapGraphErrorsSpectrum[graphspecname]->SetPoint(graphpoint,fenergybin[ipoint],fflux[ipoint]);
				//fMapGraphErrorsSpectrum[graphspecname]->SetPointError(graphpoint,0.,0.,fflux[ipoint]-ffluxsigmaminus[ipoint],ffluxsigmaplus[ipoint]-fflux[ipoint]);
				switch(fPlotStyle) {
				case Default:
				case Paper:
					fMapGraphErrorsSpectrum[graphspecname]->SetPointError(graphpoint,
																		  0.,
																		  0.,
																		  fflux[ipoint]-ffluxsigmaminus[ipoint],
																		  ffluxsigmaplus[ipoint]-fflux[ipoint]);
					break;
				case UserFriendly:
					fMapGraphErrorsSpectrum[graphspecname]->SetPointError(graphpoint,
																		  emean-binemin,
																		  binemax-emean,
																		  fflux[ipoint]-ffluxsigmaminus[ipoint],
																		  ffluxsigmaplus[ipoint]-fflux[ipoint]);
				}
				DEBUG_OUT_L(2) << "Spectrum : adding normal points en=" << fenergybin[ipoint] << " flux=" << fflux[ipoint] << std::endl;
				graphpoint++;
			}
			else { // upperlimits
				fMapArrowsSpectrum[arrowspecname].push_back(new TArrow(fenergybin[ipoint],fflux3sigmaplus[ipoint],
																	   fenergybin[ipoint],(fflux3sigmaplus[ipoint])*fscalearrowfactor,0.02,"|-|>"));
				DEBUG_OUT_L(2) << "Spectrum : adding upperlimit en= " << fenergybin[ipoint] << " upperlimit = " << fflux3sigmaplus[ipoint] << " bottom arrow = " << 
					(fflux3sigmaplus[ipoint]-fflux[ipoint])*fscalearrowfactor << std::endl;
				fMapArrowsSpectrum[arrowspecname].back()->SetLineWidth(1);
				fMapArrowsSpectrum[arrowspecname].back()->SetLineColor(kBlack);
				fMapArrowsSpectrum[arrowspecname].back()->SetFillColor(kBlack);
			}
		}
	}
	else {
		WARN_OUT << "Can't plot spectrum, something's wrong with number of points for " << hypo.GetName() << std::endl;
	}
  
	// Residuals
  
	// graph
	std::string graphresname = MakeKeyGraphErrorsResiduals(hypo.GetName());
	fMapGraphErrorsResiduals[graphresname] = new TGraphAsymmErrors();
  
	fMapGraphErrorsResiduals[graphresname]->SetName(graphresname.c_str());
	fMapGraphErrorsResiduals[graphresname]->SetMarkerColor(kRed);
	fMapGraphErrorsResiduals[graphresname]->SetMarkerStyle(20);
	fMapGraphErrorsResiduals[graphresname]->SetMarkerSize(0.8);
	fMapGraphErrorsResiduals[graphresname]->SetLineWidth(1);
	fMapGraphErrorsResiduals[graphresname]->SetLineColor(kBlack);

	// JLK ADD FOR ADA
	std::string graphresname_on = MakeKeyGraphErrorsResidualsOn(hypo.GetName());
	fMapGraphErrorsResidualsOn[graphresname_on] = new TGraphAsymmErrors();
  
	fMapGraphErrorsResidualsOn[graphresname_on]->SetName(graphresname_on.c_str());
	fMapGraphErrorsResidualsOn[graphresname_on]->SetMarkerColor(kBlue);
	fMapGraphErrorsResidualsOn[graphresname_on]->SetMarkerStyle(20);
	fMapGraphErrorsResidualsOn[graphresname_on]->SetMarkerSize(0.8);
	fMapGraphErrorsResidualsOn[graphresname_on]->SetLineWidth(1);
	fMapGraphErrorsResidualsOn[graphresname_on]->SetLineColor(kBlack);
	
	std::string graphresname_off = MakeKeyGraphErrorsResidualsOff(hypo.GetName());
	fMapGraphErrorsResidualsOff[graphresname_off] = new TGraphAsymmErrors();
  
	fMapGraphErrorsResidualsOff[graphresname_off]->SetName(graphresname_off.c_str());
	fMapGraphErrorsResidualsOff[graphresname_off]->SetMarkerColor(kGreen-3);
	fMapGraphErrorsResidualsOff[graphresname_off]->SetMarkerStyle(20);
	fMapGraphErrorsResidualsOff[graphresname_off]->SetMarkerSize(0.8);
	fMapGraphErrorsResidualsOff[graphresname_off]->SetLineWidth(1);
	fMapGraphErrorsResidualsOff[graphresname_off]->SetLineColor(kBlack);
	
	// arrows
	std::string arrowresname = MakeKeyArrowsResiduals(hypo.GetName());

	// fill the graph
	if(fenergybin.size()>0 && fresiduals.size()>0 && fresidualssigmaplus.size()>0 && fresidualssigmaminus.size()>0 
	   && fresiduals3sigmaplus.size()>0 && fresiduals3sigmaminus.size()>0) {
    
		unsigned int graphpoint(0);
		for(unsigned int ipoint(0); ipoint<fenergybin.size(); ipoint++) {
      
			if(fresiduals[ipoint]>0. && fresidualssigmaminus[ipoint]!=0) { // JLK
				double binemin(InfoBand.ebin[ipoint+firstbin].GetEmin()); // taking into account first bin from bands
				double binemax(InfoBand.ebin[ipoint+firstbin].GetEmax()); // taking into account first bin from bands
				double emean = fenergybin[ipoint];
				fMapGraphErrorsResiduals[graphresname]->SetPoint(graphpoint,fenergybin[ipoint],fresiduals[ipoint]-1);
				fMapGraphErrorsResidualsOn[graphresname_on]->SetPoint(graphpoint,fenergybin[ipoint],fresiduals_on[ipoint]-1);
				fMapGraphErrorsResidualsOff[graphresname_off]->SetPoint(graphpoint,fenergybin[ipoint],fresiduals_off[ipoint]-1);
				switch(fPlotStyle) {
				case Default:
				case Paper:
					fMapGraphErrorsResiduals[graphresname]->SetPointError(graphpoint,
																		  0.,
																		  0.,
																		  fresiduals[ipoint]-fresidualssigmaminus[ipoint],
																		  fresidualssigmaplus[ipoint]-fresiduals[ipoint]);
					fMapGraphErrorsResidualsOn[graphresname_on]->SetPointError(graphpoint,
																			   0.,
																			   0.,
																			   0.,
																			   0.);
					fMapGraphErrorsResidualsOff[graphresname_off]->SetPointError(graphpoint,
																				 0.,
																				 0.,
																				 0.,
																				 0.);
					break;
				case UserFriendly:
					fMapGraphErrorsResiduals[graphresname]->SetPointError(graphpoint,
																		  emean-binemin,
																		  binemax-emean,
																		  fresiduals[ipoint]-fresidualssigmaminus[ipoint],
																		  fresidualssigmaplus[ipoint]-fresiduals[ipoint]);
					fMapGraphErrorsResidualsOn[graphresname_on]->SetPointError(graphpoint,
																			   emean-binemin,
																			   binemax-emean,
																			   0.,
																			   0.);
					fMapGraphErrorsResidualsOff[graphresname_off]->SetPointError(graphpoint,
																				 emean-binemin,
																				 binemax-emean,
																				 0.,
																				 0.);

				}
				graphpoint++;
				DEBUG_OUT_L(2) << "Residuals : adding residuals points, point " << graphpoint 
							   << " : binemin=" << binemin 
							   << " ebin=" << fenergybin[ipoint] 
							   << " binemax=" << binemax 
							   << " res=" << fresiduals[ipoint]-1 << std::endl;
			}
			else { // upper limits, no upper limit on residuals JLK
	
			}
      
		}
    
	}
	else {
		WARNING << "Can't plot residuals, something's wrong with the number of points" << hypo.GetName() << std::endl;
	}

}

/**
 * \brief Clean TGraphAsymmErrors
 */
void START::PlotFactory::CleanGraphSpectrum()
{

	for(std::map<std::string,std::vector<TArrow*> >::iterator arrow=fMapArrowsResiduals.begin(); arrow!=fMapArrowsResiduals.end(); arrow++) {
		for(unsigned int i(0); i<fMapArrowsResiduals[arrow->first].size(); i++) {
			if(arrow->second[i]!=0) delete arrow->second[i];
			arrow->second[i]=0;
		}
		fMapArrowsResiduals[arrow->first].clear();
	}
	fMapArrowsResiduals.clear();

	for(std::map<std::string,std::vector<TArrow*> >::iterator arrow=fMapArrowsSpectrum.begin(); arrow!=fMapArrowsSpectrum.end(); arrow++) {
		for(unsigned int i(0); i<fMapArrowsSpectrum[arrow->first].size(); i++) {
			if(arrow->second[i]!=0) delete arrow->second[i];
			arrow->second[i]=0;
		}
		fMapArrowsSpectrum[arrow->first].clear();
	}
	fMapArrowsSpectrum.clear();

	for(std::map<std::string,TGraphAsymmErrors*>::iterator graph=fMapGraphErrorsSpectrum.begin(); graph!=fMapGraphErrorsSpectrum.end(); graph++) {
		if(graph->second!=0) delete graph->second;
		graph->second=0;
	}
	fMapGraphErrorsSpectrum.clear();

	for(std::map<std::string,TGraphAsymmErrors*>::iterator graph=fMapGraphErrorsResiduals.begin(); graph!=fMapGraphErrorsResiduals.end(); graph++) {
		if(graph->second!=0) delete graph->second;
		graph->second=0;
	}
	fMapGraphErrorsResiduals.clear();
}

/**
 * \brief Init PolyLine for the butterfly
 */
void START::PlotFactory::InitPolyLines(Hypothesis &hypo)
{

	unsigned int nbpoints = 400;
  
	// get energy range
	std::pair<double,double> erange = GetButterflyEnergyRange(hypo);
	double emin(0.), emax(0.);
	emin = erange.first;
	emax = erange.second;

	double log10emin = TMath::Log10(emin);
	double log10emax = TMath::Log10(emax);

	// butterfly contour's points
	std::vector<double> x;
	std::vector<double> y;

	// butterfly's vectors
	std::pair<std::vector<double>,std::vector<double> > pairvectors;

	std::string polylinename = MakeKeyPolyLineSpectrum(hypo.GetName());
	std::string polylinevectorsname = MakeKeyPolyLineVectorsSpectrum(hypo.GetName());

	switch(fButterfly) {

	case NoButterfly :

		break;

	case LinearButterfly : 
    
		ComputeLinearButterfly(hypo,x,y,nbpoints,log10emin,log10emax);

		pairvectors = std::make_pair(x,y);

		fMapButterFlyVectors[polylinevectorsname] = pairvectors;

		fMapPolyLineSpectrum[polylinename] = new TPolyLine(x.size());
		fMapPolyLineSpectrum[polylinename]->SetFillColor(kGreen-3);
		fMapPolyLineSpectrum[polylinename]->SetLineColor(1);
		fMapPolyLineSpectrum[polylinename]->SetLineWidth(1);
    
		// Fill butterfly
		for(unsigned int ipoint(0); ipoint<x.size(); ipoint++) fMapPolyLineSpectrum[polylinename]->SetPoint(ipoint,x[ipoint],y[ipoint]);
    
		// Copy it in Hypothesis
		hypo.AddButterfly(Hypothesis::LinearButterfly,*fMapPolyLineSpectrum[polylinename]);

		break;

	case LogarithmButterfly :

		ComputeLogarithmButterfly(hypo,x,y,nbpoints,log10emin,log10emax);

		pairvectors = std::make_pair(x,y);

		fMapButterFlyVectors[polylinevectorsname] = pairvectors;

		fMapPolyLineSpectrum[polylinename] = new TPolyLine(x.size());
		fMapPolyLineSpectrum[polylinename]->SetFillColor(kGreen-3);
		fMapPolyLineSpectrum[polylinename]->SetLineColor(1);
		fMapPolyLineSpectrum[polylinename]->SetLineWidth(1);
    
		// Fill butterfly
		for(unsigned int ipoint(0); ipoint<x.size(); ipoint++) fMapPolyLineSpectrum[polylinename]->SetPoint(ipoint,x[ipoint],y[ipoint]);

		// Copy it in Hypothesis
		hypo.AddButterfly(Hypothesis::LogarithmButterfly,*fMapPolyLineSpectrum[polylinename]);

		break;

	case CausticButterfly :

		ComputeCausticButterfly(hypo,log10emin,log10emax,nbpoints);

		break;

	case ContoursButterfly :

		ComputeContoursButterfly(hypo,x,y,log10emin,log10emax,nbpoints);

		pairvectors = std::make_pair(x,y);

		fMapButterFlyVectors[polylinevectorsname] = pairvectors;

		fMapPolyLineSpectrum[polylinename] = new TPolyLine(x.size());
		fMapPolyLineSpectrum[polylinename]->SetFillColor(kGreen-3);
		fMapPolyLineSpectrum[polylinename]->SetLineColor(1);
		fMapPolyLineSpectrum[polylinename]->SetLineWidth(1);
    
		// Fill butterfly
		for(unsigned int ipoint(0); ipoint<x.size(); ipoint++) fMapPolyLineSpectrum[polylinename]->SetPoint(ipoint,x[ipoint],y[ipoint]);

		// Copy it in Hypothesis
		hypo.AddButterfly(Hypothesis::ContoursButterfly,*fMapPolyLineSpectrum[polylinename]);

		break;

	default : 

		WARNING << "No butterfly options specified so you won't have butterfly!" << std::endl;

	}

}

/**
 * \brief Clean memory
 */
void START::PlotFactory::CleanPolyLines()
{
  
	for(std::map<std::string,TPolyLine*>::iterator polyline=fMapPolyLineSpectrum.begin(); polyline!=fMapPolyLineSpectrum.end(); polyline++) {
		if(polyline->second!=0) delete polyline->second;
		polyline->second=0;
	}
	fMapPolyLineSpectrum.clear();
  
	fMapButterFlyVectors.clear();

}

void START::PlotFactory::CleanGraphArraySpectrum() {

	for(std::map<std::string,std::vector<TGraph*> >::iterator map=fMapGraphArraySpectrum.begin(); map!=fMapGraphArraySpectrum.end(); map++) {
		for(unsigned int igraph(0); igraph !=map->second.size(); igraph++) {
			if(map->second[igraph]!=0) delete map->second[igraph];
			map->second[igraph]=0;
		}
	}
	fMapGraphArraySpectrum.clear();

}

/**
 * \brief compute butterfly with covariance of the flux at first order
 */
void START::PlotFactory::ComputeLinearButterfly(Hypothesis &hypo,std::vector<double> &x, std::vector<double> &y,
												unsigned int nbpoints, double log10emin, double log10emax) {

	x.clear(); y.clear();

	double binsize((log10emax-log10emin)/nbpoints);

	for(double ien=log10emin; ien<=log10emax*1.001; ien+=binsize) {
    
		double cov = hypo.GetFluxFitParams(TMath::Power(10,ien))+hypo.GetSigmaFlux(TMath::Power(10,ien));
		if(cov<0) {
			WARN_OUT << "Negative quantity with butterfly (phi+cov)!!!! Should be impossible!! " << std::endl;
      
			std::cout << "energy = " << TMath::Power(10,ien) << "flux = " << hypo.GetFluxFitParams(TMath::Power(10,ien)) 
					  << " covariance = " << hypo.GetSigmaFlux(TMath::Power(10,ien)) 
					  << " x=" << x.back() << " y=" << y.back() << std::endl;
			continue;
		}
		x.push_back(TMath::Power(10,ien));
		y.push_back(cov);
    
	}
  
	for(double ien=log10emax; ien>=log10emin*1.001; ien-=binsize) {
    
		double cov = hypo.GetFluxFitParams(TMath::Power(10,ien))-hypo.GetSigmaFlux(TMath::Power(10,ien));
		if(cov<=0) {
			if(fverbose) {
				WARN_OUT << "Negative quantity with butterfly (phi-cov)!!!! This point is set to zero : " << std::endl;
	
				std::cout << "energy = " << TMath::Power(10,ien) << "flux = " << hypo.GetFluxFitParams(TMath::Power(10,ien)) 
						  << " covariance = " << hypo.GetSigmaFlux(TMath::Power(10,ien)) 
						  << " x=" << x.back() << " y=" << y.back() << std::endl;
			}
		}
		x.push_back(TMath::Power(10,ien));
		y.push_back(cov);
    
	}

	// we close the polyline
	x.push_back(x[0]);
	y.push_back(y[0]);

}

/**
 * \brief compute butterfly with covariance of logarithm of the flux at first order
 */
void START::PlotFactory::ComputeLogarithmButterfly(Hypothesis &hypo,std::vector<double> &x, std::vector<double> &y,
												   unsigned int nbpoints, double log10emin, double log10emax) {

	x.clear(); y.clear();

	double binsize((log10emax-log10emin)/nbpoints);

	for(double ien=log10emin; ien<=log10emax*1.001; ien+=binsize) {
    
		double cov = 
			hypo.GetFluxFitParams(TMath::Power(10,ien))*TMath::Exp(+hypo.GetSigmaFlux(TMath::Power(10,ien))/hypo.GetFluxFitParams(TMath::Power(10,ien)));
		if(cov<0) {
			WARN_OUT << "Negative quantity with butterfly (phi+cov)!!!! Normally impossible! You should ctr+z; kill -9 %." << std::endl;
      
			std::cout << "energy = " << TMath::Power(10,ien) << "flux = " << hypo.GetFluxFitParams(TMath::Power(10,ien)) 
					  << " covariance = " << hypo.GetSigmaFlux(TMath::Power(10,ien)) 
					  << " x=" << x.back() << " y=" << y.back() << std::endl;
      
			continue;
		}
		x.push_back(TMath::Power(10,ien));
		y.push_back(cov);
    
	}
  
	for(double ien=log10emax; ien>=log10emin*1.001; ien-=binsize) {
    
		double cov = 
			hypo.GetFluxFitParams(TMath::Power(10,ien))*TMath::Exp(-hypo.GetSigmaFlux(TMath::Power(10,ien))/hypo.GetFluxFitParams(TMath::Power(10,ien)));
		if(cov<=0) {
			WARN_OUT << "Negative quantity with butterfly (phi-cov)!!!! Kill this point : " << std::endl;
      
			std::cout << "energy = " << TMath::Power(10,ien) << "flux = " << hypo.GetFluxFitParams(TMath::Power(10,ien)) 
					  << " covariance = " << hypo.GetSigmaFlux(TMath::Power(10,ien)) 
					  << " x=" << x.back() << " y=" << y.back() << std::endl;
			continue;
		}
		x.push_back(TMath::Power(10,ien));
		y.push_back(cov);
    
	}

	// we close the polyline
	x.push_back(x[0]);
	y.push_back(y[0]);

}

/**
 * \brief compute butterfly with contours.
 */
void START::PlotFactory::ComputeCausticButterfly(Hypothesis &hypo,double log10emin, double log10emax, unsigned int nbpoints) {

	if(hypo.GetParametersNb()>2) {
		WARNING<< "Butterfly with contours is only implemented for hypothesis with at most 2 parameters" << std::endl;
		INFO << "So you'll have nothing!" << std::endl;
		return;
	}

	// 1 sigma
	std::vector <std::pair <std::vector <std::pair <double,double> >, std::pair<int,int> > > vectorcontour1 = hypo.GetContourSigma1();
  
	if(vectorcontour1.size()==0) {
		WARNING << "Contours at 1 sigma are not computed OR don't seem really nice, so you won't have butterfly with contours..." << std::endl;
		return;
	}

	// JLK : Un peu hardcode pour le moment sachant que l'on ne sait pas faire les contours a D>2...
	// donc en gros ca ne marchera que pour la pwl... sooréé

	unsigned int index(999);

	std::vector<std::string> paramnames= hypo.GetParametersNames();

	bool paramreversed=false;

	// we look for phi0 and gamma for now...
	for(unsigned int icont(0); icont<vectorcontour1.size(); icont++) {

		if(paramnames[vectorcontour1[icont].second.first]=="phi0" && paramnames[vectorcontour1[icont].second.second]=="Gamma") {
			index=icont;
		} 
		else if(paramnames[vectorcontour1[icont].second.first]=="Gamma" && paramnames[vectorcontour1[icont].second.second]=="phi0") {
			index=icont;
			paramreversed=true;
		}
	}

	if(index==999) {
		WARNING << "I can't deal with your contours... Sorry no contour butterfly" << std::endl;
		return;
	}

	std::vector<std::pair<double,double> > pointscontour = vectorcontour1[index].first; // get contours points 

	unsigned int points = pointscontour.size();

	std::vector<TGraph*> TGraphArray;

	DEBUG_OUT_L(2) << "hypothesis = " << hypo.GetName() << std::endl;

	for(unsigned int igraph(0); igraph<points; igraph++) {
		TString graphname="igraphbutterflycontour_";
		graphname+=hypo.GetName();
		graphname+="_";
		graphname+=igraph;
		TGraphArray.push_back(new TGraph());
		TGraphArray.back()->SetName(graphname.Data());
		TGraphArray.back()->SetLineColor(kGreen-3);
		TGraphArray.back()->SetLineWidth(1);
		DEBUG_OUT_L(2) << graphname << std::endl;
		std::vector<double> paramsforhypo;
		if(!paramreversed) {
			paramsforhypo.push_back(pointscontour[igraph].first);
			paramsforhypo.push_back(pointscontour[igraph].second);
		}
		else {
			paramsforhypo.push_back(pointscontour[igraph].second);
			paramsforhypo.push_back(pointscontour[igraph].first);
		}

		hypo.SetParameters(paramsforhypo);
    
		double binsize((log10emax-log10emin)/nbpoints);
		unsigned int ipoint(0);
    
		for(double ien(log10emin); ien<=log10emax; ien+=binsize) {
			TGraphArray.back()->SetPoint(ipoint,TMath::Power(10,ien),hypo.GetFlux(TMath::Power(10,ien)));
			DEBUG_OUT_L(2) << "ipoint=" << ipoint << "x=" << TMath::Power(10,ien) << " y=" << hypo.GetFlux(TMath::Power(10,ien)) << std::endl;
			ipoint++;
		}
    
	}

	std::string mapname = MakeKeyGraphArraySpectrum(hypo.GetName());

	for(unsigned int igraph(0); igraph<TGraphArray.size(); igraph++) {
		fMapGraphArraySpectrum[mapname].push_back(new TGraph(*TGraphArray[igraph]));
	}
  
	// clean
	for(unsigned int igraph(0); igraph<TGraphArray.size(); igraph++) {
		delete TGraphArray[igraph];
		TGraphArray[igraph] = 0;
	}

}

/**
 * \brief compute butterfly in a "caustic" way.
 */
void START::PlotFactory::ComputeContoursButterfly(Hypothesis &hypo,std::vector<double> &x, std::vector<double> &y,
												  double log10emin, double log10emax, unsigned int nbpoints)
{

	x.clear(); y.clear();

	if(hypo.GetParametersNb()>2) {
		WARNING<< "Butterfly with contours is only implemented for hypothesis with at most 2 parameters" << std::endl;
		INFO << "So you'll have nothing!" << std::endl;
		return;
	}

	// 1 sigma
	std::vector <std::pair <std::vector <std::pair <double,double> >, std::pair<int,int> > > vectorcontour1 = hypo.GetContourSigma1();
	DEBUG_OUT << "vectorcontour1.size() = " << vectorcontour1.size() << std::endl;
	if(vectorcontour1.size()==0) {
		WARNING << "Contours at 1 sigma are not computed OR don't seem really nice, so you won't have butterfly with contours..." << std::endl;
		return;
	}

	// JLK : Un peu hardcode pour le moment sachant que l'on ne sait pas faire les contours a D>2...
	// donc en gros ca ne marchera que pour la pwl... sooréé

	unsigned int index(999);

	std::vector<std::string> paramnames= hypo.GetParametersNames();

	bool paramreversed=false;

	// we look for phi0 and gamma for now...
	for(unsigned int icont(0); icont<vectorcontour1.size(); icont++) {

		if(paramnames[vectorcontour1[icont].second.first]=="phi0" && paramnames[vectorcontour1[icont].second.second]=="Gamma") {
			index=icont;
		} 
		else if(paramnames[vectorcontour1[icont].second.first]=="Gamma" && paramnames[vectorcontour1[icont].second.second]=="phi0") {
			index=icont;
			paramreversed=true;
		}
	}

	if(index==999) {
		WARNING << "I can't deal with your contours... Sorry no contour butterfly" << std::endl;
		return;
	}

	std::vector<std::pair<double,double> > pointscontour = vectorcontour1[index].first; // get contours points 

	if(pointscontour.size()<3) {
		WARNING << "Contours points are too low (" << pointscontour.size() << "), so no butterfly" << std::endl;
		return;
	}

	unsigned int npoints = pointscontour.size();

	double binsize((log10emax-log10emin)/nbpoints);

	for(double ien=log10emin; ien<=log10emax*1.001; ien+=binsize) {
    
		double butt(0);
		std::vector<double> flux;

		for(unsigned int ipoint(0); ipoint<npoints; ipoint++) {

			std::vector<double> paramsforhypo;
			if(!paramreversed) {
				paramsforhypo.push_back(pointscontour[ipoint].first);
				paramsforhypo.push_back(pointscontour[ipoint].second);
			}
			else {
				paramsforhypo.push_back(pointscontour[ipoint].second);
				paramsforhypo.push_back(pointscontour[ipoint].first);
			}
			hypo.SetParameters(paramsforhypo);

			flux.push_back(hypo.GetFlux(TMath::Power(10,ien)));

		}

		std::vector<double>::const_iterator it;
		it = max_element(flux.begin(), flux.end());
		butt = *it;

		x.push_back(TMath::Power(10,ien));
		y.push_back(butt);
    
	}
  
	for(double ien=log10emax; ien>=log10emin*1.001; ien-=binsize) {
    
		double butt(0);
		std::vector<double> flux;

		for(unsigned int ipoint(0); ipoint<npoints; ipoint++) {

			std::vector<double> paramsforhypo;
			if(!paramreversed) {
				paramsforhypo.push_back(pointscontour[ipoint].first);
				paramsforhypo.push_back(pointscontour[ipoint].second);
			}
			else {
				paramsforhypo.push_back(pointscontour[ipoint].second);
				paramsforhypo.push_back(pointscontour[ipoint].first);
			}
			hypo.SetParameters(paramsforhypo);

			flux.push_back(hypo.GetFlux(TMath::Power(10,ien)));

		}

		std::vector<double>::const_iterator it;
		it = min_element(flux.begin(), flux.end());
		butt = *it;

		x.push_back(TMath::Power(10,ien));
		y.push_back(butt);
    
	}

	// we close the polyline
	x.push_back(x[0]);
	y.push_back(y[0]);

}


/**
 * \brief Initialize Canvas
 */
//void START::PlotFactory::InitTF1Residuals(const Hypothesis &hypo)
void START::PlotFactory::InitTF1Residuals(const Hypothesis &hypo)
{

	std::pair<double,double> erange = GetPlotEnergyRange(hypo);
	double emin(erange.first), emax(erange.second);
  
	std::string tf1name = MakeKeyTF1Residuals(hypo.GetName());
	fMapTF1Residuals[tf1name] = new TF1(tf1name.c_str(),"0",emin,emax);
	fMapTF1Residuals[tf1name]->SetLineColor(kBlack);
	fMapTF1Residuals[tf1name]->SetLineWidth(1);

}

/**
 * \brief Initialize spectrum TF1 (fit)
 */
void START::PlotFactory::InitTF1FitSpectrum(const Hypothesis &hypo)
{

	std::pair<double,double> erange = GetPlotEnergyRange(hypo);
	double emin(erange.first), emax(erange.second);  
  
	Hypothesis *hypoclone = hypo.clone(); // JLK : memory leak but i haven't a clean way

	hypoclone->SetParameters(hypo.GetFittedParameters());

	std::string tf1name = MakeKeyTF1FitSpectrum(hypo.GetName());
	fMapTF1FitSpectrum[tf1name] = new TF1(tf1name.c_str(),hypoclone,emin,emax,0,"Hypothesis");
	fMapTF1FitSpectrum[tf1name]->SetLineColor(kBlue);
	fMapTF1FitSpectrum[tf1name]->SetLineWidth(1);
	fMapTF1FitSpectrum[tf1name]->SetRange(hypo.GetMinimizationEnergyRange().first,hypo.GetMinimizationEnergyRange().second);
	fMapTF1FitSpectrum[tf1name]->SetNpx(500); // fix number of points for graphical representation
}

/**
 * \brief Clean memory
 */
void START::PlotFactory::CleanTF1FitSpectrum()
{
  
	for(std::map<std::string,TF1*>::iterator tf1=fMapTF1FitSpectrum.begin(); tf1!=fMapTF1FitSpectrum.end(); tf1++) {
		if(tf1->second!=0) delete tf1->second;
		tf1->second=0;
	}
	fMapTF1FitSpectrum.clear();
}

/**
 * \brief Clean memory
 */
void START::PlotFactory::CleanTF1Residuals()
{
  
	for(std::map<std::string,TF1*>::iterator tf1=fMapTF1Residuals.begin(); tf1!=fMapTF1Residuals.end(); tf1++) {
		if(tf1->second!=0) delete tf1->second;
		tf1->second=0;
	}
	fMapTF1Residuals.clear();
}


/**
 * \brief Draw ROOT's objetcs
 */
void START::PlotFactory::DrawSpectrumAndResiduals(const Hypothesis &hypo)
{
	
    switch(fButterfly) {
    case NoButterfly:
		DEBUG_OUT << "butterfly = nope" << std::endl;
		break;
    case LinearButterfly:
		DEBUG_OUT << "butterfly = linear" << std::endl;
		break;
    case LogarithmButterfly:
		DEBUG_OUT << "butterfly = logarithm" << std::endl;
		break;
    case ContoursButterfly:
		DEBUG_OUT << "butterfly = contours" << std::endl;
		break;
    case CausticButterfly:
		DEBUG_OUT << "butterfly = caustic" << std::endl;
    }

	// draw maincanvas
  
	std::string canvasname = MakeKeyCanvasSpectrum(hypo.GetName());
	DEBUG_OUT << "canvas = " << canvasname << std::endl;
	fMapCanvasSpectrumAndResiduals[canvasname]->Draw();
	fMapCanvasSpectrumAndResiduals[canvasname]->cd();
  
	//// draw Spectrum
  
	std::string padspecname = MakeKeyPadSpectrum(hypo.GetName());
	DEBUG_OUT << "pad = " << padspecname << std::endl;
	fMapPadSpectrum[padspecname]->Draw(""); // draw spectrum's pad
	fMapPadSpectrum[padspecname]->cd();
  
	std::string th2specname = MakeKeyTH2Spectrum(hypo.GetName());
	DEBUG_OUT << "th2 = " << th2specname << std::endl;
	fMapTH2Spectrum[th2specname]->Draw(""); // draw spectrum's pad
  

	// butterfly

	if(fButterfly==LinearButterfly || fButterfly==LogarithmButterfly || fButterfly==ContoursButterfly) {
		std::string polyname = MakeKeyPolyLineSpectrum(hypo.GetName());
		DEBUG_OUT << "polyline = " << polyname << std::endl;
		if(fMapPolyLineSpectrum[polyname]!=0) { // if there is no butterfly, we don't draw it...
			fMapPolyLineSpectrum[polyname]->Draw("f"); // draw spectrum's butterfly
			fMapPolyLineSpectrum[polyname]->Draw(""); // draw spectrum's butterfly
		}
	}
	else if(fButterfly==CausticButterfly) {
		std::string grapharrayname = MakeKeyGraphArraySpectrum(hypo.GetName());
		DEBUG_OUT << "grapharray = " << grapharrayname << std::endl;
		DEBUG_OUT << "map size for bitterfly contour = " << fMapGraphArraySpectrum[grapharrayname].size() << std::endl;
		for(unsigned int igraph(0); igraph < fMapGraphArraySpectrum[grapharrayname].size(); igraph++) {
			fMapGraphArraySpectrum[grapharrayname][igraph]->Draw("L"); // draw spectrum butterfly
			DEBUG_OUT << "Drawnig TGraph " << igraph << " for butterfly" << std::endl;
		}
	}
	else {
		// no butterfly
	}

	// draw fit
	fMapPadSpectrum[padspecname]->cd();
	std::string tf1namespec = MakeKeyTF1FitSpectrum(hypo.GetName());
	DEBUG_OUT << "tf1 = " << tf1namespec << std::endl;
	if(fMapTF1FitSpectrum[tf1namespec]!=0 && fdrawfit) fMapTF1FitSpectrum[tf1namespec]->Draw("same"); // draw TF1

	// draw points
	std::string graphspecname = MakeKeyGraphErrorsSpectrum(hypo.GetName());
	DEBUG_OUT << "graph = " << graphspecname << std::endl;
	fMapGraphErrorsSpectrum[graphspecname]->Draw("P same"); // draw spectrum's graph

	// draw upper limits
	std::string arrowspecname = MakeKeyArrowsSpectrum(hypo.GetName());
	DEBUG_OUT << "arrows = " << arrowspecname << std::endl;
	for(unsigned int i(0); i<fMapArrowsSpectrum[arrowspecname].size(); i++) {
		if(fMapArrowsSpectrum[arrowspecname][i]!=0) fMapArrowsSpectrum[arrowspecname][i]->Draw(""); // draw spectrum's upperlimit
	}
  
	//draw Residuals
  
	fMapCanvasSpectrumAndResiduals[canvasname]->cd();
  
	std::string padresname = MakeKeyPadResiduals(hypo.GetName());
	DEBUG_OUT << "padresname = " << padresname << std::endl;
	fMapPadResiduals[padresname]->Draw(""); // draw residuals' pad
	fMapPadResiduals[padresname]->cd();
  
	std::string th2resname = MakeKeyTH2Residuals(hypo.GetName());
	DEBUG_OUT << "th2 = " << th2resname << std::endl;
	fMapTH2Residuals[th2resname]->Draw(""); // draw residual's th2

	// JLK add for ADA
	fMapPadResiduals[padresname]->cd();
	std::string graphresname = MakeKeyGraphErrorsResiduals(hypo.GetName());
	std::string graphresname_on = MakeKeyGraphErrorsResidualsOn(hypo.GetName());
	std::string graphresname_off = MakeKeyGraphErrorsResidualsOff(hypo.GetName());
	TLegend *LegRes = 0;
	LegRes = new TLegend(0.9,0.8,1.,1.);
	LegRes->SetLineColor(0);
	LegRes->SetFillStyle(1001);
	// blabla
	switch (fResidualsStyle) {
	case Excess:
		fMapGraphErrorsResiduals[graphresname]->Draw("P"); // draw residual's graph
		break;
	case ExcessOff:
		fMapGraphErrorsResiduals[graphresname]->Draw("P"); // draw residual's graph
		fMapGraphErrorsResidualsOff[graphresname_off]->Draw("P same"); // draw residual's graph
		LegRes->AddEntry(fMapGraphErrorsResiduals[graphresname],"Excess","p");
		LegRes->AddEntry(fMapGraphErrorsResidualsOff[graphresname_off],"OFF","p");
		LegRes->Draw("same");
		break;
	case OnOff:
		fMapGraphErrorsResidualsOn[graphresname_on]->Draw("P"); // draw residual's graph
		fMapGraphErrorsResidualsOff[graphresname_off]->Draw("P same"); // draw residual's graph
		LegRes->AddEntry(fMapGraphErrorsResidualsOn[graphresname_on],"ON","p");
		LegRes->AddEntry(fMapGraphErrorsResidualsOff[graphresname_off],"OFF","p");
		LegRes->Draw("same");
		break;
	}
	// std::string graphresname = MakeKeyGraphErrorsResiduals(hypo.GetName());
	// DEBUG_OUT << "graph = " << graphresname << std::endl;
	// fMapGraphErrorsResiduals[graphresname]->Draw("P"); // draw residual's graph
  
	fMapPadResiduals[padresname]->cd();
	std::string tf1name = MakeKeyTF1Residuals(hypo.GetName());
	DEBUG_OUT << "tf1 = " << tf1name << std::endl;
	fMapTF1Residuals[tf1name]->Draw("same"); // draw TF1
  
	std::string arrowresname = MakeKeyArrowsResiduals(hypo.GetName()); // JLK : we do not draw upperlimits on residuals
	for(unsigned int i(0); i<fMapArrowsResiduals[arrowresname].size(); i++) {
		//if(fMapArrowsResiduals[arrowresname][i]!=0) fMapArrowsResiduals[arrowresname][i]->Draw(""); // draw residual's arrows
	}
  
	// draw TPaveText

	switch(fPlotStyle) {
	case Paper:
		break;
	case Default:
	case UserFriendly:
		fMapPadSpectrum[padspecname]->cd();
		std::string pavetextname = MakeKeyPaveTextSpectrum(hypo.GetName());
		if(fMapPaveTextSpectrum[pavetextname]!=0) fMapPaveTextSpectrum[pavetextname]->Draw();
		break;
	}

	fMapCanvasSpectrumAndResiduals[canvasname]->cd();

	TLatex *drawsourcename = new TLatex();

	switch(fPlotStyle) {
	case Paper:
		break;
	case Default:
	case UserFriendly:
		fMapPadSpectrum[padspecname]->cd();
		drawsourcename->SetNDC();
		drawsourcename->SetTextSize(0.037);
		drawsourcename->SetTextColor(1);
		drawsourcename->SetTextFont(52);
		drawsourcename->DrawLatex(0.1,0.96,fsourcename.c_str());
		break;
	}


}

/**
 * \brief Function to plot an hypothesis.
 *
 * \param hypo Hypothesis
 */
void START::PlotFactory::PlotHypothesis(Hypothesis &hypo) {
  
	if(hypo.GetSpectralType()!=Hypothesis::Differential) return; // JLK : we don't do plots for integrated or energyflux

	InitCanvasSpectrum(hypo);
	InitAndFillGraphs(hypo);
	InitPolyLines(hypo);
	InitPads(hypo);
	InitTH2(hypo);
	InitTF1Residuals(hypo);
	if(fdrawfit) InitTF1FitSpectrum(hypo);
	InitPaveTextSpectrum(hypo);
	DrawSpectrumAndResiduals(hypo);

}

void START::PlotFactory::PlotLikelihoodScans() {

	CleanCanvasLikelihoodScan();
	CleanGraphScanLikelihood();

	for(std::vector<Hypothesis*>::iterator hypo=fHypothesisArray.begin(); hypo!=fHypothesisArray.end(); hypo++) {

		if((*hypo)->GetConvergence()) {

			InitCanvasLikelihoodScan(**hypo);
			InitGraphLikelihoodScan(**hypo);      
			DrawLikelihoodScans(**hypo);

		}
		else {
			WARNING << "Skip hypothesis " << (*hypo)->GetName() << " because it didn't converge!" << std::endl;
			continue;
		}

	}

	INFO << "Plotting Likelihood scans... ok" << "\033[0m" << std::endl;
}

/**
 * \brief Draw ROOT's objetcs from likelihood scan
 */
void START::PlotFactory::DrawLikelihoodScans(const Hypothesis &hypo) {

	for(unsigned int iscan(0); iscan<hypo.GetScansLikelihood().size(); iscan++) {
    
		DEBUG_OUT << "Enter in loop for scans" << std::endl;

		std::string canvasname = MakeKeyCanvasLikelihoodScan(hypo.GetName(),hypo.GetParametersNames()[hypo.GetScansLikelihood()[iscan].first]);
		fMapCanvasLikelihoodScan[canvasname]->Draw();

		fMapCanvasLikelihoodScan[canvasname]->cd();

		std::string graphname = MakeKeyGraphLikelihoodScan(hypo.GetName(),hypo.GetParametersNames()[hypo.GetScansLikelihood()[iscan].first]);
		fMapGraphLikelihoodScan[graphname]->Draw("AP");
		std::string xtitle = hypo.GetLatexParametersNames()[hypo.GetScansLikelihood()[iscan].first];
		if(hypo.GetParametersUnits()[hypo.GetScansLikelihood()[iscan].first]!="") {
			xtitle += " (";
			xtitle += hypo.GetParametersUnits()[hypo.GetScansLikelihood()[iscan].first];
			xtitle += ")";
		}
		switch(fPlotStyle) {
		case Default:
		case Paper:
		case UserFriendly:
			fMapGraphLikelihoodScan[graphname]->GetHistogram()->GetXaxis()->SetTitle(xtitle.c_str());
			fMapGraphLikelihoodScan[graphname]->GetHistogram()->GetXaxis()->CenterTitle();
			fMapGraphLikelihoodScan[graphname]->GetHistogram()->GetYaxis()->SetLabelSize(0.03);
      
			fMapGraphLikelihoodScan[graphname]->GetHistogram()->GetYaxis()->SetTitle("log L");
			fMapGraphLikelihoodScan[graphname]->GetHistogram()->GetYaxis()->CenterTitle();
			fMapGraphLikelihoodScan[graphname]->GetHistogram()->GetYaxis()->SetLabelSize(0.02);
		}

		TMarker *Marker = 0;
		Marker = new TMarker();
		Marker->SetX(hypo.GetFittedParameters()[hypo.GetScansLikelihood()[iscan].first]);
		Marker->SetY(hypo.GetMaximumLikelihood());
		Marker->SetMarkerStyle(5);
		Marker->SetMarkerColor(kBlue);
		Marker->SetMarkerSize(1.5);

		Marker->Draw();

	}

}

/**
 * \brief Init Likelihood scan TCanvas
 */
void START::PlotFactory::InitCanvasLikelihoodScan(const Hypothesis &hypo) {

	int canvaswidth(800), canvasheight(700); 

	DEBUG_OUT << "size scan = " << hypo.GetScansLikelihood().size() << std::endl;

	for(unsigned int iscan(0); iscan<hypo.GetScansLikelihood().size(); iscan++) {
    
		DEBUG_OUT << "Canvas" << std::endl;

		std::string canvasname = MakeKeyCanvasLikelihoodScan(hypo.GetName(),hypo.GetParametersNames()[hypo.GetScansLikelihood()[iscan].first]);
		fMapCanvasLikelihoodScan[canvasname] = new TCanvas(canvasname.c_str(),canvasname.c_str(),canvaswidth,canvasheight);
		fMapCanvasLikelihoodScan[canvasname]->SetGridy();
		fMapCanvasLikelihoodScan[canvasname]->SetGridx();
	}

}

/**
 * \brief Clean memory
 */
void START::PlotFactory::CleanCanvasLikelihoodScan()
{

	for(std::map<std::string,TCanvas*>::iterator canvas=fMapCanvasLikelihoodScan.begin(); canvas!=fMapCanvasLikelihoodScan.end(); canvas++) {
		if(canvas->second!=0) delete canvas->second;
		canvas->second=0;
	}
	fMapCanvasLikelihoodScan.clear();

}

/**
 * \brief Init Likelihood scan TGraph
 */
void START::PlotFactory::InitGraphLikelihoodScan(const Hypothesis &hypo) {

	DEBUG_OUT << "Enter" << std::endl;

	for(unsigned int iscan(0); iscan<hypo.GetScansLikelihood().size(); iscan++) {
    
		std::vector<double> x, y;

		for(unsigned int ipoint(0); ipoint<hypo.GetScansLikelihood()[iscan].second.size(); ipoint++) {
			x.push_back(hypo.GetScansLikelihood()[iscan].second[ipoint].first);
			y.push_back(hypo.GetScansLikelihood()[iscan].second[ipoint].second);
		}

		if(x.size()!=y.size()) continue;

		std::string graphname = MakeKeyGraphLikelihoodScan(hypo.GetName(),hypo.GetParametersNames()[hypo.GetScansLikelihood()[iscan].first]);
		fMapGraphLikelihoodScan[graphname] = new TGraph(x.size());
		fMapGraphLikelihoodScan[graphname]->SetName(graphname.c_str());
		std::string graphtitle = "Likelihood profile for parameter ";
		graphtitle+=hypo.GetLatexParametersNames()[hypo.GetScansLikelihood()[iscan].first];
		fMapGraphLikelihoodScan[graphname]->SetTitle(graphtitle.c_str());
		fMapGraphLikelihoodScan[graphname]->SetMarkerColor(kRed);
		fMapGraphLikelihoodScan[graphname]->SetMarkerStyle(20);
		fMapGraphLikelihoodScan[graphname]->SetMarkerSize(0.6);
		fMapGraphLikelihoodScan[graphname]->SetLineWidth(1);
		fMapGraphLikelihoodScan[graphname]->SetLineColor(kBlack);

		DEBUG_OUT << "Scan for param " << hypo.GetParametersNames()[hypo.GetScansLikelihood()[iscan].first] << " :"<< std::endl;

		for(unsigned int ipoint(0); ipoint<x.size(); ipoint++) {
			fMapGraphLikelihoodScan[graphname]->SetPoint(ipoint,x[ipoint],-y[ipoint]);
			DEBUG_OUT << "point " << ipoint << " x=" << x[ipoint] << " y=" << y[ipoint] << std::endl;
		}

	}

}

/**
 * \brief Clean memory
 */
void START::PlotFactory::CleanGraphScanLikelihood()
{
	for(std::map<std::string,TGraph*>::iterator graph=fMapGraphLikelihoodScan.begin(); graph!=fMapGraphLikelihoodScan.end(); graph++) {
		if(graph->second!=0) delete graph->second;
		graph->second=0;
	}
	fMapGraphLikelihoodScan.clear();
}

/**
 * \brief Routine to rebin the hypothesis, to make a nice plot (does not touch the real data)
 *
 * Principle : Create a copy of the hypothesis, change the BandVector by a summarize band and rebin given a minimum significance per bin
 */
void START::PlotFactory::RebinHypothesis(Hypothesis &hypo, double sigrebin,ButterflyOption butterfly, bool drawfit) {

	fButterfly = butterfly;
	fdrawfit = drawfit;


	if(hypo.GetSpectralType()!=Hypothesis::Differential) return; // JLK we don't rebin EnergyFlux or Integrated

	if(hypo.GetConvergence()) {

		INFO << "Rebinning hypothesis " << hypo.GetName() << " (" << sigrebin << " sigma)" << std::endl;
    
		// Sanity check
		if (hypo.GetBandArray().size()==0) {
			WARNING << "Can't rebin anything since the bandarray has no size" << std::endl;
			return;
		}
    
		// Copy the first element of band array (just to have a starting Band with the good format)
		Band BandInfo(hypo.GetBandArray()[0]); // Information from the initial bin
		BandInfo.ClearBandInfo();
		BandInfo.SetKeepBand(1);
    
		// We copy the hypothesis to have the rebin hypothesis
		Hypothesis *hypo_rebin = hypo.clone();
    
		std::ostringstream oss;
		oss << sigrebin << "sigmaperbin";
    
		TString hypo_newname_str_suffix = "_rebin_";
		hypo_newname_str_suffix+=TString(oss.str().c_str());
		TString hypo_newname_str = hypo.GetName()+hypo_newname_str_suffix;
		oss.str("");

		hypo_rebin->SetName(hypo_newname_str);
		hypo_rebin->SetTitle(hypo_newname_str);
		hypo_rebin->GetBandArray().clear(); // Security Reset
    
		// We're going to stack the different band inside one band
		//std::vector<Band> hypo_bandarray = hypo.GetBandArray();
		// VIM : I asolutely don't understand why I have to do this... I must admit that I don't really like 
		//it (maybe there is a C++ limitation that I am not aware of...). If I don't do that, it will complain... 
		//In fact the problem appear in the iterator below. But I don't understant.........
    
		BandInfo.AddInfoFromBands(hypo.GetBandArray(),true);
    
		// VIM :   We're going to tag which bin I have to group
		Band BandInfoRebin(BandInfo); // VIM : This should be improve, it is a little bit annoying to use this trick
		std::vector<Band> BandRebinArray;
		BandRebinArray.push_back(BandInfoRebin);
    
		BandRebinArray[0].ebin.clear(); // We remove the ebin vector, to replace it by the rebin vector we're going to compute.
    
		//VIM : The principle is that we compute the significance of one the first bin. If this does not fulfill 
		//the sigma criteria, then we combine this bin with the next one. And so 
		//on until the criteria is fulfilled or the end of the non rebin ebin array is reached.
		//If the next bin in the process is rejected (KeepBin==0) then we stop the rebin, even if the criteria is 
		//not obtained (this should never happen in fact, this is a sanity check).
    
		int bin_first=0;
		int bin_second=0;
		while (bin_first<(int)BandInfo.ebin.size()) {
      
			if (BandInfo.ebin[bin_first].GetKeepBin()==0) {
				++bin_first;
				continue;
			}
      
			double sig=-1;
      
			EnergyBin RegroupBin = EnergyBin();
			RegroupBin.SetEmin( BandInfo.ebin[bin_first].GetEmin() );
			RegroupBin.SetKeepBin(1);
      
			// VIM : I take data from the current bin, if this is not enought to fullfill te sigrebin condition, then I use the next bin, etc...
			bin_second = bin_first;
      
			while(sig<sigrebin && bin_second<(int)BandInfo.ebin.size() && BandInfo.ebin[bin_second].GetKeepBin()!=0) {
	
				RegroupBin.AddInfoFromEBin(BandInfo.ebin[bin_second]);
				RegroupBin.SetEmax(BandInfo.ebin[bin_second].GetEmax());
	
				if (RegroupBin.GetAlpha()!=0) {
					sig = STARTUtils::LiMaSignificance(int(RegroupBin.GetOn()),int(RegroupBin.GetOff()),RegroupBin.GetAlpha());
				}
				else {
					sig=0.;
				}
	
				++bin_second;
			}
      
			// VIM : We complete the information of the bin (emid and emean), then we add this bin to the Band.
			RegroupBin.SetEmid();
			double emean_bin = hypo_rebin->GetMeanBinEnergy(RegroupBin.GetEmin(),RegroupBin.GetEmax());
			RegroupBin.SetEmean(emean_bin);
			BandRebinArray[0].ebin.push_back(RegroupBin);
      
			bin_first = bin_second;
		}

		BandRebinArray[0].SetKeepBand(1);
		hypo_rebin->SetBandArray(BandRebinArray);
    
		INFO_OUT << "Band after Rebinning :" << std::endl;
		BandRebinArray[0].Print();
    
		// We compute the residuals
		ResidualsFactory *residu_rebin = new ResidualsFactory(*hypo_rebin);
		residu_rebin->ComputeResidualsRolke();
		delete residu_rebin;
    
		fRebinnedHypothesisArray.push_back(hypo_rebin); // we add a component in the rebinhypo vector
    
		// We copy the rebinned residuals and add it in hypothesis
		Residuals Res = hypo_rebin->GetMapResiduals()[Hypothesis::MakeKeyResiduals()]; // as it is a unrebinned
		hypo.AddResiduals(Res,sigrebin);

		// we plot the hypothesis
		PlotHypothesis(*hypo_rebin);
    
	}
	else {
		INFO << "Skip rebin for hypothesis " << hypo.GetName() << " because it didn't converge!" << std::endl;
	}

}
  

/**
 * \brief Kill bad points wich have nothing to do on the spectrum and residuals plot 
 *
 * If we have an upper limit AND the assiociated 3sigma+ irrevelant on the residual, 
 * or no residuals, or a crazy point we kill the point.
 *
 * \param quiet don't write badpoints if it's true (the function is used one time for 
 * TH2 initialization and a second time for the TGraph and we need only once the 
 * informations)
 */
void START::PlotFactory::BadPointsKiller(std::vector<double> &energybin, 
										 std::vector<double> &residuals, 
										 std::vector<double> &residualssigmaplus, 
										 std::vector<double> &residualssigmaminus,
										 std::vector<double> &residuals3sigmaplus, 
										 std::vector<double> &residuals3sigmaminus, 
										 std::vector<double> &flux, 
										 std::vector<double> &fluxsigmaplus, 
										 std::vector<double> &fluxsigmaminus, 
										 std::vector<double> &flux3sigmaplus, 
										 std::vector<double> &flux3sigmaminus,
										 bool quiet)
{

	std::vector<double> energybin_tmp;
	std::vector<double> residuals_tmp; 
	std::vector<double> residualssigmaplus_tmp; 
	std::vector<double> residualssigmaminus_tmp; 
	std::vector<double> residuals3sigmaplus_tmp; 
	std::vector<double> residuals3sigmaminus_tmp; 
	std::vector<double> flux_tmp;
	std::vector<double> fluxsigmaplus_tmp; 
	std::vector<double> fluxsigmaminus_tmp; 
	std::vector<double> flux3sigmaplus_tmp; 
	std::vector<double> flux3sigmaminus_tmp;

	double limitpoint(50.);

	for(unsigned int ipoint(0); ipoint<energybin.size(); ipoint++) {
		if(residuals[ipoint]==0) { // JLK : no residuals so we delete those points
			if(!quiet) {
				WARN_OUT << "This point is killed! (no residuals)" << std::endl;
				INFO_OUT << "e_mean" << " " << "residual_exp" << " " << "residual_sigma_minus" << " " 
						 << "residual_sigma_plus" << std::endl;
				INFO_OUT << "1sigma : " << energybin[ipoint] << " " << residuals[ipoint] -1 << " " 
						 << residualssigmaminus[ipoint] -1 << " " 
						 << residualssigmaplus[ipoint] -1 << " " << std::endl;
				INFO_OUT << "3sigma : " << energybin[ipoint] << " " << residuals[ipoint] -1 << " " 
						 << residuals3sigmaminus[ipoint] -1 << " " 
						 << residuals3sigmaplus[ipoint] -1 << " " << std::endl;
			}
			continue; 
		}
		else if( (residualssigmaminus[ipoint]==0) && (TMath::Abs(residuals3sigmaplus[ipoint]-1.)>limitpoint) ) { // JLK : kill crazy upperlimit
			if(!quiet) {
				WARN_OUT << "This point is killed! (upper limit with useless stat (JLK : TBC))" << std::endl;
				INFO_OUT << "e_mean" << " " << "residual_exp" << " " << "residual_sigma_minus" << " " 
						 << "residual_sigma_plus" << std::endl;
				INFO_OUT << "1sigma : " << energybin[ipoint] << " " << residuals[ipoint] -1 << " " << residualssigmaminus[ipoint] -1 << " " 
						 << residualssigmaplus[ipoint] -1 << " " << std::endl;
				INFO_OUT << "3sigma : " << energybin[ipoint] << " " << residuals[ipoint] -1 << " " << residuals3sigmaminus[ipoint] -1 << " " 
						 << residuals3sigmaplus[ipoint] -1 << " " << std::endl;
			}
			continue; 
		}
		else if((fresidualssigmaplus[ipoint]-1.)>limitpoint || (1.-fresidualssigmaplus[ipoint])>limitpoint) { // JLK : kill crazy point
			if(!quiet) {
				WARN_OUT << "This point is killed! (point with useless stat (JLK : TBC))" << std::endl;
				INFO_OUT << "e_mean" << " " << "residual_exp" << " " << "residual_sigma_minus" << " " 
						 << "residual_sigma_plus" << std::endl;
				INFO_OUT << "1sigma : " << energybin[ipoint] << " " << residuals[ipoint] -1 << " " << residualssigmaminus[ipoint] -1 << " " 
						 << residualssigmaplus[ipoint] -1 << " " << std::endl;
				INFO_OUT << "3sigma : " << energybin[ipoint] << " " << residuals[ipoint] -1 << " " << residuals3sigmaminus[ipoint] -1 << " " 
						 << residuals3sigmaplus[ipoint] -1 << " " << std::endl;
			}
			continue; 
		}
		else {
			energybin_tmp.push_back(energybin[ipoint]);
			residuals_tmp.push_back(residuals[ipoint]);
			residualssigmaplus_tmp.push_back(residualssigmaplus[ipoint]);
			residualssigmaminus_tmp.push_back(residualssigmaminus[ipoint]);
			residuals3sigmaplus_tmp.push_back(residuals3sigmaplus[ipoint]);
			residuals3sigmaminus_tmp.push_back(residuals3sigmaminus[ipoint]);
			flux_tmp.push_back(flux[ipoint]);
			fluxsigmaplus_tmp.push_back(fluxsigmaplus[ipoint]);
			fluxsigmaminus_tmp.push_back(fluxsigmaminus[ipoint]);
			flux3sigmaplus_tmp.push_back(flux3sigmaplus[ipoint]);
			flux3sigmaminus_tmp.push_back(flux3sigmaminus[ipoint]);
		}
	}
 
	energybin = energybin_tmp;
	residuals = residuals_tmp; 
	residualssigmaplus = residualssigmaplus_tmp; 
	residualssigmaminus = residualssigmaminus_tmp; 
	residuals3sigmaplus = residuals3sigmaplus_tmp; 
	residuals3sigmaminus = residuals3sigmaminus_tmp; 
	flux = flux_tmp;
	fluxsigmaplus = fluxsigmaplus_tmp; 
	fluxsigmaminus = fluxsigmaminus_tmp; 
	flux3sigmaplus = flux3sigmaplus_tmp; 
	flux3sigmaminus = flux3sigmaminus_tmp;

	//DEBUG_OUT << "energybin.size()=" << energybin.size() << std::endl;

	DEBUG_OUT << "START::PlotFactory::BadPointsKiller : job done." << std::endl;
}



/**
 * \brief Determine the TH2's residuals height
 *
 * By default the height is 1 but if there if the sigma + is greater than 1
 * we build a greater TH2
 */
double START::PlotFactory::GetPlotResidualsMaximum(const std::vector<double> &residualssigmaminus, 
												   const std::vector<double> &residualssigmaplus)
{
  
	double limitresiduals(1.);

	bool isupperlimit(false);
	std::vector<double> resplus,resminus;
	for(unsigned int ipoint(0); ipoint<residualssigmaplus.size(); ipoint++) {
		if(residualssigmaminus[ipoint]!=0) {
			resplus.push_back(residualssigmaplus[ipoint]);
			resminus.push_back(residualssigmaminus[ipoint]);
		}
		if(isupperlimit) continue;
		if(residualssigmaminus[ipoint]==0) isupperlimit=true;
	}

	std::vector<double>::const_iterator it1,it2;

	isupperlimit=false; // JLK : no upper limits on residuals. If one day it changes erase this line

	if(!isupperlimit) {
		it1 = max_element(resplus.begin(),resplus.end());
		it2 = max_element(resminus.begin(),resminus.end());
	}
	if(isupperlimit) {
		it1 = max_element(residualssigmaplus.begin(),residualssigmaplus.end()); // JLK no upper limit on residuals
		it2 = max_element(residualssigmaplus.begin(),residualssigmaplus.end()); // JLK no upper limit on residuals
	}
  
	for(unsigned int ires(0); ires<resplus.size(); ires++) {
		DEBUG_OUT << "res-=" << resminus[ires] -1 << " res+=" << resplus[ires]-1 << std::endl;
	}

	double limit1 = TMath::Abs(*it1-1.);
	double limit2 = TMath::Abs(*it2-1.);

	DEBUG_OUT << "limit1 = " << limit1 << " limit2=" << limit2 << std::endl;

	if(limit1>1. && limit1>=limit2) limitresiduals=limit1+0.1;
	else if(limit2>1. && limit2>=limit2 && !isupperlimit) limitresiduals=limit2+0.1;
	else limitresiduals = 1.;

	DEBUG_OUT << " START::PlotFactory::GetPlotResidualsMaximum : limit = " << limitresiduals << std::endl;

	return limitresiduals;
  
}

/**
 * \brief Get energy range for spectrum and residuals (for the TH2s)
 */
std::pair<double,double> START::PlotFactory::GetPlotEnergyRange(const Hypothesis &hypothesis)
{
  
	const Hypothesis *hypo = &hypothesis;

	std::vector<double> binmin, binmax;

	std::vector<Band> BandArray = hypo->GetBandArray(); 

	for(std::vector<Band>::const_iterator band = BandArray.begin();band!=BandArray.end();++band) {
    
		if(band->GetKeepBand()==0) continue; // We are interested only in the selected bands

		for (std::vector<EnergyBin>::const_iterator bin = band->ebin.begin();bin!=band->ebin.end();++bin) {

			if (bin->GetKeepBin()==0) continue; // skip energy bins below threshold

			binmin.push_back(bin->GetEmin());
			binmax.push_back(bin->GetEmax());

		}
	}

	std::pair<double,double> energyrange;

	double shift(0.4);

	std::vector<double>::const_iterator it1, it2;
	it1 = min_element(binmin.begin(),binmin.end());
	it2 = max_element(binmax.begin(),binmax.end());

	double emin(*it1), emax(*it2);
	emin=emin*TMath::Exp(-shift);
	emax=emax*TMath::Exp(+shift);

	energyrange = std::make_pair(emin,emax);

	DEBUG_OUT << "PlotEnergyRange : emin = " << emin << "emax = " << emax << std::endl;

	return energyrange;
}


/**
 * \brief Get energy range for the butterfly (fit range)
 */
std::pair<double,double> START::PlotFactory::GetButterflyEnergyRange(const Hypothesis &hypothesis)
{

	const Hypothesis *hypo = &hypothesis;

	std::vector<double> binmin, binmax;
	std::vector<Band> BandArray = hypo->GetBandArray(); // JLK : Je passe par une affectation sinon la boucle sur les bins est infinie. prob de pointeurs qui me depasse

	//for(std::vector<Band>::const_iterator band = hypo->GetBandArray().begin();band!=hypo->GetBandArray().end();++band) {
	for(std::vector<Band>::const_iterator band = BandArray.begin();band!=BandArray.end();++band) {
    
		if(band->GetKeepBand()==0) continue; // We are interested only in the selected bands

		for (std::vector<EnergyBin>::const_iterator bin = band->ebin.begin();bin!=band->ebin.end();++bin) {

			if (bin->GetKeepBin() == 0) continue; // skip energy bins below threshold

			binmin.push_back(bin->GetEmin());
			binmax.push_back(bin->GetEmax());

		}
	}

	std::pair<double,double> energyrange;

	std::vector<double>::const_iterator it1, it2;
	it1 = min_element(binmin.begin(),binmin.end());
	it2 = max_element(binmax.begin(),binmax.end());

	double emin(*it1), emax(*it2);

	energyrange = std::make_pair(emin,emax);

	DEBUG_OUT << "ButterflyEnergyRange : emin = " << emin << " emax = " << emax << std::endl;
  
	return energyrange;
}

/**
 * \brief Get flux range for plot
 */
std::pair<double,double> START::PlotFactory::GetPlotFluxRange(Hypothesis &hypothesis,
															  const std::vector<double> &energybin,
															  const std::vector<double> &flux, 
															  const std::vector<double> &fluxsigmaminus,
															  const std::vector<double> &fluxsigmaplus,
															  const std::vector<double> &flux3sigmaplus)
{

	DEBUG_OUT_L(2) << "Flux points :" << std::endl;

	for(unsigned int ipoint(0); ipoint<flux.size(); ipoint++) {
		DEBUG_OUT_L(2) << ipoint << " flux = " << flux[ipoint] 
					   << " fluxsigmaminus = " << fluxsigmaminus[ipoint]
					   << " fluxsigmaplus = " << fluxsigmaplus[ipoint] << std::endl;
	}

	// We first look at flux points

	bool isupperlimit(false);

	std::vector<unsigned int> indexupperlimit;

	// get upperlimits index
	for(unsigned int ipoint(0); ipoint<flux.size(); ipoint++) {
		if(fluxsigmaminus[ipoint]==0 || flux[ipoint]<=0.) indexupperlimit.push_back(ipoint);
	}

	if(indexupperlimit.size()>0) {
		isupperlimit=true;
		DEBUG_OUT << "Upper limit : " << std::endl;
		for(unsigned int ip(0); ip<indexupperlimit.size(); ip++) 
			DEBUG_OUT << "energy=" << energybin[indexupperlimit[ip]] << " flux=" << flux[indexupperlimit[ip]] 
					  << " fluxsigmaminus=" << fluxsigmaminus[indexupperlimit[ip]]
					  << " flux3sigmaplus=" << flux3sigmaplus[indexupperlimit[ip]] 
					  << " flux3sigmaplusbottomarrow=" << flux3sigmaplus[indexupperlimit[ip]]*fscalearrowfactor << std::endl;
	}

	std::vector<double> vfluxmin;
	std::vector<double> vfluxmax;

	std::vector<double> vfluxminupper;
	std::vector<double> vfluxmaxupper;

	// no upperlimit
	for(unsigned int i(0); i<flux.size(); i++) {
		if(flux[i]>0. && fluxsigmaminus[i]!=0) {
			vfluxmax.push_back(fluxsigmaplus[i]);
			vfluxmin.push_back(fluxsigmaminus[i]);
		}
	}

	DEBUG_OUT <<"after no upper limits, vfluxmin.size() = " << vfluxmin.size() << std::endl;

	// we add upper limits
	if(isupperlimit) {
		for(unsigned int upper(0); upper<indexupperlimit.size(); upper++) {
			vfluxmax.push_back(flux3sigmaplus[indexupperlimit[upper]]);
			vfluxmin.push_back(flux3sigmaplus[indexupperlimit[upper]]*fscalearrowfactor);
		}
	}

	DEBUG_OUT <<"after upper limits, vfluxmin.size() = " << vfluxmin.size() << std::endl;

	DEBUG_OUT_L(2) << "Points from flux : " << std::endl;

	for(unsigned int ipoint(0); ipoint<vfluxmax.size(); ipoint++) {
		DEBUG_OUT_L(2) << ipoint << " fluxmin = " << vfluxmin[ipoint] << " fluxmax = " << vfluxmax[ipoint] << std::endl;
	}

	double fluxminpoints;
	double fluxmaxpoints;

	std::vector<double>::const_iterator itpoints1, itpoints2;
	itpoints1 = min_element(vfluxmin.begin(), vfluxmin.end());
	itpoints2 = max_element(vfluxmax.begin(), vfluxmax.end());
	fluxminpoints = *itpoints1;
	fluxmaxpoints = *itpoints2;

	DEBUG_OUT_L(2) << "found fluxmin = " << fluxminpoints <<  " fluxmax = " << fluxmaxpoints << std::endl;

	//We look at butterfly

	double fluxminbutt(999);
	double fluxmaxbutt(999);
	if(fButterfly==LinearButterfly) { // if butterfly is negative we set flux limits with LogarithmButterfly

		std::string polylinevectorsname = MakeKeyPolyLineVectorsSpectrum(hypothesis.GetName());
		std::pair<std::vector<double>,std::vector<double> > butterflyvectors;
		butterflyvectors = fMapButterFlyVectors[polylinevectorsname];
		std::vector<double> x = butterflyvectors.first;
		std::vector<double> y = butterflyvectors.second;
    
		// we look for negative value
		bool negativevalue = false;
		for(unsigned int ipoint(0); ipoint<y.size(); ipoint++) {
			if(y[ipoint]<0) {
				negativevalue=true;
				break;
			}
		}

		DEBUG_OUT_L(2) << "Points from butterfly : " << std::endl;
    
		for(unsigned int ipoint(0); ipoint<x.size(); ipoint++) {
			DEBUG_OUT_L(2) << ipoint << " x = " << x[ipoint] << " y = " << y[ipoint] << std::endl;
		}

		if(negativevalue) { // we take limits from LogarithmButterfly
			// get energy range
			std::pair<double,double> erange = GetButterflyEnergyRange(hypothesis);
			double emin(0.), emax(0.);
			emin = erange.first;
			emax = erange.second;
      
			double log10emin = TMath::Log10(emin);
			double log10emax = TMath::Log10(emax);
      
			unsigned int nbpoints(100);

			x.clear(); y.clear();
			ComputeLogarithmButterfly(hypothesis,x,y,nbpoints,log10emin,log10emax);
		}

		if(x.size()>0 && y.size()>0) {
			std::vector<double>::const_iterator itbutt1, itbutt2;
			itbutt1 = min_element(y.begin(), y.end());
			itbutt2 = max_element(y.begin(), y.end());
			fluxminbutt = *itbutt1;
			fluxmaxbutt = *itbutt2;
		}

		DEBUG_OUT_L(2) << "found fluxminbutt = " << fluxminbutt <<  " fluxmaxbutt = " << fluxmaxbutt << std::endl;
	}
	if(fButterfly==LogarithmButterfly || fButterfly==ContoursButterfly) {

		std::string polylinevectorsname = MakeKeyPolyLineVectorsSpectrum(hypothesis.GetName());
		std::pair<std::vector<double>,std::vector<double> > butterflyvectors;
		butterflyvectors = fMapButterFlyVectors[polylinevectorsname];
		std::vector<double> x = butterflyvectors.first;
		std::vector<double> y = butterflyvectors.second;
    
		DEBUG_OUT_L(2) << "Points from butterfly : " << std::endl;
    
		for(unsigned int ipoint(0); ipoint<x.size(); ipoint++) {
			DEBUG_OUT_L(2) << ipoint << " x = " << x[ipoint] << " y = " << y[ipoint] << std::endl;
		}

		if(x.size()>0 && y.size()>0) {
			std::vector<double>::const_iterator itbutt1, itbutt2;
			itbutt1 = min_element(y.begin(), y.end());
			itbutt2 = max_element(y.begin(), y.end());
			fluxminbutt = *itbutt1;
			fluxmaxbutt = *itbutt2;
		}

		DEBUG_OUT_L(2) << "found fluxminbutt = " << fluxminbutt <<  " fluxmaxbutt = " << fluxmaxbutt << std::endl;
	}
	else if (fButterfly==CausticButterfly){

		std::string grapharrayname = MakeKeyGraphArraySpectrum(hypothesis.GetName());
		std::vector<double> y; //will contain the flux from butterfly graph contours
		y.clear();

		DEBUG_OUT_L(2) << "Points from caustic butterfly (contours) : " << std::endl;

		for(unsigned int igraph(0); igraph<fMapGraphArraySpectrum[grapharrayname].size();igraph++) {
      
			DEBUG_OUT << "igraph = " << igraph << " nb points =" << fMapGraphArraySpectrum[grapharrayname][igraph]->GetN() << std::endl;

			for(int ipoint(0); ipoint<fMapGraphArraySpectrum[grapharrayname][igraph]->GetN(); ipoint++) {
				double a(0),b(0);
				fMapGraphArraySpectrum[grapharrayname][igraph]->GetPoint(ipoint,a,b);
				y.push_back(b);
				DEBUG_OUT_L(2) << ipoint << " x = " << a << " y = " << b << std::endl;
			}

		}

		if(y.size()>0) {
			std::vector<double>::const_iterator itbutt1, itbutt2;
			itbutt1 = min_element(y.begin(), y.end());
			itbutt2 = max_element(y.begin(), y.end());
			fluxminbutt = *itbutt1;
			fluxmaxbutt = *itbutt2;
		}

	}

	DEBUG_OUT << "fluxminpoint = " << fluxminpoints << " fluxmaxpoint = " << fluxmaxpoints << std::endl;
	DEBUG_OUT << "fluxminbutt = " << fluxminbutt << " fluxmaxbutt = " << fluxmaxbutt << std::endl;

	std::pair<double,double> fluxrange;
	double fluxmin, fluxmax;
  
	if(fluxminbutt!=999 && fluxmaxbutt!=999) {
		if(fluxminpoints>fluxminbutt) fluxmin = fluxminbutt;
		else fluxmin = fluxminpoints;
    
		if(fluxmaxpoints>fluxmaxbutt) fluxmax = fluxmaxpoints;
		else fluxmax = fluxmaxbutt;
	}
	else {
		fluxmin = fluxminpoints;
		fluxmax = fluxmaxpoints;
	}

	// shift
	double shift(1.);
	fluxmin*=TMath::Power(10.,-shift);
	fluxmax*=TMath::Power(10.,+shift);
  
	fluxrange = std::make_pair(fluxmin,fluxmax);

	if(isupperlimit) {
		DEBUG_OUT_L(2) << "Upperlimit!! " << std::endl;
		for(unsigned int upper(0); upper<indexupperlimit.size(); upper++) {
			DEBUG_OUT_L(2) << "energy = " << energybin[upper] << std::endl;
		}
	}

	DEBUG_OUT << " START::PlotFactory::GetPlotFluxRange : fluxmin = " << fluxmin << " fluxmax = " << fluxmax << std::endl;

	return fluxrange;
}


/**
 * \brief Make and draw contours for all pair of fitted parameters 
 *
 * All the ROOT's objects are initialized here and graph with contours
 * at 1, 2 and 3 sigmas are Filled here.
 */
void START::PlotFactory::PlotContours()
{

	CleanCanvasContours();
	CleanMultiGraphs();

	for(std::vector<Hypothesis*>::iterator hypo=fHypothesisArray.begin(); hypo!=fHypothesisArray.end(); hypo++) {
    
		if((*hypo)->GetConvergence()) {

			std::vector<std::string> paramname = (*hypo)->GetParametersNames();
			std::vector<std::string> paramunits = (*hypo)->GetParametersUnits();

			std::vector <std::pair <std::vector <std::pair <double,double> >, std::pair<int,int> > > vectorcontour1 = (*hypo)->GetContourSigma1();
			std::vector <std::pair <std::vector <std::pair <double,double> >, std::pair<int,int> > > vectorcontour2 = (*hypo)->GetContourSigma2();
			std::vector <std::pair <std::vector <std::pair <double,double> >, std::pair<int,int> > > vectorcontour3 = (*hypo)->GetContourSigma3();

			bool skipvec1 = true;
			bool skipvec2 = true;
			bool skipvec3 = true;

			if(vectorcontour1.size()>0) skipvec1=false;
			if(vectorcontour2.size()>0) skipvec2=false;
			if(vectorcontour2.size()>0) skipvec3=false;

			if(skipvec1 && skipvec2 && skipvec3) WARN_OUT << "Contours are not computed so they can't be plot!" << std::endl;

			int nbcanvas(0); // number of canvas

			if(!skipvec1) {
				nbcanvas = vectorcontour1.size();
			}
			else if(!skipvec2) {
				nbcanvas = vectorcontour2.size();
			} 
			else if(!skipvec3) {
				nbcanvas = vectorcontour3.size();
			} 

			for(int icont(0); icont<nbcanvas; icont++) {
  
				std::vector<std::pair<double,double> > pointscontour1;
				if(!skipvec1) pointscontour1 = vectorcontour1[icont].first; // get contours points 
				std::vector<std::pair<double,double> > pointscontour2;
				if(!skipvec2) pointscontour2 = vectorcontour2[icont].first; // get contours points 
				std::vector<std::pair<double,double> > pointscontour3;
				if(!skipvec3) pointscontour3 = vectorcontour3[icont].first; // get contours points 
	
				bool skippoints1 = true;
				bool skippoints2 = true;
				bool skippoints3 = true;

				if(!skipvec1 && pointscontour1.size()>3) skippoints1=false;
				if(!skipvec2 && pointscontour2.size()>3) skippoints2=false;
				if(!skipvec3 && pointscontour3.size()>3) skippoints3=false;

				if(skippoints1) INFO << "No contours at 1 sigma" << std::endl;
				if(skippoints2) INFO << "No contours at 2 sigma" << std::endl;
				if(skippoints3) INFO << "No contours at 3 sigma" << std::endl;
	
				std::pair<int,int> param1_12;
				if(!skippoints1) param1_12 = vectorcontour1[icont].second; // get the indices of the two parameters contours
				std::pair<int,int> param2_12;
				if(!skippoints2) param2_12 = vectorcontour2[icont].second; // get the indices of the two parameters contours
				std::pair<int,int> param3_12;
				if(!skippoints3) param3_12 = vectorcontour3[icont].second; // get the indices of the two parameters contours
	
				/* We build the canvas */

				std::string canvasname="";
				if(!skippoints1) {
					canvasname+=MakeKeyCanvasContours((*hypo)->GetName(),paramname[param1_12.first],paramname[param1_12.second]);
				}
				else if(!skippoints2) {
					canvasname+=MakeKeyCanvasContours((*hypo)->GetName(),paramname[param2_12.first],paramname[param2_12.second]);
				}
				else if(!skippoints2) {
					canvasname+=MakeKeyCanvasContours((*hypo)->GetName(),paramname[param3_12.first],paramname[param3_12.second]);
				}
				else continue;

				fMapCanvasContours[canvasname] = new TCanvas(canvasname.c_str(),canvasname.c_str());
				std::map<std::string,TCanvas*>::iterator canvas = fMapCanvasContours.find(canvasname);

				// We build the names of the TGraphs 

				TGraph *grcontour1sigma = new TGraph();
				TGraph *grcontour2sigma = new TGraph();
				TGraph *grcontour3sigma = new TGraph();

				TString grname1 = "Contours1sigma";
				grname1+=(*hypo)->GetName();
				grname1+="_";
				if(!skippoints1) {
					grname1+=paramname[param1_12.first];
					grname1+="_";
					grname1+=paramname[param1_12.second];
				}
				grcontour1sigma->SetName(grname1);
				TString grname2 = "Contours2sigma";
				grname2+=(*hypo)->GetName();
				grname2+="_";
				if(!skippoints2) {
					grname2+=paramname[param2_12.first];
					grname2+="_";
					grname2+=paramname[param2_12.second];
				}
				grcontour2sigma->SetName(grname2);

				TString grname3 = "Contours3sigma";
				grname3+=(*hypo)->GetName();
				grname3+="_";
				if(!skippoints1) {
					grname3+=paramname[param3_12.first];
					grname3+="_";
					grname3+=paramname[param3_12.second];
				}
				grcontour3sigma->SetName(grname3);

				DEBUG_OUT << "graph 1 sigma's name : " << grname1 << std::endl;
				DEBUG_OUT << "graph 2 sigma's name : " << grname2 << std::endl;
				DEBUG_OUT << "graph 3 sigma's name : " << grname3 << std::endl;


				// Multigraph 

				std::string multigrname = "";
				if(!skippoints1) {
					multigrname+=MakeKeyMultiGraphContours((*hypo)->GetName(),paramname[param1_12.first],paramname[param1_12.second]);
				}
				else if(!skippoints2) {
					multigrname+=MakeKeyMultiGraphContours((*hypo)->GetName(),paramname[param2_12.first],paramname[param2_12.second]);
				}
				else if(!skippoints2) {
					multigrname+=MakeKeyMultiGraphContours((*hypo)->GetName(),paramname[param3_12.first],paramname[param3_12.second]);
				}
				else continue;
				fMapMultiGraphContours[multigrname] = new TMultiGraph();
				std::map<std::string,TMultiGraph*>::iterator multigr = fMapMultiGraphContours.find(multigrname);

				multigr->second->SetName(multigrname.c_str());

				DEBUG_OUT << "multigraph's name : " << multigrname << std::endl;

				// We build the name of the axis 

				std::string xtitle;
				std::string ytitle;

				std::vector<std::string> latexparamname = (*hypo)->GetLatexParametersNames();

				if(!skippoints1) {
					xtitle = latexparamname[param1_12.first];
					if(paramunits[param1_12.first].size()> 0) xtitle+=" (";
					xtitle+=paramunits[param1_12.first];
					if(paramunits[param1_12.first].size()> 0) xtitle+=")";
					ytitle = latexparamname[param1_12.second];
					if(paramunits[param1_12.second].size()> 0) ytitle+=" (";
					ytitle+=paramunits[param1_12.second];
					if(paramunits[param1_12.second].size()> 0) ytitle+=")";
				}
				else if(!skippoints2) {
					xtitle = latexparamname[param2_12.first];
					if(paramunits[param2_12.first].size()> 0) xtitle+=" (";
					xtitle+=paramunits[param2_12.first];
					if(paramunits[param2_12.first].size()> 0) xtitle+=")";
					ytitle = latexparamname[param2_12.second];
					if(paramunits[param2_12.second].size()> 0) ytitle+=" (";
					ytitle+=paramunits[param2_12.second];
					if(paramunits[param2_12.second].size()> 0) ytitle+=")";
				}
				else if(!skippoints3) {
					xtitle = latexparamname[param3_12.first];
					if(paramunits[param3_12.first].size()> 0) xtitle+=" (";
					xtitle+=paramunits[param3_12.first];
					if(paramunits[param3_12.first].size()> 0) xtitle+=")";
					ytitle = latexparamname[param3_12.second];
					if(paramunits[param3_12.second].size()> 0) ytitle+=" (";
					ytitle+=paramunits[param3_12.second];
					if(paramunits[param3_12.second].size()> 0) ytitle+=")";
				}
				else continue;

				DEBUG_OUT << "xtitle = " << xtitle << std::endl;
				DEBUG_OUT << "ytitle = " << ytitle << std::endl;

				// We build the titles of the TGraphs 

				TString grtitle = "Confidence intervals ";
				TString multigrtitle = "Confidence intervals "; 

				if(!skippoints1) {
					grtitle+=paramname[param1_12.first];
					grtitle+="/";
					grtitle+=paramname[param1_12.second];
					multigrtitle+=paramname[param1_12.first];
					multigrtitle+="/";
					multigrtitle+=paramname[param1_12.second];
				}
				else if(!skippoints2) {
					grtitle+=paramname[param2_12.first];
					grtitle+="/";
					grtitle+=paramname[param2_12.second];
					multigrtitle+=paramname[param2_12.first];
					multigrtitle+="/";
					multigrtitle+=paramname[param2_12.second];
				}
				else if(!skippoints2) {
					grtitle+=paramname[param3_12.first];
					grtitle+="/";
					grtitle+=paramname[param3_12.second];
					multigrtitle+=paramname[param3_12.first];
					multigrtitle+="/";
					multigrtitle+=paramname[param3_12.second];
				}
				else continue;

				DEBUG_OUT << "graph 1 sigma's title : " << grtitle << std::endl;
				DEBUG_OUT << "graph 2 sigma's title : " << grtitle << std::endl;
				DEBUG_OUT << "graph 3 sigma's title : " << grtitle << std::endl;

				grcontour1sigma->SetTitle(grtitle);
				grcontour2sigma->SetTitle(grtitle);
				grcontour3sigma->SetTitle(grtitle);

				multigr->second->SetTitle(multigrtitle);

				/* Fill Tgraph with contours*/
				DEBUG_OUT_L(2) << "Contours points sigma 1" << std::endl;
				if(!skippoints1) {
					for(unsigned int ipoint(0); ipoint<pointscontour1.size(); ipoint++) {
						double x1 = pointscontour1[ipoint].first;
						double y1 = pointscontour1[ipoint].second;
						DEBUG_OUT_L(2) << "x1 = " << x1 << " | y1 = " << y1 << std::endl;
						grcontour1sigma->SetPoint(ipoint,x1,y1);
					}
				}
				DEBUG_OUT_L(2) << "Contours points sigma 2" << std::endl;
				if(!skippoints2) {
					for(unsigned int ipoint(0); ipoint<pointscontour2.size(); ipoint++) {
						double x2 = pointscontour2[ipoint].first;
						double y2 = pointscontour2[ipoint].second;
						DEBUG_OUT_L(2) << "x2 = " << x2 << " | y2 = " << y2 << std::endl;
						grcontour2sigma->SetPoint(ipoint,x2,y2);
					} 
				}
				DEBUG_OUT_L(2) << "Contours points sigma 3" << std::endl;
				if(!skippoints3) {
					for(unsigned int ipoint(0); ipoint<pointscontour3.size(); ipoint++) {
						double x3 = pointscontour3[ipoint].first;
						double y3 = pointscontour3[ipoint].second;
						DEBUG_OUT_L(2) << "x3 = " << x3 << " | y3 = " << y3 << std::endl;
						grcontour3sigma->SetPoint(ipoint,x3,y3);
					}
				}
				// We close the contours
				if(!skippoints1) grcontour1sigma->SetPoint(pointscontour1.size(),pointscontour1[0].first,pointscontour1[0].second);
				if(!skippoints2) grcontour2sigma->SetPoint(pointscontour2.size(),pointscontour2[0].first,pointscontour2[0].second);
				if(!skippoints3) grcontour3sigma->SetPoint(pointscontour3.size(),pointscontour3[0].first,pointscontour3[0].second);
	
				/* Set line color */
	
				if(!skippoints1) grcontour1sigma->SetLineColor(kRed);
				if(!skippoints2) grcontour2sigma->SetLineColor(kBlue);
				if(!skippoints3) grcontour3sigma->SetLineColor(kGreen);

				TLegend *legmulti = new TLegend(0.6,0.8,0.89,0.89);
				legmulti->SetLineColor(0);
				//legmulti->SetLineStyle(0);
				//legmulti->SetFillColor(0);
				legmulti->SetFillStyle(0);

				if(!skippoints1) legmulti->AddEntry(grcontour1sigma,"68.27% confidence level","l");
				if(!skippoints2) legmulti->AddEntry(grcontour2sigma,"95.45% confidence level","l");
				if(!skippoints3) legmulti->AddEntry(grcontour3sigma,"99.73% confidence level","l");

				if(!skippoints1) multigr->second->Add(grcontour1sigma);
				if(!skippoints2) multigr->second->Add(grcontour2sigma);
				if(!skippoints3) multigr->second->Add(grcontour3sigma);

				TMarker *Marker = 0;
				if(!skippoints1 || !skippoints2 || !skippoints3) {
					unsigned int param1(0),param2(0);
					if(!skippoints1) {
						param1=param1_12.first;
						param2=param1_12.second;
					}
					else if(!skippoints2) {
						param1=param2_12.first;
						param2=param2_12.second;
					}
					else if(!skippoints3) {
						param1=param3_12.first;
						param2=param3_12.second;
					}
					Marker = new TMarker();
					Marker->SetX((*hypo)->GetFittedParameters()[param1]);
					Marker->SetY((*hypo)->GetFittedParameters()[param2]);
					Marker->SetMarkerStyle(5);
					Marker->SetMarkerColor(kBlue);
					Marker->SetMarkerSize(1.5);
				}

				canvas->second->cd();

				multigr->second->Draw("AC");
				if(Marker!=0) Marker->Draw();

				multigr->second->GetXaxis()->SetTitle(xtitle.c_str());
				multigr->second->GetXaxis()->SetLabelSize(0.03);
				multigr->second->GetXaxis()->CenterTitle();
				multigr->second->GetYaxis()->SetTitle(ytitle.c_str());
				multigr->second->GetYaxis()->SetLabelSize(0.03);
				multigr->second->GetYaxis()->CenterTitle();

				canvas->second->cd();

				legmulti->Draw("same");

				canvas->second->Update();
				/*
				  delete grcontour1sigma; grcontour1sigma=0;
				  delete grcontour2sigma; grcontour3sigma=0;
				  delete grcontour3sigma; grcontour3sigma=0;
				  delete legmulti; legmulti=0;
				*/
			} // end loop canvas

		}
		else {
			WARNING << "Skip hypothesis " << (*hypo)->GetName() << " because it didn't converge!" << std::endl;
			continue;
		}
    
	} // end loop hypothesis

	INFO << "Plotting contours... ok" << std::endl;
}


/**
 * \brief Clean canvas's contours
 */
void START::PlotFactory::CleanCanvasContours()
{
  
	for(std::map<std::string,TCanvas*>::iterator canvas=fMapCanvasContours.begin(); canvas!=fMapCanvasContours.end(); canvas++) {
		if(canvas->second!=0) delete canvas->second;
		canvas->second=0;
	}
	fMapCanvasContours.clear();
}

/**
 * \brief Clean TGraphAsymmErrors
 */
void START::PlotFactory::CleanMultiGraphs()
{

	for(std::map<std::string,TMultiGraph*>::iterator graph=fMapMultiGraphContours.begin(); graph!=fMapMultiGraphContours.end(); graph++) {
		if(graph->second!=0) delete graph->second;
		graph->second=0;
	}
	fMapMultiGraphContours.clear();
}

/**
 * \brief Canvas spectrum key
 */
std::string START::PlotFactory::MakeKeyCanvasSpectrum(std::string hypothesisname)
{
	std::string key;
	key = "TCanvas_SpectrumAndResiduals";
	key+="_";
	switch(fButterfly) {
	case NoButterfly:
		key+="NoButterfly";
		break;
	case LinearButterfly:
		key+="LinearButterfly";
		break;
	case LogarithmButterfly:
		key+="LogarithmButterfly";
		break;
	case ContoursButterfly:
		key+="ContoursButterfly";
		break;
	case CausticButterfly:
		key+="CausticButterfly";
		break;
	}
	key+="_";
	key+=hypothesisname;

	return key;
}


/**
 * \brief Canvas spectrum key
 */
std::string START::PlotFactory::MakeKeyCanvasContours(std::string hypothesisname, std::string param1, std::string param2)
{
	std::string key;
	key = "TCanvas_Contours";
	key+="_";
	key+=hypothesisname;
	key+="_";
	key+=param1;
	key+="_";
	key+=param2;

	return key;
}

/**
 * \brief Pad's spectrum key
 */
std::string START::PlotFactory::MakeKeyPadSpectrum(std::string hypothesisname)
{
	std::string key;
	key = "TPad_Spectrum";
	key+="_";
	switch(fButterfly) {
	case NoButterfly:
		key+="NoButterfly";
		break;
	case LinearButterfly:
		key+="LinearButterfly";
		break;
	case LogarithmButterfly:
		key+="LogarithmButterfly";
		break;
	case ContoursButterfly:
		key+="ContoursButterfly";
		break;
	case CausticButterfly:
		key+="CausticButterfly";
		break;
	}
	key+="_";
	key+=hypothesisname;

	return key;

}

/**
 * \brief Pad's residuals key
 */
std::string START::PlotFactory::MakeKeyPadResiduals(std::string hypothesisname)
{

	std::string key;
	key = "TPad_Residuals";
	key+="_";
	switch(fButterfly) {
	case NoButterfly:
		key+="NoButterfly";
		break;
	case LinearButterfly:
		key+="LinearButterfly";
		break;
	case LogarithmButterfly:
		key+="LogarithmButterfly";
		break;
	case ContoursButterfly:
		key+="ContoursButterfly";
		break;
	case CausticButterfly:
		key+="CausticButterfly";
		break;
	}
	key+="_";
	key+=hypothesisname;

	return key;

}

/**
 * \brief TH2's spectrum key
 */
std::string START::PlotFactory::MakeKeyTH2Spectrum(std::string hypothesisname)
{

	std::string key;
	key = "TH2_Spectrum";
	key+="_";
	switch(fButterfly) {
	case NoButterfly:
		key+="NoButterfly";
		break;
	case LinearButterfly:
		key+="LinearButterfly";
		break;
	case LogarithmButterfly:
		key+="LogarithmButterfly";
		break;
	case ContoursButterfly:
		key+="ContoursButterfly";
		break;
	case CausticButterfly:
		key+="CausticButterfly";
		break;
	}
	key+="_";
	key+=hypothesisname;

	return key;

}

/**
 * \brief TH2's residuals key
 */
std::string START::PlotFactory::MakeKeyTH2Residuals(std::string hypothesisname)
{

	std::string key;
	key = "TH2_Residuals";
	key+="_";
	switch(fButterfly) {
	case NoButterfly:
		key+="NoButterfly";
		break;
	case LinearButterfly:
		key+="LinearButterfly";
		break;
	case LogarithmButterfly:
		key+="LogarithmButterfly";
		break;
	case ContoursButterfly:
		key+="ContoursButterfly";
		break;
	case CausticButterfly:
		key+="CausticButterfly";
		break;
	}
	key+="_";
	key+=hypothesisname;

	return key;

}

/**
 * \brief TGraphAsymmError's spectrum key
 */
std::string START::PlotFactory::MakeKeyGraphErrorsSpectrum(std::string hypothesisname)
{

	std::string key;
	key = "TGraphAssymErrors_Spectrum";
	key+="_";
	switch(fButterfly) {
	case NoButterfly:
		key+="NoButterfly";
		break;
	case LinearButterfly:
		key+="LinearButterfly";
		break;
	case LogarithmButterfly:
		key+="LogarithmButterfly";
		break;
	case ContoursButterfly:
		key+="ContoursButterfly";
		break;
	case CausticButterfly:
		key+="CausticButterfly";
		break;
	}
	key+="_";
	key+=hypothesisname;

	return key;

}

/**
 * \brief TGraphAsymmError's residuals key
 */
std::string START::PlotFactory::MakeKeyGraphErrorsResiduals(std::string hypothesisname)
{

	std::string key;
	key = "TGraphAssymErrors_Residuals";
	key+="_";
	switch(fButterfly) {
	case NoButterfly:
		key+="NoButterfly";
		break;
	case LinearButterfly:
		key+="LinearButterfly";
		break;
	case LogarithmButterfly:
		key+="LogarithmButterfly";
		break;
	case ContoursButterfly:
		key+="ContoursButterfly";
		break;
	case CausticButterfly:
		key+="CausticButterfly";
		break;
	}
	key+="_";
	key+=hypothesisname;

	return key;

}

/**
 * \brief TGraphAsymmError's residuals key
 */
std::string START::PlotFactory::MakeKeyGraphErrorsResidualsOn(std::string hypothesisname)
{

	std::string key;
	key = "TGraphAssymErrors_ResidualsOn";
	key+="_";
	switch(fButterfly) {
	case NoButterfly:
		key+="NoButterfly";
		break;
	case LinearButterfly:
		key+="LinearButterfly";
		break;
	case LogarithmButterfly:
		key+="LogarithmButterfly";
		break;
	case ContoursButterfly:
		key+="ContoursButterfly";
		break;
	case CausticButterfly:
		key+="CausticButterfly";
		break;
	}
	key+="_";
	key+=hypothesisname;

	return key;

}

/**
 * \brief TGraphAsymmError's residuals key
 */
std::string START::PlotFactory::MakeKeyGraphErrorsResidualsOff(std::string hypothesisname)
{

	std::string key;
	key = "TGraphAssymErrors_ResidualsOff";
	key+="_";
	switch(fButterfly) {
	case NoButterfly:
		key+="NoButterfly";
		break;
	case LinearButterfly:
		key+="LinearButterfly";
		break;
	case LogarithmButterfly:
		key+="LogarithmButterfly";
		break;
	case ContoursButterfly:
		key+="ContoursButterfly";
		break;
	case CausticButterfly:
		key+="CausticButterfly";
		break;
	}
	key+="_";
	key+=hypothesisname;

	return key;

}


/**
 * \brief TGraphAsymmError's residuals key
 */
std::string START::PlotFactory::MakeKeyMultiGraphContours(std::string hypothesisname, std::string param1, std::string param2)
{

	std::string key;
	key = "TMultiGraph_Contours";
	key+="_";
	key+=hypothesisname;
	key+="_";
	key+=param1;
	key+="_";
	key+=param2;
	return key;

}

/**
 * \brief TArrow's residuals key
 */
std::string START::PlotFactory::MakeKeyArrowsResiduals(std::string hypothesisname)
{

	std::string key;
	key = "TArrow_Residuals";
	key+="_";
	switch(fButterfly) {
	case NoButterfly:
		key+="NoButterfly";
		break;
	case LinearButterfly:
		key+="LinearButterfly";
		break;
	case LogarithmButterfly:
		key+="LogarithmButterfly";
		break;
	case ContoursButterfly:
		key+="ContoursButterfly";
		break;
	case CausticButterfly:
		key+="CausticButterfly";
		break;
	}
	key+="_";
	key+=hypothesisname;

	return key;

}

/**
 * \brief TArrow's spectrum key
 */
std::string START::PlotFactory::MakeKeyArrowsSpectrum(std::string hypothesisname)
{

	std::string key;
	key = "TArrow_Spectrum";
	key+="_";
	switch(fButterfly) {
	case NoButterfly:
		key+="NoButterfly";
		break;
	case LinearButterfly:
		key+="LinearButterfly";
		break;
	case LogarithmButterfly:
		key+="LogarithmButterfly";
		break;
	case ContoursButterfly:
		key+="ContoursButterfly";
		break;
	case CausticButterfly:
		key+="CausticButterfly";
		break;
	}
	key+="_";
	key+=hypothesisname;

	return key;

}

/**
 * \brief TPolyLine's spectrum key
 */
std::string START::PlotFactory::MakeKeyPolyLineSpectrum(std::string hypothesisname)
{

	std::string key;
	key = "TPolyLine_Spectrum";
	key+="_";
	switch(fButterfly) {
	case NoButterfly:
		key+="NoButterfly";
		break;
	case LinearButterfly:
		key+="LinearButterfly";
		break;
	case LogarithmButterfly:
		key+="LogarithmButterfly";
		break;
	case ContoursButterfly:
		key+="ContoursButterfly";
		break;
	case CausticButterfly:
		key+="CausticButterfly";
		break;
	}
	key+="_";
	key+=hypothesisname;

	return key;

}

/**
 * \brief TPolyLine's vectors spectrum key
 */
std::string START::PlotFactory::MakeKeyPolyLineVectorsSpectrum(std::string hypothesisname)
{

	std::string key;
	key = "TPolyLineVectors_Spectrum";
	key+="_";
	switch(fButterfly) {
	case NoButterfly:
		key+="NoButterfly";
		break;
	case LinearButterfly:
		key+="LinearButterfly";
		break;
	case LogarithmButterfly:
		key+="LogarithmButterfly";
		break;
	case ContoursButterfly:
		key+="ContoursButterfly";
		break;
	case CausticButterfly:
		key+="CausticButterfly";
		break;
	}
	key+="_";
	key+=hypothesisname;

	return key;

}

/**
 * \brief TGraph array spectrum key (butterfly with contours)
 */
std::string START::PlotFactory::MakeKeyGraphArraySpectrum(std::string hypothesisname) {

	std::string key;
	key = "TGraphArrayContoursSpectrum";
	key+="_";
	switch(fButterfly) {
	case NoButterfly:
		key+="NoButterfly";
		break;
	case LinearButterfly:
		key+="LinearButterfly";
		break;
	case LogarithmButterfly:
		key+="LogarithmButterfly";
		break;
	case ContoursButterfly:
		key+="ContoursButterfly";
		break;
	case CausticButterfly:
		key+="CausticButterfly";
		break;
	}
	key+="_";
	key+=hypothesisname;

	return key;

}

/**
 * \brief TF1's residuals residuals key
 */
std::string START::PlotFactory::MakeKeyTF1Residuals(std::string hypothesisname)
{

	std::string key;
	key = "TF1_Residuals";
	key+="_";
	switch(fButterfly) {
	case NoButterfly:
		key+="NoButterfly";
		break;
	case LinearButterfly:
		key+="LinearButterfly";
		break;
	case LogarithmButterfly:
		key+="LogarithmButterfly";
		break;
	case ContoursButterfly:
		key+="ContoursButterfly";
		break;
	case CausticButterfly:
		key+="CausticButterfly";
		break;
	}
	key+="_";
	key+=hypothesisname;

	return key;

}

/**
 * \brief canvas's likelihood key
 */
std::string START::PlotFactory::MakeKeyCanvasLikelihoodScan(std::string hypothesisname,std::string param) {

	std::string key;
	key = "TCanvas_LikelihoodScan";
	key+="_";
	key+=hypothesisname;
	key+="_";
	key+=param;

	return key;

}

/**
 * \brief Graph's likelihood key
 */
std::string START::PlotFactory::MakeKeyGraphLikelihoodScan(std::string hypothesisname,std::string param) {

	std::string key;
	key = "TGraph_LikelihoodScan";
	key+="_";
	key+=hypothesisname;
	key+="_";
	key+=param;

	return key;

}

/**
 * \brief canvas's light curves key
 */
std::string START::PlotFactory::MakeKeyCanvasLightCurves(std::string hypothesisname,std::string lctype) {

	std::string key;
	key = "TCanvas_LightCurves";
	key+="_";
	key+= lctype;
	key+="_";
	switch(fLightCurvesErrors) {
	case Gaussian:
		key+="GaussianErrors";
		break;
	case Rolke:
		key+="RolkeErrors";
	}
	key+="_";
	key+=hypothesisname;

	return key;

}

/**
 * \brief Pad's residuals key
 */
std::string START::PlotFactory::MakeKeyPadLightCurves(std::string hypothesisname,std::string lctype)
{

	std::string key;
	key = "TPad_LightCurves";
	key+="_";
	key+=lctype;
	key+="_";
	switch(fLightCurvesErrors) {
	case Gaussian:
		key+="GaussianErrors";
		break;
	case Rolke:
		key+="RolkeErrors";
	}
	key+="_";
	key+=hypothesisname;

	return key;

}

/**
 * \brief TH2's LC key
 */
std::string START::PlotFactory::MakeKeyTH2LightCurves(std::string hypothesisname,std::string lctype)
{

	std::string key;
	key = "TH2_LightCurves";
	key+="_";
	key+= lctype;
	key+="_";
	switch(fLightCurvesErrors) {
	case Gaussian:
		key+="GaussianErrors";
		break;
	case Rolke:
		key+="RolkeErrors";
	}
	key+="_";
	key+=hypothesisname;

	return key;

}

/**
 * \brief Graph's light curves key
 */
std::string START::PlotFactory::MakeKeyGraphLightCurves(std::string hypothesisname,std::string lctype) {

	std::string key;
	key = "TGraph_LightCurves";
	key+="_";
	key+=lctype;
	key+="_";
	switch(fLightCurvesErrors) {
	case Gaussian:
		key+="GaussianErrors";
		break;
	case Rolke:
		key+="RolkeErrors";
	}
	key+="_";
	key+=hypothesisname;

	return key;

}

/**
 * \brief TF1's mean flux light curve
 */
std::string START::PlotFactory::MakeKeyTF1MeanFluxLightCurves(std::string hypothesisname,std::string lctype) {

	std::string key;
	key = "TF1_MeanFluxLC";
	key+="_";
	key+=lctype;
	key+="_";
	switch(fLightCurvesErrors) {
	case Gaussian:
		key+="GaussianErrors";
		break;
	case Rolke:
		key+="RolkeErrors";
	}
	key+="_";
	key+=hypothesisname;

	return key;
}

/**
 * \brief TF1's zero flux light curve
 */
std::string START::PlotFactory::MakeKeyTF1ZeroFluxLightCurves(std::string hypothesisname,std::string lctype) {

	std::string key;
	key = "TF1_ZeroFluxLC";
	key+="_";
	key+=lctype;
	key+="_";
	switch(fLightCurvesErrors) {
	case Gaussian:
		key+="GaussianErrors";
		break;
	case Rolke:
		key+="RolkeErrors";
	}
	key+="_";
	key+=hypothesisname;

	return key;
}

/**
 * \brief canvas's light curves key
 */
std::string START::PlotFactory::MakeKeyArrowsLightCurves(std::string hypothesisname,std::string lctype) {

	std::string key;
	key = "TArrows_LightCurves";
	key+="_";
	key+= lctype;
	key+="_";
	switch(fLightCurvesErrors) {
	case Gaussian:
		key+="GaussianErrors";
		break;
	case Rolke:
		key+="RolkeErrors";
	}
	key+="_";
	key+=hypothesisname;

	return key;

}

/**
 * \brief TF1's spectrum key
 */
std::string START::PlotFactory::MakeKeyTF1FitSpectrum(std::string hypothesisname) {

	std::string key;
	key = "TF1_Spectrum";
	key+="_";
	key+=hypothesisname;

	return key;
}

/**
 * \brief PaveText's spectrum key
 */
std::string START::PlotFactory::MakeKeyPaveTextSpectrum(std::string hypothesisname) {

	std::string key;
	key = "PaveText_Spectrum";
	key+="_";
	switch(fButterfly) {
	case NoButterfly:
		key+="NoButterfly";
		break;
	case LinearButterfly:
		key+="LinearButterfly";
		break;
	case LogarithmButterfly:
		key+="LogarithmButterfly";
		break;
	case ContoursButterfly:
		key+="ContoursButterfly";
		break;
	case CausticButterfly:
		key+="CausticButterfly";
		break;
	}
	key+="_";
	key+=hypothesisname;

	return key;
}

/**
 * \brief Plot main canvas with all MWL attributs
 */
void START::PlotFactory::PlotMultiWaveLength(MultiWaveLengthFactory &Mwl)
{

	DEBUG_OUT << "START!" << std::endl;

	CleanCanvasMwl();
	CleanPadMwl();
	CleanTH2Mwl();
	CleanGraphErrorsMwl();
	CleanTF1Mwl();
	CleanArrowsMwl();
	CleanButterflyMwl();

	InitCanvasMwl(Mwl);
	InitPadMwl(Mwl);
	InitTH2Mwl(Mwl);
	InitGraphErrorsMwl(Mwl);
	InitTF1Mwl(Mwl);
	InitButterflyMwl(Mwl);

	DrawMultiWaveLength(Mwl);

	DEBUG_OUT << "END!" << std::endl;

}

/**
 * \brief Draw main canvas with all MWL attributs
 */
void START::PlotFactory::DrawMultiWaveLength(MultiWaveLengthFactory &Mwl) {

	DEBUG_OUT << "START!" << std::endl;

	// TCanvas
	std::string canvasname = MakeKeyCanvasMwl(Mwl.GetName());
	fMapCanvasMwl[canvasname]->Draw();
	fMapCanvasMwl[canvasname]->cd();
	DEBUG_OUT_L(2) << canvasname << " " << fMapCanvasMwl[canvasname] << std::endl;

	// TPad
	std::string padname = MakeKeyPadMwl(Mwl.GetName());
	fMapPadMwl[padname]->Draw();
	fMapPadMwl[padname]->cd();
	DEBUG_OUT_L(2) << padname << " " << fMapPadMwl[padname] << std::endl;

	// TH2
	std::string th2name = MakeKeyTH2Mwl(Mwl.GetName());
	fMapTH2Mwl[th2name]->DrawCopy();

	DEBUG_OUT_L(2) << th2name << " " << fMapTH2Mwl[th2name] << std::endl;

  
	// Butterflies
	for(std::map<std::string,TPolyLine*>::iterator butt=fMapButterflyMwl.begin();
		butt!=fMapButterflyMwl.end(); ++butt) {
		fMapPadMwl[padname]->cd();
		if(butt->second!=0) {
			DEBUG_OUT_L(2) << butt->first << " " << butt->second->GetN() << std::endl;
			butt->second->Draw("f");
			butt->second->Draw();
		}
	}

	// TF1
	for(std::map<std::string,TF1*>::iterator tf1=fMapTF1Mwl.begin();
		tf1!=fMapTF1Mwl.end(); ++tf1) {
		fMapPadMwl[padname]->cd();
		if(tf1->second!=0) {
			tf1->second->Draw("same");
			DEBUG_OUT_L(2) << tf1->first << " " << tf1->second->Eval(0.5) << std::endl;
		}
	}

	// Points
	for(std::map<std::string,TGraphAsymmErrors*>::iterator graph=fMapGraphErrorsMwl.begin();
		graph!=fMapGraphErrorsMwl.end(); ++graph) {
		fMapPadMwl[padname]->cd();
		if(graph->second!=0) {
			graph->second->Draw("P");
			DEBUG_OUT_L(2) << graph->first << " " << graph->second->GetN() << std::endl;
		}
	}

	// UpperLimits
	for(std::map<std::string,std::vector<TArrow*> >::iterator arrow=fMapArrowsMwl.begin(); arrow!=fMapArrowsMwl.end(); ++arrow) {
		for(unsigned int i(0); i<arrow->second.size(); i++) {
			if(arrow->second[i]!=0) {
				fMapPadMwl[padname]->cd();
				arrow->second[i]->Draw();
				DEBUG_OUT_L(2) << "arrow " << i+1 << std::endl;
			}
		}

	}

	// Draw source name
	fMapPadMwl[padname]->cd();
	std::ostringstream plotname;
	plotname << fsourcename;
	TLatex *drawsourcename = new TLatex();
	drawsourcename->SetNDC();
	drawsourcename->SetTextSize(0.032);
	drawsourcename->SetTextColor(1);
	drawsourcename->SetTextFont(52);
	switch(fPlotStyle) {
	case Paper:
		break;
	case Default:
	case UserFriendly:
		drawsourcename->DrawLatex(0.10,0.9,plotname.str().c_str());
	}

	// Legend
	fMapPadMwl[padname]->cd();  
	if(fLegendMwl) fLegendMwl->Draw();

	//fMapCanvasMwl[canvasname]->Update();
	fMapCanvasMwl[canvasname]->ForceUpdate();
  
	DEBUG_OUT << "END!" << std::endl;

}

/**
 * \brief Init Mwl
 */
void START::PlotFactory::InitCanvasMwl(MultiWaveLengthFactory &Mwl) {

	DEBUG_OUT << "START!" << std::endl;

	int canvaswidth(0), canvasheight(0); 
  
	switch(fPlotStyle) {
	case Default:
	case Paper:
	case UserFriendly:
		canvaswidth=800;
		canvasheight=600;
	}

	std::string canvasname = MakeKeyCanvasMwl(Mwl.GetName());
	fMapCanvasMwl[canvasname] = new TCanvas(canvasname.c_str(),canvasname.c_str(),canvaswidth,canvasheight);

	DEBUG_OUT << "END!" << std::endl;
}

/**
 * \brief Init Mwl
 */
void START::PlotFactory::InitPadMwl(MultiWaveLengthFactory &Mwl) {
	DEBUG_OUT << "START!" << std::endl;

	// Definition of spectrum's and residuals' pad size

	double x1mwl(0.),x2mwl(0.),y1mwl(0.),y2mwl(0.);

	switch(fPlotStyle) {
	case Default:
	case Paper:
	case UserFriendly:
    
		x1mwl=0.01;
		y1mwl=0.01;
		x2mwl=0.99;
		y2mwl=0.99;
		/*
		  x1mwl=0.1;
		  y1mwl=0.1;
		  x2mwl=0.9;
		  y2mwl=0.9;
		*/
	}

	std::string padname = MakeKeyPadMwl(Mwl.GetName());
	fMapPadMwl[padname] = new TPad(padname.c_str(),padname.c_str(),x1mwl,y1mwl,x2mwl,y2mwl,0);
	fMapPadMwl[padname]->SetLogx();
	fMapPadMwl[padname]->SetLogy();
	fMapPadMwl[padname]->SetGridx();
	fMapPadMwl[padname]->SetGridy();
	fMapPadMwl[padname]->SetFillColor(0);

	switch(fPlotStyle) {
	case Default:
	case Paper:
	case UserFriendly:
		//fMapPadMwl[padname]->SetTopMargin(0.1);
		//fMapPadMwl[padname]->SetBottomMargin(0.15);
		fMapPadMwl[padname]->SetRightMargin(0.02);
		DEBUG_OUT_L(2) << "hello" << std::endl;
	}

	DEBUG_OUT << "END!" << std::endl;
}

/**
 * \brief Init Mwl
 */
void START::PlotFactory::InitTH2Mwl(MultiWaveLengthFactory &Mwl) {

	DEBUG_OUT << "START!" << std::endl;

	std::string th2name = MakeKeyTH2Mwl(Mwl.GetName());
	std::pair<double,double> energyrange = Mwl.GetEnergyRange();
	std::pair<double,double> fluxrange = Mwl.GetFluxRange();

	double emin(energyrange.first), emax(energyrange.second);
	double fmin(fluxrange.first), fmax(fluxrange.second);

	//double fshift((fmax-fmin)*0.1);
	double fshift(0.4);
	double eshift(0.1);

	emin*=TMath::Power(10.,-eshift);
	emax*=TMath::Power(10.,eshift);
	fmin*=TMath::Power(10.,-fshift);
	fmax*=TMath::Power(10.,+fshift);

	DEBUG_OUT << "emin=" << emin << " emax=" << emax << std::endl;
	DEBUG_OUT << "fmin=" << fmin << " fmax=" << fmax << std::endl;

	fMapTH2Mwl[th2name] = new TH2D(th2name.c_str(),
								   th2name.c_str(),
								   100,emin,emax,
								   100,fmin,fmax);
	fMapTH2Mwl[th2name]->SetStats(kFALSE);
	fMapTH2Mwl[th2name]->SetTitle("");

	std::ostringstream xaxisname;
	MultiWaveLengthFactory::EnergyUnits EUnits = Mwl.GetEnergyUnits();
	switch(EUnits) {
	case MultiWaveLengthFactory::TeV:
		xaxisname << "Energy (TeV)";
		break;
	case MultiWaveLengthFactory::GeV:
		xaxisname << "Energy (GeV)";
		break;
	case MultiWaveLengthFactory::MeV:
		xaxisname << "Energy (MeV)";
		break;
	case MultiWaveLengthFactory::keV:
		xaxisname << "Energy (keV)";
		break;
	case MultiWaveLengthFactory::Hz:
		xaxisname << "#nu (Hz)";
	}

	std::ostringstream yaxisname;
	MultiWaveLengthFactory::FluxUnits FUnits = Mwl.GetFluxUnits();
	switch(FUnits) {
	case MultiWaveLengthFactory::vFv:
		yaxisname << "#nu F_{#nu}" << " " << "(erg.cm^{-2}.s^{-1})";
		break;
	case MultiWaveLengthFactory::ESquareFlux_TeV:
		yaxisname << "E^{2} #Phi(E)" << " " << "(TeV.cm^{-2}.s^{-1})";
		break;
	case MultiWaveLengthFactory::Differential_TeV:
		yaxisname << "dN/dE" << " " << "(cm^{-2}.s^{-1}.TeV^{-1})";
		break;
	case MultiWaveLengthFactory::Differential_GeV:
		yaxisname << "dN/dE" << " " << "(cm^{-2}.s^{-1}.GeV^{-1})";
		break;
	case MultiWaveLengthFactory::Differential_MeV:
		yaxisname << "dN/dE" << " " << "(cm^{-2}.s^{-1}.MeV^{-1})";
		break;
	case MultiWaveLengthFactory::Differential_keV:
		yaxisname << "dN/dE" << " " << "(cm^{-2}.s^{-1}.keV^{-1})";
	}

	switch(fPlotStyle) {
	case Default:
	case Paper:
	case UserFriendly:
		fMapTH2Mwl[th2name]->GetXaxis()->SetTitle(xaxisname.str().c_str());
		fMapTH2Mwl[th2name]->GetYaxis()->SetTitle(yaxisname.str().c_str());

		fMapTH2Mwl[th2name]->GetXaxis()->SetTitleSize(0.04);
		fMapTH2Mwl[th2name]->GetXaxis()->SetTitleOffset(1.20);
		fMapTH2Mwl[th2name]->GetXaxis()->SetLabelSize(0.03);
		fMapTH2Mwl[th2name]->GetXaxis()->SetLabelOffset(0.002);
		//fMapTH2Mwl[th2name]->GetXaxis()->SetMoreLogLabels();
    
		fMapTH2Mwl[th2name]->GetYaxis()->SetTitleSize(0.04);
		fMapTH2Mwl[th2name]->GetYaxis()->SetTitleOffset(1.12);
		fMapTH2Mwl[th2name]->GetYaxis()->SetLabelSize(0.03);
		fMapTH2Mwl[th2name]->GetXaxis()->SetLabelOffset(0.002);    
	}

	DEBUG_OUT << "END!" << std::endl;

}

/**
 * \brief Init Mwl
 */
void START::PlotFactory::InitGraphErrorsMwl(MultiWaveLengthFactory &Mwl) {

	DEBUG_OUT << "START!" << std::endl;

	std::map<std::string,Residuals> map = Mwl.GetMapResiduals();
  
	for(std::map<std::string,Residuals>::iterator CRes=map.begin(); CRes!=map.end(); ++CRes) {

		Residuals Res = CRes->second;

		Color_t MarkerColor = Res.GetMarkerColor();
		Size_t MarkerSize = Res.GetMarkerSize();
		Style_t MarkerStyle = Res.GetMarkerStyle();
		Color_t LineColor = Res.GetLineColor();
		Width_t LineWidth = Res.GetLineWidth();


		// TGraph
		std::string graphname = MakeKeyGraphErrorsMwl(Mwl.GetName(),CRes->first);
		fMapGraphErrorsMwl[graphname] = new TGraphAsymmErrors();
		fMapGraphErrorsMwl[graphname]->SetName(graphname.c_str());
		fMapGraphErrorsMwl[graphname]->SetMarkerSize(MarkerSize);
		fMapGraphErrorsMwl[graphname]->SetMarkerColor(MarkerColor);
		fMapGraphErrorsMwl[graphname]->SetMarkerStyle(MarkerStyle);
		fMapGraphErrorsMwl[graphname]->SetLineColor(LineColor);
		fMapGraphErrorsMwl[graphname]->SetLineWidth(LineWidth);

		// UpperLimits
		std::string arrowsname = MakeKeyArrowsMwl(Mwl.GetName(),CRes->first);
		fMapArrowsMwl[arrowsname].clear();

		// Attributs
		std::vector<bool> IsUpperLimits = Res.GetIsUpperLimits();
		std::vector<double> Energy = Res.GetEnergy();
		std::vector<double> EnergyMinus = Res.GetEnergyMinus();
		std::vector<double> EnergyPlus = Res.GetEnergyPlus();
		std::vector<double> Flux = Res.GetFlux();
		std::vector<double> FluxSigmaMinus = Res.GetFluxSigmaMinus();
		std::vector<double> FluxSigmaPlus = Res.GetFluxSigmaPlus();
		std::vector<double> Flux3SigmaPlus = Res.GetFlux3SigmaPlus();

		for(unsigned int i(0); i<Energy.size(); i++) {
      
			if(IsUpperLimits[i]==false) { // no uppper limit

				if(DEBUG>1) {
					std::cout << IsUpperLimits[i] << " " << Energy[i] << " " << EnergyMinus[i] << " " << EnergyPlus[i] << " "
							  << Flux[i] << " " << FluxSigmaMinus[i] << " " << FluxSigmaPlus[i] << " " << FluxSigmaPlus[i]
							  << std::endl;
				}

				fMapGraphErrorsMwl[graphname]->SetPoint(i,
														Energy[i],
														Flux[i]);

				if(FluxSigmaMinus[i]>0. && FluxSigmaPlus[i]>0. && EnergyMinus[i]>0. && EnergyPlus[i]>0.) { // if errors
					fMapGraphErrorsMwl[graphname]->SetPointError(i,
																 Energy[i]-EnergyMinus[i],
																 EnergyPlus[i]-Energy[i],
																 Flux[i]-FluxSigmaMinus[i],
																 FluxSigmaPlus[i]-FluxSigmaMinus[i]);
					DEBUG_OUT_L(2) << "CASE 1" << std::endl;
				}
				else if(FluxSigmaMinus[i]>0. && FluxSigmaPlus[i]>0.) { // if errors
					fMapGraphErrorsMwl[graphname]->SetPointError(i,
																 0.,
																 0.,
																 Flux[i]-FluxSigmaMinus[i],
																 FluxSigmaPlus[i]-FluxSigmaMinus[i]);
					DEBUG_OUT_L(2) << "CASE 2" << std::endl;
				}
				else if(EnergyMinus[i]>0. && EnergyPlus[i]>0.){
					fMapGraphErrorsMwl[graphname]->SetPointError(i,
																 Energy[i]-EnergyMinus[i],
																 EnergyPlus[i]-Energy[i],
																 0.,
																 0.);
					DEBUG_OUT_L(2) << "CASE 3" << std::endl;
				}
			}
			else {
				fMapArrowsMwl[arrowsname].push_back(new TArrow(Energy[i],
															   Flux3SigmaPlus[i],
															   Energy[i],
															   Flux3SigmaPlus[i]*fscalearrowfactor,
															   0.02,
															   "|-|>"));
				fMapArrowsMwl[arrowsname].back()->SetLineWidth(LineWidth);
				fMapArrowsMwl[arrowsname].back()->SetLineColor(MarkerColor);
				fMapArrowsMwl[arrowsname].back()->SetFillColor(MarkerColor);
				DEBUG_OUT_L(2) << "CASE 4" << std::endl;
			}

		}

	}

	DEBUG_OUT << "END!" << std::endl;

}


/**
 * \brief Init TF1 Mwl
 */
void START::PlotFactory::InitTF1Mwl(MultiWaveLengthFactory &Mwl) {

	DEBUG_OUT << "START!" << std::endl;

	std::map<std::string,TF1*> map = Mwl.GetMapTF1();
  
	for(std::map<std::string,TF1*>::iterator tf1=map.begin(); tf1!=map.end(); ++tf1) {

		std::string tf1name = MakeKeyTF1Mwl(Mwl.GetName(),tf1->first);
		fMapTF1Mwl[tf1name] = new TF1(*(tf1->second));
		fMapTF1Mwl[tf1name]->SetName(tf1name.c_str());
		fMapTF1Mwl[tf1name]->SetTitle(tf1name.c_str());
	}

	DEBUG_OUT << "END!" << std::endl;

}

/**
 * \brief Init Butterfly Mwl
 */
void START::PlotFactory::InitButterflyMwl(MultiWaveLengthFactory &Mwl) {

	DEBUG_OUT << "START!" << std::endl;

	std::map<std::string,TPolyLine*> map = Mwl.GetMapButterfly();
  
	for(std::map<std::string,TPolyLine*>::iterator butt=map.begin(); butt!=map.end(); ++butt) {

		std::string buttname = MakeKeyButterflyMwl(Mwl.GetName(),butt->first);
		fMapButterflyMwl[buttname] = new TPolyLine(*(butt->second));
	}

	DEBUG_OUT << "END!" << std::endl;

}

/**
 * \brief Init Butterfly Mwl
 */
void START::PlotFactory::InitLegendMwl(MultiWaveLengthFactory &Mwl) {

	DEBUG_OUT << "START!" << std::endl;

	std::map<std::string,TPolyLine*> map = Mwl.GetMapButterfly();
  
	for(std::map<std::string,TPolyLine*>::iterator butt=map.begin(); butt!=map.end(); ++butt) {

	}

	DEBUG_OUT << "END!" << std::endl;

}

/**
 * \brief Clean Canvas Mwl
 */
void START::PlotFactory::CleanCanvasMwl() {
	for(std::map<std::string,TCanvas*>::iterator obj=fMapCanvasMwl.begin(); obj!=fMapCanvasMwl.end(); ++obj) {
		if(obj->second!=0) delete obj->second;
		obj->second=0;
	}
	fMapCanvasMwl.clear();
}


/**
 * \brief Clean TPad Mwl
 */
void START::PlotFactory::CleanPadMwl() {
	for(std::map<std::string,TPad*>::iterator obj=fMapPadMwl.begin(); obj!=fMapPadMwl.end(); ++obj) {
		if(obj->second!=0) delete obj->second;
		obj->second=0;
	}
	fMapPadMwl.clear();
}


/**
 * \brief Clean TH2 Mwl
 */
void START::PlotFactory::CleanTH2Mwl() {
	for(std::map<std::string,TH2*>::iterator obj=fMapTH2Mwl.begin(); obj!=fMapTH2Mwl.end(); ++obj) {
		if(obj->second!=0) delete obj->second;
		obj->second=0;
	}
	fMapTH2Mwl.clear();
}


/**
 * \brief Clean Graph Mwl
 */
void START::PlotFactory::CleanGraphErrorsMwl() {
	for(std::map<std::string,TGraphAsymmErrors*>::iterator obj=fMapGraphErrorsMwl.begin(); obj!=fMapGraphErrorsMwl.end(); ++obj) {
		if(obj->second!=0) delete obj->second;
		obj->second=0;
	}
	fMapGraphErrorsMwl.clear();
}

/**
 * \brief Clean Arrows Mwl
 */
void START::PlotFactory::CleanArrowsMwl() {
	for(std::map<std::string,std::vector<TArrow*> >::iterator obj=fMapArrowsMwl.begin(); obj!=fMapArrowsMwl.end(); ++obj) {
		for(unsigned int i(0); i<obj->second.size(); i++) {
			if(obj->second[i]!=0) delete obj->second[i];
			obj->second[i]=0;
		}
		obj->second.clear();
	}
	fMapArrowsMwl.clear();
}

/**
 * \brief Clean TF1 Mwl
 */
void START::PlotFactory::CleanTF1Mwl() {
	for(std::map<std::string,TF1*>::iterator obj=fMapTF1Mwl.begin(); obj!=fMapTF1Mwl.end(); ++obj) {
		if(obj->second!=0) delete obj->second;
		obj->second=0;
	}
	fMapTF1Mwl.clear();
}

/**
 * \brief Clean Butterfly Mwl
 */
void START::PlotFactory::CleanButterflyMwl() {
	for(std::map<std::string,TPolyLine*>::iterator obj=fMapButterflyMwl.begin(); obj!=fMapButterflyMwl.end(); ++obj) {
		if(obj->second!=0) delete obj->second;
		obj->second=0;
	}
	fMapButterflyMwl.clear();
}

/**
 * \brief Make key Canvas Mwl
 */
std::string START::PlotFactory::MakeKeyCanvasMwl(std::string mwlname) {
	std::ostringstream key;
	key << mwlname.c_str() << "_Canvas"; 
	return key.str();
}

/**
 * \brief Make key Pad Mwl
 */
std::string START::PlotFactory::MakeKeyPadMwl(std::string mwlname) {
	std::ostringstream key;
	key << mwlname.c_str() << "_Pad"; 
	return key.str();
}

/**
 * \brief Make key TH2 Mwl
 */
std::string START::PlotFactory::MakeKeyTH2Mwl(std::string mwlname) {
	std::ostringstream key;
	key << mwlname.c_str() << "_TH2"; 
	return key.str();
}

/**
 * \brief Make key TGraphErrors Mwl
 */
std::string START::PlotFactory::MakeKeyGraphErrorsMwl(std::string mwlname,std::string graphname) {
	std::ostringstream key;
	key << mwlname.c_str() << "_GraphErrors_" << graphname.c_str(); 
	return key.str();
}

/**
 * \brief Make key Arrows  Mwl
 */
std::string START::PlotFactory::MakeKeyArrowsMwl(std::string mwlname,std::string graphname) {
	std::ostringstream key;
	key << mwlname.c_str() << "_Arrows_" << graphname.c_str(); 
	return key.str();
}

/**
 * \brief Make key TF1  Mwl
 */
std::string START::PlotFactory::MakeKeyTF1Mwl(std::string mwlname,std::string tf1name) {
	std::ostringstream key;
	key << mwlname.c_str() << "_TF1_" << tf1name.c_str(); 
	return key.str();
}

/**
 * \brief Make key Butterfly  Mwl
 */
std::string START::PlotFactory::MakeKeyButterflyMwl(std::string mwlname,std::string butterfly) {
	std::ostringstream key;
	key << mwlname.c_str() << "_Butterfly_" << butterfly.c_str(); 
	return key.str();
}

/**
 * \brief Add legend for points
 * \param Mwl MultiWaveLengthFactory object used for MWL data creation
 * \param name Name given to the points in MultiWaveLengthFactory
 * \param label Text associated with points
 * \param options Options for legend
 * Options can take the following values :
 * <ul> 
 * <li> "L" draw line associated with TAttLine if obj inherits from TAttLine </li>
 * <li> "P" draw polymarker associated with TAttMarker if obj inherits from TAttMarker </li>
 * <li> "F" draw a box with fill associated wit TAttFill if obj inherits TAttFill </li>
 * <li> "E" draw vertical error bar if option "L" is also specified </li>
 * </ul> 
 */
void START::PlotFactory::AddPointsLegendEntryMwl(MultiWaveLengthFactory &Mwl,std::string name, std::string label, std::string options) {

	DEBUG_OUT << "START!" << std::endl;

	if(fLegendMwl==0) {
		WARN_OUT << "You should call START::PlotFactory::SetMwlLegendAttributs before using this function" << std::endl;
		INFO_OUT << "Legend will be initialized with default attributs" << std::endl;
		SetMwlLegendAttributs();
	} 

	for(std::map<std::string,TGraphAsymmErrors*>::iterator graph=fMapGraphErrorsMwl.begin(); 
		graph!=fMapGraphErrorsMwl.end(); ++graph) {
		std::string iname = MakeKeyGraphErrorsMwl(Mwl.GetName(),name);
		DEBUG_OUT << "Add Entry" << iname << "?" << std::endl;
		if(graph->first==iname && graph->second!=0) {
			DEBUG_OUT << "Entry" << iname << " ok!" << std::endl;
			fLegendMwl->AddEntry(fMapGraphErrorsMwl[iname],label.c_str(),options.c_str());
		}
	}

	std::string canvasname = MakeKeyCanvasMwl(Mwl.GetName());
	fMapCanvasMwl[canvasname]->ForceUpdate();

	DEBUG_OUT << "END!" << std::endl;

}

/**
 * \brief Add legend for Upperlimits
 * \param Mwl MultiWaveLengthFactory object used for MWL data creation
 * \param name Name given to the points in MultiWaveLengthFactory
 * \param label Text associated with points
 * \param options Options for legend
 * Options can take the following values :
 * <ul> 
 * <li> "L" draw line associated with TAttLine if obj inherits from TAttLine </li>
 * <li> "P" draw polymarker associated with TAttMarker if obj inherits from TAttMarker </li>
 * <li> "F" draw a box with fill associated wit TAttFill if obj inherits TAttFill </li>
 * <li> "E" draw vertical error bar if option "L" is also specified </li>
 * </ul> 
 */
void START::PlotFactory::AddUpperLimitsLegendEntryMwl(MultiWaveLengthFactory &Mwl,std::string name, std::string label, std::string options) {

	if(fLegendMwl==0) {
		WARN_OUT << "You should call START::PlotFactory::SetMwlLegendAttributs before using this function" << std::endl;
		INFO_OUT << "Legend will be initialized with default attributs" << std::endl;
		SetMwlLegendAttributs();
	} 

	for(std::map<std::string,std::vector<TArrow*> >::iterator graph=fMapArrowsMwl.begin(); 
		graph!=fMapArrowsMwl.end(); ++graph) {
		std::string iname = MakeKeyArrowsMwl(Mwl.GetName(),name);
		if(graph->second.size()>0) {
			if(graph->first==iname && graph->second[0]!=0) {
				fLegendMwl->AddEntry(fMapArrowsMwl[iname][0],label.c_str(),options.c_str());
				break;
			}
		}
	}

	std::string canvasname = MakeKeyCanvasMwl(Mwl.GetName());
	fMapCanvasMwl[canvasname]->ForceUpdate();

}

/**
 * \brief Add legend for TF1
 * \param Mwl MultiWaveLengthFactory object used for MWL data creation
 * \param name Name given to the points in MultiWaveLengthFactory
 * \param label Text associated with points
 * \param options Options for legend
 * Options can take the following values :
 * <ul> 
 * <li> "L" draw line associated with TAttLine if obj inherits from TAttLine </li>
 * <li> "P" draw polymarker associated with TAttMarker if obj inherits from TAttMarker </li>
 * <li> "F" draw a box with fill associated wit TAttFill if obj inherits TAttFill </li>
 * <li> "E" draw vertical error bar if option "L" is also specified </li>
 * </ul> 
 */
void START::PlotFactory::AddTF1LegendEntryMwl(MultiWaveLengthFactory &Mwl,std::string name, std::string label, std::string options) {

	if(fLegendMwl==0) {
		WARN_OUT << "You should call START::PlotFactory::SetMwlLegendAttributs before using this function" << std::endl;
		INFO_OUT << "Legend will be initialized with default attributs" << std::endl;
		SetMwlLegendAttributs();
	} 

	DEBUG_OUT << "size map tf1 " << fMapTF1Mwl.size() << std::endl;

	for(std::map<std::string,TF1*>::iterator graph=fMapTF1Mwl.begin(); 
		graph!=fMapTF1Mwl.end(); ++graph) {
		std::string iname = MakeKeyTF1Mwl(Mwl.GetName(),name);
		DEBUG_OUT << "Add Entry" << iname << "?" << std::endl;
		if(graph->first==iname && graph->second!=0) {
			fLegendMwl->AddEntry(fMapTF1Mwl[iname],label.c_str(),options.c_str());
			DEBUG_OUT << "Entry" << iname << " ok!" << std::endl;
		}
	}

	std::string canvasname = MakeKeyCanvasMwl(Mwl.GetName());
	fMapCanvasMwl[canvasname]->ForceUpdate();

}

/**
 * \brief Add legend for Butterfly
 * \param Mwl MultiWaveLengthFactory object used for MWL data creation
 * \param name Name given to the points in MultiWaveLengthFactory
 * \param label Text associated with points
 * \param options Options for legend
 * Options can take the following values :
 * <ul> 
 * <li> "L" draw line associated with TAttLine if obj inherits from TAttLine </li>
 * <li> "P" draw polymarker associated with TAttMarker if obj inherits from TAttMarker </li>
 * <li> "F" draw a box with fill associated wit TAttFill if obj inherits TAttFill </li>
 * <li> "E" draw vertical error bar if option "L" is also specified </li>
 * </ul> 
 */
void START::PlotFactory::AddButterflyLegendEntryMwl(MultiWaveLengthFactory &Mwl,std::string name, std::string label, std::string options) {

	if(fLegendMwl==0) {
		WARN_OUT << "You should call START::PlotFactory::SetMwlLegendAttributs before using this function" << std::endl;
		INFO_OUT << "Legend will be initialized with default attributs" << std::endl;
		SetMwlLegendAttributs();
	} 

	for(std::map<std::string,TPolyLine*>::iterator graph=fMapButterflyMwl.begin(); 
		graph!=fMapButterflyMwl.end(); ++graph) {
		std::string iname = MakeKeyButterflyMwl(Mwl.GetName(),name);
		if(graph->first==iname && graph->second!=0) {
			fLegendMwl->AddEntry(fMapButterflyMwl[iname],label.c_str(),options.c_str());
		}
	}

	std::string canvasname = MakeKeyCanvasMwl(Mwl.GetName());
	fMapCanvasMwl[canvasname]->ForceUpdate();

}


void START::PlotFactory::SetMwlLegendAttributs(double x1, double y1, double x2, double y2, std::string title) {
	fLegendMwl = new TLegend(x1,y1,x2,y2,title.c_str());
	fLegendMwl->SetFillColor(0);
	fLegendMwl->SetFillStyle(1001);
	fLegendMwl->SetLineColor(1);
}

void START::PlotFactory::UpdateCanvasMwl(MultiWaveLengthFactory &Mwl) {
	std::string canvasname = MakeKeyCanvasMwl(Mwl.GetName());
	fMapCanvasMwl[canvasname]->ForceUpdate();
}


