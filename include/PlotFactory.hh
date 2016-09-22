#ifndef _PLOTFACTORY_
#define _PLOTFACTORY_

// STL
#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <map>

// ROOT
#include <TObject.h>
#include <TString.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH2.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TGraphPainter.h>
#include <TMultiGraph.h>
#include <TGraphAsymmErrors.h>
#include <TPad.h>
#include <TArrow.h>
#include <TPaveText.h>
#include <TPolyLine.h>

// START
#include "MultiWaveLengthFactory.hh"

namespace START {
	class Hypothesis;

	/**
	 * 
	 * \brief Plot spectrum, residuals, confidence intervals, likelihood scans and light curve
	 * for an Hypothesis or a vector of Hypothesis.
	 *
	 * There are six user's functions : PlotContours, PlotSpectrumAndResiduals, RebinHypothesis
	 * PlotLikelihoodScans, PlotLightCurves and SaveObjectsInROOTFile.
	 *
	 * <ul>
	 * <li> PlotContours is used to draw contours at 68%, 95% and 99% confidence intervals </li>
	 * <li> PlotSpectrumAndResiduals is used to draw the hypothesis'spectrum with butterfly and its residuals.
	 * There is four options for the butterfly : </li>
	 *    <ul>
	 *    <li> NoButterfly : no butterfly is plotted </li>
	 *    <li> LinearButterfly : butterfly is plotted with covariance of first order </li>
	 *    <li> LogButterfly : butterfly is plotted with covariance of first order (Fermi way) </li>
	 *    <li> ContoursButterfly : butterfly is plotted with contours. They must be computed before. </li>
	 *    <li> CausticButterfly : causted is plotted with contours (with TGraph). They must be computed before. </li>
	 *    </ul>
	 * <li> PlotLikelihoodScans is used to draw likelihood scans </li>
	 * <li> RebinHypothesis is used to re-compute and draw again the points used in the spectrum and the residuals's plot </li>
	 * <li> SaveObjectsInROOTFile saves ROOT objects (plots, graph, legend...) </li>
	 * </ul>
	 *
	 * \author HAP-Fr team
	 */
	class PlotFactory : public TObject
	{

	public:

		/**
		 * \brief Options for the butterfly's drawing on the spectrum canvas.
		 * \param NoButterfly butterfly is not drawn
		 * \param LinearButterfly is drawn with the covaricance of the flux at 1st order
		 * \param LogarithmButterfly butterfly is drawn with the log10 covaricance of the flux at 1st order
		 * \param ContoursButterfly butterfly is drawn with the contours which have to be computed before
		 */
		typedef enum {NoButterfly, LinearButterfly, LogarithmButterfly, ContoursButterfly, CausticButterfly} ButterflyOption;

		/**
		 * \brief Options for the plot style.
		 * \param Paper
		 * \param FriendlyUser
		 */
		typedef enum {Paper, UserFriendly, Default} PlotStyle;

		/**
		 * \brief Options for the plot style.
		 * \param Paper
		 * \param FriendlyUser
		 */
		typedef enum {Excess, ExcessOff, OnOff} ResidualsStyle;
	  
	  
		/**
		 * \brief Options for light curve time axis
		 */
		typedef enum {MJD, YearMonth, YearMonthDay, HourMinuteSecond} LightCurvesTimeAxis;

		/**
		 * \brief Options for light curve errors handling
		 * \param Gaussian Errors handled with simple errors propagation
		 * \param Rolke Errors handled with Rolke
		 */
		typedef enum {Gaussian,Rolke} LightCurvesErrors;

		PlotFactory(std::string sourcename="", PlotStyle Style=UserFriendly,
					ResidualsStyle ResStyle = Excess, bool verbose=false);
		PlotFactory(std::vector<Hypothesis*> &hypo, std::string sourcename="",
					PlotStyle Style=Default, ResidualsStyle ResStyle=Excess, bool verbose=false); // constructor to analyse a set of hypothesis
		PlotFactory(Hypothesis &hypo, std::string sourcename="",
					PlotStyle Style=UserFriendly, ResidualsStyle ResStyle=Excess, bool verbose=false); // constructor to analyse one hypothesis

		~PlotFactory(); // destructor

		void PlotContours(); // Plot contours at 1, 2 and 3 sigmas 

		void PlotSpectrumAndResiduals(ButterflyOption butterfly=NoButterfly, bool drawfit=false); // Plot spectrum, residuals and butterfly

		void PlotLikelihoodScans(); // Plot Likelihood scans

		void PlotLightCurves(LightCurvesTimeAxis TimeAxisOption=MJD, LightCurvesErrors ErrorsOption=Rolke); // Plot LightCurves

		void RebinHypothesis(Hypothesis &hypo, double rebining=2., ButterflyOption butterfly=NoButterfly, bool drawfit=false); // rebin one hypothesis

		void SaveObjectsInROOTFile(TString rootfile); // Save objects in a root file

		void SetMinimalInfosOnPad(bool IsMinimalInfosOnPad) {fminimalinfosonpad=IsMinimalInfosOnPad;}; 

		void PlotMultiWaveLength(MultiWaveLengthFactory &Mwl);
		void AddPointsLegendEntryMwl(MultiWaveLengthFactory &Mwl,std::string name, std::string label, std::string options);
		void AddButterflyLegendEntryMwl(MultiWaveLengthFactory &Mwl,std::string name, std::string label, std::string options);
		void AddTF1LegendEntryMwl(MultiWaveLengthFactory &Mwl,std::string name, std::string label, std::string options);
		void AddUpperLimitsLegendEntryMwl(MultiWaveLengthFactory &Mwl,std::string name, std::string label, std::string options);
		void UpdateCanvasMwl(MultiWaveLengthFactory& Mwl);
		void SetMwlLegendAttributs(double x1=0.7, double y1=0.7, double x2=0.9, double y2=0.9, std::string title="");

	private:

		/*Private functions*/

		void SaveSpectrumAndResidualsROOTObjects();
		void SaveContoursROOTObjects();
		void SaveLikelihoodScansROOTObjects();
		void SaveLightCurvesROOTObjects();

		void DrawSpectrumAndResiduals(const Hypothesis &hypo); // Draw objects
		void DrawLikelihoodScans(const Hypothesis &hypo); // Draw Likelihood scans
		void DrawLightCurves(const Hypothesis &hypo,std::string lctype); // Draw light curves 

		// Return the TH2's residuals height
		double GetPlotResidualsMaximum(const std::vector<double> &residualssigmaminus, const std::vector<double> &residualssigmaplus);

		// Set flux range for plot
		std::pair<double,double> GetPlotFluxRange(Hypothesis &hypothesis, const std::vector<double> &energybin, 
												  const std::vector<double> &flux, const std::vector<double> &fluxsigmaminus,
												  const std::vector<double> &fluxsigmaplus, const std::vector<double> &flux3sigmaplus);

		std::pair<double,double> GetPlotEnergyRange(const Hypothesis &hypothesis); //< Set energy range for TH2's

		std::pair<double,double> GetLightCurvesPlotTimeRange(const Hypothesis &hypothesis,std::string lctype); // Set time range for TH2's
		std::pair<double,double> GetLightCurvesPlotFluxRange(const Hypothesis &hypothesis,std::string lctype); // Set flux range for TH2's

		std::pair<double,double> GetButterflyEnergyRange(const Hypothesis &hypothesis); // Set energy range for butterfly plot

		void BadPointsKiller(std::vector<double> &energybin, std::vector<double> &residuals, std::vector<double> &residualssigmaplus, 
							 std::vector<double> &residualssigmaminus,std::vector<double> &residuals3sigmaplus, 
							 std::vector<double> &residuals3sigmaminus, std::vector<double> &flux, std::vector<double> &fluxsigmaplus,
							 std::vector<double> &fluxsigmaminus, std::vector<double> &flux3sigmaplus, std::vector<double> &flux3sigmaminus,
							 bool quiet); // Kill bad points on the spectrum and residuals plot

		void InitCanvasSpectrum(const Hypothesis &hypo); // build TCanvas
		void InitPads(const Hypothesis &hypo); // build TPads
		void InitTH2(Hypothesis &hypo); // build TH2
		void InitAndFillGraphs(const Hypothesis &hypo); // build and fill TGraphAsymmErrors
		void InitPolyLines(Hypothesis &hypo); // build butterfly
		void InitTF1Residuals(const Hypothesis &hypo); /// build equation y=0 for residuals
		void InitTF1FitSpectrum(const Hypothesis &hypo); /// build hypothesis fit
		void InitCanvasLikelihoodScan(const Hypothesis &hypo); // build TCanvas
		void InitGraphLikelihoodScan(const Hypothesis &hypo); // build TCanvas
		void InitPaveTextSpectrum(Hypothesis &hypo); // build TPaveText

		void InitGraphLightCurves(const Hypothesis &hypo,std::string lctype); // build TCanvas
		void InitCanvasLightCurves(const Hypothesis &hypo,std::string lctype); // build TCanvas
		void InitTH2LightCurves(Hypothesis &hypo,std::string lctype);
		void InitPadLightCurves(Hypothesis &hypo,std::string lctype);
		void InitTF1MeanFluxLightCurves(Hypothesis &hypo,std::string lctype); // build TF1
		void InitTF1ZeroFluxLightCurves(Hypothesis &hypo,std::string lctype); // build TF1

		void CleanCanvasSpectrum(); // delete TCanvas
		void CleanCanvasContours(); // delete TCanvas
		void CleanCanvasLikelihoodScan(); // delete TCanvas
		void CleanPads(); // delete TPads
		void CleanTH2(); // delete TH2
		void CleanGraphSpectrum(); // delete TGraphAsymmErrors
		void CleanMultiGraphs(); // delete TMultiGraphs
		void CleanPolyLines(); // delete TPolyLine (butterfly)
		void CleanTF1Residuals(); // delete TF1
		void CleanTF1FitSpectrum(); // delete TF1
		void CleanGraphScanLikelihood(); // delete TMultiGraphs
		void CleanGraphArraySpectrum(); // delete array of graph for contours
		void CleanPaveTextSpectrum(); // delete pavetext

		void CleanGraphLightCurves(); // delete TMultiGraphs
		void CleanCanvasLightCurves(); // delete TCanvas
		void CleanTH2LightCurves(); // delete TH2
		void CleanPadLightCurves(); // delete Pad
		void CleanTF1MeanFluxLightCurves(); // delete tf1
		void CleanTF1ZeroFluxLightCurves(); // delete tf1

		void PlotHypothesis(Hypothesis &hypo); // function to plot a hypothesis.

		void ComputeLinearButterfly(Hypothesis &hypo,std::vector<double> &x, std::vector<double> &y,
									unsigned int nbpoints, double log10emin, double log10emax);
		void ComputeLogarithmButterfly(Hypothesis &hypo,std::vector<double> &x, std::vector<double> &y,
									   unsigned int nbpoints, double log10emin, double log10emax);
		void ComputeCausticButterfly(Hypothesis &hypo,double log10emin, double log10emax, unsigned int nbpoints);
		void ComputeContoursButterfly(Hypothesis &hypo,std::vector<double> &x, std::vector<double> &y,
									  double log10emin, double log10emax, unsigned int nbpoints);


		// Key makers to play with maps

		std::string MakeKeyCanvasSpectrum(std::string hypothesisname);
		std::string MakeKeyCanvasContours(std::string hypothesisname, std::string param1, std::string param2);

		std::string MakeKeyPadSpectrum(std::string hypothesisname);
		std::string MakeKeyPadResiduals(std::string hypothesisname);

		std::string MakeKeyPaveTextSpectrum(std::string hypothesisname);

		std::string MakeKeyTH2Spectrum(std::string hypothesisname);
		std::string MakeKeyTH2Residuals(std::string hypothesisname);

		std::string MakeKeyGraphErrorsSpectrum(std::string hypothesisname);
		std::string MakeKeyGraphErrorsResiduals(std::string hypothesisname);
		std::string MakeKeyGraphErrorsResidualsOn(std::string hypothesisname);
		std::string MakeKeyGraphErrorsResidualsOff(std::string hypothesisname);

		std::string MakeKeyMultiGraphContours(std::string hypothesisname, std::string param1, std::string param2);

		std::string MakeKeyArrowsResiduals(std::string hypothesisname);
		std::string MakeKeyArrowsSpectrum(std::string hypothesisname);

		std::string MakeKeyPolyLineSpectrum(std::string hypothesisname);
		std::string MakeKeyPolyLineVectorsSpectrum(std::string hypothesisname);
		std::string MakeKeyGraphArraySpectrum(std::string hypothesisname);

		std::string MakeKeyTF1Residuals(std::string hypothesisname);

		std::string MakeKeyTF1FitSpectrum(std::string hypothesisname);

		std::string MakeKeyCanvasLikelihoodScan(std::string hypothesisname,std::string param);
		std::string MakeKeyGraphLikelihoodScan(std::string hypothesisname,std::string param);

		std::string MakeKeyCanvasLightCurves(std::string hypothesisname,std::string lctype);
		std::string MakeKeyGraphLightCurves(std::string hypothesisname,std::string lctype);
		std::string MakeKeyTF1MeanFluxLightCurves(std::string hypothesisname,std::string lctype);
		std::string MakeKeyTF1ZeroFluxLightCurves(std::string hypothesisname,std::string lctype);
		std::string MakeKeyTH2LightCurves(std::string hypothesisname,std::string lctype);
		std::string MakeKeyPadLightCurves(std::string hypothesisname,std::string lctype);
		std::string MakeKeyArrowsLightCurves(std::string hypothesisname,std::string lctype);


		/* Private attributs */

		bool fverbose;

		std::string fsourcename; ///< source name

		ButterflyOption fButterfly;
		PlotStyle fPlotStyle;
		ResidualsStyle fResidualsStyle;
		LightCurvesTimeAxis fLightCurvesTimeAxis;
		LightCurvesErrors fLightCurvesErrors;
		bool fdrawfit;

		bool fminimalinfosonpad;

		TFile *fSaveRootFile; ///< ROOT savefile

		double fscalearrowfactor; ///< arrow scaling for upperlimit

		// TPad height and width
		double fx1padspec; ///< spectrum pad x1
		double fy1padspec; ///< spectrum pad y1
		double fx2padspec; ///< spectrum pad x2
		double fy2padspec; ///< spectrum pad y2
		double fx1padres; ///< residulas pad x1
		double fy1padres; ///< residulas pad y1
		double fx2padres; ///< residulas pad x2
		double fy2padres; ///< residulas pad y2
		double fx1padlc; ///< lc pad x1
		double fy1padlc; ///< lc pad y1
		double fx2padlc; ///< lc pad x2
		double fy2padlc; ///< lc pad y2
		// Pad : 
		//          (x2,y2)
		// (x1,y1)

		//Residuals and spectrum's points temporary

		std::vector<double> fenergybin; ///< current energy
		std::vector<double> fresiduals; ///< current residuals
		std::vector<double> fresiduals_on; ///< current residuals
		std::vector<double> fresiduals_off; ///< current residuals
		std::vector<double> fresidualssigmaplus; ///< current residuals 1 sigma plus 
		std::vector<double> fresidualssigmaminus; ///< current residuals 1 sigma minus 
		std::vector<double> fresiduals3sigmaplus; ///< current residuals 3 sigma plus 
		std::vector<double> fresiduals3sigmaminus; ///< current residuals 3 sigma plus 
		std::vector<double> fflux; ///< current flux 
		std::vector<double> ffluxsigmaplus; ///< current flux 1 sigma plus 
		std::vector<double> ffluxsigmaminus; ///< current flux 1 sigma minus 
		std::vector<double> fflux3sigmaplus; ///< current flux 3 sigma plus 
		std::vector<double> fflux3sigmaminus; ///< current flux 3 sigma minus 

		// useful objects 

		std::vector<Hypothesis*> fHypothesisArray; ///< hypothesis
		std::vector<Hypothesis*> fRebinnedHypothesisArray; ///< rebinned hypothesis

		// saved objects, one element in a map is defined by the name of an hypothesis

		std::map<std::string,TCanvas*> fMapCanvasContours; ///< contours' TCanvas
		std::map<std::string,TCanvas*> fMapCanvasSpectrumAndResiduals; ///< spectrum's, residuals' and butterfly's TCanvas

		std::map<std::string,TPaveText*> fMapPaveTextSpectrum; ///< map of TPaveText for main spectrum canvas if needed

		std::map<std::string,TPad*> fMapPadSpectrum; ///< spectrum's pad
		std::map<std::string,TPad*> fMapPadResiduals; ///< residuals's pad

		std::map<std::string,TH2D*> fMapTH2Spectrum; ///< spectrum's TH2F
		std::map<std::string,TH2D*> fMapTH2Residuals; ///< residuals's TH2F

		std::map<std::string,TGraphAsymmErrors*> fMapGraphErrorsSpectrum; ///< spectrum's TGraphAsymmErrors
		std::map<std::string,TGraphAsymmErrors*> fMapGraphErrorsResiduals; ///< residuals's TGraphAsymmErrors
		std::map<std::string,TGraphAsymmErrors*> fMapGraphErrorsResidualsOn; ///< residuals's TGraphAsymmErrors
		std::map<std::string,TGraphAsymmErrors*> fMapGraphErrorsResidualsOff; ///< residuals's TGraphAsymmErrors

		std::map<std::string,TMultiGraph*> fMapMultiGraphContours; ///< contours' TMultiGraph

		std::map<std::string,std::vector<TArrow*> > fMapArrowsSpectrum; ///< spectrum's lower limits arrows
		std::map<std::string,std::vector<TArrow*> > fMapArrowsResiduals; ///< spectrum's lower limits arrows

		std::map<std::string,TPolyLine*> fMapPolyLineSpectrum; ///< spectrum's butterflies 
		std::map<std::string,std::vector<TGraph*> > fMapGraphArraySpectrum; ///< spectrum's butterflies 

		std::map<std::string,TF1*> fMapTF1Residuals; ///< spectrum's butterflies 

		std::map<std::string,TF1*> fMapTF1FitSpectrum; ///< spectrum's fit

		std::map<std::string,TCanvas*> fMapCanvasLikelihoodScan; ///< scan canvas
		std::map<std::string,TGraph*> fMapGraphLikelihoodScan; ///< scan graphs

		std::map<std::string,std::pair<std::vector<double>,std::vector<double> > > fMapButterFlyVectors; ///< vector used to buil butterfly

		std::map<std::string,TCanvas*> fMapCanvasLightCurves; ///< light curves canvas
		std::map<std::string,TGraphAsymmErrors*> fMapGraphLightCurves; ///< light curves graphs
		std::map<std::string,TH2D*> fMapTH2LightCurves; ///< th2 LC
		std::map<std::string,TPad*> fMapPadLightCurves; ///< pad light curve
		std::map<std::string,TF1*> fMapTF1MeanFluxLightCurves; ///< mean flux TF1 line
		std::map<std::string,TF1*> fMapTF1ZeroFluxLightCurves; ///< mean flux TF1 line
		std::map<std::string,std::vector<TArrow*> > fMapArrowsLightCurves; ///< lc's upper limits

		// MWL
		std::map<std::string,TCanvas*> fMapCanvasMwl;
		std::map<std::string,TPad*> fMapPadMwl;
		std::map<std::string,TH2*> fMapTH2Mwl;
		std::map<std::string,TGraphAsymmErrors*> fMapGraphErrorsMwl;
		std::map<std::string,std::vector<TArrow*> > fMapArrowsMwl;
		std::map<std::string,TF1*> fMapTF1Mwl;
		std::map<std::string,TPolyLine*> fMapButterflyMwl;

		TLegend *fLegendMwl;

		void DrawMultiWaveLength(MultiWaveLengthFactory &Mwl);

		void InitCanvasMwl(MultiWaveLengthFactory &Mwl);
		void InitPadMwl(MultiWaveLengthFactory &Mwl);
		void InitTH2Mwl(MultiWaveLengthFactory &Mwl);
		void InitGraphErrorsMwl(MultiWaveLengthFactory &Mwl);
		void InitTF1Mwl(MultiWaveLengthFactory &Mwl);
		void InitButterflyMwl(MultiWaveLengthFactory &Mwl);
		void InitLegendMwl(MultiWaveLengthFactory &Mwl);

		void CleanCanvasMwl();
		void CleanPadMwl();
		void CleanTH2Mwl();
		void CleanGraphErrorsMwl();
		void CleanTF1Mwl();
		void CleanButterflyMwl();
		void CleanArrowsMwl();

		std::string MakeKeyCanvasMwl(std::string mwlname);
		std::string MakeKeyPadMwl(std::string mwlname);
		std::string MakeKeyTH2Mwl(std::string mwlname);
		std::string MakeKeyGraphErrorsMwl(std::string mwlname,std::string graphname);
		std::string MakeKeyTF1Mwl(std::string mwlname,std::string tf1name);
		std::string MakeKeyArrowsMwl(std::string mwlname,std::string tf1name);
		std::string MakeKeyButterflyMwl(std::string mwlname,std::string butterflyname);


#ifdef CONSTRUCT_LIBRARY
	public:
		ClassDef(START::PlotFactory,1);
#endif

	};
}
#endif
