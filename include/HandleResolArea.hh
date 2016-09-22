#ifndef _STORERESOLAREA_
#define _STORERESOLAREA_

// STL
#include <iostream>
#include <map>
#include <vector>
#include <utility>

// ROOT
#include <TString.h>
#include <TObject.h>
#include <TMultiGraph.h>
#include <TCanvas.h>

// START
#include "MonteCarlo.hh"
#include "Band.hh"
#include "STARTUtils.hh"
#include "Config.hh"

namespace START {

	/**
	 * \brief Stock in bands effective area, resolutions and biais stored
	 * in IRF files. Additional infos are then stock in Band and EnergyBin.
	 * 
	 * Histograms are stored in maps and their keys are defined by efficacity, 
	 * zenith angle and offset, then interpolated areas, biais and resolution
	 * are stored in bands to speed up the computations of the expected excess.
	 *
	 * The energy threshold is computed and then used as a flag for each bin
	 * in energy to determine if the bins can be used in the minimization.
	 *
	 * \author HAP-Fr team
	 */
	class HandleResolArea : public TObject
	{
  
	public :

		typedef enum {AreaAndBiaisResolMethod=0, AreaMaxMethod} ESafeThresholdMethod;

		/**
		 * \brief Criteria to accept an EnergyBin using the mc energy threshold
		 * \param Safe the first EnergyBin is accepted if the emin is greater than EthMC
		 * \param Explorer the first EnergyBin is accepted if Emin is closer to EthMC than Emax
		 * \param Insane the EnergyBin is accepted in any case
		 */
		typedef enum {Safe,Explorer,Insane} EnergyBinThresholdCondition;
  
		// Constructors :
		HandleResolArea() {}; //Needed by ROOT
		//HandleResolArea(TString configname, TString configtype, TString mcprod, bool forceuseofdistribution=false, bool readnow=false);
		//HandleResolArea(TString configname, TString mcprod, bool forceuseofdistribution=false, bool readnow=false);
  		HandleResolArea(const START::Config &Conf,
						START::STARTUtils::IrfOpt Type,
						bool readnow=false); //
		virtual ~HandleResolArea();
  
		// set energy threshold condition
		void SetESafeThresholdCondition(ESafeThresholdMethod emeth, double areamin_or_areacondition, double resolfactor=2.);
		void SetEnergyBinThresholdCondition(EnergyBinThresholdCondition Condition) {fEnergyBinThresholdCondition=Condition;};

		// copy in bands infos
		void SetAndUpdateBands(std::vector<Band> &BandArray);

		// set use of distribution
		void SetForceUseDistribution(bool fdistrib=true) { fForceUseDistribution = fdistrib;};

		std::map<TString, TH1F *> GetMapHistoResol() const {return fmapHistoResol;}; ///< Get histograms resolution
		std::map<TString, TH1F *> GetMapHistoBiais() const {return fmapHistoBiais;}; ///< Get histograms biais
		std::map<TString, TH1F *> GetMapHistoArea() const {return fmapHistoArea;}  ///< Get histograms area
		std::map<TString, std::pair<std::vector<double >, std::vector<double > > > GetMapVectorsResol() const {return fmapVectorsResol;} ///<Get vectors energy resolution
		std::map<TString, std::pair<std::vector<double >, std::vector<double > > > GetMapVectorsBiais() const {return fmapVectorsBiais;} ///<Get vectors energy biais
		std::map<TString, std::pair<std::vector<double >, std::vector<double > > > GetMapVectorsArea() const {return fmapVectorsArea;} ///<Get vectors energy area

		std::map<TString,std::pair<double,double> > GetPairResolMeanSigma() const {return fpairResolMeanSigma;} //JLK je sais pas a quoi ça sert ça va degager
		std::map<TString,double> GetMapAreaValues() const {return fmapAreaValues;} //JLK je sais pas a quoi ça sert ça va degager
		std::map<TString,double> GetMapResolValues() const {return fmapResolValues;} //JLK je sais pas a quoi ça sert ça va degager
		std::map<TString,double> GetMapBiaisValues() const {return fmapBiaisValues;} //JLK je sais pas a quoi ça sert ça va degager

		// Make keys for the maps with eff, zen, and offset
		TString MakeKeyWeight(double eff, double zen, double off);

		// Make keys for the maps with eff, zen, offset and energy
		TString MakeKeyWeightEnergy(double eff, double zen, double off, double ene);
  
		// Make keys for the maps with eff, zen, offset and tel
		TString MakeKeyRun(double eff, double zen, double off, TString tel);
  
		// Make keys for the maps with eff, zen, offset, energy and tel
		TString MakeKeyRunEnergy(double eff, double zen, double off, double en, TString tel);
  
		// Make a string which contain the MC type
		TString MakeStringTelCode(int telcode);

		void DrawIRF(double eff, double zen, double off, unsigned int telcode);
		void DrawAllIRF();

		void DrawIRFBandDiagnostic(Band &DrawBand);

		void SetMcProductionName(TString mcname) {fMcProductionName=mcname;};


		//Ces fonctions doivent etre public pour que je puisse y avoir acces des reprojband
		// Copy Interpolators in bands   
		int CopyInterpolatorsInBands(std::vector<Band> &InputBandArray);
		// Copy Interpolators in bins    
		int CopyInterpolatorsInBins(std::vector<Band> &InputBandArray);
		// Determine value of the first MC energy bin with valid area, resol and bias
		void SetInBandsFirstEmcBin(std::vector<Band> &InputBandArray); 
		// Copy MC's energy threshold in bands
		void SetInBandsMCEnergyThreshold(std::vector<Band> &InputBandArray);
		// Set flag in bins if they are used
		void SetInBinsKeepBin(std::vector<Band> &InputBandArray);
		// Set effective area in bins
		void SetInBinsEffectiveArea(std::vector<Band> &InputBandArray);
  

	private :

		// Copy and assignement are a priori not necessary thus not allowed
		HandleResolArea(HandleResolArea const &StoreCopy); // Copy constructor
		HandleResolArea &operator=(HandleResolArea const& StoreCopy);

		// read and store histos in maps
		int ReadAndStore();

		// Build two vectors x and y whith the values contained in the histogram
		// stored in the map "Map" with the key "Key"
		void MakeVectorsFromMap(std::map<TString, TH1F *> Map, TString Key,
								std::vector<double> &x, std::vector<double> &y) const;
  
		// Make vector of area or resolution for a band with a multidimensional interpolation
		// x : energy and y : resolution, biais or area
		void MakeVectorForBand(Band myBand,
							   std::map<TString,std::pair<std::vector<double >,std::vector<double > > > mapVector,
							   std::vector<double> &x,std::vector<double> &y);
  
		// Make Interpolated function for area, resolution and biais (used in the likelihood)
		void MakeSmoothFunction(const std::vector<double> x, const std::vector<double> y,
								std::vector<double> &xinter, std::vector<double> &yinter);
  
		// Copy Vectors of area or resolution in bands
		int CopyVectorsInBands();
  
		// Copy Interpolators in bands
		int CopyInterpolatorsInBands();

		// Copy Interpolators in bins
		int CopyInterpolatorsInBins();

		// Determine value of the first MC energy bin with valid area, resol and bias
		void SetInBandsFirstEmcBin();
  
		//double FindZeroInBandArray(std::vector<double> x, std::vector<double> y, bool &zeroexist);
		static double FindZeroInBandArray(std::vector<double> x, std::vector<double> y, int &zeroexist);

		// Copy MC's energy threshold in bands

		void SetInBandsMCEnergyThreshold();
		static void ComputeSafeThresholdFromAreaAndResolFactor(Band &InputBand, double areamin, double resolfactor);
		static void ComputeSafeThresholdFromAreaMaxFraction(Band &InputBand, double areamaxfraction);

 
		// Set flag in bins if they are used
		void SetInBinsKeepBin();
    
		// Set effective area in bins
		void SetInBinsEffectiveArea();
    
		void ClearTablesVector(); // VIM : Function to clear the vector fill in readandstore
		void FillBandsWithDistributionBoolean(); // VIM : Function to say to the band vector that the distribution should be use instead of sigma and mean.
		void MakeWeightedInstrumentDistributionForBand(Band myBand,std::map<double,std::pair< std::vector<double>,std::vector<double> > > &DistributionInstrumentMapTable);  

		// allows to determine which are the needed MC ranges to do the interpolations
		void GetNeededMcRange(std::vector<Band> const &BandArray,
							  double &mineff, double &maxeff,
							  double &minoff, double &maxoff,
							  double &minzen, double &maxzen);
		
	private :

		double FindEffMinMC(Band myBand);
		double FindEffMaxMC(Band myBand);
		double FindOffMinMC(Band myBand);
		double FindOffMaxMC(Band myBand);
		double FindZenMinMC(Band myBand);
		double FindZenMaxMC(Band myBand);

		TString fMcProductionName;

		std::vector<Band> *fBandArray; //! Band array

		bool fForceUseDistribution; // VIM : This variable can be usefull for distinguish between fullcontainment and plike hypo !
		bool fIsDataSummaryFromDynamicAlloc;
		bool fIsInstrumentTableAlreadyLoad; //VIM : This variable is in case you want to load the table asap in the constructor


		TString fConfigName; ///< name of the configuration
		TString fConfigType; ///< type of the configuration

		//Min & max values of MC values to look for
		double fEffMin, fEffMax;
		double fOffMin, fOffMax;
		double fZenMin, fZenMax;

		// Map of histos for different off, zen and eff
		// containing the associated histogramms of resolution,
		// biais and effective area 
		std::map<TString, TH1F *> fmapHistoResol;
		std::map<TString, TH1F *> fmapHistoBiais;
		std::map<TString, TH1F *> fmapHistoArea; 
		std::map<TString, TH1F *> fmapHistoDistrib; // VIM For the full distribution of R(E,E')
  
		// Map of pairs which contain the mean and the sigma of
		// the resolution
		std::map<TString,std::pair<double,double> > fpairResolMeanSigma;

		// Maps which contains area,resol and biais values
		std::map<TString,double> fmapAreaValues;
		std::map<TString,double> fmapBiaisValues;
		std::map<TString,double> fmapResolValues;

		// Maps which contains vectors energy and values of histo
		std::map<TString, std::pair<std::vector<double >, std::vector<double > > > fmapVectorsResol;
		std::map<TString, std::pair<std::vector<double >, std::vector<double > > > fmapVectorsBiais;  
		std::map<TString, std::pair<std::vector<double >, std::vector<double > > > fmapVectorsArea;



		// For the spectrum method for safethreshold 
		double fAreaMin; ///< AreaMin used for the condition of the safe threshold, en Hectares.
		double fResolFactor; ///< For the safethreshol, the Biais should be inferior to fResolFactor * Resol ! 

		ESafeThresholdMethod fSafeThresholdMethod;
		EnergyBinThresholdCondition fEnergyBinThresholdCondition;

		// For the MaxArea Method
		double fFractionMaxArea; ///< Fraction of the Maximum area to cut to define the safe threshold 
  
		std::map<TString,TCanvas*> fMapIRFBiaisCanvas;
		std::map<TString,TCanvas*> fMapIRFResolutionCanvas;
		std::map<TString,TCanvas*> fMapIRFEffectiveAreaCanvas;

		std::map<TString,TGraph*> fMapIRFBiaisGraph;
		std::map<TString,TGraph*> fMapIRFResolutionGraph;
		std::map<TString,TGraph*> fMapIRFEffectiveAreaGraph;

		std::map<TString,TCanvas*> fMapBandIRFEffectiveAreaCanvas;
		std::map<TString,TCanvas*> fMapBandIRResolutionCanvas;
		std::map<TString,TCanvas*> fMapBandIRFBiaisCanvas;

		std::map<TString,TMultiGraph*> fMapBandIRFBiaisMultiGraph;
		std::map<TString,TMultiGraph*> fMapBandIRFEffectiveAreaMultiGraph;
		std::map<TString,TMultiGraph*> fMapBandIRFResolutionMultiGraph;


		TString MakeKeyIRFBiaisCanvas(double eff, double zen, double off, unsigned int telcode);
		TString MakeKeyIRFResolutionCanvas(double eff, double zen, double off, unsigned int telcode);
		TString MakeKeyIRFEffectiveAreaCanvas(double eff, double zen, double off, unsigned int telcode);

		TString MakeKeyIRFBiaisGraph(double eff, double zen, double off, unsigned int telcode);
		TString MakeKeyIRFResolutionGraph(double eff, double zen, double off, unsigned int telcode);
		TString MakeKeyIRFEffectiveAreaGraph(double eff, double zen, double off, unsigned int telcode);

		TString MakeKeyIRFBandEffectiveAreaCanvas(Band &BandToDraw);
		TString MakeKeyIRFBandResolutionCanvas(Band &BandToDraw);
		TString MakeKeyIRFBandBiaisCanvas(Band &BandToDraw);

		TString MakeKeyIRFBandBiaisMultiGraph(Band &BandToDraw);
		TString MakeKeyIRFBandResolutionMultiGraph(Band &BandToDraw);
		TString MakeKeyIRFBandEffectiveAreaMultiGraph(Band &BandToDraw);

		void BuildIRFBiaisCanvas(double eff, double zen, double off, unsigned int telcode);
		void BuildIRFResolutionCanvas(double eff, double zen, double off, unsigned int telcode);
		void BuildIRFEffectiveAreaCanvas(double eff, double zen, double off, unsigned int telcode);

		void BuildBandEffectiveAreaCanvas(Band &BandToDraw);
		void BuildBandResolutionCanvas(Band &BandToDraw);
		void BuildBandBiaisCanvas(Band &BandToDraw);

		TString MakeBandHeader(Band &BandToDraw);

		bool buildirfindividualscanvas; ///< build effective area, resolution and biais canvas

		START::MonteCarlo *fMc;
		
#ifdef CONSTRUCT_LIBRARY
	public:
		ClassDef(START::HandleResolArea,1);
#endif
	};
}
#endif
