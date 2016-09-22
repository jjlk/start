#ifndef _BAND_
#define _BAND_

// STL
#include <iostream>
#include <vector>
#include <map>
#include <sstream>

// ROOT
#include <TObject.h>
#include <Math/InterpolationTypes.h>
#include <Math/Interpolator.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TDirectory.h>

// START
#include "EnergyBin.hh"

namespace START {
	/**
	 * \brief Store EnergyBin in a vector and contain global informations
	 * about the band : zenith angle, offset, efficiency, time, etc...
	 *
	 * \warning The code expect that the first vector in fDistributionInstrumentMapTable which represent the range of ln(ereco/etrue) is the same range for each true energy
	 *
	 * \todo : Maybe create a function to initialize the interpolated vector (fVectorInterArea, fVectorInterResolution, etc...)
	 * \author HAP-Fr team
	 */

	class Band : public TNamed
	{
  
	public:
  
		Band(double zenon=0.0, double offset=0.0, double efficiency=0.0, double alpha=0.0, 
			 double livetime=0.0, unsigned int telcode=0);
		~Band(); 
		Band(Band const &BandCopy); // Copy constructor
		START::Band &operator=(START::Band const &BandCopy); // assignment operator

		virtual void Print(Option_t *option="") const; // *MENU*
		void PrintBand(std::ostream &os=std::cout) const;

		std::vector<EnergyBin> ebin; // we keep it public

		inline unsigned int   GetTelCode() const { return fTelCode;}
		inline double GetZenON()   const  { return fzenon;}
		inline double GetZenOFF()  const  { return fzenoff;}
		inline double GetOffset()  const  { return foffset;}
		inline double GetAzimuth() const {return fazimuth;}
		inline double GetEff()     const { return fefficiency;}
		inline double GetAlphaRun()const { return falpharun;}
		inline double GetLiveTime()  const { return flivetime;}
		inline unsigned int GetNbRun() const{ return fnbrun;}
		inline unsigned int GetKeepBand() const {return fkeepband;}
		inline double GetRunStartTime() const {return frunstarttime;}
		inline double GetRunEndTime() const {return frunendtime;}
		inline double GetLiveTimeFraction() const {return flivetimefraction;};

		// inline std::vector<double> GetVectorEnergy() const {return fVectorEnergy;}
		// inline std::vector<double> GetVectorArea() const {return fVectorArea;}
		// inline std::vector<double> GetVectorResolution() const {return fVectorResolution;}
		// inline std::vector<double> GetVectorBiais() const {return fVectorBiais;}
		// inline std::vector<double> GetVectorInterEnergy() const {return fVectorInterEnergy;}
		// inline std::vector<double> GetVectorInterArea() const {return fVectorInterArea;}
		// inline std::vector<double> GetVectorInterResolution() const {return fVectorInterResolution;}
		// inline std::vector<double> GetVectorInterBiais() const {return ffVectorInterBiais;}

		inline std::vector<double> &GetVectorEnergy()  {return fVectorEnergy;}
		inline std::vector<double> &GetVectorArea()  {return fVectorArea;}
		inline std::vector<double> &GetVectorResolution()  {return fVectorResolution;}
		inline std::vector<double> &GetVectorBiais()  {return fVectorBiais;}
		inline std::vector<double> &GetVectorInterEnergy()  {return fVectorInterEnergy;}
		inline std::vector<double> &GetVectorInterArea()  {return fVectorInterArea;}
		inline std::vector<double> &GetVectorInterResolution()  {return fVectorInterResolution;}
		inline std::vector<double> &GetVectorInterBiais()  {return fVectorInterBiais;}


		inline double GetEthMC() const {return fEthMC;}
		inline int GetFirstEmcBinNumber() const {return fFirstEmcBinNumber;}
		double GetFirstEmcVal() const 
		{if (fUseInstrumentEnergyDistribution) {return fFirstEmcValDistrib;} else {return fFirstEmcValFit;}} // VIM
  
		void SetTelCode(unsigned int TelCode)  { fTelCode = TelCode;}
		void SetZenON(double zenon)    { fzenon = zenon;}
		void SetZenOFF(double zenoff)   { fzenoff = zenoff;}
		void SetOffset(double offset)   { foffset = offset;}
		void SetAzimuth(double azimuth)  { fazimuth = azimuth;}
		void SetEff(double efficiency)      { fefficiency = efficiency;}
		void SetAlphaRun(double alpharun) { falpharun = alpharun;}
		void SetLiveTime(double livetime)  { flivetime = livetime;}
		void SetNbRun(unsigned int nbrun) { fnbrun = nbrun;}
		void SetKeepBand(unsigned int keepband) { fkeepband = keepband; }
		void SetRunStartTime(double runstarttime) { frunstarttime = runstarttime; }
		void SetRunEndTime(double runendtime) { frunendtime = runendtime; }
		void SetLiveTimeFraction(double livetimefraction) {flivetimefraction=livetimefraction;};
  
		void SetVectorEnergy(std::vector<double> VectorEnergy) {fVectorEnergy=VectorEnergy;}
		void SetVectorArea(std::vector<double> VectorArea) {fVectorArea=VectorArea;}
		void SetVectorResolution(std::vector<double> VectorResolution) {fVectorResolution=VectorResolution;}
		void SetVectorBiais(std::vector<double> VectorBiais) {fVectorBiais=VectorBiais;}
		void SetVectorInterEnergy(std::vector<double> VectorInterEnergy) {fVectorInterEnergy=VectorInterEnergy;}
		void SetVectorInterArea(std::vector<double> VectorInterArea) {fVectorInterArea=VectorInterArea;}
		void SetVectorInterResolution(std::vector<double> VectorInterResolution) {fVectorInterResolution=VectorInterResolution;}
		void SetVectorInterBiais(std::vector<double> VectorInterBiais) {fVectorInterBiais=VectorInterBiais;}

		void SetEthMC(double EthMC) {fEthMC = EthMC;}
		void SetFirstEmcBinNumber(int FirstEmcBinNumber) {fFirstEmcBinNumber = FirstEmcBinNumber;}

		// Return the area, biais or resolution for energy E
		double GetInterpolatedArea(double E);
		double GetInterpolatedBiais(double E);
		double GetInterpolatedResolution(double E);

		// Set the interpolators for area, biais and resolution
		void SetGSLInterpolatorForArea(std::vector<double> &x, std::vector<double> &y);
		void SetGSLInterpolatorForBiais(std::vector<double> &x, std::vector<double> &y);
		void SetGSLInterpolatorForResolution(std::vector<double> &x, std::vector<double> &y);
  
		// VIM :
		void SetUseOfInstrumentEnergyDistribution(bool usedistrib=true) {fUseInstrumentEnergyDistribution = usedistrib;}
		bool GetUseOfInstrumentEnergyDistribution() {return fUseInstrumentEnergyDistribution;}
  
		void SetFirstEmcValDistrib(float FirstEmcVal) {fFirstEmcValDistrib = FirstEmcVal;}
		void SetFirstEmcValFit(float FirstEmcVal) {fFirstEmcValFit = FirstEmcVal;}
		double GetFirstEmcValDistrib() {return fFirstEmcValDistrib;}
		double GetFirstEmcValFit() {return fFirstEmcValFit;}
  
  
		std::map< double, std::pair< std::vector<double>, std::vector<double> > > GetDistributionVectorTable() {return fDistributionInstrumentMapTable;}
		std::map< double, ROOT::Math::Interpolator* > GetDistributionInterpTable() {return fDistributionInstrumentInterpTable;}
  
		void SetDistributionVectorTable(  std::map<double, std::pair< std::vector<double>, std::vector<double> > > &distribtable) {
			fDistributionInstrumentMapTable = distribtable;
		}
  
		void InitDistributionInterpTable();
		void DeleteInstrumentInterpTable();
		double GetEResolProbabilityFromDistribution(double Ereco, double Etrue);
		double GetEResolProbabilityFromFit(double Ereco, double Etrue);
		double GetEResolProbability(double Ereco, double Etrue);

		void ClearBandInfo();
		void AddInfoFromBands(const std::vector<START::Band> &VecBand, bool IsUseEthresh=false);
		void InitEbinAlpha(double alpharun);
		void InitEbinLiveTime(double livetime);
  
		double GetNOnTot(bool IsUseEthresh=true) const;
		double GetNOffTot(bool IsUseEthresh=true) const;
		double GetExcessTot(bool IsUseEthresh=true) const;
		double GetSignificanceTot(bool IsUseEthresh=true) const;
		double GetNSthTot(bool IsUseEthresh=true) const;
		double GetNOffFitted(bool IsUseEthresh=true) const;
  
		void DisplayArea(Bool_t save=false); // *MENU*
		void DisplayEnergyResolutionForTrueEnergy(Double_t Energy_True); // *MENU*
		void DisplayEnergyResolutionForRecoEnergy(Double_t Energy_Reco); // *MENU*
		void DisplayEnergyResolutionForTrueOrRecoEnergy(Double_t Energy, Bool_t IsTrueEnergy);
		void DisplayEnergyResolution(Bool_t save=false); // *MENU*
    
		// Simple function to have a standard init of the bins
		void InitEnergyBins(Double_t emin, Double_t emax, Int_t nbins);
		EnergyBin* GetEnergyBin(Double_t ener);

		std::string ConstructBandName() const;
		std::string ConstructBandTitle() const;
		void UpdateNameAndTitle();
		void InitInterpolator(); // *MENU*
  
	private:
		std::string ConstructBandNameOrTitle(Bool_t constructname) const;

		unsigned int fnbrun; ///< number of the run
		unsigned int fTelCode; ///< telcode of the run
		double fzenon; ///< ON zenith of the run
		double fzenoff; ///< OFF zenith of the run
		double foffset; ///< offset of the run
		double fazimuth; ///< mean azimuth of the run
		double falpharun; ///< normalization ON/OFF region of the run
		double fefficiency; ///< efficiency of the run
		double flivetime; ///< livetime of the run
		double flivetimefraction; ///< time efficiency of the run
		unsigned int fkeepband; ///< 1 the band is kept, 0 if not
		double frunstarttime; ///< time of first event in MJD
		double frunendtime; ///< time of last event in MJD

		std::vector<double> fVectorEnergy;
		std::vector<double> fVectorArea;
		std::vector<double> fVectorResolution;
		std::vector<double> fVectorBiais;
		std::vector<double> fVectorInterEnergy; 
		std::vector<double> fVectorInterArea; 
		std::vector<double> fVectorInterResolution;
		std::vector<double> fVectorInterBiais;
		std::map< double, std::pair< std::vector<double>, std::vector<double> > > fDistributionInstrumentMapTable; ///< Map that contain the  RMF. Map Key : True Energy, the  value is a pair of vector : first vector of ln(Ereco/Etrue) and second, the value associated. Through the code, it is expected that the vector have the  same size and range 
		std::map< double, ROOT::Math::Interpolator* > fDistributionInstrumentInterpTable; //!
  
		double fEthMC; ///< energy threshold
		int fFirstEmcBinNumber; ///< bin number in vector MonteCarlo::GetEnergy() for which 
		///< area, resol and bias are available.
		///< In the estimation of the expected number of gammas, the
		///< integral on true energy must avoid going below this energy.
			  
		double fFirstEmcValFit; ///< corresponding MC energy value (point like)
		double fFirstEmcValDistrib; ///< corresponding MC energy value (etends)
  
		bool fUseInstrumentEnergyDistribution; ///< In order to tell to the program to use the distribution
  
		ROOT::Math::Interpolation::Type fInterpolTypeForBand; ///< Interpolation type
  
		// Interpolator of area, resolution and biais
		ROOT::Math::Interpolator *fInterpolArea; //!  ///< Area interpolator
		ROOT::Math::Interpolator *fInterpolBiais; //!  ///< Biais interpolator
		ROOT::Math::Interpolator *fInterpolResol; //! ///< Resolution interpolator

		// JLK Function to get a Band stack from a vector of Band
		// from VIM implementation (grouping of bands)
		// IRF won't be available
		//Band StackBandFromVectorOfBands();
  
#ifdef CONSTRUCT_LIBRARY
	public:
		ClassDef(START::Band,1);
#endif
  

	};

}

#endif
