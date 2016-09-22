#ifndef _BANDSFACTORY_
#define _BANDSFACTORY_

// STL
#include <iostream>
#include <vector>

// ROOT
#include <TObject.h>
#include <TString.h>

// START
#include "HandleResolArea.hh"
#include "Config.hh"
#include "Band.hh"
#include "STARTUtils.hh"

namespace START {

	/**
	 * \brief Stock data in a vector of Band
	 *
	 * \author HAP-Fr team
	 */
	class BandsFactory : public TObject
	{

	public:

		BandsFactory() {}; //Needed by ROOT
		BandsFactory(const Config &config,
					 START::STARTUtils::IrfOpt Type,
					 bool verbose=false);
		~BandsFactory();

		int MakeBands(std::vector<Band> &BandArray);

		void PrintBands(const std::vector<Band> &BandArray) const;
		void PrintSummaryBand(const std::vector<Band> &BandArray) const;
		void FillBandInfo(std::vector<Band> &BandArray, std::map<TString,std::vector<int> >  &InfoReprojArray,bool LowStat,bool VeryLowStat,Config myConfig) const;
		//void ReprojBands(std::vector<Band> &BandArray, std::map<TString, std::vector<int> >  &InfoReprojArray, std::vector<Band> &ReprojArray,TString configname, HandleResolArea &HandleIRF);
		//void RebinEnergy(const std::vector<Band> &BandArray, std::vector<Band> &BandRebinArray, double sigrebin,double MinE, HandleResolArea &HandleIRF); 
		void SetFitEnergyRangeFlagInEnergyBins(std::vector<Band> &BandArray) const;
		void CheckAndSetIfBandsAreInMCLimits(std::vector<Band> &BandArray) const;
		bool CheckForSelectedBands(const std::vector<Band> &BandArray) const;

		bool CheckIfEventInTimeWindow(double currenttime);
		bool CheckIfEntireRunIncludedInSingleTimeWindow(double TimeFirstEvt,double TimeLastEvt);
		void Change_N_ON(std::vector<Band> &BandArray, double n_on) const;
  
	private:

		//Band SummaryBand;

		Config *fConfig;

		bool fverbose;

		std::vector<int> fRunList; // run list
		std::vector<unsigned int> fRunTelType; // telcode list 
		std::vector<std::string> fInputFileList; // input file list

		double fTimeToAdd;
		double fMiddleZen;
		int fBinNumber;
		double fBinWidth;
		double fEmin, fEmax;
		int fNbBands;

		START::STARTUtils::IrfOpt fIrfOpt; // JLK Add
		
#ifdef CONSTRUCT_LIBRARY
	public:
		ClassDef(START::BandsFactory,1);
#endif

	};
}
#endif
