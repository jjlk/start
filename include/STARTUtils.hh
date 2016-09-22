#ifndef _STARTUTILS_
#define _STARTUTILS_

// STL
#include <vector>
#include <utility>
#include <map>


// ROOT
#include <TObject.h>
#include <TString.h>
#include <TGraph.h>

// START

namespace START {
	/**
	 * \brief Contains functions usefull in different classes
	 *
	 * Basically you put here a function in the case you don't know where to put it
	 *
	 * \author HAP-Fr team
	 */
	class STARTUtils : public TObject
	{
  
	public: 
		virtual ~STARTUtils() {};

		/**
		 * \brief Options for the IRFs.
		 * \param HapFr H.E.S.S. stereo mode
		 * \param HapFr H.E.S.S. hybrid mode
		 * \param HapFr H.E.S.S. mono mode
		 */
		typedef enum {HapFr_Hess_Stereo, HapFr_Hess_Hybrid, HapFr_Hess_Mono} IrfOpt;

		/**
		 * \brief Options for the handling of bands during the minimization
		 * \param None: each band is treated individually
		 * \param None: the Bands are merged in one band (no interpolation between IRF is done)
		 */	  
	        typedef enum {NoGrouping, StackInOneBand, StackInZenBands} BandsGroupingTypeForMinimization;
		
		static double LiMaSignificance(int non, int noff, double alpha); // return Li-Ma significance

		static double GetExcess(int non, int noff, double alpha); // return excess

		static double GetExpectedOff(int non, int noff, double alpha, double ns); // return expected Off

		static double GetExpectedOn(int non, int noff, double alpha, double ns); // return expected On

		static double GetNormalizationConstante(); // return constante used during the minimization to avoid precisions issues

		static double GetTeVToErg(); // return convertion factor from TeV to ergs
		static double GetTeVToMeV(); // return convertion factor from TeV to MeV
		static double GetTeVToGeV(); // return convertion factor from TeV to GeV
		static double GetTeVTokeV(); // return convertion factor from TeV to keV
		static double GetTeVToHz(); // return convertion factor from TeV to Hz

		static TString GetRunROOTFileName(TString config, unsigned int run, unsigned int telcode); // return root's run file name

		static TString GetRunInfoTreeName() {return "RunInfoTree";}
		static TString GetOnEventsTreeName() {return "EventsTree_BgMakerOn";}
		static TString GetOffEventsTreeName() {return "EventsTree_BgMakerOff";}
		static TString GetAllEventsTreeName() {return "EventsTree_AllSelected";}

		static double GetMJDFromSashSeconds(double seconds); // return MJD time from SASH seconds
		static double GetSashSecondsFromMJD(double mjd); // return SASH seconds from MJD
		static double GetLightCurveMJDReferenceTime() {return 51544.;}; // return 01/01/2000 in MJD

		static int GetIntegerExposantFromDecimalNumber(double number);

		static std::vector<Int_t> GetTelescopeListFromPattern(UInt_t TelPattern);
    
		static std::map<Int_t,TGraph*> GetContourFromSTARTFile(std::string filename, std::string canname);

		static std::vector<std::pair<double,double> > GetNewMoonMJDPeriods();


#ifdef CONSTRUCT_LIBRARY
	public:
		ClassDef(START::STARTUtils,1);
#endif
  
	};
}  
#endif
