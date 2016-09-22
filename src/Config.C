// STL
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>

// ROOT
#include <TString.h>

// START
#include "Config.hh"

// Utilities
#define DEBUG 0
#include "debugging.hh"
#include "StringTools.hh"

#define INFO std::cout << INFOCOLOR << "Config> " << RESETCOLOR
#define WARNING std::cout << WARNINGCOLOR << "Config> " << RESETCOLOR

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::Config)
#endif


/**
 * \brief Constructor for ROOT
 */
START::Config::Config() 
:fFileConfig(""),
	fFileList(""),
	fUserCodeVersion(""),
	fUserUseInputFileList(false),
	fUserRunListsFolderAddress(""),
	fUserRunListFileBaseName(""),
	fUserInstrumentResponseAddress(""),
	fUserRootFilesFolderAddress(""),
	fUserAnalysisConfig(""),
	fUserAnalysisConfigType(""),
	fUserForceAreaAnalysisConfig(""),
	fUserMcProductionName(""),
	fUserOutputFolderName(""),
	fUserMaxCosZenBinWidth(0.),
	fUserZenithMax(0.),
	fUserOffsetMax(0.),
	fUserERangeMin(0.025),
	fUserERangeMax(125),
	fUserERangeBinsNumber(30),
	fUserFitEMin(-1),
	fUserFitEMax(-1),
	fUserAreaMin(4.),
	fUserSkipHighThreshRuns(0.),
	fUserAutomaticFit(0),
	fUserIsMc(false),
	fUserTmin(0.),
	fUserVerbose(false),
	fUserEref(1.),
	fUserDrawDataSummaryPlots(false),
	fUserDrawPlots(true),
    fUserResidualsStyle(true),
	fUserPlotStyle(0),
	fUserFitResultsDrawing(false),
	fUserSourceName(""),
	fUserHypothesisComparison(false),
	fUserComputeContours(0),
	fUserContoursNumberOfPoints(30),
	fUserComputeScanLikelihood(0),
	fUserMinimizationLevel(0),
	fUserLightCurveTimeInterval(0),
	fUserLightCurveTimeAxisUnits(0),
    fUserIrfOpt("HapFr_Hess_Stereo"),
    fUserBandsGroupingType("NoGrouping"),
	fUserOutputRootFileName("StartOutput.root")
{

	fUserHypothesisIntegrationEnergyRange = std::make_pair(0.2,120.);
	fUserLightCurveTimeRange = std::make_pair(0.,0.);
	fUserLightCurveIntegratedFluxEnergyRange = std::make_pair(0.,0.);
	fUserLightCurveUserTimeIntervals.clear();
	fUserAddHypothesis.clear();
	fUserButterfly.clear();
	fUserHypothesisRebin.clear();
	fUserLightCurveTimeCuttingType.clear();
	fUserLightCurveErrorsHandling.clear();
}

/**
 * \brief Constructor
 *
 * \param config file (path/name)
 * \param readlistnow Parameter to require that the list is immediately read (false is usefull when you're doing other analyis than HAP-Fr)
 */
START::Config::Config(TString file_config, Bool_t readlistnow) 
	:fFileConfig(file_config),
	 fFileList(""),
	 fUserCodeVersion(""),
	 fUserUseInputFileList(false),
	 fUserRunListsFolderAddress(""),
	 fUserRunListFileBaseName(""),
	 fUserInstrumentResponseAddress(""),
	 fUserRootFilesFolderAddress(""),
	 fUserAnalysisConfig(""),
	 fUserAnalysisConfigType(""),
	 fUserForceAreaAnalysisConfig(""),
	 fUserMcProductionName(""),
	 fUserOutputFolderName(""),
	 fUserMaxCosZenBinWidth(0.),
	 fUserZenithMax(0.),
	 fUserOffsetMax(0.),
	 fUserERangeMin(0.),
	 fUserERangeMax(0.),
	 fUserERangeBinsNumber(0),
	 fUserFitEMin(-1),
	 fUserFitEMax(-1),
	 fUserAreaMin(4.),
	 fBandZenMax(0.),
	 fBandOffMax(0.),
	 fBandEffMax(0.),
	 fBandZenMin(0.),
	 fBandOffMin(0.),
	 fBandEffMin(0.),
	 fBandZenbin_UseReproj(0),
	 fBandOffbin_UseReproj(0),
	 fBandEffbin_UseReproj(0),
	 fBandZenbin_SmallStat(0),
	 fBandOffbin_SmallStat(0),
	 fBandEffbin_SmallStat(0),
	 fBandZenbin_VerySmallStat(0),
	 fBandOffbin_VerySmallStat(0),
	 fBandEffbin_VerySmallStat(0),
	 fUserSkipHighThreshRuns(0.),
	 fUserAutomaticFit(0),
	 fUserIsMc(false),
	 fUserTmin(0.),
	 fUserVerbose(false),
	 fUserEref(1.),
	 fUserDrawDataSummaryPlots(false),
	 fUserDrawPlots(true),
	 fUserResidualsStyle(0),
	 fUserPlotStyle(0),
	 fUserFitResultsDrawing(false),
	 fUserSourceName(""),
	 fUserHypothesisComparison(false),
	 fUserComputeContours(0),
	 fUserContoursNumberOfPoints(30),
	 fUserComputeScanLikelihood(0),
	 fUserMinimizationLevel(0),
	 fUserLightCurveTimeInterval(0),
	 fUserLightCurveTimeAxisUnits(0),
	 fUserOutputRootFileName("StartOutput.root")
{

	fUserHypothesisIntegrationEnergyRange = std::make_pair(0.2,120.);
	fUserLightCurveTimeRange = std::make_pair(0.,0.);
	fUserLightCurveIntegratedFluxEnergyRange = std::make_pair(0.,0.);
	fUserAddHypothesis.clear();
	fUserButterfly.clear();
	fUserHypothesisRebin.clear();
	fUserLightCurveUserTimeIntervals.clear();
	fUserLightCurveTimeCuttingType.clear();
	fUserLightCurveErrorsHandling.clear();

	if(ReadConfig(fFileConfig)==1) INFO << "Reading config... ok" << std::endl;
	else {
		WARNING << "Reading config " << fFileConfig << "... Failed" << std::endl;
		exit(EXIT_FAILURE);
	}

	if (readlistnow) {
		// We build the run's list
		TString listfilename = fUserRunListsFolderAddress+"/"+fUserRunListFileBaseName+".list";
    
		if(ReadList(listfilename)==1) { 
			INFO << "Reading list... ok" << std::endl;
		}
		else {
			WARNING << "Reading list : " << fFileList << "... Failed" << std::endl;
			exit(EXIT_FAILURE);
		}
    
	}

	PrintConfig();
	if (fUserUseInputFileList) { PrintInputFileList(); }
	else { PrintRunList(); }
}

/**
 * \brief Destructor
 */
START::Config::~Config()
{

}

/**
 * \brief Read the configuration file and set the members variables
 */
int START::Config::ReadConfig(TString file_config) {
  
	fFileConfig=file_config; 
  
	std::ifstream configfile(fFileConfig.Data());
	if (!configfile.is_open()){
		std::cout << "Problem opening config file " << fFileConfig.Data() << std::endl;
		return 0;
	}

	std::string input;

	while (1){     

		std::getline(configfile, input, '\n');
    
		if (input.empty()) { //skip empty lines
			if(configfile.eof()) break;
		}
    
		// if (input.find("###")==0) //found ### at col 0
		//   continue;

		// VIM : Add a mecanism to extract everything before the Sharp symbol that will be used as a Comment symbol
		// TODO : Need to add a mechanism to see if a line is only composed by white space
		unsigned long ksharp = input.find("#");
		input = input.substr(0,ksharp);
		if (input.empty()) {
			// At this step, this means that the full line is commented !
			continue;
		}

    
		// JLK : we skipped delimiters
		if(input.find("[User]")==0 || input.find("[Data]")==0 || input.find("[Analysis]")==0
		   || input.find("[Hypothesis]")==0 || input.find("[MC]")==0 || input.find("[Plots]")==0
		   || input.find("[LightCurve]")==0) {
			continue;
		}
      
		DEBUG_OUT << input << std::endl;
    
		//SP: il faut un peu ameliorer la methode pour prendre en compte la
		//possibilite que l'utilisateur mette des espaces apres sa variable (apres la 2eme colonne)
    
		if (input.find("UserUseInputFileList") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserUseInputFileList=(bool)atoi(input.c_str());
		}

		//deuxieme possibilite pour faute d'orthographe
		if ( (input.find("UserRunListsFolderAddress") != std::string::npos) || 
			 (input.find("UserRunListsFolderAdress") != std::string::npos) ) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserRunListsFolderAddress=input;
		}

		if (input.find("UserRunListFileBaseName") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserRunListFileBaseName=input;
		}
    
		if (input.find("UserAnalysisConfig") != std::string::npos
			&& input.find("UserAnalysisConfigType") == std::string::npos) { //MPA: otherwise it gets overwritten by this option!
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserAnalysisConfig=input;
		}

		if (input.find("UserAnalysisConfigType") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserAnalysisConfigType=input;
		}

		//ADA add forcedConfig
		if (input.find("UserForceAreaAnalysisConfig") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserForceAreaAnalysisConfig=input;
		}

		if (input.find("UserRootFilesFolderAddress") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserRootFilesFolderAddress=input;
		}

		if (input.find("UserOutputFolderName") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserOutputFolderName=input;
		}

		// name of the source
		if (input.find("UserSourceName") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserSourceName=input;
		}

		// name of the output root file
		if (input.find("UserOutputRootFileName") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserOutputRootFileName=input;
		}

		// JLK add
		if (input.find("UserIrfOpt") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserIrfOpt=input;
		}

		// JLK add
		if (input.find("UserBandsGroupingTypeForMinimization") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserBandsGroupingType=input;
		}

		
		if (input.find("UserMcProductionName") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserMcProductionName=input;
		}

		if (input.find("UserIsMc") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserIsMc=(Bool_t)atoi(input.c_str());
		}

		if (input.find("UserTmin") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserTmin=atof(input.c_str());
			if(fUserTmin<=0.) {
				WARNING << "UserTmin (in seconds) should be positive and at least > 0 !!! ==> EXIT" << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		// Zenith binning
		if (input.find("UserMaxCosZenBinWidth") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserMaxCosZenBinWidth=atof(input.c_str());
			if(fUserMaxCosZenBinWidth<=0.) {
				WARNING << "UserMaxCosZenBinWidth should be positive and at least > 0 !!! ==> EXIT" << std::endl;
				INFO << "UserMaxCosZenBinWidth = " << fUserMaxCosZenBinWidth << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		// MJD window
		if (input.find("UserMJDwindow") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			int strsize=input.size();
			int coma=input.find(",");
			fUserMJDwindow.push_back(std::make_pair<double,double>(atof(input.substr(0,coma).c_str()),atof(input.substr(coma+1,strsize-coma).c_str())));
		}

		// Integration energy range for integrated hypothesis
		if (input.find("UserHypothesisIntegrationEnergyRange") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			int strsize=input.size();
			int coma=input.find(",");
			fUserHypothesisIntegrationEnergyRange = std::make_pair(atof(input.substr(0,coma).c_str()),atof(input.substr(coma+1,strsize-coma).c_str()));
		}

		// Maximal offset
		if (input.find("UserOffsetMax") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserOffsetMax=atof(input.c_str());
			if(fUserOffsetMax<=0.) {
				WARNING << "UserOffsetMax (in degrees) should be positive and at least > 0 !!! ==> EXIT" << std::endl;
				INFO << "UserOffsetMax = " << fUserOffsetMax << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		// Maximal zenith angle
		if (input.find("UserZenithMax") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserZenithMax=atof(input.c_str());
			if(fUserZenithMax<=0.) {
				WARNING << "UserZenithMax (in degrees) should be positive and at least > 0 !!! ==> EXIT" << std::endl;
				INFO << "UserZenithMax = " << fUserZenithMax << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		// Energy range
		if (input.find("UserERangeMin") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserERangeMin=atof(input.c_str());
			if(fUserERangeMin<=0.) {
				WARNING << "UserERangeMin (in TeV) should be positive and at least > 0 !!! ==> EXIT" << std::endl;
				INFO << "UserERangeMin = " << fUserERangeMin << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		// Energy range
		if (input.find("UserERangeMax") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserERangeMax=atof(input.c_str());
			if(fUserERangeMax<=0. || fUserERangeMax<fUserERangeMin) {
				WARNING << "UserERangeMax (in TeV) should be positive and greater than UserERangeMin !!! ==> EXIT" << std::endl;
				INFO << "UserERangeMin = " << fUserERangeMin << std::endl;
				INFO << "UserERangeMax = " << fUserERangeMax << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		// Energy range
		if (input.find("UserERangeBinsNumber") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserERangeBinsNumber=atoi(input.c_str());
			if(fUserERangeBinsNumber<=0) {
				WARNING << "UserERangeBinsNumber should be positive and at least > 0 !!! ==> EXIT" << std::endl;
				INFO << "UserERangeBinsNumber = " << fUserERangeBinsNumber << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		//Selection of energy range for spectrum fit
		if (input.find("UserFitEMin") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserFitEMin=atof(input.c_str());
			if(fUserFitEMin<=0) {
				WARNING << "UserFitEMin (in TeV) should be positive and at least > 0 !!! ==> EXIT" << std::endl;
				INFO << "UserFitEMin = " << fUserFitEMin << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		//Selection of energy range for spectrum fit
		if (input.find("UserFitEMax") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserFitEMax=atof(input.c_str());
			if(fUserFitEMax<=0. || fUserFitEMax<fUserFitEMin) {
				WARNING << "UserFitEMax (in TeV) should be positive and greater than UserFitEMax !!! ==> EXIT" << std::endl;
				INFO << "UserFitEMin = " << fUserFitEMin << std::endl;
				INFO << "UserFitEMax = " << fUserFitEMax << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		// Minimum Area for safe threshhold
		if (input.find("UserAreaMin") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserAreaMin=atof(input.c_str());
			if(fUserAreaMin<=0) {
				WARNING << "UserAreaMin (in hectare) should be positive and at least > 0!!! ==> EXIT" << std::endl;
				INFO << "UserAreaMin = " << fUserAreaMin << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		//LJ
		// Maximum zenith band value for the reprojection
		if (input.find("BandZenMax") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fBandZenMax=atof(input.c_str());
			if(fBandZenMax<0) {
				WARNING << "BandZenMax should be positive!!! ==> EXIT" << std::endl;
				INFO << "BandZenMax = " << fBandZenMax << std::endl;
				exit(EXIT_FAILURE);
			}
			if(fBandZenMax>70) {
				WARNING << "BandZenMax should be inferior to 70degrees!!! ==> EXIT" << std::endl;
				INFO << "BandZenMax = " << fBandZenMax << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		// Maximum offset band value for the reprojection
		if (input.find("BandOffMax") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fBandOffMax=atof(input.c_str());
			if(fBandOffMax<0) {
				WARNING << "BandOffMax should be positive!!! ==> EXIT" << std::endl;
				INFO << "BandOffMax = " << fBandOffMax << std::endl;
				exit(EXIT_FAILURE);
			}
			if(fBandOffMax>2.5) {
				WARNING << "BandOffMax should be inferior to 2.5!!! ==> EXIT" << std::endl;
				INFO << "BandOffMax = " << fBandOffMax << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		// Maximum efficacite band value for the reprojection
		if (input.find("BandEffMax") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fBandEffMax=atof(input.c_str());
			if(fBandEffMax<0) {
				WARNING << "BandEffMax should be positive!!! ==> EXIT" << std::endl;
				INFO << "BandEffMax = " << fBandEffMax << std::endl;
				exit(EXIT_FAILURE);
			}
			if(fBandEffMax>100) {
				WARNING << "BandEffMax should be inferior to 100%!!! ==> EXIT" << std::endl;					    
				INFO << "BandEffMax = " << fBandEffMax << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		//Minimum zenith band value for the reprojection
		if (input.find("BandZenMin") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fBandZenMin=atof(input.c_str());
			if(fBandZenMin<0) {
				WARNING << "BandZenMin should be positive!!! ==> EXIT" << std::endl;
				INFO << "BandZenMin = " << fBandZenMin << std::endl;
				exit(EXIT_FAILURE);
			}
			if(fBandZenMin>70) {
				WARNING << "BandZenMin should be inferior to 70degrees!!! ==> EXIT" << std::endl;
				INFO << "BandZenMin = " << fBandZenMin << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		// Minimum offset band value for the reprojection
		if (input.find("BandOffMin") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fBandOffMin=atof(input.c_str());
			if(fBandOffMin<0) {
				WARNING << "BandOffMin should be positive!!! ==> EXIT" << std::endl;
				INFO << "BandOffMin = " << fBandOffMin << std::endl;
				exit(EXIT_FAILURE);
			}
			if(fBandOffMin>2.5) {
				WARNING << "BandOffMin should be inferior to 2.5!!! ==> EXIT" << std::endl;
				INFO << "BandOffMin = " << fBandOffMin << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		// Minimum efficacite band value for the reprojection
		if (input.find("BandEffMin") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fBandEffMin=atof(input.c_str());
			if(fBandEffMin<0) {
				WARNING << "BandEffMin should be positive!!! ==> EXIT" << std::endl;
				INFO << "BandEffMin = " << fBandEffMin << std::endl;
				exit(EXIT_FAILURE);
			}
			if(fBandEffMin>100) {
				WARNING << "BandEffMin should be inferior to 100%!!! ==> EXIT" << std::endl;
				INFO << "BandEffMin = " << fBandEffMin << std::endl;
				exit(EXIT_FAILURE);
			}
		}
    
		//Number of zenith bin for the first level of reprojection
		if (input.find("BandZenbin_UseReproj") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fBandZenbin_UseReproj=atoi(input.c_str());
			if(fBandZenbin_UseReproj<0) {
				WARNING << "BandZenbin_UseReproj should be positive!!! ==> EXIT" << std::endl;
				INFO << "BandZenbin_UseReproj = " << fBandZenbin_UseReproj << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		//Number of offset bin for the first level of reprojection
		if (input.find("BandOffbin_UseReproj") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fBandOffbin_UseReproj=atoi(input.c_str());
			if(fBandOffbin_UseReproj<0) {
				WARNING << "BandOffbin_UseReproj should be positive!!! ==> EXIT" << std::endl;
				INFO << "BandOffbin_UseReproj = " << fBandOffbin_UseReproj << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		//Number of efficiency bin for the first level of reprojection
		if (input.find("BandEffbin_UseReproj") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fBandEffbin_UseReproj=atoi(input.c_str());
			if(fBandEffbin_UseReproj<0) {
				WARNING << "BandEffbin_UseReproj should be positive!!! ==> EXIT" << std::endl;
				INFO << "BandEffbin_UseReproj = " << fBandEffbin_UseReproj << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		//Number of zenith bin for the second level of reprojection
		if (input.find("BandZenbin_SmallStat") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fBandZenbin_SmallStat=atoi(input.c_str());
			if(fBandZenbin_SmallStat<0) {
				WARNING << "BandZenbin_SmallStat should be positive!!! ==> EXIT" << std::endl;
				INFO << "BandZenbin_SmallStat = " << fBandZenbin_SmallStat << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		//Number of offset bin for the second level of reprojection
		if (input.find("BandOffbin_SmallStat") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fBandOffbin_SmallStat=atoi(input.c_str());
			if(fBandOffbin_SmallStat<0) {
				WARNING << "BandOffbin_SmallStat should be positive!!! ==> EXIT" << std::endl;
				INFO << "BandOffbin_SmallStat = " << fBandOffbin_SmallStat << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		//Number of efficiency bin for the second level of reprojection
		if (input.find("BandEffbin_SmallStat") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fBandEffbin_SmallStat=atoi(input.c_str());
			if(fBandEffbin_SmallStat<0) {
				WARNING << "BandEffbin_SmallStat should be positive!!! ==> EXIT" << std::endl;
				INFO << "BandEffbin_SmallStat = " << fBandEffbin_SmallStat << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		//Number of zenith bin for the third level of reprojection
		if (input.find("BandZenbin_VerySmallStat") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fBandZenbin_VerySmallStat=atoi(input.c_str());
			if(fBandZenbin_VerySmallStat<0) {
				WARNING << "BandZenbin_VerySmallStat should be positive!!! ==> EXIT" << std::endl;
				INFO << "BandZenbin_VerySmallStat = " << fBandZenbin_VerySmallStat << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		//Number of offset bin for the third level of reprojection
		if (input.find("BandOffbin_VerySmallStat") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fBandOffbin_VerySmallStat=atoi(input.c_str());
			if(fBandOffbin_VerySmallStat<0) {
				WARNING << "BandOffbin_VerySmallStat should be positive!!! ==> EXIT" << std::endl;
				INFO << "BandOffbin_VerySmallStat = " << fBandOffbin_VerySmallStat << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		//Number of efficiency bin for the third level of reprojection
		if (input.find("BandEffbin_VerySmallStat") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fBandEffbin_VerySmallStat=atoi(input.c_str());
			if(fBandEffbin_VerySmallStat<0) {
				WARNING << "BandEffbin_VerySmallStat should be positive!!! ==> EXIT" << std::endl;
				INFO << "BandEffbin_VerySmallStat = " << fBandEffbin_VerySmallStat << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		//LJ
		// Light curve time cutting type
		if (input.find("UserLightCurveTimeCuttingType") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			unsigned int type(0);
			type = atoi(input.c_str());
			fUserLightCurveTimeCuttingType.push_back(type);
		}

		// Light curve time range
		if (input.find("UserLightCurveTimeRange") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			int strsize=input.size();
			int coma=input.find(",");
			fUserLightCurveTimeRange=std::make_pair<double,double>(atof(input.substr(0,coma).c_str()),
																   atof(input.substr(coma+1,strsize-coma).c_str()));
		}

		// Light curve time interval
		if (input.find("UserLightCurveTimeInterval") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			double interval(atof(input.c_str()));
			fUserLightCurveTimeInterval=interval;
		}

		// Light curve integrated flux energy range
		if (input.find("UserLightCurveIntegratedFluxEnergyRange") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			int strsize=input.size();
			int coma=input.find(",");
			fUserLightCurveIntegratedFluxEnergyRange=std::make_pair<double,double>(atof(input.substr(0,coma).c_str()),
																				   atof(input.substr(coma+1,strsize-coma).c_str()));
		}

		// Light curve integrated flux energy range
		if (input.find("UserLightCurveUserTimeIntervals") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			int strsize=input.size();
			int coma=input.find(",");
			fUserLightCurveUserTimeIntervals.push_back(std::make_pair<double,double>(atof(input.substr(0,coma).c_str()),
																					 atof(input.substr(coma+1,strsize-coma).c_str())));
		}

		if (input.find("UserLightCurveTimeAxisUnits") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserLightCurveTimeAxisUnits=(unsigned int)atoi(input.c_str());
		}

		// Light curve time cutting type
		if (input.find("UserLightCurveErrorsHandling") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			unsigned int type(0);
			type = atoi(input.c_str());
			fUserLightCurveErrorsHandling.push_back(type);
		}

		//By default, the runs or bands with energy threshold higher than the UserFitEMin ares skipped
		//We can force to keep these runs by setting to 1 the UserSkipHighThreshRuns parameter
		if (input.find("UserSkipHighThreshRuns") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserSkipHighThreshRuns=atof(input.c_str());
		}



		// Eref
		if (input.find("UserEref") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserEref=atof(input.c_str());
		}

		// Rebin
		if (input.find("UserHypothesisRebin") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserHypothesisRebin.push_back(atof(input.c_str()));
		}

		// Selection of Hypothesis
		if (input.find("UserAddHypothesis") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			unsigned int hypo(0);
			hypo = atoi(input.c_str());
			fUserAddHypothesis.push_back(hypo);
		}

		//Selection of butterfly
		if (input.find("UserButterfly") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserButterfly.push_back(atoi(input.c_str()));
		}

		// Hypothesis comparison
		if (input.find("UserHypothesisComparison") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserHypothesisComparison=(bool)atoi(input.c_str());
		}
    
		if (input.find("UserDrawDataSummaryPlots") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserDrawDataSummaryPlots=(bool)atoi(input.c_str());
		}

		if (input.find("UserDrawPlots") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserDrawPlots=(bool)atoi(input.c_str());
		}

		if (input.find("UserResidualsStyle") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserResidualsStyle=atoi(input.c_str());
		}
		
		if (input.find("UserPlotStyle") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserPlotStyle=(unsigned int)atoi(input.c_str());
		}

		if (input.find("UserFitResultsDrawing") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserFitResultsDrawing=(bool)atoi(input.c_str());
		}

		if (input.find("UserComputeContours") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserComputeContours=(unsigned int)atoi(input.c_str());
		}

		if (input.find("UserContoursNumberOfPoints") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserContoursNumberOfPoints=(unsigned int)atoi(input.c_str());
		}

		if (input.find("UserComputeScanLikelihood") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserComputeScanLikelihood=(unsigned int)atoi(input.c_str());
		}

		if (input.find("UserMinimizationLevel") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserMinimizationLevel=(unsigned int)atoi(input.c_str());
		}

		//energy range automaticaly calculated in spectrum
		if (input.find("UserAutomaticFit") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserAutomaticFit=atoi(input.c_str());
		}
    
		if (input.find("UserCodeVersion") != std::string::npos) {
			int n=input.rfind(" ");       //search the beginning of the user substring
			input.erase(0,n+1);           //let in the string only the user substring
			fUserCodeVersion=input;
		}
    
	}
 
	configfile.close();

	if(DEBUG) PrintConfig();

	return 1;

}

 
/**
 * \brief Print infos about configuration
 */

void START::Config::PrintConfig() {
	INFO  << "Reading Config from file " << fFileConfig << " :" << std::endl;

	// User
	//std::cout << "\033[34;01m" << "[User]" << "\033[0m" << std::endl;
	INFO << "[User]" << std::endl;
	std::cout << "UserOutputFolderName =                   " << fUserOutputFolderName << std::endl;
	std::cout << "UserUseInputFileList =                   " << fUserUseInputFileList << std::endl;
	std::cout << "UserRunListsFolderAddress =              " << fUserRunListsFolderAddress << std::endl;
	std::cout << "UserRunListFileBaseName =                " << fUserRunListFileBaseName << std::endl; 
	std::cout << "UserOutputRootFileName =                 " << fUserOutputRootFileName << std::endl;
	// Data
	//std::cout << "\033[34;01m" << "[Data]" << "\033[0m" << std::endl
	INFO << "[Data]" << std::endl;
	std::cout << "UserRootFilesFolderAddress =             " << fUserRootFilesFolderAddress << std::endl; 
	std::cout << "UserOffsetMax =                          " << fUserOffsetMax << std::endl;
	std::cout << "UserZenithMax =                          " << fUserZenithMax << std::endl;
	std::cout << "UserMaxCosZenBinWidth =                  " << fUserMaxCosZenBinWidth << std::endl; 
	std::cout << "UserTmin =                               " << fUserTmin << std::endl;
	if(fUserMJDwindow.size()>0) {
		for(unsigned int imjd(0); imjd<fUserMJDwindow.size(); imjd++) {
			std::cout << "MJD window :                             " << "[" << fUserMJDwindow[imjd].first
					  << " - " << fUserMJDwindow[imjd].second << "]" << std::endl;
		}
	}
	else {
		std::cout << "UserMJDwindow =                          " << "No UserMJDwindow specified" << std::endl;
	}

	// Analysis
	//std::cout << "\033[34;01m" << "[Analysis]" << "\033[0m" << std::endl;
	INFO << "[Analysis]" << std::endl;
	std::cout << "UserAnalysisConfig =                     " << fUserAnalysisConfig << std::endl;
	std::cout << "UserAnalysisConfigType =                 " << fUserAnalysisConfigType << std::endl;
	std::cout << "UserForceAreaAnalysisConfig =            " << fUserForceAreaAnalysisConfig << std::endl;
	std::cout << "UserAreaMin =                            " << fUserAreaMin << std::endl;
	std::cout << "UserERangeBinsNumber =                   " << fUserERangeBinsNumber << std::endl;   
	std::cout << "UserERangeMin =                          " << fUserERangeMin << std::endl;
	std::cout << "UserERangeMax =                          " << fUserERangeMax << std::endl;

	if(fUserFitEMin!=-1 && fUserFitEMax!=-1) {
		std::cout << "UserFitEMin =                            " << fUserFitEMin << std::endl;
		std::cout << "UserFitEMax =                            " << fUserFitEMax << std::endl;
	}
	else {
		std::cout << "UserFitEMin =                            " << "not specified" << std::endl;
		std::cout << "UserFitEMax =                            " << "not specified" << std::endl;
	}

	std::cout << "BandZenMax =                             " << fBandZenMax << std::endl;
	std::cout << "BandOffMax =                             " << fBandOffMax << std::endl;
	std::cout << "BandEffMax =                             " << fBandEffMax << std::endl;
	std::cout << "BandZenMin =                             " << fBandZenMin << std::endl;
	std::cout << "BandOffMin =                             " << fBandOffMin << std::endl;
	std::cout << "BandEffMin =                             " << fBandEffMin << std::endl;
	std::cout << "BandZenbin_UseReproj =                             " << fBandZenbin_UseReproj << std::endl;
	std::cout << "BandOffbin_UseReproj =                             " << fBandOffbin_UseReproj << std::endl;
	std::cout << "BandEffbin_UseReproj =                             " << fBandEffbin_UseReproj << std::endl;
	std::cout << "BandZenbin_SmallStat =                             " << fBandZenbin_SmallStat << std::endl;
	std::cout << "BandOffbin_SmallStat =                             " << fBandOffbin_SmallStat << std::endl;
	std::cout << "BandEffbin_SmallStat =                             " << fBandEffbin_SmallStat << std::endl;
	std::cout << "BandZenbin_VerySmallStat =                             " << fBandZenbin_VerySmallStat << std::endl;
	std::cout << "BandOffbin_VerySmallStat =                             " << fBandOffbin_VerySmallStat << std::endl;
	std::cout << "BandEffbin_VerySmallStat =                             " << fBandEffbin_VerySmallStat << std::endl;
	std::cout << "UserComputeContours =                    " << fUserComputeContours << std::endl; 
	std::cout << "UserContoursNumberOfPoints =             " << fUserContoursNumberOfPoints << std::endl; 
	std::cout << "UserComputeScanLikelihood =              " << fUserComputeScanLikelihood << std::endl; 
	std::cout << "UserMinimizationLevel =                  " << fUserMinimizationLevel << std::endl; 
	std::cout << "UserSkipHighThreshRuns =                 " << fUserSkipHighThreshRuns << std::endl; 

	// LightCurve
	//std::cout << "\033[34;01m" << "[LightCurve]" << "\033[0m" << std::endl;
	INFO << "[LightCurve]" << std::endl;
	if(fUserLightCurveTimeCuttingType.size()==0)
		std::cout << "UserLightCurveTimeCuttingType            " << "no light curve" << std::endl;
	else {
		for(unsigned int type(0); type<fUserLightCurveTimeCuttingType.size(); type++) {
			std::cout << "UserLightCurveTimeCuttingType            " << fUserLightCurveTimeCuttingType[type] << std::endl;
			if(fUserLightCurveTimeCuttingType[type]==10 && fUserLightCurveUserTimeIntervals.size()==0) {
				std::cout << "fUserLightCurveUserTimeIntervals         " << "[" << "not specified" << std::endl;
				WARNING << "Troubles! ==> EXIT" << std::endl;
				exit(EXIT_FAILURE);
			}
			if(fUserLightCurveUserTimeIntervals.size()>0 && fUserLightCurveTimeCuttingType[type]==10) {
				for(std::vector<std::pair<double,double> >::iterator time=fUserLightCurveUserTimeIntervals.begin(); 
					time!=fUserLightCurveUserTimeIntervals.end(); ++time) {
					std::cout << "fUserLightCurveUserTimeIntervals         " << "[" << time->first
							  << ":" << time->second << "]" << std::endl;
				}
			}
		}
		std::cout << "UserLightCurveTimeRange                  " << "[" << fUserLightCurveTimeRange.first
				  << ":" << fUserLightCurveTimeRange.second << "]" << std::endl;
		std::cout << "UserLightCurveTimeInterval               " << fUserLightCurveTimeInterval << std::endl;
		std::cout << "UserLightCurveIntegratedFluxEnergyRange  " << "[" << fUserLightCurveIntegratedFluxEnergyRange.first
				  << ":" << fUserLightCurveIntegratedFluxEnergyRange.second << "]" << std::endl;
		std::cout << "UserLightCurveTimeAxisUnits              " << fUserLightCurveTimeAxisUnits << std::endl;
		for(unsigned int err(0); err<fUserLightCurveErrorsHandling.size(); err++) {
			std::cout << "UserLightCurveErrorsHandling             " << fUserLightCurveErrorsHandling[err] << std::endl;	
		}
	}
	// MC
	//std::cout << "\033[34;01m" << "[MC]" << "\033[0m" << std::endl;
	INFO << "[MC]" << std::endl;
	std::cout << "UserIrfOpt =                             " << fUserIrfOpt << std::endl;
	std::cout << "UserMcProductionName =                   " << fUserMcProductionName << std::endl;  
	std::cout << "UserIsMc =                               " << fUserIsMc << std::endl; 

	// hypothesis
	//std::cout << "\033[34;01m" << "[Hypothesis]" << "\033[0m" << std::endl;
	INFO << "[Hypothesis]" << std::endl;
	std::string pwl = "PowerLaw";
	std::string exp = "ExpoCutOffPowerLaw";
	std::string log = "LogParabolic";
	std::string bpwl = "BrokenPowerLaw";
	std::string sbpwl = "SmoothBrokenPowerLaw";
	std::string sexp = "SuperExpoCutOffPowerLaw";
	for(unsigned int ihyp(0); ihyp<fUserAddHypothesis.size(); ihyp++) {
		if(fUserAddHypothesis[ihyp]==1)   std::cout << "UserAddHypothesis =                      " << pwl << std::endl;
		if(fUserAddHypothesis[ihyp]==2)   std::cout << "UserAddHypothesis =                      " << exp << std::endl;
		if(fUserAddHypothesis[ihyp]==3)   std::cout << "UserAddHypothesis =                      " << log << std::endl;
		if(fUserAddHypothesis[ihyp]==4)   std::cout << "UserAddHypothesis =                      " << bpwl << std::endl;
		if(fUserAddHypothesis[ihyp]==5)   std::cout << "UserAddHypothesis =                      " << sexp << std::endl;
		if(fUserAddHypothesis[ihyp]==6)   std::cout << "UserAddHypothesis =                      " << sbpwl << std::endl;
	}

	std::cout << "UserEref =                               " << fUserEref << std::endl; 
	std::cout << "UserHypothesisIntegrationEnergyRange =   " << "[" << fUserHypothesisIntegrationEnergyRange.first 
			  << ";" << fUserHypothesisIntegrationEnergyRange.second << "]" << std::endl;
	std::cout << "UserHypothesisComparison =               " << fUserHypothesisComparison << std::endl; 

	// plots
	//std::cout << "\033[34;01m" << "[Plots]" << "\033[0m" << std::endl;
	INFO << "[Plots]" << std::endl;
	if(fUserDrawPlots==0) 
		std::cout << "UserDrawPlots =                          " << "no plots" << std::endl;
	else {
		std::cout << "UserDrawPlots =                          " << fUserDrawPlots << std::endl; 
		std::cout << "UserPlotStyle =                          " << fUserPlotStyle << std::endl;
		std::cout << "UserResidualsStyle =                     " << fUserResidualsStyle << std::endl; 
		std::cout << "UserSourceName =                         " << fUserSourceName << std::endl; 
		for(unsigned int ibutt(0); ibutt<fUserButterfly.size(); ibutt++) {
			std::cout << "UserButterfly =                          " << fUserButterfly[ibutt] << std::endl; 
		}
		for(unsigned int ireb(0); ireb<fUserHypothesisRebin.size(); ireb++) {
			std::cout << "UserHypothesisRebin =                    " << fUserHypothesisRebin[ireb] << std::endl; 
		}
		std::cout << "UserFitResultsDrawing =                  " << fUserFitResultsDrawing << std::endl; 
		std::cout << "UserDrawDataSummaryPlots =               " << fUserDrawDataSummaryPlots << std::endl; 
	}

	// keep those : ?
	std::cout << "JLK :  Keep these options?" << std::endl;
	std::cout << "UserCodeVersion =                        " << fUserCodeVersion << std::endl; 
	std::cout << "UserInstrumentResponseAddress =          " << fUserInstrumentResponseAddress << std::endl;

}

/**
 * \brief Read the run numbers and the associated telcode 
 * and stock them in a vector. Alternatively read a file list.
 */

int START::Config::ReadList(TString listname)
{
	fFileList = listname;

	std::ifstream runlist(listname.Data());
	if (!runlist.is_open()){
		std::cout << "Problem opening list file " << listname.Data() << std::endl;
		return 0;
	}
  
	std::string input;

	if (fUserUseInputFileList) {

		//read input file list and fill run list with a counter
		int count = 0;
		while (!runlist.eof()){     
			std::getline(runlist, input, '\n');

			if (input[0]=='#') continue; //excludes commented lines
			if (input=="") continue; //excludes empty line

			fInputFileList.push_back(input);
			count++;
			fRunList.push_back(count);
			fAllowedTelList.push_back(0); //MPA: TODO: try to avoid this line!!!
		}

	}
	else {

		//read run list with tel pattern
		while (1){
			std::getline(runlist, input, '\n');
			if((input.find("*") != std::string::npos) || (input.find("h") != std::string::npos)) {
				runlist >> input;
				continue;
			}

			if(runlist.eof()) break;
			if (input[0]=='#') continue; //excludes commented lines
			if (input=="") continue; //excludes empty line
			std::vector<std::string> inlist = split(input.c_str(),' ');

			if (inlist.size()==0) continue;
			int run = atoi(inlist[0].c_str());
			if(run == 0) continue;
			int allowedtels = 0;
			if(inlist.size()>1) allowedtels = atoi(inlist[1].c_str());

			if(allowedtels != 0) {
				fRunList.push_back(run);
				fAllowedTelList.push_back(allowedtels);
			}
			else {
				//std::cout << "Rejecting Run " << run
				//                     << ", all telescopes failed run selection." << std::endl;
			}
		}

	}
  
	runlist.close();
  
	return 1;
}

/**
 * \brief Print the runs to process
 */
void START::Config::PrintRunList() {

	INFO << fRunList.size() << " runs to process : " << std::endl;
	for(unsigned int irun(0); irun<fRunList.size(); irun++) 
		std::cout << "run " << fRunList[irun] << " with telcode " << fAllowedTelList[irun] << std::endl;

}

/**
 * \brief Print the input files to process
 */
void START::Config::PrintInputFileList() {

	INFO << fInputFileList.size() << " input files to process : " << std::endl;
	for(unsigned int ifile(0); ifile<fInputFileList.size(); ifile++) 
		std::cout << "file " << fInputFileList[ifile] << std::endl;

}

