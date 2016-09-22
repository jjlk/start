#ifndef _MINIMIZEFACTORY_
#define _MINIMIZEFACTORY_

// STL
#include <vector>
#include <string>
#include <map>
#include <utility>

// ROOT
#include <TObject.h>
#include <TStopwatch.h>

// Minuit2
#include <Minuit2/FunctionMinimum.h>
#include "Minuit2/MnMinimize.h"

// START
#include "STARTUtils.hh"
#include "Config.hh"

namespace START {
	class Band;
	class Hypothesis;


	/**
	 * \brief Used to minimize one or more Hypothesis.
	 *
	 * Minimization is done with MINUIT2.
	 * You can minimize one hypothesis or a vector of hypothesis, 
	 * depending on which constructor is called.
	 *
	 * \warning copy constructor and assignment operator are not necessary 
	 * thus not implemented.
	 *
	 * \warning At this step, the we store the BandArray inside the hypothesis. 
	 * Thus it became different than the one the user gave.... 
	 * (the theoritical signal is filled inside the hypothesis)

	 * \author HAP-Fr team
	 */
	class MinimizeFactory : public TObject
	{

	public :

		/**
		 * \brief Options for the minimization for differential hypothesis.
		 * \param Standard : minimization consists of computed fitted parameters for hypothesis and integrated flux
		 * \param Minos : minimization consists of computed fitted parameters for hypothesis and by doing Minos analysis
		 * \param Light : minimization consists of computed fitted parameters for hypothesis
		 * \param Medium : minimization consists of computed fitted parameters for hypothesis, integrated flux and energyflux
		 * \param Full : minimization consists of computed fitted parameters and doing Minos for hypothesis, then computed integrated and energy flux
		 */
		typedef enum {Standard, Minos, Light, Medium, Full} MinimizationType;

		MinimizeFactory(); // needed by ROOT
		MinimizeFactory(Hypothesis &hypo, Config &Configuration, bool verbose=false); // defaut constructor, need the data and the hypothesis
		MinimizeFactory(std::vector<START::Hypothesis*> &hypo, Config &Configuration, bool verbose=false); //constructor, need the data and the vector of hypothesis*

		~MinimizeFactory(); // destructor

		void SetBandsGroupingType(const START::STARTUtils::BandsGroupingTypeForMinimization Type) {fBandsGroupingType = Type;};

		void MakeMinimization(std::vector<START::Band> const &BandArray, MinimizationType type=Light); // launch minimization

		void MakeContours(std::vector<START::Band> const &BandArray,std::string hyponame, unsigned int param1, 
						  unsigned int param2,unsigned int npoints=20); // make contours for 2 parameters
		void MakeContours(std::vector<START::Band> const &BandArray,std::string hyponame, std::string param1, 
						  std::string param2,unsigned int npoints); // make contours for 2 parameters

		void MakeAllContours(std::vector<START::Band> const &BandArray, int npoints=20); // make contours for all parameters
  
		void ScanLikelihoodAfterMinimization(Hypothesis &hypo,std::vector<START::Band> &BandArray, unsigned param,
											 unsigned int npoints=40, double low=-999, double high=-999); // scan region

		void ScanLikelihoodAfterMinimization(Hypothesis &hypo,std::vector<START::Band> &BandArray,std::string param,
											 unsigned int npoints=40, double low=-999, double high=-999); // scan region
  
		void ScanLikelihoodForAllHypothesisAfterMinimization(std::vector<START::Band> &BandArray, unsigned int npoints=40); // scan region for all param and hypo

		void ComputeAndPrintLikelihoodRatio(); // compute and print likelihood ratio

		void CompareHypothesisByMinimization(std::vector<START::Band> &BandArray, Hypothesis &hypo1, Hypothesis &hypo2); // compare hypothesis with minimization

		void CompareAllHypothesisByMinimization(std::vector<START::Band> &BandArray); // compare all hypothesis with minimization

		void PrintMinimizationSummary(); // print summary of the minimiation

		std::pair<double,double> ComputeMinimizationEnergyRange(const std::vector<START::Band> BandArray); // get minimization energy range

	private :

		void ComputeMinosErrors(const std::vector<START::Band> &BandArray, Hypothesis &hypo, ROOT::Minuit2::FunctionMinimum &min);
		void ComputeIntegratedFlux(const std::vector<START::Band> &BandArray, Hypothesis &hypo);
		void ComputeEnergyFlux(const std::vector<START::Band> &BandArray, Hypothesis &hypo);

		void CopyExcessInBandsHypothesis(Hypothesis *hypo, std::vector<START::Band> const &BandArray);
  
		void PrintMinimizationCarateristics(Hypothesis &hypo); ///< print minimization caracteristics such as edm, iterations' number and minimum
		void PrintFittedParametersWithErrors(Hypothesis &hypo); ///< print fitted parameters with their errors
		void PrintCovarianceMatrix(Hypothesis &hypo); ///< print covariance matrix
		void PrintDecorrelationEnergy(Hypothesis &hypo); ///< print covariance matrix

		ROOT::Minuit2::MnUserParameterState GetUserParameterForMinuit2(const Hypothesis &hypo);

		std::vector<std::pair<bool,bool> > GetVectorPairFixedParamRedoMinimization(const Hypothesis &hypo);

		void ComputeContours(const std::vector<START::Band> &BandArray, Hypothesis &hypothesis, 
							 ROOT::Minuit2::FunctionMinimum &funcminimum, 
							 unsigned int param1, unsigned int param2, 
							 unsigned int npoints); ///< Compute contours at 1, 2 and 3 sigma

		bool IsThereAnyNaNInFittedParameters(std::vector<double> fittedparams, std::vector<double> fittedparamerrors);

		void PrintResultsAndAddInfosInHypothesis(Hypothesis &hypo,
												 const ROOT::Minuit2::FunctionMinimum &MinimumFromMINUIT2,
												 const ROOT::Minuit2::MnUserCovariance &CovFromMINUIT2,
												 TStopwatch &TimeFromROOT);
		void PrintResultsAndAddInfosInHypothesisWithoutFitting(Hypothesis &hypo,
															   const ROOT::Minuit2::MnMinimize &MigradSimplex,
															   TStopwatch &TimeFromROOT);


		MinimizationType fMinimizationType;

		bool fverbose; ///< verbose level

		/* caracteristics of the minimization */
		unsigned int fmaxcall; ///< maximum number of iterations
		double ftolerance; ///< EDM_min for a good fit is at least 0.001*0.5*tolerance
		double frequiredEDM; ///< required EDM
		unsigned int fstrategy; ///< strategy, 0 for low precision and 2 for high precision

		/* Contours */
  
		std::vector<std::pair<std::vector<std::pair<double,double> >, std::pair<int,int> > > fcontoursigma1; ///< contours at 1 sigma (points,param1,param2)
		std::vector<std::pair<std::vector<std::pair<double,double> >, std::pair<int,int> > > fcontoursigma2; ///< contours at 1 sigma (points,param1,param2)
		std::vector<std::pair<std::vector<std::pair<double,double> >, std::pair<int,int> > > fcontoursigma3; ///< contours at 1 sigma (points,param1,param2)

		std::map<std::string,std::vector<std::pair<std::vector<std::pair<double,double> >, std::pair<int,int> > > > fMapContours1Sigma;
		std::map<std::string,std::vector<std::pair<std::vector<std::pair<double,double> >, std::pair<int,int> > > > fMapContours2Sigma;
		std::map<std::string,std::vector<std::pair<std::vector<std::pair<double,double> >, std::pair<int,int> > > > fMapContours3Sigma;

		/* Scan regions */

		std::vector<std::pair<unsigned int, std::vector<std::pair<double,double> > > > fregionscans; ///< scans (param,points)

		std::map<std::string, std::vector<std::pair<unsigned int, std::vector<std::pair<double,double> > > > > fMapScanLikeliHood;

		/* Objects */

		std::vector<std::pair<START::Hypothesis*,ROOT::Minuit2::FunctionMinimum*> > fPairHypothesisFunctionMinimum; ///< vector of pair containing hypothesis and associated function minimum

		// Type of grouping
		START::STARTUtils::BandsGroupingTypeForMinimization fBandsGroupingType;

		//ADA Give access to the configuration file 
		Config *fConfig;

#ifdef CONSTRUCT_LIBRARY
	public:
		ClassDef(START::MinimizeFactory,1);
#endif
	};
}
#endif

