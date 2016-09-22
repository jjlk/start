// STL
#include <iostream>
#include <iomanip>

// ROOT
#include <TMath.h>

// MINUIT2
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnSimplex.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMachinePrecision.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnScan.h"

// START
#include "Hypothesis.hh"
#include "FCNLikelihood.hh"
#include "MinimizeFactory.hh"
#include "STARTUtils.hh"
#include "Band.hh"
#include "SumHypothesis.hh"
#include "DataSummary.hh"

// Utilities
//#define DEBUG 1000
#include "debugging.hh"

#define INFO std::cout << INFOCOLOR << "MinimizeFactory> " << RESETCOLOR
#define WARNING std::cout << WARNINGCOLOR << "MinimizeFactory> " << RESETCOLOR

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::MinimizeFactory)
#endif


//ADA
extern int gErrorIgnoreLevel;


/**
 * \brief Constructor needed by ROOT
 *
 */
START::MinimizeFactory::MinimizeFactory()
	:fmaxcall(200000),
	 ftolerance(0.01),
	 frequiredEDM(0.),
	 fstrategy(2)
{

}

/**
 * \brief Constructor to minimize on hypothesis
 */
START::MinimizeFactory::MinimizeFactory(Hypothesis &hypo, Config &Configuration, bool verbose)
	:fverbose(verbose),
	 fmaxcall(200000),
	 ftolerance(0.01),
	 frequiredEDM(0.),
	 fstrategy(2)
{

        fConfig = &Configuration;
	// We add a composant in the vector of the hypothesis
	// and the functionminimum is allocated later
	fPairHypothesisFunctionMinimum.clear();
	std::pair<Hypothesis*,ROOT::Minuit2::FunctionMinimum*> init;
	init.first = &hypo;
	init.second = 0;
	fPairHypothesisFunctionMinimum.push_back(init);

	fcontoursigma1.clear();
	fcontoursigma2.clear();
	fcontoursigma3.clear();

}

/**
 * \brief Constructor to minimize a vector of hypothesis hypothesis
 */
START::MinimizeFactory::MinimizeFactory(std::vector<Hypothesis*> &hypo, Config &Configuration, bool verbose)
	:fverbose(verbose),
	 fmaxcall(200000),
	 ftolerance(0.01),
	 //ftolerance(1),
	 frequiredEDM(0.),
	 fstrategy(2)
{

        fConfig = &Configuration;

	// We add composants in the vector pair : hypothesis, the functionminimum is allocate later
	fPairHypothesisFunctionMinimum.clear();

	for(std::vector<Hypothesis*>::iterator it=hypo.begin(); it!=hypo.end(); ++it) {
		std::pair<Hypothesis*,ROOT::Minuit2::FunctionMinimum*> init;
		init.first = &(**it);
		init.second = 0;
		fPairHypothesisFunctionMinimum.push_back(init);
	}

	fcontoursigma1.clear();
	fcontoursigma2.clear();
	fcontoursigma3.clear();

}

/**
 * \brief Destructor
 */
START::MinimizeFactory::~MinimizeFactory()
{

	for(std::vector<std::pair<Hypothesis*,ROOT::Minuit2::FunctionMinimum*> >::iterator hypo=fPairHypothesisFunctionMinimum.begin(); 
		hypo!=fPairHypothesisFunctionMinimum.end(); ++hypo) { // loop on hypothesis
		if((*hypo).second!=0) delete (*hypo).second;
		(*hypo).second = 0;
	}
  
}

/**
 * \brief Minimize!
 *
 * \param BandArray vector of bands
 *
 */
void START::MinimizeFactory::MakeMinimization(std::vector<Band> const &BandArray, MinimizationType type)
{

	// define type of minimization

	fMinimizationType=type;

	// Compute minimization energy range

	std::pair<double,double> minimizationenergyrange = ComputeMinimizationEnergyRange(BandArray);

	INFO << "Starting minimization between " << minimizationenergyrange.first 
		 << " and " << minimizationenergyrange.second << " TeV!" << std::endl;

	for(std::vector<std::pair<Hypothesis*,ROOT::Minuit2::FunctionMinimum*> >::iterator hypo=fPairHypothesisFunctionMinimum.begin(); 
		hypo!=fPairHypothesisFunctionMinimum.end(); hypo++) { // loop on hypothesis

		// copy minimization energy range in hypothesis

		(*hypo).first->SetMinimizationEnergyRange(minimizationenergyrange.first,minimizationenergyrange.second);    

		/* Initialisation of parameters */

		//ADA 
		gErrorIgnoreLevel = 1001;
		ROOT::Minuit2::MnUserParameterState Parameters;

		std::vector<double> param_init = (*hypo).first->GetParameters();
		std::vector<double> paramerr_init = (*hypo).first->GetParametersErrors();

		Parameters = GetUserParameterForMinuit2(*(*hypo).first);

		// Minimize minimizer : call MIGRAD then if it fails call SIMPLEX and finally
		// call MIGRAD again
    
		ROOT::Minuit2::MnStrategy Strategy(fstrategy); // set the strategy 0,1,2

		FCNLikelihood *FCN = new FCNLikelihood(BandArray,(*(*hypo).first), *fConfig, fverbose); // likelihood definition
		// JLK add Band regrouping
		FCN->SetBandsGroupingType(fBandsGroupingType);

		ROOT::Minuit2::MnMinimize MigradSimplex(*FCN,Parameters,Strategy); // m

		// If there is fixed parameters :

		std::vector<std::pair<unsigned int, std::pair<double, bool> > > fixedparameters = (*hypo).first->GetFixedParameters();

		if(fixedparameters.size()>0) {

			for(unsigned int ifix(0); ifix<fixedparameters.size(); ifix++) { // We fix the parameters
	
				unsigned int fixedparam = fixedparameters[ifix].first; 
				double fixedvalue = fixedparameters[ifix].second.first;
				bool redominimization = fixedparameters[ifix].second.second;
	
				MigradSimplex.SetValue(fixedparam,fixedvalue);
				MigradSimplex.Fix(fixedparam);

				if(redominimization) { // We minimize if necessary

					ROOT::Minuit2::FunctionMinimum FixMinimumFromMINUIT2 = MigradSimplex(fmaxcall,ftolerance);
					const ROOT::Minuit2::MnUserCovariance CovFixedFromMINUIT2 = FixMinimumFromMINUIT2.UserCovariance();
					DEBUG_OUT << "Minimum from MINUIT2 : " << FixMinimumFromMINUIT2 << std::endl;
					DEBUG_OUT << "Covariance from MINUIT2" << CovFixedFromMINUIT2 << std::endl;
					MigradSimplex.Release(fixedparam); // We release the parameters

				}

			}

		}
		std::vector<std::pair<unsigned int, std::pair<double, double> > > limitedparameters = (*hypo).first->GetLimitedParameters();
    
		if(limitedparameters.size()>0) {
      
			for(unsigned int ilimit(0); ilimit<limitedparameters.size(); ilimit++) { // We fix the parameters
	
				unsigned int limitedparam = limitedparameters[ilimit].first; 
				double limitedvaluemin = limitedparameters[ilimit].second.first;
				double limitedvaluemax = limitedparameters[ilimit].second.second;
	
				MigradSimplex.SetLimits(limitedparam,limitedvaluemin, limitedvaluemax);
	
			}
		}
    

		std::cout << "****************************************************************" << std::endl;

		// start main minimization

		TStopwatch TimeFromROOT;
		TimeFromROOT.Start();
		//ADA Hack to get the result with fixed params  
		if(fixedparameters.size()==(*hypo).first->GetParametersNb()) {
			INFO << "All parameters Fixed ..............." << std::endl;
			CopyExcessInBandsHypothesis((*hypo).first,BandArray);
			INFO << "No Minimization But Excess copied in bands" << std::endl;
			PrintResultsAndAddInfosInHypothesisWithoutFitting(*(*hypo).first, MigradSimplex,TimeFromROOT);
			return;
		}
		ROOT::Minuit2::FunctionMinimum MinimumFromMINUIT2 = MigradSimplex(fmaxcall,ftolerance);
    
		DEBUG_OUT << "Minimum from Minuit : " << MinimumFromMINUIT2 << std::endl;
    
		(*hypo).second = new ROOT::Minuit2::FunctionMinimum(MinimumFromMINUIT2);
    
		TimeFromROOT.Stop();
		const ROOT::Minuit2::MnUserCovariance CovFromMINUIT2 = (*hypo).second->UserCovariance();
    
		DEBUG_OUT << "DEBUG_OUT <<Covariance from MINUIT2 : " << CovFromMINUIT2 << std::endl;
    
		PrintResultsAndAddInfosInHypothesis(*(*hypo).first,MinimumFromMINUIT2,CovFromMINUIT2,TimeFromROOT);
		DEBUG_OUT << "apres PrintResultsAndAddInfosInHypothesis" << std::endl;
		// Lea hack
		delete FCN; FCN=0;
		if(!(*hypo).first->GetIsMinimumValid() || !(*hypo).first->GetIsCovarianceValid() || !(*hypo).first->GetAreFittedParametersValid())
			continue;
	
		//delete FCN; FCN=0;

		if((*hypo).first->GetSpectralType()==Hypothesis::Differential) {

			switch(fMinimizationType) {
			case Light:
				break;
			case Medium:
				ComputeIntegratedFlux(BandArray,*(*hypo).first);
				ComputeEnergyFlux(BandArray,*(*hypo).first);
				break;
			case Full:
				ComputeMinosErrors(BandArray,*(*hypo).first,*(*hypo).second);
				ComputeIntegratedFlux(BandArray,*(*hypo).first);
				ComputeEnergyFlux(BandArray,*(*hypo).first);
				break;
			case Standard:
				ComputeIntegratedFlux(BandArray,*(*hypo).first);
				break;
			case Minos:
				ComputeMinosErrors(BandArray,*(*hypo).first,*(*hypo).second);
			}


		}
		DEBUG_OUT << "after swich fMinimizationType " << std::endl;

		// copy excess in bands 
		CopyExcessInBandsHypothesis((*hypo).first,BandArray);

	} // end loop on hypothesis

	// Print likelihood ratios

	std::cout << "****************************************************************" << std::endl;

	ComputeAndPrintLikelihoodRatio();

	std::cout << "****************************************************************" << std::endl;

	INFO << "Minimization... ok" << std::endl;

}

ROOT::Minuit2::MnUserParameterState START::MinimizeFactory::GetUserParameterForMinuit2(const Hypothesis &hypo) {

	std::vector<std::string> paramfitname = hypo.GetParametersNames();
  
	std::vector<std::pair<unsigned int,double> > normalizedparam = 
		hypo.GetVectorNormalizedParameters();
  
	/* Initialisation of parameters */
  
	ROOT::Minuit2::MnUserParameterState Parameters;
  
	std::vector<double> param_init = hypo.GetParameters();
	std::vector<double> paramerr_init = hypo.GetParametersErrors();
  
	// We normalize the components of the flux which have a lower limit (necessary)
	// All parameters such as phi0 in flux=phi0*E^-gamma must have a lowerlimit because :
	// 1) the likelihood can't have a negative expected excess (it's our choice)
	// 2) We normalize them in a way which make it easy for MINUIT2 (precision's matters)
	// The non-normalized/havenolimit parameters can't have a lower limit because MINUIT2 don't like it
	// (internals changes of variables in MINUIT2) and it's again our choice to forbid it ant at least,
	// to minimize those effects.
  
	const double normalization = STARTUtils::GetNormalizationConstante();
  
	DEBUG_OUT << "Before normalization" << std::endl;

	for(unsigned int ipar(0); ipar<hypo.GetParametersNb(); ipar++) { // loop on all parameters
    
		DEBUG_OUT << paramfitname[ipar] << " :" << std::endl;

		if(normalizedparam.size()>0) {
			// we check if the parameter have to be normalized and we stock the indice and the lower limit value
			unsigned int paramtonormalize(666);
			double lowerlimitvalue(666);
			for(unsigned int ilow(0); ilow<normalizedparam.size(); ilow++) { //loop on the parameter which have a lower limit
				if(normalizedparam[ilow].first==ipar) {
					paramtonormalize = normalizedparam[ilow].first;
					lowerlimitvalue = normalizedparam[ilow].second;
					break;
				}
			}
      
			if(paramtonormalize==666 && lowerlimitvalue==666) { // we don't normalize
				Parameters.Add(paramfitname[ipar],param_init[ipar],paramerr_init[ipar]);
				DEBUG_OUT << "No lower limit on parameter " << paramfitname[ipar] 
						  << " with value " << param_init[ipar] << std::endl;
			}
			else { // we normalize
				Parameters.Add(paramfitname[ipar],param_init[ipar]/normalization,paramerr_init[ipar]/normalization);
				Parameters.SetLowerLimit(paramtonormalize,lowerlimitvalue/normalization);
				DEBUG_OUT << "Lower limit on parameter " << paramtonormalize
						  << " with value " << param_init[ipar] 
						  << " and lower limit " << lowerlimitvalue << std::endl;
			}


		}
		else { // if there is no normalized parameters, we go here (exotic law)
			DEBUG_OUT << "In condition \"no normalized parameters\" " << std::endl;
			Parameters.Add(paramfitname[ipar],param_init[ipar],paramerr_init[ipar]);
			DEBUG_OUT << "No lower limit on parameter " << paramfitname[ipar] << std::endl;
		}

	}

	DEBUG_OUT << "After normalization" << std::endl;

	DEBUG_OUT << Parameters << std::endl;

	return Parameters;

}

/**
 * \brief Build a pair of vector containing two boolean (isparamfixed,redominimization)
 */
std::vector<std::pair<bool,bool> > START::MinimizeFactory::GetVectorPairFixedParamRedoMinimization(const Hypothesis &hypo) {

	std::vector<std::pair<unsigned int, std::pair<double, bool> > > fixedparameters =
		hypo.GetFixedParameters();

	// containing info on re-minimization if any, size of numbers of hypothesis's parameters
	std::vector<std::pair<bool,bool> > isfixedisredominimization;
  
	// filling the vector with true if reminimization is done for param i      
	for(unsigned int iparam(0); iparam<hypo.GetParametersNb(); iparam++) {
    
		isfixedisredominimization.push_back(std::make_pair(false,true));
    
		for(unsigned int ifix(0); ifix<fixedparameters.size(); ifix++) {
      
			unsigned int fixedparam = fixedparameters[ifix].first; 
			bool redo = fixedparameters[ifix].second.second;	
      
			if(fixedparam==iparam && redo) {
				isfixedisredominimization.back().first=true;
				isfixedisredominimization.back().second=true;
	
			}
			else if(fixedparam==iparam && !redo) {
				isfixedisredominimization.back().first=true;
				isfixedisredominimization.back().second=false;
			}
		}
    
	}
  
	for(unsigned int iparam(0); iparam<isfixedisredominimization.size(); iparam++) {
		DEBUG_OUT << "iparam=" << iparam << " isfixed=" << isfixedisredominimization[iparam].first
				  << " redominimization=" << isfixedisredominimization[iparam].second 
				  << std::endl;
	}    

	return isfixedisredominimization;

}

/**
 * \brief Add fitted parameters in hypothesis and print results
 * \param hypo hypothesis
 * \param CovFromMINUIT2 covriance from MINUIT2
 * \param MinimumFromMINUIT2 minimum from MINUIT2
 * \param TimeFromROOT TStopwatch class containing time
 * \param isfixedisredominimization pair of vector containing two boolean (isfixedparam,redominimization)
 */
void START::MinimizeFactory::PrintResultsAndAddInfosInHypothesis(Hypothesis &hypo,
																 const ROOT::Minuit2::FunctionMinimum &MinimumFromMINUIT2,
																 const ROOT::Minuit2::MnUserCovariance &CovFromMINUIT2,
																 TStopwatch &TimeFromROOT)
{

	bool isminimumvalid = MinimumFromMINUIT2.IsValid();
	bool iscovariancevalid = MinimumFromMINUIT2.HasValidCovariance();
	bool arefittedparametersvalid = MinimumFromMINUIT2.HasValidParameters();
	bool convergence(false);
	if(isminimumvalid==true && arefittedparametersvalid==true && iscovariancevalid==true)
		convergence = true;
	else
		convergence = false;

	double fcnminimum = MinimumFromMINUIT2.Fval();
	double edmvalue = MinimumFromMINUIT2.Edm();
	double iteration = MinimumFromMINUIT2.NFcn();
  
	frequiredEDM = 0.001*0.5*ftolerance;
  
	std::vector<double> paramfit;
	std::vector<double> paramfiterr;
	std::vector<std::vector<double> > covariancefit;
	INFO << "frequiredEDM = " << frequiredEDM << std::endl;
	INFO << "Results of the minimization for " << hypo.GetName() << " (Eref=" << hypo.GetEref() << " TeV";

	if(hypo.GetSpectralType()==Hypothesis::Integrated || hypo.GetSpectralType()==Hypothesis::EnergyFlux) {
		std::cout << ", integration range is between " 
				  << hypo.GetIntegratedFitEnergyRange().first
				  << " and "
				  << hypo.GetIntegratedFitEnergyRange().second 
				  << " TeV";
	}

	std::cout << ") : " << std::endl;
  
	// JLK: Trap, can have the convergence status and NaN in fitted parameters or errors... ==> Counter-trap!
	if(IsThereAnyNaNInFittedParameters(MinimumFromMINUIT2.UserState().Params(),MinimumFromMINUIT2.UserState().Errors())) {
		WARNING << "NaN detected in fitted parameters! Your hypothesis will be tagged as non-convergent!" << std::endl;
		arefittedparametersvalid=false;
		iscovariancevalid=false;
		isminimumvalid=false;
		convergence = false;
	}
  
	std::vector<std::pair<unsigned int,double> > normalizedparam = 
		hypo.GetVectorNormalizedParameters();

	const double normalization = STARTUtils::GetNormalizationConstante();

	std::vector<std::pair<bool,bool> > isfixedisredominimization = GetVectorPairFixedParamRedoMinimization(hypo);
  
	/* set caracteristics of the minimization */
	hypo.SetMaximumLikelihood(-fcnminimum);
	hypo.SetEDM(edmvalue);
	hypo.SetIteration(iteration);
	hypo.SetIsMinimumValid(isminimumvalid);
	hypo.SetIsCovarianceValid(iscovariancevalid);
	hypo.SetAreFittedParametersValid(arefittedparametersvalid);
	hypo.SetConvergence(convergence);

	if(convergence) {
    
		INFO << "MINUIT2 has converged!" << std::endl;
    
		if(edmvalue>frequiredEDM) {
			WARNING << "EDM is greater than required EDM " << "(" << frequiredEDM << ")" << std::endl;
		}
    
		/* set fitted parameters and their errors */
		paramfit = MinimumFromMINUIT2.UserState().Params();
		paramfiterr = MinimumFromMINUIT2.UserState().Errors();
    
		/* Denormalization of the parameters */
		for(unsigned int ilow(0); ilow<normalizedparam.size(); ilow++) {
			paramfit[normalizedparam[ilow].first] = paramfit[normalizedparam[ilow].first]*normalization;
			paramfiterr[normalizedparam[ilow].first] = paramfiterr[normalizedparam[ilow].first]*normalization;
		}
    
		if(hypo.GetSpectralType()==Hypothesis::EnergyFlux) {
			paramfit[0]=paramfit[0]*STARTUtils::GetTeVToErg();
			paramfiterr[0]=paramfiterr[0]*STARTUtils::GetTeVToErg();
		}

		hypo.SetFittedParameters(paramfit);
		hypo.SetFittedParametersErrors(paramfiterr);

		/* get the covariance matrix */
    
		covariancefit.resize(hypo.GetParametersNb(),std::vector<double>(hypo.GetParametersNb(),0.));
    
		unsigned int icov(0); // covariance indice
    
		for(unsigned int ipar(0); ipar<hypo.GetParametersNb(); ipar++) {
      
			DEBUG_OUT_L(2) << "##############ipar = " << ipar << std::endl;
      
			unsigned int jcov(0); // covariance indice
      
			for(unsigned int jpar(0); jpar<hypo.GetParametersNb(); jpar++) {
	
				DEBUG_OUT_L(2) << "ipar = " << ipar
							   << " fixed=" << isfixedisredominimization[ipar].first
							   << " isredo=" << isfixedisredominimization[ipar].second << std::endl;
				DEBUG_OUT_L(2) << "jpar = " << jpar
							   << " fixed=" << isfixedisredominimization[jpar].first
							   << " isredo=" << isfixedisredominimization[jpar].second << std::endl;
	
				if( (isfixedisredominimization[ipar].first==false && isfixedisredominimization[ipar].second==true)
					|| (isfixedisredominimization[ipar].first==false && isfixedisredominimization[ipar].second==false)
					|| (isfixedisredominimization[ipar].first==true && isfixedisredominimization[ipar].second==true) ) { //param not fixed or fixed with redo 
	  
					DEBUG_OUT_L(2) << "**b1covi = " << icov << " covj=" << jcov << std::endl;   
	  
					if( (isfixedisredominimization[jpar].first==false && isfixedisredominimization[jpar].second==true)
						|| (isfixedisredominimization[jpar].first==false && isfixedisredominimization[jpar].second==false)
						|| (isfixedisredominimization[jpar].first==true && isfixedisredominimization[jpar].second==true) ) { //param not fixed or fixed with redo
	    
						covariancefit[ipar][jpar]=CovFromMINUIT2(icov,jcov);
						jcov++;   
	    
					}
					else if(isfixedisredominimization[jpar].first==true && isfixedisredominimization[jpar].second==false) { // param fixed with no redo
	    
						covariancefit[ipar][jpar]=0.;
	    
					}
	  
					DEBUG_OUT_L(2) << "**a1covi = " << icov << " covj=" << jcov << std::endl;   
	  
				}
				else if(isfixedisredominimization[ipar].first==true && isfixedisredominimization[ipar].second==false) { //param fixed with no redo
					DEBUG_OUT_L(2) << "**4bcovi = " << icov << " covj=" << jcov << std::endl;
	  
					covariancefit[ipar][jpar]=0.;
	  
					DEBUG_OUT_L(2) << "**4acovi = " << icov << " covj=" << jcov << std::endl;   
				}
	
				if(normalizedparam.size()>0) { // denormalization
	  
					for(unsigned int ilow(0); ilow<normalizedparam.size(); ilow++) {
						if(ipar==normalizedparam[ilow].first) covariancefit[ipar][jpar]=covariancefit[ipar][jpar]*normalization;
						if(jpar==normalizedparam[ilow].first) covariancefit[ipar][jpar]=covariancefit[ipar][jpar]*normalization;
					}
	  
				}
	
			} // loop parameter j
      
			if(isfixedisredominimization[ipar].second==true) icov++; // we loop on cov on parameter at a time
      
		} // loop parameter i
    
		// Converting tev to erg for EnergyFlux hypothesis
		if(hypo.GetSpectralType()==Hypothesis::EnergyFlux) {

			for(unsigned int ipar(0); ipar<covariancefit.size(); ipar++) {
				for(unsigned int jpar(0); jpar<covariancefit.size(); jpar++) {
					if(ipar==0) covariancefit[ipar][jpar]=covariancefit[ipar][jpar]*STARTUtils::GetTeVToErg();
					if(jpar==0) covariancefit[ipar][jpar]=covariancefit[ipar][jpar]*STARTUtils::GetTeVToErg();
				}
			}
      
		}

		hypo.SetFittedCovarianceMatrix(covariancefit);

		// Get the time of the minimization 
		INFO << "Time to make the minimization : " 
			 << TimeFromROOT.RealTime() 
			 << "s" << std::endl;
    
	}
	else if(!arefittedparametersvalid && !isminimumvalid && !iscovariancevalid) {
		WARNING << "MINUIT2 has not converged!" << std::endl;
		WARNING << "Fitted parameters are not valid!" << std::endl;
		WARNING << "Function minimum is wrong!" << std::endl;
		WARNING << "Matrix covariance is wrong!" << std::endl;
		INFO_OUT << "Change something (number of bins or intial values for example)" << std::endl;
		paramfit.assign(hypo.GetParametersNb(),0.);
		paramfiterr.assign(hypo.GetParametersNb(),0.);
		covariancefit.resize(hypo.GetParametersNb(),std::vector<double>(hypo.GetParametersNb(),0.));
		hypo.SetFittedParameters(paramfit);
		hypo.SetFittedParametersErrors(paramfiterr);
		hypo.SetFittedCovarianceMatrix(covariancefit);
		hypo.SetIsMinimumValid(isminimumvalid);
		hypo.SetIsCovarianceValid(iscovariancevalid);
		hypo.SetAreFittedParametersValid(arefittedparametersvalid);
		hypo.SetConvergence(false);
	}
	else if(arefittedparametersvalid==false) {
		WARNING << "MINUIT2 has not converged!" << std::endl;
		WARNING << "Fitted parameters are not valid!" << std::endl;
		INFO_OUT << "Change something (number of bins or intial values for example)" << std::endl;
		paramfit.assign(hypo.GetParametersNb(),0.);
		paramfiterr.assign(hypo.GetParametersNb(),0.);
		covariancefit.resize(hypo.GetParametersNb(),std::vector<double>(hypo.GetParametersNb(),0.));
		hypo.SetIsMinimumValid(isminimumvalid);
		hypo.SetIsCovarianceValid(iscovariancevalid);
		hypo.SetAreFittedParametersValid(arefittedparametersvalid);
		hypo.SetConvergence(false);
	}
	else if(isminimumvalid==false){
    
		WARNING << "MINUIT2 has not converged!" << std::endl;
		WARNING << "Function minimum is wrong!" << std::endl;
		INFO_OUT << "Change something (number of bins or intial values for example)" << std::endl;
    
		/* get the fitted values and errors */
		paramfit = MinimumFromMINUIT2.UserState().Params();
		paramfiterr = MinimumFromMINUIT2.UserState().Errors();
    
		for(unsigned int ilow(0); ilow<normalizedparam.size(); ilow++) {
			paramfit[normalizedparam[ilow].first] = paramfit[normalizedparam[ilow].first]*normalization;
			paramfiterr[normalizedparam[ilow].first] = paramfiterr[normalizedparam[ilow].first]*normalization;
			DEBUG_OUT << "value = " << MinimumFromMINUIT2.UserState().Params()[normalizedparam[ilow].first] << std::endl;
		}
    
		hypo.SetFittedParameters(paramfit);
		hypo.SetFittedParametersErrors(paramfiterr);

		covariancefit.resize(hypo.GetParametersNb(),std::vector<double>(hypo.GetParametersNb(),0.));
		hypo.SetFittedCovarianceMatrix(covariancefit);
    
		/* Get the time of the minimization */
		INFO << "Time to make the minimization : " 
			 << TimeFromROOT.RealTime() 
			 << "s" << std::endl;

	}
	else if(iscovariancevalid==false){
		WARNING << "MINUIT2 has not converged!" << std::endl;
		WARNING << "Matrix covariance is wrong!" << std::endl;
		INFO_OUT << "Change something (number of bins or intial values for example)" << std::endl;
    
		/* get the fitted values and errors */
		paramfit = MinimumFromMINUIT2.UserState().Params();
		paramfiterr = MinimumFromMINUIT2.UserState().Errors();
    
		// denormalization
		for(unsigned int ilow(0); ilow<normalizedparam.size(); ilow++) {
			paramfit[normalizedparam[ilow].first] = paramfit[normalizedparam[ilow].first]*normalization;
			paramfiterr[normalizedparam[ilow].first] = paramfiterr[normalizedparam[ilow].first]*normalization;
			DEBUG_OUT << "value = " << MinimumFromMINUIT2.UserState().Params()[normalizedparam[ilow].first] << std::endl;
		}

		hypo.SetFittedParameters(paramfit);
		hypo.SetFittedParametersErrors(paramfiterr);
    
		covariancefit.resize(hypo.GetParametersNb(),std::vector<double>(hypo.GetParametersNb(),0.));
		hypo.SetFittedCovarianceMatrix(covariancefit);

		/* Get the time of the minimization */
		INFO << "Time to make the minimization : " 
			 << TimeFromROOT.RealTime() 
			 << " seconds" << std::endl;
    
	}

    PrintMinimizationCarateristics(hypo); // caracteristics
    if(arefittedparametersvalid) PrintFittedParametersWithErrors(hypo); // parameters + flux
    if(iscovariancevalid) PrintCovarianceMatrix(hypo); // covariance
    if(hypo.GetConvergence() && hypo.GetSpectralType()==Hypothesis::Differential 
       && hypo.GetParametersNb()>1) // JLK hack pour la comparaison, pardon... 
		PrintDecorrelationEnergy(hypo); // decorrelation energy + flux

}



/**
 * \brief Add fitted parameters in hypothesis and print results
 * \param hypo hypothesis
 * \param CovFromMINUIT2 covriance from MINUIT2
 * \param MinimumFromMINUIT2 minimum from MINUIT2
 * \param TimeFromROOT TStopwatch class containing time
 * \param isfixedisredominimization pair of vector containing two boolean (isfixedparam,redominimization)
 */
void START::MinimizeFactory::PrintResultsAndAddInfosInHypothesisWithoutFitting(Hypothesis &hypo,
																			   const ROOT::Minuit2::MnMinimize &MigradSimplex,
																			   TStopwatch &TimeFromROOT)
{

	bool isminimumvalid = true;
	bool iscovariancevalid = true;
	bool arefittedparametersvalid = true;
	bool convergence(true);

	double fcnminimum = 0; //MinimumFromMINUIT2.Fval();
	double edmvalue = 1e-5; //MinimumFromMINUIT2.Edm();
	double iteration = 1; //MinimumFromMINUIT2.NFcn();
  
	frequiredEDM = 0.001*0.5*ftolerance;
  
	std::vector<double> paramfit;
	std::vector<double> paramfiterr;
	std::vector<std::vector<double> > covariancefit;
	//  INFO << "frequiredEDM = " << frequiredEDM << std::endl;
	INFO << "Results of the minimization for " << hypo.GetName() << " (Eref=" << hypo.GetEref() << " TeV";

	if(hypo.GetSpectralType()==Hypothesis::Integrated || hypo.GetSpectralType()==Hypothesis::EnergyFlux) {
		std::cout << ", integration range is between " 
				  << hypo.GetIntegratedFitEnergyRange().first
				  << " and "
				  << hypo.GetIntegratedFitEnergyRange().second 
				  << " TeV";
	}

	std::cout << ") : " << std::endl;
  
	std::vector<std::pair<unsigned int,double> > normalizedparam = 
		hypo.GetVectorNormalizedParameters();

	const double normalization = STARTUtils::GetNormalizationConstante();

	std::vector<std::pair<bool,bool> > isfixedisredominimization = GetVectorPairFixedParamRedoMinimization(hypo);
  
	/* set caracteristics of the minimization */
	hypo.SetMaximumLikelihood(-fcnminimum);
	hypo.SetEDM(edmvalue);
	hypo.SetIteration(iteration);
	hypo.SetIsMinimumValid(isminimumvalid);
	hypo.SetIsCovarianceValid(iscovariancevalid);
	hypo.SetAreFittedParametersValid(arefittedparametersvalid);
	hypo.SetConvergence(convergence);

	INFO << "Using MINUIT2 without fitting !" << std::endl;
    
    /* set fitted parameters and their errors */
    paramfit = MigradSimplex.Params();
    paramfiterr = MigradSimplex.Errors();
    
    /* Denormalization of the parameters */
    for(unsigned int ilow(0); ilow<normalizedparam.size(); ilow++) {
		paramfit[normalizedparam[ilow].first] = paramfit[normalizedparam[ilow].first]*normalization;
		paramfiterr[normalizedparam[ilow].first] = paramfiterr[normalizedparam[ilow].first]*normalization;
    }
    
    if(hypo.GetSpectralType()==Hypothesis::EnergyFlux) {
		paramfit[0]=paramfit[0]*STARTUtils::GetTeVToErg();
		paramfiterr[0]=paramfiterr[0]*STARTUtils::GetTeVToErg();
    }

    hypo.SetFittedParameters(paramfit);
    hypo.SetFittedParametersErrors(paramfiterr);
    
    covariancefit.resize(hypo.GetParametersNb(),std::vector<double>(hypo.GetParametersNb(),0.));
    
    unsigned int icov(0); // covariance indice
    
    hypo.SetFittedCovarianceMatrix(covariancefit);
    
    // Get the time of the minimization 
    INFO << "Time to make the minimization : " 
		 << TimeFromROOT.RealTime() 
		 << "s" << std::endl;
    
}
/**
 * \brief Compute errors using Minos
 */
void START::MinimizeFactory::ComputeMinosErrors(const std::vector<Band> &BandArray, Hypothesis &hypo, ROOT::Minuit2::FunctionMinimum &min) {

        FCNLikelihood FCN(BandArray,hypo,*fConfig, fverbose);

	ROOT::Minuit2::MnMinos MinosErrors(FCN,min,fstrategy);

	//INFO_OUT << "Starting Minos analysis :" << std::endl;

	std::vector<std::pair<double,double> > minoserrors;

	for(unsigned ipar(0); ipar<hypo.GetParametersNb(); ipar++) {
		minoserrors.push_back(MinosErrors(ipar));
	}

	// denormalization
	std::vector<std::pair<unsigned int,double> > normalizedparam = 
		hypo.GetVectorNormalizedParameters();
	double normalization = STARTUtils::GetNormalizationConstante();
	for(unsigned int ilow(0); ilow<normalizedparam.size(); ilow++) {
		minoserrors[normalizedparam[ilow].first].first = minoserrors[normalizedparam[ilow].first].first*normalization;
		minoserrors[normalizedparam[ilow].first].second = minoserrors[normalizedparam[ilow].first].second*normalization;
	}

	INFO << "Parameters with Minos errors :" << std::endl;
  
	for(unsigned int i(0); i<hypo.GetParametersNb(); i++) {
		std::cout << hypo.GetParametersNames()[i] << " = ("  << hypo.GetFittedParameters()[i] 
				  << " - " << -minoserrors[i].first 
				  << " + " << minoserrors[i].second 
				  << ") " << hypo.GetParametersUnits()[i] << std::endl;
	}

	hypo.SetMinosErrors(minoserrors);

}

/**
 * \brief Compute integrated flux
 */
void START::MinimizeFactory::ComputeIntegratedFlux(const std::vector<Band> &BandArray, Hypothesis &hypo) {

	if(hypo.GetConvergence()) {
		std::cout << "****************************************************************" << std::endl;
	}
	else {
		WARNING << hypo.GetName() << " didn't converge so we won't minimize Integrated Hypothesis!" << std::endl;
		return;
	}

	Hypothesis *hypoIF = hypo.clone();

	std::string hypoIFname = hypo.GetName();
	hypoIFname+="_Integrated";

	hypoIF->SetName(hypoIFname.c_str());
	hypoIF->SetSpectralType(Hypothesis::Integrated);
	hypoIF->InitParametersCaracteristics();
	hypoIF->SetParameters(hypo.GetFittedParameters());
	hypoIF->SetParametersErrors(hypo.GetParametersErrors());
	hypoIF->GetFixedParameters().clear();

	ROOT::Minuit2::MnStrategy Strategy(fstrategy); // set the strategy 0,1,2
  
	FCNLikelihood FCN(BandArray,*hypoIF,*fConfig,fverbose);

	ROOT::Minuit2::MnUserParameterState Parameters = GetUserParameterForMinuit2(*hypoIF);

	TStopwatch TimeFromROOT;
	TimeFromROOT.Start();

	ROOT::Minuit2::MnMinimize MigradSimplex(FCN,Parameters,Strategy);

  
	// Fix parameters for the previous fit
	std::vector<std::pair<unsigned int,double> > normalizedparam = hypoIF->GetVectorNormalizedParameters();

	for(unsigned int ipar(0); ipar<hypoIF->GetParametersNb(); ipar++) {
		bool normpar(false);
		for(std::vector<std::pair<unsigned int,double> >::const_iterator it=normalizedparam.begin(); it!=normalizedparam.end(); ++it) {
			if(it->first==ipar) { 
				normpar=true;
				break;
			}
		}
		if(!normpar) hypoIF->FixParameter(ipar,hypo.GetFittedParameters()[ipar],false);
	}

	std::vector<std::pair<unsigned int, std::pair<double, bool> > > fixedparameters = hypoIF->GetFixedParameters();
  
	if(fixedparameters.size()>0) {
    
		for(unsigned int ifix(0); ifix<fixedparameters.size(); ifix++) { // We fix the parameters
      
			unsigned int fixedparam = fixedparameters[ifix].first; 
			double fixedvalue = fixedparameters[ifix].second.first;
      
			MigradSimplex.SetValue(fixedparam,fixedvalue);
			MigradSimplex.Fix(fixedparam);

		}
	}
  
	std::vector<std::pair<unsigned int, std::pair<double, double> > > limitedparameters = hypoIF->GetLimitedParameters();
  
	if(limitedparameters.size()>0) {
    
		for(unsigned int ilimit(0); ilimit<limitedparameters.size(); ilimit++) { // We fix the parameters
      
			unsigned int limitedparam = limitedparameters[ilimit].first; 
			double limitedvaluemin = limitedparameters[ilimit].second.first;
			double limitedvaluemax = limitedparameters[ilimit].second.second;
      
			MigradSimplex.SetLimits(limitedparam,limitedvaluemin, limitedvaluemax);

		}
	}
  

	ROOT::Minuit2::FunctionMinimum MinimumFromMINUIT2 = MigradSimplex(fmaxcall,ftolerance);

	ROOT::Minuit2::MnUserCovariance CovFromMINUIT2 = MinimumFromMINUIT2.UserCovariance();

	TimeFromROOT.Stop();

	PrintResultsAndAddInfosInHypothesis(*hypoIF,MinimumFromMINUIT2,CovFromMINUIT2,TimeFromROOT);

	if(hypoIF->GetConvergence())
		hypo.SetFittedIntegratedFlux(hypoIF->GetFittedParameters()[0],hypoIF->GetFittedParametersErrors()[0]);

	delete hypoIF;
	hypoIF=0;

}

/**
 * \brief Compute integrated energy flux
 */
void START::MinimizeFactory::ComputeEnergyFlux(const std::vector<Band> &BandArray, Hypothesis &hypo) {

	if(hypo.GetConvergence()) {
		std::cout << "****************************************************************" << std::endl;
	}
	else {
		WARNING << hypo.GetName() << " didn't converge so we won't minimize EnergyFlux Hypothesis!" << std::endl;
		return;
	}

	Hypothesis *hypoEF = hypo.clone();

	std::string hypoEFname = hypo.GetName();
	hypoEFname+="_EnergyFlux";

	hypoEF->SetName(hypoEFname.c_str());
	hypoEF->SetSpectralType(Hypothesis::EnergyFlux);
	hypoEF->InitParametersCaracteristics();
	hypoEF->SetParameters(hypo.GetFittedParameters());
	hypoEF->SetParametersErrors(hypo.GetParametersErrors());
	hypoEF->GetFixedParameters().clear();
	hypoEF->GetLimitedParameters().clear();

	ROOT::Minuit2::MnStrategy Strategy(fstrategy); // set the strategy 0,1,2
  
	FCNLikelihood FCN(BandArray,*hypoEF,*fConfig,fverbose);

	ROOT::Minuit2::MnUserParameterState Parameters = GetUserParameterForMinuit2(*hypoEF);

	TStopwatch TimeFromROOT;
	TimeFromROOT.Start();

	ROOT::Minuit2::MnMinimize MigradSimplex(FCN,Parameters,Strategy);

	// Fix parameters from the previous fit
	std::vector<std::pair<unsigned int,double> > normalizedparam = hypoEF->GetVectorNormalizedParameters();

	for(unsigned int ipar(0); ipar<hypoEF->GetParametersNb(); ipar++) {
		bool normpar(false);
		for(std::vector<std::pair<unsigned int,double> >::const_iterator it=normalizedparam.begin(); it!=normalizedparam.end(); ++it) {
			if(it->first==ipar) { 
				normpar=true;
				break;
			}
		}
		if(!normpar) hypoEF->FixParameter(ipar,hypo.GetFittedParameters()[ipar],false);
	}

	std::vector<std::pair<unsigned int, std::pair<double, bool> > > fixedparameters = hypoEF->GetFixedParameters();
  
	if(fixedparameters.size()>0) {
    
		for(unsigned int ifix(0); ifix<fixedparameters.size(); ifix++) { // We fix the parameters
      
			unsigned int fixedparam = fixedparameters[ifix].first; 
			double fixedvalue = fixedparameters[ifix].second.first;
      
			MigradSimplex.SetValue(fixedparam,fixedvalue);
			MigradSimplex.Fix(fixedparam);

		}
	}

	std::vector<std::pair<unsigned int, std::pair<double, double> > > limitedparameters = hypoEF->GetLimitedParameters();
  
	if(limitedparameters.size()>0) {
    
		for(unsigned int ilimit(0); ilimit<limitedparameters.size(); ilimit++) { // We fix the parameters
      
			unsigned int limitedparam = limitedparameters[ilimit].first; 
			double limitedvaluemin = limitedparameters[ilimit].second.first;
			double limitedvaluemax = limitedparameters[ilimit].second.second;
      
			MigradSimplex.SetLimits(limitedparam,limitedvaluemin, limitedvaluemax);

		}
	}

	//

	ROOT::Minuit2::FunctionMinimum MinimumFromMINUIT2 = MigradSimplex(fmaxcall,ftolerance);

	ROOT::Minuit2::MnUserCovariance CovFromMINUIT2 = MinimumFromMINUIT2.UserCovariance();

	TimeFromROOT.Stop();

	PrintResultsAndAddInfosInHypothesis(*hypoEF,MinimumFromMINUIT2,CovFromMINUIT2,TimeFromROOT);
	if(hypoEF->GetConvergence())
		hypo.SetFittedEnergyFlux(hypoEF->GetFittedParameters()[0],hypoEF->GetFittedParametersErrors()[0]);

	delete hypoEF;
	hypoEF=0;

}

/**
 * \brief Print Likelihood ratios for all hypothesis stored in vectors
 */
void START::MinimizeFactory::ComputeAndPrintLikelihoodRatio() {

	if(fPairHypothesisFunctionMinimum.size()>1) {

		INFO << "Comparison between hypothesis : " << std::endl;
    
		for(std::vector<std::pair<Hypothesis*,ROOT::Minuit2::FunctionMinimum*> >::iterator hypo1=fPairHypothesisFunctionMinimum.begin(); 
			hypo1!=fPairHypothesisFunctionMinimum.end(); hypo1++) { // loop on hypothesis
      
			for(std::vector<std::pair<Hypothesis*,ROOT::Minuit2::FunctionMinimum*> >::iterator hypo2=fPairHypothesisFunctionMinimum.begin(); 
				hypo2!=fPairHypothesisFunctionMinimum.end(); hypo2++) { // loop on hypothesis
	
				if((*hypo1).first->GetSpectralType()!=(*hypo2).first->GetSpectralType()) continue;

				if((*hypo1).first->GetConvergence() && (*hypo2).first->GetConvergence()) {
	  
					if((*hypo1).first->GetName()==(*hypo2).first->GetName()) continue; // skip same hypothesis
	  
					std::cout << "   lambda=-2Log(L_{" << (*hypo1).first->GetName() << "}/L_{" << (*hypo2).first->GetName() << "}) = "
							  << -2*((*hypo1).first->GetMaximumLikelihood()-(*hypo2).first->GetMaximumLikelihood())
							  << std::endl;
	  
				}
				else WARNING << "Can't compute the ratio of Likelihood : "
							 << "L{" << (*hypo1).first->GetName() << "/" 
							 << (*hypo2).first->GetName() << "}"
							 << std::endl;
	
			}
      
		}
    
	}
  
}

/**
 * \brief Make countours for two parameters at 1, 2 and 3 sigma
 * \param BandArray data
 * \param hyponame name of the hypothesis
 * \param param1 parameter1 of contour interval
 * \param param2 parameter2 of contour interval
 * \param npoints number of contour's points
 *
 */
void START::MinimizeFactory::MakeContours(std::vector<Band> const &BandArray,std::string hyponame, unsigned int param1, 
										  unsigned int param2,unsigned int npoints)
{
	for(std::vector<std::pair<Hypothesis*,ROOT::Minuit2::FunctionMinimum*> >::iterator hypo=fPairHypothesisFunctionMinimum.begin();
		hypo!=fPairHypothesisFunctionMinimum.end(); hypo++) {

		if((*hypo).first->GetName()!=hyponame) continue; // we are only interested in hyponame hypothesis

		if((*hypo).first->GetConvergence()) {

			std::vector<std::string> paramnames = (*hypo).first->GetParametersNames();

			INFO << "Starting computing contours " << paramnames[param1] 
				 << "/"  << paramnames[param2] << " for hypothesis " 
				 << (*hypo).first->GetName() << "!" << std::endl;

			ComputeContours(BandArray,*((*hypo).first),*((*hypo).second),param1,param2,npoints);

			(*hypo).first->GetContourSigma1().clear();
			(*hypo).first->GetContourSigma2().clear();
			(*hypo).first->GetContourSigma3().clear();

			(*hypo).first->SetContourSigma1(fMapContours1Sigma[(*hypo).first->GetName()]);
			(*hypo).first->SetContourSigma2(fMapContours2Sigma[(*hypo).first->GetName()]);
			(*hypo).first->SetContourSigma3(fMapContours3Sigma[(*hypo).first->GetName()]);

			INFO << "Computing contours for " << (*hypo).first->GetName() << "... ok" << std::endl;

		}
		else {
			WARNING << "Skip hypothesis " << (*hypo).first->GetName() << " because it didn't converge" << std::endl;
			continue;
		}

	} // loop on hypothesis

}

/**
 * \brief Make countours for two parameters at 1, 2 and 3 sigma
 *
 * \param BandArray data
 * \param hyponame name of the hypothesis
 * \param param1 parameter1 of contour interval
 * \param param2 parameter2 of contour interval
 * \param npoints number of contour's points
 *
 */
void START::MinimizeFactory::MakeContours(std::vector<Band> const &BandArray,std::string hyponame, std::string param1, 
										  std::string param2,unsigned int npoints) {

	for(std::vector<std::pair<Hypothesis*,ROOT::Minuit2::FunctionMinimum*> >::iterator hypo=fPairHypothesisFunctionMinimum.begin();
		hypo!=fPairHypothesisFunctionMinimum.end(); hypo++) {
    
		if((*hypo).first->GetName()!=hyponame) continue; // we are only interested in hyponame hypothesis
    
		std::vector<std::string> paramnames = (*hypo).first->GetParametersNames();

		unsigned par1(999), par2(999);
    
		for(unsigned int i(0); i<paramnames.size(); i++) {
			if(paramnames[i]==param1) {
				par1=i;
				break;
			}
		}

		for(unsigned int i(0); i<paramnames.size(); i++) {
			if(paramnames[i]==param2) {
				par2=i;
				break;
			}
		}

		if(par1!=999 && par2!=999)
			MakeContours(BandArray,hyponame,par1,par2,npoints);
		else {
			WARNING << "You ask for contours beteween param " 
					<< param1 << " and " << param2 << ", which probably doesn't exist for hypothesis " 
					<< (*hypo).first->GetName() << "." << std::endl;
			INFO << "Minimization will proceed without it..." << std::endl;
		}
	}

}

/**
 * \brief Make countours for every couples of parameter at 1, 2 and 3 sigma
 * \param BandArray data
 * \param npoints number of contour's points
 *
 * \warning this function may take a very very long time to compute contours
 * so use it only if you are sure about you minimization's results
 *
 */
void START::MinimizeFactory::MakeAllContours(std::vector<Band> const &BandArray, int npoints)
{

	for(std::vector<std::pair<Hypothesis*,ROOT::Minuit2::FunctionMinimum*> >::iterator hypo=fPairHypothesisFunctionMinimum.begin();
		hypo!=fPairHypothesisFunctionMinimum.end(); hypo++) {

		if((*hypo).first->GetConvergence()) {

			std::vector<std::string> paramnames = (*hypo).first->GetParametersNames();

			fcontoursigma1.clear();
			fcontoursigma2.clear();
			fcontoursigma3.clear();
      
			for(unsigned int ipar(0); ipar<(*hypo).first->GetParametersNb(); ipar++) {
	
				for(unsigned int jpar(0); jpar<(*hypo).first->GetParametersNb(); jpar++) {
	  
					if(jpar<=ipar) continue; // n(n-1)/2 contours for n parameters

					INFO << "Starting computing contours " << paramnames[ipar] 
						 << "/"  << paramnames[jpar] << " for hypothesis " 
						 << (*hypo).first->GetName() << "!" << std::endl;

					ComputeContours(BandArray,*((*hypo).first),*((*hypo).second),ipar,jpar,npoints);
				}
			}
 
			(*hypo).first->GetContourSigma1().clear();
			(*hypo).first->GetContourSigma2().clear();
			(*hypo).first->GetContourSigma3().clear();

			(*hypo).first->SetContourSigma1(fMapContours1Sigma[(*hypo).first->GetName()]);
			(*hypo).first->SetContourSigma2(fMapContours2Sigma[(*hypo).first->GetName()]);
			(*hypo).first->SetContourSigma3(fMapContours3Sigma[(*hypo).first->GetName()]);

		}
		else {
			WARNING << "Skip hypothesis " << (*hypo).first->GetName() << " because it didn't converge" << std::endl;
			continue;
		}

	} // loop on hypothesis

	INFO << "Computing contours... ok" << std::endl;

}

/**
 * \brief Compute contours at 1, 2 and 3 sigma for param1 and param2
 */
void START::MinimizeFactory::ComputeContours(const std::vector<Band> &BandArray, Hypothesis &hypothesis, 
											 ROOT::Minuit2::FunctionMinimum &funcminimum,
											 unsigned int param1, unsigned int param2, unsigned int npoints)
{

	DEBUG_OUT << "hypothesis : " << hypothesis.GetName() << std::endl;

	DEBUG_OUT << "Minimum from minuit2 :" << std::endl;
	DEBUG_OUT << funcminimum << std::endl;

	Hypothesis *hypo = &hypothesis;
	ROOT::Minuit2::FunctionMinimum *funcmin = &funcminimum;

	ROOT::Minuit2::MnStrategy Strategy(fstrategy); // set the strategy 0,1,2

	INFO_OUT << "1 sigma" << std::endl; // 1 sigma contours

	std::vector<std::pair<double,double> > contour1;

	FCNLikelihood *FCN = new FCNLikelihood(BandArray,*hypo,*fConfig,fverbose);
	FCN->SetErrorDef(2.3/2.); // 68.27% confidence level

	ROOT::Minuit2::MnContours *Contours = new ROOT::Minuit2::MnContours(*FCN,*funcmin,Strategy);
	contour1 = (*Contours)(param1,param2,npoints);
	DEBUG_OUT << "contour1.size()" << fMapContours1Sigma[hypothesis.GetName()].size() << std::endl;

	delete FCN; FCN = 0;
	delete Contours; Contours = 0;

	INFO_OUT << "2 sigma" << std::endl;

	std::vector<std::pair<double,double> > contour2;

	FCN = new FCNLikelihood(BandArray,*hypo,*fConfig,fverbose); 
	FCN->SetErrorDef(6.19/2.); // 95.45% confidence level 
	Contours = new ROOT::Minuit2::MnContours(*FCN,*funcmin,Strategy);
	contour2 = (*Contours)(param1,param2,npoints);
	DEBUG_OUT << "contour2.size()" << fMapContours2Sigma[hypothesis.GetName()].size() << std::endl;

	delete FCN; FCN = 0;
	delete Contours; Contours = 0;

	INFO_OUT << "3 sigma" << std::endl;

	std::vector<std::pair<double,double> > contour3;

	FCN = new FCNLikelihood(BandArray,*hypo,*fConfig,fverbose); 
	FCN->SetErrorDef(11.83/2.); // 99.73% confidence level 
	Contours = new ROOT::Minuit2::MnContours(*FCN,*funcmin,Strategy);
	contour3 = (*Contours)(param1,param2,npoints);
	DEBUG_OUT << "contour3.size()" << fMapContours3Sigma[hypothesis.GetName()].size() << std::endl;

	delete FCN; FCN = 0;
	delete Contours; Contours = 0;

	// JLK : sometimes contours have different size so we have to denormalized for each contours

	std::vector<std::pair<unsigned int,double> > normalizedparam = hypo->GetVectorNormalizedParameters();

	DEBUG_OUT << "Denormalization : " << std::endl;

	const double normalization = STARTUtils::GetNormalizationConstante();

	for(unsigned int icont(0); icont<contour1.size(); icont++) { // denormalization of the parameters
		for(unsigned int ilow(0); ilow<normalizedparam.size(); ilow++) {
			if(param1==normalizedparam[ilow].first) contour1[icont].first=contour1[icont].first*normalization;
			if(param2==normalizedparam[ilow].first) contour1[icont].second=contour1[icont].second*normalization;
		}
	}

	for(unsigned int icont(0); icont<contour2.size(); icont++) { // denormalization of the parameters
		for(unsigned int ilow(0); ilow<normalizedparam.size(); ilow++) {
			if(param1==normalizedparam[ilow].first) contour2[icont].first=contour2[icont].first*normalization;
			if(param2==normalizedparam[ilow].first) contour2[icont].second=contour2[icont].second*normalization;
		}
	}

	for(unsigned int icont(0); icont<contour3.size(); icont++) { // denormalization of the parameters
		for(unsigned int ilow(0); ilow<normalizedparam.size(); ilow++) {
			if(param1==normalizedparam[ilow].first) contour3[icont].first=contour3[icont].first*normalization;
			if(param2==normalizedparam[ilow].first) contour3[icont].second=contour3[icont].second*normalization;
		}
	}

	DEBUG_OUT << "Copy contours into maps" << std::endl;
	fMapContours1Sigma[hypothesis.GetName()].push_back(std::make_pair(contour1,std::make_pair(param1,param2)));
	fMapContours2Sigma[hypothesis.GetName()].push_back(std::make_pair(contour2,std::make_pair(param1,param2)));
	fMapContours3Sigma[hypothesis.GetName()].push_back(std::make_pair(contour3,std::make_pair(param1,param2)));

}

/**
 * \brief Scan likelihood domain for one parameter (other parameters fixed)
 *
 * \param hypo hypothesis
 * \param BandArray data
 * \param param parameter of scans
 * \param low lower limit of param
 * \param high higher limit of param
 * \param npoints number of contour's points
 *
 * If low and high are not specified, scan will be done at 2 sigma from the minimum
 */
void START::MinimizeFactory::ScanLikelihoodAfterMinimization(Hypothesis &hypo,std::vector<Band> &BandArray,unsigned param,
															 unsigned int npoints, double low,double high) {

	hypo.GetScansLikelihood().clear();
 
	FCNLikelihood *FCN = new FCNLikelihood(BandArray,hypo,*fConfig,fverbose);

	const double normalization = STARTUtils::GetNormalizationConstante();
    std::vector<std::pair<unsigned int,double> > normalizedparam = hypo.GetVectorNormalizedParameters();

	std::vector<double> fittedparams = hypo.GetFittedParameters();
	std::vector<double> fittedparamserrors = hypo.GetFittedParametersErrors();

	for(unsigned int ipar(0); ipar<fittedparams.size(); ipar++) {
		for(unsigned int ilow(0); ilow<normalizedparam.size(); ilow++) { //loop on the parameter which have a lower limit
			if(normalizedparam[ilow].first==ipar) { // if the parameter has a lower limit we normalize
				fittedparams[ipar]=fittedparams[ipar]/normalization;
				fittedparamserrors[ipar]=fittedparamserrors[ipar]/normalization;
			}
		}
	}

	ROOT::Minuit2::MnScan LikelihoodBuilder(*FCN,fittedparams,fittedparamserrors,fstrategy);
  
	std::vector<std::pair<double,double> > scan;

	if(low==-999 || high==-999) {
		INFO << "Scanning likelihood region between for parameter " << hypo.GetParametersNames()[param] <<  std::endl;
		scan = LikelihoodBuilder.Scan(param,npoints);
	}
	else {
		INFO << "Scanning likelihood region between " << low << " and " << high << " for parameter " << param <<  std::endl;
		scan = LikelihoodBuilder.Scan(param,npoints,low,high);
	}

	// de-normalization
	for(std::vector<std::pair<unsigned int,double> >::const_iterator norm=normalizedparam.begin(); 
		norm!=normalizedparam.end(); ++norm) {
		if((*norm).first==param) {
			for(std::vector<std::pair<double,double> >::iterator tonorm=scan.begin(); tonorm!=scan.end();
				++tonorm) {
				(*tonorm).first*=STARTUtils::GetNormalizationConstante();
			}
		}
	}

	fMapScanLikeliHood[hypo.GetName()].push_back(std::make_pair(param,scan));

	for(unsigned int ipoint(0); ipoint<fMapScanLikeliHood[hypo.GetName()].back().second.size(); ipoint++) {
		DEBUG_OUT << "point " << " 1 : " << fMapScanLikeliHood[hypo.GetName()].back().second[ipoint].first << " " 
				  << fMapScanLikeliHood[hypo.GetName()].back().second[ipoint].second << std::endl;
	}

	DEBUG_OUT << "number of scans" << fMapScanLikeliHood[hypo.GetName()].size() << std::endl;
  
	hypo.SetScansLikelihood(fMapScanLikeliHood[hypo.GetName()]);
	DEBUG_OUT << "number of scans hypo" << hypo.GetScansLikelihood().size() << std::endl;
	delete FCN; FCN=0;
  
}

/**
 * \brief Scan likelihood domain for one parameter (other parameters fixed)
 *
 * \param hypo hypothesis
 * \param BandArray data
 * \param param parameter of scans
 * \param low lower limit of param
 * \param high higher limit of param
 * \param npoints number of contour's points
 *
 * If low and high are not specified, scan will be done at 2 sigma from the minimum
 */
void START::MinimizeFactory::ScanLikelihoodAfterMinimization(Hypothesis &hypo,std::vector<Band> &BandArray,std::string param,
															 unsigned int npoints, double low,double high) {

	std::vector<std::string> paramnames = hypo.GetParametersNames();
  
	unsigned par(999);
  
	for(unsigned int i(0); i<paramnames.size(); i++) {
		if(paramnames[i]==param) {
			par=i;
			break;
		}
	}

	if(par!=999)
		ScanLikelihoodAfterMinimization(hypo,BandArray,par,npoints,low,high);
	else {
		WARNING << "You ask for scanning likelihood region beteween " 
				<< low << " and " << high << " for param " << param <<  ", which probably doesn't exist for hypothesis " 
				<< hypo.GetName() << "." << std::endl;
		INFO_OUT << "Minimization will proceed without it..." << std::endl;
	}
  
}

/**
 * \brief Scan likelihood domain for all parameters and all hypothesis
 *
 * \param BandArray data
 * \param npoints number of scan's points
 *
 * Scan will be done at 2 sigma from the minimum
 */
void START::MinimizeFactory::ScanLikelihoodForAllHypothesisAfterMinimization(std::vector<Band> &BandArray, unsigned int npoints) {
    
	for(std::vector<std::pair<Hypothesis*,ROOT::Minuit2::FunctionMinimum*> >::iterator hypo1=fPairHypothesisFunctionMinimum.begin(); 
		hypo1!=fPairHypothesisFunctionMinimum.end(); hypo1++) { // loop on hypothesis
    
		INFO << "Scanning hypothesis " << (*hypo1).first->GetName() << "!" << std::endl;      
    
		for(unsigned int ipar(0); ipar<(*hypo1).first->GetParametersNb(); ipar++) {

			ScanLikelihoodAfterMinimization(*(*hypo1).first,BandArray,ipar,npoints);

		}

      
	}
  

}

/**
 * \brief Print minimization caracteristics such as edm, iterations' number and minimum
 * for current hypothesis
 */
void START::MinimizeFactory::PrintMinimizationCarateristics(Hypothesis &hypo)
{

	INFO << "FCN : " << std::endl;
	std::cout << "   maximum likelihood (log) = "  << std::setprecision(12)<< hypo.GetMaximumLikelihood() << std::endl;
	std::cout << "   EDM = " << hypo.GetEDM() << " (required EDM = " << frequiredEDM << ")" <<  std::endl;
	std::cout << "   number of iterations = " << hypo.GetIteration() << std::endl;
}

/**
 * \brief Print fitted parameters with their errors for current hypothesis
 */
void START::MinimizeFactory::PrintFittedParametersWithErrors(Hypothesis &hypo)
{

	INFO << "Parameters : " << std::endl;

	for(unsigned int i(0); i<hypo.GetParametersNb(); i++) {
		std::cout << "   " << hypo.GetParametersNames()[i] << " = ("  << hypo.GetFittedParameters()[i] 
				  << " +/- " << hypo.GetFittedParametersErrors()[i] 
				  << ") " << hypo.GetParametersUnits()[i] << std::endl;
	}

	double eref = 1.;
	std::cout << "   phi(" << eref << " TeV) = (" << hypo.GetFluxFitParams(eref)
			  << " +/- " << hypo.GetSigmaFlux(eref) << ") " 
			  << hypo.GetFluxUnits() << std::endl; 

}

/**
 * \brief Print covariance matrix for current hypothesis
 */
void START::MinimizeFactory::PrintCovarianceMatrix(Hypothesis &hypo)
{

	INFO << "Covariance matrix : " << std::endl;

	for(unsigned int i(0); i<hypo.GetParametersNb(); i++) {
		std::cout << "   ";
		for(unsigned int j(0); j<hypo.GetParametersNb(); j++) {
			std::cout << hypo.GetFittedCovarianceMatrix()[i][j] <<  "  ";
		}
		std::cout << std::endl;
	}

}

/**
 * \brief Print covariance matrix
 */
void START::MinimizeFactory::PrintDecorrelationEnergy(Hypothesis &hypo)
{

	INFO << "Decorrelation energy : " << std::endl;
	double edecclin = hypo.GetLinearDecorrelationEnergy();
	INFO_OUT << "Decorelation energy found from first order flux covariance : " << std::endl;
	std::cout << "   Ed = " << edecclin << " TeV" << std::endl;
	std::cout << "   phi(" << edecclin << " TeV) = (" << hypo.GetFluxFitParams(edecclin)
			  << " +/- " << hypo.GetSigmaFlux(edecclin) << ") " 
			  << hypo.GetFluxUnits() << std::endl;
	double edeclog = hypo.GetLogarithmDecorrelationEnergy();
	INFO_OUT << "Decorelation energy found from first order log flux covariance : " << std::endl;
	std::cout << "   Ed = " << edeclog << " TeV" << std::endl;
	std::cout << "   phi(" << edeclog << " TeV) = (" << hypo.GetFluxFitParams(edeclog)
			  << " +/- " << hypo.GetSigmaFlux(edeclog) << ") " 
			  << hypo.GetFluxUnits() << std::endl;
}

/**
 * \brief Print covariance, fitted parameters and minimization's caracteristics
 * for all hypothesis
 */
void START::MinimizeFactory::PrintMinimizationSummary()
{

}

/**
 * \brief Loop on fitted parameters and detects if there is NaN. If there is, set to false 
 * hypothesis's parameters which reveals the good quality of the fit
 */
bool START::MinimizeFactory::IsThereAnyNaNInFittedParameters(std::vector<double> fittedparams, std::vector<double> fittedparamerrors)
{

	for(unsigned int ipar(0); ipar<fittedparams.size(); ipar++) {
		if(TMath::IsNaN(fittedparams[ipar])) return true;
	}

	for(unsigned int ipar(0); ipar<fittedparamerrors.size(); ipar++) {
		if(TMath::IsNaN(fittedparamerrors[ipar])) return true;
	}

	return false;
}

/**
 *
 * \brief Store informations in bands vector member of hypothesis :
 * </ul>
 * </li> expected signal in each band and each energy bin <li/>
 * </li> mean energy of each energy bin given the fitted law <li/>
 * </li> theoretical differential flux (dN/dE) for the energy bin mean energy <li/>
 * <ul/>
 */
void START::MinimizeFactory::CopyExcessInBandsHypothesis(Hypothesis *hypo, std::vector<Band> const &BandArray)
{
	ComputeResults *CompRes = NULL;
  
	if(hypo->GetConvergence()) {
		CompRes = new ComputeResults(BandArray,*hypo);
    
		std::vector<Band> BandArray_copy = BandArray;
    
		//Check if BandArray exists
		if (!(BandArray_copy.size()>0)) {
			WARN_OUT <<"BandArray is empty. Exit StoreExpectedExcess"<<std::endl;
			exit(EXIT_FAILURE);
		}
    
		// Clean the Sth in bins
		for(std::vector<Band>::iterator band = BandArray_copy.begin();band!=BandArray_copy.end();++band) {
			for (std::vector<EnergyBin>::iterator bin = band->ebin.begin();bin!=band->ebin.end();++bin) {
				bin->SetSth(0.);
			}
		}
    
		for(std::vector<Band>::iterator band = BandArray_copy.begin();band!=BandArray_copy.end();++band) {
      
			if(band->GetKeepBand()==0) continue; // We are interested only in the selected bands
      
			for (std::vector<EnergyBin>::iterator bin = band->ebin.begin();bin!=band->ebin.end();++bin) {
	
				if (bin->GetKeepBin()==0) continue; // skip energy bins below threshold
	
				// expected excess
				double S = CompRes->FunctionExpectedExcess((&*band),(&*bin),hypo->GetFittedParameters());
				bin->SetSth(S);
				double expected_off = STARTUtils::GetExpectedOff(TMath::Nint(bin->GetOn()),
																 TMath::Nint(bin->GetOff()),
																 band->GetAlphaRun(),S);
				// VIM : Maybe use the Alpha value inside the EnergyBin instead here.
				double expected_on = S + band->GetAlphaRun()*expected_off;

				bin->SetOnFitted(expected_on);
				bin->SetOffFitted(expected_off);
	  
				double emin = bin->GetEmin();
				double emax = bin->GetEmax();
				double emean = hypo->GetMeanBinEnergy(emin,emax);
	
				bin->SetEmean(emean);
	
				DEBUG_OUT_L(2) << "emean=" << emean << " Sth=" << S << " Sexp=" << bin->GetOn() << "exected_on = " 
							   << expected_on << " expected_off = " << expected_off << " noff = " << bin->GetOff() << std::endl;
	
			} //loop on bin
      
		} // loop on band
    
		if(CompRes!=0) delete CompRes;
		CompRes=0;
    
		hypo->SetBandArray(BandArray_copy); // We stock a BandArray into the hypothesis
		DEBUG_OUT_L(2) << "After SetBandArray(BandArray_copy)" << std::endl;    
		DataSummary DataSum(BandArray_copy);
		DEBUG_OUT_L(2) << "After DataSummary DataSum(BandArray_copy);" << std::endl; 
		DataSum.PrintSummaryBand(BandArray_copy,true);

	}
	DEBUG_OUT_L(2) << "CopyExcessInBandsHypothesis Done" << std::endl;    
}

/**
 * \brief Compare hypothesis by minimizing :
 * \f[
 * \phi(E) =  \chi \phi_{H_1}(E) + (1-\chi) \phi_{H_2}(E)
 * \f]
 */
void START::MinimizeFactory::CompareHypothesisByMinimization(std::vector<Band> &BandArray, Hypothesis &hypo1, Hypothesis &hypo2)
{

	if(hypo1.GetSpectralType()!=hypo2.GetSpectralType()) {
		WARNING << "We can't compare hypothesis with different spectral type!" << std::endl;
		return;
	}
	else INFO << "Comparing hypothesis " << hypo1.GetName() << " and " << hypo2.GetName() << "!" << std::endl;

	TString name = "Comparison_";
	name+=hypo1.GetName();
	name+="_";
	name+=hypo2.GetName();

	SumHypothesis HypoComp(name,hypo1,hypo2,true);

	ROOT::Minuit2::MnStrategy Strategy(fstrategy); // set the strategy 0,1,2
  
	FCNLikelihood FCN(BandArray,HypoComp,*fConfig,fverbose);

	ROOT::Minuit2::MnUserParameterState Parameters = GetUserParameterForMinuit2(HypoComp);
	Parameters.SetLimits(0,0.,1.); // we limit x parameter from 0 to 1
	TStopwatch TimeFromROOT;
	TimeFromROOT.Start();

	ROOT::Minuit2::MnMinimize MigradSimplex(FCN,Parameters,Strategy);

	ROOT::Minuit2::FunctionMinimum MinimumFromMINUIT2 = MigradSimplex(fmaxcall,ftolerance);

	ROOT::Minuit2::MnUserCovariance CovFromMINUIT2 = MinimumFromMINUIT2.UserCovariance();

	TimeFromROOT.Stop();

	PrintResultsAndAddInfosInHypothesis(HypoComp,MinimumFromMINUIT2,CovFromMINUIT2,TimeFromROOT);

	std::cout << "****************************************************************" << std::endl;

}

/**
 * \brief Compare all hypothesis by minimizing :
 * \f[
 * \phi(E) =  \chi \phi_{H_1}(E) + (1-\chi) \phi_{H_2}(E)
 * \f]
 *
 * \param BandArray vector of Band
 */
void START::MinimizeFactory::CompareAllHypothesisByMinimization(std::vector<Band> &BandArray)
{

	std::cout << "****************************************************************" << std::endl;

	for(std::vector<std::pair<Hypothesis*,ROOT::Minuit2::FunctionMinimum*> >::iterator hypo1=fPairHypothesisFunctionMinimum.begin(); 
		hypo1!=fPairHypothesisFunctionMinimum.end(); hypo1++) { // loop on hypothesis
    
		for(std::vector<std::pair<Hypothesis*,ROOT::Minuit2::FunctionMinimum*> >::iterator hypo2=fPairHypothesisFunctionMinimum.begin(); 
			hypo2!=fPairHypothesisFunctionMinimum.end(); hypo2++) { // loop on hypothesis
      
			if((*hypo1).first->GetSpectralType()!=(*hypo2).first->GetSpectralType()) continue;
      
			if((*hypo1).first->GetConvergence() && (*hypo2).first->GetConvergence()) {
	
				if((*hypo1).first->GetName()==(*hypo2).first->GetName()) continue;
	
				CompareHypothesisByMinimization(BandArray,*(*hypo1).first,*(*hypo2).first);	

			}

		}

	}

	INFO << "Comparison between hypothesis... ok" << std::endl;

}

/**
 * \brief Compute minimization energy range
 *
 * \param BandArray vector of band
 */
std::pair<double,double> START::MinimizeFactory::ComputeMinimizationEnergyRange(const std::vector<Band> BandArray)
{

	std::vector<double> vemin, vemax;

	for(std::vector<Band>::const_iterator band=BandArray.begin(); band!=BandArray.end(); ++band) {

		if(band->GetKeepBand()==0) continue;

		for(std::vector<EnergyBin>::const_iterator bin=band->ebin.begin(); bin!=band->ebin.end(); ++bin) {

			if(bin->GetKeepBin()==0) continue;

			vemin.push_back(bin->GetEmin());
			vemax.push_back(bin->GetEmax());

		}
    
	}

	std::vector<double>::const_iterator it1, it2;

	it1 = min_element(vemin.begin(),vemin.end());
	it2 = max_element(vemax.begin(),vemax.end());
	double eth(0.), emax(0);

	if(vemin.size()!=0 && vemax.size()!=0) {
		eth = *it1;
		emax = *it2;
	}

	std::pair<double,double> erange = std::make_pair(eth,emax);

	return erange;

}
