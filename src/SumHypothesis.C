// STL
#include <vector>

// ROOT
#include <TMath.h>

// START
#include "SumHypothesis.hh"

// Utilities
#define DEBUG 0
#include "debugging.hh"

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::SumHypothesis)
#endif

/**
 * \brief Default constructor needed by root
 */
START::SumHypothesis::SumHypothesis()
  :Hypothesis()
{

}

/**
 * \brief Constructor to build a new hypothesis which is the sum of two hypothesis
 * or to build a weighted hypothesis build with the two hypothesis
 * \param name hypothesis' name
 * \param hypo1 first hypothesis \f$H_1\f$
 * \param hypo1 second hypothesis \f$H_2\f$
 * \param comparison if false we sum \f$H_1\f$ and \f$H_2\f$ hypothesis : 
 * \f[
 * \phi(E) = \phi_{H_1}(E) + \phi_{H_2}(E)
 * \f]
 * and if not, we build the hypothesis : 
 * \f[
 * \phi(E) =  \chi \phi_{H_1}(E) + (1-\chi) \phi_{H_2}(E)
 * \f]
 * used to compare hypothesis by minimization
 *
 * In case of hypothesis' sum, we build a new hypothesis with parameters
 * ordered from \f$H_1\f$ to \f$H_2\f$, and add the hypothesis' names to their names
 * (e.g. for an hypothesis PowerLaw named PWL, the first parameter will be named "phi0_PWL")
 */
START::SumHypothesis::SumHypothesis(TString name, const Hypothesis &hypo1, const Hypothesis &hypo2, bool comparison)
  :Hypothesis(name)
{

  fiscomparison=comparison;

  fHypo1 = hypo1.clone();
  fHypo2 = hypo2.clone();

  // Check spectral type
  if(hypo1.GetSpectralType()!=hypo2.GetSpectralType()) {
    WARN_OUT << "We can't add or compare hypothesis with different spectral type!" << std::endl;
    INFO_OUT << "You'll have hypothesis with " << fHypo1->GetName() << "spectral type..." << std::endl;
    fHypo2->SetSpectralType(fHypo1->GetSpectralType());
  }

  fSpectralType = fHypo1->GetSpectralType();


  // Check reference energy
  if(fHypo1->GetEref()!=fHypo2->GetEref()) {
    WARN_OUT << "Summing hypothesis with different reference energies!" << std::endl;
    INFO_OUT << fHypo1->GetName() << " : Eref = " << fHypo1->GetEref() << std::endl;
    INFO_OUT << fHypo2->GetName() << " : Eref = " << fHypo2->GetEref() << std::endl;
    INFO_OUT << "We set reference energy at 1. TeV..." << std::endl;
  }

  feref = fHypo1->GetEref();

  // Check integration energy range
  if(fHypo1->GetIntegratedFitEnergyRange().first!=fHypo2->GetIntegratedFitEnergyRange().first ||
     fHypo1->GetIntegratedFitEnergyRange().second!=fHypo2->GetIntegratedFitEnergyRange().second) {
    WARN_OUT << "Summing hypothesis with different integrated energy range!" << std::endl;
    INFO_OUT << fHypo1->GetName() << " : E1 = " << fHypo1->GetIntegratedFitEnergyRange().first 
	     << " and E2 = " << fHypo1->GetIntegratedFitEnergyRange().second << std::endl;
    INFO_OUT << fHypo2->GetName() << " : E1 = " << fHypo2->GetIntegratedFitEnergyRange().first 
	     << " and E2 = " << fHypo2->GetIntegratedFitEnergyRange().second << std::endl;
    INFO_OUT << "We set integrated energy range between "<< fHypo1->GetIntegratedFitEnergyRange().first 
	     << " and " << fHypo1->GetIntegratedFitEnergyRange().second << " TeV..." << std::endl;
    SetIntegratedFitEnergyRange(fHypo1->GetIntegratedFitEnergyRange().first,fHypo1->GetIntegratedFitEnergyRange().second); 
  }

  if(fiscomparison) {
    fnumberofparameters=1; // JLK
    fparam.resize(fnumberofparameters,0.);
    fparamerr.resize(fnumberofparameters,0.);
    fparamname.resize(fnumberofparameters,"");
    flatexparamname.resize(fnumberofparameters,"");
    fcovariancefit.resize(fnumberofparameters,std::vector<double >(fnumberofparameters,0.));
    fparamunits.resize(fnumberofparameters,"");
    ffirstderivatives.resize(fnumberofparameters,0.);
  }
  else {
    fnumberofparameters=fHypo1->GetParametersNb()+fHypo2->GetParametersNb();
  }

  fcovariancefit.resize(fnumberofparameters,std::vector<double >(fnumberofparameters,0.));

  if(fiscomparison) {
    fparam[0] = 0.5;
    fparamerr[0] = 0.;
  }
  else {

    for(unsigned int ipar(0); ipar<fHypo1->GetParametersNb(); ipar++) {
      fparam.push_back(fHypo1->GetParameters()[ipar]);
      fparamerr.push_back(fHypo1->GetParametersErrors()[ipar]);
    }
    for(unsigned int ipar(0); ipar<fHypo2->GetParametersNb(); ipar++) {
      fparam.push_back(fHypo2->GetParameters()[ipar]);
      fparamerr.push_back(fHypo2->GetParametersErrors()[ipar]);
    }

    // copy of normalized parameters
    fnormalizedparam.clear();

    std::vector<std::pair<unsigned int,double> > normparam1 = fHypo1->GetVectorNormalizedParameters();
    std::vector<std::pair<unsigned int,double> > normparam2 = fHypo2->GetVectorNormalizedParameters();

    for(unsigned int inorm(0); inorm<normparam1.size(); inorm++) {
      SetLowerLimitOnNormalizationParameter(normparam1[inorm].first,normparam1[inorm].second);
    }

    for(unsigned int inorm(0); inorm<normparam2.size(); inorm++) {
      SetLowerLimitOnNormalizationParameter(fHypo1->GetParametersNb()+normparam2[inorm].first,normparam2[inorm].second);
    }

    // copy of fixed parameters
    ffixedparameter.clear();
    std::vector<std::pair<unsigned int, std::pair<double, bool> > > fixedparam1 = fHypo1->GetFixedParameters();
    std::vector<std::pair<unsigned int, std::pair<double, bool> > > fixedparam2 = fHypo2->GetFixedParameters();

    for(unsigned int ifix(0); ifix<fixedparam1.size(); ifix++) {
      FixParameter(fixedparam1[ifix].first,fixedparam1[ifix].second.first,fixedparam1[ifix].second.second);
    }

    for(unsigned int ifix(0); ifix<fixedparam2.size(); ifix++) {
      FixParameter(fixedparam2[ifix].first,fixedparam2[ifix].second.first,fixedparam2[ifix].second.second);
    }

  }

  InitParametersCaracteristics();

  flatexformula = fHypo1->GetLatexFormula();
  flatexformula+=" + ";
  flatexformula = fHypo2->GetLatexFormula();

  fmodelname = fHypo1->GetModelName();
  fmodelname = " + ";
  fmodelname = fHypo2->GetModelName();

}

/**
 * \brief Destructor
 */
START::SumHypothesis::~SumHypothesis()
{
  if(fIntegrantFlux!=0) delete fIntegrantFlux;
  fIntegrantFlux=0;
  if(fIntegrantFluxTimesE!=0) delete fIntegrantFluxTimesE;
  fIntegrantFluxTimesE=0;
  if(fFunctorIntegrantFlux!=0) delete fFunctorIntegrantFlux;
  fFunctorIntegrantFlux=0;
  if(fFunctorIntegrantFluxTimesE!=0) delete fFunctorIntegrantFluxTimesE;
  fFunctorIntegrantFluxTimesE=0;

  if(fIntegrantFitFlux!=0) delete fIntegrantFitFlux;
  fIntegrantFitFlux=0;
  if(fIntegrantFitFluxTimesE!=0) delete fIntegrantFitFluxTimesE;
  fIntegrantFitFluxTimesE=0;
  if(fFunctorIntegrantFitFlux!=0) delete fFunctorIntegrantFitFlux;
  fFunctorIntegrantFitFlux=0;
  if(fFunctorIntegrantFitFluxTimesE!=0) delete fFunctorIntegrantFitFluxTimesE;
  fFunctorIntegrantFitFluxTimesE=0;

  /* Butterflies */
  for(std::map<std::string,TPolyLine*>::iterator butt=fMapButterfly.begin(); butt!=fMapButterfly.end(); ++butt) {
    if(butt->second!=0) delete butt->second;
    butt->second = 0;
  }
  fMapButterfly.clear();

  if(fHypo1!=0) delete fHypo1;
  fHypo1=0;
  if(fHypo2!=0) delete fHypo2;
  fHypo2=0;
}

/**
 * \brief Copy Constructor
 */
START::SumHypothesis::SumHypothesis(SumHypothesis const &StoreCopy)
  :Hypothesis(StoreCopy.fName)
{

  /* MINUIT2's parameters */
  fmaximumlikelihood = StoreCopy.fmaximumlikelihood; // maximum of the likelihood sor the fitted parameters
  fedm = StoreCopy.fedm; // normalized distance to the minimum
  fiteration = StoreCopy.fiteration; // number of iteration fro the minimization
  fisminimumvalid = StoreCopy.fisminimumvalid; // if true the minimum is valid
  fiscovariancevalid = StoreCopy.fiscovariancevalid; // if true the covariance is valid
  farefittedparametersvalid = StoreCopy.farefittedparametersvalid; // if true fitted parameters are valid
  fconvergence = StoreCopy.fconvergence;

  /* Parameters */
  fnumberofparameters = StoreCopy.fnumberofparameters; // number of parameters of the hypothesis
  feref = StoreCopy.feref; // reference energy
  fparamunits = StoreCopy.fparamunits; // strings with parameters' units
  // std::string fname; // name of the hypothesis
  fparam = StoreCopy.fparam; // parameters of the hypothesis 
  fparamfit = StoreCopy.fparamfit; // fitted parameters of the hypothesis 
  fparamerr = StoreCopy.fparamerr; // parameters' errors 
  fparamerrfit = StoreCopy.fparamerrfit; // fitted parameters' errors 
  fcovariancefit = StoreCopy.fcovariancefit; // covariance matrix
  fparamname = StoreCopy.fparamname; // parameters' names
  flatexparamname = StoreCopy.flatexparamname; // parameters' names
  fnormalizedparam = StoreCopy.fnormalizedparam; // number of the normalized parameter and the lower limit
  ffixedparameter = StoreCopy.ffixedparameter;
  ffirstderivatives = StoreCopy.ffirstderivatives; // first derivatives of the flux for parameters 0, 1, 2 .. N
  flatexformula = StoreCopy.flatexformula; // latex formula
  fintegratedfitenergyrange = StoreCopy.fintegratedfitenergyrange;
  ffittedintegratedflux = StoreCopy.ffittedintegratedflux;
  ffittedenergyflux = StoreCopy.ffittedenergyflux;
  fmodelname = StoreCopy.fmodelname;
  fminimizationenergyrange = StoreCopy.fminimizationenergyrange;
  flineardecorrelationenergy = StoreCopy.flineardecorrelationenergy; // linear decorrelation energy
  flogarithmdecorrelationenergy = StoreCopy.flogarithmdecorrelationenergy; // logarithm decorrelation energy
  fcontoursdecorrelationenergy = StoreCopy.fcontoursdecorrelationenergy; // contours decorrelation energy

  /* Contours */
  fcontoursigma1 = StoreCopy.fcontoursigma1; // contours sigma1
  fcontoursigma2 = StoreCopy.fcontoursigma2; // contours sigma2
  fcontoursigma3 = StoreCopy.fcontoursigma3; // contours sigma2

  /* Likelihood scans */
  fscanslikelihood=StoreCopy.fscanslikelihood;

  /* Data with excess and mean energy in bin */
  fBandArray = StoreCopy.fBandArray; // Bands containing data, bin emean, sth for residuals and plotfactory

  /* Data binned in time */
  fMapTimeBinVector = StoreCopy.fMapTimeBinVector;
  fLightCurveIntegratedFluxEnergyRange = StoreCopy.fLightCurveIntegratedFluxEnergyRange;

  /* Residuals */
  fresiduals = StoreCopy.fresiduals; // residuals
  fresidualssigmaplus = StoreCopy.fresidualssigmaplus; // residuals error + at 1 sigma
  fresidualssigmaminus = StoreCopy.fresidualssigmaminus; // residuals error - at 1 sigma
  fresiduals3sigmaplus = StoreCopy.fresiduals3sigmaplus; // residuals error + at 3 sigma
  fresiduals3sigmaminus = StoreCopy.fresiduals3sigmaminus; // residuals error - at 3 sigma
  fflux = StoreCopy.fflux; // experimental flux
  ffluxsigmaplus = StoreCopy.ffluxsigmaplus; // error + on experimental flux at 1 sigma
  ffluxsigmaminus = StoreCopy.ffluxsigmaminus; // error - on experimental flux at 1 sigma
  fflux3sigmaplus = StoreCopy.fflux3sigmaplus; // error + on experimental flux at 3 sigma
  fflux3sigmaminus = StoreCopy.fflux3sigmaminus; // error - on experimental flux at 3 sigma
  femean = StoreCopy.femean; // mean energy bin

  fMapResiduals = StoreCopy.fMapResiduals;

  /* Butterflies */
  for(std::map<std::string,TPolyLine*>::const_iterator butt=StoreCopy.fMapButterfly.begin(); 
      butt!=StoreCopy.fMapButterfly.end(); ++butt) {
    if(butt->second!=0) fMapButterfly[butt->first] = new TPolyLine(*(butt->second));
  }

  /* Integrator */
  
  frelativeprecision = StoreCopy.frelativeprecision;
  fabsoluteprecision = StoreCopy.fabsoluteprecision;
  fIntegrationType = StoreCopy.fIntegrationType;
   
  /* Spectral type */
  fSpectralType = StoreCopy.fSpectralType;

}

/**
 * \brief Function to make a copy of this class (necessary because polymorphism of Hypothesis)
 */
START::SumHypothesis* START::SumHypothesis::clone() const
{
  return new SumHypothesis(*this);
}

/**
 * \brief Return the flux for a given energy x used for Minuit2's callings
 */
double START::SumHypothesis::GetFlux(double x) const
{
  if(fiscomparison) {
    return (fparam[0]*fHypo1->GetFluxFitParams(x)+(1.-fparam[0])*fHypo2->GetFluxFitParams(x));
  }
  else {
    std::vector<double> paramhypo1;
    for(unsigned int ipar(0); ipar<fHypo1->GetParametersNb(); ipar++) paramhypo1.push_back(fparam[ipar]);
    fHypo1->SetParameters(paramhypo1);

    std::vector<double> paramhypo2;
    for(unsigned int ipar(fHypo1->GetParametersNb()); ipar<fHypo2->GetParametersNb(); ipar++) paramhypo2.push_back(fparam[ipar]);
    fHypo2->SetParameters(paramhypo2);

    return (fHypo1->GetFlux(x)+fHypo2->GetFlux(x));
  }
}


/**
 * \brief Define the first derivatives of the flux.
 * \param x energy
 */
void START::SumHypothesis::DerivativesFormulae(double x)
{

  ffirstderivatives.clear();

  if(!fiscomparison) {
    
    std::vector<double> paramfithypo1;
    for(unsigned int ipar(0); ipar<fHypo1->GetParametersNb(); ipar++) 
      paramfithypo1.push_back(fparamfit[ipar]);
    fHypo1->SetFittedParameters(paramfithypo1);
    fHypo1->DerivativesFormulae(x);

    std::vector<double> paramfithypo2;
    for(unsigned int ipar(fHypo1->GetParametersNb()); ipar<fHypo1->GetParametersNb()+fHypo2->GetParametersNb(); ipar++) 
      paramfithypo2.push_back(fparamfit[ipar]);
    fHypo2->SetFittedParameters(paramfithypo2);
    fHypo2->DerivativesFormulae(x);

    for(unsigned int ipar(0); ipar<fHypo1->GetParametersNb(); ipar++) {
      ffirstderivatives.push_back(fHypo1->GetFirstDerivatives()[ipar]);
    }
    for(unsigned int ipar(0); ipar<fHypo2->GetParametersNb(); ipar++) {
      ffirstderivatives.push_back(fHypo2->GetFirstDerivatives()[ipar]);
    }

  }
  

}

/**
 * \brief Initialization of parameter's name and units
 *
 * Parameters will be named with original parameters names plus the name of the hypothesis
 */
void START::SumHypothesis::InitParametersCaracteristics()
{

  if(fiscomparison) {
    fparamname[0] = "X";
    flatexparamname[0] = "#mchi";
    fparamunits[0] = ""; 
  }
  else {

    for(unsigned int ipar(0); ipar<fHypo1->GetParametersNb(); ipar++) {
      std::string name(fHypo1->GetParametersNames()[ipar]), latexname(fHypo1->GetLatexParametersNames()[ipar]), units(fHypo1->GetParametersUnits()[ipar]);
      name+="_";
      name+=fHypo1->GetName();
      fparamname.push_back(name);
      flatexparamname.push_back(latexname);
      fparamunits.push_back(units);
    }

    for(unsigned int ipar(0); ipar<fHypo2->GetParametersNb(); ipar++) {
      std::string name(fHypo2->GetParametersNames()[ipar]), latexname(fHypo2->GetLatexParametersNames()[ipar]), units(fHypo2->GetParametersUnits()[ipar]);
      name+="_";
      name+=fHypo2->GetName();
      fparamname.push_back(name);
      flatexparamname.push_back(latexname);
      fparamunits.push_back(units);
    }

  }

}
