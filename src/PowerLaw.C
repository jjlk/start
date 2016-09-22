// STL
#include <vector>

// ROOT
#include <TMath.h>

// START
#include "PowerLaw.hh"

// Utilities
#define DEBUG 0
#include "debugging.hh"

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::PowerLaw)
#endif

/**
 * \brief Default constructor needed by root
 */
START::PowerLaw::PowerLaw()
  :Hypothesis()
{
}

/**
 * \brief Constructor
 */
START::PowerLaw::PowerLaw(TString name, SpectralType type, double phi0,double gamma, double phi0err, double gammaerr)
  :Hypothesis(name)
{

  fSpectralType = type;

  fnumberofparameters=2;

  feref = 1.;

  fparam.resize(fnumberofparameters,0.);
  fparamerr.resize(fnumberofparameters,0.);
  fparamname.resize(fnumberofparameters,"");
  flatexparamname.resize(fnumberofparameters,"");
  fcovariancefit.resize(fnumberofparameters,std::vector<double >(fnumberofparameters,0.));
  fparamunits.resize(fnumberofparameters,"");
  ffirstderivatives.resize(fnumberofparameters,0.);

  fparam[0] = phi0;
  fparam[1] = gamma;

  fparamerr[0] = phi0err;
  fparamerr[1] = gammaerr;

  InitParametersCaracteristics();

  fmodelname = "PowerLaw";

  flatexformula = "#phi_{0} #left(#frac{E}{E_{Ref}} #right)^{-#Gamma}";

  SetLowerLimitOnNormalizationParameter(0,TMath::Power(10,-18));

  SetIntegratedFitEnergyRange(1.,10.); // user should change this for integrated flux

}

/**
 * \brief Destructor
 */
START::PowerLaw::~PowerLaw()
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
}

/**
 * \brief Copy Constructor
 */
START::PowerLaw::PowerLaw(PowerLaw const &StoreCopy)
//  :Hypothesis(StoreCopy.fName+TString("_copy"))
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
START::PowerLaw* START::PowerLaw::clone() const
{
  return new PowerLaw(*this);
}

/**
 * \brief Return the flux for a given energy x used for Minuit2's callings
 */
double START::PowerLaw::GetFlux(double x) const
{
  double a,b,c;
  
  switch(fSpectralType) {
  case Differential:
    return (fparam[0]*TMath::Power((x/feref),-fparam[1]));
    break;
  case Integrated:
    a = 1.-fparam[1];
    b = TMath::Power(fintegratedfitenergyrange.second,a);
    c = TMath::Power(fintegratedfitenergyrange.first,a);
    return (fparam[0]*a/(b-c)*TMath::Power(x,-fparam[1]));
    break;
  case EnergyFlux:
    a = 2.-fparam[1];
    b = TMath::Power(fintegratedfitenergyrange.second,a);
    c = TMath::Power(fintegratedfitenergyrange.first,a);
    return (fparam[0]*a/(b-c)*TMath::Power(x,-fparam[1]));
    break;
  default:
    std::cout << "Unknown Spectral Type !" << std::endl;
    return 0.;
    break;
  }
}

/**
 * \brief Return the flux for a given energy x with fitted parameter
 */
/*
double START::PowerLaw::GetFluxFitParams(double x) const
{
  fparam = fparamfit;
  return GetFlux(x);
}
*/
/**
 * \brief Return mean value for bin by taking into acount the hypothesis.
 * \param x1 lower edge bin energy 
 * \param x1 higher edge bin energy 
 */
double START::PowerLaw::GetMeanBinEnergy(double x1, double x2) const
{

  double gm1 =  fparamfit[1]-1;
  double gm2 =  fparamfit[1]-2;
  double diff1= TMath::Power(x2,-gm1) - TMath::Power(x1,-gm1);
  double diff2= TMath::Power(x2,-gm2) - TMath::Power(x1,-gm2);
  double energy_mean = (gm1/gm2)*(diff2/diff1);

  return energy_mean;

}

/**
 * \brief Define the first derivatives of the flux.
 * \param x energy
 */
void START::PowerLaw::DerivativesFormulae(double x)
{
  
  ffirstderivatives[0] = TMath::Power((x/feref),-fparamfit[1]); // dphi/dphi0
  ffirstderivatives[1] = -fparamfit[0]*TMath::Log(x/feref)*TMath::Power((x/feref),-fparamfit[1]); // dphi/dgamma
  
}

/**
 * \brief Initialization of parameter's name and units
 */
void START::PowerLaw::InitParametersCaracteristics() {

  switch(fSpectralType) {
  case Differential:
    fparamname[0] = "phi0";
    fparamname[1] = "Gamma";
    flatexparamname[0] = "#phi_{0}";
    flatexparamname[1] = "#Gamma";
    fparamunits[0] = "cm^{-2}.s^{-1}.TeV^{-1}"; 
    fparamunits[1] = ""; 
    break;
  case Integrated:
    fparamname[0] = "IF";
    fparamname[1] = "Gamma";
    flatexparamname[0] = "I_{F}";
    flatexparamname[1] = "#Gamma";
    fparamunits[0] = "cm^{-2}.s^{-1}"; 
    fparamunits[1] = ""; 
    break;
  case EnergyFlux:
    fparamname[0] = "EF";
    fparamname[1] = "Gamma";
    flatexparamname[0] = "E_{F}";
    flatexparamname[1] = "#Gamma";
    fparamunits[0] = "erg.cm^{-2}.s^{-1}"; 
    fparamunits[1] = ""; 
  }

}
