// STL
#include <vector>

// ROOT
#include <TMath.h>

// START
#include "SmoothBrokenPowerLaw.hh"

// Utilities
#define DEBUG 1
#include "debugging.hh"

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::SmoothBrokenPowerLaw)
#endif

/**
 * \brief Default constructor needed by root
 */
START::SmoothBrokenPowerLaw::SmoothBrokenPowerLaw()
:Hypothesis()
{
  fFunctorFlux=0; fFunctorFluxTimesEnergy=0;
  fFluxIntegrant=0; fFluxTimesEnergyIntegrant=0;
}

/**
 * \brief Constructor. Parameters are initial values for the fit
 */
START::SmoothBrokenPowerLaw::SmoothBrokenPowerLaw(TString name, SpectralType type, double phi0, double ebreak,
					   double beta, double gamma1, double gamma2, double phi0err, 
					   double ebreakerr, double betaerr, double gamma1err,
					   double gamma2err)
  :Hypothesis(name)
{
  
  fSpectralType = type;

  fnumberofparameters=5;

  feref = 1.;

  fparam.resize(fnumberofparameters,0.);
  fparamerr.resize(fnumberofparameters,0.);
  fparamname.resize(fnumberofparameters,"");
  flatexparamname.resize(fnumberofparameters,"");
  fcovariancefit.resize(fnumberofparameters,std::vector<double >(fnumberofparameters,0.));
  fparamunits.resize(fnumberofparameters,"");
  ffirstderivatives.resize(fnumberofparameters,0.);

  fparam[0] = phi0;
  fparam[1] = ebreak;
  fparam[2] = beta;
  fparam[3] = gamma1;
  fparam[4] = gamma2;

  fparamerr[0] = phi0err;
  fparamerr[1] = ebreakerr;
  fparamerr[2] = betaerr;
  fparamerr[3] = gamma1err; 
  fparamerr[4] = gamma2err; 

  InitParametersCaracteristics();

  fmodelname = "SmoothBrokenPowerLaw";

  flatexformula = 
    "#phi_{0} #left (#frac{E}{E_{Ref}} #right)^{-#Gamma_{1}} #left (1 + (#frac{E}{E_{break}} #right)^{-#frac{#Gamma_{2}-#Gamma_{1}}{#beta}} #right)^{-#beta}";

  SetLowerLimitOnNormalizationParameter(0,TMath::Power(10,-18));

  SetIntegratedFitEnergyRange(1.,10.); // user should change this to integrate flux

  // Integrator for integrated flux

  fFunctorFlux=0; fFunctorFluxTimesEnergy=0;
  fFluxIntegrant=0; fFluxTimesEnergyIntegrant=0;

  fFunctorFlux = new  ROOT::Math::Functor1D(this,&START::SmoothBrokenPowerLaw::FluxFormula);
  fFunctorFluxTimesEnergy = new  ROOT::Math::Functor1D(this,&START::SmoothBrokenPowerLaw::FluxTimesEnergy);

  fFluxIntegrant = new ROOT::Math::GSLIntegrator(fIntegrationType,fabsoluteprecision,frelativeprecision);
  fFluxIntegrant->SetFunction(*fFunctorFlux);

  fFluxTimesEnergyIntegrant = new ROOT::Math::GSLIntegrator(fIntegrationType,fabsoluteprecision,frelativeprecision);
  fFluxTimesEnergyIntegrant->SetFunction(*fFunctorFluxTimesEnergy);

}

/**
 * \brief Function to make a copy of this class (necessary because polymorphism of Hypothesis)
 *
 */
START::SmoothBrokenPowerLaw* START::SmoothBrokenPowerLaw::clone() const
{
  return new SmoothBrokenPowerLaw(*this);
}

/**
 * \brief Destructor
 */
START::SmoothBrokenPowerLaw::~SmoothBrokenPowerLaw()
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

  if(fFunctorFlux!=0) delete fFunctorFlux;
  fFunctorFlux=0;
  if(fFunctorFluxTimesEnergy!=0) delete fFunctorFluxTimesEnergy;
  fFunctorFluxTimesEnergy=0;
  if(fFluxIntegrant!=0) delete fFluxIntegrant;
  fFluxIntegrant=0;
  if(fFluxTimesEnergyIntegrant!=0) delete fFluxTimesEnergyIntegrant;
  fFluxTimesEnergyIntegrant=0;

  /* Butterflies */
  for(std::map<std::string,TPolyLine*>::iterator butt=fMapButterfly.begin(); butt!=fMapButterfly.end(); ++butt) {
    if(butt->second!=0) delete butt->second;
    butt->second = 0;
  }
  fMapButterfly.clear();

}

/**
 * \brief Copy constructor
 */
START::SmoothBrokenPowerLaw::SmoothBrokenPowerLaw(SmoothBrokenPowerLaw const &StoreCopy)
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

  /* Pointers*/
  fFunctorFlux=0;
  fFunctorFluxTimesEnergy=0;
  fFluxIntegrant=0;
  fFluxTimesEnergyIntegrant=0;

  fFunctorFlux=0; fFunctorFluxTimesEnergy=0;
  fFluxIntegrant=0; fFluxTimesEnergyIntegrant=0;

  fFunctorFlux = new  ROOT::Math::Functor1D(this,&START::SmoothBrokenPowerLaw::FluxFormula);
  fFunctorFluxTimesEnergy = new  ROOT::Math::Functor1D(this,&START::SmoothBrokenPowerLaw::FluxTimesEnergy);

  fFluxIntegrant = new ROOT::Math::GSLIntegrator(fIntegrationType,fabsoluteprecision,frelativeprecision);
  fFluxIntegrant->SetFunction(*fFunctorFlux);

  fFluxTimesEnergyIntegrant = new ROOT::Math::GSLIntegrator(fIntegrationType,fabsoluteprecision,frelativeprecision);
  fFluxTimesEnergyIntegrant->SetFunction(*fFunctorFluxTimesEnergy);

}

/**
 * brief Definition of the flux used in GetFlux divided by normalization.
 * We need this function in order to compute the integrated flux
 */
double START::SmoothBrokenPowerLaw::FluxFormula(double x) const {

  return(TMath::Power(x/feref,-fparam[3])*TMath::Power(1+TMath::Power(x/fparam[1],(fparam[4]-fparam[3])/fparam[2]),-fparam[2]));

}

/**
 * \brief Definition of the flux time the energy divides by normalization used in GetFlux.
 * We need this function in order to compute the energy flux
 */
double START::SmoothBrokenPowerLaw::FluxTimesEnergy(double x) const {
  return FluxFormula(x)*x;
}
	      
/**
 * \brief Return the flux for a given energy x.
 *
 * Used for Minuit2's callings.
 *
 * We use numerical integration for integrated flux
 */
double START::SmoothBrokenPowerLaw::GetFlux(double x) const
{
  double integrant(0);
  switch(fSpectralType) {
  case Differential:
    return (fparam[0]*FluxFormula(x));
    break;
  case Integrated:
    integrant = fFluxIntegrant->Integral(fintegratedfitenergyrange.first,fintegratedfitenergyrange.second);
    return (fparam[0]/integrant*FluxFormula(x));
    break;
  case EnergyFlux:
    integrant = fFluxTimesEnergyIntegrant->Integral(fintegratedfitenergyrange.first,fintegratedfitenergyrange.second);
    return (fparam[0]/integrant*FluxFormula(x));
    break;
  default:
    std::cout << "Unknown Spectral Type !" << std::endl;
    return 0.;
    break;
  }

}

/**
 * \brief Define the first derivatives of the flux.
 */
void START::SmoothBrokenPowerLaw::DerivativesFormulae(double x)
{
  double phi0(fparamfit[0]);
  double ebreak(fparamfit[1]);
  double beta(fparamfit[2]);
  double gamma1(fparamfit[3]);
  double gamma2(fparamfit[4]);
  double xebreak(x/ebreak);
  double xeref(x/feref);

  double pwl(TMath::Power(xeref,-gamma1));
  double index((gamma2-gamma1)/beta);
  double corr1(TMath::Power(xebreak,index));
  double corr2(1.+corr1);
  double globalcorr1(TMath::Power(corr1,-beta));
  double globalcorr2(TMath::Power(corr2,-beta));
  
  ffirstderivatives[0] = GetFluxFitParams(x)/fparamfit[0]; // dphi/dphi0
  ffirstderivatives[1] = phi0*index/beta*pwl*TMath::Power(xebreak,index-1.)*TMath::Power(corr2,-beta-1.); // dphi/dEbreak
  ffirstderivatives[2] = phi0*pwl*globalcorr2*(index*corr1*TMath::Log(xeref)/corr2-TMath::Log(corr2)); // dphi/dbeta
  ffirstderivatives[3] = phi0*pwl*(-globalcorr2*TMath::Log(xeref)+corr1*TMath::Power(corr2,-beta-1.)*TMath::Log(xebreak)); // dphi/dgamma1
  ffirstderivatives[4] = -phi0*pwl*corr1*TMath::Power(corr2,-beta-1.)*TMath::Log(xebreak);

}


/**
 * \brief Initialization of parameter's name and units
 */
void START::SmoothBrokenPowerLaw::InitParametersCaracteristics() {

  switch(fSpectralType) {
  case Differential:
    fparamname[0] = "phi0";
    fparamname[1] = "Ebreak";
    fparamname[2] = "beta";
    fparamname[3] = "Gamma1";
    fparamname[4] = "Gamma2";
    flatexparamname[0] = "#phi_{0}";
    flatexparamname[1] = "E_{break}";
    flatexparamname[2] = "#beta";
    flatexparamname[4] = "#Gamma1";
    flatexparamname[4] = "#Gamma2";
    fparamunits[0] = "cm^{-2}.s^{-1}.TeV^{-1}"; 
    fparamunits[1] = "TeV"; 
    fparamunits[2] = ""; 
    fparamunits[3] = ""; 
    fparamunits[4] = ""; 
    break;
  case Integrated:
    fparamname[0] = "IF";
    fparamname[1] = "Ebreak";
    fparamname[2] = "beta";
    fparamname[3] = "Gamma1";
    fparamname[4] = "Gamma2";
    flatexparamname[0] = "I_{F}";
    flatexparamname[1] = "E_{break}";
    flatexparamname[2] = "#beta";
    flatexparamname[3] = "#Gamma1";
    flatexparamname[4] = "#Gamma2";
    fparamunits[0] = "cm^{-2}.s^{-1}"; 
    fparamunits[1] = "TeV"; 
    fparamunits[2] = ""; 
    fparamunits[3] = ""; 
    fparamunits[4] = ""; 
    break;
  case EnergyFlux:
    fparamname[0] = "EF";
    fparamname[1] = "Ebreak";
    fparamname[2] = "beta";
    fparamname[3] = "Gamma1";
    fparamname[4] = "Gamma2";
    flatexparamname[0] = "E_{F}";
    flatexparamname[1] = "E_{break}";
    flatexparamname[2] = "#beta";
    flatexparamname[3] = "#Gamma1";
    flatexparamname[4] = "#Gamma2";
    fparamunits[0] = "erg.cm^{-2}.s^{-1}"; 
    fparamunits[1] = "TeV"; 
    fparamunits[2] = ""; 
    fparamunits[3] = "";  
    fparamunits[4] = "";  
  }

}
