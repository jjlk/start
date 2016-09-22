// STL
#include <iostream>
#include <iomanip>
#include <utility>
#include <string>
#include <sstream>

// ROOT
#include <TMath.h>
#include <Math/Derivator.h>

// START
#include "Hypothesis.hh"
//#include "GSLError.h"

#define DEBUG 0
#include "debugging.hh"

#define INFO std::cout << INFOCOLOR << "Hypothesis> " << RESETCOLOR
#define WARNING std::cout << WARNINGCOLOR << "Hypothesis> " << RESETCOLOR

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::Hypothesis)
#endif

/**
 * \brief Constructor
 *
 * Default needed by root (no instanciation inside this constructor, should not be used by the user)
 */
START::Hypothesis::Hypothesis()
:TNamed(),
  fIntegrantFlux(0),
  fIntegrantFluxTimesE(0),
  fFunctorIntegrantFlux(0),
  fFunctorIntegrantFluxTimesE(0),
  fIntegrantFitFlux(0),
  fIntegrantFitFluxTimesE(0),
  fFunctorIntegrantFitFlux(0),
  fFunctorIntegrantFitFluxTimesE(0)
{
  
}

/**
 * \brief Constructor
 *
 * Abstract class so it can't be instantiate!
 */
START::Hypothesis::Hypothesis(TString hypothesisname)
  :TNamed(hypothesisname,hypothesisname),
   fnumberofparameters(0),
   feref(1.),
   fmaximumlikelihood(0.),
   fedm(0.),
   fiteration(0),
   fisminimumvalid(false),
   fiscovariancevalid(false),
   farefittedparametersvalid(false),
   fconvergence(false),
   fIntegrantFlux(0),
   fIntegrantFluxTimesE(0),
   fFunctorIntegrantFlux(0),
   fFunctorIntegrantFluxTimesE(0),
   fIntegrantFitFlux(0),
   fIntegrantFitFluxTimesE(0),
   fFunctorIntegrantFitFlux(0),
   fFunctorIntegrantFitFluxTimesE(0)
{
  fparam.clear();
  fparamfit.clear();
  fparamerr.clear();
  fparamerrfit.clear();
  fcovariancefit.clear();
  fparamname.clear();
  fnormalizedparam.clear();
  fparamunits.clear();

  fcontoursigma1.clear();
  fcontoursigma2.clear();
  fcontoursigma3.clear();

  fBandArray.clear();

  fresiduals.clear();
  fresidualssigmaplus.clear();
  fresidualssigmaminus.clear();
  fresiduals3sigmaplus.clear();
  fresiduals3sigmaminus.clear();
  fflux.clear();
  ffluxsigmaplus.clear();
  ffluxsigmaminus.clear();
  fflux3sigmaplus.clear();
  fflux3sigmaminus.clear();
  femean.clear();

  ffixedparameter.clear();
  flimitedparameter.clear();
  fnormalizedparam.clear();

  ffirstderivatives.clear();

  fIntegrationType = ROOT::Math::IntegrationOneDim::kADAPTIVE; // type  of integration
  fabsoluteprecision = 0.; 
  frelativeprecision = 0.0001;

  ffittedintegratedflux = std::make_pair(-1.,-1.);
  ffittedenergyflux = std::make_pair(-1.,-1.);

  flineardecorrelationenergy = -1;
  flogarithmdecorrelationenergy = -1;
  fcontoursdecorrelationenergy = -1;

  /* Functor to put functions into Root objects*/

  fFunctorIntegrantFlux = new ROOT::Math::Functor1D(this,&START::Hypothesis::GetFlux);
  fFunctorIntegrantFluxTimesE = new ROOT::Math::Functor1D(this,&START::Hypothesis::GetFluxTimesEnergy);

  fFunctorIntegrantFitFlux = new ROOT::Math::Functor1D(this,&START::Hypothesis::GetFluxFitParams); // fit param
  fFunctorIntegrantFitFluxTimesE = new ROOT::Math::Functor1D(this,&START::Hypothesis::GetFluxFitParamsTimesEnergy);

  /* Integral */
  fIntegrantFlux = new ROOT::Math::GSLIntegrator(fIntegrationType,fabsoluteprecision,frelativeprecision); 
  fIntegrantFlux->SetFunction(*fFunctorIntegrantFlux);

  fIntegrantFluxTimesE = new ROOT::Math::GSLIntegrator(fIntegrationType,fabsoluteprecision,frelativeprecision); 
  fIntegrantFluxTimesE->SetFunction(*fFunctorIntegrantFluxTimesE);

  fIntegrantFitFlux = new ROOT::Math::GSLIntegrator(fIntegrationType,fabsoluteprecision,frelativeprecision); // fit params
  fIntegrantFitFlux->SetFunction(*fFunctorIntegrantFitFlux);

  fIntegrantFitFluxTimesE = new ROOT::Math::GSLIntegrator(fIntegrationType,fabsoluteprecision,frelativeprecision); // fit params
  fIntegrantFitFluxTimesE->SetFunction(*fFunctorIntegrantFitFluxTimesE);

  /* Root finder */
  fRootFinderType = ROOT::Math::RootFinder::kGSL_BRENT;

  /* Butterflies */
  for(std::map<std::string,TPolyLine*>::iterator butt=fMapButterfly.begin(); butt!=fMapButterfly.end(); ++butt) {
    butt->second = 0;
  }
  fMapButterfly.clear();

}

/*
 * \brief Destructor
 */
START::Hypothesis::~Hypothesis()
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
 * \brief Return the flux for a given energy x with fitted parameters.
 *
 * We use numerical integration for integrated flux
 */
double START::Hypothesis::GetFluxFitParams(double x) const
{
  fparam = fparamfit;
  return GetFlux(x);
}

/**
 * \brief Function used to put an Hypothesis derived based class in
 * an TF1 object. To do so, build your TF1 as follows : 
 * \code
 * Hypothesis *hypo = new PowerLaw("PWL");
 * hypo->SetParameters(hypo->GetFittedParameters()); // or a vector
 * TF1 *f = new TF1("TF1PowerLaw",hypo,emin,emax,0,"Hypothesis");
 * \endcode
 *
 */
double START::Hypothesis::operator() (double *x, double *p)
{
  double energy = x[0];
  if(p!=0) std::cout << "strange..." << std::endl;
  return GetFlux(energy);
}

/**
 * \brief Return variance of flux for energy x
 */
double START::Hypothesis::GetSigmaFlux(double x)
{

  double sigma(0.);

  DerivativesFormulae(x); // derivatives at energy x

  // variance of the flux = sum_i dphi/dxi|fit*cov(xi,xi) + 2*sum_i sum_j dphi/dxi|fit dphi/dxj|fit*cov(xi,xj)

  for(unsigned int ipar(0); ipar<fnumberofparameters; ipar++) {

    double tmplin = ffirstderivatives[ipar]*ffirstderivatives[ipar]*fcovariancefit[ipar][ipar];
    sigma+=tmplin;

  }

  for(int unsigned ipar1(0); ipar1<fnumberofparameters; ipar1++) {

    for(int unsigned ipar2(ipar1+1); ipar2<fnumberofparameters; ipar2++) {

      double tmpcross = 2*ffirstderivatives[ipar1]*ffirstderivatives[ipar2]*fcovariancefit[ipar1][ipar2];
      sigma+=tmpcross;

    }

  }

  return TMath::Sqrt(sigma);

}

/**
 * \brief Return the value of \f$ log10(\Phi(E)+\sigma(E)) - log10(\Phi(E)+\sigma(E)) \f$ at energy \f$ E \f$
 */
double START::Hypothesis::GetLinearWidthSigma(double x)
{

  double sigma = GetSigmaFlux(x);
  double flux = GetFluxFitParams(x);
  double value = TMath::Log10(flux+sigma)-TMath::Log10(flux-sigma);

  if(TMath::IsNaN(value)) {
    DEBUG_OUT << "Quantity log10(flux+sigma) - log10(flux-sigma) = " << value << " is undefined!!" << std::endl;
    DEBUG_OUT << "Flux = "<< GetFluxFitParams(x) << " Sigma = " << GetSigmaFlux(x) << std::endl;
    return TMath::Log10(flux+sigma); // if we return a constant, the derivative will be zero for 
                                     // a while and the decorrelation energy will be misdefined.
                                     // We might have the good value with some luck with this primitive technic...
  }

  return (value);

}

/**
 * \brief Return the value of \f$ \Phi(E) exp(+\frac{\sigma(E)}{\Phi(E)}) - \Phi(E) exp(-\frac{\sigma(E)}{\Phi(E)}) \f$
 */
double START::Hypothesis::GetLogarithmWidthSigma(double x)
{

  double flux = GetFluxFitParams(x);
  double sigma = GetSigmaFlux(x);
  double expo = TMath::Exp(sigma/flux);
  double value = TMath::Log10(flux*expo)-TMath::Log10(flux/expo);

  if(TMath::IsNaN(value)) {
    // JLK : shouldn't happened
    WARNING << "Quantity Phi(E)*exp(+sigma(E)/Phi(E)) - Phi(E)*exp(-sigma/Phi(E)) = " << value << " is undefined!!" << std::endl;
    INFO << "Flux = "<< flux << " Sigma = " << sigma << std::endl;
    return 1.;
  }

  return (value);

}

/**
 * \brief Return the value of width of the butterfly defined by the contours at a fixed energy
 * If the contours are not computed, we return START::Hypothesis::GetLinearWidthSigma
 *
 * \warning this function only works for PowerLaw hypothesis since we don't have multidimensional contours
 * (for the moment)
 */
double START::Hypothesis::GetContoursWidthSigma(double x)
{

  if(fnumberofparameters>2) {
    INFO << "Hypothesis "<< GetName() << "not suited to have contour's butterfly so you will have linear decorrelation energ" << std::endl;
    return GetLinearWidthSigma(x);
  }

  // 1 sigma
  if(fcontoursigma1.size()==0) {
    INFO << "No contours, so you will have linear decorrelation energy" << std::endl;
    return GetLinearWidthSigma(x);
  }

  // JLK : Un peu hardcode pour le moment sachant que l'on ne sait pas faire les contours a D>2...
  // donc en gros ca ne marchera que pour la pwl... sooréé

  unsigned int index(999);
  bool paramreversed=false;

  // we look for phi0 and gamma for now...
  for(unsigned int icont(0); icont<fcontoursigma1.size(); icont++) {

    if(fparamname[fcontoursigma1[icont].second.first]=="phi0" && fparamname[fcontoursigma1[icont].second.second]=="Gamma") {
      index=icont;
    } 
    else if(fparamname[fcontoursigma1[icont].second.first]=="Gamma" && fparamname[fcontoursigma1[icont].second.second]=="phi0") {
      index=icont;
      paramreversed=true;
    }
  }

  if(index==999) {
    INFO << "Don't recognize parameters so you will have linear decorrelation energy" << std::endl;
    return GetLinearWidthSigma(x);
  }

  std::vector<std::pair<double,double> > pointscontour = fcontoursigma1[index].first; // get contours points 

  if(pointscontour.size()==0) {
    INFO << "Contours are empty so you will have linear decorelation energy" << std::endl;
    return GetLinearWidthSigma(x);
  }

  unsigned int npoints = pointscontour.size();
  std::vector<double> flux;

  for(unsigned int ipoint(0); ipoint<npoints; ipoint++) {
    
    std::vector<double> paramsforhypo;
    if(!paramreversed) {
      paramsforhypo.push_back(pointscontour[ipoint].first);
      paramsforhypo.push_back(pointscontour[ipoint].second);
    }
    else {
      paramsforhypo.push_back(pointscontour[ipoint].second);
      paramsforhypo.push_back(pointscontour[ipoint].first);
    }
    SetParameters(paramsforhypo);
    flux.push_back(GetFlux(x));
    
  }
  double valuemin(0),valuemax(0);
  std::vector<double>::const_iterator itmax, itmin;
  itmax = max_element(flux.begin(), flux.end());
  valuemax = *itmax;
  itmin = min_element(flux.begin(), flux.end());
  valuemin = *itmin;  

  return (TMath::Log10(valuemax)-TMath::Log10(valuemin));

}

/**
 * \brief Return the value of the first derivative of START::Hypothesis::GetLinearWidthSigma
 *
 * \param x energy
 */
double START::Hypothesis::GetFirstDerivativeLinearWidthSigma(double x)
{
  double dsigma(0.);

  ROOT::Math::Functor1D FunctorSigmaFlux(this,&START::Hypothesis::GetLinearWidthSigma);

  ROOT::Math::Derivator SigmaFluxDerivative(FunctorSigmaFlux);

  dsigma = SigmaFluxDerivative.Eval(x);

  return dsigma;

}

/**
 * \brief Return the value of the first derivative of START::Hypothesis::GetLogarithmWidthSigma
 *
 * \param x energy
 */
double START::Hypothesis::GetFirstDerivativeLogarithmWidthSigma(double x)
{
  double dsigma(0.);

  ROOT::Math::Functor1D FunctorSigmaFlux(this,&START::Hypothesis::GetLogarithmWidthSigma);

  ROOT::Math::Derivator SigmaFluxDerivative(FunctorSigmaFlux);

  dsigma = SigmaFluxDerivative.Eval(x);

  return dsigma;

}

/**
 * \brief Return the value of the first derivative of START::Hypothesis::GetContoursWidthSigma
 *
 * \param x energy
 */
double START::Hypothesis::GetFirstDerivativeContoursWidthSigma(double x)
{
  double dsigma(0.);

  ROOT::Math::Functor1D FunctorSigmaFlux(this,&START::Hypothesis::GetContoursWidthSigma);

  ROOT::Math::Derivator SigmaFluxDerivative(FunctorSigmaFlux);

  dsigma = SigmaFluxDerivative.Eval(x);

  return dsigma;

}

/**
 * \brief Return decorrelation energy from flux covariance
 */
double START::Hypothesis::GetLinearDecorrelationEnergy() 
{

  if(flineardecorrelationenergy==-1) 
    flineardecorrelationenergy = FindLinearDecorrelationEnergy();

  return flineardecorrelationenergy;
}

/**
 * \brief Return decorrelation energy from the log of the flux covariance
 */
double START::Hypothesis::GetLogarithmDecorrelationEnergy()
{

  if(flogarithmdecorrelationenergy==-1)
    flogarithmdecorrelationenergy = FindLogarithmDecorrelationEnergy();

  return flogarithmdecorrelationenergy;

}

/**
 * \brief Return decorrelation energy from contours
 */
double START::Hypothesis::GetContoursDecorrelationEnergy()
{

  if(fcontoursdecorrelationenergy==-1) {
    if(GetModelName()=="PowerLaw" && fcontoursigma1.size()>0)
      fcontoursdecorrelationenergy = FindContoursDecorrelationEnergy();
    else {
      WARNING << "Contours are note computed or are in bad shape or you are no using a PowerLaw hypothesis!" << std::endl;
      INFO << "So you will have decorrelation energy with flux covariance at first order" << std::endl;
    }
  }
  return fcontoursdecorrelationenergy;

}

/**
 * \brief Find decorrelation energy with log flux variance
 * A test is done to see if quantity flux - sigma becomes negative
 * an if it is negative, we take the one computed with log of the flux 
 * of the covariance
 */
double START::Hypothesis::FindLinearDecorrelationEnergy() 
{

  double emin(0.), emax(0.); // range of energy for the root finder
  bool findNaN(false); // if the value flux - sigma is negative it will be true
  double npoint(1000.);
  double estart(TMath::Log10(fminimizationenergyrange.first)), estop(TMath::Log10(fminimizationenergyrange.second));
  double ebinsize((estop-estart)/npoint);

  DEBUG_OUT << " minimization range : " << fminimizationenergyrange.first << " " << fminimizationenergyrange.second << std::endl;

  double testvalue(0); // value for emax if a negative value for flux - sigma is found

  for(double ien(estart); ien<estop*1.001; ien+=ebinsize) {
    double sigma = GetSigmaFlux(TMath::Power(10.,ien));
    double flux = GetFluxFitParams(TMath::Power(10.,ien));
    DEBUG_OUT << "energy = " << TMath::Power(10.,ien) << " flux = " << flux << " sigma = " << sigma << std::endl;
    if(flux-sigma<0) {
      DEBUG_OUT << "find it (top)!!" << std::endl;
      testvalue = TMath::Power(10.,ien)-TMath::Power(10.,ebinsize);
      findNaN=true;
      break;
    }
  }

  if(!findNaN) { // everything ok
    emin = fminimizationenergyrange.first;
    emax = fminimizationenergyrange.second;
  }
  else { // covariance is greater than flux
    emin = fminimizationenergyrange.first;
    emax = testvalue;
    WARNING << "Be carefull with inconsistencies!" << std::endl;
    INFO << "Negative value detected for flux-sigma wich will lead to NaN values," << " so decorrelation energy search will be done between " 
	 << emin << " and " << emax << " TeV." << std::endl;
  }

  DEBUG_OUT << " root finder range : " << emin << " " << emax << std::endl;

  ROOT::Math::Functor1D FunctorDerivativeSigmaFlux(this,&START::Hypothesis::GetFirstDerivativeLinearWidthSigma);

  ROOT::Math::RootFinder DecorrelationEnergyFinder(fRootFinderType);
  DecorrelationEnergyFinder.SetFunction(FunctorDerivativeSigmaFlux,emin,emax); 
  DecorrelationEnergyFinder.Solve();

  double edec(DecorrelationEnergyFinder.Root());

  // JLK if these conditions are met, decorrelation energy is probably wrong so we compute it in log
  if(edec==0. || edec==emin || edec==emax || edec<0.) {
    WARNING << "Decorrelation energy is probably wrong : E_d = " << edec << " TeV." << std::endl;
    INFO << "We compute decorrelation energy with START::Hypothesis::FindLogarithmDecorrelationEnergy..." << std::endl;
    edec = FindLogarithmDecorrelationEnergy();
    INFO << "New decorralation energy is : " << edec << " TeV." << std::endl;
  }

  return edec;

}

/**
 * \brief Find decorrelation energy with log flux variance
 */
double START::Hypothesis::FindLogarithmDecorrelationEnergy() 
{

  ROOT::Math::Functor1D FunctorDerivativeSigmaFlux(this,&START::Hypothesis::GetFirstDerivativeLogarithmWidthSigma);

  ROOT::Math::RootFinder DecorrelationEnergyFinder(fRootFinderType);
  DecorrelationEnergyFinder.SetFunction(FunctorDerivativeSigmaFlux,fminimizationenergyrange.first,fminimizationenergyrange.second); 
  DecorrelationEnergyFinder.Solve();

  return DecorrelationEnergyFinder.Root();

}

/**
 * \brief Find decorrelation energy with contours
 */
double START::Hypothesis::FindContoursDecorrelationEnergy() 
{

  ROOT::Math::Functor1D FunctorDerivativeSigmaFlux(this,&START::Hypothesis::GetFirstDerivativeContoursWidthSigma);

  ROOT::Math::RootFinder DecorrelationEnergyFinder(fRootFinderType);
  DecorrelationEnergyFinder.SetFunction(FunctorDerivativeSigmaFlux,fminimizationenergyrange.first,fminimizationenergyrange.second); 
  DecorrelationEnergyFinder.Solve();

  return DecorrelationEnergyFinder.Root();

}

/**
 * \brief Print parameters of the flux
 */
void START::Hypothesis::PrintParameters(std::ostream &os) const
{

  for(unsigned int i(0); i<fnumberofparameters; i++) {
    os << fparamname[i] << " = ("  << fparam[i] << " +/- " << fparamerr[i] 
	      << ")" << fparamunits[i]  << std::endl;
  }

}

/**
 * \brief Print fitted parameters of the flux
 */
void START::Hypothesis::PrintFittedParameters(std::ostream &os) const
{

  for(unsigned int i(0); i<fnumberofparameters; i++) {
    os << fparamname[i] << " = ("  << fparamfit[i] << " +/- " << fparamerrfit[i] 
	      << ")" << fparamunits[i]  << std::endl;
  }

}

/**
 * \brief Print covariance matrix
 */
void START::Hypothesis::PrintFittedCovarianceMatrix(std::ostream &os) const
{

  for(unsigned int i(0); i<fnumberofparameters; i++) {
    for(unsigned int j(0); j<fnumberofparameters; j++) {
      os << fcovariancefit[i][j] <<  "  ";
    }
    os << std::endl;
  }

}

/**
 * \brief Print minimization caracteristics such as edm, iterations' number and minimum of the likelihood
 */
void START::Hypothesis::PrintMinimizationCarateristics(std::ostream &os) const
{
  if(GetConvergence())
    os << "Minimization worked!" << std::endl;
  else os << "Minimization didn't work!" << std::endl;
  os << "Energy range [" << fminimizationenergyrange.first << ";" 
     << fminimizationenergyrange.second << "] TeV"<< std::endl;
  os << "Inteagrated flux energy range [" << fintegratedfitenergyrange.first << ";" 
     << fintegratedfitenergyrange.second << "] TeV"<< std::endl;
  os << "ERef=" << feref << " TeV" << std::endl;
  os << "Convergence = " << GetConvergence() << std::endl;
  os << "IsMinimumValid = " << GetIsMinimumValid() << std::endl;
  os << "IsCovarianceValid = " << GetIsCovarianceValid() << std::endl;
  os << "AreFittedParametersValid = " << GetAreFittedParametersValid() << std::endl;
  os << "maximum likelihood = " << std::setprecision(12) << fmaximumlikelihood << std::endl;
  os << "EDM = " << fedm << " (required EDM = " << 0.001*0.5*0.1 << ")" <<  std::endl;
  os << "number of iterations = " << fiteration << std::endl;
}

/**
 * \brief Print Bands summary
 */
void START::Hypothesis::PrintSummaryBand(std::ostream &os, bool UseEThreshold) const
{
  if(fBandArray.size()==0) return;
  Band SummaryBand(fBandArray[0]);
  SummaryBand.ClearBandInfo();
  SummaryBand.AddInfoFromBands(fBandArray,UseEThreshold);
  SummaryBand.PrintBand(os);
}

/**
 * \brief Print residuals, flux and errors at 1 and 3 sigma
 */
void START::Hypothesis::PrintResiduals(std::ostream &os) const
{
  for(std::map<std::string,Residuals>::const_iterator MapRes=fMapResiduals.begin(); 
      MapRes!=fMapResiduals.end(); ++MapRes) {
    os << "###Residuals "<< MapRes->first << ":"<< std::endl;
    MapRes->second.PrintResiduals(os);
    os << std::endl;
  }
}

/**
 * \brief Print butetrflies
 */
void START::Hypothesis::PrintButterflies(std::ostream &os) const // print residuals
{
  for(std::map<std::string,TPolyLine*>::const_iterator butt=fMapButterfly.begin(); 
      butt!=fMapButterfly.end(); ++butt) {
    os << "###" << butt->first << ":" << std::endl;
    for(int ipoint(0); ipoint<butt->second->GetN(); ipoint++) 
      os << butt->second->GetX()[ipoint] << " " << butt->second->GetY()[ipoint] << std::endl;
    os << std::endl;
  }
}

/**
 *
 */
void START::Hypothesis::PrintContours(std::ostream &os) const
{

  if(fcontoursigma1.size()==0) {
    os << "###No 1 sigma contours" << std::endl;
  }
  else {
    for(std::vector<std::pair<std::vector<std::pair<double,double> >, std::pair<int,int> > >::const_iterator cont=fcontoursigma1.begin(); 
	cont!=fcontoursigma1.end(); cont++) {
      std::vector<std::pair<double,double> > points = cont->first;
      std::pair<int,int> param = cont->second;
      os << "####1 sigma contours for " << fparamname[param.first] << "/" << fparamname[param.second] << ":" << std::endl;
      for(std::vector<std::pair<double,double> >::const_iterator xy=points.begin(); xy!=points.end(); xy++)
	os << xy->first << " " << xy->second << std::endl;
    }
  }

  if(fcontoursigma2.size()==0) {
    os << "###No 2 sigma contours" << std::endl;
  }
  else {
    for(std::vector<std::pair<std::vector<std::pair<double,double> >, std::pair<int,int> > >::const_iterator cont=fcontoursigma2.begin(); 
	cont!=fcontoursigma2.end(); cont++) {
      std::vector<std::pair<double,double> > points = cont->first;
      std::pair<int,int> param = cont->second;
      os << "####2 sigma contours for " << fparamname[param.first] << "/" << fparamname[param.second] << ":" << std::endl;
      for(std::vector<std::pair<double,double> >::const_iterator xy=points.begin(); xy!=points.end(); xy++)
	os << xy->first << " " << xy->second << std::endl;
    }
  }

  if(fcontoursigma3.size()==0) {
    os << "###No 3 sigma contours" << std::endl;
  }
  else {
    for(std::vector<std::pair<std::vector<std::pair<double,double> >, std::pair<int,int> > >::const_iterator cont=fcontoursigma3.begin(); 
	cont!=fcontoursigma3.end(); cont++) {
      std::vector<std::pair<double,double> > points = cont->first;
      std::pair<int,int> param = cont->second;
      os << "####3 sigma contours for " << fparamname[param.first] << "/" << fparamname[param.second] << ":" << std::endl;
      for(std::vector<std::pair<double,double> >::const_iterator xy=points.begin(); xy!=points.end(); xy++)
	os << xy->first << " " << xy->second << std::endl;
    }
  }

}

/**
 * \brief Print hypothesis' parameters
 */
void START::Hypothesis::Print(Option_t *option) const
{
  PrintHypothesis(std::cout);
}

/**
 * \brief Print hypothesis' parameters
 * \param os file or screen
 */
void START::Hypothesis::PrintHypothesis(std::ostream &os) const
{
  os << "Hypothesis " << GetName()
     << " for model " << fmodelname << std::endl; 
  os << std::endl;
  os << "###Minimization caracteristics : " << std::endl;
  PrintMinimizationCarateristics(os);
  os << std::endl;
  os << "###Fitted parameters : " << std::endl;
  PrintFittedParameters(os);
  os << std::endl;
  os << "###Fitted Covariance : " << std::endl;
  PrintFittedCovarianceMatrix(os);
  os << std::endl;
  os << "###Decorelation energy (JLK Find a way to determine which on it is):" << std::endl;
  os << "E_dec (Linear)=" << flineardecorrelationenergy << std::endl;
  os << "E_dec (Logarithm)=" << flogarithmdecorrelationenergy << std::endl;
  os << "E_dec (Contours)=" << fcontoursdecorrelationenergy << std::endl;
  os << std::endl;
  os << "###Integrated flux:" << ffittedintegratedflux.first << "+/-" << ffittedintegratedflux.second << std::endl;
  os << "###Energy flux:" << ffittedenergyflux.first << "+/-" << ffittedenergyflux.second << std::endl;
  
  if(GetConvergence()) {
    os << "###SummaryBand (All infos) : " << std::endl;
    PrintSummaryBand(os,false);
    os << std::endl;
    os << "###SummaryBand (above treshold) : " << std::endl;
    PrintSummaryBand(os,true);
    os << std::endl;
    PrintResiduals(os);
    os << std::endl;
    PrintContours(os);
    os << std::endl;
    for(std::map<std::string,TimeBinVector>::const_iterator MapTBin=fMapTimeBinVector.begin(); 
	MapTBin!=fMapTimeBinVector.end(); ++MapTBin) {
      os << "###Data binned in time "<< MapTBin->first << " [" << fLightCurveIntegratedFluxEnergyRange.first 
	 <<";" << fLightCurveIntegratedFluxEnergyRange.second << "] (TeV):"<< std::endl;
      MapTBin->second.PrintTimeBinVector(os);
      os << std::endl;
    }
    PrintButterflies(os);
  }
}

/**
 * \brief Set a lower limit on the parameter 'param' with value 'value'
 *
 * We normalize the components of the flux which have a lower limit (necessary)
 * All parameters such as phi0 in flux=phi0*E^-gamma must have a lowerlimit because :
 * <ul> 
 *    <li> the likelihood can't have a negative expected excess (it's our choice) </li>
 *    <li> we normalize them in a way which make it easy for MINUIT2 (precision's matters) </li>
 * </ul>
 *  The non-normalized/havenolimit parameters can't have limit(s) because MINUIT2 don't like it
 *  (internals changes of variables in MINUIT2) and it's again our choice to forbid it ant at least,
 *  to minimize those effects.
 *
 * \param param : parameter to normalize
 * \param value : lower limit of the parameter
 *
 * \todo add fixed parameters
 */
void START::Hypothesis::SetLowerLimitOnNormalizationParameter(unsigned int param, double value)
{
  fnormalizedparam.push_back(std::make_pair(param,value));
}

/**
 * \brief Return mean energy bin by taking account the hypothesis
 *
 * \param x1 min bin energy
 * \param x2 max bin energy
 */ 
double START::Hypothesis::GetMeanBinEnergy(double x1,double x2) const
{
  double meanenergy(0.);
  double a(0.), b(0.);

  a=fIntegrantFitFlux->Integral(x1,x2);
  b=fIntegrantFitFluxTimesE->Integral(x1,x2);

  meanenergy = b/a;

  return meanenergy;
}

/**
 * \brief Get integral of flux 
 *
 * \param x1 energy low
 * \param x1 energy high
 */
double START::Hypothesis::GetFluxIntegral(double x1, double x2) const 
{
  if(fIntegrantFlux!=0) return fIntegrantFlux->Integral(x1,x2);
  else return 1.;
}

/**
 * \brief Get integral of energy times flux 
 *
 * \param x1 energy low
 * \param x1 energy high
 */
double START::Hypothesis::GetFluxTimesEnergyIntegral(double x1, double x2) const 
{
  if(fIntegrantFluxTimesE!=0) return fIntegrantFluxTimesE->Integral(x1,x2);
  else return 1.;
}

/**
 * \brief Get integral of flux for fitted parameters
 *
 * \param x1 energy low
 * \param x1 energy high
 */
double START::Hypothesis::GetFluxIntegralFitParams(double x1, double x2) const 
{
  if(fIntegrantFitFlux!=0) return fIntegrantFitFlux->Integral(x1,x2);
  else return 1.;
}

/**
 * \brief Get integral of energy times flux  for fitted parameters
 *
 * \param x1 energy low
 * \param x1 energy high
 */
double START::Hypothesis::GetFluxTimesEnergyIntegralFitParams(double x1, double x2) const 
{
  if(fIntegrantFitFluxTimesE!=0) return fIntegrantFitFluxTimesE->Integral(x1,x2);
  else return 1.;
}

/**
 * \brief return flux times energy for fitted parameters
 * \param x energy
 */
double START::Hypothesis::GetFluxFitParamsTimesEnergy(double x) const
{
  return (GetFluxFitParams(x)*x);
}

/**
 * \brief return flux times energy
 * \param x energy
 */
double START::Hypothesis::GetFluxTimesEnergy(double x) const
{
  return (GetFlux(x)*x);
}

/**
 * \brief Fix a parameter for the minimization.
 * \param param Index of the parameter to fix
 * \param value Value of the parameter
 * \param redominimization if true, the minimization is done an other time with the release 
 * of the parameter
 */
void START::Hypothesis::FixParameter(unsigned int param, double value, bool redominimization) {
  if(param>=fnumberofparameters) {
    WARNING << "You ask to fix parameter " 
	    << param << ", which doesn't exist for hypothesis " << GetName() << "." << std::endl;
    INFO << "Minization will proceed without it..." << std::endl;
  }
  else
    ffixedparameter.push_back(std::make_pair(param,std::make_pair(value,redominimization)));
  DEBUG_OUT << "param fixed = " << param << std::endl;
  
  if(ffixedparameter.size()==fnumberofparameters) {
    WARNING << "You fixed all the parameter of the hypothesis " << GetName() << " !" << std::endl;
    INFO << "I don't know how MINUIT2 can deal with that "
	 << "so your fixed parameters we'll be erased..." << std::endl;
    //ffixedparameter.clear();
  }
  

}

/**
 * \brief Fix a parameter for the minimization.
 * \param param Name of the parameter to fix
 * \param value Value of the parameter
 * \param redominimization If true, the minimization is done an other time with the release 
 * of the parameter
 */
void START::Hypothesis::FixParameter(std::string param, double value, bool redominimization) {
  unsigned int indice(999999);
  for(unsigned int ipar(0); ipar<fnumberofparameters; ipar++) {
    if(param==fparamname[ipar]) {
      indice=ipar;
      FixParameter(indice,value,redominimization);
      break;
    }
  }
  
  if(indice==999999) {
    WARNING << "You asked to fix the parameter " 
	     << param << ", which doesn't exist for hypothesis " << GetName() << "." << std::endl;
    INFO << "Minization will proceed without it..." << std::endl;
  }

}


/**
 * \brief Limit a parameter for the minimization.
 * \param param Index of the parameter to fix
 * \param value Value of the parameter
 * \param redominimization if true, the minimization is done an other time with the release 
 * of the parameter
 */
void START::Hypothesis::LimitParameter(unsigned int param, double valuemin, double valuemax) {
  if(param>=fnumberofparameters) {
    WARNING << "You ask to limit parameter " 
	    << param << ", which doesn't exist for hypothesis " << GetName() << "." << std::endl;
    INFO << "Minization will proceed without it..." << std::endl;
  }
  else
    flimitedparameter.push_back(std::make_pair(param,std::make_pair(valuemin,valuemax)));
  DEBUG_OUT << "param limited = " << param << std::endl;

}

/**
 * \brief Return energy flux units
 */
std::string START::Hypothesis::GetEnergyFluxUnits() const {

  std::string units("");

  switch(fSpectralType) {
  case Differential:
  case Integrated:
  case EnergyFlux:
    units = "erg.cm^{-2}.s^{-1}";
  }

  return units;

}

/**
 * \brief Return integrated flux units
 */
std::string START::Hypothesis::GetIntegratedFluxUnits() const {

  std::string units("");

  switch(fSpectralType) {
  case Differential:
  case Integrated:
  case EnergyFlux:
    units = "cm^{-2}.s^{-1}";
  }

  return units;

}

/**
 * \brief Return flux units
 */
std::string START::Hypothesis::GetFluxUnits() const {

  std::string units("");

  switch(fSpectralType) {
  case Differential:
  case Integrated:
  case EnergyFlux:
    units = "cm^{-2}.s^{-1}.TeV^{-1}";
  }

  return units;

}

/**
 * \brief Return key for residuals
 * \param sigma number of sigma's residuals. Zero by default (unrebinned)
 */
std::string START::Hypothesis::MakeKeyResiduals(double sigma) {
  std::ostringstream key;
  if(sigma==0.) {
    key << "UnRebinned";
  }
  else {
    key.precision(2);
    key << sigma << "Sigmas_Rebinned";
  }
  return key.str();
}

/**
 * \brief Return key for residuals
 * \param TypeOfButterfly ButterflyType value
 */
std::string START::Hypothesis::MakeKeyButterfly(ButterflyType TypeOfButterfly) {
  std::ostringstream key;
  key << "Butterfly_";
  switch(TypeOfButterfly) {
  case LinearButterfly:
    key << "Linear";
    break;
  case LogarithmButterfly:
    key << "Logarithm";
    break;
  case ContoursButterfly:
    key << "Contours";
  }
  return key.str();
}
