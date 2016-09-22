#ifndef _LOGPARABOLIC_
#define _LOGPARABOLIC_

#include "Hypothesis.hh"

namespace START {

  /**
   * \brief LogParabolic Hypothesis with 3 parameters \f$\phi_0\f$, \f$\alpha\f$ and \f$\beta\f$ :
   * \f[
   * \Phi(E) = \phi_{0} \left (\frac{E}{E_{Ref}} \right)^{-\alpha -\beta \log \frac{E}{E_{Ref}}}
   * \f]
   */
  class LogParabolic : public Hypothesis
  {

  public:

    LogParabolic(); // default constructor

    LogParabolic(TString name, SpectralType type=Differential, double phi0=1.e-12,double alpha=2.5,
		 double beta=0., double phi0err=0., double alphaerr=0., double betaerr=0.); //values for parameters and errors

    LogParabolic(LogParabolic const &StoreCopy); // Copy constructor

    virtual ~LogParabolic(); // destructor

    virtual LogParabolic* clone() const; // VIM : Tricks to create a copy constructor
  
    virtual double GetFlux(double x) const; // return the value of the flux for current parameters

    virtual void DerivativesFormulae(double x); // Used for the butterfly

    virtual void InitParametersCaracteristics(); // Initialization of parameters names and units

  private:

    //Definition of the flux time the energy used in GetFlux.
    //We need this function in order to compute the integrated flux
    double FluxFormula(double x) const;
  
    //Definition of the flux time the energy used in GetFlux.
    //We need this function in order to compute the energy flux
    double FluxTimesEnergy(double x) const;

    // members needed to compute the integrals
    ROOT::Math::Functor1D *fFunctorFlux, *fFunctorFluxTimesEnergy; //!
    ROOT::Math::GSLIntegrator *fFluxIntegrant, *fFluxTimesEnergyIntegrant; //!

#ifdef CONSTRUCT_LIBRARY
  public:
    ClassDef(START::LogParabolic,1);
#endif
  };
}
#endif
