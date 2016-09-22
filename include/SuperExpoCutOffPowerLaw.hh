#ifndef _SUPEREXPOCUTOFFPOWERLAW_
#define _SUPEREXPOCUTOFFPOWERLAW_

#include "Hypothesis.hh"
namespace START {
  /**
   * \brief  PowerLaw with a super exponential cut-off Hypothesis with 4 parameters \f$\phi_0\f$, 
   * \f$\Gamma_{1}\f$, \f$\Gamma_{2}\f$ and \f$\beta\f$ :
   * \f[
   * \Phi(E) = \phi_0 \left (\frac{E}{E_{Ref}} \right)^{-\Gamma_{1}} \exp{[-(\beta E)^{\Gamma_{2}}]}
   * \f]
   */
  class SuperExpoCutOffPowerLaw : public Hypothesis
  {

  public:

    SuperExpoCutOffPowerLaw(); // default constructor

    SuperExpoCutOffPowerLaw(TString name, SpectralType type=Differential, double phi0=1.e-12,
			    double gamma1=2.5, double gamma2=0.1, double beta=0.1, // constructor used for minimization with initial 
			    double phi0err=0., double gamma1err=0., double gamma2err=0.,  // values for parameters and errors
			    double betaerr=0.);

    SuperExpoCutOffPowerLaw(SuperExpoCutOffPowerLaw const &StoreCopy); // Copy constructor

    virtual ~SuperExpoCutOffPowerLaw(); // destructor

    virtual SuperExpoCutOffPowerLaw* clone() const; // VIM : Tricks to create a copy constructor
  
    virtual double GetFlux(double x) const; // return the value of the flux for current parameters

    virtual void DerivativesFormulae(double x); // Used for the butterfly. 

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
    ClassDef(START::SuperExpoCutOffPowerLaw,1);
#endif
  };
}
#endif
