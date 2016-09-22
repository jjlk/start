#ifndef _SMOOTHBROKENPOWERLAW_
#define _SMOOTHBROKENPOWERLAW_

#include "Hypothesis.hh"
namespace START {
  /**
   * \brief Smooth BrokenPowerLaw Hypothesis with 5 parameters \f$\phi_0\f$, \f$E_{break}\f$, \f$\beta\f$, \f$\Gamma_1\f$ and \f$\Gamma_2\f$ :
   * \f[
   * \Phi(E) = \phi_{0} \left (\frac{E}{E_{Ref}} \right)^{-\Gamma_{1}} \left (1 + \left (\frac{E}{E_{break}} \right)^{\frac{\Gamma_{2}-\Gamma_{1}}{\beta}} \right)^{-\beta}
   * \f]
   */
  class SmoothBrokenPowerLaw : public Hypothesis
  {

  public:

    SmoothBrokenPowerLaw(); // default constructor

    SmoothBrokenPowerLaw(TString name, SpectralType type=Differential, double phi0=1.e-12, double ebreak=1.,
			 double beta=0.1, double gamma1=2., double gamma2=2.5, double phi0err=0.,
			 double ebreakerr=0., double betaerr=0., double gamma1err=0., double gamma2err=0.); 

    SmoothBrokenPowerLaw(SmoothBrokenPowerLaw const &StoreCopy); // Copy constructor

    virtual ~SmoothBrokenPowerLaw(); // destructor

    virtual SmoothBrokenPowerLaw* clone() const; // VIM : Tricks to create a copy constructor
  
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
    ClassDef(START::SmoothBrokenPowerLaw,1);
#endif
  };
}
#endif
