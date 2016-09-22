#ifndef _EXPOCUTOFFPOWERLAW_
#define _EXPOCUTOFFPOWERLAW_

#include "Hypothesis.hh"

namespace START {
  /**
   * \brief  PowerLaw with exponential cut-off Hypothesis with 3 parameters \f$\phi_0\f$, \f$\Gamma\f$ and \f$\beta\f$ :
   * \f[
   * \Phi(E) = \phi_0 \left (\frac{E}{E_{Ref}} \right)^{-\Gamma} \exp{(-\beta E)}
   * \f]
   */
  class ExpoCutOffPowerLaw : public Hypothesis
  {

  public:

    ExpoCutOffPowerLaw(); // default constructor

    ExpoCutOffPowerLaw(TString name, SpectralType type=Differential, double phi0=1.e-12,double gamma=2.5,double beta=0.1, // constructor used for minimization with initial 
		       double phi0err=0., double gammaerr=0., double betaerr=0.); //values for parameters and errors

    ExpoCutOffPowerLaw(ExpoCutOffPowerLaw const &StoreCopy); // Copy constructor

    virtual ~ExpoCutOffPowerLaw(); // destructor

    virtual ExpoCutOffPowerLaw* clone() const; // VIM : Tricks to create a copy constructor
  
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
    ClassDef(START::ExpoCutOffPowerLaw,1);
#endif
  };
}
#endif
