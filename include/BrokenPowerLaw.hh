#ifndef _BROKENPOWERLAW_
#define _BROKENPOWERLAW_

#include "Hypothesis.hh"

namespace START {
  /**
   * \brief Broken powerLaw Hypothesis with 4 parameters \f$\phi_0\f$, \f$E_{break}\f$ \f$\Gamma1\f$ and \f$\Gamma2\f$ :
   * \f[
   * \Phi(E) = \left\{
   *           \begin{array}{ll}
   *              \phi_0 \left (\frac{E}{E_{Ref}} \right)^{-\Gamma_1} & \qquad \mathrm{if}\quad E < E_{break}  \\
   *              \phi_0 \left (\frac{E_{break}}{E_{Ref}} \right)^{\Gamma2-\Gamma1}  \left (\frac{E}{E_{Ref}} \right)^{-\Gamma_2} & \qquad \mathrm{if}\quad E\geq E_{break} \\
   *           \end{array}
   *           \right.
   * \f]
   */
  class BrokenPowerLaw : public Hypothesis
  {

  public:

    BrokenPowerLaw(); // default constructor

    BrokenPowerLaw(TString name, SpectralType type=Differential, double phi0=1.e-12,double ebreak=10.,double gamma1=1., 
		   double gamma2=1., double phi0err=0.,double ebreakerr=0.,double gamma1err=0., double gamma2err=0.); 

    BrokenPowerLaw(BrokenPowerLaw const &StoreCopy); // Copy constructor

    virtual ~BrokenPowerLaw(); // destructor

    virtual BrokenPowerLaw* clone() const; // VIM : Tricks to create a copy constructor
  
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
    ClassDef(START::BrokenPowerLaw,1);
#endif
  };
}
#endif
