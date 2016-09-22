#ifndef _POWERLAW_
#define _POWERLAW_

#include "Hypothesis.hh"

namespace START {
  /**
   * \brief Power Law Hypothesis with 2 parameters \f$\phi_{0}\f$ and \f$\Gamma\f$ :
   * \f[
   * \Phi(E) =  \phi_0 \left (\frac{E}{E_{Ref}} \right)^{-\Gamma}
   * \f]
   */
  class PowerLaw : public Hypothesis
  {
  
  public:
    PowerLaw(); // default constructor

    PowerLaw(TString name, SpectralType type=Differential, double phi0=1.e-12, double gamma=2.5,double phi0err=0., double gammaerr=0.);
    // constructor used for minimization with initial values for parameters and errors

    PowerLaw(PowerLaw const &StoreCopy); // Copy constructor

    virtual ~PowerLaw(); // destructor

    virtual PowerLaw* clone() const; // VIM : Tricks to create a copy constructor
  
    virtual double GetFlux(double x) const; // return the value of the flux for current parameters

    //virtual double GetFluxFitParams(double x) const; // return the value of the flux for fitted parameters

    virtual double GetMeanBinEnergy(double x1,double x2) const; // Compute the mean energy bin for energy x1 and x2

    virtual void DerivativesFormulae(double x); //  Used for the butterfly. return a vector composed with first derivatives of flux for param 0,1,... n at energy x

    virtual void InitParametersCaracteristics(); // Initialization of parameters names and units

    /* Added for convenience */

    double GetPhi0() const {return fparam[0];}; ///< Get the flux normalization
    double GetPhi0Err() const {return fparamerr[0];}; ///< Get the flux normalization's error
    double GetGamma() const {return fparam[1];}; ///< Get the spectral index
    double GetGammaErr() const {return fparamerr[1];}; ///< Get the spectral index's error

    void SetPhi0(double phi0) {fparam[0]=phi0;}; ///< Set the flux normalization
    void SetPhi0err(double phi0err) {fparamerr[0]=phi0err;};  ///< Set the flux normalization's error
    void SetGamma(double gamma) {fparam[1]=gamma;}; ///< Set the spectral index
    void SetGammaerr(double gammaerr) {fparam[1]=gammaerr;}; ///< Set the spectral index's error

#ifdef CONSTRUCT_LIBRARY
  public:
    ClassDef(START::PowerLaw,1);
#endif
  };

}
#endif
