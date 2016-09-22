#ifndef _SUM3HYPOTHESIS_
#define _SUM3HYPOTHESIS_

#include "Hypothesis.hh"

namespace START {
  /**
   * \brief Used to sum hypothesis \f$H_1\f$ and \f$H_2\f$ by building the Hypothesis : 
   * \f[
   * \Phi(E) = \phi_{H_1}(E) + \phi_{H_2}(E)
   * \f]
   * or to compare two Hypothesis \f$H_1\f$ and \f$H_2\f$ by minimizing
   * \f[
   * \Phi(E) =  \chi \phi_{H_1}(E) + (1-\chi) \phi_{H_2}(E)
   * \f]
   */
  class Sum3Hypothesis : public Hypothesis
  {
  
  public:
    Sum3Hypothesis(); // default constructor

    Sum3Hypothesis(TString name, const Hypothesis &hypo1, const Hypothesis &hypo2, const Hypothesis &hypo3, bool comparison=false); // constructor taking two hypothesis

    Sum3Hypothesis(Sum3Hypothesis const &StoreCopy); // Copy constructor

    virtual ~Sum3Hypothesis(); // destructor

    virtual Sum3Hypothesis* clone() const; // Tricks to create a copy constructor
  
    virtual double GetFlux(double x) const; // return the value of the flux for current parameters

    virtual void DerivativesFormulae(double x); // 

    virtual void InitParametersCaracteristics(); // Initialization of parameters names and units

  private:
    Hypothesis *fHypo1, *fHypo2, *fHypo3;
    bool fiscomparison;

#ifdef CONSTRUCT_LIBRARY
  public:
    ClassDef(START::Sum3Hypothesis,1);
#endif
  };

}
#endif
