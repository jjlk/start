#ifndef _SUMHYPOTHESIS_
#define _SUMHYPOTHESIS_

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
  class SumHypothesis : public Hypothesis
  {
  
  public:
    SumHypothesis(); // default constructor

    SumHypothesis(TString name, const Hypothesis &hypo1, const Hypothesis &hypo2, bool comparison=false); // constructor taking two hypothesis

    SumHypothesis(SumHypothesis const &StoreCopy); // Copy constructor

    virtual ~SumHypothesis(); // destructor

    virtual SumHypothesis* clone() const; // Tricks to create a copy constructor
  
    virtual double GetFlux(double x) const; // return the value of the flux for current parameters

    virtual void DerivativesFormulae(double x); // 

    virtual void InitParametersCaracteristics(); // Initialization of parameters names and units

  private:
    Hypothesis *fHypo1, *fHypo2;
    bool fiscomparison;

#ifdef CONSTRUCT_LIBRARY
  public:
    ClassDef(START::SumHypothesis,1);
#endif
  };

}
#endif
