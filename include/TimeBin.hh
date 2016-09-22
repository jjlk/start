#ifndef _TIMEBIN_
#define _TIMEBIN_

// STL
#include <iostream>

// ROOT
#include <TObject.h>

// START
#include "STARTUtils.hh"

namespace START {
  /**
   * \brief Stock data of one bin in time (time, excess, On, Off, and integrated flux)
   *
   * \author HAP-Fr team
   */

  class TimeBin : public TObject
  {

  public:

    TimeBin(double TimeMin=0., double TimeMax=0., double intflux=0.);
    ~TimeBin();

    TimeBin(TimeBin const &TimeBinCopy);
    TimeBin &operator=(TimeBin const &TimeBinCopy);
  
    // Print functions
    virtual void Print(Option_t *option="") const;
    void PrintTimeBin(std::ostream &os=std::cout) const;

    // Set functions
    void SetTimeMin(double timemin) {fTimeMin=timemin;};
    void SetTimeMean(double timemean) {fTimeMean=timemean;};
    void SetTimeMax(double timemax) {fTimeMax=timemax;};
    void SetLiveTime(double livetime) {fLiveTime=livetime;};
    void SetExpectedOn(double expectedon) {fExpectedOn=expectedon;};
    void SetExpectedOff(double expectedoff) {fExpectedOff=expectedoff;};
    void SetExpectedExcess(double expectedexcess) {fExpectedExcess=expectedexcess;};
    void SetOn(double on) {fOn=on;};
    void SetOff(double off) {fOff=off;};
    void SetExcess(double excess) {fExcess=excess;};
    void SetAlpha(double alpha) {fAlpha=alpha;}; // mean alpha
    void SetIntegratedFlux(double integratedflux) {fIntegratedFlux=integratedflux;};
    void SetIntegratedFluxError(double integratedfluxerror) {fIntegratedFluxError=integratedfluxerror;};
    void SetIntegratedFluxError1SigmaMinus(double integratedfluxerror) {fIntegratedFluxError1SigmaMinus=integratedfluxerror;};
    void SetIntegratedFluxError1SigmaPlus(double integratedfluxerror) {fIntegratedFluxError1SigmaPlus=integratedfluxerror;};
    void SetIntegratedFluxError3SigmaMinus(double integratedfluxerror) {fIntegratedFluxError3SigmaMinus=integratedfluxerror;};
    void SetIntegratedFluxError3SigmaPlus(double integratedfluxerror) {fIntegratedFluxError3SigmaPlus=integratedfluxerror;};
    void SetIsUpperLimit(bool isupper) {fIsUpperLimit=isupper;};

    // Get functions
    double GetTimeMin() const {return fTimeMin;};
    double GetTimeMean() const {return fTimeMean;};
    double GetTimeMax() const {return fTimeMax;};
    double GetLiveTime() const {return fLiveTime;};
    double GetExpectedOn() const {return fExpectedOn;};
    double GetExpectedOff() const {return fExpectedOff;};
    double GetExpectedExcess() const {return fExpectedExcess;};
    double GetOn() const {return fOn;};
    double GetOff() const {return fOff;};
    double GetExcess() const {return fExcess;};
    double GetAlpha() const {return fAlpha;};
    double GetSignificance() const {return STARTUtils::LiMaSignificance(fOn,fOff,fAlpha);};
    double GetIntegratedFlux() const {return fIntegratedFlux;};
    double GetIntegratedFluxError() const {return fIntegratedFluxError;};
    double GetIntegratedFluxError1SigmaMinus() const {return fIntegratedFluxError1SigmaMinus;};
    double GetIntegratedFluxError1SigmaPlus() const {return fIntegratedFluxError1SigmaPlus;};
    double GetIntegratedFluxError3SigmaMinus() const {return fIntegratedFluxError3SigmaMinus;};
    double GetIntegratedFluxError3SigmaPlus() const {return fIntegratedFluxError3SigmaPlus;};
    bool GetIsUpperLimit() const {return fIsUpperLimit;};

  private:
    double fTimeMin;
    double fTimeMax;
    double fTimeMean;
    double fLiveTime;
    double fExpectedOn;
    double fExpectedOff;
    double fExpectedExcess;
    double fOn;
    double fOff;
    double fExcess;
    double fAlpha;
    double fIntegratedFlux;
    double fIntegratedFluxError; ///< classic error propagation
    double fIntegratedFluxError1SigmaMinus; ///< Rolke
    double fIntegratedFluxError1SigmaPlus; ///< Rolke
    double fIntegratedFluxError3SigmaMinus; ///< Rolke
    double fIntegratedFluxError3SigmaPlus; ///< Rolke
    bool fIsUpperLimit; ///< Rolke

#ifdef CONSTRUCT_LIBRARY
  public:
    ClassDef(START::TimeBin,1);
#endif
   
  };
}
#endif
