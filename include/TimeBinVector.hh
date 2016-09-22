#ifndef _TIMEBINVECTOR_
#define _TIMEBINVECTOR_

// STL
#include <iostream>
#include <vector>

// ROOT
#include <TObject.h>

// START
#include "TimeBin.hh"

namespace START {
  /**
   * \brief Store several TimeBin in a vector and contain global informations
   * about the total timed binned data
   *
   * \author HAP-Fr team
   */

  class TimeBinVector : public TObject
  {

  public:

    TimeBinVector();
    TimeBinVector(std::vector<START::TimeBin> const &TimeBinArray);

    ~TimeBinVector(); 
    TimeBinVector(TimeBinVector const &TimeBinVectorCopy); // Copy constructor
    TimeBinVector &operator=(TimeBinVector const &TimeBinVectorCopy); // assignment operator

    std::vector<START::TimeBin> tbin;

    virtual void Print(Option_t *option="") const;
    void PrintTimeBinVector(std::ostream &os=std::cout) const;

    // Set functions
    void SetLiveTime(double livetime) {fLiveTime=livetime;};
    void SetExpectedOn(double expectedon) {fExpectedOn=expectedon;};
    void SetExpectedOff(double expectedoff) {fExpectedOff=expectedoff;};
    void SetExpectedExcess(double expectedexcess) {fExpectedExcess=expectedexcess;};
    void SetOn(double measuredon) {fOn=measuredon;};
    void SetOff(double measuredoff) {fOff=measuredoff;};
    void SetExcess(double measuredexcess) {fExcess=measuredexcess;};

    // Get functions
    double GetLiveTime() const;
    double GetExpectedOn() const;
    double GetExpectedOff() const;
    double GetExpectedExcess() const;
    double GetOn() const;
    double GetOff() const;
    double GetExcess() const;

  private:
    double fLiveTime;
    double fExpectedOn;
    double fExpectedOff;
    double fExpectedExcess;
    double fOn;
    double fOff;
    double fExcess;

#ifdef CONSTRUCT_LIBRARY
  public:
    ClassDef(START::TimeBinVector,1);
#endif
  
  };

}

#endif
