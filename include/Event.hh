#ifndef _EVENT_
#define _EVENT_

// STL
#include <iostream>

// ROOT
#include <TObject.h>

namespace START {
  /**
   * \brief Stock event information
   *
   * \author HAP-Fr team
   */

  class Event : public TObject
  {

  public:

    Event(bool ison=false, double energy=0., double timmjd=0.);
    ~Event();

    Event(Event const &EventCopy);
    Event &operator=(Event const &EventCopy);
  
    // Print functions
    virtual void Print(Option_t *option="") const;
    void PrintEvent(std::ostream &os=std::cout) const;

    // Set functions
    void SetTimeMJD(double timemjd) {fTimeMJD=timemjd;};
    void SetEnergy(double energy) {fEnergy=energy;};
    void SetIsOn(bool ison) {fIsOn=ison; fIsOff=!ison;};
    void SetIsOff(bool isoff) {fIsOff=isoff; fIsOn=!isoff;};

    // Get functions
    double GetEnergy() const {return fEnergy;};
    double GetTimeMJD() const {return fTimeMJD;};
    bool GetIsOn() const {return fIsOn;};
    bool GetIsOff() const {return fIsOff;};

  private:

    double fEnergy;
    double fTimeMJD;
    bool fIsOn;
    bool fIsOff;

#ifdef CONSTRUCT_LIBRARY
  public:
    ClassDef(START::Event,1);
#endif
   
  };
}
#endif
