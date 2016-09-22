#ifndef _ENERGYBIN_
#define _ENERGYBIN_

// STL
#include <iostream>
#include <vector>
#include <utility>
#include <cstdlib>

// ROOT
#include <TObject.h>
#include "Math/InterpolationTypes.h"
#include "Math/Interpolator.h"

// START
#include "Event.hh"

namespace START {
  /**
   * \brief Stock data of one bin in energy (on, off, emin, emax, emid, time, significance...)
   *
   * \author HAP-Fr team
   */

  class EnergyBin : public TObject
  {

  public:

    EnergyBin();
    EnergyBin(double emin, double emax);
    ~EnergyBin();
  
    EnergyBin(EnergyBin const &EnergyBinCopy);
    EnergyBin &operator=(EnergyBin const &EnergyBinCopy);

    virtual void Print(Option_t *option="") const;
    void PrintEnergyBin(std::ostream &os=std::cout) const;

    double GetInterpolatedPartialIntegral(double E);

    void ClearEnergyBinInfo();
    void AddInfoFromEBin(const EnergyBin &AnotherEBin);
			      
    inline unsigned int    GetKeepBin() const    {return fkeepbin;} ///< Get fkeepbin
    inline double GetEmin() const       {return fe_min;} ///< Get fe_min
    inline double GetEmax() const       {return fe_max;} ///< Get fe_max
    inline double GetEmid() const       {return fe_mid;} ///< Get fe_mid
    inline double GetEmean() const      {return fe_mean;} ///< Get fe_mean
    inline double GetOn() const         {return fon;} ///< Get fon
    inline double GetOff() const        {return foff;} ///< Get foff
    inline double GetAcceff() const     {return facceff;} ///< Get facceff
    inline double GetSth() const        {return fsth;} ///< Get fsth
    inline double GetOnFitted() const   {return fon_fitted;} ///< Get fon_fitted
    inline double GetOffFitted() const  {return foff_fitted;} ///< Get foff_fitted
    inline double GetAlpha() const      {return falpha;} ///< Get falpha
    inline double GetLiveTime() const   {return flivetime;} ///< Get flivetime
    std::pair<std::vector<double>, std::vector<double> > GetPartialIntegral() const {return fpartialintegral;} ///< Get fpartialintegral
    std::vector<START::Event> GetEvents() const {return fEvents;}

    void SetKeepBin(int keepbin)         {fkeepbin = keepbin;} ///< Set fkeepbin
    void SetEmin(double e_min)           {fe_min = e_min;} ///< Set fe_min
    void SetEmax(double e_max)           {fe_max = e_max;} ///< Set fe_max
    void SetEmid()                       {fe_mid = 0.5*(fe_min+fe_max);} ///< Set fe_mid
    void SetEmid(double e_mid)           {fe_mid = e_mid;} ///< Set fe_mid
    void SetEmean(double e_mean)         {fe_mean = e_mean;} ///< Set fe_mean
    void SetOn(double on)                {fon = on;} ///< Set fon
    void SetOff(double off)              {foff = off;} ///< Set foff
    void SetAcceff(double acceff)        {facceff = acceff;} ///< Set facceff
    void SetSth(double sth)              {fsth = sth;} ///< Set fsth
    void SetOnFitted(double on_fitted)   {fon_fitted = on_fitted;} ///< Set fon_fitted
    void SetOffFitted(double off_fitted) {foff_fitted = off_fitted;} ///< Set goff_fitted
    void SetAlpha(double alpha)          {falpha = alpha;} ///< Set falpha
    void SetLiveTime(double livetime)    {flivetime = livetime;} ///< Set the livetime in the bin
    void SetPartialIntegral(std::pair< std::vector<double>, std::vector<double> > partialintegral) {fpartialintegral=partialintegral;} ///< Set fpartialintegral

    void SetGSLInterpolatorForPartialIntegral(std::pair<std::vector<double>, std::vector<double> > partialintegral); ///< Set GSL interpolator

    void SetEvents(std::vector<START::Event> Events) {fEvents = Events;}; ///< Set Events
    void AddEvent(START::Event Evt) {fEvents.push_back(Evt);}; ///< Add an Event

  private:

    unsigned int fkeepbin; ///< Equal 1 if bin is used and 0 otherwise
    double fe_min; ///< Minimal limit of the bin in energy
    double fe_mid; ///< Center of the bin's energy
    double fe_max; ///< Superior limit of the bin in energy 
    double fe_mean; ///< Center of the energy bin after by taking into account the spectral fit
    double fon; ///< Number of ON event
    double foff; ///< Number of OFF event
    double facceff; ///< Ratio of Time on over Time oFF
    double fsth; ///< Number of expected excess
    double falpha; ///< Alpha (usefull in the case you want to stack the result, because alpha will change from bin to bin)
    double fon_fitted;  ///< fitted ON
    double foff_fitted; ///< fitted OFF
    double flivetime;   ///< livetime of the band
    // JLK add this : 
    /*
      fsigma; ///< significance
      fexcess; ///< excess
    */
  
    ROOT::Math::Interpolation::Type fInterpolTypeForBin; ///< Type of interpolation
  
    std::pair<std::vector<double>, std::vector<double> > fpartialintegral; ///< contain energy and integrated resolution
  
    ROOT::Math::Interpolator *fInterpolPartialIntegral; //! ///< GSL interpolator

    std::vector<START::Event> fEvents; ///< Contains Evt information

#ifdef CONSTRUCT_LIBRARY
  public:
    ClassDef(START::EnergyBin,1);
#endif
   
  };
}
#endif
