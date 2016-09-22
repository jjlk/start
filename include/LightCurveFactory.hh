#ifndef _LIGHTCURVEFACTORY_
#define _LIGHTCURVEFACTORY_

// STL
#include <vector>
#include <string>

// ROOT
#include <TObject.h>

// START
#include "Config.hh"
#include "TimeBinVector.hh"

namespace START {

  
class Hypothesis;
class Band;

  /**
   * \brief Fill TimeBinVector
   *
   * \author HAP-Fr team
   */

  class LightCurveFactory : public TObject
  {

  public:

    /**
     * \brief different binning
     */
    typedef enum {RunByRun, MinuteByMinute, HourByHour, NightByNight, DayByDay, 
		  WeekByWeek, MonthByMonth, YearByYear, GivenTimeInterval, 
		  UserTimeIntervals, PeriodByPeriod} TimeCuttingType;

    LightCurveFactory() {}; // Needed by ROOT
    LightCurveFactory(std::vector<START::Hypothesis*> HypothesisArray, Config &Configuration, bool verbose=false);
    LightCurveFactory(Hypothesis &Hypo, Config &Configuration, bool verbose=false);

    ~LightCurveFactory();

    void SetAndBuildLightCurveData(TimeCuttingType TimeCutting=RunByRun);
    void BuildGivenTimeIntervalLightCurve(Hypothesis &hypo, TimeBinVector &TDataContainer);

    std::vector<std::pair<double, double> > GetAllTimeIntervals(const std::vector<START::Band> &BandArray);
  
  private:

    std::vector<START::Hypothesis*> fHypothesisArray; ///< Hypothesis array
    Config *fConfig;

    bool fverbose; ///< verbosity level

    TimeCuttingType fTimeCutting; ///< time type cutting

    std::string fTimeCuttingName;

#ifdef CONSTRUCT_LIBRARY
  public:
    ClassDef(START::LightCurveFactory,1);
#endif
  };
}
#endif
