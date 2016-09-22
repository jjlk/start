//STL
#include <sstream>

// ROOT

// START
#include "TimeBin.hh"

// Utilities
#define DEBUG 0
#include "debugging.hh"

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::TimeBin)
#endif

/**
 * \brief Constructor
 */
START::TimeBin::TimeBin(double tmin, double tmax, double intflux)
   :fTimeMin(tmin),
   fTimeMax(tmax),
   fTimeMean((tmin+tmax)/2.),
   fLiveTime(0.),
   fExpectedOn(0.),
   fExpectedOff(0.),
   fExpectedExcess(0.),
   fOn(0.),
   fOff(0.),
   fExcess(0.),
   fAlpha(0.),
   fIntegratedFlux(intflux),
   fIntegratedFluxError(0.),
   fIntegratedFluxError1SigmaMinus(0.),
   fIntegratedFluxError1SigmaPlus(0.),
   fIntegratedFluxError3SigmaMinus(0.),
   fIntegratedFluxError3SigmaPlus(0.),
   fIsUpperLimit(false)
{

}

/**
 * \brief Destructor
 */
START::TimeBin::~TimeBin()
{

}

/**
 * \brief Copy constructor
 */
START::TimeBin::TimeBin(TimeBin const &TimeBinCopy)
  :TObject(TimeBinCopy),
   fTimeMin(TimeBinCopy.fTimeMin),
   fTimeMax(TimeBinCopy.fTimeMax),
   fTimeMean(TimeBinCopy.fTimeMean),
   fLiveTime(TimeBinCopy.fLiveTime),
   fExpectedOn(TimeBinCopy.fExpectedOn),
   fExpectedOff(TimeBinCopy.fExpectedOff),
   fExpectedExcess(TimeBinCopy.fExpectedExcess),
   fOn(TimeBinCopy.fOn),
   fOff(TimeBinCopy.fOff),
   fExcess(TimeBinCopy.fExcess),
   fAlpha(TimeBinCopy.fAlpha),
   fIntegratedFlux(TimeBinCopy.fIntegratedFlux),
   fIntegratedFluxError(TimeBinCopy.fIntegratedFluxError),
   fIntegratedFluxError1SigmaMinus(TimeBinCopy.fIntegratedFluxError1SigmaMinus),
   fIntegratedFluxError1SigmaPlus(TimeBinCopy.fIntegratedFluxError1SigmaPlus),
   fIntegratedFluxError3SigmaMinus(TimeBinCopy.fIntegratedFluxError3SigmaMinus),
   fIntegratedFluxError3SigmaPlus(TimeBinCopy.fIntegratedFluxError3SigmaPlus),
   fIsUpperLimit(TimeBinCopy.fIsUpperLimit)
{

}

/**
 * \brief Assignment operator
 */
START::TimeBin &START::TimeBin::operator=(TimeBin const &TimeBinCopy)
{
  if(this != &TimeBinCopy) {
    fTimeMin = TimeBinCopy.fTimeMin;
    fTimeMax = TimeBinCopy.fTimeMax;
    fTimeMean = TimeBinCopy.fTimeMean;
    fLiveTime = TimeBinCopy.fLiveTime;
    fExpectedOn = TimeBinCopy.fExpectedOn;
    fExpectedOff = TimeBinCopy.fExpectedOff;
    fExpectedExcess = TimeBinCopy.fExpectedExcess;
    fOn = TimeBinCopy.fOn;
    fOff = TimeBinCopy.fOff;
    fExcess = TimeBinCopy.fExcess;
    fAlpha = TimeBinCopy.fAlpha;
    fIntegratedFlux = TimeBinCopy.fIntegratedFlux;
    fIntegratedFluxError = TimeBinCopy.fIntegratedFluxError;
    fIntegratedFluxError1SigmaMinus = TimeBinCopy.fIntegratedFluxError1SigmaMinus;
    fIntegratedFluxError1SigmaPlus = TimeBinCopy.fIntegratedFluxError1SigmaPlus;
    fIntegratedFluxError3SigmaMinus = TimeBinCopy.fIntegratedFluxError3SigmaMinus;
    fIntegratedFluxError3SigmaPlus = TimeBinCopy.fIntegratedFluxError3SigmaPlus;
    fIsUpperLimit = TimeBinCopy.fIsUpperLimit;
  }
  return (*this);
}

/**
 * \brief Print TimeBin info
 */  
void START::TimeBin::Print(Option_t *option) const
{
  PrintTimeBin();
}

/**
 * \brief Print TimeBin info
 */
void START::TimeBin::PrintTimeBin(std::ostream &os) const
{
  std::ostringstream binheader;
  binheader.setf(std::ios::fixed);
  binheader.width(10);
  binheader.precision(2);
  binheader << fTimeMin;
  binheader.width(10);
  binheader.precision(2);
  binheader << fTimeMean;
  binheader.width(10);
  binheader.precision(2);
  binheader << fTimeMax;
  binheader.width(10);
  binheader.precision(2);
  binheader << fLiveTime;
  binheader.width(10);
  binheader.precision(2);
  binheader << fOn;
  binheader.width(10);
  binheader.precision(2);
  binheader << fOff;
  binheader.width(10);
  binheader.precision(3);
  binheader << fAlpha;
  binheader.width(10);
  binheader.precision(2);
  binheader << fExcess;
  binheader.width(10);
  binheader.precision(2);
  binheader << STARTUtils::LiMaSignificance(fOn,fOff,fAlpha);
  binheader.width(10);
  binheader.precision(2);
  binheader << fExpectedOn;
  binheader.width(10);
  binheader.precision(2);
  binheader << fExpectedOff;
  binheader.width(11);
  binheader.precision(2);
  binheader << fExpectedExcess;
  binheader.width(14);
  binheader.precision(3);
  binheader << std::scientific << fIntegratedFlux;
  binheader.width(14);
  binheader.precision(3);
  binheader << std::scientific << fIntegratedFluxError;
  binheader.width(14);
  binheader.precision(3);
  binheader << std::scientific << fIntegratedFluxError1SigmaMinus;
  binheader.width(14);
  binheader.precision(3);
  binheader << std::scientific << fIntegratedFluxError1SigmaPlus;
  binheader.width(14);
  binheader.precision(3);
  binheader << std::scientific << fIntegratedFluxError3SigmaMinus;
  binheader.width(14);
  binheader.precision(3);
  binheader << std::scientific << fIntegratedFluxError3SigmaPlus;
  // JLK error
  os << binheader.str().c_str() << std::endl;
}


