//STL
#include <sstream>

// ROOT

// START
#include "TimeBinVector.hh"

// Utilities
#define DEBUG 0
#include "debugging.hh"

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::TimeBinVector)
#endif

/**
 * \brief Constructor
 */
START::TimeBinVector::TimeBinVector()
   :fLiveTime(0.),
   fExpectedOn(0.),
   fExpectedOff(0.),
   fExpectedExcess(0.),
   fOn(0.),
   fOff(0.),
   fExcess(0.)
{
  tbin.clear();
}

/**
 * \brief Constructor
 */
START::TimeBinVector::TimeBinVector(std::vector<TimeBin> const &TimeBinArray)
   :fLiveTime(0.),
   fExpectedOn(0.),
   fExpectedOff(0.),
   fExpectedExcess(0.),
   fOn(0.),
   fOff(0.),
   fExcess(0.)
{
  tbin.clear();
  tbin = TimeBinArray;

  for(std::vector<TimeBin>::const_iterator bin=TimeBinArray.begin(); bin!=TimeBinArray.end(); ++bin) {
    fLiveTime+=bin->GetLiveTime();
    fExpectedOn+=bin->GetExpectedOn();
    fExpectedOff+=bin->GetExpectedOff();
    fExpectedExcess+=bin->GetExpectedExcess();
    fOn+=bin->GetOn();
    fOff+=bin->GetOff();
    fExcess+=bin->GetExcess();
  }

}

/**
 * \brief Destructor
 */
START::TimeBinVector::~TimeBinVector()
{

}

/**
 * \brief Copy constructor
 */
START::TimeBinVector::TimeBinVector(TimeBinVector const &TimeBinVectorCopy)
  :TObject(TimeBinVectorCopy),
   fLiveTime(TimeBinVectorCopy.fLiveTime),
   fExpectedOn(TimeBinVectorCopy.fExpectedOn),
   fExpectedOff(TimeBinVectorCopy.fExpectedOff),
   fExpectedExcess(TimeBinVectorCopy.fExpectedExcess),
   fOn(TimeBinVectorCopy.fOn),
   fOff(TimeBinVectorCopy.fOff),
   fExcess(TimeBinVectorCopy.fExcess)
{
  tbin = TimeBinVectorCopy.tbin;
}

/**
 * \brief Assignment operator
 */
START::TimeBinVector &START::TimeBinVector::operator=(TimeBinVector const &TimeBinVectorCopy)
{
  if(this != &TimeBinVectorCopy) {
    fLiveTime = TimeBinVectorCopy.fLiveTime;
    fExpectedOn = TimeBinVectorCopy.fExpectedOn;
    fExpectedOff = TimeBinVectorCopy.fExpectedOff;
    fExpectedExcess = TimeBinVectorCopy.fExpectedExcess;
    fOn = TimeBinVectorCopy.fOn;
    fOff = TimeBinVectorCopy.fOff;
    fExcess = TimeBinVectorCopy.fExcess;
    tbin = TimeBinVectorCopy.tbin;
  }
  return (*this);
}

/**
 * \brief Print TimeBinVector info
 */  
void START::TimeBinVector::Print(Option_t *option) const
{
  PrintTimeBinVector();
}

/**
 * \brief Print TimeBinVector info
 */
void START::TimeBinVector::PrintTimeBinVector(std::ostream &os) const
{
  std::ostringstream binvectorheader;
  binvectorheader.setf(std::ios::fixed);
  binvectorheader << "TotalLiveTime(h) ";
  binvectorheader.precision(2);
  binvectorheader << GetLiveTime() << " ";
  binvectorheader << "TotalOn ";
  binvectorheader.precision(2);
  binvectorheader << GetOn() << " ";
  binvectorheader << "TotalOff ";
  binvectorheader.precision(2);
  binvectorheader << GetOff() << " ";
  binvectorheader << "TotalExcess ";
  binvectorheader.precision(2);
  binvectorheader << GetExcess() << " ";
  binvectorheader << "TotalExpOn ";
  binvectorheader.precision(2);
  binvectorheader << GetExpectedOn() << " ";
  binvectorheader << "TotalExpOff ";
  binvectorheader.precision(2);
  binvectorheader << GetExpectedOff() << " ";
  binvectorheader << "TotalExpExcess ";
  binvectorheader.precision(2);
  binvectorheader << GetExpectedExcess() << " ";
  os << binvectorheader.str().c_str() << std::endl;

  std::ostringstream binheader;
  binheader.setf(std::ios::fixed);
  binheader.width(10);
  binheader << " Tmin ";
  binheader.width(10);
  binheader << " Tmean ";
  binheader.width(10);
  binheader << " Tmax ";
  binheader.width(12);
  binheader << " LiveTime ";
  binheader.width(8);
  binheader << " ON ";
  binheader.width(10);
  binheader << " Off ";
  binheader.width(11);
  binheader << " Alpha ";
  binheader.width(10);
  binheader << " Excess ";
  binheader.width(10);
  binheader << " Sigma ";
  binheader.width(10);
  binheader << " ExpOn ";
  binheader.width(10);
  binheader << " ExpOff ";
  binheader.width(13);
  binheader << " ExpExcess ";
  binheader.width(10);
  binheader << " Flux ";
  binheader.width(15);
  binheader << " GaussErr ";
  binheader.width(16);
  binheader << " RolkeErr1- ";
  binheader.width(14);
  binheader << " RolkeErr1+ ";
  binheader.width(14);
  binheader << " RolkeErr3- ";
  binheader.width(14);
  binheader << " RolkeErr3+ ";
  os << binheader.str().c_str() << std::endl;

  for(std::vector<TimeBin>::const_iterator bin=tbin.begin(); bin!=tbin.end(); ++bin) {
    bin->PrintTimeBin(os);
  }
}

double START::TimeBinVector::GetLiveTime() const {
  double livetime(0.);
  for(std::vector<TimeBin>::const_iterator bin=tbin.begin(); bin!=tbin.end(); ++bin) {
    livetime+=bin->GetLiveTime();
  }
  return livetime;
}

double START::TimeBinVector::GetExpectedOn() const {
  double expectedon(0.);
  for(std::vector<TimeBin>::const_iterator bin=tbin.begin(); bin!=tbin.end(); ++bin) {
    expectedon+=bin->GetExpectedOn();
  }
  return expectedon;
}

double START::TimeBinVector::GetExpectedOff() const {
  double expectedoff(0.);
  for(std::vector<TimeBin>::const_iterator bin=tbin.begin(); bin!=tbin.end(); ++bin) {
    expectedoff+=bin->GetExpectedOff();
  }
  return expectedoff;
}

double START::TimeBinVector::GetExpectedExcess() const {
  double expectedexcess(0.);
  for(std::vector<TimeBin>::const_iterator bin=tbin.begin(); bin!=tbin.end(); ++bin) {
    expectedexcess+=bin->GetExpectedExcess();
  }
  return expectedexcess;
}

double START::TimeBinVector::GetOn() const {
  double on(0.);
  for(std::vector<TimeBin>::const_iterator bin=tbin.begin(); bin!=tbin.end(); ++bin) {
    on+=bin->GetOn();
  }
  return on;
}

double START::TimeBinVector::GetOff() const {
  double off(0.);
  for(std::vector<TimeBin>::const_iterator bin=tbin.begin(); bin!=tbin.end(); ++bin) {
    off+=bin->GetOff();
  }
  return off;
}

double START::TimeBinVector::GetExcess() const {
  double excess(0.);
  for(std::vector<TimeBin>::const_iterator bin=tbin.begin(); bin!=tbin.end(); ++bin) {
    excess+=bin->GetExcess();
  }
  return excess;
}


