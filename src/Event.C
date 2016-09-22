//STL
#include <sstream>

// ROOT

// START
#include "Event.hh"

// Utilities
#define DEBUG 0
#include "debugging.hh"

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::Event)
#endif

/**
 * \brief Constructor
 */
START::Event::Event(bool ison, double energy, double timmjd)
  :fTimeMJD(timmjd),
  fEnergy(energy),
  fIsOn(ison)
{
  fIsOff = !fIsOn;
}

/**
 * \brief Destructor
 */
START::Event::~Event()
{

}

/**
 * \brief Copy constructor
 */
START::Event::Event(Event const &EventCopy)
  :TObject(EventCopy),
   fTimeMJD(EventCopy.fTimeMJD),
   fEnergy(EventCopy.fEnergy),
   fIsOn(EventCopy.fIsOn),
   fIsOff(EventCopy.fIsOff)
{

}

/**
 * \brief Assignment operator
 */
START::Event &START::Event::operator=(Event const &EventCopy)
{
  if(this != &EventCopy) {
    fTimeMJD = EventCopy.fTimeMJD;
    fEnergy = EventCopy.fEnergy;
    fIsOn = EventCopy.fIsOn;
    fIsOff = EventCopy.fIsOff;
  }
  return (*this);
}

/**
 * \brief Print Event info
 */  
void START::Event::Print(Option_t *option) const
{
  PrintEvent();
}

/**
 * \brief Print Event info
 */
void START::Event::PrintEvent(std::ostream &os) const
{
  std::ostringstream evtheader;
  if(fIsOn)
    evtheader << "OnEvt: ";
  else
    evtheader << "OffEvt: ";
  evtheader << "Energy=";
  evtheader.precision(3); evtheader.setf(std::ios::fixed);
  evtheader << fEnergy << " TeV, TimeMJD=";
  evtheader.precision(3); evtheader.setf(std::ios::fixed);
  evtheader<< fTimeMJD;
  os << evtheader.str().c_str() << std::endl;
}


