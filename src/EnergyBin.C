// STL
#include <iomanip>
#include <sstream>

// START
#include "EnergyBin.hh"
#include "STARTUtils.hh"

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::EnergyBin)
#endif

/**
 * \brief Constructor
 */
START::EnergyBin::EnergyBin()
:fInterpolPartialIntegral(0)
{

  fkeepbin = 1;
  fe_min = 0; 
  fe_mid = 0; 
  fe_max = 0; 
  fe_mean=0;
  fon = 0; 
  foff = 0; 
  facceff = 0;
  fsth = 0; 
  falpha = 0;
  fon_fitted=0;
  foff_fitted=0;
  flivetime=0.;
  fInterpolTypeForBin = ROOT::Math::Interpolation::kLINEAR;
  fEvents.clear();
}

/**
 * \brief Constructor
 */
START::EnergyBin::EnergyBin(double emin, double emax)
:fInterpolPartialIntegral(0)
{

  fkeepbin = 1;
  fe_min = emin; 
  fe_mid = 0.5*(emax+emin); 
  fe_max = emax; 
  fe_mean = fe_mid;
  fon = 0; 
  foff = 0; 
  facceff = 0;
  fsth = 0; 
  falpha = 0;
  fon_fitted=0;
  foff_fitted=0;
  flivetime=0.;
  fInterpolTypeForBin = ROOT::Math::Interpolation::kLINEAR;
  fEvents.clear();
}

/**
 * \brief Destructor
 */
START::EnergyBin::~EnergyBin()
{

  if (fInterpolPartialIntegral!=0) delete fInterpolPartialIntegral; 
  fInterpolPartialIntegral = 0;

}

/**
 * \brief Copy constructor
 */
START::EnergyBin::EnergyBin(EnergyBin const &EnergyBinCopy)
  :TObject(EnergyBinCopy),
   fkeepbin(EnergyBinCopy.fkeepbin),
   fe_min(EnergyBinCopy.fe_min),
   fe_mid(EnergyBinCopy.fe_mid),
   fe_max(EnergyBinCopy.fe_max),
   fe_mean(EnergyBinCopy.fe_mean),
   fon(EnergyBinCopy.fon),
   foff(EnergyBinCopy.foff),
   facceff(EnergyBinCopy.facceff),
   fsth(EnergyBinCopy.fsth),
   falpha(EnergyBinCopy.falpha),
   fon_fitted(EnergyBinCopy.fon_fitted),
   foff_fitted(EnergyBinCopy.foff_fitted),
   flivetime(EnergyBinCopy.flivetime),
   fInterpolTypeForBin(EnergyBinCopy.fInterpolTypeForBin),
   fpartialintegral(EnergyBinCopy.fpartialintegral),
   fInterpolPartialIntegral(0), // JLK
   fEvents(EnergyBinCopy.fEvents)
{

  if (fpartialintegral.first.size()>4 && fpartialintegral.second.size()>4) {
    fInterpolPartialIntegral = new ROOT::Math::Interpolator(fpartialintegral.first,fpartialintegral.second,fInterpolTypeForBin);
  }
  else {
    fInterpolPartialIntegral=0;
  }
}

/**
 * \brief Assignment operator
 */
START::EnergyBin &START::EnergyBin::operator=(EnergyBin const &EnergyBinCopy)
{
  if(this != &EnergyBinCopy) {

    fe_min=EnergyBinCopy.fe_min;
    fe_mid=EnergyBinCopy.fe_mid;
    fe_max=EnergyBinCopy.fe_max;
    fe_mean=EnergyBinCopy.fe_mean;
    fkeepbin=EnergyBinCopy.fkeepbin;
    fon=EnergyBinCopy.fon;
    foff=EnergyBinCopy.foff;
    facceff=EnergyBinCopy.facceff;
    fsth=EnergyBinCopy.fsth;
    falpha=EnergyBinCopy.falpha;
    flivetime=EnergyBinCopy.flivetime;
    fpartialintegral=EnergyBinCopy.fpartialintegral;
    fInterpolTypeForBin=EnergyBinCopy.fInterpolTypeForBin;
    fon_fitted=EnergyBinCopy.fon_fitted;
    foff_fitted=EnergyBinCopy.foff_fitted;    
    fInterpolPartialIntegral=0; // JLK
    fEvents=EnergyBinCopy.fEvents;

    if (fpartialintegral.first.size()>4 && fpartialintegral.second.size()>4 ) {
      fInterpolPartialIntegral = new ROOT::Math::Interpolator(fpartialintegral.first,fpartialintegral.second,fInterpolTypeForBin);
    }
    else {
      fInterpolPartialIntegral=0;
    }

  }

  return (*this);
}

void START::EnergyBin::Print(Option_t *) const {

  PrintEnergyBin();

}

/**
 * \brief Print caracteristics of the energy bin
 *
 * Print : 
 * <ul> 
 *    <li> inferior limit of the bin </li>
 *    <li> mean of the bin </li>
 *    <li> superior limit of the bin </li>
 *    <li> number of ON events </li>
 *    <li> number of OFF events </li>
 *    <li> \f$\alpha\f$ defined as the ratio of \f$T_{ON}\f$ and \f$T_{OFF}\f$ </li>
 *    <li> number of excess events defined by \f$ N_{ON} - \alpha N_{OFF}\f$ </li>
 *    <li> Li-Ma significance of the bin </li>
 *    <li> live time of the bin </li>
 *    <li> effective area at mean energy bin </li>
 *    <li> expected number of ON events </li>
 *    <li> expected number of OFF events </li>
 *    <li> expected number of excess events </li>
 * </ul>
 */
void START::EnergyBin::PrintEnergyBin(std::ostream &os) const {

  std::ostringstream binheader;
  binheader.setf(std::ios::fixed);
  binheader.width(1);
  binheader.precision(0);
  binheader << fkeepbin;
  binheader.width(10);
  binheader.precision(3);
  binheader << fe_min;
  binheader.width(10);
  binheader << fe_mean;
  binheader.width(10);
  binheader << fe_max;
  binheader.precision(2);
  binheader.width(10);
  binheader << fon;
  binheader.width(10);
  binheader << foff;
  binheader.width(10);
  binheader << falpha;
  binheader.width(10);
  binheader << (fon-falpha*foff);
  binheader.width(10);
  binheader << STARTUtils::LiMaSignificance(fon,foff,falpha);
  binheader.width(10);
  binheader << flivetime;
  binheader.width(10);
  binheader << facceff;
  binheader.width(10);
  binheader << fsth+falpha*foff_fitted;
  binheader.width(10);
  binheader << foff_fitted;
  binheader.width(10);
  binheader << fsth;

  os << binheader.str().c_str() << std::endl;

}

/**
 * \brief Set interpolator of integrated resolutionn in bin
 */
void START::EnergyBin::SetGSLInterpolatorForPartialIntegral(std::pair<std::vector<double>, std::vector<double> > partialintegral)
{
  fpartialintegral = partialintegral;
  
  if (fInterpolPartialIntegral!=0) delete fInterpolPartialIntegral;
  fInterpolPartialIntegral=0;
  
  if (fpartialintegral.first.size()>4 && fpartialintegral.second.size()>4) {
    fInterpolPartialIntegral = new ROOT::Math::Interpolator(fpartialintegral.first,fpartialintegral.second,fInterpolTypeForBin);
  }
  else {
    fInterpolPartialIntegral = 0;
  }
}

/**
 * \brief Return the integrated resolution for true energy
 *
 *\param E true energy
 */
double START::EnergyBin::GetInterpolatedPartialIntegral(double E) {
  if(fInterpolPartialIntegral==0) {
    std::cout << "You have to initialize Interpolators before to use them!!!! ==> EXIT" << std::endl;
    exit(EXIT_FAILURE);
  }
  return (fInterpolPartialIntegral->Eval(E));
}

/**
 * \brief Function to remove all the informations from a bin, keep only the energy range
 * (fpartialintegral & fInterpolPartialIntegral are not reset)
 * 
 */
void START::EnergyBin::ClearEnergyBinInfo() {

  fkeepbin=0;
  fon=0;
  foff=0;
  facceff=0;
  fsth=0;
  falpha=0;
  fon_fitted=0;
  foff_fitted=0;
  flivetime=0.; 
  fEvents.clear();

}


/**
 * \brief Function to add information from another EnergyBin.
 *
 * \warning Does not handle the emin emid and emax values...
 *
 * \warning I still don't like the treatment I made for alpha....
 *
 * \todo A better handling of te alpha averaging thing........ (see also in Band.C and ResidualsFactory.C). Actually, I dislike how I did it !
 */

void START::EnergyBin::AddInfoFromEBin(const EnergyBin &AnotherEBin) {
  double alphamean = 0.;

  if (foff==0. || AnotherEBin.GetOff()==0.) { // VIM : This value, or the expected value ????
    alphamean = (falpha*flivetime + AnotherEBin.GetAlpha()*AnotherEBin.GetLiveTime())/(flivetime+AnotherEBin.GetLiveTime());
  }
  else {
    alphamean = ( falpha*foff + AnotherEBin.GetAlpha()*AnotherEBin.GetOff() )/ (foff + AnotherEBin.GetOff());
  }

  fon = fon+AnotherEBin.GetOn(); 
  foff = foff+AnotherEBin.GetOff(); 
  fsth = fsth+AnotherEBin.GetSth(); 
  fon_fitted=fon_fitted+AnotherEBin.GetOnFitted();
  foff_fitted=foff_fitted+AnotherEBin.GetOffFitted();
  flivetime = flivetime+AnotherEBin.GetLiveTime();
  falpha = alphamean;
  facceff = facceff + AnotherEBin.GetAcceff();
  fEvents = AnotherEBin.GetEvents();
}
