// STL
#include <vector>

// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TMath.h>
#include <TRolke.h>

// START
#include "LightCurveFactory.hh"
#include "Hypothesis.hh"
#include "Band.hh"
#include "EnergyBin.hh"
#include "STARTUtils.hh"
#include "ComputeResults.hh"
#include "EnergyBin.hh"
#include "Event.hh"

// Utilities
#define DEBUG 0
#include "debugging.hh"

#define INFO std::cout << INFOCOLOR << "LightCurveFactory> " << RESETCOLOR
#define WARNING std::cout << WARNINGCOLOR << "LightCurveFactory> " << RESETCOLOR

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::LightCurveFactory)
#endif

/**
 * \brief Constructor for a vector of Hypothesis
 * \param HypothesisArray vector of hypothesis
 * \param verbose verbose level
 */
START::LightCurveFactory::LightCurveFactory(std::vector<Hypothesis*> HypothesisArray, Config &Configuration, bool verbose)
{
  fHypothesisArray.clear();
  for(std::vector<Hypothesis*>::iterator hypo=HypothesisArray.begin(); hypo!=HypothesisArray.end(); ++hypo) {
    fHypothesisArray.push_back(&(**hypo));
  }

  fConfig = new Config(Configuration);

  fverbose=verbose;
}

/**
 * \brief Constructor for one Hypothesis
 * \param Hypo Hypothesis
 * \param verbose verbose level
 */
START::LightCurveFactory::LightCurveFactory(Hypothesis &Hypo, Config &Configuration, bool verbose)
{

  fHypothesisArray.clear();
  fHypothesisArray.push_back(&Hypo);

  fConfig = new Config(Configuration);

  fverbose=verbose;
}

/**
 * \brief Destructor
 */
START::LightCurveFactory::~LightCurveFactory()
{
  if(fConfig!=0) delete fConfig;
  fConfig=0;
}

void START::LightCurveFactory::SetAndBuildLightCurveData(TimeCuttingType TimeCutting)
{

  fTimeCutting = TimeCutting;

  switch(fTimeCutting) {
  case RunByRun:
    fTimeCuttingName="RunByRun";
    break;
  case MinuteByMinute:
    fTimeCuttingName="MinuteByMinute";
    break;
  case HourByHour:
    fTimeCuttingName="HourByHour";
    break;
  case NightByNight:
    fTimeCuttingName="NightByNight";
    break;
  case DayByDay:
    fTimeCuttingName="DayByDay";
    break;
  case WeekByWeek:
    fTimeCuttingName="WeekByWeek";
    break;
  case MonthByMonth:
    fTimeCuttingName="MonthByMonth";
    break;
  case YearByYear:
    fTimeCuttingName="YearByYear";
    break;
  case GivenTimeInterval:
    fTimeCuttingName="GivenTimeInterval";
    fTimeCuttingName+="(";
    fTimeCuttingName+=fConfig->GetUserLightCurveTimeInterval();
    fTimeCuttingName+="s)";
    break;
  case UserTimeIntervals:
    fTimeCuttingName="UserTimeIntervals";
    break;
  case PeriodByPeriod:
    fTimeCuttingName="PeriodByPeriod";
    break;
  default:
    WARN_OUT << "Nothing to do here ==> EXIT!" << std::endl;
    exit(EXIT_FAILURE);
  }

  DEBUG_OUT << "Light curve users parameters :" << std::endl;
  DEBUG_OUT << "TimeCuttingType = " << fTimeCuttingName << std::endl;
  DEBUG_OUT << "TLCmin = " << fConfig->GetUserLightCurveTimeRange().first
	    << ", TLCmax = " << fConfig->GetUserLightCurveTimeRange().second << std::endl;
  DEBUG_OUT << "TimeInterval = " << fConfig->GetUserLightCurveTimeInterval() << std::endl;
  DEBUG_OUT << "Integrated flux energy range = [" << fConfig->GetUserLightCurveIntegratedFluxEnergyRange().first
	    << ":" << fConfig->GetUserLightCurveIntegratedFluxEnergyRange().second << "]" << std::endl;

  for(std::vector<Hypothesis*>::iterator hypos=fHypothesisArray.begin(); hypos!=fHypothesisArray.end(); hypos++) {
    
    Hypothesis *hypo = &**hypos;

    if(hypo->GetConvergence()) {
      
      DEBUG_OUT << "LC for hypothesis " << hypo->GetName() << std::endl;
      
      hypo->SetLightCurveIntegratedFluxEnergyRange(fConfig->GetUserLightCurveIntegratedFluxEnergyRange().first,
						   fConfig->GetUserLightCurveIntegratedFluxEnergyRange().second);

      TimeBinVector TimeBinnedData;
      
      switch(fTimeCutting) {
      case RunByRun:
      case MinuteByMinute:
      case HourByHour:
      case NightByNight:
      case DayByDay:
      case WeekByWeek:
      case MonthByMonth:
      case YearByYear:
      case GivenTimeInterval:
      case UserTimeIntervals:
      case PeriodByPeriod:
	BuildGivenTimeIntervalLightCurve(*hypo,TimeBinnedData);
	break;
      default:
	WARN_OUT << "Nothing to do here ==> EXIT!" << std::endl;
	exit(EXIT_FAILURE);
      }
      INFO << "Data binned in time for " << hypo->GetName() << ": " << std::endl;
      TimeBinnedData.Print();
      hypo->AddTimeBinVector(fTimeCuttingName,TimeBinnedData);
    }
    else {
      WARNING << "Skip hypothesis " << hypo->GetName() << " because it didn't converge!" << std::endl;
      continue;
    }
    
    
  } // end loop on hypothesis
  
  INFO << "Building "<< fTimeCuttingName <<" LightCurve... ok" << std::endl;
}


/**
 * \brief Return time intervals
 */
std::vector<std::pair<double, double> > START::LightCurveFactory::GetAllTimeIntervals(const std::vector<Band> &BandArray) {

  double tmin(fConfig->GetUserLightCurveTimeRange().first);
  double tmax(fConfig->GetUserLightCurveTimeRange().second);
  double tbinsize(0.);

  switch(fTimeCutting) { // in seconds
  case RunByRun:
    break; // no bin size
  case MinuteByMinute:
    tbinsize=60.;
    break;
  case HourByHour:
    tbinsize=60.*60.;
    break;
  case NightByNight:
    tbinsize=60.*60.*12.;
    break;
  case DayByDay:
    tbinsize=60.*60.*24.;
    break;
  case WeekByWeek:
    tbinsize=60.*60.*24.*7.;
    break;
  case MonthByMonth:
    tbinsize=60.*60.*24.*30.;
    break;
  case YearByYear:
    tbinsize=60.*60.*24.*365.25;
    break;
  case GivenTimeInterval:
    tbinsize=fConfig->GetUserLightCurveTimeInterval();
    break;
  case UserTimeIntervals:
    // Nothing to do, intervals are given by user
    break;
  case PeriodByPeriod:
    // Nothing to do, intervals are given with a specified function
    break;
  default:
    WARN_OUT << "Nothing to do here ==> EXIT!" << std::endl;
    exit(EXIT_FAILURE);
  }

  tbinsize/=(3600.*24.); // convert seconds in day

  std::vector<std::pair<double, double> > TimeIntervals;
  std::vector<std::pair<double, double> > TimeIntervalsTmp;

  switch(fTimeCutting) { // in seconds
  case RunByRun:
    
    for(std::vector<Band>::const_iterator band=BandArray.begin(); band!=BandArray.end(); ++band) {
      if(band->GetKeepBand()==0) continue;
      if(band->GetRunStartTime()>=tmin && band->GetRunEndTime()<=tmax)
	TimeIntervals.push_back(std::make_pair(band->GetRunStartTime(),band->GetRunEndTime()));
    }
    break;
  case MinuteByMinute:
  case HourByHour:
  case NightByNight:
  case DayByDay:
  case WeekByWeek:
  case MonthByMonth:
  case YearByYear:
  case GivenTimeInterval:
    for(double time(tmin); time<tmax*1.001; time+=tbinsize) {

      bool keepinterval(false);
      for(std::vector<Band>::const_iterator band=BandArray.begin(); band!=BandArray.end(); ++band) {
	if(band->GetKeepBand()==0) continue;
	// we add intervals only if there is intersection with at least one run
	if( (time<=band->GetRunStartTime() && (time+tbinsize)>=band->GetRunEndTime()) || // run included in interval
	    (time<=band->GetRunStartTime() && (time+tbinsize)<=band->GetRunEndTime()) || // tbin_min lower than runstart and tbin_max lower than runend
	    (time>=band->GetRunStartTime() && (time+tbinsize)>=band->GetRunEndTime()) || // tbin_min greater than runstart and tbin_max greater than runend
	    (time>=band->GetRunStartTime() && (time+tbinsize)<=band->GetRunEndTime()))  { // interval included in run
	  keepinterval=true;
	  break;
	}
      }

      if(keepinterval) TimeIntervals.push_back(std::make_pair(time,time+tbinsize));

    }

    break;
  case UserTimeIntervals:
    TimeIntervals = fConfig->GetUserLightCurveUserTimeIntervals();
    break;
  case PeriodByPeriod:
    TimeIntervalsTmp = STARTUtils::GetNewMoonMJDPeriods();

    for(std::vector<std::pair<double,double> >::const_iterator Period=TimeIntervalsTmp.begin();
	Period!=TimeIntervalsTmp.end(); ++Period) {

      bool keepinterval(false);

      double timestart(Period->first), timeend(Period->second);
      
      if(timestart<tmin) continue;
      if(timeend>tmax) continue;
      
      for(std::vector<Band>::const_iterator band=BandArray.begin(); band!=BandArray.end(); ++band) {
	if(band->GetKeepBand()==0) continue;

	// we add intervals only if there is intersection with at least one run
	if( (timestart<=band->GetRunStartTime() && timeend>=band->GetRunEndTime()) || // run included in interval
	    (timestart<=band->GetRunStartTime() && timeend<=band->GetRunEndTime()) || // tbin_min lower than runstart and tbin_max lower than runend
	    (timestart>=band->GetRunStartTime() && timeend>=band->GetRunEndTime()) || // tbin_min greater than runstart and tbin_max greater than runend
	    (timestart>=band->GetRunStartTime() && timeend<=band->GetRunEndTime()))  { // interval included in run
	  keepinterval=true;
	  break;
	}
      }
      
      if(keepinterval) TimeIntervals.push_back(std::make_pair(Period->first,Period->second));

    }


    break;
  default:
    WARN_OUT << "Nothing to do here ==> EXIT!" << std::endl;
    exit(EXIT_FAILURE);
  }

  DEBUG_OUT << "tmin=" << tmin << " tmax=" << tmax << " tbinsize=" << tbinsize << std::endl;

  if(TimeIntervals.size()>0) {
  DEBUG_OUT << "Intervals size(" << TimeIntervals.size() << "), first([" << TimeIntervals.front().first << ":" 
	    << TimeIntervals.front().second << "])"
	    << " , last([" << TimeIntervals.back().first << ":" << TimeIntervals.back().second << "])" << std::endl;
  }
  else {
    WARN_OUT << "Null intervals size ==>EXIT!" << std::endl;
    //exit(EXIT_FAILURE);
  }

  return TimeIntervals;
}


/**
 * \brief Build light curve with a time interval defined by the user
 */
void START::LightCurveFactory::BuildGivenTimeIntervalLightCurve(Hypothesis &hypo, TimeBinVector &TDataContainer)
{

  std::vector<Band> BandArray;
  BandArray = hypo.GetBandArray();

  std::vector<std::pair<double, double> > TimeIntervals = GetAllTimeIntervals(BandArray);

  for(std::vector<std::pair<double, double> >::const_iterator intervals=TimeIntervals.begin(); 
      intervals!=TimeIntervals.end(); ++intervals) {
    
    DEBUG_OUT_L(2) << "TimeInterval [" << intervals->first << ":" << intervals->second << "]" << std::endl;
    
    // declaration of time bin attributs
    double TBinOn(0.), TBinOff(0.), TBinExpOn(0.), TBinExpOff(0.), TBinLivetime(0.);
    double TBinExcess(0.), TBinExpExcess(0.), TBinIntegratedFlux(0.), TBinIntegratedFluxError(0.);    
    double TBinFactorErrorIntegaredFlux(0.);
    double TBinIntegratedFluxError1SigmaPlus(0.), TBinIntegratedFluxError1SigmaMinus(0.);
    double TBinIntegratedFluxError3SigmaPlus(0.), TBinIntegratedFluxError3SigmaMinus(0.);
    double alpha_mean(0.);
    double alpha_mean_backup(0.);

    for(std::vector<Band>::iterator band=BandArray.begin(); band!=BandArray.end(); ++band) {
      
      if(band->GetKeepBand()==0) continue; // we skipp not selected bands
      
      DEBUG_OUT_L(2) << "BandInterval [" << band->GetRunStartTime() << ":" << band->GetRunEndTime() << "]" << std::endl;

      // we skipp the band if they are not in the interval
      if( (intervals->first<band->GetRunStartTime() && intervals->second<band->GetRunStartTime()) || // interval not in run
	  (intervals->first>band->GetRunEndTime() && intervals->first>band->GetRunEndTime()) ) { // interval not in run
	
	DEBUG_OUT_L(2) << "Band " << band->GetNbRun() << " [" << band->GetRunStartTime() 
		       << ":" << band->GetRunEndTime() <<"] is not in MJD range (timeintervals=" 
		       << intervals->first << "," << intervals->second << ")" << std::endl;
	continue;
      }
      DEBUG_OUT_L(2) << "Band " << band->GetNbRun() << " [" << band->GetRunStartTime() 
		     << ":" << band->GetRunEndTime() <<"] is in MJD range (timeintervals=" 
		     << intervals->first << "," << intervals->second << ")!" << std::endl;
      
      // We build a temporary vector of energy bins which will only contains
      // ON and OFF numbers, and IRF previously set for integrated flux computation
      std::vector<EnergyBin> TmpEnergyBinVector = band->ebin;
      for(std::vector<EnergyBin>::iterator bin=TmpEnergyBinVector.begin(); bin!=TmpEnergyBinVector.end(); ++bin) {
				bin->SetOn(0);
				bin->SetOff(0);
				bin->SetOnFitted(0);
				bin->SetOffFitted(0);
      }

      // Filling ON and OFF for vector of EnergyBin
      int acceptedon(0);
      // ON events loop
      
      for(std::vector<EnergyBin>::iterator bin=TmpEnergyBinVector.begin(); bin!=TmpEnergyBinVector.end(); ++bin) {

	if(bin->GetKeepBin()==0) continue;

	std::vector<START::Event> EventsList = bin->GetEvents();

	for(std::vector<START::Event>::const_iterator OnEvt=EventsList.begin(); OnEvt!=EventsList.end(); ++OnEvt) {

	  if(OnEvt->GetIsOff()) continue; // we skipp off events

	  if(OnEvt->GetTimeMJD()>=intervals->first && OnEvt->GetTimeMJD()<=intervals->second) {
			/*
	    double ebinmin(bin->GetEmin());
	    double ebinmax(bin->GetEmax());

	    // we skipp energy bins not in the range of the energy range light curve
	    if(fConfig->GetUserLightCurveIntegratedFluxEnergyRange().first>ebinmax || 
	       fConfig->GetUserLightCurveIntegratedFluxEnergyRange().second<ebinmin)
	      continue;
	  
	    // If Ebinmin>EminLC then Ebinmin=EminLC
	    // If Ebinmax>EmaxLC then Ebinmax=EmaxLC

	    if(fConfig->GetUserLightCurveIntegratedFluxEnergyRange().first>ebinmin && 
	       fConfig->GetUserLightCurveIntegratedFluxEnergyRange().first<ebinmax)
	      ebinmin = fConfig->GetUserLightCurveIntegratedFluxEnergyRange().first;
	    if(fConfig->GetUserLightCurveIntegratedFluxEnergyRange().second>ebinmin && 
	       fConfig->GetUserLightCurveIntegratedFluxEnergyRange().second<ebinmax)
	      ebinmax =  fConfig->GetUserLightCurveIntegratedFluxEnergyRange().second;
	  
	    if(OnEvt->GetEnergy()>=ebinmin && OnEvt->GetEnergy()<ebinmax) {
	      bin->SetOn(bin->GetOn()+1);
	      acceptedon++;
	    }
			*/
	      bin->SetOn(bin->GetOn()+1);
	      acceptedon++;


	  } 

	} // loop on Event

      } // loop on TimeBin
      
      // OFF events loop on
      int acceptedoff(0);

      for(std::vector<EnergyBin>::iterator bin=TmpEnergyBinVector.begin(); bin!=TmpEnergyBinVector.end(); ++bin) {

	if(bin->GetKeepBin()==0) continue;

	std::vector<START::Event> EventsList = bin->GetEvents();

	for(std::vector<START::Event>::const_iterator OffEvt=EventsList.begin(); OffEvt!=EventsList.end(); ++OffEvt) {

	  if(OffEvt->GetIsOn()) continue; // we skipp On events

	  if(OffEvt->GetTimeMJD()>=intervals->first && OffEvt->GetTimeMJD()<=intervals->second) {

			/*
	    double ebinmin(bin->GetEmin());
	    double ebinmax(bin->GetEmax());

	    // we skipp energy bins not in the range of the energy range light curve
	    if(fConfig->GetUserLightCurveIntegratedFluxEnergyRange().first>ebinmax || 
	       fConfig->GetUserLightCurveIntegratedFluxEnergyRange().second<ebinmin)
	      continue;

	    // If Ebinmin>EminLC then Ebinmin=EminLC
	    // If Ebinmax>EmaxLC then Ebinmax=EmaxLC
	    if(fConfig->GetUserLightCurveIntegratedFluxEnergyRange().first>ebinmin && 
	       fConfig->GetUserLightCurveIntegratedFluxEnergyRange().first<ebinmax)
	      ebinmin = fConfig->GetUserLightCurveIntegratedFluxEnergyRange().first;
	    if(fConfig->GetUserLightCurveIntegratedFluxEnergyRange().second>ebinmin && 
	       fConfig->GetUserLightCurveIntegratedFluxEnergyRange().second<ebinmax)
	      ebinmax =  fConfig->GetUserLightCurveIntegratedFluxEnergyRange().second;

	    if(OffEvt->GetEnergy()>=ebinmin && OffEvt->GetEnergy()<ebinmax) {
	      bin->SetOff(bin->GetOff()+1);
	      acceptedoff++;
	    }
			*/

	      bin->SetOff(bin->GetOff()+1);
	      acceptedoff++;


	  }
	}
      }
      
      DEBUG_OUT_L(2) << " BandOn(" << band->GetNOnTot(true) <<"), BandOff(" << band->GetNOffTot(true) << "),"
		     << " acceptedOn(" << acceptedon << ") acceptedOff(" << acceptedoff << ")" << std::endl;

      // Livetime of interval
      double LivetimeToAdd(0.);
      if(intervals->first>=band->GetRunStartTime() && intervals->second<=band->GetRunEndTime()) {
	LivetimeToAdd=(intervals->second-intervals->first); // interval included in run
	DEBUG_OUT_L(2) << "Case 1" << std::endl;
      }
      else if(intervals->first>=band->GetRunStartTime() && intervals->second>=band->GetRunEndTime()) {
	LivetimeToAdd=(band->GetRunEndTime()-intervals->first); // above and in run
	DEBUG_OUT_L(2) << "Case 2" << std::endl;
      }
      else if(intervals->first<=band->GetRunStartTime() && intervals->second<=band->GetRunEndTime()) {
	LivetimeToAdd=(intervals->second-band->GetRunStartTime()); // below and in run
	DEBUG_OUT_L(2) << "Case 3" << std::endl;
      }
      else if(intervals->first<=band->GetRunStartTime() && intervals->second>=band->GetRunEndTime()) {
	LivetimeToAdd=(band->GetRunEndTime()-band->GetRunStartTime());
	DEBUG_OUT_L(2) << "Case 4" << std::endl; // run included in interval
      }
      else { // not in run, shouldn't be here
	WARN_OUT << "We should'nt be here ==> EXIT!" << std::endl;
	std::cout << "Band : " << band->GetRunStartTime() << " " << band->GetRunEndTime() << std::endl;
	std::cout << "Intervals : " << intervals->first << " " << intervals->second << std::endl;
	exit(EXIT_FAILURE);
      }
      LivetimeToAdd*=24.; // MJD to hours
      LivetimeToAdd*=band->GetLiveTimeFraction();
      DEBUG_OUT_L(2) << "BandLiveTime(" << band->GetLiveTime() << "), LivetimeToAdd(" << LivetimeToAdd << ")" << std::endl;

      // Compute integrated flux
      double SumOn(0.), SumOff(0.), SumExpExcess(0.);
      ComputeResults *Calculate = new ComputeResults(BandArray,hypo);
      for(std::vector<EnergyBin>::iterator bin=TmpEnergyBinVector.begin(); bin!=TmpEnergyBinVector.end(); ++bin) {
	if(bin->GetKeepBin()==0) continue; // we skip bins below threshold
	// JLK : appliquer interval en energie aux energies mesurees.
	double Sth = Calculate->FunctionExpectedExcess((&*band),(&*bin),
																								 hypo.GetFittedParameters());
	Sth/=band->GetLiveTime(); // Interval Livetime will be used instead of band's
	Sth*=LivetimeToAdd;

	SumOn+=bin->GetOn();
	SumOff+=bin->GetOff();
	SumExpExcess+=Sth;


      } // end EnergyBins loop 
      
      alpha_mean+=band->GetAlphaRun()*SumOff;
      alpha_mean_backup+=band->GetAlphaRun()*LivetimeToAdd;

      TBinOn+=SumOn;
      TBinOff+=SumOff;
      TBinExpExcess+=SumExpExcess;
      TBinExcess+=STARTUtils::GetExcess(SumOn,SumOff,band->GetAlphaRun());
      TBinExpOn=STARTUtils::GetExpectedOn(SumOn,SumOff,band->GetAlphaRun(),SumExpExcess);
      TBinExpOff=STARTUtils::GetExpectedOff(SumOn,SumOff,band->GetAlphaRun(),SumExpExcess);
      TBinLivetime+=LivetimeToAdd;
      TBinFactorErrorIntegaredFlux+=(SumOn+TMath::Power(band->GetAlphaRun(),2.)*SumOff);

      DEBUG_OUT << "SumOn=" << SumOn << " SumOff=" << SumOff << " Alpha=" << band->GetAlphaRun() 
		<< " Excess=" << STARTUtils::GetExcess(SumOn,SumOff,band->GetAlphaRun()) 
		<< " ExpExcess=" << SumExpExcess << std::endl;
      DEBUG_OUT << "TBinOn=" << TBinOn << " TBinOff=" << TBinOff 
		<< " TBinExcess=" << TBinExcess 
		<< " TBinExpExcess=" << TBinExpExcess << std::endl;

      // Clean
      delete Calculate; Calculate=0;
      
    } // end Bands loop

    if(TBinLivetime==0) continue;

    double thflux = hypo.GetFluxIntegralFitParams(fConfig->GetUserLightCurveIntegratedFluxEnergyRange().first,
						  fConfig->GetUserLightCurveIntegratedFluxEnergyRange().second);

    TBinIntegratedFlux=TBinExcess;
    TBinIntegratedFlux/=TBinExpExcess;
    TBinIntegratedFlux*=thflux;

    TBinIntegratedFluxError=TMath::Sqrt(TBinFactorErrorIntegaredFlux);
    TBinIntegratedFluxError/=TBinExpExcess;
    TBinIntegratedFluxError*=thflux;

    if(TBinOff>0) alpha_mean/=TBinOff;
    if(TBinLivetime>0) alpha_mean_backup/=TBinLivetime;
    if(alpha_mean==0) alpha_mean=alpha_mean_backup;

    TRolke ConfInter1Sigma;
    ConfInter1Sigma.SetCL(0.6827); // 1 sigma    
    // VIM : Set Bounding method to true (Rolke 2005) --> yield a bit more overcoverage when 
    // significant negative excess (but give results more according to our feeling)
    ConfInter1Sigma.SetBounding(true);
    double minns68;
    double maxns68;
    ConfInter1Sigma.SetPoissonBkgKnownEff((int)TBinOn,(int)TBinOff,1.0/alpha_mean,1.0);    
    ConfInter1Sigma.GetLimits(minns68,maxns68);
	  
    TRolke ConfInter3Sigma;
    ConfInter3Sigma.SetCL(0.9973); // 3 sigma
    // VIM : Set Bounding method to true (Rolke 2005) --> yield a bit more overcoverage when 
    // significant negative excess (but give results more according to our feeling)
    ConfInter3Sigma.SetBounding(true);
    double minns99;
    double maxns99;
    ConfInter3Sigma.SetPoissonBkgKnownEff((int)TBinOn,(int)TBinOff,1.0/alpha_mean,1.0);    
    ConfInter3Sigma.GetLimits(minns99,maxns99);   

    TBinIntegratedFluxError1SigmaPlus = maxns68*thflux/TBinExpExcess;
    TBinIntegratedFluxError1SigmaMinus = minns68*thflux/TBinExpExcess;

    TBinIntegratedFluxError3SigmaPlus = maxns99*thflux/TBinExpExcess;
    TBinIntegratedFluxError3SigmaMinus = minns99*thflux/TBinExpExcess;

    TDataContainer.tbin.push_back(TimeBin(intervals->first,intervals->second));
    TDataContainer.tbin.back().SetOn(TBinOn);
    TDataContainer.tbin.back().SetOff(TBinOff);
    TDataContainer.tbin.back().SetAlpha(alpha_mean);
    TDataContainer.tbin.back().SetExcess(TBinExcess);
    TDataContainer.tbin.back().SetExpectedOn(TBinExpOn);
    TDataContainer.tbin.back().SetExpectedOff(TBinExpOff);
    TDataContainer.tbin.back().SetExpectedExcess(TBinExpExcess);
    TDataContainer.tbin.back().SetIntegratedFlux(TBinIntegratedFlux);
    TDataContainer.tbin.back().SetIntegratedFluxError(TBinIntegratedFluxError);
    TDataContainer.tbin.back().SetLiveTime(TBinLivetime);
    TDataContainer.tbin.back().SetIntegratedFluxError1SigmaMinus(TBinIntegratedFluxError1SigmaMinus);
    TDataContainer.tbin.back().SetIntegratedFluxError1SigmaPlus(TBinIntegratedFluxError1SigmaPlus);
    TDataContainer.tbin.back().SetIntegratedFluxError3SigmaMinus(TBinIntegratedFluxError3SigmaMinus);
    TDataContainer.tbin.back().SetIntegratedFluxError3SigmaPlus(TBinIntegratedFluxError3SigmaPlus);
    if(TBinIntegratedFlux<0 || TBinIntegratedFluxError1SigmaMinus==0)
      TDataContainer.tbin.back().SetIsUpperLimit(true);
    else
      TDataContainer.tbin.back().SetIsUpperLimit(false);

    if(DEBUG) {
      TDataContainer.Print();
    }

  } // end time intervals loop
  
}
