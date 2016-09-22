// STL
#include <iostream>
#include <vector>

// ROOT
#include <TTree.h>
#include <TMath.h>
#include <TFile.h>
#include <TBranch.h>
#include <TROOT.h>

// START
#include "BandsFactory.hh"
#include "Config.hh"
#include "EnergyBin.hh"
#include "Band.hh"
#include "MonteCarlo.hh"
#include "STARTUtils.hh"
#include "Event.hh"
#include "ComputeResults.hh"
#include "HandleResolArea.hh"
// Utilities
#define DEBUG 0
#include "debugging.hh"

#define INFO std::cout << INFOCOLOR << "BandsFactory> " << RESETCOLOR
#define WARNING std::cout << WARNINGCOLOR << "BandsFactory> " << RESETCOLOR

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::BandsFactory)
#endif

/**
 * \brief Constructor
 *
 */
START::BandsFactory::BandsFactory(const Config &config,
				  START::STARTUtils::IrfOpt Type,
				  bool verbose)
{

  gROOT->cd();
	
  fConfig = new Config(config);
  fIrfOpt = Type;

  fRunList = fConfig->GetRunList();
  fRunTelType = fConfig->GetAllowedTelList();
  fInputFileList = fConfig->GetInputFileList();

  fTimeToAdd=40587; //SP: conversion entre temps SASH et temps MJD, j'aimerais qu'un oeil externe vérifie si ok

  fBinNumber = fConfig->GetUserERangeBinsNumber();

  fEmax = TMath::Log10(fConfig->GetUserERangeMax());
  fEmin = TMath::Log10(fConfig->GetUserERangeMin());

  fBinWidth = (fEmax - fEmin)/fBinNumber;

  fNbBands = -1; // will be used to access element of bandarrays
  // not sexy at all but it works
}

/**
 * destructor
 */
START::BandsFactory::~BandsFactory()
{
  if(fConfig!=0) delete fConfig;
  fConfig=0;
}


/**
 * \brief Read informations of the runs contained in rootfiles
 * and stock them in the bojects bands
 *
 * \param vector of bands
 */
int START::BandsFactory::MakeBands(std::vector<Band> &BandArray)
{
  BandArray.clear(); // JLK un minimum quand meme...
  
  gROOT->Reset();

  for(unsigned int irun(0); irun<fRunList.size(); irun++) {

    /*******************************************************************************/
    /************************* Open File of the run ********************************/
    /*******************************************************************************/

    TString treepath = "";
    if (fConfig->GetUserUseInputFileList()) {
      treepath+=fInputFileList[irun];
    }
    else {
      treepath+=fConfig->GetUserRootFilesFolderAddress();
      treepath+="/";

      treepath+=STARTUtils::GetRunROOTFileName(fConfig->GetUserAnalysisConfig(),fRunList[irun],fRunTelType[irun]);
      /*
	treepath+="/run";
	treepath += fRunList[irun];
	treepath += "_";
	treepath += fConfig->GetUserAnalysisConfig();
	// JLK : need to be improved by using MonteCarlo class
	// If I give a telcode to MC function it gives me the right string to add
	if(fRunTelType[irun]==14 || fRunTelType[irun]==22 || fRunTelType[irun]==26 || fRunTelType[irun]==28) {
	treepath+="_3Tel";
	}
	treepath+=".root";
      */
    }

    TFile *treefile = new TFile(treepath,"READ");
    
    if(!treefile->IsOpen()) {
      WARNING << "Failed to open file :" << treepath << std::endl;
      continue;
    }

    /*******************************************************************************/
    /******** Get Global Parameter of the run, i.e., muon eficiency, alpha *********/
    /*******************************************************************************/

    TTree *RunInfo_TTree = (TTree *)treefile->Get(STARTUtils::GetRunInfoTreeName());
    if(!RunInfo_TTree){
      WARNING << "Run " << fRunList[irun] << " RunInfo_TTree does not exist => SKIP RUN" << std::endl;
      RunInfo_TTree=0;
      continue;
    }

    RunInfo_TTree->SetDirectory(treefile);
  
    TBranch *alpha_branch = RunInfo_TTree->GetBranch("fAlpha");
    TBranch *muonefi_branch = RunInfo_TTree->GetBranch("fMuonEff");
    TBranch *duration_branch = RunInfo_TTree->GetBranch("fDuration");
    TBranch *livetime_branch = RunInfo_TTree->GetBranch("fLiveTime");

    Double_t AlphaRUN;
    Double_t MuonEffRUN;
    Double_t DurationRUN;
    Double_t LiveTimeRUN;

    alpha_branch->SetAddress(&AlphaRUN);
    muonefi_branch->SetAddress(&MuonEffRUN);
    duration_branch->SetAddress(&DurationRUN);
    livetime_branch->SetAddress(&LiveTimeRUN);

    alpha_branch->GetEntry();
    muonefi_branch->GetEntry();
    duration_branch->GetEntry();
    livetime_branch->GetEntry();

    DEBUG_OUT << fRunList[irun] 
	      << " : AlphaRUN = " << AlphaRUN 
	      << " MuonEffRUN=" << MuonEffRUN 
	      << " DurationRUN=" << DurationRUN 
	      << " LiveTimeRUN=" << LiveTimeRUN << std::endl;

    /* 
       SP: below is a way to recover the LiveTimeFraction which comes from DB RunQuality.
       Have noticed that DurationRUN (which comes from DB Run_end) is generally 2 or 3% larger than
       the time between the first and last event in the AllSelected TTree, and this seems to be also true
       if DurationRUN is compared to event times in the RawData file.
       For this reason, it is better to take a ScanDurationTime which comes from Events (from AllEvents TTree to have enough statistics)
       and re-determine the LiveTime using this ScanDurationTime and the LiveTimeFraction.
    */
    double LiveTimeFraction = LiveTimeRUN/DurationRUN;

    //ADA Hack for Mono
    if(fIrfOpt==START::STARTUtils::HapFr_Hess_Mono)
      LiveTimeFraction = 1.07;

    if(LiveTimeFraction<0.8) {
      WARNING << "Run "<<fRunList[irun]
	      <<" Very Big DEAD TIME " << LiveTimeFraction <<" (check if DB is correct)"<<std::endl;
    }

    if(MuonEffRUN<0.5) {
      WARNING << "Run "<<fRunList[irun]
	      <<" Very Low MUON EFF. " << MuonEffRUN <<std::endl;
    }

    /*******************************************************************************/
    /************************ Prepare to read BgMakerOff TTree *********************/
    /*******************************************************************************/

    TTree *EventsTTree_BgMakerOff = (TTree*) treefile->Get(STARTUtils::GetOffEventsTreeName());
    if(!EventsTTree_BgMakerOff){
      WARNING << "Run " << fRunList[irun] << " EventsTTree_BgMakerOff does not exist => SKIP RUN" << std::endl;
      EventsTTree_BgMakerOff=0;
      continue;
    }
    EventsTTree_BgMakerOff->SetDirectory(treefile);

    Int_t nbentriesOff = EventsTTree_BgMakerOff->GetEntries();
    TBranch *offtime_branch_secs = EventsTTree_BgMakerOff->GetBranch("fEvtTimeSec");
    TBranch *offtime_branch_nano = EventsTTree_BgMakerOff->GetBranch("fEvtTimeNano");
    TBranch *offzen_branch = EventsTTree_BgMakerOff->GetBranch("fEvtZen");
    TBranch *Energy_branch_OFF = EventsTTree_BgMakerOff->GetBranch("fEvtEnergy");
    TBranch *EvtOffset_OFF = EventsTTree_BgMakerOff->GetBranch("fEvtOffset");

    Int_t OffEvtTimeSec;
    Int_t OffEvtTimeNano;
    Double_t OffZen;
    Double_t OffEvtEnergy;
    Double_t OffEvtOffset;

    offtime_branch_secs->SetAddress(&OffEvtTimeSec);
    offtime_branch_nano->SetAddress(&OffEvtTimeNano);
    offzen_branch->SetAddress(&OffZen);
    Energy_branch_OFF->SetAddress(&OffEvtEnergy);
    EvtOffset_OFF->SetAddress(&OffEvtOffset);


    
    // JLK : SP add :

    /*****************************************************************************************/
    /****************** Determine useful times from events individual times ******************/
    /*****************************************************************************************/

    Double_t TimeFirstEvt=0.;
    Double_t TimeLastEvt=0.;
    Double_t ScanDurationTime=0.;
    Double_t MiddleRunTime=0.;

    if(!fConfig->GetUserIsMc()) {
      TTree *EventsTTree_AllSelected = 0;
      //For Mono : lots of events so use the On Events TTree
      if(fIrfOpt==START::STARTUtils::HapFr_Hess_Mono)
	EventsTTree_AllSelected = (TTree *)treefile->Get(STARTUtils::GetOnEventsTreeName());
      else
	EventsTTree_AllSelected = (TTree *)treefile->Get(STARTUtils::GetAllEventsTreeName());
      if(!EventsTTree_AllSelected){
	WARNING << "Run " << fRunList[irun] << " EventsTTree_AllSelected does not exist => SKIP RUN" << std::endl;
	EventsTTree_AllSelected=0;
	continue;
      }
      EventsTTree_AllSelected->SetDirectory(treefile);
		  
      Int_t nbentriesAll = EventsTTree_AllSelected->GetEntries();
      
      TBranch *alltime_branch_secs = EventsTTree_AllSelected->GetBranch("fEvtTimeSec");
      TBranch *alltime_branch_nano = EventsTTree_AllSelected->GetBranch("fEvtTimeNano");
      TBranch *allzen_branch = EventsTTree_AllSelected->GetBranch("fEvtZen");
      
      Int_t AllEvtTimeSec;
      Int_t AllEvtTimeNano;
      Double_t AllEvtZen;
      
      alltime_branch_secs->SetAddress(&AllEvtTimeSec);
      alltime_branch_nano->SetAddress(&AllEvtTimeNano);
      allzen_branch->SetAddress(&AllEvtZen);
 
      // Get first event time
      alltime_branch_secs->GetEntry(0);
      alltime_branch_nano->GetEntry(0);
      TimeFirstEvt = (Double_t)AllEvtTimeSec+ ((Double_t)AllEvtTimeNano)/1e9;
      
      // Get last event time 
      alltime_branch_secs->GetEntry(nbentriesAll-1);
      alltime_branch_nano->GetEntry(nbentriesAll-1);
      TimeLastEvt = (Double_t)AllEvtTimeSec+ ((Double_t)AllEvtTimeNano)/1e9;
      
      ScanDurationTime = TimeLastEvt - TimeFirstEvt;
      
      // Check if events are correctly ordered
      Double_t prev=0;
      int nbadevents(0);
      for(int i=0;i<nbentriesAll;i++) {
	alltime_branch_secs->GetEntry(i);
	alltime_branch_nano->GetEntry(i);
	allzen_branch->GetEntry(i);
	Double_t TimeEvt = (Double_t)AllEvtTimeSec+ ((Double_t)AllEvtTimeNano)/1.e9;  
	if (i!=1 && TimeEvt<prev) {
	  nbadevents++;
	  WARNING << "Time order failure : "<<std::endl;
	  INFO << "Event (in AllSelected) # "<<i<<" has time "<<TimeEvt<<" which is smaller than previous event :"<<prev<<std::endl;
	  INFO << "Should check if run is ok" << std::endl;
	  if(nbadevents > 10){
	    WARNING << "More than 10 time order failure, Exiting ...  "<<std::endl;
	    exit(EXIT_FAILURE);
	  }
	}
	prev=TimeEvt;
      }
      
      
      // Determine time of the middle of the run 
      MiddleRunTime = 0.5*(TimeFirstEvt+TimeLastEvt);
      
      if (fverbose) {
	INFO << "Scan of AllSelected" << std::endl;
	std::cout << "- Events correctly ordered in time" << std::endl;
	std::cout << "- Number of events: " << nbentriesAll << std::endl;
	std::cout << "- First event time: " <<  TimeFirstEvt << std::endl;
	std::cout << "- Last event time: " <<  TimeLastEvt << std::endl;
	std::cout << "- Recovered observation time: " << ScanDurationTime << std::endl;
	std::cout << "- Duration time from DB (via RunInfo TTree): " << DurationRUN << std::endl;
	std::cout << "- Live time from DB  (via RunInfo TTree): " << LiveTimeRUN << std::endl;
      }



    }

    /*******************************************************************************/
    /************************* Prepare to read BgMakerOn TTree *********************/
    /*******************************************************************************/

    TTree *EventsTTree_BgMakerOn = (TTree *)treefile->Get(STARTUtils::GetOnEventsTreeName());
    if(!EventsTTree_BgMakerOn){
      WARNING << "Run " << fRunList[irun] << " EventsTTree_BgMakerOn does not exist => SKIP RUN" << std::endl;
      EventsTTree_BgMakerOn=0;
      continue;
    }
    EventsTTree_BgMakerOn->SetDirectory(treefile);

    Int_t nbentriesOn = EventsTTree_BgMakerOn->GetEntries();

    TBranch *ontime_branch_secs = EventsTTree_BgMakerOn->GetBranch("fEvtTimeSec");
    TBranch *ontime_branch_nano = EventsTTree_BgMakerOn->GetBranch("fEvtTimeNano");
    TBranch *onzen_branch = EventsTTree_BgMakerOn->GetBranch("fEvtZen");
    TBranch *Energy_branch_ON = EventsTTree_BgMakerOn->GetBranch("fEvtEnergy");
    TBranch *EvtOffset_ON  = EventsTTree_BgMakerOn->GetBranch("fEvtOffset");

    Int_t OnEvtTimeSec;
    Int_t OnEvtTimeNano;
    Double_t OnZen;
    Double_t OnEvtEnergy;
    Double_t OnEvtOffset;

    ontime_branch_secs->SetAddress(&OnEvtTimeSec);
    ontime_branch_nano->SetAddress(&OnEvtTimeNano);
    onzen_branch->SetAddress(&OnZen);
    Energy_branch_ON->SetAddress(&OnEvtEnergy);    
    EvtOffset_ON->SetAddress(&OnEvtOffset);
     
    //ADA 06/08/2016    : moved the following block here so as nbentriesOn is defined
    /*******************************************************************************/
    /**************************** Skip runs if problems ****************************/
    /*******************************************************************************/
    /* For hard config analysis, it happens than Non=0 */

    if(!fConfig->GetUserIsMc()) {
      //ADA 08/06/2016
      // if ((int)nbentriesOff==0 ||
      // 	  (int)Energy_branch_OFF->GetEntries()==0 ||
      // 	  (int)EvtOffset_OFF->GetEntries()==0 || 
      // 	  (int)offzen_branch->GetEntries()==0) {
      if ((int)nbentriesOff==0 && (int)nbentriesOn==0){
	WARNING << "Run " << fRunList[irun] << " Problem with the number of events => SKIP RUN" << std::endl;
	continue;
      }
    }


    /*******************************************************************************/
    /****************** Determine Zbands strategy & zbands edges  ******************/
    /*******************************************************************************/
    double ZenOnMin=90.; /*initialisation*/
    double ZenOnMax=0.;
    double ZenOffMin=90.;
    double ZenOffMax=0.;
    bool UseTwoZenBands=false;
    int RunNumberOfZBands=1;

    if(!fConfig->GetUserIsMc()) {
      for(int i=0;i<nbentriesOn;i++) {
	onzen_branch->GetEntry(i);	  	  
	if (OnZen<ZenOnMin)ZenOnMin=OnZen;
	if (OnZen>ZenOnMax)ZenOnMax=OnZen;
      }
      for(int i=0;i<nbentriesOff;i++) {
	offzen_branch->GetEntry(i);	  
	if (OffZen<ZenOffMin)ZenOffMin=OffZen;
	if (OffZen>ZenOffMax)ZenOffMax=OffZen;
      }
      
      double DeltaCosZenOn = TMath::Cos(ZenOnMin*TMath::Pi()/180.) - TMath::Cos(ZenOnMax*TMath::Pi()/180.);
      
      if (DeltaCosZenOn>(fConfig->GetUserMaxCosZenBinWidth())) {
	/* The following condition corresponds to what is supported below*/
	if(fConfig->GetUserMJDwindow().size()==0 || CheckIfEntireRunIncludedInSingleTimeWindow(TimeFirstEvt,TimeLastEvt)) {
	  UseTwoZenBands=true;
	  RunNumberOfZBands=2;
	  //if (fConfig->GetUserVerbose()) {
	  INFO << "Run " << fRunList[irun] << " (Zmin,Zmax)=(" << ZenOffMin << "," << ZenOffMax 
	       << ") : we'll use 2 zen bands because delta cos(zen) = "<< DeltaCosZenOn<<std::endl;
	  //}
	}
      }
      
    }


    /**************************************************************************/
    /****       Initialization of Bands objects & BinEnergy Objects        ****/
    /****       All attributs have value 0 except bin properties.          ****/
    /**************************************************************************/

    
    for(int zbin=0;zbin<RunNumberOfZBands;zbin++) {

      fNbBands++; // increase number of band
      BandArray.push_back(Band());
      BandArray.back().SetNbRun(fRunList[irun]);
      for(int enbin=0;enbin<fBinNumber;enbin++) {
	BandArray.back().ebin.push_back(EnergyBin());
	BandArray.back().ebin.back().SetEmin(TMath::Power(10,fEmin+enbin*fBinWidth));
	BandArray.back().ebin.back().SetEmax(TMath::Power(10,fEmin+(enbin+1)*fBinWidth));
	BandArray.back().ebin.back().SetEmid((BandArray.back().ebin.back().GetEmin()+(BandArray.back().ebin.back().GetEmax()))/2.); //SP: a quoi ça sert ? moyenne linéaire...
      }

    }

    /****************************************************/
    /*************** Determine livetimes ****************/
    /****************************************************/

    if(!fConfig->GetUserIsMc()) {
      /* If no time windows or if entire run included in a single time window */
      if(fConfig->GetUserMJDwindow().size()==0 || CheckIfEntireRunIncludedInSingleTimeWindow(TimeFirstEvt,TimeLastEvt)) {

	Double_t RunLiveTime = ScanDurationTime*LiveTimeFraction ;
	DEBUG_OUT << "JLK: " << " RunLiveTime=" << RunLiveTime << std::endl;
	if (UseTwoZenBands==false) {	
	  BandArray[fNbBands].SetLiveTime(RunLiveTime);
	  BandArray[fNbBands].SetLiveTimeFraction(LiveTimeFraction);
	  BandArray[fNbBands].InitEbinLiveTime(RunLiveTime);
	  // JLK
	  BandArray[fNbBands].SetRunStartTime(STARTUtils::GetMJDFromSashSeconds(TimeFirstEvt));
	  BandArray[fNbBands].SetRunEndTime(STARTUtils::GetMJDFromSashSeconds(TimeLastEvt));
	} else {
	  double starttime = STARTUtils::GetMJDFromSashSeconds(TimeFirstEvt);
	  double endtime = STARTUtils::GetMJDFromSashSeconds(TimeLastEvt);
	  BandArray[fNbBands-1].SetLiveTime(0.5*RunLiveTime);
	  BandArray[fNbBands-1].SetLiveTimeFraction(LiveTimeFraction);
	  BandArray[fNbBands-1].InitEbinLiveTime(0.5*RunLiveTime);
	  // JLK
	  BandArray[fNbBands-1].SetRunStartTime(starttime);
	  BandArray[fNbBands-1].SetRunEndTime(starttime+0.5*(endtime-starttime));

	  BandArray[fNbBands].SetLiveTime(0.5*RunLiveTime);
	  BandArray[fNbBands].SetLiveTimeFraction(LiveTimeFraction);
	  BandArray[fNbBands].InitEbinLiveTime(0.5*RunLiveTime);
	  // JLK
	  BandArray[fNbBands].SetRunStartTime(starttime+0.5*(endtime-starttime));
	  BandArray[fNbBands].SetRunEndTime(endtime);
	}

      }
      else {
	// Run is cut in one or more time windows 
	if (UseTwoZenBands==false) {
	  Double_t ObsTime=0;
	  for (int i=0;i<(int)fConfig->GetUserMJDwindow().size();i++) {
	    
	    // Get Max(begin_run,begin_window) 
	    Double_t max = std::max(STARTUtils::GetMJDFromSashSeconds(TimeFirstEvt),fConfig->GetUserMJDwindow()[i].first);
	    DEBUG_OUT << "JLK: TimeFirstEvt=" << TimeFirstEvt << " fConfig->GetUserMJDwindow()[i].first=" << fConfig->GetUserMJDwindow()[i].first 
		      <<  "StartutilsTimeLastEvt"<<std::endl;
	    // Get Min(end_run,end_window) 
	    Double_t min = std::min(STARTUtils::GetMJDFromSashSeconds(TimeLastEvt),fConfig->GetUserMJDwindow()[i].second);
	
	    if(max>min) {
	      ObsTime += (max-min);
	    } else {
	      INFO <<"Time window ["<<fConfig->GetUserMJDwindow()[i].first<<","<<fConfig->GetUserMJDwindow()[i].second<<"] is not in run "<<fRunList[irun]<<std::endl;
	    }
	  }

	  ObsTime/=3600.;// in hours
	  DEBUG_OUT << "JLK: " << " ObsTime=" << ObsTime << " LiveTimeFraction=" << LiveTimeFraction << std::endl;
	  BandArray[fNbBands].SetLiveTime(ObsTime*LiveTimeFraction);
	  BandArray[fNbBands].InitEbinLiveTime(ObsTime*LiveTimeFraction);
	  // JLK how to treat multiple windows for runstarttime and runendtime??
	  BandArray[fNbBands].SetRunStartTime(STARTUtils::GetMJDFromSashSeconds(TimeFirstEvt)); // JLK please check
	  BandArray[fNbBands].SetRunEndTime(STARTUtils::GetMJDFromSashSeconds(TimeLastEvt)); // JLK please check
	} else {
	  WARNING << "Cut of run into 2 bands not supported if there is time windows AND run not entirely included in a single time window" << std::endl;
	  INFO << "UseTwoZenBands should have been set to false before. Exit." <<std::endl;
	  exit(EXIT_FAILURE);
	}
      }

    } else {
      //SP: for MC (check if RunInfoTTree filled or LiveTimeRUN overwritten)
      BandArray[fNbBands].SetLiveTime(LiveTimeRUN);
      BandArray[fNbBands].InitEbinLiveTime(LiveTimeRUN);
      if(LiveTimeRUN==0.) {
	WARNING "MC with no livetime, assume a run with 1680s" <<std::endl;
	double livetimemc(1680./60.);
	BandArray[fNbBands].SetLiveTime(livetimemc);
	BandArray[fNbBands].InitEbinLiveTime(livetimemc);
      }
    }


    /***************************************************/
    /********  Fill arrays for BgMakerOn events  *******/
    /***************************************************/
    for(int i=0;i<nbentriesOn;i++) {
      
      ontime_branch_secs->GetEntry(i);
      ontime_branch_nano->GetEntry(i);
      onzen_branch->GetEntry(i);	  
      Energy_branch_ON->GetEntry(i);
      EvtOffset_ON->GetEntry(i);

      Double_t CurrentTimeOn = (Double_t)OnEvtTimeSec+ ((Double_t)OnEvtTimeNano)/1e9;

      //std::cout<<"xxx "<<CurrentTimeOn-1.07679e+09<<" "<<OnZen<<std::endl;

      if (!fConfig->GetUserIsMc() && (CurrentTimeOn<TimeFirstEvt || CurrentTimeOn>TimeLastEvt)) {
	WARNING << "Something wrong with time in run :" << fRunList[irun] 
		<< ". Fix it. (CurrentTimeOn:" << CurrentTimeOn 
		<< " TimeFirstEvt: " << TimeFirstEvt << " TimeLastEvt: " 
		<< TimeLastEvt << " --> TimeFirstEvt-CurrentTimeOn = " 
		<< TimeFirstEvt-CurrentTimeOn << " CurrentTimeOn-TimeLastEvt = " 
		<< CurrentTimeOn-TimeLastEvt << ")" <<std::endl;
	exit(EXIT_FAILURE);
      }

      Int_t Zbin=0; // Define the band to use
      if (UseTwoZenBands==true && CurrentTimeOn>MiddleRunTime) Zbin=1;	    
      
      // When needed, we check if the events are in the required time windows
      if (fConfig->GetUserMJDwindow().size()>0 && !CheckIfEventInTimeWindow(CurrentTimeOn)) 
	continue; //Skip if without time window
       
      if (OnEvtEnergy<fConfig->GetUserERangeMin())
	continue; //Skip events with energy lower than requested
      if (OnEvtEnergy>fConfig->GetUserERangeMax())
	continue; //Skip events with energy higher than requested
      
      Int_t Ebin = (Int_t)((TMath::Log10(OnEvtEnergy) - TMath::Log10(fConfig->GetUserERangeMin()))/fBinWidth);
      
      Double_t MOffsetON = (Double_t)OnEvtOffset;
      Double_t MZenON = (Double_t)OnZen;

      // Event information for light curves, firt argument tells that we deal with an On event
      START::Event *OnEvt = new START::Event(true,OnEvtEnergy,STARTUtils::GetMJDFromSashSeconds(OnEvtTimeSec));
      //std::cout << "JLK: Energy=" << OnEvtEnergy << " Time=" << STARTUtils::GetMJDFromSashSeconds(OnEvtTimeSec) << std::endl;
      //OnEvt->Print();
      //SP: il y a-t-il une raison pour qu'on n'utilise pas les offset moyen associés à BgMaker ON et OFF ? Vincent ?      
      //SP: Eventuels problemes liés à des stats limitées.
      if(UseTwoZenBands==false) {
	BandArray[fNbBands].ebin[Ebin].SetOn(BandArray[fNbBands].ebin[Ebin].GetOn()+1);
	BandArray[fNbBands].ebin[Ebin].AddEvent(*OnEvt);
	BandArray[fNbBands].SetOffset(BandArray[fNbBands].GetOffset()+MOffsetON);
	BandArray[fNbBands].SetZenON(BandArray[fNbBands].GetZenON()+MZenON);
      }
      else {
	if(Zbin==0) {
	  BandArray[fNbBands-1].ebin[Ebin].SetOn(BandArray[fNbBands-1].ebin[Ebin].GetOn()+1);
	  BandArray[fNbBands-1].ebin[Ebin].AddEvent(*OnEvt);
	  BandArray[fNbBands-1].SetOffset(BandArray[fNbBands-1].GetOffset()+MOffsetON);
	  BandArray[fNbBands-1].SetZenON(BandArray[fNbBands-1].GetZenON()+MZenON);
	}
	if(Zbin==1) {
	  BandArray[fNbBands].ebin[Ebin].SetOn(BandArray[fNbBands].ebin[Ebin].GetOn()+1);
	  BandArray[fNbBands].ebin[Ebin].AddEvent(*OnEvt);
	  BandArray[fNbBands].SetOffset(BandArray[fNbBands].GetOffset()+MOffsetON);
	  BandArray[fNbBands].SetZenON(BandArray[fNbBands].GetZenON()+MZenON);
	}
      }

      delete OnEvt; OnEvt = 0;

    } // loop on ON entries
    

    /****************************************************/
    /******** Fill arrays for BgMakerOff events  ********/
    /****************************************************/

    //SP: vérifier à quoi ça sert. En regionbg le offset moyen on et off devrait être le même...
    double offsetOFF1(0.); // JLK : it will be the mean offset for the band 
    double offsetOFF2(0.); // event in case of zero stat for on events 

    for(int i=0;i<nbentriesOff;i++)
      {
	offtime_branch_secs->GetEntry(i);
	offtime_branch_nano->GetEntry(i);
	offzen_branch->GetEntry(i);
	Energy_branch_OFF->GetEntry(i);
	EvtOffset_OFF->GetEntry(i);

	Double_t CurrentTimeOff = (Double_t)OffEvtTimeSec+ ((Double_t)OffEvtTimeNano)/1e9;

	/* When needed, we check if the events are in the required time windows  */
	if (fConfig->GetUserMJDwindow().size()>0 && !CheckIfEventInTimeWindow(CurrentTimeOff))
	  continue;

	Int_t Zbin=0; // Define the band to use
	if (UseTwoZenBands==true && CurrentTimeOff>MiddleRunTime) Zbin=1;	    
 
	if (OffEvtEnergy<fConfig->GetUserERangeMin())
	  continue; //Skip events with energy lower than requested
	if (OffEvtEnergy>fConfig->GetUserERangeMax())
	  continue; //Skip events with energy higher than requested

	Int_t Ebin = (Int_t)((TMath::Log10(OffEvtEnergy) - TMath::Log10(fConfig->GetUserERangeMin()))/fBinWidth);
	Double_t MZenOFF = (Double_t)OffZen;	      

	// Event information for light curves, firt argument tells that we deal with an Off event
	START::Event *OffEvt = new START::Event(false,OffEvtEnergy,STARTUtils::GetMJDFromSashSeconds(OffEvtTimeSec));
	//std::cout << "JLK: Energy=" << OffEvtEnergy << " Time=" << STARTUtils::GetMJDFromSashSeconds(OffEvtTimeSec) << std::endl;
	//OffEvt->Print();
	if(UseTwoZenBands==false) {
	  BandArray[fNbBands].ebin[Ebin].SetOff(BandArray[fNbBands].ebin[Ebin].GetOff()+1);
	  BandArray[fNbBands].ebin[Ebin].AddEvent(*OffEvt);
	  BandArray[fNbBands].SetZenOFF(BandArray[fNbBands].GetZenOFF()+MZenOFF);
	  offsetOFF1+=OffEvtOffset; // we fill the off offset
	}
	else {
	  if(Zbin==0) {
	    BandArray[fNbBands-1].ebin[Ebin].SetOff(BandArray[fNbBands-1].ebin[Ebin].GetOff()+1);
	    BandArray[fNbBands-1].ebin[Ebin].AddEvent(*OffEvt);
	    BandArray[fNbBands-1].SetZenOFF(BandArray[fNbBands-1].GetZenOFF()+MZenOFF);
	    offsetOFF1+=OffEvtOffset; // we fill the off offset
	  }
	  if(Zbin==1) {
	    BandArray[fNbBands].ebin[Ebin].SetOff(BandArray[fNbBands].ebin[Ebin].GetOff()+1);
	    BandArray[fNbBands].ebin[Ebin].AddEvent(*OffEvt);
	    BandArray[fNbBands].SetZenOFF(BandArray[fNbBands].GetZenOFF()+MZenOFF);
	    offsetOFF2+=OffEvtOffset; // we fill the off offset
	  }
	}
	
	delete OffEvt; OffEvt = 0;
      } // loop on OFF entries
    
    /*****************************************************************************/
    /******** Determine mean values over energies for each zband (1 or 2) ********/
    /********    Skip runs with time or offset outside user requirements  ********/
    /*****************************************************************************/

    for(int zbin=0;zbin<RunNumberOfZBands;zbin++) {

      TString BandOrRun = "";
      if (RunNumberOfZBands==1) BandOrRun = "RUN";
      if (RunNumberOfZBands==2) BandOrRun = "BAND";
      
      // Check if time of the run is inside the limits fixed by the user
            
      //SP ICI. DEUX CHOSES: on peut faire probablement plus court et simple
      //SP: dans Bands: ajouter fDuration ? 
      if(UseTwoZenBands==false) {
	if(BandArray[fNbBands].GetLiveTime()<(fConfig->GetUserTmin()*LiveTimeFraction)) {
	  WARNING << "Run " << fRunList[irun] << ": TimeAtZOn = " 
		  << (BandArray[fNbBands].GetLiveTime()/LiveTimeFraction) << " < UserTmin = " << fConfig->GetUserTmin()
		  << " => SKIP" << BandOrRun << std::endl;
	  BandArray[fNbBands].SetKeepBand(0);
	  continue;
	}
      }
      else {
	if(zbin==0) {
	  if(BandArray[fNbBands-1].GetLiveTime()<(fConfig->GetUserTmin()*LiveTimeFraction)) {
	    WARNING << "Run " << fRunList[irun] << " TimeAtZOn = " 
		    << (BandArray[fNbBands-1].GetLiveTime()/LiveTimeFraction) << " < UserTmin = " << fConfig->GetUserTmin()
		    << " => SKIP" << BandOrRun << std::endl;
	    BandArray[fNbBands-1].SetKeepBand(0);
	    continue;
	  }
	}
	else if(zbin==1) {
	  if(BandArray[fNbBands].GetLiveTime()<(fConfig->GetUserTmin()*LiveTimeFraction)) {
	    WARNING << "Run " << fRunList[irun] << " TimeAtZOn = " 
		    << (BandArray[fNbBands].GetLiveTime()/LiveTimeFraction) << " < UserTmin = " << fConfig->GetUserTmin()
		    << " => SKIP" << BandOrRun << std::endl;
	    BandArray[fNbBands].SetKeepBand(0);
	    continue;
	  }
	}
      }
      
      // Set livetime TON in hours in band
      if(UseTwoZenBands==false) {
	BandArray[fNbBands].SetLiveTime(BandArray[fNbBands].GetLiveTime()/3600.);
	BandArray[fNbBands].InitEbinLiveTime(BandArray[fNbBands].GetLiveTime());	
      }
      else {
	if(zbin==0) {
	  BandArray[fNbBands-1].SetLiveTime(BandArray[fNbBands-1].GetLiveTime()/3600.);
	  BandArray[fNbBands-1].InitEbinLiveTime(BandArray[fNbBands-1].GetLiveTime());
	}
	if(zbin==1) {
	  BandArray[fNbBands].SetLiveTime(BandArray[fNbBands].GetLiveTime()/3600.);
	  BandArray[fNbBands].InitEbinLiveTime(BandArray[fNbBands].GetLiveTime());
	}
      }
      
      // Set alpha in band
      
      if(UseTwoZenBands==false) {
	BandArray[fNbBands].SetAlphaRun(AlphaRUN);
	BandArray[fNbBands].InitEbinAlpha(AlphaRUN);
      }
      else {
	if(zbin==0) {
	  BandArray[fNbBands-1].SetAlphaRun(AlphaRUN);
	  BandArray[fNbBands-1].InitEbinAlpha(AlphaRUN);
	}  
	if(zbin==1) {
	  BandArray[fNbBands].SetAlphaRun(AlphaRUN);
	  BandArray[fNbBands].InitEbinAlpha(AlphaRUN);
	}
      }
      
      // Set efficacity in band
      
      if(UseTwoZenBands==false) {
	BandArray[fNbBands].SetEff(MuonEffRUN*100);
      }
      else {
	if(zbin==0) BandArray[fNbBands-1].SetEff(MuonEffRUN*100);
	if(zbin==1) BandArray[fNbBands].SetEff(MuonEffRUN*100);
      }
      
      // Set telcode in band
      
      if(UseTwoZenBands==false) {
	BandArray[fNbBands].SetTelCode(fRunTelType[irun]);
      }
      else {
	if(zbin==0) BandArray[fNbBands-1].SetTelCode(fRunTelType[irun]);
	if(zbin==1) BandArray[fNbBands].SetTelCode(fRunTelType[irun]);
      }


      double nbOnBand1(0);
      double nbOnBand2(0);     // Number of events ON and OFF 
      double nbOffBand1(0);    // in band 1 or 2
      double nbOffBand2(0);

      // Compute number of ON and OFF in the bands
	
      if(UseTwoZenBands==false) {
	for(int ebin=0;ebin<fBinNumber;ebin++) {
	  nbOnBand1+=BandArray[fNbBands].ebin[ebin].GetOn();
	  nbOffBand1+=BandArray[fNbBands].ebin[ebin].GetOff();
	}
      }
      else {
	if(zbin==0) {
	  for(int ebin=0;ebin<fBinNumber;ebin++) {
	    nbOnBand1+=BandArray[fNbBands-1].ebin[ebin].GetOn();
	    nbOffBand1+=BandArray[fNbBands-1].ebin[ebin].GetOff();
	  }
	}
	else if(zbin==1) {
	  for(int ebin=0;ebin<fBinNumber;ebin++) {
	    nbOnBand2+=BandArray[fNbBands].ebin[ebin].GetOn();
	    nbOffBand2+=BandArray[fNbBands].ebin[ebin].GetOff();
	  }
	}
      }
      
      // Fill zenith angle of the band
      
      if(UseTwoZenBands==false) {
	if(nbOffBand1>0) BandArray[fNbBands].SetZenOFF((double)(BandArray[fNbBands].GetZenOFF()/nbOffBand1));
	else if(!fConfig->GetUserIsMc()) {
	  //ADA 08/06/2016 
	  if(!nbOnBand1>0){
	    WARNING <<"Run "<<fRunList[irun]<<" NoffBand=" << nbOffBand1 <<" and NOnBand = " << nbOnBand1<< "=> SKIP "<<BandOrRun <<std::endl;
	    BandArray[fNbBands].SetKeepBand(0);
	    continue;
	  }
	}
      }
      else {
	if(zbin==0 && nbOffBand1>0) BandArray[fNbBands-1].SetZenOFF((BandArray[fNbBands-1].GetZenOFF()/nbOffBand1));
	else if(!fConfig->GetUserIsMc() && zbin==0) {
	  //ADA 08/06/2016 
	  if(!nbOnBand1>0){
	    WARNING <<"Run "<<fRunList[irun]<<" NoffBand=" << nbOffBand1 <<" and NOnBand = " << nbOnBand1<< "=> SKIP "<<BandOrRun <<std::endl;
	    BandArray[fNbBands-1].SetKeepBand(0);
	    continue;
	  }
	}
	if(zbin==1 && nbOffBand2>0) BandArray[fNbBands].SetZenOFF((double)(BandArray[fNbBands].GetZenOFF()/nbOffBand2));
	else if(!fConfig->GetUserIsMc() && zbin==1) {
	  //ADA 08/06/2016 
	  if(!nbOnBand2>0){
	    WARNING <<"Run "<<fRunList[irun]<<" NoffBand2=" << nbOffBand2 <<" and NOnBand2 = " << nbOnBand2<< "=> SKIP "<<BandOrRun <<std::endl;
	    BandArray[fNbBands].SetKeepBand(0);
	    continue;
	  }
	}
      }
      // Fill offset of the band
      
      if(UseTwoZenBands==false) {
	if(nbOnBand1>0) {
	  BandArray[fNbBands].SetOffset((double)(BandArray[fNbBands].GetOffset()/nbOnBand1));
	  BandArray[fNbBands].SetZenON((double)(BandArray[fNbBands].GetZenON()/nbOnBand1));
	}
	else if(nbOffBand1>0) {
	  WARNING << "Zero ON events for band so we are using offset from OFF events" <<std::endl;
	  BandArray[fNbBands].SetOffset(offsetOFF1/nbOffBand1);
	  BandArray[fNbBands].SetZenON(BandArray[fNbBands].GetZenOFF());
	}
	else {
	  WARNING << "Zero ON and OFF events so we are using an offset of 0.5 for band" << std::endl;
	  BandArray[fNbBands].SetOffset(0.5); //SP: ATTENTION PROVISOIRE AGN... pas bon ça, y revenir
	  BandArray[fNbBands].SetZenON(BandArray[fNbBands].GetZenOFF());
	}
      }
      else {
	if(zbin==0 && nbOnBand1>0) { //The offset should be taken from run ?
	  BandArray[fNbBands-1].SetOffset((double)(BandArray[fNbBands-1].GetOffset()/nbOnBand1));
	  BandArray[fNbBands-1].SetZenON((double)(BandArray[fNbBands-1].GetZenON()/nbOnBand1));
	}
	else if(zbin==0 && nbOffBand1>0) {
	  WARNING << "Zero ON events for band so we are using offset from OFF events" << std::endl;
	  BandArray[fNbBands-1].SetOffset(offsetOFF1/nbOffBand1);
	  BandArray[fNbBands-1].SetZenON(BandArray[fNbBands-1].GetZenOFF());
	}
	else if(zbin==1 && nbOnBand2>0) {
	  BandArray[fNbBands].SetOffset((double)(BandArray[fNbBands].GetOffset()/nbOnBand2));
	  BandArray[fNbBands].SetZenON((double)(BandArray[fNbBands].GetZenON()/nbOnBand2));
	}
	else if(zbin==1 && nbOffBand2>0) {
	  WARNING << "Zero ON events for band so we are using offset from OFF events" << std::endl;
	  BandArray[fNbBands].SetOffset(offsetOFF2/nbOffBand2);
	  BandArray[fNbBands].SetZenON(BandArray[fNbBands].GetZenOFF());
	}
	else {
	  if(zbin==0) { 
	    WARNING << "Zero ON and OFF events so we are using an offset of 0.5 for band" << std::endl;
	    BandArray[fNbBands-1].SetOffset(0.5); //SP: ATTENTION PROVISOIRE AGN
	    BandArray[fNbBands-1].SetZenON(BandArray[fNbBands-1].GetZenOFF());
	  }
	  if(zbin==1) {
	    WARNING << "Zero ON and OFF events so we are using an offset of 0.5 for band" << std::endl;
	    BandArray[fNbBands].SetOffset((double)(0.5)); //SP: ATTENTION PROVISOIRE AGN
	    BandArray[fNbBands].SetZenON(BandArray[fNbBands-1].GetZenOFF());
	  }
	}
      }
      
      // Check if zenith angle of the band is lower than the one fixed by the user
      
      if(UseTwoZenBands==false) {
	if(fConfig->GetUserZenithMax()>0 && BandArray[fNbBands].GetZenON()>fConfig->GetUserZenithMax()) {
	  WARNING <<"Run "<<fRunList[irun]<<" zenith > zenmax : "
		  <<BandArray[fNbBands].GetZenON()<<">"<<fConfig->GetUserZenithMax()<<" => SKIP "<<BandOrRun <<std::endl;
	  BandArray[fNbBands].SetKeepBand(0);
	  continue;
	}
      }
      else {
	if(zbin==0) {
	  if(fConfig->GetUserZenithMax()>0 && BandArray[fNbBands-1].GetZenON()>fConfig->GetUserZenithMax()) {
	    WARNING <<"Run "<<fRunList[irun]<<" zenith > zenmax : "
		    <<BandArray[fNbBands-1].GetZenON()<<">"<<fConfig->GetUserZenithMax()<<" => SKIP "<<BandOrRun <<std::endl;
	    BandArray[fNbBands-1].SetKeepBand(0);
	    continue;
	  }
	}
	else if(zbin==1) {
	  if(fConfig->GetUserZenithMax()>0 && BandArray[fNbBands].GetZenON()>fConfig->GetUserZenithMax()) {
	    WARNING <<"Run "<<fRunList[irun]<<": zenith > zenmax : "
		    <<BandArray[fNbBands].GetZenON()<<">"<<fConfig->GetUserZenithMax()<<" => SKIP "<<BandOrRun <<std::endl;
	    BandArray[fNbBands].SetKeepBand(0);
	    continue;
	  }
	}
      }
      
      // Check if offset of the band is inside the limits fixed by the user
      
      if(UseTwoZenBands==false) {
	if(BandArray[fNbBands].GetOffset() > fConfig->GetUserOffsetMax()) {
	  WARNING <<"Run "<<fRunList[irun]<<" offset > offmax : "
		  << BandArray[fNbBands].GetOffset() <<">"<<fConfig->GetUserOffsetMax()<<" => SKIP "<<BandOrRun <<std::endl;
	  BandArray[fNbBands].SetKeepBand(0);
	  continue;
	}
      }
      else {
	if(zbin==0) {
	  if(BandArray[fNbBands-1].GetOffset() > fConfig->GetUserOffsetMax()) {
	    WARNING <<"Run "<<fRunList[irun]<<" offset > offmax : "
		    << BandArray[fNbBands-1].GetOffset() <<">"<<fConfig->GetUserOffsetMax()<<" => SKIP "<<BandOrRun <<std::endl;
	    BandArray[fNbBands-1].SetKeepBand(0);
	    continue;
	  }
	}
	else if(zbin==1) {
	  if(BandArray[fNbBands].GetOffset() > fConfig->GetUserOffsetMax()) {
	    WARNING <<"Run "<<fRunList[irun]<<" offset > offmax : "
		    << BandArray[fNbBands].GetOffset() <<">"<<fConfig->GetUserOffsetMax()<<" => SKIP "<<BandOrRun << std::endl;
	    BandArray[fNbBands].SetKeepBand(0);
	    continue;
	  }
	}
      }
      
    }
    
    delete EventsTTree_BgMakerOn;
    EventsTTree_BgMakerOn = 0;
    delete EventsTTree_BgMakerOff;
    EventsTTree_BgMakerOff=0;
    delete RunInfo_TTree;
    RunInfo_TTree=0;
    
    treefile->Close();

  } // loop on run
  
  if(DEBUG>0) {
    for(std::vector<Band>::const_iterator band=BandArray.begin(); band!=BandArray.end(); ++band)
      band->Print();
  }
  
  if(BandArray.size()==0) {
    WARNING << "Reading data from tree... FAILED" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  CheckAndSetIfBandsAreInMCLimits(BandArray); // check if band's parameters are in MC range
  
  SetFitEnergyRangeFlagInEnergyBins(BandArray); // set userfit emin and userfit emax flags in bins

  if(!CheckForSelectedBands(BandArray)) {
    INFO << "All bands are rejected!" << std::endl;
    WARNING << "Reading data from tree... FAILED" << std::endl;
    exit(EXIT_FAILURE);
  }
  else INFO << "Reading data from tree... ok" << std::endl;

  for (std::vector<Band>::iterator it_band = BandArray.begin(); it_band!=BandArray.end();++it_band) {
    it_band->UpdateNameAndTitle(); // Put the correct name for each band !
  }

  return 1; 
							   							   
}

/**
 * \brief Check if there is at least one accepeted band. Return true if there 
 * is no accepted band.
 *
 * \param Array of band
 */
bool START::BandsFactory::CheckForSelectedBands(const std::vector<Band> &BandArray) const {

  int accepetedbands(0);
  
  for(unsigned int iband(0); iband<BandArray.size(); iband++) {
    if(BandArray[iband].GetKeepBand()==1) accepetedbands++;
  }
  DEBUG_OUT << "acceptedbands=" << accepetedbands << std::endl;
  if(accepetedbands<1) return false;
  else return true;

  return true;
}

/**
 * \brief We check if the band's parameters are within the limits of the MC values.
 * If not, we skipped the bands
 *
 * \param Array of band
 */
void START::BandsFactory::CheckAndSetIfBandsAreInMCLimits(std::vector<Band> &BandArray) const {

  MonteCarlo MC(fIrfOpt);
  
  for(unsigned int iband(0); iband<BandArray.size(); iband++) {
    
    // offset
    
    if(BandArray[iband].GetOffset()<MC.GetOffset().front() || BandArray[iband].GetOffset()>MC.GetOffset().back()) {
      BandArray[iband].SetKeepBand(0);
      if(BandArray[iband].GetOffset()<MC.GetOffset().front()) {
	WARNING << "Run " << BandArray[iband].GetNbRun()
		<<" offset < offset min MC : "
		<< BandArray[iband].GetOffset() <<"<"<< MC.GetOffset().front()
		<<" => SKIP "<< "Run" <<std::endl;
      }
      else {
	WARNING << "Run " << BandArray[iband].GetNbRun()
		<<"  offset > offset max MC : "
		<< BandArray[iband].GetOffset() <<">"<< MC.GetOffset().back()
		<<" => SKIP "<< "Run" <<std::endl;
      }
    }

    // zenith

    if(BandArray[iband].GetZenON()<MC.GetZenith().front() || BandArray[iband].GetZenON()>MC.GetZenith().back()) {
      BandArray[iband].SetKeepBand(0);
      if(BandArray[iband].GetZenON()<MC.GetZenith().front()) {
	WARNING << "Run " << BandArray[iband].GetNbRun()
		<<" zenith < zenith min MC : "
		<< BandArray[iband].GetZenON() <<"<"<< MC.GetZenith().front()
		<<" => SKIP "<< "Run" <<std::endl;
      }
      else {
	WARNING << "Run " << BandArray[iband].GetNbRun()
		<<" zenith > zenith max MC : "
		<< BandArray[iband].GetZenON() <<">"<< MC.GetZenith().back()
		<<" => SKIP "<< "Run" <<std::endl;
      }
    }

    for(unsigned int i(0); i<BandArray.size(); i++) {
      if(BandArray[i].GetEff()>45. && BandArray[i].GetEff()<50.) {
	WARNING << "JLK : the efficiency of the band is " << BandArray[i].GetEff() 
		<< " so while we wait for Santiago and Vincent to fix HandleResolArea, \n"
		<< "we put the efficiency to 50. Sorry..." << std::endl; //SP revoir ???
	
	BandArray[i].SetEff(50.);	
      }
    }
    
    // efficiency
    if(MC.GetEfficiency().size()  >= 2){
      if(BandArray[iband].GetEff()<MC.GetEfficiency().front() || BandArray[iband].GetEff()>MC.GetEfficiency().back()) {
      
	BandArray[iband].SetKeepBand(0);

	if(BandArray[iband].GetEff()<MC.GetEfficiency().front()) {
	  WARNING << "Run " << BandArray[iband].GetNbRun()
		  <<" efficiency < efficiency min MC : "
		  << BandArray[iband].GetEff() <<"<"<< MC.GetEfficiency().front()
		  <<" => SKIP "<< "Run" << std::endl;
	}
	else {
	  WARNING << "Run " << BandArray[iband].GetNbRun()
		  <<" efficiency > efficiency max MC : "
		  << BandArray[iband].GetEff() <<">"<< MC.GetEfficiency().back()
		  <<" => SKIP "<< "Run" <<std::endl;
	}
      }
    }

  }
  

}

/**
 * \brief if energy range if specified by the user, we keep only energybins
 * which satisfies efitmin<ebinmean<efitmax
 *
 * \param BandArray vector of bands
 *
 * \todo JLK : Change condition to ebinmin>= efitmin && ebinmax<= efitmax?
 */
void START::BandsFactory::SetFitEnergyRangeFlagInEnergyBins(std::vector<Band> &BandArray) const {

  double fitenergymin = fConfig->GetUserFitEMin();
  double fitenergymax = fConfig->GetUserFitEMax();  

  DEBUG_OUT << "fitenergymin=" << fitenergymin << std::endl;
  DEBUG_OUT << "fitenergymax=" << fitenergymax << std::endl;

  for(std::vector<Band>::iterator band = BandArray.begin(); band!=BandArray.end(); ++band) {

    if(fitenergymin==-1 && fitenergymax==-1) continue; // no fit energy range specified by the user

    if(band->GetKeepBand()==0) continue; // interested in selected bands

    for(std::vector<EnergyBin>::iterator bin = band->ebin.begin(); bin!=band->ebin.end(); ++bin) {

      if(bin->GetEmid()>=fitenergymin && bin->GetEmid()<=fitenergymax) bin->SetKeepBin(1); // JLK change to emin-emax from bin?
      else bin->SetKeepBin(0);

      DEBUG_OUT << "emean=" << bin->GetEmid() << " keep=" << bin->GetKeepBin() << std::endl;

    }

  }

}

/**
 * \brief Print All the infos contained in the bands
 *
 * \todo : Add the possibility to have the result above the threshold energy
 */
void START::BandsFactory::PrintBands(const std::vector<Band> &BandArray) const {
  int nbrun = fConfig->GetRunList().size();
  int nbbands = BandArray.size();
  int acceptedbands(0);
  double number_on_used(0);
  double number_off_used(0);
  double number_on(0);
  double number_off(0);

  for(unsigned int iband(0); iband<BandArray.size(); iband++) {

    if(BandArray[iband].GetKeepBand()==1) {
      acceptedbands++;
    }
    
    for(unsigned int ibin(0); ibin<BandArray[iband].ebin.size(); ibin++) {

      number_on+=BandArray[iband].ebin[ibin].GetOn();
      number_off+=BandArray[iband].ebin[ibin].GetOff();

      if(BandArray[iband].ebin[ibin].GetKeepBin() == 0) continue;
      
      if(BandArray[iband].GetKeepBand()==1) {
	number_on_used+=BandArray[iband].ebin[ibin].GetOn();
	number_off_used+=BandArray[iband].ebin[ibin].GetOff();
      }

    }

  }
  
  INFO << "Bands caracteristics : " << std::endl;
  for(int i(0); i<nbbands; i++) {
    BandArray[i].Print();
  }

  PrintSummaryBand(BandArray);

  INFO << "Info on bands :" << std::endl;
  std::cout << "Number of runs = " << nbrun << std::endl;
  std::cout << "Number of bands = " << nbbands << std::endl;
  std::cout << "Number of accepted bands = " << acceptedbands << std::endl;
  std::cout << "Number of ON events = " << number_on << std::endl;
  std::cout << "Number of OFF events = " << number_off << std::endl;

}


/** AL : Reproj the band Array in Zen,alpha,theta,ect....
 *\ for very low statistics
 */

void START::BandsFactory::FillBandInfo(std::vector<Band> &BandArray,std::map<TString, std::vector<int> > &InfoReprojArray,bool LowStat,bool VeryLowStat, Config myConfig)const
{

  //std::cout<<" *************Fill BandInfo and reproj !  "<<std::endl;

  double Eff_max=myConfig.GetBandEffMax();
  double Eff_min=myConfig.GetBandEffMin();
  int Eff_bin=myConfig.GetBandEffbin_UseReproj();
  double Eff_bin_size=(Eff_max-Eff_min)/Eff_bin;
 
  double pi=3.14159265359;
  //G0P1 -BG1 (400 excess):
  double Zen_max=cos(myConfig.GetBandZenMin()*pi/180.);//0deg
  double Zen_min=cos(myConfig.GetBandZenMax()*pi/180.);//70deg -> tot interval = 0.66 
  int Zen_bin=myConfig.GetBandZenbin_UseReproj();
  double Zen_bin_size=(Zen_max-Zen_min)/Zen_bin;

  double theta_max=myConfig.GetBandOffMax();
  double theta_min=myConfig.GetBandOffMin();
  int theta_bin=myConfig.GetBandOffbin_UseReproj();
  double theta_bin_size=(theta_max-theta_min)/theta_bin;
 
  if(LowStat){
    //std::cout<<"Low Stat !!!  -> LARGE BANDS !!!! "<<std::endl;
       
      
    Eff_bin=myConfig.GetBandEffbin_SmallStat();
    Eff_bin_size=(Eff_max-Eff_min)/Eff_bin;
    

    Zen_bin=myConfig.GetBandZenbin_SmallStat();
    Zen_bin_size=(Zen_max-Zen_min)/Zen_bin;
 

    theta_bin=myConfig.GetBandOffbin_SmallStat();
    theta_bin_size=(theta_max-theta_min)/theta_bin;


  }

  if(VeryLowStat){
    //std::cout<<"VEry Low Stat !!!  -> VERY LARGE BANDS !!!! "<<std::endl;
       
      
    Eff_bin=myConfig.GetBandEffbin_VerySmallStat();
    Eff_bin_size=(Eff_max-Eff_min)/Eff_bin;
    

    Zen_bin=myConfig.GetBandZenbin_VerySmallStat();;
    Zen_bin_size=(Zen_max-Zen_min)/Zen_bin;
 

    theta_bin=myConfig.GetBandOffbin_VerySmallStat();
    theta_bin_size=(theta_max-theta_min)/theta_bin;

  }
  
  int z_out=0.;
  int e_out=0.;
  int t_out=0.;
  
  
  //it est accrement� de 1 a chaque run
  int it=0;
  for(std::vector<Band>::const_iterator band=BandArray.begin(); band!=BandArray.end(); ++band) 
    //for(int iband=0; iband< BandArray.size();iband++)
    {

      //commente car ca risque de decaller les indices it pour se reperer dans Bandarray
      //   if(band->GetKeepBand()==0)
      // 	{
      //	  std::cout<<"Band->GetKeepBand is false !   NbRun =  "<<band->GetNbRun() <<"   AlphaRun = "<<band->GetAlphaRun()<<std::endl;
      //	  continue;
      //	}
      
      //      std::cout<<"band = "<<band<<"entier it = "<<it<<std::endl;
      
      for(int e=0;e<Eff_bin;e++)
	{
	  double Eff=Eff_min+(Eff_bin_size/2.0)+ (e*Eff_bin_size);
	  //std::cout<<" band->GetEff() = "<<band->GetEff() <<"   Eff = "<<Eff <<"   e = "<<e<<std::endl;
	  if((band->GetEff()>(Eff -(Eff_bin_size/2.0))) && (band->GetEff()<(Eff +(Eff_bin_size/2.0))))	
	    { 
	      e_out=e;
	      break;	
	    }
	}
     
      
      for(int z=0;z<Zen_bin;z++)
	{
	  double Zen=Zen_min+(Zen_bin_size/2.0) + (z*Zen_bin_size);
	  //std::cout<<" band->GetZenON() = "<<band->GetZenON() <<"   Zen = "<<Zen <<"   z = "<<z<<std::endl;
	
	 
	  double CosZen = TMath::Cos(band->GetZenON()*TMath::Pi()/180.);
	  
	  if((CosZen>(Zen -(Zen_bin_size/2.0))) && (CosZen<(Zen +(Zen_bin_size/2.0))))	
	    { 
	      z_out=z;
	      break;	
	    }	  
	}
      
      for(int t=0;t<theta_bin;t++)  
	{
	  double Theta=theta_min+(theta_bin_size/2.0)+(t*theta_bin_size);
	  //std::cout<<" band->GetOffset() = "<<band->GetOffset() <<"   Theta = "<<Theta <<"   t = "<<t<<std::endl;
	  if((band->GetOffset()>(Theta -(theta_bin_size/2.0))) && (band->GetOffset()<(Theta +(theta_bin_size/2.0))))	
	    { 
	      t_out=t;
	      break;	
	    }
	}


	
     
      std::map<TString, std::vector<int> >::iterator im;
      //defini la key pour la map
      std::ostringstream name;
      name << "Band_Eff" << e_out << "_CosZen" <<z_out << "_Offset" <<t_out  ;

      //std::cout<<"name = "<<name.str().c_str() <<"  band->GetEff() = "<<band->GetEff() <<"  band->GetZenON() = "<<band->GetZenON() <<"  band->GetOffset() = "<<band->GetOffset() <<std::endl;

      TString key(name.str().c_str());
      //Et ben la fonction map::find renvoie l'it�rateur end() dans ce cas, end() qui renvoie un it�rateur vers le dernier �l�ment de la s�quence
      im=InfoReprojArray.find(key);
      
   
      //si im est differend du .end() de Info reproj, il ajoute a la bande en eff,zen et offset deja existantes, le run. pusk back ajoute une ligne au vector avec le nombre=numero du run a ajoute dans la bande.
      //sinon si c'est egal au dernier element de Inforeproj, il cree une nouvelle key string associee au nom de la abnde en eff,zen et off qui n existait pas encore avant et lui attribue le numero dun run 
      if(im!=InfoReprojArray.end())
	{
	  if(it>0)
	    //std::cout<<"(*im).first = "<<(*im).first <<"  (*im).second.size() = "<<(*im).second.size() <<std::endl; 
	    InfoReprojArray[(*im).first].push_back(it);
	}
      else 
	InfoReprojArray[key].push_back(it);
      
      
      it++;
    }
  
  //Check if infoArray is ok:
  //std::cout<<"InfoReprojArray.size() = "<<InfoReprojArray.size() <<"   it = "<<it <<std::endl;
  /*for( std::map<TString, std::vector<int> >::iterator in =InfoReprojArray.begin(); in !=InfoReprojArray.end(); ++in)
    {
    //std::cout<<"(*in).first = "<<(*in).first <<" taille du vecteur i  (*in).second.size()= "<<(*in).second.size() <<std::endl;
    for(int n=0;n<(*in).second.size();n++)
    //std::cout<<"  element "<<n<<" du vecteur in :  = "<<(*in).second[n]<<std::endl; 
    }*/


 
}

/*


  void START::BandsFactory::ReprojBands(std::vector<Band> &BandArray,
  std::map<TString, std::vector<int> > &InfoReprojArray,
  std::vector<Band> &ReprojArray,
  TString configname,
  HandleResolArea &HandleIRF)
  {
  int nband=0;
  //loop on InfoReprojArray keys
  for(std::map<TString, std::vector<int> > ::iterator i=InfoReprojArray.begin(); i!=InfoReprojArray.end(); i++) 
  {
  //std::cout<<" ***************************************** bande numero : "<<nband <<std::endl;  
  Band *InterBand = new Band(BandArray[0]);
  InterBand->ClearBandInfo();
  InterBand->AddInfoFromSelectedBands(BandArray,(*i).second,configname,false);

      
  if((InterBand->GetNbRun()) !=0)
  {
  //InterBand->Print();
  ReprojArray.push_back(*InterBand);
  }
      
  delete InterBand;
  InterBand=0;
  nband++;
      
  }  
  std::cout << "avant copyinterpolator in band" << std::endl;
  if (HandleIRF.CopyInterpolatorsInBands(ReprojArray)==-1) {
  WARNING << "Assignment of GSLInterpolators in Band... FAILED" << std::endl;
  exit(EXIT_FAILURE);
  }
  else {
  INFO << "Assignment of GSLInterpolators in Band... ok" << std::endl;
  }
  HandleIRF.SetInBandsFirstEmcBin(ReprojArray);
  HandleIRF.SetInBandsMCEnergyThreshold(ReprojArray);
  HandleIRF.SetInBinsKeepBin(ReprojArray);
  HandleIRF.SetInBinsEffectiveArea(ReprojArray);
  //calcul de f(E)
  ComputeResults *CompRes = new ComputeResults(ReprojArray);
  std::cout << "avant makepartialintegrand" << std::endl;
  if(CompRes->MakeVectorPartialIntegral(ReprojArray)==-1) {
  WARNING << "Copying effective resolution in Bins... FAILED" << std::endl;
  exit(EXIT_FAILURE);
  }
  else {
  INFO << "Copying effective resolution in Bins... ok" << std::endl;
  }
  std::cout << "apres makepartialintegrand" << std::endl;
  delete CompRes;
  std::cout <<"avant copyinterpolator in bin" << std::endl;
  if (HandleIRF.CopyInterpolatorsInBins(ReprojArray)==-1) {
  WARNING << "Assignment of GSLInterpolators in Bins... FAILED" << std::endl;
  exit(EXIT_FAILURE);
  }
  else {
  INFO << "Assignment of GSLInterpolators in EnergyBins... ok" << std::endl;
  }
  std::cout <<"apres copyinterpolator in bin" << std::endl;

  }

  void START::BandsFactory::RebinEnergy(const std::vector<Band> &BandArray,
  std::vector<Band> &BandRebinArray,
  double sigrebin,
  double MinE,
  HandleResolArea &HandleIRF) {
 
  //BandRebinArray is a copy of BandArray
  //Then we clear its Ebin
  //Then we copy the rebinned Energy bins for Each Band of the Array


  //std::cout<<" ********* Entering RebinEnergy function ************ "<<std::endl;


  // Sanity check
  if (BandArray.size()==0) {
  WARNING << "Can't rebin anything since the bandarray has no size" << std::endl;
  return;
  }

  Band *SummaryBand = new Band(BandArray[0]);
  SummaryBand->ClearBandInfo();
  SummaryBand->AddInfoFromBands(BandArray,false);

   
  Band *BandSummaryRebin = new Band(*SummaryBand);
  BandSummaryRebin->ebin.clear();

  for (unsigned int iband=0;iband<BandArray.size();++iband)
  BandRebinArray.push_back(BandArray[iband]);
    
  for (unsigned int iband=0;iband<BandArray.size();++iband)
  BandRebinArray[iband].ebin.clear(); 
    
   
  //Rebin BandSummaryRebin first
    
  int bin_first=0;
  int bin_second=0;

  

  while (bin_first<(int)SummaryBand->ebin.size()) {

  //std::cout<<"loop over energy bin : ibin = "<< bin_first<<std::endl;
      
  if (SummaryBand->ebin[bin_first].GetKeepBin()==0) {
  ++bin_first;
  //std::cout<<" --> dont keep bin "<<std::endl;
  continue;
  }


  int bin1=bin_first;
  while (bin1<(int)SummaryBand->ebin.size())
  {
  if(SummaryBand->ebin[bin1].GetKeepBin()==1) 
  ++bin1;
  else
  break;
  }

  //std::cout<<"bin max with 1 = "<<bin1 <<std::endl;

  double sig=-1;
      
  EnergyBin RegroupBin = EnergyBin();
  RegroupBin.SetEmin( SummaryBand->ebin[bin_first].GetEmin() );
  RegroupBin.SetKeepBin(1);
      
  std::vector<EnergyBin> RegroupBin_Array;
  for (unsigned int iband=0;iband<BandArray.size();++iband) 
  { 
  RegroupBin_Array.push_back(EnergyBin());
  RegroupBin_Array[iband].SetEmin(BandArray[iband].ebin[bin_first].GetEmin() );
  RegroupBin_Array[iband].SetKeepBin(1);
  }
      
  bin_second = bin_first;
   
      
  while(( sig<sigrebin && bin_second<(int)SummaryBand->ebin.size() && SummaryBand->ebin[bin_second].GetKeepBin()!=0) ||
  (bin_second<(int)SummaryBand->ebin.size() && (bin1 - bin_second)<(bin_second - bin_first) && SummaryBand->ebin[bin_second].GetKeepBin()!=0)) {	
  //std::cout<<" sig  = "<<sig <<" for a sig threshold of : "<<sigrebin <<"-->  bin_first = "<< bin_first<<"    bin_second = "<<bin_second <<std::endl;
	

  RegroupBin.AddInfoFromEBin(SummaryBand->ebin[bin_second]);
  RegroupBin.SetEmax(SummaryBand->ebin[bin_second].GetEmax());
	
  for (unsigned int iband=0;iband<BandArray.size();++iband){ 
  RegroupBin_Array[iband].AddInfoFromEBin(BandArray[iband].ebin[bin_second]);
  RegroupBin_Array[iband].SetEmax(BandArray[iband].ebin[bin_second].GetEmax());
  }
	
  if (RegroupBin.GetAlpha()!=0) {
  sig = STARTUtils::LiMaSignificance(int(RegroupBin.GetOn()),int(RegroupBin.GetOff()),RegroupBin.GetAlpha());
  }
  else {
  sig=0.;
  }
  //std::cout<<"in the loop : sig for the new bin is = "<<sig  <<" sigma"<<std::endl;
	
  ++bin_second;

  //no rebin before 5 TeV 
  if((SummaryBand->ebin[bin_second].GetEmin())<MinE && bin_second!=bin_first)
  break;
  }
      
  // VIM : We complete the information of the bin (emid and emean), then we add this bin to the Band.
  RegroupBin.SetEmid();
  double emean_bin =RegroupBin.GetEmid();
  RegroupBin.SetEmean(emean_bin);
  BandSummaryRebin->ebin.push_back(RegroupBin);
      
      
  for (unsigned int iband=0;iband<BandArray.size();++iband) { 
  RegroupBin_Array[iband].SetEmid();
  double Emean=RegroupBin_Array[iband].GetEmid();
  //std::cout<<"iband = "<<iband<<"   Emean = "<<Emean <<std::endl;
  RegroupBin_Array[iband].SetEmean(Emean);
  BandRebinArray[iband].ebin.push_back(RegroupBin_Array[iband]); 
  }
      
  bin_first = bin_second;
  }
    
  BandSummaryRebin->SetKeepBand(1);
  for (unsigned int iband=0;iband<BandArray.size();++iband) { 
  if(BandArray[iband].GetKeepBand()==0) continue;
  //if(iband==8) BandRebinArray[iband].SetKeepBand(0);
  BandRebinArray[iband].SetKeepBand(1);
  //std::cout<<"BandRebinArray[iband].ebin.size()"<<BandRebinArray[iband].ebin.size() <<std::endl;
  }
    
  std::cout << "avant copyinterpolator in band" << std::endl;
  if (HandleIRF.CopyInterpolatorsInBands(BandRebinArray)==-1) {
  WARNING << "Assignment of GSLInterpolators in Band... FAILED" << std::endl;
  exit(EXIT_FAILURE);
  }
  else {
  INFO << "Assignment of GSLInterpolators in Band... ok" << std::endl;
  }
  HandleIRF.SetInBandsFirstEmcBin(BandRebinArray);
  HandleIRF.SetInBandsMCEnergyThreshold(BandRebinArray);
  HandleIRF.SetInBinsKeepBin(BandRebinArray);
  HandleIRF.SetInBinsEffectiveArea(BandRebinArray);
  //calcul de f(E)
  ComputeResults *CompRes = new ComputeResults(BandRebinArray);
  std::cout << "avant makepartialintegrand" << std::endl;
  if(CompRes->MakeVectorPartialIntegral(BandRebinArray)==-1) {
  WARNING << "Copying effective resolution in Bins... FAILED" << std::endl;
  exit(EXIT_FAILURE);
  }
  else {
  INFO << "Copying effective resolution in Bins... ok" << std::endl;
  }
  std::cout << "apres makepartialintegrand" << std::endl;
  delete CompRes;
  std::cout <<"avant copyinterpolator in bin" << std::endl;
  if (HandleIRF.CopyInterpolatorsInBins(BandRebinArray)==-1) {
  WARNING << "Assignment of GSLInterpolators in Bins... FAILED" << std::endl;
  exit(EXIT_FAILURE);
  }
  else {
  INFO << "Assignment of GSLInterpolators in EnergyBins... ok" << std::endl;
  }
  std::cout <<"apres copyinterpolator in bin" << std::endl;
    
  }
*/
 
/**
 * \brief Print information of bands such as ON, OFF, exess, time... in one band
 * \param BandArray vector of band
 */
void START::BandsFactory::PrintSummaryBand(const std::vector<Band> &BandArray) const
{

  Band *SummaryBand = new Band(BandArray[0]);
  SummaryBand->ClearBandInfo();
  SummaryBand->AddInfoFromBands(BandArray,false);

  INFO << "SummaryBand : " << std::endl;
  SummaryBand->Print();

  delete SummaryBand;
  SummaryBand=0;

}

/**
 * \brief Check if event's time is in time window defined by the user
 * \param currenttime time of current event
 */
bool START::BandsFactory::CheckIfEventInTimeWindow(double currenttime) {

  //double CurrentMJD = fTimeToAdd + currenttime/(24.*3600.);
  double CurrentMJD = STARTUtils::GetMJDFromSashSeconds(currenttime); // JLK
  for(unsigned int i=0;i<fConfig->GetUserMJDwindow().size();i++) {
    if(CurrentMJD>fConfig->GetUserMJDwindow()[i].first && CurrentMJD<=fConfig->GetUserMJDwindow()[i].second) {
      //if(fConfig->GetUserVerbose()) printf("currenttime=%f  fTimeToAdd=%f        => MJD=%f\n",currenttime,fTimeToAdd,CurrentMJD);
      // To much lines for verbose mode
      DEBUG_OUT_L(3) << "currenttime=" << currenttime << " TimeToAdd=" << fTimeToAdd << " ==>CurrentMJD=" << CurrentMJD << std::endl;
      return true;
    }
  }
  return false;
}

/**
 * \brief Check if run's time is in time window defined by the user
 * \param TimeFirstEvt time of first event
 * \param TimeLastEvt time of last event
 */
bool START::BandsFactory::CheckIfEntireRunIncludedInSingleTimeWindow(double TimeFirstEvt, double TimeLastEvt) {

  bool IsIncluded=false;

  //double TimeFirstEvtMJD = fTimeToAdd + TimeFirstEvt/(24.*3600.);
  //double TimeLastEvtMJD = fTimeToAdd + TimeLastEvt/(24.*3600.);
  double TimeFirstEvtMJD = STARTUtils::GetMJDFromSashSeconds(TimeFirstEvt); // JLK
  double TimeLastEvtMJD = STARTUtils::GetMJDFromSashSeconds(TimeLastEvt); // JLK

  for(unsigned int i=0;i<fConfig->GetUserMJDwindow().size();i++) {
    if(TimeFirstEvtMJD>=fConfig->GetUserMJDwindow()[i].first && TimeLastEvtMJD<=fConfig->GetUserMJDwindow()[i].second) {
      if(fverbose) std::cout<<"Run included in window #"<<(i+1)<<std::endl;
      IsIncluded=true;
    }
  }

  return IsIncluded;
}


void START::BandsFactory::Change_N_ON(std::vector<Band> &BandArray, double n_on) const {
  for(unsigned int iband(0); iband<BandArray.size(); iband++) {
    
    for(unsigned int ibin(0); ibin<BandArray[iband].ebin.size(); ibin++) {

      BandArray[iband].ebin[ibin].SetOn(n_on);
	
    }

  }

}  
  
