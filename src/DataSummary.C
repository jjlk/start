// STL
#include <iostream>
#include <algorithm>
#include <utility>

// ROOT
#include <TF1.h>

// START
#include "DataSummary.hh"
#include "Band.hh"
#include "STARTUtils.hh"

// Utilities
#define DEBUG 0
#include "debugging.hh"

#define INFO std::cout << INFOCOLOR << "DataSummary> " << RESETCOLOR
#define WARNING std::cout << WARNINGCOLOR << "DataSummary> " << RESETCOLOR

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::DataSummary)
#endif

/**
 * \brief Default constructor
 */
START::DataSummary::DataSummary()
:fCanvasEfficiency(0),
	fCanvasOffset(0),
	fCanvasZenithON(0),
	fCanvasZenithOFF(0),
	fCanvasEnergyON(0),
	fCanvasEnergyOFF(0),
	fCanvasLiveTime(0),
	fCanvasAlpha(0),
	fCanvasSignificanceVsLivetime(0),
	fCanvasOffEventsVsLivetime(0),
	fCanvasOnEventsVsLivetime(0),
	fCanvasExcessVsLivetime(0),
	fCanvasOffsetVsLivetime(0),
	fCanvasTimeMJD(0),
	fHistoEfficiency(0),
	fHistoOffset(0),
	fHistoZenithON(0),
	fHistoZenithOFF(0),
	fHistoEnergyON(0),
	fHistoEnergyOFF(0),
	fHistoLiveTime(0),
	fHistoAlpha(0),
	fGraphSignificanceVsLivetime(0),
	fGraphOffEventsVsLivetime(0),
	fGraphOnEventsVsLivetime(0),
	fGraphExcessVsLivetime(0),
	fGraphOffsetVsLivetime(0),
	fHistoTimeMJD(0)
{

}

/**
 * \brief Constructor which compute by default the pairs of MC parameters
 * for all bands and the telescope needed.
 */
START::DataSummary::DataSummary(std::vector<Band> const &BandArray)
	:fCanvasEfficiency(0),
	 fCanvasOffset(0),
	 fCanvasZenithON(0),
	 fCanvasZenithOFF(0),
	 fCanvasEnergyON(0),
	 fCanvasEnergyOFF(0),
	 fCanvasLiveTime(0),
	 fCanvasAlpha(0),
	 fCanvasSignificanceVsLivetime(0),
	 fCanvasOffEventsVsLivetime(0),
	 fCanvasOnEventsVsLivetime(0),
	 fCanvasExcessVsLivetime(0),
	 fCanvasOffsetVsLivetime(0),
	 fCanvasTimeMJD(0),
	 fHistoEfficiency(0),
	 fHistoOffset(0),
	 fHistoZenithON(0),
	 fHistoZenithOFF(0),
	 fHistoEnergyON(0),
	 fHistoEnergyOFF(0),
	 fHistoLiveTime(0),
	 fHistoAlpha(0),
	 fGraphSignificanceVsLivetime(0),
	 fGraphOffEventsVsLivetime(0),
	 fGraphOnEventsVsLivetime(0),
	 fGraphExcessVsLivetime(0),
	 fGraphOffsetVsLivetime(0),
	 fHistoTimeMJD(0)
{
	//MakePairs(BandArray);
	//DetermineTeltype(BandArray);
}
/**
 * \brief Copy Constructor
 * \warning root objects are not copied
 */
START::DataSummary::DataSummary(START::DataSummary const &DataSumCopy)
	:TObject(DataSumCopy),
	 fPairEfficiency(DataSumCopy.fPairEfficiency),
	 fPairZenith(DataSumCopy.fPairZenith),
	 fPairOffset(DataSumCopy.fPairOffset),
	 fEfficiency(DataSumCopy.fEfficiency),
	 fZenith(DataSumCopy.fZenith),
	 fOffset(DataSumCopy.fOffset),
	 fWhichTelType(DataSumCopy.fWhichTelType)
{
  
	fCanvasEfficiency = 0;
	fCanvasOffset = 0;
	fCanvasZenithON = 0;
	fCanvasZenithOFF = 0;
	fCanvasEnergyON = 0;
	fCanvasEnergyOFF = 0;
	fCanvasLiveTime = 0;
	fCanvasAlpha = 0;
	fCanvasSignificanceVsLivetime = 0;
	fCanvasTimeMJD = 0;
	fHistoEfficiency = 0;
	fHistoOffset = 0;
	fHistoZenithON = 0;
	fHistoZenithOFF = 0;
	fHistoEnergyON = 0;
	fHistoEnergyOFF = 0;
	fHistoLiveTime = 0;
	fHistoAlpha = 0;
	fGraphSignificanceVsLivetime = 0;
	fHistoTimeMJD = 0;
}

/**
 * \brief Assignment operator
 * \warning root objects are not copied
 */
START::DataSummary& START::DataSummary::operator=(START::DataSummary const &DataSumCopy)
{
	if(this != &DataSumCopy) {
		fPairEfficiency=DataSumCopy.fPairEfficiency;
		fPairZenith=DataSumCopy.fPairZenith;
		fPairOffset=DataSumCopy.fPairOffset;
		fZenith=DataSumCopy.fZenith;
		fOffset=DataSumCopy.fOffset;
		fEfficiency=DataSumCopy.fEfficiency;
		fWhichTelType=DataSumCopy.fWhichTelType;

		fCanvasEfficiency = 0;
		fCanvasOffset = 0;
		fCanvasZenithON = 0;
		fCanvasZenithOFF = 0;
		fCanvasEnergyON = 0;
		fCanvasEnergyOFF = 0;
		fCanvasLiveTime = 0;
		fCanvasAlpha = 0;
		fCanvasSignificanceVsLivetime = 0;
		fCanvasTimeMJD = 0;
		fHistoEfficiency = 0;
		fHistoOffset = 0;
		fHistoZenithON = 0;
		fHistoZenithOFF = 0;
		fHistoEnergyON = 0;
		fHistoEnergyOFF = 0;
		fHistoLiveTime = 0;
		fHistoAlpha = 0;
		fGraphSignificanceVsLivetime = 0;
		fHistoTimeMJD = 0;
	}

	return (*this);
}

/**
 * \brief Destructor
 */
START::DataSummary::~DataSummary()
{

	// JLK : I don't know how to manage the memory and the convinience.
	// So I don't delete those pointers
	CleanROOTObjetcs();
}

/**
 * \brief Delete ROOT' objects
 */
void START::DataSummary::CleanROOTObjetcs()
{

	//CleanCanvas();
	//CleanTH1();
	//CleanGraph();

}

/**
 * \brief Delete TCanvas
 */
void START::DataSummary::CleanCanvas()
{
	if(fCanvasEfficiency!=0) delete fCanvasEfficiency;
	fCanvasEfficiency=0;
	if(fCanvasOffset!=0) delete fCanvasOffset;
	fCanvasOffset=0;
	if(fCanvasZenithON!=0) delete fCanvasZenithON;
	fCanvasZenithON=0;
	if(fCanvasZenithOFF!=0) delete fCanvasZenithOFF;
	fCanvasZenithOFF=0;
	if(fCanvasEnergyON!=0) delete fCanvasEnergyON;
	fCanvasEnergyON=0;
	if(fCanvasEnergyOFF!=0) delete fCanvasEnergyOFF;
	fCanvasEnergyOFF=0;
	if(fCanvasLiveTime!=0) delete fCanvasLiveTime;
	fCanvasLiveTime=0;
	if(fCanvasAlpha!=0) delete fCanvasAlpha;
	fCanvasAlpha=0;
	if(fCanvasSignificanceVsLivetime!=0) delete fCanvasSignificanceVsLivetime;
	fCanvasSignificanceVsLivetime=0;
	if(fCanvasOffEventsVsLivetime!=0) delete fCanvasOffEventsVsLivetime;
	fCanvasOffEventsVsLivetime=0;
	if(fCanvasOnEventsVsLivetime!=0) delete fCanvasOnEventsVsLivetime;
	fCanvasOnEventsVsLivetime=0;
	if(fCanvasExcessVsLivetime!=0) delete fCanvasExcessVsLivetime;
	fCanvasExcessVsLivetime=0;
	if(fCanvasOffsetVsLivetime!=0) delete fCanvasOffsetVsLivetime;
	fCanvasOffsetVsLivetime=0;
	if(fCanvasTimeMJD!=0) delete fCanvasTimeMJD;
	fCanvasTimeMJD = 0;
}

/**
 * \brief Delete TH1
 */
void START::DataSummary::CleanTH1()
{

	if(fHistoEfficiency!=0) delete fHistoEfficiency;
	fHistoEfficiency=0;
	if(fHistoOffset!=0) delete fHistoOffset;
	fHistoOffset=0;
	if(fHistoZenithON!=0) delete fHistoZenithON;
	fHistoZenithON=0;
	if(fHistoZenithOFF!=0) delete fHistoZenithOFF;
	fHistoZenithOFF=0;
	if(fHistoEnergyON!=0) delete fHistoEnergyON;
	fHistoEnergyON=0;
	if(fHistoEnergyOFF!=0) delete fHistoEnergyOFF;
	fHistoEnergyOFF=0;
	if(fHistoAlpha!=0) delete fHistoAlpha;
	fHistoAlpha=0;
	if(fHistoLiveTime!=0) delete fHistoLiveTime;
	fHistoLiveTime=0;
	if(fHistoTimeMJD!=0) delete fHistoTimeMJD;
	fHistoTimeMJD = 0;

}

/**
 * \brief Delete TH1
 */
void START::DataSummary::CleanGraph()
{

	if(fGraphSignificanceVsLivetime!=0) delete fGraphSignificanceVsLivetime;
	fGraphSignificanceVsLivetime=0;
	if(fGraphOffEventsVsLivetime!=0) delete fGraphOffEventsVsLivetime;
	fGraphOffEventsVsLivetime=0;
	if(fGraphOnEventsVsLivetime!=0) delete fGraphOnEventsVsLivetime;
	fGraphOnEventsVsLivetime=0;
	if(fGraphExcessVsLivetime!=0) delete fGraphExcessVsLivetime;
	fGraphExcessVsLivetime=0;
	if(fGraphOffsetVsLivetime!=0) delete fGraphOffsetVsLivetime;
	fGraphOffsetVsLivetime=0;

}


/**
 * \brief Build the MC min and max values of eff, zen and off
 * used to get the effective areas, resolution and Biais
 * from the MC DSTs.
 */
void START::DataSummary::MakePairs(std::vector<Band> const &BandArray)
{
	/*
	MonteCarlo MC;

	double minoff(0), maxoff(0);
	double mineff(0), maxeff(0);
	double minzen(0), maxzen(0);

	for(unsigned int i(0); i<BandArray.size(); i++) {
		if(BandArray[i].GetKeepBand()==0) continue;
		fOffset.push_back(BandArray[i].GetOffset());
		fEfficiency.push_back(BandArray[i].GetEff());
		fZenith.push_back(BandArray[i].GetZenON());
	}

	// Find the min and max values of the bands parameters off, zen
	// and eff

	std::vector<double>::const_iterator it1, it2, it3, it4, it5, it6;
	it1 = min_element(fOffset.begin(), fOffset.end());
	minoff = *it1;
	it2 = max_element(fOffset.begin(), fOffset.end());
	maxoff = *it2;
	//std::cout << *it1 << std::endl;
	//std::cout << *it2 << std::endl;
	it3 = min_element(fEfficiency.begin(), fEfficiency.end());
	mineff = *it3;
	it4 = max_element(fEfficiency.begin(), fEfficiency.end());
	maxeff = *it4;

	it5 = min_element(fZenith.begin(), fZenith.end());
	minzen = *it5;
	it6 = max_element(fZenith.begin(), fZenith.end());
	maxzen = *it6;

	// Determine the MC values 

	double minoffMC(0), maxoffMC(0);
	double mineffMC(0), maxeffMC(0);
	double minzenMC(0), maxzenMC(0);

	// Determine efficiencies pair

	for(unsigned int i(0); i<MC.GetEfficiency().size()-1; i++) {
		if(mineff > MC.GetEfficiency()[i] && mineff < MC.GetEfficiency()[i+1]) mineffMC = MC.GetEfficiency()[i]; 
	}

	for(unsigned int i(0); i<MC.GetEfficiency().size()-1; i++) {
		if(maxeff > MC.GetEfficiency()[i] && maxeff < MC.GetEfficiency()[i+1]) maxeffMC = MC.GetEfficiency()[i+1]; 
	}

	for(unsigned int i(0); i<MC.GetEfficiency().size(); i++) DEBUG_OUT << MC.GetEfficiency()[i] << " ";
	DEBUG_OUT << std::endl;
	for(unsigned int i(0); i<fEfficiency.size(); i++) DEBUG_OUT << fEfficiency[i] << " ";
	DEBUG_OUT << std::endl;
	DEBUG_OUT << "mineff " << mineff << " maxeff " << maxeff << std::endl;  
	DEBUG_OUT << "mineffMC " << mineffMC << " maxeffMC " << maxeffMC << std::endl;


	fPairEfficiency = std::make_pair(mineffMC,maxeffMC);

	// Determine offset pair

	for(unsigned int i(0); i<MC.GetOffset().size()-1; i++) {
		if(minoff > MC.GetOffset()[i] && minoff < MC.GetOffset()[i+1]){ 
			minoffMC = MC.GetOffset()[i]; 
		}
	}

	for(unsigned int i(0); i<MC.GetOffset().size()-1; i++) {
		if(maxoff > MC.GetOffset()[i] && maxoff < MC.GetOffset()[i+1]){ 
			maxoffMC = MC.GetOffset()[i+1]; 
		}
	}
	for(unsigned int i(0); i<MC.GetOffset().size(); i++) DEBUG_OUT << MC.GetOffset()[i] << " ";
	DEBUG_OUT << std::endl;
	for(unsigned int i(0); i<fOffset.size(); i++) DEBUG_OUT << fOffset[i] << " ";
	DEBUG_OUT << std::endl;
	DEBUG_OUT << "minoff " << minoff << " maxoff " << maxoff << std::endl; 
	DEBUG_OUT << "minoffMC " << minoffMC << " maxoffMC " << maxoffMC << std::endl;


	fPairOffset = std::make_pair(minoffMC,maxoffMC);

	// Determine zenith pair

	for(unsigned int i(0); i<MC.GetZenith().size()-1; i++) {
		if(minzen > MC.GetZenith()[i] && minzen > MC.GetZenith()[i+1]) minzenMC = MC.GetZenith()[i]; 
	}

	for(unsigned int i(0); i<MC.GetZenith().size()-1; i++) {
		if(maxzen > MC.GetZenith()[i] && maxzen < MC.GetZenith()[i+1]) maxzenMC = MC.GetZenith()[i+1]; 
	}


	for(unsigned int i(0); i<MC.GetZenith().size(); i++) DEBUG_OUT << MC.GetZenith()[i] << " ";
	DEBUG_OUT << std::endl;
	for(unsigned int i(0); i<fZenith.size(); i++) DEBUG_OUT << fZenith[i] << " ";
	DEBUG_OUT << std::endl;
	DEBUG_OUT << "minzen " << minzen << " maxzen " << maxzen << std::endl;
	DEBUG_OUT << "minzenMC " << minzenMC << " maxzenMC " << maxzenMC << std::endl;


	fPairZenith = std::make_pair(minzenMC,maxzenMC);

	*/
}

/**
 * \brief Used to know what kink of telescopes is used in bands
 */
void START::DataSummary::DetermineTeltype(std::vector<Band> const &BandArray) {

	/*
	  MonteCarlo MC;
	  std::vector<TString> telMC = MC.GetTelType();
	  int u(0), v(0); // used to leave the loop
 
	  for(unsigned int i(0); i<BandArray.size(); i++) {

	  if(BandArray[i].GetKeepBand()==0) continue; // JLK we are only interested in selected bands


	  if(BandArray[i].GetTelCode()==30 && u==0) {
      fWhichTelType.push_back(telMC[0]);
      u++;
	  }
	  else if(BandArray[i].GetTelCode()!=30 && v==0) {
      v++;
      fWhichTelType.push_back(telMC[1]);
	  }
	  }
	*/
}

/**
 * \brief Print all the infos contained in the bands
 */
void START::DataSummary::PrintDataInfo(const std::vector<Band> &BandArray) const {

	std::vector<double> vemin, vemax, vethmc;

	INFO << "Bands carateristics : " << std::endl;
	for(std::vector<Band>::const_iterator band=BandArray.begin(); band!=BandArray.end(); ++band) {
		if(band->GetKeepBand()==0) continue; // we look for selected data only
		band->Print();
		vethmc.push_back(band->GetEthMC());
		for(std::vector<EnergyBin>::const_iterator bin=band->ebin.begin(); bin!=band->ebin.end(); ++bin) {
			if(bin->GetKeepBin()==0) continue;
			vemin.push_back(bin->GetEmin());
			vemax.push_back(bin->GetEmax());
		}
	}

	// energy treshold
	std::vector<double>::const_iterator it1, it2;
	it1 = min_element(vemin.begin(),vemin.end());
	it2 = max_element(vemax.begin(),vemax.end());
	double eth(0.), emax(0.);
	if(vemin.size()!=0 && vemax.size()!=0) {
		eth = *it1;
		emax = *it2;
	}

	std::vector<double>::const_iterator itethmc;
	itethmc = min_element(vethmc.begin(),vethmc.end());
	double ethmc(0.);
	if(vethmc.size()!=0) ethmc = *itethmc;

	Band SummaryBand(BandArray[0]);
	SummaryBand.ClearBandInfo();
	SummaryBand.AddInfoFromBands(BandArray,true);
  
	double sumlivetime = SummaryBand.GetLiveTime();

	double sumsignificance = SummaryBand.GetSignificanceTot(true);
	double sumnon = SummaryBand.GetNOnTot(true);
	double sumnoff = SummaryBand.GetNOffTot(true);
	double sumexcess = SummaryBand.GetExcessTot(true);

	double sumsignificanceall = SummaryBand.GetSignificanceTot(false);
	double sumnonall = SummaryBand.GetNOnTot(false);
	double sumnoffall = SummaryBand.GetNOffTot(false);
	double sumexcessall = SummaryBand.GetExcessTot(false);

	PrintSummaryBand(BandArray,false);
	PrintSummaryBand(BandArray,true);
	INFO << "Info on bands :" << std::endl;
	std::cout << "Number of bands = " << BandArray.size() << std::endl;
	std::cout << "Livetime = " << sumlivetime << " hours" << std::endl;
	std::cout << "MC energy treshold is " << ethmc << " TeV" << std::endl;
	std::cout << "Energy treshold is " << eth << " TeV" << std::endl;
	std::cout << "Energy max is " << emax << " TeV" << std::endl;
	INFO_OUT << "All infos :" << std::endl;
	std::cout << "Number of ON events = " << sumnonall << std::endl;
	std::cout << "Number of OFF events = " << sumnoffall << std::endl;
	std::cout << "Excess = " << sumexcessall << std::endl;
	std::cout << "Significance = " << sumsignificanceall << std::endl;
	INFO_OUT << "*Above threshold :" << std::endl;
	std::cout << "Number of ON events = " << sumnon << " (" 
			  << (sumnon/sumnonall)*100. << "%)" << std::endl;
	std::cout << "Number of OFF events = " << sumnoff << " (" 
			  << (sumnoffall/sumnoffall)*100. << "%)" << std::endl;
	std::cout << "Excess = " << sumexcess << std::endl;
	std::cout << "Significance = " << sumsignificance << std::endl;

}

void START::DataSummary::PrintDataInfoInFile(const std::vector<Band> &BandArray,std::ostream &os) const
{

	std::vector<double> vemin, vemax, vethmc;

	//INFO << "Bands carateristics : " << std::endl;
	for(std::vector<Band>::const_iterator band=BandArray.begin(); band!=BandArray.end(); ++band) {
		if(band->GetKeepBand()==0) continue; // we look for selected data only
		//band->Print();
		vethmc.push_back(band->GetEthMC());
		for(std::vector<EnergyBin>::const_iterator bin=band->ebin.begin(); bin!=band->ebin.end(); ++bin) {
			if(bin->GetKeepBin()==0) continue;
			vemin.push_back(bin->GetEmin());
			vemax.push_back(bin->GetEmax());
		}
	}

	// energy treshold
	std::vector<double>::const_iterator it1, it2;
	it1 = min_element(vemin.begin(),vemin.end());
	it2 = max_element(vemax.begin(),vemax.end());
	double eth(0.), emax(0.);
	if(vemin.size()!=0 && vemax.size()!=0) {
		eth = *it1;
		emax = *it2;
	}

	std::vector<double>::const_iterator itethmc;
	itethmc = min_element(vethmc.begin(),vethmc.end());
	double ethmc(0.);
	if(vethmc.size()!=0) ethmc = *itethmc;

	Band SummaryBand(BandArray[0]);
	SummaryBand.ClearBandInfo();
	SummaryBand.AddInfoFromBands(BandArray,true);
  
	double sumlivetime = SummaryBand.GetLiveTime();

	double sumsignificance = SummaryBand.GetSignificanceTot(true);
	double sumnon = SummaryBand.GetNOnTot(true);
	double sumnoff = SummaryBand.GetNOffTot(true);
	double sumexcess = SummaryBand.GetExcessTot(true);

	double sumsignificanceall = SummaryBand.GetSignificanceTot(false);
	double sumnonall = SummaryBand.GetNOnTot(false);
	double sumnoffall = SummaryBand.GetNOffTot(false);
	double sumexcessall = SummaryBand.GetExcessTot(false);

	os << "Info from bands : " << std::endl;
	os << "Number of bands = " << BandArray.size() << std::endl;
	os << "Livetime = " << sumlivetime << " hours" << std::endl;
	os << "MC energy treshold is " << ethmc << " TeV" << std::endl;
	os << "Energy treshold is " << eth << " TeV" << std::endl;
	os << "Energy max is " << emax << " TeV" << std::endl;
	os << "All infos :" << std::endl;
	os << "Number of ON events = " << sumnonall << std::endl;
	os << "Number of OFF events = " << sumnoffall << std::endl;
	os << "Excess = " << sumexcessall << std::endl;
	os << "Significance = " << sumsignificanceall << std::endl;
	os << "*Above threshold :" << std::endl;
	os << "Number of ON events = " << sumnon << " (" 
	   << (sumnon/sumnonall)*100. << "%)" << std::endl;
	os << "Number of OFF events = " << sumnoff << " (" 
	   << (sumnoffall/sumnoffall)*100. << "%)" << std::endl;
	os << "Excess = " << sumexcess << std::endl;
	os << "Significance = " << sumsignificance << std::endl;
}



/**
 * \brief Print infos of all bands in one
 */
void START::DataSummary::PrintSummaryBand(const std::vector<Band> &BandArray, bool UseEThreshold) const
{

	Band *SummaryBand = new Band(BandArray[0]);
	SummaryBand->ClearBandInfo();
	SummaryBand->AddInfoFromBands(BandArray,UseEThreshold);
	SummaryBand->SetKeepBand(1);

	if (UseEThreshold) {
		INFO << "SummaryBand (Above Ethreshold) : " << std::endl;
	}
	else {
		INFO << "SummaryBand (AllInfo) :" << std::endl;
	}
	SummaryBand->Print();

	delete SummaryBand;
	SummaryBand=NULL;

}

/**
 * \brief Init and fill summary histograms
 */
void START::DataSummary::InitAndFillSummaryHistograms(std::vector<Band> const &BandArray) {

	//CleanCanvas();
	//CleanTH1();
	//CleanGraph();

	/*
	  MonteCarlo MC;

	  unsigned int bineff = MC.GetEfficiency().size();
	  double bineffmin = MC.GetEfficiency().front();
	  double bineffmax = MC.GetEfficiency().back();
	  unsigned int binoff = MC.GetOffset().size();
	  double binoffmin = MC.GetOffset().front();
	  double binoffmax = MC.GetOffset().back();
	  unsigned int binzen = MC.GetZenith().size();
	  double binzenmin = MC.GetZenith().front();
	  double binzenmax = MC.GetZenith().back();
	*/

	// Band parameters distribution

	unsigned int nbin(50);

	fHistoEfficiency = new TH1D("HistoEfficiency","Efficiency distribution from selected data;Efficiency(%)",nbin,30,100);
	fHistoEfficiency->SetFillColor(kGreen-3);
	fHistoOffset = new TH1D("HistoOffset","Offset distribution from selected data;Offset",nbin,0.,3.);
	fHistoOffset->SetFillColor(kGreen-3);
	fHistoZenithON = new TH1D("HistoZenithON","Zenith ON distribution from selected data;ZenithON",nbin,0.,80.);
	fHistoZenithON->SetFillColor(kGreen-3);
	fHistoZenithOFF = new TH1D("HistoZenithOff","Zenith OFF distribution from selected data;ZenithOFF",nbin,0.,80.);
	fHistoZenithOFF->SetFillColor(kGreen-3);
	fHistoLiveTime = new TH1D("HistoLiveTime","Livetime distribution from selected data;Livetime(s)",nbin,0.,2000);
	fHistoLiveTime->SetFillColor(kGreen-3);
	fHistoAlpha = new TH1D("HistoAlpha","Alpha (Ton/Toff) distribution from selected data;Alpha",nbin,0.,1.);
	fHistoAlpha->SetFillColor(kGreen-3);
	fHistoTimeMJD = new TH1D("HistoTimeMJD","Time (MJD) distribution;Time(MJD)",nbin,50000,60000);
	fHistoTimeMJD->SetFillColor(kGreen-3);

	double emin(BandArray[0].ebin.front().GetEmin());
	double emax(BandArray[0].ebin.back().GetEmax());
	fHistoEnergyON = new TH1D("HistoEnergyON","Energy of ON events distribution from selected data;Energy ON evts",10000,emin,emax);
	fHistoEnergyON->SetFillColor(kGreen-3);
	fHistoEnergyOFF = new TH1D("HistoEnergyOFF","Energy of OFF events distribution from selected data;Energy Off evts",10000,emin,emax);
	fHistoEnergyOFF->SetFillColor(kGreen-3);

	for(std::vector<Band>::const_iterator band = BandArray.begin(); band!=BandArray.end(); ++band) {

		if(band->GetKeepBand()==0) continue;

		// Fill Efficiency
		fHistoEfficiency->Fill(band->GetEff());
		// Fill Offset
		fHistoOffset->Fill(band->GetOffset());
		// Fill Zenith on
		fHistoZenithON->Fill(band->GetZenON());
		// Fill Zenith off
		fHistoZenithOFF->Fill(band->GetZenOFF());
		// Fill LiveTime (seconds)
		fHistoLiveTime->Fill(band->GetLiveTime()*3600.);
		// Fill Alpha
		fHistoAlpha->Fill(band->GetAlphaRun());
		// Fill Time (MJD)
		fHistoTimeMJD->Fill(band->GetRunStartTime());

		for(std::vector<EnergyBin>::const_iterator bin = band->ebin.begin(); bin!=band->ebin.end(); bin++) {
     
			if(bin->GetKeepBin()==0) continue;

			unsigned int non = bin->GetOn();
			unsigned int noff = bin->GetOff();

			// Fill ON energy
			for(unsigned int ion(0); ion<non; ion++) fHistoEnergyON->Fill(bin->GetEmid());
			// Fill OFF energy
			for(unsigned int ioff(0); ioff<noff; ioff++) fHistoEnergyOFF->Fill(bin->GetEmid());

		}

	}

	// significance distribution
  
	fGraphSignificanceVsLivetime = new TGraph();
	fGraphSignificanceVsLivetime->SetName("GraphSignificanceVsLivetime");
	fGraphSignificanceVsLivetime->SetTitle("Li-Ma significance vs cumulated livetime");
	fGraphSignificanceVsLivetime->SetFillColor(0);
	fGraphSignificanceVsLivetime->SetMarkerColor(1);
	fGraphSignificanceVsLivetime->SetMarkerStyle(6);

	fGraphOffEventsVsLivetime = new TGraph();
	fGraphOffEventsVsLivetime->SetName("GraphOffEventsVsLivetime");
	fGraphOffEventsVsLivetime->SetTitle("#OFF events vs cumulated livetime");
	fGraphOffEventsVsLivetime->SetFillColor(0);
	fGraphOffEventsVsLivetime->SetMarkerColor(1);
	fGraphOffEventsVsLivetime->SetMarkerStyle(6);

	fGraphOnEventsVsLivetime = new TGraph();
	fGraphOnEventsVsLivetime->SetName("GraphOnEventsVsLivetime");
	fGraphOnEventsVsLivetime->SetTitle("#ON events vs cumulated livetime");
	fGraphOnEventsVsLivetime->SetFillColor(0);
	fGraphOnEventsVsLivetime->SetMarkerColor(1);
	fGraphOnEventsVsLivetime->SetMarkerStyle(6);

	fGraphExcessVsLivetime = new TGraph();
	fGraphExcessVsLivetime->SetName("GraphExcessVsLivetime");
	fGraphExcessVsLivetime->SetTitle("Excess vs cumulated livetime");
	fGraphExcessVsLivetime->SetFillColor(0);
	fGraphExcessVsLivetime->SetMarkerColor(1);
	fGraphExcessVsLivetime->SetMarkerStyle(6);

	fGraphOffsetVsLivetime = new TGraph();
	fGraphOffsetVsLivetime->SetName("GraphOffsetVsLivetime");
	fGraphOffsetVsLivetime->SetTitle("Offset vs cumulated livetime");
	fGraphOffsetVsLivetime->SetFillColor(0);
	fGraphOffsetVsLivetime->SetMarkerStyle(6);

	std::vector<Band> TmpBand;
	unsigned int irun(0);

	for(std::vector<Band>::const_iterator band = BandArray.begin(); band!=BandArray.end(); ++band) {

		if(band->GetKeepBand()==0) continue;
    
		TmpBand.push_back(*band);
    
		Band SummaryBand (TmpBand[0]);
		SummaryBand.ClearBandInfo();
		SummaryBand.AddInfoFromBands(TmpBand,true);

		double offset = band->GetOffset();

		double sumlivetime = SummaryBand.GetLiveTime();
		double sumsignificance = SummaryBand.GetSignificanceTot(true);
		double sumnon = SummaryBand.GetNOnTot(true);
		double sumnoff = SummaryBand.GetNOffTot(true);
		double sumexcess = SummaryBand.GetExcessTot(true);

		fGraphSignificanceVsLivetime->SetPoint(irun,sumlivetime,sumsignificance);
		fGraphOffEventsVsLivetime->SetPoint(irun,sumlivetime,sumnoff);
		fGraphOnEventsVsLivetime->SetPoint(irun,sumlivetime,sumnon);
		fGraphExcessVsLivetime->SetPoint(irun,sumlivetime,sumexcess);
		fGraphOffsetVsLivetime->SetPoint(irun,sumlivetime,offset);
		irun++;

	}

}

/**
 * \brief Draw summary histograms
 */
void START::DataSummary::DrawSummaryHistograms(std::vector<Band> const &BandArray) {

	InitAndFillSummaryHistograms(BandArray);

	fCanvasEfficiency = new TCanvas("CanvasEfficiency","CanvasEfficiency");
	fCanvasEfficiency->cd();
	fHistoEfficiency->Draw();

	fCanvasOffset = new TCanvas("CanvasOffset","CanvasOffset");
	fCanvasOffset->cd();
	fHistoOffset->Draw();

	fCanvasZenithON = new TCanvas("CanvasZenithON","CanvasZenithON");
	fCanvasZenithON->cd();
	fHistoZenithON->Draw();

	fCanvasZenithOFF = new TCanvas("CanvasZenithOFF","CanvasZenithOFF");
	fCanvasZenithOFF->cd();
	fHistoZenithOFF->Draw();

	fCanvasLiveTime = new TCanvas("CanvasLiveTime","CanvasLiveTime");
	fCanvasLiveTime->cd();
	fHistoLiveTime->Draw();

	fCanvasAlpha = new TCanvas("CanvasAlpha","CanvasAlpha");
	fCanvasAlpha->cd();
	fHistoAlpha->Draw();

	fCanvasTimeMJD = new TCanvas("CanvasTimeMJD","CanvasTimeMJD");
	fCanvasTimeMJD->cd();
	fHistoTimeMJD->Draw();

	// energy
	fCanvasEnergyON = new TCanvas("CanvasEnergyON","CanvasEnergyON");
	fCanvasEnergyON->cd();
	fCanvasEnergyON->SetGridx();
	fCanvasEnergyON->SetGridy();
	fCanvasEnergyON->SetLogx();
	fCanvasEnergyON->SetLogy();
	fHistoEnergyON->Draw();

	fCanvasEnergyOFF = new TCanvas("CanvasEnergyOFF","CanvasEnergyOFF");
	fCanvasEnergyOFF->cd();
	fCanvasEnergyOFF->SetGridx();
	fCanvasEnergyOFF->SetGridy();
	fCanvasEnergyOFF->SetLogx();
	fCanvasEnergyOFF->SetLogy();
	fHistoEnergyOFF->Draw();

	fCanvasSignificanceVsLivetime = new TCanvas("CanvasSignificanceVsLivetime","CanvasSignificanceVsLivetime");
	fCanvasSignificanceVsLivetime->cd();
	double x2limit(0.),y2limit(0.);
	fGraphSignificanceVsLivetime->GetPoint(fGraphSignificanceVsLivetime->GetN(),x2limit,y2limit);
	TF1 *sigmafit = new TF1("sigmafit","[0]*sqrt(x)",0,x2limit);
	sigmafit->SetParName(0,"norm");
	sigmafit->SetParameter(0,1.);
	sigmafit->SetLineColor(kGreen-3);
	sigmafit->SetLineStyle(2);
	sigmafit->SetNpx(1000);
	fCanvasSignificanceVsLivetime->cd();
	fGraphSignificanceVsLivetime->SetMarkerStyle(20);
	fGraphSignificanceVsLivetime->SetMarkerSize(0.6);
	fGraphSignificanceVsLivetime->Draw("AP");
	fGraphSignificanceVsLivetime->Fit(sigmafit,"R");
	double x1limit(0.),y1limit(0.);
	std::vector<double> minpoint;
	for(unsigned int ipoint(0); ipoint<(unsigned int)fGraphSignificanceVsLivetime->GetN(); ipoint++) {
		fGraphSignificanceVsLivetime->GetPoint(ipoint,x1limit,y1limit);
		minpoint.push_back(y1limit);
	}
	std::vector<double>::iterator it;
	it = min_element(minpoint.begin(),minpoint.end());
	y1limit=*it;
	if(y1limit<0) fGraphSignificanceVsLivetime->GetHistogram()->SetMinimum(y1limit-0.3);
	else fGraphSignificanceVsLivetime->GetHistogram()->SetMinimum(0.);
	//fGraphSignificanceVsLivetime->Fit(sigmafit,"QR");
	fGraphSignificanceVsLivetime->GetHistogram()->GetXaxis()->SetTitle("Livetime (hours)");
	fGraphSignificanceVsLivetime->GetHistogram()->GetYaxis()->SetTitle("Significance");
	//sigmafit->Draw("same");

	fCanvasOffEventsVsLivetime = new TCanvas("CanvasOffEventsVsLivetime","CanvasOffEventsVsLivetime");
	fCanvasOffEventsVsLivetime->cd();
	fGraphOffEventsVsLivetime->SetMarkerStyle(20);
	fGraphOffEventsVsLivetime->SetMarkerSize(0.6);
	fGraphOffEventsVsLivetime->Draw("AP");
	fGraphOffEventsVsLivetime->GetHistogram()->GetXaxis()->SetTitle("Livetime (hours)");
	fGraphOffEventsVsLivetime->GetHistogram()->GetYaxis()->SetTitle("#OFF events");
	fGraphOffEventsVsLivetime->GetHistogram()->GetYaxis()->SetLabelSize(0.03);

	fCanvasOnEventsVsLivetime = new TCanvas("CanvasOnEventsVsLivetime","CanvasOnEventsVsLivetime");
	fCanvasOnEventsVsLivetime->cd();
	fGraphOnEventsVsLivetime->SetMarkerStyle(20);
	fGraphOnEventsVsLivetime->SetMarkerSize(0.6);
	fGraphOnEventsVsLivetime->Draw("AP");
	fGraphOnEventsVsLivetime->GetHistogram()->GetXaxis()->SetTitle("Livetime (hours)");
	fGraphOnEventsVsLivetime->GetHistogram()->GetYaxis()->SetTitle("#ON events");
	fGraphOnEventsVsLivetime->GetHistogram()->GetYaxis()->SetLabelSize(0.03);

	fCanvasExcessVsLivetime = new TCanvas("CanvasExcessVsLivetime","CanvasExcessVsLivetime");
	fCanvasExcessVsLivetime->cd();
	fGraphExcessVsLivetime->SetMarkerStyle(20);
	fGraphExcessVsLivetime->SetMarkerSize(0.6);
	fGraphExcessVsLivetime->Draw("AP");
	fGraphExcessVsLivetime->GetHistogram()->GetXaxis()->SetTitle("Livetime (hours)");
	fGraphExcessVsLivetime->GetHistogram()->GetYaxis()->SetTitle("Excess");

	fCanvasOffsetVsLivetime = new TCanvas("CanvasOffsetVsLivetime","CanvasOffsetVsLivetime");
	fCanvasOffsetVsLivetime->cd();
	fGraphOffsetVsLivetime->SetMarkerStyle(20);
	fGraphOffsetVsLivetime->SetMarkerSize(0.6);
	fGraphOffsetVsLivetime->Draw("AP");
	fGraphOffsetVsLivetime->GetHistogram()->GetXaxis()->SetTitle("Livetime (hours)");
	fGraphOffsetVsLivetime->GetHistogram()->GetYaxis()->SetTitle("Offset (deg.)");
}
