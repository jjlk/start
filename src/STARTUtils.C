//STL
#include <iostream>
#include <string>
#include <sstream>
#include <map>

// ROOT
#include <TFile.h>
#include <TROOT.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TIterator.h>
#include <TList.h>

// START
#include "STARTUtils.hh"

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::STARTUtils);
#endif

/**
 * \brief Return Li and Ma significance
 *
 * \param non Number of ON events
 * \param noff Number of OFF events
 * \param alpha ratio of TON and TOFF
 */
double START::STARTUtils::LiMaSignificance(int non, int noff, double alpha) {

  double alpha_p1 = alpha+1.;
  double non_d  = (double)non;
  double noff_d = (double)noff;
  double ntot_d = non_d+noff_d;

  double part1 = 0.;
  if (non) {
    part1 = non_d * TMath::Log( (alpha_p1/alpha) * (non_d/ntot_d) );
  }

  double part2 = 0.;
  if (noff) {
    part2 = noff_d * TMath::Log(alpha_p1*noff_d/ntot_d);
  }
  
  double sigsq = 2.*(part1+part2);
  double sign = 1.;

  if ((non_d-alpha*noff_d)<0.) {
    sign = -1.;
  }
  return (sign*TMath::Sqrt(sigsq));
}

/**
 * \brief Return excess
 */
double START::STARTUtils::GetExcess(int non, int noff, double alpha)
{
  return (non-noff*alpha);
}


/**
 * \brief Function to get the expected background value.
 * Taken from the phd thesis of Piron, p. 204
 *
 * \param non Number of On events
 * \param noff Number of Off events
 * \param alpha Normalisation between on and off region
 * \param ns Expected signal
 */
double START::STARTUtils::GetExpectedOff(int non, int noff, double alpha, double ns) {
  double non_d = double(non);
  double noff_d = double(noff);
  double alpha_p1 = alpha+1.;
  double a = alpha*(non_d+noff_d)-(alpha_p1*ns);
  double b = a*a + 4.*alpha*alpha_p1*noff_d*ns;
  double pbar;
  
  if ( (noff==0) && (a<=0.) ) {
    pbar = 0.;
  }
  else {
    pbar = (a+TMath::Sqrt(b))/(2.*alpha*alpha_p1);
  }
  return pbar;
}

/**
 * \brief Function to get the expected On value.
 *
 * \param non Number of On events
 * \param noff Number of Off events
 * \param alpha Normalisation between on and off region
 * \param ns Expected signal
 */
double START::STARTUtils::GetExpectedOn(int non, int noff, double alpha, double ns)
{
  double pbar = GetExpectedOff(non,noff,alpha,ns);
  return (ns+alpha*pbar);
}

/**
 * \brief Return constante used during the minimization to avoid precisions issues
 */
double START::STARTUtils::GetNormalizationConstante() {
  return 1.e-11;
}

/**
 * \brief Return convertion factor from TeV to ergs
 */
double START::STARTUtils::GetTeVToErg() {
  return TMath::Qe()*1.e19; 
}

/**
 * \brief Return convertion factor from TeV to GeV
 */
double START::STARTUtils::GetTeVToGeV() {
  return TMath::Power(10.,3.); 
}

/**
 * \brief Return convertion factor from TeV to MeV
 */
double START::STARTUtils::GetTeVToMeV() {
  return TMath::Power(10.,6.); 
}

/**
 * \brief Return convertion factor from TeV to keV
 */
double START::STARTUtils::GetTeVTokeV() {
  return TMath::Power(10.,9.); 
}

/**
 * \brief Return convertion factor from TeV to Hz
 */
double START::STARTUtils::GetTeVToHz() {
  return (TMath::Qe()/TMath::H()*1.e12); 
}

/**
 * \brief Return MJD time from SASH seconds
 */
double START::STARTUtils::GetMJDFromSashSeconds(double seconds)
{
  double currentMJD = 40587.;
  currentMJD+=seconds/(24.*3600.);

  return currentMJD;
}

/**
 * \brief Return Sash seconds fomr MJD
 */
double START::STARTUtils::GetSashSecondsFromMJD(double mjd)
{
  double sashseconds = mjd-40587.;
  sashseconds=sashseconds*24.*3600.;

  return sashseconds;
}

/**
 * \brief Return run's root file name
 */
TString START::STARTUtils::GetRunROOTFileName(TString config, unsigned int run, unsigned int telcode)
{

  TString name("run");
  name+=run;
  name+="_";
  name+=config;
  if(telcode==14 || telcode==22 || telcode==26 || telcode==28) {
    name+="_3Tel";
  }
  name+=".root";

  return name;

}

/**
 * \brief Get base 10 exponent for a given number
 * @return Base10 exponent of the number (ex : number = 1.23e45 --> 45);
 */
int START::STARTUtils::GetIntegerExposantFromDecimalNumber(double number)
{
  return  ( (number!=0.) ? Int_t(TMath::Floor(TMath::Log10(TMath::Abs(number)))) : 0 );
}

/*
  int START::STARTUtils::GetIntegerExposantFromDecimalNumber(double number)
{

  int min(-20), max(20), exponent(0);
  
  if(number<0.) number=-number;

  bool IsNumberGreaterThan1(false);
  if(number>1.) IsNumberGreaterThan1=true;
  else if(number<1.) IsNumberGreaterThan1=false;
  else return 0;

  if(!IsNumberGreaterThan1) {
    //std::cout << "<1" << std::endl;
    exponent=max;
    double tmp(11.);
    while(tmp>=10.) {
      tmp=number*TMath::Power(10.,exponent);
      //std::cout << "exponent=" << exponent << " number=" << tmp << std::endl;
      exponent--;
      if(exponent<=min) break;
    };
    
    exponent++;

  }
  else {

    //std::cout << ">1" << std::endl;
    exponent=min;
    double tmp(-1.);
    while(tmp<1.) {
      tmp=number*TMath::Power(10.,exponent);
      //std::cout << "exponent=" << exponent << " number=" << tmp << std::endl;
      exponent++;
      if(exponent>=max) break;
    };

    exponent--;

  }
   
  //std::cout << "exponent=" << -exponent << std::endl;
  
  return (-exponent);
}
*/

std::vector<Int_t> START::STARTUtils::GetTelescopeListFromPattern(UInt_t TelPattern) {
  std::vector<Int_t> tellist;
  UInt_t TelNumber = 1;
  UInt_t TelId=0;
  
  while (TelNumber<TelPattern) {
    TelNumber = (TelNumber << 1);
    ++TelId;
    //std::cout << "TelNumber = " << TelNumber << " CT" << TelId << std::endl;
    if (TelPattern & TelNumber) {
      // This tel is in the Pattern
      tellist.push_back(TelId);
    }
  }

  return tellist;
}

/**
 * \brief Function to retrieve the contours from the root file at the end of the START fitting procedure
 * \param filename Path to the output ROOT file of START
 * \param cannam Name of the canvas that hold the contours
 *
 * \return map where the key is the number of sigma  of the contour, and the object is the contour in a TGraph
 **/
std::map<Int_t,TGraph*> START::STARTUtils::GetContourFromSTARTFile(std::string filename, std::string canname) {

  std::map<Int_t,TGraph*> map_cont;
  
  TFile *file = NULL;
  file = new TFile(filename.c_str());
  gROOT->cd();

  if (file==NULL) {
    std::cout << "Can't open the file : " << filename << " ! Please Check ! " << std::endl;
    return map_cont;
  }

  TCanvas *can = NULL;
  can = (TCanvas*)file->Get(canname.c_str());
  if (0==can) {
    std::cout << "The canvas :  " << canname << " has not been found in the file ! Please Check ! " << std::endl;
    return map_cont;
  }
  
  TList *list = NULL;
  list = can->GetListOfPrimitives();
  
  if (list==NULL) {
    std::cout << "Can't find the list of object inside the canvas! " << std::endl;
    return map_cont;
  }
  
  TIterator *iter = NULL;
  iter = list->MakeIterator();
  iter->Reset();
  
  for (int i=0;i<list->GetSize();++i) {
    TObject* obj = iter->Next();
    if (!obj) continue;
    // Contours are stored in a TMultiGraph:
    TMultiGraph *grObj = NULL;
    grObj = dynamic_cast<TMultiGraph*>(obj);
    
    if (grObj==NULL) {  
      continue;
    }
    else {
      //std::cout << "list entry : " << i << " We got a graph : " << std::endl;
      std::cout << "GraphName = " << grObj->GetName() << " GraphTitle = " << grObj->GetTitle() << std::endl;
      // Now Loop on the Graph List inside the TMultiGraph to retrieve the graphs;
      TList *grlist = grObj->GetListOfGraphs();
      std::cout << "grlist = " << grlist << std::endl;
      std::cout << "grlist->GetSize() = " << std::endl;
      std::cout << grlist->GetSize() << std::endl;
      TIterator *griter = 0;
      griter = grlist->MakeIterator();
      std::cout << "griter = " << griter << std::endl;
      griter->Reset();
      for (int j=0;j<grlist->GetSize();++j)
	{
	  std::cout << "ici" << std::endl;
	  TGraph *cur_gr = 0;
	  cur_gr = dynamic_cast< TGraph* > ( griter->Next() );
	  std::cout << "la" << std::endl;
	  if (!cur_gr) {
	    std::cout << "strange ! "<< std::endl;
	    continue;
	  }
	  //std::cout << " cur_gr->GetName() = " << cur_gr->GetName() << "  cur_gr->GetTitle() = " << cur_gr->GetTitle() << std::endl;
	  cur_gr->SetLineColor(1);
	  cur_gr->SetMarkerColor(1);
	  cur_gr->SetMarkerStyle(8);
	  // Get the name of the graph and retrieve the Sigma Number :
	  std::string sigmanumber = cur_gr->GetName();
	  //Contours2sigmaPWL_phi0_Gamma
	  std::string contname("Contours");
	  // Purpose : A START Graph contours name is like : Contours2sigmaPWL_phi0_Gamma.
	  // So I remove from "sigma" to the end.
	  // Then I remove Contours
	  // The only remaining string should be the number of sigma
	  std::cout << "sigmanumber = " << sigmanumber << std::endl;
	  sigmanumber = sigmanumber.erase( sigmanumber.find("sigma"), std::string::npos );
	  std::cout << "sigmanumber = " << sigmanumber << std::endl;
	  sigmanumber = sigmanumber.substr( sigmanumber.find(contname) + contname.length(), std::string::npos);
	  std::cout << "sigmanumber = " << sigmanumber << std::endl;
	  std::istringstream istream_nsigma( sigmanumber );
	  Int_t nsigma;
	  if (!(istream_nsigma >> nsigma)) {
	    std::cout << "Can't stream sigmanumber = " << sigmanumber << std::endl;
	    continue;
	  }
	  std::pair<std::map<Int_t,TGraph*>::iterator,Bool_t> ret;
	  ret = map_cont.insert( std::pair<Int_t,TGraph*>(nsigma,cur_gr) );
	  if (!ret.second) {
	    std::cout << "Can't Add the Graph : " << cur_gr->GetName() << " : An entry Already exist for the Significance : " << nsigma << std::endl;
	  }
	}
    }
  }
  return map_cont;
}


/**
 * \brief Return vector of pair containing start and end of full moon periods
 * from 1950 to 2060 =)
 *
 * Took from http://wise-obs.tau.ac.il/~eran/Wise/wise_calen.html
 */
std::vector<std::pair<double,double> > START::STARTUtils::GetNewMoonMJDPeriods()
{
  std::vector<std::pair<double,double> > MJDPeriod;
  MJDPeriod.push_back(std::make_pair(33285.33,33314.93));
  MJDPeriod.push_back(std::make_pair(33314.93,33344.44));
  MJDPeriod.push_back(std::make_pair(33344.44,33373.87));
  MJDPeriod.push_back(std::make_pair(33373.87,33403.22));
  MJDPeriod.push_back(std::make_pair(33403.22,33432.53));
  MJDPeriod.push_back(std::make_pair(33432.53,33461.83));
  MJDPeriod.push_back(std::make_pair(33461.83,33491.18));
  MJDPeriod.push_back(std::make_pair(33491.18,33520.62));
  MJDPeriod.push_back(std::make_pair(33520.62,33550.18));
  MJDPeriod.push_back(std::make_pair(33550.18,33579.87));
  MJDPeriod.push_back(std::make_pair(33579.87,33609.64));
  MJDPeriod.push_back(std::make_pair(33609.64,33639.43));
  MJDPeriod.push_back(std::make_pair(33639.43,33669.20));
  MJDPeriod.push_back(std::make_pair(33669.20,33698.88));
  MJDPeriod.push_back(std::make_pair(33698.88,33728.45));
  MJDPeriod.push_back(std::make_pair(33728.45,33757.90));
  MJDPeriod.push_back(std::make_pair(33757.90,33787.24));
  MJDPeriod.push_back(std::make_pair(33787.24,33816.52));
  MJDPeriod.push_back(std::make_pair(33816.52,33845.80));
  MJDPeriod.push_back(std::make_pair(33845.80,33875.12));
  MJDPeriod.push_back(std::make_pair(33875.12,33904.53));
  MJDPeriod.push_back(std::make_pair(33904.53,33934.03));
  MJDPeriod.push_back(std::make_pair(33934.03,33963.66));
  MJDPeriod.push_back(std::make_pair(33963.66,33993.40));
  MJDPeriod.push_back(std::make_pair(33993.40,34023.21));
  MJDPeriod.push_back(std::make_pair(34023.21,34053.02));
  MJDPeriod.push_back(std::make_pair(34053.02,34082.76));
  MJDPeriod.push_back(std::make_pair(34082.76,34112.37));
  MJDPeriod.push_back(std::make_pair(34112.37,34141.85));
  MJDPeriod.push_back(std::make_pair(34141.85,34171.21));
  MJDPeriod.push_back(std::make_pair(34171.21,34200.52));
  MJDPeriod.push_back(std::make_pair(34200.52,34229.82));
  MJDPeriod.push_back(std::make_pair(34229.82,34259.14));
  MJDPeriod.push_back(std::make_pair(34259.14,34288.51));
  MJDPeriod.push_back(std::make_pair(34288.51,34317.96));
  MJDPeriod.push_back(std::make_pair(34317.96,34347.53));
  MJDPeriod.push_back(std::make_pair(34347.53,34377.21));
  MJDPeriod.push_back(std::make_pair(34377.21,34406.99));
  MJDPeriod.push_back(std::make_pair(34406.99,34436.79));
  MJDPeriod.push_back(std::make_pair(34436.79,34466.54));
  MJDPeriod.push_back(std::make_pair(34466.54,34496.18));
  MJDPeriod.push_back(std::make_pair(34496.18,34525.71));
  MJDPeriod.push_back(std::make_pair(34525.71,34555.15));
  MJDPeriod.push_back(std::make_pair(34555.15,34584.52));
  MJDPeriod.push_back(std::make_pair(34584.52,34613.85));
  MJDPeriod.push_back(std::make_pair(34613.85,34643.18));
  MJDPeriod.push_back(std::make_pair(34643.18,34672.54));
  MJDPeriod.push_back(std::make_pair(34672.54,34701.97));
  MJDPeriod.push_back(std::make_pair(34701.97,34731.49));
  MJDPeriod.push_back(std::make_pair(34731.49,34761.11));
  MJDPeriod.push_back(std::make_pair(34761.11,34790.80));
  MJDPeriod.push_back(std::make_pair(34790.80,34820.53));
  MJDPeriod.push_back(std::make_pair(34820.53,34850.24));
  MJDPeriod.push_back(std::make_pair(34850.24,34879.91));
  MJDPeriod.push_back(std::make_pair(34879.91,34909.50));
  MJDPeriod.push_back(std::make_pair(34909.50,34939.02));
  MJDPeriod.push_back(std::make_pair(34939.02,34968.46));
  MJDPeriod.push_back(std::make_pair(34968.46,34997.85));
  MJDPeriod.push_back(std::make_pair(34997.85,35027.22));
  MJDPeriod.push_back(std::make_pair(35027.22,35056.60));
  MJDPeriod.push_back(std::make_pair(35056.60,35086.04));
  MJDPeriod.push_back(std::make_pair(35086.04,35115.53));
  MJDPeriod.push_back(std::make_pair(35115.53,35145.07));
  MJDPeriod.push_back(std::make_pair(35145.07,35174.65));
  MJDPeriod.push_back(std::make_pair(35174.65,35204.27));
  MJDPeriod.push_back(std::make_pair(35204.27,35233.93));
  MJDPeriod.push_back(std::make_pair(35233.93,35263.59));
  MJDPeriod.push_back(std::make_pair(35263.59,35293.23));
  MJDPeriod.push_back(std::make_pair(35293.23,35322.81));
  MJDPeriod.push_back(std::make_pair(35322.81,35352.33));
  MJDPeriod.push_back(std::make_pair(35352.33,35381.80));
  MJDPeriod.push_back(std::make_pair(35381.80,35411.25));
  MJDPeriod.push_back(std::make_pair(35411.25,35440.70));
  MJDPeriod.push_back(std::make_pair(35440.70,35470.16));
  MJDPeriod.push_back(std::make_pair(35470.16,35499.61));
  MJDPeriod.push_back(std::make_pair(35499.61,35529.07));
  MJDPeriod.push_back(std::make_pair(35529.07,35558.55));
  MJDPeriod.push_back(std::make_pair(35558.55,35588.07));
  MJDPeriod.push_back(std::make_pair(35588.07,35617.64));
  MJDPeriod.push_back(std::make_pair(35617.64,35647.26));
  MJDPeriod.push_back(std::make_pair(35647.26,35676.90));
  MJDPeriod.push_back(std::make_pair(35676.90,35706.53));
  MJDPeriod.push_back(std::make_pair(35706.53,35736.14));
  MJDPeriod.push_back(std::make_pair(35736.14,35765.73));
  MJDPeriod.push_back(std::make_pair(35765.73,35795.28));
  MJDPeriod.push_back(std::make_pair(35795.28,35824.80));
  MJDPeriod.push_back(std::make_pair(35824.80,35854.26));
  MJDPeriod.push_back(std::make_pair(35854.26,35883.69));
  MJDPeriod.push_back(std::make_pair(35883.69,35913.10));
  MJDPeriod.push_back(std::make_pair(35913.10,35942.51));
  MJDPeriod.push_back(std::make_pair(35942.51,35971.94));
  MJDPeriod.push_back(std::make_pair(35971.94,36001.42));
  MJDPeriod.push_back(std::make_pair(36001.42,36030.95));
  MJDPeriod.push_back(std::make_pair(36030.95,36060.55));
  MJDPeriod.push_back(std::make_pair(36060.55,36090.21));
  MJDPeriod.push_back(std::make_pair(36090.21,36119.91));
  MJDPeriod.push_back(std::make_pair(36119.91,36149.61));
  MJDPeriod.push_back(std::make_pair(36149.61,36179.26));
  MJDPeriod.push_back(std::make_pair(36179.26,36208.84));
  MJDPeriod.push_back(std::make_pair(36208.84,36238.34));
  MJDPeriod.push_back(std::make_pair(36238.34,36267.77));
  MJDPeriod.push_back(std::make_pair(36267.77,36297.16));
  MJDPeriod.push_back(std::make_pair(36297.16,36326.52));
  MJDPeriod.push_back(std::make_pair(36326.52,36355.87));
  MJDPeriod.push_back(std::make_pair(36355.87,36385.25));
  MJDPeriod.push_back(std::make_pair(36385.25,36414.70));
  MJDPeriod.push_back(std::make_pair(36414.70,36444.25));
  MJDPeriod.push_back(std::make_pair(36444.25,36473.91));
  MJDPeriod.push_back(std::make_pair(36473.91,36503.65));
  MJDPeriod.push_back(std::make_pair(36503.65,36533.43));
  MJDPeriod.push_back(std::make_pair(36533.43,36563.16));
  MJDPeriod.push_back(std::make_pair(36563.16,36592.82));
  MJDPeriod.push_back(std::make_pair(36592.82,36622.37));
  MJDPeriod.push_back(std::make_pair(36622.37,36651.84));
  MJDPeriod.push_back(std::make_pair(36651.84,36681.22));
  MJDPeriod.push_back(std::make_pair(36681.22,36710.54));
  MJDPeriod.push_back(std::make_pair(36710.54,36739.83));
  MJDPeriod.push_back(std::make_pair(36739.83,36769.15));
  MJDPeriod.push_back(std::make_pair(36769.15,36798.53));
  MJDPeriod.push_back(std::make_pair(36798.53,36828.04));
  MJDPeriod.push_back(std::make_pair(36828.04,36857.67));
  MJDPeriod.push_back(std::make_pair(36857.67,36887.40));
  MJDPeriod.push_back(std::make_pair(36887.40,36917.20));
  MJDPeriod.push_back(std::make_pair(36917.20,36946.99));
  MJDPeriod.push_back(std::make_pair(36946.99,36976.73));
  MJDPeriod.push_back(std::make_pair(36976.73,37006.35));
  MJDPeriod.push_back(std::make_pair(37006.35,37035.85));
  MJDPeriod.push_back(std::make_pair(37035.85,37065.24));
  MJDPeriod.push_back(std::make_pair(37065.24,37094.54));
  MJDPeriod.push_back(std::make_pair(37094.54,37123.82));
  MJDPeriod.push_back(std::make_pair(37123.82,37153.11));
  MJDPeriod.push_back(std::make_pair(37153.11,37182.47));
  MJDPeriod.push_back(std::make_pair(37182.47,37211.93));
  MJDPeriod.push_back(std::make_pair(37211.93,37241.50));
  MJDPeriod.push_back(std::make_pair(37241.50,37271.18));
  MJDPeriod.push_back(std::make_pair(37271.18,37300.96));
  MJDPeriod.push_back(std::make_pair(37300.96,37330.78));
  MJDPeriod.push_back(std::make_pair(37330.78,37360.57));
  MJDPeriod.push_back(std::make_pair(37360.57,37390.24));
  MJDPeriod.push_back(std::make_pair(37390.24,37419.78));
  MJDPeriod.push_back(std::make_pair(37419.78,37449.19));
  MJDPeriod.push_back(std::make_pair(37449.19,37478.53));
  MJDPeriod.push_back(std::make_pair(37478.53,37507.83));
  MJDPeriod.push_back(std::make_pair(37507.83,37537.13));
  MJDPeriod.push_back(std::make_pair(37537.13,37566.48));
  MJDPeriod.push_back(std::make_pair(37566.48,37595.90));
  MJDPeriod.push_back(std::make_pair(37595.90,37625.41));
  MJDPeriod.push_back(std::make_pair(37625.41,37655.03));
  MJDPeriod.push_back(std::make_pair(37655.03,37684.76));
  MJDPeriod.push_back(std::make_pair(37684.76,37714.55));
  MJDPeriod.push_back(std::make_pair(37714.55,37744.33));
  MJDPeriod.push_back(std::make_pair(37744.33,37774.02));
  MJDPeriod.push_back(std::make_pair(37774.02,37803.61));
  MJDPeriod.push_back(std::make_pair(37803.61,37833.09));
  MJDPeriod.push_back(std::make_pair(37833.09,37862.49));
  MJDPeriod.push_back(std::make_pair(37862.49,37891.84));
  MJDPeriod.push_back(std::make_pair(37891.84,37921.18));
  MJDPeriod.push_back(std::make_pair(37921.18,37950.52));
  MJDPeriod.push_back(std::make_pair(37950.52,37979.92));
  MJDPeriod.push_back(std::make_pair(37979.92,38009.39));
  MJDPeriod.push_back(std::make_pair(38009.39,38038.96));
  MJDPeriod.push_back(std::make_pair(38038.96,38068.62));
  MJDPeriod.push_back(std::make_pair(38068.62,38098.33));
  MJDPeriod.push_back(std::make_pair(38098.33,38128.04));
  MJDPeriod.push_back(std::make_pair(38128.04,38157.72));
  MJDPeriod.push_back(std::make_pair(38157.72,38187.35));
  MJDPeriod.push_back(std::make_pair(38187.35,38216.91));
  MJDPeriod.push_back(std::make_pair(38216.91,38246.40));
  MJDPeriod.push_back(std::make_pair(38246.40,38275.82));
  MJDPeriod.push_back(std::make_pair(38275.82,38305.20));
  MJDPeriod.push_back(std::make_pair(38305.20,38334.58));
  MJDPeriod.push_back(std::make_pair(38334.58,38364.00));
  MJDPeriod.push_back(std::make_pair(38364.00,38393.46));
  MJDPeriod.push_back(std::make_pair(38393.46,38422.97));
  MJDPeriod.push_back(std::make_pair(38422.97,38452.53));
  MJDPeriod.push_back(std::make_pair(38452.53,38482.12));
  MJDPeriod.push_back(std::make_pair(38482.12,38511.74));
  MJDPeriod.push_back(std::make_pair(38511.74,38541.40));
  MJDPeriod.push_back(std::make_pair(38541.40,38571.05));
  MJDPeriod.push_back(std::make_pair(38571.05,38600.67));
  MJDPeriod.push_back(std::make_pair(38600.67,38630.23));
  MJDPeriod.push_back(std::make_pair(38630.23,38659.73));
  MJDPeriod.push_back(std::make_pair(38659.73,38689.20));
  MJDPeriod.push_back(std::make_pair(38689.20,38718.66));
  MJDPeriod.push_back(std::make_pair(38718.66,38748.11));
  MJDPeriod.push_back(std::make_pair(38748.11,38777.57));
  MJDPeriod.push_back(std::make_pair(38777.57,38807.02));
  MJDPeriod.push_back(std::make_pair(38807.02,38836.48));
  MJDPeriod.push_back(std::make_pair(38836.48,38865.96));
  MJDPeriod.push_back(std::make_pair(38865.96,38895.50));
  MJDPeriod.push_back(std::make_pair(38895.50,38925.08));
  MJDPeriod.push_back(std::make_pair(38925.08,38954.71));
  MJDPeriod.push_back(std::make_pair(38954.71,38984.35));
  MJDPeriod.push_back(std::make_pair(38984.35,39013.98));
  MJDPeriod.push_back(std::make_pair(39013.98,39043.59));
  MJDPeriod.push_back(std::make_pair(39043.59,39073.18));
  MJDPeriod.push_back(std::make_pair(39073.18,39102.72));
  MJDPeriod.push_back(std::make_pair(39102.72,39132.22));
  MJDPeriod.push_back(std::make_pair(39132.22,39161.67));
  MJDPeriod.push_back(std::make_pair(39161.67,39191.07));
  MJDPeriod.push_back(std::make_pair(39191.07,39220.47));
  MJDPeriod.push_back(std::make_pair(39220.47,39249.88));
  MJDPeriod.push_back(std::make_pair(39249.88,39279.32));
  MJDPeriod.push_back(std::make_pair(39279.32,39308.82));
  MJDPeriod.push_back(std::make_pair(39308.82,39338.38));
  MJDPeriod.push_back(std::make_pair(39338.38,39368.01));
  MJDPeriod.push_back(std::make_pair(39368.01,39397.70));
  MJDPeriod.push_back(std::make_pair(39397.70,39427.42));
  MJDPeriod.push_back(std::make_pair(39427.42,39457.11));
  MJDPeriod.push_back(std::make_pair(39457.11,39486.74));
  MJDPeriod.push_back(std::make_pair(39486.74,39516.28));
  MJDPeriod.push_back(std::make_pair(39516.28,39545.74));
  MJDPeriod.push_back(std::make_pair(39545.74,39575.14));
  MJDPeriod.push_back(std::make_pair(39575.14,39604.50));
  MJDPeriod.push_back(std::make_pair(39604.50,39633.85));
  MJDPeriod.push_back(std::make_pair(39633.85,39663.21));
  MJDPeriod.push_back(std::make_pair(39663.21,39692.61));
  MJDPeriod.push_back(std::make_pair(39692.61,39722.10));
  MJDPeriod.push_back(std::make_pair(39722.10,39751.71));
  MJDPeriod.push_back(std::make_pair(39751.71,39781.43));
  MJDPeriod.push_back(std::make_pair(39781.43,39811.20));
  MJDPeriod.push_back(std::make_pair(39811.20,39840.97));
  MJDPeriod.push_back(std::make_pair(39840.97,39870.68));
  MJDPeriod.push_back(std::make_pair(39870.68,39900.28));
  MJDPeriod.push_back(std::make_pair(39900.28,39929.79));
  MJDPeriod.push_back(std::make_pair(39929.79,39959.20));
  MJDPeriod.push_back(std::make_pair(39959.20,39988.55));
  MJDPeriod.push_back(std::make_pair(39988.55,40017.84));
  MJDPeriod.push_back(std::make_pair(40017.84,40047.14));
  MJDPeriod.push_back(std::make_pair(40047.14,40076.48));
  MJDPeriod.push_back(std::make_pair(40076.48,40105.92));
  MJDPeriod.push_back(std::make_pair(40105.92,40135.49));
  MJDPeriod.push_back(std::make_pair(40135.49,40165.18));
  MJDPeriod.push_back(std::make_pair(40165.18,40194.96));
  MJDPeriod.push_back(std::make_pair(40194.96,40224.77));
  MJDPeriod.push_back(std::make_pair(40224.77,40254.54));
  MJDPeriod.push_back(std::make_pair(40254.54,40284.22));
  MJDPeriod.push_back(std::make_pair(40284.22,40313.78));
  MJDPeriod.push_back(std::make_pair(40313.78,40343.22));
  MJDPeriod.push_back(std::make_pair(40343.22,40372.56));
  MJDPeriod.push_back(std::make_pair(40372.56,40401.84));
  MJDPeriod.push_back(std::make_pair(40401.84,40431.12));
  MJDPeriod.push_back(std::make_pair(40431.12,40460.44));
  MJDPeriod.push_back(std::make_pair(40460.44,40489.85));
  MJDPeriod.push_back(std::make_pair(40489.85,40519.36));
  MJDPeriod.push_back(std::make_pair(40519.36,40549.00));
  MJDPeriod.push_back(std::make_pair(40549.00,40578.73));
  MJDPeriod.push_back(std::make_pair(40578.73,40608.54));
  MJDPeriod.push_back(std::make_pair(40608.54,40638.35));
  MJDPeriod.push_back(std::make_pair(40638.35,40668.08));
  MJDPeriod.push_back(std::make_pair(40668.08,40697.68));
  MJDPeriod.push_back(std::make_pair(40697.68,40727.15));
  MJDPeriod.push_back(std::make_pair(40727.15,40756.52));
  MJDPeriod.push_back(std::make_pair(40756.52,40785.83));
  MJDPeriod.push_back(std::make_pair(40785.83,40815.14));
  MJDPeriod.push_back(std::make_pair(40815.14,40844.47));
  MJDPeriod.push_back(std::make_pair(40844.47,40873.85));
  MJDPeriod.push_back(std::make_pair(40873.85,40903.31));
  MJDPeriod.push_back(std::make_pair(40903.31,40932.88));
  MJDPeriod.push_back(std::make_pair(40932.88,40962.56));
  MJDPeriod.push_back(std::make_pair(40962.56,40992.32));
  MJDPeriod.push_back(std::make_pair(40992.32,41022.11));
  MJDPeriod.push_back(std::make_pair(41022.11,41051.84));
  MJDPeriod.push_back(std::make_pair(41051.84,41081.48));
  MJDPeriod.push_back(std::make_pair(41081.48,41111.00));
  MJDPeriod.push_back(std::make_pair(41111.00,41140.44));
  MJDPeriod.push_back(std::make_pair(41140.44,41169.82));
  MJDPeriod.push_back(std::make_pair(41169.82,41199.17));
  MJDPeriod.push_back(std::make_pair(41199.17,41228.51));
  MJDPeriod.push_back(std::make_pair(41228.51,41257.89));
  MJDPeriod.push_back(std::make_pair(41257.89,41287.33));
  MJDPeriod.push_back(std::make_pair(41287.33,41316.85));
  MJDPeriod.push_back(std::make_pair(41316.85,41346.46));
  MJDPeriod.push_back(std::make_pair(41346.46,41376.13));
  MJDPeriod.push_back(std::make_pair(41376.13,41405.84));
  MJDPeriod.push_back(std::make_pair(41405.84,41435.53));
  MJDPeriod.push_back(std::make_pair(41435.53,41465.19));
  MJDPeriod.push_back(std::make_pair(41465.19,41494.78));
  MJDPeriod.push_back(std::make_pair(41494.78,41524.31));
  MJDPeriod.push_back(std::make_pair(41524.31,41553.77));
  MJDPeriod.push_back(std::make_pair(41553.77,41583.17));
  MJDPeriod.push_back(std::make_pair(41583.17,41612.56));
  MJDPeriod.push_back(std::make_pair(41612.56,41641.96));
  MJDPeriod.push_back(std::make_pair(41641.96,41671.41));
  MJDPeriod.push_back(std::make_pair(41671.41,41700.89));
  MJDPeriod.push_back(std::make_pair(41700.89,41730.42));
  MJDPeriod.push_back(std::make_pair(41730.42,41759.98));
  MJDPeriod.push_back(std::make_pair(41759.98,41789.58));
  MJDPeriod.push_back(std::make_pair(41789.58,41819.21));
  MJDPeriod.push_back(std::make_pair(41819.21,41848.86));
  MJDPeriod.push_back(std::make_pair(41848.86,41878.50));
  MJDPeriod.push_back(std::make_pair(41878.50,41908.10));
  MJDPeriod.push_back(std::make_pair(41908.10,41937.64));
  MJDPeriod.push_back(std::make_pair(41937.64,41967.13));
  MJDPeriod.push_back(std::make_pair(41967.13,41996.60));
  MJDPeriod.push_back(std::make_pair(41996.60,42026.07));
  MJDPeriod.push_back(std::make_pair(42026.07,42055.53));
  MJDPeriod.push_back(std::make_pair(42055.53,42084.98));
  MJDPeriod.push_back(std::make_pair(42084.98,42114.42));
  MJDPeriod.push_back(std::make_pair(42114.42,42143.88));
  MJDPeriod.push_back(std::make_pair(42143.88,42173.37));
  MJDPeriod.push_back(std::make_pair(42173.37,42202.92));
  MJDPeriod.push_back(std::make_pair(42202.92,42232.53));
  MJDPeriod.push_back(std::make_pair(42232.53,42262.16));
  MJDPeriod.push_back(std::make_pair(42262.16,42291.81));
  MJDPeriod.push_back(std::make_pair(42291.81,42321.44));
  MJDPeriod.push_back(std::make_pair(42321.44,42351.06));
  MJDPeriod.push_back(std::make_pair(42351.06,42380.63));
  MJDPeriod.push_back(std::make_pair(42380.63,42410.16));
  MJDPeriod.push_back(std::make_pair(42410.16,42439.63));
  MJDPeriod.push_back(std::make_pair(42439.63,42469.05));
  MJDPeriod.push_back(std::make_pair(42469.05,42498.44));
  MJDPeriod.push_back(std::make_pair(42498.44,42527.83));
  MJDPeriod.push_back(std::make_pair(42527.83,42557.24));
  MJDPeriod.push_back(std::make_pair(42557.24,42586.70));
  MJDPeriod.push_back(std::make_pair(42586.70,42616.23));
  MJDPeriod.push_back(std::make_pair(42616.23,42645.82));
  MJDPeriod.push_back(std::make_pair(42645.82,42675.49));
  MJDPeriod.push_back(std::make_pair(42675.49,42705.21));
  MJDPeriod.push_back(std::make_pair(42705.21,42734.94));
  MJDPeriod.push_back(std::make_pair(42734.94,42764.61));
  MJDPeriod.push_back(std::make_pair(42764.61,42794.20));
  MJDPeriod.push_back(std::make_pair(42794.20,42823.70));
  MJDPeriod.push_back(std::make_pair(42823.70,42853.12));
  MJDPeriod.push_back(std::make_pair(42853.12,42882.49));
  MJDPeriod.push_back(std::make_pair(42882.49,42911.84));
  MJDPeriod.push_back(std::make_pair(42911.84,42941.18));
  MJDPeriod.push_back(std::make_pair(42941.18,42970.55));
  MJDPeriod.push_back(std::make_pair(42970.55,42999.99));
  MJDPeriod.push_back(std::make_pair(42999.99,43029.54));
  MJDPeriod.push_back(std::make_pair(43029.54,43059.21));
  MJDPeriod.push_back(std::make_pair(43059.21,43088.97));
  MJDPeriod.push_back(std::make_pair(43088.97,43118.76));
  MJDPeriod.push_back(std::make_pair(43118.76,43148.51));
  MJDPeriod.push_back(std::make_pair(43148.51,43178.17));
  MJDPeriod.push_back(std::make_pair(43178.17,43207.72));
  MJDPeriod.push_back(std::make_pair(43207.72,43237.17));
  MJDPeriod.push_back(std::make_pair(43237.17,43266.54));
  MJDPeriod.push_back(std::make_pair(43266.54,43295.86));
  MJDPeriod.push_back(std::make_pair(43295.86,43325.14));
  MJDPeriod.push_back(std::make_pair(43325.14,43354.45));
  MJDPeriod.push_back(std::make_pair(43354.45,43383.84));
  MJDPeriod.push_back(std::make_pair(43383.84,43413.34));
  MJDPeriod.push_back(std::make_pair(43413.34,43442.98));
  MJDPeriod.push_back(std::make_pair(43442.98,43472.73));
  MJDPeriod.push_back(std::make_pair(43472.73,43502.53));
  MJDPeriod.push_back(std::make_pair(43502.53,43532.33));
  MJDPeriod.push_back(std::make_pair(43532.33,43562.06));
  MJDPeriod.push_back(std::make_pair(43562.06,43591.68));
  MJDPeriod.push_back(std::make_pair(43591.68,43621.18));
  MJDPeriod.push_back(std::make_pair(43621.18,43650.55));
  MJDPeriod.push_back(std::make_pair(43650.55,43679.86));
  MJDPeriod.push_back(std::make_pair(43679.86,43709.13));
  MJDPeriod.push_back(std::make_pair(43709.13,43738.43));
  MJDPeriod.push_back(std::make_pair(43738.43,43767.79));
  MJDPeriod.push_back(std::make_pair(43767.79,43797.26));
  MJDPeriod.push_back(std::make_pair(43797.26,43826.83));
  MJDPeriod.push_back(std::make_pair(43826.83,43856.52));
  MJDPeriod.push_back(std::make_pair(43856.52,43886.30));
  MJDPeriod.push_back(std::make_pair(43886.30,43916.11));
  MJDPeriod.push_back(std::make_pair(43916.11,43945.89));
  MJDPeriod.push_back(std::make_pair(43945.89,43975.55));
  MJDPeriod.push_back(std::make_pair(43975.55,44005.09));
  MJDPeriod.push_back(std::make_pair(44005.09,44034.50));
  MJDPeriod.push_back(std::make_pair(44034.50,44063.83));
  MJDPeriod.push_back(std::make_pair(44063.83,44093.14));
  MJDPeriod.push_back(std::make_pair(44093.14,44122.46));
  MJDPeriod.push_back(std::make_pair(44122.46,44151.82));
  MJDPeriod.push_back(std::make_pair(44151.82,44181.24));
  MJDPeriod.push_back(std::make_pair(44181.24,44210.76));
  MJDPeriod.push_back(std::make_pair(44210.76,44240.38));
  MJDPeriod.push_back(std::make_pair(44240.38,44270.10));
  MJDPeriod.push_back(std::make_pair(44270.10,44299.88));
  MJDPeriod.push_back(std::make_pair(44299.88,44329.64));
  MJDPeriod.push_back(std::make_pair(44329.64,44359.32));
  MJDPeriod.push_back(std::make_pair(44359.32,44388.90));
  MJDPeriod.push_back(std::make_pair(44388.90,44418.38));
  MJDPeriod.push_back(std::make_pair(44418.38,44447.79));
  MJDPeriod.push_back(std::make_pair(44447.79,44477.16));
  MJDPeriod.push_back(std::make_pair(44477.16,44506.51));
  MJDPeriod.push_back(std::make_pair(44506.51,44535.87));
  MJDPeriod.push_back(std::make_pair(44535.87,44565.28));
  MJDPeriod.push_back(std::make_pair(44565.28,44594.76));
  MJDPeriod.push_back(std::make_pair(44594.76,44624.32));
  MJDPeriod.push_back(std::make_pair(44624.32,44653.96));
  MJDPeriod.push_back(std::make_pair(44653.96,44683.64));
  MJDPeriod.push_back(std::make_pair(44683.64,44713.33));
  MJDPeriod.push_back(std::make_pair(44713.33,44743.00));
  MJDPeriod.push_back(std::make_pair(44743.00,44772.63));
  MJDPeriod.push_back(std::make_pair(44772.63,44802.20));
  MJDPeriod.push_back(std::make_pair(44802.20,44831.69));
  MJDPeriod.push_back(std::make_pair(44831.69,44861.13));
  MJDPeriod.push_back(std::make_pair(44861.13,44890.54));
  MJDPeriod.push_back(std::make_pair(44890.54,44919.94));
  MJDPeriod.push_back(std::make_pair(44919.94,44949.36));
  MJDPeriod.push_back(std::make_pair(44949.36,44978.83));
  MJDPeriod.push_back(std::make_pair(44978.83,45008.33));
  MJDPeriod.push_back(std::make_pair(45008.33,45037.87));
  MJDPeriod.push_back(std::make_pair(45037.87,45067.43));
  MJDPeriod.push_back(std::make_pair(45067.43,45097.03));
  MJDPeriod.push_back(std::make_pair(45097.03,45126.67));
  MJDPeriod.push_back(std::make_pair(45126.67,45156.31));
  MJDPeriod.push_back(std::make_pair(45156.31,45185.94));
  MJDPeriod.push_back(std::make_pair(45185.94,45215.52));
  MJDPeriod.push_back(std::make_pair(45215.52,45245.05));
  MJDPeriod.push_back(std::make_pair(45245.05,45274.54));
  MJDPeriod.push_back(std::make_pair(45274.54,45304.02));
  MJDPeriod.push_back(std::make_pair(45304.02,45333.48));
  MJDPeriod.push_back(std::make_pair(45333.48,45362.94));
  MJDPeriod.push_back(std::make_pair(45362.94,45392.37));
  MJDPeriod.push_back(std::make_pair(45392.37,45421.81));
  MJDPeriod.push_back(std::make_pair(45421.81,45451.27));
  MJDPeriod.push_back(std::make_pair(45451.27,45480.78));
  MJDPeriod.push_back(std::make_pair(45480.78,45510.36));
  MJDPeriod.push_back(std::make_pair(45510.36,45539.98));
  MJDPeriod.push_back(std::make_pair(45539.98,45569.63));
  MJDPeriod.push_back(std::make_pair(45569.63,45599.28));
  MJDPeriod.push_back(std::make_pair(45599.28,45628.91));
  MJDPeriod.push_back(std::make_pair(45628.91,45658.52));
  MJDPeriod.push_back(std::make_pair(45658.52,45688.08));
  MJDPeriod.push_back(std::make_pair(45688.08,45717.59));
  MJDPeriod.push_back(std::make_pair(45717.59,45747.03));
  MJDPeriod.push_back(std::make_pair(45747.03,45776.42));
  MJDPeriod.push_back(std::make_pair(45776.42,45805.80));
  MJDPeriod.push_back(std::make_pair(45805.80,45835.19));
  MJDPeriod.push_back(std::make_pair(45835.19,45864.61));
  MJDPeriod.push_back(std::make_pair(45864.61,45894.10));
  MJDPeriod.push_back(std::make_pair(45894.10,45923.66));
  MJDPeriod.push_back(std::make_pair(45923.66,45953.29));
  MJDPeriod.push_back(std::make_pair(45953.29,45983.00));
  MJDPeriod.push_back(std::make_pair(45983.00,46012.74));
  MJDPeriod.push_back(std::make_pair(46012.74,46042.46));
  MJDPeriod.push_back(std::make_pair(46042.46,46072.10));
  MJDPeriod.push_back(std::make_pair(46072.10,46101.64));
  MJDPeriod.push_back(std::make_pair(46101.64,46131.09));
  MJDPeriod.push_back(std::make_pair(46131.09,46160.48));
  MJDPeriod.push_back(std::make_pair(46160.48,46189.83));
  MJDPeriod.push_back(std::make_pair(46189.83,46219.16));
  MJDPeriod.push_back(std::make_pair(46219.16,46248.51));
  MJDPeriod.push_back(std::make_pair(46248.51,46277.90));
  MJDPeriod.push_back(std::make_pair(46277.90,46307.39));
  MJDPeriod.push_back(std::make_pair(46307.39,46337.01));
  MJDPeriod.push_back(std::make_pair(46337.01,46366.74));
  MJDPeriod.push_back(std::make_pair(46366.74,46396.53));
  MJDPeriod.push_back(std::make_pair(46396.53,46426.31));
  MJDPeriod.push_back(std::make_pair(46426.31,46456.02));
  MJDPeriod.push_back(std::make_pair(46456.02,46485.63));
  MJDPeriod.push_back(std::make_pair(46485.63,46515.13));
  MJDPeriod.push_back(std::make_pair(46515.13,46544.53));
  MJDPeriod.push_back(std::make_pair(46544.53,46573.86));
  MJDPeriod.push_back(std::make_pair(46573.86,46603.15));
  MJDPeriod.push_back(std::make_pair(46603.15,46632.44));
  MJDPeriod.push_back(std::make_pair(46632.44,46661.79));
  MJDPeriod.push_back(std::make_pair(46661.79,46691.23));
  MJDPeriod.push_back(std::make_pair(46691.23,46720.81));
  MJDPeriod.push_back(std::make_pair(46720.81,46750.51));
  MJDPeriod.push_back(std::make_pair(46750.51,46780.30));
  MJDPeriod.push_back(std::make_pair(46780.30,46810.11));
  MJDPeriod.push_back(std::make_pair(46810.11,46839.87));
  MJDPeriod.push_back(std::make_pair(46839.87,46869.55));
  MJDPeriod.push_back(std::make_pair(46869.55,46899.11));
  MJDPeriod.push_back(std::make_pair(46899.11,46928.54));
  MJDPeriod.push_back(std::make_pair(46928.54,46957.87));
  MJDPeriod.push_back(std::make_pair(46957.87,46987.15));
  MJDPeriod.push_back(std::make_pair(46987.15,47016.43));
  MJDPeriod.push_back(std::make_pair(47016.43,47045.76));
  MJDPeriod.push_back(std::make_pair(47045.76,47075.18));
  MJDPeriod.push_back(std::make_pair(47075.18,47104.70));
  MJDPeriod.push_back(std::make_pair(47104.70,47134.33));
  MJDPeriod.push_back(std::make_pair(47134.33,47164.07));
  MJDPeriod.push_back(std::make_pair(47164.07,47193.87));
  MJDPeriod.push_back(std::make_pair(47193.87,47223.67));
  MJDPeriod.push_back(std::make_pair(47223.67,47253.39));
  MJDPeriod.push_back(std::make_pair(47253.39,47282.99));
  MJDPeriod.push_back(std::make_pair(47282.99,47312.46));
  MJDPeriod.push_back(std::make_pair(47312.46,47341.82));
  MJDPeriod.push_back(std::make_pair(47341.82,47371.14));
  MJDPeriod.push_back(std::make_pair(47371.14,47400.46));
  MJDPeriod.push_back(std::make_pair(47400.46,47429.80));
  MJDPeriod.push_back(std::make_pair(47429.80,47459.19));
  MJDPeriod.push_back(std::make_pair(47459.19,47488.66));
  MJDPeriod.push_back(std::make_pair(47488.66,47518.23));
  MJDPeriod.push_back(std::make_pair(47518.23,47547.90));
  MJDPeriod.push_back(std::make_pair(47547.90,47577.65));
  MJDPeriod.push_back(std::make_pair(47577.65,47607.42));
  MJDPeriod.push_back(std::make_pair(47607.42,47637.14));
  MJDPeriod.push_back(std::make_pair(47637.14,47666.76));
  MJDPeriod.push_back(std::make_pair(47666.76,47696.29));
  MJDPeriod.push_back(std::make_pair(47696.29,47725.74));
  MJDPeriod.push_back(std::make_pair(47725.74,47755.13));
  MJDPeriod.push_back(std::make_pair(47755.13,47784.49));
  MJDPeriod.push_back(std::make_pair(47784.49,47813.86));
  MJDPeriod.push_back(std::make_pair(47813.86,47843.24));
  MJDPeriod.push_back(std::make_pair(47843.24,47872.69));
  MJDPeriod.push_back(std::make_pair(47872.69,47902.21));
  MJDPeriod.push_back(std::make_pair(47902.21,47931.80));
  MJDPeriod.push_back(std::make_pair(47931.80,47961.46));
  MJDPeriod.push_back(std::make_pair(47961.46,47991.14));
  MJDPeriod.push_back(std::make_pair(47991.14,48020.81));
  MJDPeriod.push_back(std::make_pair(48020.81,48050.46));
  MJDPeriod.push_back(std::make_pair(48050.46,48080.06));
  MJDPeriod.push_back(std::make_pair(48080.06,48109.60));
  MJDPeriod.push_back(std::make_pair(48109.60,48139.07));
  MJDPeriod.push_back(std::make_pair(48139.07,48168.50));
  MJDPeriod.push_back(std::make_pair(48168.50,48197.91));
  MJDPeriod.push_back(std::make_pair(48197.91,48227.33));
  MJDPeriod.push_back(std::make_pair(48227.33,48256.78));
  MJDPeriod.push_back(std::make_pair(48256.78,48286.26));
  MJDPeriod.push_back(std::make_pair(48286.26,48315.77));
  MJDPeriod.push_back(std::make_pair(48315.77,48345.30));
  MJDPeriod.push_back(std::make_pair(48345.30,48374.87));
  MJDPeriod.push_back(std::make_pair(48374.87,48404.48));
  MJDPeriod.push_back(std::make_pair(48404.48,48434.12));
  MJDPeriod.push_back(std::make_pair(48434.12,48463.77));
  MJDPeriod.push_back(std::make_pair(48463.77,48493.38));
  MJDPeriod.push_back(std::make_pair(48493.38,48522.95));
  MJDPeriod.push_back(std::make_pair(48522.95,48552.47));
  MJDPeriod.push_back(std::make_pair(48552.47,48581.96));
  MJDPeriod.push_back(std::make_pair(48581.96,48611.43));
  MJDPeriod.push_back(std::make_pair(48611.43,48640.90));
  MJDPeriod.push_back(std::make_pair(48640.90,48670.34));
  MJDPeriod.push_back(std::make_pair(48670.34,48699.76));
  MJDPeriod.push_back(std::make_pair(48699.76,48729.20));
  MJDPeriod.push_back(std::make_pair(48729.20,48758.67));
  MJDPeriod.push_back(std::make_pair(48758.67,48788.20));
  MJDPeriod.push_back(std::make_pair(48788.20,48817.80));
  MJDPeriod.push_back(std::make_pair(48817.80,48847.44));
  MJDPeriod.push_back(std::make_pair(48847.44,48877.10));
  MJDPeriod.push_back(std::make_pair(48877.10,48906.75));
  MJDPeriod.push_back(std::make_pair(48906.75,48936.39));
  MJDPeriod.push_back(std::make_pair(48936.39,48965.99));
  MJDPeriod.push_back(std::make_pair(48965.99,48995.53));
  MJDPeriod.push_back(std::make_pair(48995.53,49025.00));
  MJDPeriod.push_back(std::make_pair(49025.00,49054.41));
  MJDPeriod.push_back(std::make_pair(49054.41,49083.78));
  MJDPeriod.push_back(std::make_pair(49083.78,49113.15));
  MJDPeriod.push_back(std::make_pair(49113.15,49142.54));
  MJDPeriod.push_back(std::make_pair(49142.54,49171.99));
  MJDPeriod.push_back(std::make_pair(49171.99,49201.51));
  MJDPeriod.push_back(std::make_pair(49201.51,49231.11));
  MJDPeriod.push_back(std::make_pair(49231.11,49260.79));
  MJDPeriod.push_back(std::make_pair(49260.79,49290.53));
  MJDPeriod.push_back(std::make_pair(49290.53,49320.27));
  MJDPeriod.push_back(std::make_pair(49320.27,49349.96));
  MJDPeriod.push_back(std::make_pair(49349.96,49379.56));
  MJDPeriod.push_back(std::make_pair(49379.56,49409.05));
  MJDPeriod.push_back(std::make_pair(49409.05,49438.47));
  MJDPeriod.push_back(std::make_pair(49438.47,49467.82));
  MJDPeriod.push_back(std::make_pair(49467.82,49497.15));
  MJDPeriod.push_back(std::make_pair(49497.15,49526.48));
  MJDPeriod.push_back(std::make_pair(49526.48,49555.84));
  MJDPeriod.push_back(std::make_pair(49555.84,49585.28));
  MJDPeriod.push_back(std::make_pair(49585.28,49614.83));
  MJDPeriod.push_back(std::make_pair(49614.83,49644.51));
  MJDPeriod.push_back(std::make_pair(49644.51,49674.29));
  MJDPeriod.push_back(std::make_pair(49674.29,49704.10));
  MJDPeriod.push_back(std::make_pair(49704.10,49733.85));
  MJDPeriod.push_back(std::make_pair(49733.85,49763.51));
  MJDPeriod.push_back(std::make_pair(49763.51,49793.06));
  MJDPeriod.push_back(std::make_pair(49793.06,49822.51));
  MJDPeriod.push_back(std::make_pair(49822.51,49851.87));
  MJDPeriod.push_back(std::make_pair(49851.87,49881.17));
  MJDPeriod.push_back(std::make_pair(49881.17,49910.45));
  MJDPeriod.push_back(std::make_pair(49910.45,49939.76));
  MJDPeriod.push_back(std::make_pair(49939.76,49969.15));
  MJDPeriod.push_back(std::make_pair(49969.15,49998.66));
  MJDPeriod.push_back(std::make_pair(49998.66,50028.31));
  MJDPeriod.push_back(std::make_pair(50028.31,50058.06));
  MJDPeriod.push_back(std::make_pair(50058.06,50087.87));
  MJDPeriod.push_back(std::make_pair(50087.87,50117.67));
  MJDPeriod.push_back(std::make_pair(50117.67,50147.39));
  MJDPeriod.push_back(std::make_pair(50147.39,50177.01));
  MJDPeriod.push_back(std::make_pair(50177.01,50206.49));
  MJDPeriod.push_back(std::make_pair(50206.49,50235.87));
  MJDPeriod.push_back(std::make_pair(50235.87,50265.17));
  MJDPeriod.push_back(std::make_pair(50265.17,50294.44));
  MJDPeriod.push_back(std::make_pair(50294.44,50323.75));
  MJDPeriod.push_back(std::make_pair(50323.75,50353.12));
  MJDPeriod.push_back(std::make_pair(50353.12,50382.59));
  MJDPeriod.push_back(std::make_pair(50382.59,50412.17));
  MJDPeriod.push_back(std::make_pair(50412.17,50441.86));
  MJDPeriod.push_back(std::make_pair(50441.86,50471.63));
  MJDPeriod.push_back(std::make_pair(50471.63,50501.44));
  MJDPeriod.push_back(std::make_pair(50501.44,50531.20));
  MJDPeriod.push_back(std::make_pair(50531.20,50560.86));
  MJDPeriod.push_back(std::make_pair(50560.86,50590.39));
  MJDPeriod.push_back(std::make_pair(50590.39,50619.80));
  MJDPeriod.push_back(std::make_pair(50619.80,50649.14));
  MJDPeriod.push_back(std::make_pair(50649.14,50678.46));
  MJDPeriod.push_back(std::make_pair(50678.46,50707.79));
  MJDPeriod.push_back(std::make_pair(50707.79,50737.16));
  MJDPeriod.push_back(std::make_pair(50737.16,50766.59));
  MJDPeriod.push_back(std::make_pair(50766.59,50796.11));
  MJDPeriod.push_back(std::make_pair(50796.11,50825.73));
  MJDPeriod.push_back(std::make_pair(50825.73,50855.43));
  MJDPeriod.push_back(std::make_pair(50855.43,50885.19));
  MJDPeriod.push_back(std::make_pair(50885.19,50914.93));
  MJDPeriod.push_back(std::make_pair(50914.93,50944.60));
  MJDPeriod.push_back(std::make_pair(50944.60,50974.18));
  MJDPeriod.push_back(std::make_pair(50974.18,51003.67));
  MJDPeriod.push_back(std::make_pair(51003.67,51033.09));
  MJDPeriod.push_back(std::make_pair(51033.09,51062.47));
  MJDPeriod.push_back(std::make_pair(51062.47,51091.84));
  MJDPeriod.push_back(std::make_pair(51091.84,51121.22));
  MJDPeriod.push_back(std::make_pair(51121.22,51150.64));
  MJDPeriod.push_back(std::make_pair(51150.64,51180.12));
  MJDPeriod.push_back(std::make_pair(51180.12,51209.67));
  MJDPeriod.push_back(std::make_pair(51209.67,51239.29));
  MJDPeriod.push_back(std::make_pair(51239.29,51268.95));
  MJDPeriod.push_back(std::make_pair(51268.95,51298.62));
  MJDPeriod.push_back(std::make_pair(51298.62,51328.28));
  MJDPeriod.push_back(std::make_pair(51328.28,51357.90));
  MJDPeriod.push_back(std::make_pair(51357.90,51387.48));
  MJDPeriod.push_back(std::make_pair(51387.48,51416.99));
  MJDPeriod.push_back(std::make_pair(51416.99,51446.45));
  MJDPeriod.push_back(std::make_pair(51446.45,51475.88));
  MJDPeriod.push_back(std::make_pair(51475.88,51505.30));
  MJDPeriod.push_back(std::make_pair(51505.30,51534.73));
  MJDPeriod.push_back(std::make_pair(51534.73,51564.20));
  MJDPeriod.push_back(std::make_pair(51564.20,51593.69));
  MJDPeriod.push_back(std::make_pair(51593.69,51623.20));
  MJDPeriod.push_back(std::make_pair(51623.20,51652.74));
  MJDPeriod.push_back(std::make_pair(51652.74,51682.32));
  MJDPeriod.push_back(std::make_pair(51682.32,51711.94));
  MJDPeriod.push_back(std::make_pair(51711.94,51741.58));
  MJDPeriod.push_back(std::make_pair(51741.58,51771.22));
  MJDPeriod.push_back(std::make_pair(51771.22,51800.82));
  MJDPeriod.push_back(std::make_pair(51800.82,51830.37));
  MJDPeriod.push_back(std::make_pair(51830.37,51859.89));
  MJDPeriod.push_back(std::make_pair(51859.89,51889.38));
  MJDPeriod.push_back(std::make_pair(51889.38,51918.85));
  MJDPeriod.push_back(std::make_pair(51918.85,51948.30));
  MJDPeriod.push_back(std::make_pair(51948.30,51977.73));
  MJDPeriod.push_back(std::make_pair(51977.73,52007.14));
  MJDPeriod.push_back(std::make_pair(52007.14,52036.58));
  MJDPeriod.push_back(std::make_pair(52036.58,52066.07));
  MJDPeriod.push_back(std::make_pair(52066.07,52095.63));
  MJDPeriod.push_back(std::make_pair(52095.63,52125.25));
  MJDPeriod.push_back(std::make_pair(52125.25,52154.91));
  MJDPeriod.push_back(std::make_pair(52154.91,52184.58));
  MJDPeriod.push_back(std::make_pair(52184.58,52214.24));
  MJDPeriod.push_back(std::make_pair(52214.24,52243.87));
  MJDPeriod.push_back(std::make_pair(52243.87,52273.45));
  MJDPeriod.push_back(std::make_pair(52273.45,52302.95));
  MJDPeriod.push_back(std::make_pair(52302.95,52332.39));
  MJDPeriod.push_back(std::make_pair(52332.39,52361.77));
  MJDPeriod.push_back(std::make_pair(52361.77,52391.13));
  MJDPeriod.push_back(std::make_pair(52391.13,52420.49));
  MJDPeriod.push_back(std::make_pair(52420.49,52449.91));
  MJDPeriod.push_back(std::make_pair(52449.91,52479.38));
  MJDPeriod.push_back(std::make_pair(52479.38,52508.94));
  MJDPeriod.push_back(std::make_pair(52508.94,52538.58));
  MJDPeriod.push_back(std::make_pair(52538.58,52568.31));
  MJDPeriod.push_back(std::make_pair(52568.31,52598.07));
  MJDPeriod.push_back(std::make_pair(52598.07,52627.80));
  MJDPeriod.push_back(std::make_pair(52627.80,52657.45));
  MJDPeriod.push_back(std::make_pair(52657.45,52687.00));
  MJDPeriod.push_back(std::make_pair(52687.00,52716.44));
  MJDPeriod.push_back(std::make_pair(52716.44,52745.82));
  MJDPeriod.push_back(std::make_pair(52745.82,52775.15));
  MJDPeriod.push_back(std::make_pair(52775.15,52804.47));
  MJDPeriod.push_back(std::make_pair(52804.47,52833.81));
  MJDPeriod.push_back(std::make_pair(52833.81,52863.20));
  MJDPeriod.push_back(std::make_pair(52863.20,52892.69));
  MJDPeriod.push_back(std::make_pair(52892.69,52922.31));
  MJDPeriod.push_back(std::make_pair(52922.31,52952.05));
  MJDPeriod.push_back(std::make_pair(52952.05,52981.86));
  MJDPeriod.push_back(std::make_pair(52981.86,53011.65));
  MJDPeriod.push_back(std::make_pair(53011.65,53041.37));
  MJDPeriod.push_back(std::make_pair(53041.37,53070.97));
  MJDPeriod.push_back(std::make_pair(53070.97,53100.46));
  MJDPeriod.push_back(std::make_pair(53100.46,53129.86));
  MJDPeriod.push_back(std::make_pair(53129.86,53159.18));
  MJDPeriod.push_back(std::make_pair(53159.18,53188.47));
  MJDPeriod.push_back(std::make_pair(53188.47,53217.75));
  MJDPeriod.push_back(std::make_pair(53217.75,53247.10));
  MJDPeriod.push_back(std::make_pair(53247.10,53276.55));
  MJDPeriod.push_back(std::make_pair(53276.55,53306.13));
  MJDPeriod.push_back(std::make_pair(53306.13,53335.84));
  MJDPeriod.push_back(std::make_pair(53335.84,53365.63));
  MJDPeriod.push_back(std::make_pair(53365.63,53395.44));
  MJDPeriod.push_back(std::make_pair(53395.44,53425.21));
  MJDPeriod.push_back(std::make_pair(53425.21,53454.88));
  MJDPeriod.push_back(std::make_pair(53454.88,53484.42));
  MJDPeriod.push_back(std::make_pair(53484.42,53513.85));
  MJDPeriod.push_back(std::make_pair(53513.85,53543.18));
  MJDPeriod.push_back(std::make_pair(53543.18,53572.46));
  MJDPeriod.push_back(std::make_pair(53572.46,53601.75));
  MJDPeriod.push_back(std::make_pair(53601.75,53631.08));
  MJDPeriod.push_back(std::make_pair(53631.08,53660.51));
  MJDPeriod.push_back(std::make_pair(53660.51,53690.04));
  MJDPeriod.push_back(std::make_pair(53690.04,53719.68));
  MJDPeriod.push_back(std::make_pair(53719.68,53749.41));
  MJDPeriod.push_back(std::make_pair(53749.41,53779.20));
  MJDPeriod.push_back(std::make_pair(53779.20,53808.98));
  MJDPeriod.push_back(std::make_pair(53808.98,53838.70));
  MJDPeriod.push_back(std::make_pair(53838.70,53868.29));
  MJDPeriod.push_back(std::make_pair(53868.29,53897.75));
  MJDPeriod.push_back(std::make_pair(53897.75,53927.13));
  MJDPeriod.push_back(std::make_pair(53927.13,53956.46));
  MJDPeriod.push_back(std::make_pair(53956.46,53985.78));
  MJDPeriod.push_back(std::make_pair(53985.78,54015.13));
  MJDPeriod.push_back(std::make_pair(54015.13,54044.54));
  MJDPeriod.push_back(std::make_pair(54044.54,54074.02));
  MJDPeriod.push_back(std::make_pair(54074.02,54103.58));
  MJDPeriod.push_back(std::make_pair(54103.58,54133.24));
  MJDPeriod.push_back(std::make_pair(54133.24,54162.97));
  MJDPeriod.push_back(std::make_pair(54162.97,54192.72));
  MJDPeriod.push_back(std::make_pair(54192.72,54222.42));
  MJDPeriod.push_back(std::make_pair(54222.42,54252.05));
  MJDPeriod.push_back(std::make_pair(54252.05,54281.58));
  MJDPeriod.push_back(std::make_pair(54281.58,54311.03));
  MJDPeriod.push_back(std::make_pair(54311.03,54340.44));
  MJDPeriod.push_back(std::make_pair(54340.44,54369.82));
  MJDPeriod.push_back(std::make_pair(54369.82,54399.20));
  MJDPeriod.push_back(std::make_pair(54399.20,54428.60));
  MJDPeriod.push_back(std::make_pair(54428.60,54458.05));
  MJDPeriod.push_back(std::make_pair(54458.05,54487.57));
  MJDPeriod.push_back(std::make_pair(54487.57,54517.15));
  MJDPeriod.push_back(std::make_pair(54517.15,54546.78));
  MJDPeriod.push_back(std::make_pair(54546.78,54576.43));
  MJDPeriod.push_back(std::make_pair(54576.43,54606.09));
  MJDPeriod.push_back(std::make_pair(54606.09,54635.73));
  MJDPeriod.push_back(std::make_pair(54635.73,54665.33));
  MJDPeriod.push_back(std::make_pair(54665.33,54694.89));
  MJDPeriod.push_back(std::make_pair(54694.89,54724.39));
  MJDPeriod.push_back(std::make_pair(54724.39,54753.84));
  MJDPeriod.push_back(std::make_pair(54753.84,54783.26));
  MJDPeriod.push_back(std::make_pair(54783.26,54812.69));
  MJDPeriod.push_back(std::make_pair(54812.69,54842.14));
  MJDPeriod.push_back(std::make_pair(54842.14,54871.62));
  MJDPeriod.push_back(std::make_pair(54871.62,54901.11));
  MJDPeriod.push_back(std::make_pair(54901.11,54930.62));
  MJDPeriod.push_back(std::make_pair(54930.62,54960.17));
  MJDPeriod.push_back(std::make_pair(54960.17,54989.76));
  MJDPeriod.push_back(std::make_pair(54989.76,55019.39));
  MJDPeriod.push_back(std::make_pair(55019.39,55049.04));
  MJDPeriod.push_back(std::make_pair(55049.04,55078.67));
  MJDPeriod.push_back(std::make_pair(55078.67,55108.26));
  MJDPeriod.push_back(std::make_pair(55108.26,55137.80));
  MJDPeriod.push_back(std::make_pair(55137.80,55167.31));
  MJDPeriod.push_back(std::make_pair(55167.31,55196.80));
  MJDPeriod.push_back(std::make_pair(55196.80,55226.26));
  MJDPeriod.push_back(std::make_pair(55226.26,55255.69));
  MJDPeriod.push_back(std::make_pair(55255.69,55285.10));
  MJDPeriod.push_back(std::make_pair(55285.10,55314.51));
  MJDPeriod.push_back(std::make_pair(55314.51,55343.96));
  MJDPeriod.push_back(std::make_pair(55343.96,55373.48));
  MJDPeriod.push_back(std::make_pair(55373.48,55403.07));
  MJDPeriod.push_back(std::make_pair(55403.07,55432.71));
  MJDPeriod.push_back(std::make_pair(55432.71,55462.39));
  MJDPeriod.push_back(std::make_pair(55462.39,55492.07));
  MJDPeriod.push_back(std::make_pair(55492.07,55521.73));
  MJDPeriod.push_back(std::make_pair(55521.73,55551.34));
  MJDPeriod.push_back(std::make_pair(55551.34,55580.89));
  MJDPeriod.push_back(std::make_pair(55580.89,55610.36));
  MJDPeriod.push_back(std::make_pair(55610.36,55639.76));
  MJDPeriod.push_back(std::make_pair(55639.76,55669.11));
  MJDPeriod.push_back(std::make_pair(55669.11,55698.46));
  MJDPeriod.push_back(std::make_pair(55698.46,55727.84));
  MJDPeriod.push_back(std::make_pair(55727.84,55757.28));
  MJDPeriod.push_back(std::make_pair(55757.28,55786.79));
  MJDPeriod.push_back(std::make_pair(55786.79,55816.39));
  MJDPeriod.push_back(std::make_pair(55816.39,55846.09));
  MJDPeriod.push_back(std::make_pair(55846.09,55875.85));
  MJDPeriod.push_back(std::make_pair(55875.85,55905.61));
  MJDPeriod.push_back(std::make_pair(55905.61,55935.31));
  MJDPeriod.push_back(std::make_pair(55935.31,55964.91));
  MJDPeriod.push_back(std::make_pair(55964.91,55994.40));
  MJDPeriod.push_back(std::make_pair(55994.40,56023.81));
  MJDPeriod.push_back(std::make_pair(56023.81,56053.15));
  MJDPeriod.push_back(std::make_pair(56053.15,56082.47));
  MJDPeriod.push_back(std::make_pair(56082.47,56111.79));
  MJDPeriod.push_back(std::make_pair(56111.79,56141.14));
  MJDPeriod.push_back(std::make_pair(56141.14,56170.58));
  MJDPeriod.push_back(std::make_pair(56170.58,56200.14));
  MJDPeriod.push_back(std::make_pair(56200.14,56229.83));
  MJDPeriod.push_back(std::make_pair(56229.83,56259.62));
  MJDPeriod.push_back(std::make_pair(56259.62,56289.43));
  MJDPeriod.push_back(std::make_pair(56289.43,56319.19));
  MJDPeriod.push_back(std::make_pair(56319.19,56348.85));
  MJDPeriod.push_back(std::make_pair(56348.85,56378.40));
  MJDPeriod.push_back(std::make_pair(56378.40,56407.83));
  MJDPeriod.push_back(std::make_pair(56407.83,56437.19));
  MJDPeriod.push_back(std::make_pair(56437.19,56466.48));
  MJDPeriod.push_back(std::make_pair(56466.48,56495.76));
  MJDPeriod.push_back(std::make_pair(56495.76,56525.07));
  MJDPeriod.push_back(std::make_pair(56525.07,56554.47));
  MJDPeriod.push_back(std::make_pair(56554.47,56583.98));
  MJDPeriod.push_back(std::make_pair(56583.98,56613.64));
  MJDPeriod.push_back(std::make_pair(56613.64,56643.40));
  MJDPeriod.push_back(std::make_pair(56643.40,56673.20));
  MJDPeriod.push_back(std::make_pair(56673.20,56703.00));
  MJDPeriod.push_back(std::make_pair(56703.00,56732.72));
  MJDPeriod.push_back(std::make_pair(56732.72,56762.32));
  MJDPeriod.push_back(std::make_pair(56762.32,56791.80));
  MJDPeriod.push_back(std::make_pair(56791.80,56821.18));
  MJDPeriod.push_back(std::make_pair(56821.18,56850.48));
  MJDPeriod.push_back(std::make_pair(56850.48,56879.76));
  MJDPeriod.push_back(std::make_pair(56879.76,56909.07));
  MJDPeriod.push_back(std::make_pair(56909.07,56938.45));
  MJDPeriod.push_back(std::make_pair(56938.45,56967.93));
  MJDPeriod.push_back(std::make_pair(56967.93,56997.52));
  MJDPeriod.push_back(std::make_pair(56997.52,57027.20));
  MJDPeriod.push_back(std::make_pair(57027.20,57056.97));
  MJDPeriod.push_back(std::make_pair(57056.97,57086.75));
  MJDPeriod.push_back(std::make_pair(57086.75,57116.51));
  MJDPeriod.push_back(std::make_pair(57116.51,57146.16));
  MJDPeriod.push_back(std::make_pair(57146.16,57175.68));
  MJDPeriod.push_back(std::make_pair(57175.68,57205.10));
  MJDPeriod.push_back(std::make_pair(57205.10,57234.45));
  MJDPeriod.push_back(std::make_pair(57234.45,57263.78));
  MJDPeriod.push_back(std::make_pair(57263.78,57293.12));
  MJDPeriod.push_back(std::make_pair(57293.12,57322.50));
  MJDPeriod.push_back(std::make_pair(57322.50,57351.95));
  MJDPeriod.push_back(std::make_pair(57351.95,57381.47));
  MJDPeriod.push_back(std::make_pair(57381.47,57411.07));
  MJDPeriod.push_back(std::make_pair(57411.07,57440.76));
  MJDPeriod.push_back(std::make_pair(57440.76,57470.50));
  MJDPeriod.push_back(std::make_pair(57470.50,57500.23));
  MJDPeriod.push_back(std::make_pair(57500.23,57529.89));
  MJDPeriod.push_back(std::make_pair(57529.89,57559.46));
  MJDPeriod.push_back(std::make_pair(57559.46,57588.96));
  MJDPeriod.push_back(std::make_pair(57588.96,57618.40));
  MJDPeriod.push_back(std::make_pair(57618.40,57647.80));
  MJDPeriod.push_back(std::make_pair(57647.80,57677.18));
  MJDPeriod.push_back(std::make_pair(57677.18,57706.58));
  MJDPeriod.push_back(std::make_pair(57706.58,57736.00));
  MJDPeriod.push_back(std::make_pair(57736.00,57765.48));
  MJDPeriod.push_back(std::make_pair(57765.48,57795.02));
  MJDPeriod.push_back(std::make_pair(57795.02,57824.62));
  MJDPeriod.push_back(std::make_pair(57824.62,57854.26));
  MJDPeriod.push_back(std::make_pair(57854.26,57883.91));
  MJDPeriod.push_back(std::make_pair(57883.91,57913.55));
  MJDPeriod.push_back(std::make_pair(57913.55,57943.17));
  MJDPeriod.push_back(std::make_pair(57943.17,57972.76));
  MJDPeriod.push_back(std::make_pair(57972.76,58002.30));
  MJDPeriod.push_back(std::make_pair(58002.30,58031.78));
  MJDPeriod.push_back(std::make_pair(58031.78,58061.23));
  MJDPeriod.push_back(std::make_pair(58061.23,58090.66));
  MJDPeriod.push_back(std::make_pair(58090.66,58120.10));
  MJDPeriod.push_back(std::make_pair(58120.10,58149.56));
  MJDPeriod.push_back(std::make_pair(58149.56,58179.04));
  MJDPeriod.push_back(std::make_pair(58179.04,58208.53));
  MJDPeriod.push_back(std::make_pair(58208.53,58238.04));
  MJDPeriod.push_back(std::make_pair(58238.04,58267.60));
  MJDPeriod.push_back(std::make_pair(58267.60,58297.20));
  MJDPeriod.push_back(std::make_pair(58297.20,58326.85));
  MJDPeriod.push_back(std::make_pair(58326.85,58356.50));
  MJDPeriod.push_back(std::make_pair(58356.50,58386.12));
  MJDPeriod.push_back(std::make_pair(58386.12,58415.70));
  MJDPeriod.push_back(std::make_pair(58415.70,58445.24));
  MJDPeriod.push_back(std::make_pair(58445.24,58474.74));
  MJDPeriod.push_back(std::make_pair(58474.74,58504.22));
  MJDPeriod.push_back(std::make_pair(58504.22,58533.66));
  MJDPeriod.push_back(std::make_pair(58533.66,58563.07));
  MJDPeriod.push_back(std::make_pair(58563.07,58592.47));
  MJDPeriod.push_back(std::make_pair(58592.47,58621.88));
  MJDPeriod.push_back(std::make_pair(58621.88,58651.36));
  MJDPeriod.push_back(std::make_pair(58651.36,58680.90));
  MJDPeriod.push_back(std::make_pair(58680.90,58710.52));
  MJDPeriod.push_back(std::make_pair(58710.52,58740.19));
  MJDPeriod.push_back(std::make_pair(58740.19,58769.88));
  MJDPeriod.push_back(std::make_pair(58769.88,58799.57));
  MJDPeriod.push_back(std::make_pair(58799.57,58829.22));
  MJDPeriod.push_back(std::make_pair(58829.22,58858.81));
  MJDPeriod.push_back(std::make_pair(58858.81,58888.32));
  MJDPeriod.push_back(std::make_pair(58888.32,58917.74));
  MJDPeriod.push_back(std::make_pair(58917.74,58947.11));
  MJDPeriod.push_back(std::make_pair(58947.11,58976.45));
  MJDPeriod.push_back(std::make_pair(58976.45,59005.80));
  MJDPeriod.push_back(std::make_pair(59005.80,59035.20));
  MJDPeriod.push_back(std::make_pair(59035.20,59064.67));
  MJDPeriod.push_back(std::make_pair(59064.67,59094.22));
  MJDPeriod.push_back(std::make_pair(59094.22,59123.88));
  MJDPeriod.push_back(std::make_pair(59123.88,59153.62));
  MJDPeriod.push_back(std::make_pair(59153.62,59183.40));
  MJDPeriod.push_back(std::make_pair(59183.40,59213.15));
  MJDPeriod.push_back(std::make_pair(59213.15,59242.80));
  MJDPeriod.push_back(std::make_pair(59242.80,59272.35));
  MJDPeriod.push_back(std::make_pair(59272.35,59301.78));
  MJDPeriod.push_back(std::make_pair(59301.78,59331.15));
  MJDPeriod.push_back(std::make_pair(59331.15,59360.47));
  MJDPeriod.push_back(std::make_pair(59360.47,59389.78));
  MJDPeriod.push_back(std::make_pair(59389.78,59419.11));
  MJDPeriod.push_back(std::make_pair(59419.11,59448.50));
  MJDPeriod.push_back(std::make_pair(59448.50,59478.00));
  MJDPeriod.push_back(std::make_pair(59478.00,59507.62));
  MJDPeriod.push_back(std::make_pair(59507.62,59537.37));
  MJDPeriod.push_back(std::make_pair(59537.37,59567.19));
  MJDPeriod.push_back(std::make_pair(59567.19,59596.99));
  MJDPeriod.push_back(std::make_pair(59596.99,59626.71));
  MJDPeriod.push_back(std::make_pair(59626.71,59656.31));
  MJDPeriod.push_back(std::make_pair(59656.31,59685.79));
  MJDPeriod.push_back(std::make_pair(59685.79,59715.18));
  MJDPeriod.push_back(std::make_pair(59715.18,59744.49));
  MJDPeriod.push_back(std::make_pair(59744.49,59773.78));
  MJDPeriod.push_back(std::make_pair(59773.78,59803.07));
  MJDPeriod.push_back(std::make_pair(59803.07,59832.42));
  MJDPeriod.push_back(std::make_pair(59832.42,59861.87));
  MJDPeriod.push_back(std::make_pair(59861.87,59891.46));
  MJDPeriod.push_back(std::make_pair(59891.46,59921.17));
  MJDPeriod.push_back(std::make_pair(59921.17,59950.97));
  MJDPeriod.push_back(std::make_pair(59950.97,59980.77));
  MJDPeriod.push_back(std::make_pair(59980.77,60010.53));
  MJDPeriod.push_back(std::make_pair(60010.53,60040.19));
  MJDPeriod.push_back(std::make_pair(60040.19,60069.73));
  MJDPeriod.push_back(std::make_pair(60069.73,60099.16));
  MJDPeriod.push_back(std::make_pair(60099.16,60128.49));
  MJDPeriod.push_back(std::make_pair(60128.49,60157.77));
  MJDPeriod.push_back(std::make_pair(60157.77,60187.07));
  MJDPeriod.push_back(std::make_pair(60187.07,60216.42));
  MJDPeriod.push_back(std::make_pair(60216.42,60245.85));
  MJDPeriod.push_back(std::make_pair(60245.85,60275.39));
  MJDPeriod.push_back(std::make_pair(60275.39,60305.02));
  MJDPeriod.push_back(std::make_pair(60305.02,60334.75));
  MJDPeriod.push_back(std::make_pair(60334.75,60364.52));
  MJDPeriod.push_back(std::make_pair(60364.52,60394.29));
  MJDPeriod.push_back(std::make_pair(60394.29,60423.99));
  MJDPeriod.push_back(std::make_pair(60423.99,60453.58));
  MJDPeriod.push_back(std::make_pair(60453.58,60483.05));
  MJDPeriod.push_back(std::make_pair(60483.05,60512.43));
  MJDPeriod.push_back(std::make_pair(60512.43,60541.77));
  MJDPeriod.push_back(std::make_pair(60541.77,60571.11));
  MJDPeriod.push_back(std::make_pair(60571.11,60600.48));
  MJDPeriod.push_back(std::make_pair(60600.48,60629.90));
  MJDPeriod.push_back(std::make_pair(60629.90,60659.38));
  MJDPeriod.push_back(std::make_pair(60659.38,60688.94));
  MJDPeriod.push_back(std::make_pair(60688.94,60718.58));
  MJDPeriod.push_back(std::make_pair(60718.58,60748.29));
  MJDPeriod.push_back(std::make_pair(60748.29,60778.02));
  MJDPeriod.push_back(std::make_pair(60778.02,60807.71));
  MJDPeriod.push_back(std::make_pair(60807.71,60837.32));
  MJDPeriod.push_back(std::make_pair(60837.32,60866.86));
  MJDPeriod.push_back(std::make_pair(60866.86,60896.33));
  MJDPeriod.push_back(std::make_pair(60896.33,60925.76));
  MJDPeriod.push_back(std::make_pair(60925.76,60955.16));
  MJDPeriod.push_back(std::make_pair(60955.16,60984.56));
  MJDPeriod.push_back(std::make_pair(60984.56,61013.97));
  MJDPeriod.push_back(std::make_pair(61013.97,61043.42));
  MJDPeriod.push_back(std::make_pair(61043.42,61072.92));
  MJDPeriod.push_back(std::make_pair(61072.92,61102.49));
  MJDPeriod.push_back(std::make_pair(61102.49,61132.09));
  MJDPeriod.push_back(std::make_pair(61132.09,61161.73));
  MJDPeriod.push_back(std::make_pair(61161.73,61191.37));
  MJDPeriod.push_back(std::make_pair(61191.37,61221.00));
  MJDPeriod.push_back(std::make_pair(61221.00,61250.61));
  MJDPeriod.push_back(std::make_pair(61250.61,61280.18));
  MJDPeriod.push_back(std::make_pair(61280.18,61309.70));
  MJDPeriod.push_back(std::make_pair(61309.70,61339.18));
  MJDPeriod.push_back(std::make_pair(61339.18,61368.62));
  MJDPeriod.push_back(std::make_pair(61368.62,61398.06));
  MJDPeriod.push_back(std::make_pair(61398.06,61427.51));
  MJDPeriod.push_back(std::make_pair(61427.51,61456.98));
  MJDPeriod.push_back(std::make_pair(61456.98,61486.45));
  MJDPeriod.push_back(std::make_pair(61486.45,61515.94));
  MJDPeriod.push_back(std::make_pair(61515.94,61545.46));
  MJDPeriod.push_back(std::make_pair(61545.46,61575.03));
  MJDPeriod.push_back(std::make_pair(61575.03,61604.66));
  MJDPeriod.push_back(std::make_pair(61604.66,61634.31));
  MJDPeriod.push_back(std::make_pair(61634.31,61663.96));
  MJDPeriod.push_back(std::make_pair(61663.96,61693.58));
  MJDPeriod.push_back(std::make_pair(61693.58,61723.14));
  MJDPeriod.push_back(std::make_pair(61723.14,61752.67));
  MJDPeriod.push_back(std::make_pair(61752.67,61782.17));
  MJDPeriod.push_back(std::make_pair(61782.17,61811.63));
  MJDPeriod.push_back(std::make_pair(61811.63,61841.05));
  MJDPeriod.push_back(std::make_pair(61841.05,61870.44));
  MJDPeriod.push_back(std::make_pair(61870.44,61899.83));
  MJDPeriod.push_back(std::make_pair(61899.83,61929.26));
  MJDPeriod.push_back(std::make_pair(61929.26,61958.76));
  MJDPeriod.push_back(std::make_pair(61958.76,61988.34));
  MJDPeriod.push_back(std::make_pair(61988.34,62017.99));
  MJDPeriod.push_back(std::make_pair(62017.99,62047.69));
  MJDPeriod.push_back(std::make_pair(62047.69,62077.39));
  MJDPeriod.push_back(std::make_pair(62077.39,62107.07));
  MJDPeriod.push_back(std::make_pair(62107.07,62136.70));
  MJDPeriod.push_back(std::make_pair(62136.70,62166.25));
  MJDPeriod.push_back(std::make_pair(62166.25,62195.72));
  MJDPeriod.push_back(std::make_pair(62195.72,62225.10));
  MJDPeriod.push_back(std::make_pair(62225.10,62254.44));
  MJDPeriod.push_back(std::make_pair(62254.44,62283.78));
  MJDPeriod.push_back(std::make_pair(62283.78,62313.14));
  MJDPeriod.push_back(std::make_pair(62313.14,62342.57));
  MJDPeriod.push_back(std::make_pair(62342.57,62372.08));
  MJDPeriod.push_back(std::make_pair(62372.08,62401.69));
  MJDPeriod.push_back(std::make_pair(62401.69,62431.40));
  MJDPeriod.push_back(std::make_pair(62431.40,62461.17));
  MJDPeriod.push_back(std::make_pair(62461.17,62490.95));
  MJDPeriod.push_back(std::make_pair(62490.95,62520.66));
  MJDPeriod.push_back(std::make_pair(62520.66,62550.27));
  MJDPeriod.push_back(std::make_pair(62550.27,62579.75));
  MJDPeriod.push_back(std::make_pair(62579.75,62609.14));
  MJDPeriod.push_back(std::make_pair(62609.14,62638.47));
  MJDPeriod.push_back(std::make_pair(62638.47,62667.78));
  MJDPeriod.push_back(std::make_pair(62667.78,62697.09));
  MJDPeriod.push_back(std::make_pair(62697.09,62726.45));
  MJDPeriod.push_back(std::make_pair(62726.45,62755.89));
  MJDPeriod.push_back(std::make_pair(62755.89,62785.45));
  MJDPeriod.push_back(std::make_pair(62785.45,62815.15));
  MJDPeriod.push_back(std::make_pair(62815.15,62844.95));
  MJDPeriod.push_back(std::make_pair(62844.95,62874.77));
  MJDPeriod.push_back(std::make_pair(62874.77,62904.53));
  MJDPeriod.push_back(std::make_pair(62904.53,62934.19));
  MJDPeriod.push_back(std::make_pair(62934.19,62963.72));
  MJDPeriod.push_back(std::make_pair(62963.72,62993.15));
  MJDPeriod.push_back(std::make_pair(62993.15,63022.50));
  MJDPeriod.push_back(std::make_pair(63022.50,63051.79));
  MJDPeriod.push_back(std::make_pair(63051.79,63081.07));
  MJDPeriod.push_back(std::make_pair(63081.07,63110.39));
  MJDPeriod.push_back(std::make_pair(63110.39,63139.79));
  MJDPeriod.push_back(std::make_pair(63139.79,63169.31));
  MJDPeriod.push_back(std::make_pair(63169.31,63198.97));
  MJDPeriod.push_back(std::make_pair(63198.97,63228.73));
  MJDPeriod.push_back(std::make_pair(63228.73,63258.54));
  MJDPeriod.push_back(std::make_pair(63258.54,63288.32));
  MJDPeriod.push_back(std::make_pair(63288.32,63318.03));
  MJDPeriod.push_back(std::make_pair(63318.03,63347.63));
  MJDPeriod.push_back(std::make_pair(63347.63,63377.11));
  MJDPeriod.push_back(std::make_pair(63377.11,63406.48));
  MJDPeriod.push_back(std::make_pair(63406.48,63435.79));
  MJDPeriod.push_back(std::make_pair(63435.79,63465.08));
  MJDPeriod.push_back(std::make_pair(63465.08,63494.40));
  MJDPeriod.push_back(std::make_pair(63494.40,63523.79));
  MJDPeriod.push_back(std::make_pair(63523.79,63553.28));
  MJDPeriod.push_back(std::make_pair(63553.28,63582.87));
  MJDPeriod.push_back(std::make_pair(63582.87,63612.55));
  MJDPeriod.push_back(std::make_pair(63612.55,63642.29));
  MJDPeriod.push_back(std::make_pair(63642.29,63672.07));
  MJDPeriod.push_back(std::make_pair(63672.07,63701.80));
  MJDPeriod.push_back(std::make_pair(63701.80,63731.45));
  MJDPeriod.push_back(std::make_pair(63731.45,63760.97));
  MJDPeriod.push_back(std::make_pair(63760.97,63790.40));
  MJDPeriod.push_back(std::make_pair(63790.40,63819.76));
  MJDPeriod.push_back(std::make_pair(63819.76,63849.10));
  MJDPeriod.push_back(std::make_pair(63849.10,63878.46));
  MJDPeriod.push_back(std::make_pair(63878.46,63907.86));
  MJDPeriod.push_back(std::make_pair(63907.86,63937.31));
  MJDPeriod.push_back(std::make_pair(63937.31,63966.83));
  MJDPeriod.push_back(std::make_pair(63966.83,63996.42));
  MJDPeriod.push_back(std::make_pair(63996.42,64026.09));
  MJDPeriod.push_back(std::make_pair(64026.09,64055.81));
  MJDPeriod.push_back(std::make_pair(64055.81,64085.51));
  MJDPeriod.push_back(std::make_pair(64085.51,64115.16));
  MJDPeriod.push_back(std::make_pair(64115.16,64144.74));
  MJDPeriod.push_back(std::make_pair(64144.74,64174.25));
  MJDPeriod.push_back(std::make_pair(64174.25,64203.70));
  MJDPeriod.push_back(std::make_pair(64203.70,64233.12));
  MJDPeriod.push_back(std::make_pair(64233.12,64262.53));
  MJDPeriod.push_back(std::make_pair(64262.53,64291.94));
  MJDPeriod.push_back(std::make_pair(64291.94,64321.37));
  MJDPeriod.push_back(std::make_pair(64321.37,64350.85));
  MJDPeriod.push_back(std::make_pair(64350.85,64380.37));
  MJDPeriod.push_back(std::make_pair(64380.37,64409.95));
  MJDPeriod.push_back(std::make_pair(64409.95,64439.56));
  MJDPeriod.push_back(std::make_pair(64439.56,64469.18));
  MJDPeriod.push_back(std::make_pair(64469.18,64498.82));
  MJDPeriod.push_back(std::make_pair(64498.82,64528.44));
  MJDPeriod.push_back(std::make_pair(64528.44,64558.04));
  MJDPeriod.push_back(std::make_pair(64558.04,64587.60));
  MJDPeriod.push_back(std::make_pair(64587.60,64617.11));
  MJDPeriod.push_back(std::make_pair(64617.11,64646.58));
  MJDPeriod.push_back(std::make_pair(64646.58,64676.02));
  MJDPeriod.push_back(std::make_pair(64676.02,64705.47));
  MJDPeriod.push_back(std::make_pair(64705.47,64734.92));
  MJDPeriod.push_back(std::make_pair(64734.92,64764.38));
  MJDPeriod.push_back(std::make_pair(64764.38,64793.85));
  MJDPeriod.push_back(std::make_pair(64793.85,64823.34));
  MJDPeriod.push_back(std::make_pair(64823.34,64852.88));
  MJDPeriod.push_back(std::make_pair(64852.88,64882.47));
  MJDPeriod.push_back(std::make_pair(64882.47,64912.12));
  MJDPeriod.push_back(std::make_pair(64912.12,64941.78));
  MJDPeriod.push_back(std::make_pair(64941.78,64971.43));
  MJDPeriod.push_back(std::make_pair(64971.43,65001.03));
  MJDPeriod.push_back(std::make_pair(65001.03,65030.59));
  MJDPeriod.push_back(std::make_pair(65030.59,65060.11));
  MJDPeriod.push_back(std::make_pair(65060.11,65089.59));
  MJDPeriod.push_back(std::make_pair(65089.59,65119.02));
  MJDPeriod.push_back(std::make_pair(65119.02,65148.41));
  MJDPeriod.push_back(std::make_pair(65148.41,65177.79));
  MJDPeriod.push_back(std::make_pair(65177.79,65207.18));
  MJDPeriod.push_back(std::make_pair(65207.18,65236.64));
  MJDPeriod.push_back(std::make_pair(65236.64,65266.18));
  MJDPeriod.push_back(std::make_pair(65266.18,65295.80));
  MJDPeriod.push_back(std::make_pair(65295.80,65325.48));
  MJDPeriod.push_back(std::make_pair(65325.48,65355.19));
  MJDPeriod.push_back(std::make_pair(65355.19,65384.90));
  MJDPeriod.push_back(std::make_pair(65384.90,65414.57));
  MJDPeriod.push_back(std::make_pair(65414.57,65444.17));
  MJDPeriod.push_back(std::make_pair(65444.17,65473.67));
  MJDPeriod.push_back(std::make_pair(65473.67,65503.09));
  MJDPeriod.push_back(std::make_pair(65503.09,65532.44));
  MJDPeriod.push_back(std::make_pair(65532.44,65561.77));
  MJDPeriod.push_back(std::make_pair(65561.77,65591.10));
  MJDPeriod.push_back(std::make_pair(65591.10,65620.49));
  MJDPeriod.push_back(std::make_pair(65620.49,65649.96));
  MJDPeriod.push_back(std::make_pair(65649.96,65679.52));
  MJDPeriod.push_back(std::make_pair(65679.52,65709.18));
  MJDPeriod.push_back(std::make_pair(65709.18,65738.94));
  MJDPeriod.push_back(std::make_pair(65738.94,65768.73));
  MJDPeriod.push_back(std::make_pair(65768.73,65798.49));
  MJDPeriod.push_back(std::make_pair(65798.49,65828.15));
  MJDPeriod.push_back(std::make_pair(65828.15,65857.69));
  MJDPeriod.push_back(std::make_pair(65857.69,65887.12));
  MJDPeriod.push_back(std::make_pair(65887.12,65916.47));
  MJDPeriod.push_back(std::make_pair(65916.47,65945.78));
  MJDPeriod.push_back(std::make_pair(65945.78,65975.09));
  MJDPeriod.push_back(std::make_pair(65975.09,66004.41));
  MJDPeriod.push_back(std::make_pair(66004.41,66033.81));
  MJDPeriod.push_back(std::make_pair(66033.81,66063.31));
  MJDPeriod.push_back(std::make_pair(66063.31,66092.94));
  MJDPeriod.push_back(std::make_pair(66092.94,66122.70));
  MJDPeriod.push_back(std::make_pair(66122.70,66152.53));
  MJDPeriod.push_back(std::make_pair(66152.53,66182.33));
  MJDPeriod.push_back(std::make_pair(66182.33,66212.04));
  MJDPeriod.push_back(std::make_pair(66212.04,66241.64));
  MJDPeriod.push_back(std::make_pair(66241.64,66271.11));
  MJDPeriod.push_back(std::make_pair(66271.11,66300.49));
  MJDPeriod.push_back(std::make_pair(66300.49,66329.81));
  MJDPeriod.push_back(std::make_pair(66329.81,66359.09));
  MJDPeriod.push_back(std::make_pair(66359.09,66388.38));
  MJDPeriod.push_back(std::make_pair(66388.38,66417.74));
  MJDPeriod.push_back(std::make_pair(66417.74,66447.20));
  MJDPeriod.push_back(std::make_pair(66447.20,66476.80));
  MJDPeriod.push_back(std::make_pair(66476.80,66506.51));
  MJDPeriod.push_back(std::make_pair(66506.51,66536.30));
  MJDPeriod.push_back(std::make_pair(66536.30,66566.10));
  MJDPeriod.push_back(std::make_pair(66566.10,66595.85));
  MJDPeriod.push_back(std::make_pair(66595.85,66625.50));
  MJDPeriod.push_back(std::make_pair(66625.50,66655.04));
  MJDPeriod.push_back(std::make_pair(66655.04,66684.46));
  MJDPeriod.push_back(std::make_pair(66684.46,66713.79));
  MJDPeriod.push_back(std::make_pair(66713.79,66743.09));
  MJDPeriod.push_back(std::make_pair(66743.09,66772.39));
  MJDPeriod.push_back(std::make_pair(66772.39,66801.75));
  MJDPeriod.push_back(std::make_pair(66801.75,66831.20));
  MJDPeriod.push_back(std::make_pair(66831.20,66860.74));
  MJDPeriod.push_back(std::make_pair(66860.74,66890.37));
  MJDPeriod.push_back(std::make_pair(66890.37,66920.08));
  MJDPeriod.push_back(std::make_pair(66920.08,66949.84));
  MJDPeriod.push_back(std::make_pair(66949.84,66979.60));
  MJDPeriod.push_back(std::make_pair(66979.60,67009.29));
  MJDPeriod.push_back(std::make_pair(67009.29,67038.87));
  MJDPeriod.push_back(std::make_pair(67038.87,67068.34));
  MJDPeriod.push_back(std::make_pair(67068.34,67097.73));
  MJDPeriod.push_back(std::make_pair(67097.73,67127.09));
  MJDPeriod.push_back(std::make_pair(67127.09,67156.44));
  MJDPeriod.push_back(std::make_pair(67156.44,67185.83));
  MJDPeriod.push_back(std::make_pair(67185.83,67215.25));
  MJDPeriod.push_back(std::make_pair(67215.25,67244.74));
  MJDPeriod.push_back(std::make_pair(67244.74,67274.29));
  MJDPeriod.push_back(std::make_pair(67274.29,67303.92));
  MJDPeriod.push_back(std::make_pair(67303.92,67333.60));
  MJDPeriod.push_back(std::make_pair(67333.60,67363.31));
  MJDPeriod.push_back(std::make_pair(67363.31,67392.98));
  MJDPeriod.push_back(std::make_pair(67392.98,67422.60));
  MJDPeriod.push_back(std::make_pair(67422.60,67452.14));
  MJDPeriod.push_back(std::make_pair(67452.14,67481.63));
  MJDPeriod.push_back(std::make_pair(67481.63,67511.08));
  MJDPeriod.push_back(std::make_pair(67511.08,67540.50));
  MJDPeriod.push_back(std::make_pair(67540.50,67569.91));
  MJDPeriod.push_back(std::make_pair(67569.91,67599.34));
  MJDPeriod.push_back(std::make_pair(67599.34,67628.79));
  MJDPeriod.push_back(std::make_pair(67628.79,67658.28));
  MJDPeriod.push_back(std::make_pair(67658.28,67687.82));
  MJDPeriod.push_back(std::make_pair(67687.82,67717.40));
  MJDPeriod.push_back(std::make_pair(67717.40,67747.01));
  MJDPeriod.push_back(std::make_pair(67747.01,67776.64));
  MJDPeriod.push_back(std::make_pair(67776.64,67806.27));
  MJDPeriod.push_back(std::make_pair(67806.27,67835.89));
  MJDPeriod.push_back(std::make_pair(67835.89,67865.48));
  MJDPeriod.push_back(std::make_pair(67865.48,67895.02));
  MJDPeriod.push_back(std::make_pair(67895.02,67924.52));
  MJDPeriod.push_back(std::make_pair(67924.52,67953.98));
  MJDPeriod.push_back(std::make_pair(67953.98,67983.43));
  MJDPeriod.push_back(std::make_pair(67983.43,68012.88));
  MJDPeriod.push_back(std::make_pair(68012.88,68042.33));
  MJDPeriod.push_back(std::make_pair(68042.33,68071.78));
  MJDPeriod.push_back(std::make_pair(68071.78,68101.24));
  MJDPeriod.push_back(std::make_pair(68101.24,68130.75));
  MJDPeriod.push_back(std::make_pair(68130.75,68160.30));
  MJDPeriod.push_back(std::make_pair(68160.30,68189.93));
  MJDPeriod.push_back(std::make_pair(68189.93,68219.59));
  MJDPeriod.push_back(std::make_pair(68219.59,68249.26));
  MJDPeriod.push_back(std::make_pair(68249.26,68278.90));
  MJDPeriod.push_back(std::make_pair(68278.90,68308.49));
  MJDPeriod.push_back(std::make_pair(68308.49,68338.04));
  MJDPeriod.push_back(std::make_pair(68338.04,68367.54));
  MJDPeriod.push_back(std::make_pair(68367.54,68396.99));
  MJDPeriod.push_back(std::make_pair(68396.99,68426.39));
  MJDPeriod.push_back(std::make_pair(68426.39,68455.76));
  MJDPeriod.push_back(std::make_pair(68455.76,68485.14));
  MJDPeriod.push_back(std::make_pair(68485.14,68514.55));
  MJDPeriod.push_back(std::make_pair(68514.55,68544.04));
  MJDPeriod.push_back(std::make_pair(68544.04,68573.62));
  MJDPeriod.push_back(std::make_pair(68573.62,68603.28));
  MJDPeriod.push_back(std::make_pair(68603.28,68632.99));
  MJDPeriod.push_back(std::make_pair(68632.99,68662.71));
  MJDPeriod.push_back(std::make_pair(68662.71,68692.42));
  MJDPeriod.push_back(std::make_pair(68692.42,68722.06));
  MJDPeriod.push_back(std::make_pair(68722.06,68751.61));
  MJDPeriod.push_back(std::make_pair(68751.61,68781.07));
  MJDPeriod.push_back(std::make_pair(68781.07,68810.44));
  MJDPeriod.push_back(std::make_pair(68810.44,68839.77));
  MJDPeriod.push_back(std::make_pair(68839.77,68869.09));
  MJDPeriod.push_back(std::make_pair(68869.09,68898.44));
  MJDPeriod.push_back(std::make_pair(68898.44,68927.86));
  MJDPeriod.push_back(std::make_pair(68927.86,68957.37));
  MJDPeriod.push_back(std::make_pair(68957.37,68986.99));
  MJDPeriod.push_back(std::make_pair(68986.99,69016.71));
  MJDPeriod.push_back(std::make_pair(69016.71,69046.50));
  MJDPeriod.push_back(std::make_pair(69046.50,69076.29));
  MJDPeriod.push_back(std::make_pair(69076.29,69106.01));
  MJDPeriod.push_back(std::make_pair(69106.01,69135.61));
  MJDPeriod.push_back(std::make_pair(69135.61,69165.09));
  MJDPeriod.push_back(std::make_pair(69165.09,69194.47));
  MJDPeriod.push_back(std::make_pair(69194.47,69223.79));
  MJDPeriod.push_back(std::make_pair(69223.79,69253.09));
  MJDPeriod.push_back(std::make_pair(69253.09,69282.40));
  MJDPeriod.push_back(std::make_pair(69282.40,69311.76));
  MJDPeriod.push_back(std::make_pair(69311.76,69341.20));
  MJDPeriod.push_back(std::make_pair(69341.20,69370.77));
  MJDPeriod.push_back(std::make_pair(69370.77,69400.47));
  MJDPeriod.push_back(std::make_pair(69400.47,69430.28));
  MJDPeriod.push_back(std::make_pair(69430.28,69460.10));
  MJDPeriod.push_back(std::make_pair(69460.10,69489.87));
  MJDPeriod.push_back(std::make_pair(69489.87,69519.52));
  MJDPeriod.push_back(std::make_pair(69519.52,69549.05));
  MJDPeriod.push_back(std::make_pair(69549.05,69578.47));
  MJDPeriod.push_back(std::make_pair(69578.47,69607.81));
  MJDPeriod.push_back(std::make_pair(69607.81,69637.10));
  MJDPeriod.push_back(std::make_pair(69637.10,69666.39));
  MJDPeriod.push_back(std::make_pair(69666.39,69695.71));
  MJDPeriod.push_back(std::make_pair(69695.71,69725.12));
  MJDPeriod.push_back(std::make_pair(69725.12,69754.65));
  MJDPeriod.push_back(std::make_pair(69754.65,69784.31));
  MJDPeriod.push_back(std::make_pair(69784.31,69814.07));
  MJDPeriod.push_back(std::make_pair(69814.07,69843.87));
  MJDPeriod.push_back(std::make_pair(69843.87,69873.64));
  MJDPeriod.push_back(std::make_pair(69873.64,69903.34));
  MJDPeriod.push_back(std::make_pair(69903.34,69932.94));
  MJDPeriod.push_back(std::make_pair(69932.94,69962.41));
  MJDPeriod.push_back(std::make_pair(69962.41,69991.79));
  MJDPeriod.push_back(std::make_pair(69991.79,70021.10));
  MJDPeriod.push_back(std::make_pair(70021.10,70050.40));
  MJDPeriod.push_back(std::make_pair(70050.40,70079.73));
  MJDPeriod.push_back(std::make_pair(70079.73,70109.14));
  MJDPeriod.push_back(std::make_pair(70109.14,70138.63));
  MJDPeriod.push_back(std::make_pair(70138.63,70168.22));
  MJDPeriod.push_back(std::make_pair(70168.22,70197.89));
  MJDPeriod.push_back(std::make_pair(70197.89,70227.62));
  MJDPeriod.push_back(std::make_pair(70227.62,70257.38));
  MJDPeriod.push_back(std::make_pair(70257.38,70287.10));
  MJDPeriod.push_back(std::make_pair(70287.10,70316.73));
  MJDPeriod.push_back(std::make_pair(70316.73,70346.26));
  MJDPeriod.push_back(std::make_pair(70346.26,70375.69));
  MJDPeriod.push_back(std::make_pair(70375.69,70405.07));
  MJDPeriod.push_back(std::make_pair(70405.07,70434.43));
  MJDPeriod.push_back(std::make_pair(70434.43,70463.80));
  MJDPeriod.push_back(std::make_pair(70463.80,70493.21));
  MJDPeriod.push_back(std::make_pair(70493.21,70522.67));
  MJDPeriod.push_back(std::make_pair(70522.67,70552.18));
  MJDPeriod.push_back(std::make_pair(70552.18,70581.77));
  MJDPeriod.push_back(std::make_pair(70581.77,70611.41));
  MJDPeriod.push_back(std::make_pair(70611.41,70641.10));
  MJDPeriod.push_back(std::make_pair(70641.10,70670.79));
  MJDPeriod.push_back(std::make_pair(70670.79,70700.44));
  MJDPeriod.push_back(std::make_pair(70700.44,70730.02));
  MJDPeriod.push_back(std::make_pair(70730.02,70759.54));
  MJDPeriod.push_back(std::make_pair(70759.54,70789.01));
  MJDPeriod.push_back(std::make_pair(70789.01,70818.46));
  MJDPeriod.push_back(std::make_pair(70818.46,70847.88));
  MJDPeriod.push_back(std::make_pair(70847.88,70877.30));
  MJDPeriod.push_back(std::make_pair(70877.30,70906.74));
  MJDPeriod.push_back(std::make_pair(70906.74,70936.21));
  MJDPeriod.push_back(std::make_pair(70936.21,70965.72));
  MJDPeriod.push_back(std::make_pair(70965.72,70995.27));
  MJDPeriod.push_back(std::make_pair(70995.27,71024.85));
  MJDPeriod.push_back(std::make_pair(71024.85,71054.46));
  MJDPeriod.push_back(std::make_pair(71054.46,71084.09));
  MJDPeriod.push_back(std::make_pair(71084.09,71113.71));
  MJDPeriod.push_back(std::make_pair(71113.71,71143.33));
  MJDPeriod.push_back(std::make_pair(71143.33,71172.91));
  MJDPeriod.push_back(std::make_pair(71172.91,71202.44));
  MJDPeriod.push_back(std::make_pair(71202.44,71231.93));
  MJDPeriod.push_back(std::make_pair(71231.93,71261.39));
  MJDPeriod.push_back(std::make_pair(71261.39,71290.84));
  MJDPeriod.push_back(std::make_pair(71290.84,71320.28));
  MJDPeriod.push_back(std::make_pair(71320.28,71349.72));
  MJDPeriod.push_back(std::make_pair(71349.72,71379.17));
  MJDPeriod.push_back(std::make_pair(71379.17,71408.64));
  MJDPeriod.push_back(std::make_pair(71408.64,71438.16));
  MJDPeriod.push_back(std::make_pair(71438.16,71467.74));
  MJDPeriod.push_back(std::make_pair(71467.74,71497.39));
  MJDPeriod.push_back(std::make_pair(71497.39,71527.07));
  MJDPeriod.push_back(std::make_pair(71527.07,71556.74));
  MJDPeriod.push_back(std::make_pair(71556.74,71586.37));
  MJDPeriod.push_back(std::make_pair(71586.37,71615.95));
  MJDPeriod.push_back(std::make_pair(71615.95,71645.47));
  MJDPeriod.push_back(std::make_pair(71645.47,71674.95));
  MJDPeriod.push_back(std::make_pair(71674.95,71704.37));
  MJDPeriod.push_back(std::make_pair(71704.37,71733.75));
  MJDPeriod.push_back(std::make_pair(71733.75,71763.11));
  MJDPeriod.push_back(std::make_pair(71763.11,71792.48));
  MJDPeriod.push_back(std::make_pair(71792.48,71821.93));
  MJDPeriod.push_back(std::make_pair(71821.93,71851.46));
  MJDPeriod.push_back(std::make_pair(71851.46,71881.08));
  MJDPeriod.push_back(std::make_pair(71881.08,71910.78));
  MJDPeriod.push_back(std::make_pair(71910.78,71940.51));
  MJDPeriod.push_back(std::make_pair(71940.51,71970.24));
  MJDPeriod.push_back(std::make_pair(71970.24,71999.92));
  MJDPeriod.push_back(std::make_pair(71999.92,72029.53));
  MJDPeriod.push_back(std::make_pair(72029.53,72059.03));
  MJDPeriod.push_back(std::make_pair(72059.03,72088.44));
  MJDPeriod.push_back(std::make_pair(72088.44,72117.77));
  MJDPeriod.push_back(std::make_pair(72117.77,72147.08));
  MJDPeriod.push_back(std::make_pair(72147.08,72176.41));
  MJDPeriod.push_back(std::make_pair(72176.41,72205.79));
  MJDPeriod.push_back(std::make_pair(72205.79,72235.25));
  MJDPeriod.push_back(std::make_pair(72235.25,72264.82));
  MJDPeriod.push_back(std::make_pair(72264.82,72294.49));
  MJDPeriod.push_back(std::make_pair(72294.49,72324.26));
  MJDPeriod.push_back(std::make_pair(72324.26,72354.07));
  MJDPeriod.push_back(std::make_pair(72354.07,72383.84));
  MJDPeriod.push_back(std::make_pair(72383.84,72413.50));
  MJDPeriod.push_back(std::make_pair(72413.50,72443.03));
  MJDPeriod.push_back(std::make_pair(72443.03,72472.45));
  MJDPeriod.push_back(std::make_pair(72472.45,72501.79));
  MJDPeriod.push_back(std::make_pair(72501.79,72531.10));
  MJDPeriod.push_back(std::make_pair(72531.10,72560.40));
  MJDPeriod.push_back(std::make_pair(72560.40,72589.72));
  MJDPeriod.push_back(std::make_pair(72589.72,72619.12));
  MJDPeriod.push_back(std::make_pair(72619.12,72648.63));
  MJDPeriod.push_back(std::make_pair(72648.63,72678.27));
  MJDPeriod.push_back(std::make_pair(72678.27,72708.03));
  MJDPeriod.push_back(std::make_pair(72708.03,72737.86));
  MJDPeriod.push_back(std::make_pair(72737.86,72767.66));
  MJDPeriod.push_back(std::make_pair(72767.66,72797.37));
  MJDPeriod.push_back(std::make_pair(72797.37,72826.96));
  MJDPeriod.push_back(std::make_pair(72826.96,72856.43));
  MJDPeriod.push_back(std::make_pair(72856.43,72885.80));
  MJDPeriod.push_back(std::make_pair(72885.80,72915.12));
  MJDPeriod.push_back(std::make_pair(72915.12,72944.40));
  MJDPeriod.push_back(std::make_pair(72944.40,72973.70));
  MJDPeriod.push_back(std::make_pair(72973.70,73003.07));
  MJDPeriod.push_back(std::make_pair(73003.07,73032.54));
  MJDPeriod.push_back(std::make_pair(73032.54,73062.14));
  MJDPeriod.push_back(std::make_pair(73062.14,73091.85));
  MJDPeriod.push_back(std::make_pair(73091.85,73121.63));
  MJDPeriod.push_back(std::make_pair(73121.63,73151.42));
  MJDPeriod.push_back(std::make_pair(73151.42,73181.16));
  MJDPeriod.push_back(std::make_pair(73181.16,73210.81));
  MJDPeriod.push_back(std::make_pair(73210.81,73240.34));
  MJDPeriod.push_back(std::make_pair(73240.34,73269.76));
  MJDPeriod.push_back(std::make_pair(73269.76,73299.10));
  MJDPeriod.push_back(std::make_pair(73299.10,73328.41));
  MJDPeriod.push_back(std::make_pair(73328.41,73357.72));
  MJDPeriod.push_back(std::make_pair(73357.72,73387.10));
  MJDPeriod.push_back(std::make_pair(73387.10,73416.55));
  MJDPeriod.push_back(std::make_pair(73416.55,73446.09));
  MJDPeriod.push_back(std::make_pair(73446.09,73475.72));
  MJDPeriod.push_back(std::make_pair(73475.72,73505.42));
  MJDPeriod.push_back(std::make_pair(73505.42,73535.15));
  MJDPeriod.push_back(std::make_pair(73535.15,73564.89));
  MJDPeriod.push_back(std::make_pair(73564.89,73594.57));
  MJDPeriod.push_back(std::make_pair(73594.57,73624.15));
  MJDPeriod.push_back(std::make_pair(73624.15,73653.63));
  MJDPeriod.push_back(std::make_pair(73653.63,73683.04));
  MJDPeriod.push_back(std::make_pair(73683.04,73712.41));
  MJDPeriod.push_back(std::make_pair(73712.41,73741.78));
  MJDPeriod.push_back(std::make_pair(73741.78,73771.18));
  MJDPeriod.push_back(std::make_pair(73771.18,73800.62));

  return MJDPeriod;
}
