// STL
#include <iomanip>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cmath>

// ROOT
#include <TMath.h>
#include <TLine.h>

// START
#include "STARTUtils.hh"
#include "Band.hh"

// Utilities
#define DEBUG 0
#include <debugging.hh>

#define INFO std::cout << INFOCOLOR << "Band> " << RESETCOLOR
#define WARNING std::cout << WARNINGCOLOR << "Band> " << RESETCOLOR

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::Band)
#endif

/**
 * \brief Constructor
 */
START::Band::Band(double zenon, double offset, double efficiency, double alpha, double livetime, unsigned int telcode) 
  : TNamed(),
  fInterpolArea(0),
  fInterpolBiais(0),
  fInterpolResol(0)
{
  fTelCode = telcode;
  fzenon = zenon;
  fzenoff = zenon;
  foffset = offset;
  fazimuth = 0.;
  falpharun = alpha;
  fefficiency = efficiency;
  flivetime = livetime;
  flivetimefraction = 0.;
  fnbrun = 0;
  fkeepband = 1;
  frunstarttime = 0.;
  frunendtime = 0.;
  fEthMC=0;
  fFirstEmcBinNumber=-1;
  fFirstEmcValFit=-1;  
  fFirstEmcValDistrib=-1;  
  fUseInstrumentEnergyDistribution=false;

  fVectorEnergy.clear();
  fVectorResolution.clear();
  fVectorArea.clear();
  fVectorBiais.clear();
  ebin.clear();
  fInterpolTypeForBand = ROOT::Math::Interpolation::kAKIMA;

  fInterpolArea=0;
  fInterpolBiais=0;
  fInterpolResol=0;

  SetName(ConstructBandName().c_str());
  SetTitle(ConstructBandTitle().c_str());
};

/**
 * \brief Destructor
 */
START::Band::~Band()
{

  if (fInterpolResol!=0) delete fInterpolResol; 
  fInterpolResol = 0;

  if (fInterpolBiais!=0) delete fInterpolBiais;
  fInterpolBiais = 0;

  if (fInterpolArea!=0) delete fInterpolArea;
  fInterpolArea = 0;

  DeleteInstrumentInterpTable();
  
}

/**
 * \brief Copy constructor
 */
START::Band::Band(Band const &BandCopy)
  :TNamed(BandCopy),
   fnbrun(BandCopy.fnbrun),
   fTelCode(BandCopy.fTelCode),
   fzenon(BandCopy.fzenon),
   fzenoff(BandCopy.fzenoff),
   foffset(BandCopy.foffset),
   fazimuth(BandCopy.fazimuth),
   falpharun(BandCopy.falpharun),
   fefficiency(BandCopy.fefficiency),
   flivetime(BandCopy.flivetime),
   flivetimefraction(BandCopy.flivetimefraction),
   fkeepband(BandCopy.fkeepband),
   frunstarttime(BandCopy.frunstarttime),
   frunendtime(BandCopy.frunendtime),
   fVectorEnergy(BandCopy.fVectorEnergy),
   fVectorArea(BandCopy.fVectorArea),
   fVectorResolution(BandCopy.fVectorResolution),
   fVectorBiais(BandCopy.fVectorBiais),
   fVectorInterEnergy(BandCopy.fVectorInterEnergy),
   fVectorInterArea(BandCopy.fVectorInterArea),
   fVectorInterResolution(BandCopy.fVectorInterResolution),
   fVectorInterBiais(BandCopy.fVectorInterBiais),
   fDistributionInstrumentMapTable(BandCopy.fDistributionInstrumentMapTable),
   fEthMC(BandCopy.fEthMC),
   fFirstEmcBinNumber(BandCopy.fFirstEmcBinNumber),
   fFirstEmcValFit(BandCopy.fFirstEmcValFit),
   fFirstEmcValDistrib(BandCopy.fFirstEmcValDistrib),
   fUseInstrumentEnergyDistribution(BandCopy.fUseInstrumentEnergyDistribution),
   fInterpolTypeForBand(BandCopy.fInterpolTypeForBand),
   fInterpolArea(0), // JLK
   fInterpolBiais(0), // JLK
   fInterpolResol(0) // JLK
{
  //std::cout << "START::Band> CopyConstructor Called !" << std::endl;
  ebin.clear();

  for(unsigned int ibin(0); ibin<BandCopy.ebin.size(); ibin++) ebin.push_back(BandCopy.ebin[ibin]);

  if (fVectorArea.size()>4 && fVectorEnergy.size()>4) {
    fInterpolArea  = new ROOT::Math::Interpolator(fVectorEnergy,fVectorArea,fInterpolTypeForBand);
  }
  else {
    fInterpolArea  = 0;
  }

  if (fVectorResolution.size()>4 && fVectorEnergy.size()>4) {
  }
  else {
    fInterpolResol  = 0;
  }

  if (fVectorBiais.size()>4 && fVectorEnergy.size()>4) {
    fInterpolBiais  = new ROOT::Math::Interpolator(fVectorEnergy,fVectorBiais,fInterpolTypeForBand);
  }
  else {
    fInterpolBiais  = 0;
  }
  
  InitDistributionInterpTable();  // VIM

}

/**
 * \brief Assignment operator
 */
START::Band &START::Band::operator=(START::Band const &BandCopy)
{

  if(this != &BandCopy) {
    fTelCode=BandCopy.fTelCode;
    fzenon=BandCopy.fzenon;
    fzenoff=BandCopy.fzenoff;
    foffset=BandCopy.foffset;
    fazimuth = BandCopy.fazimuth;
    falpharun=BandCopy.falpharun;
    fefficiency=BandCopy.fefficiency;
    flivetime=BandCopy.flivetime;
    flivetimefraction=BandCopy.flivetimefraction;
    fnbrun=BandCopy.fnbrun;
    fkeepband=BandCopy.fkeepband;
    frunstarttime=BandCopy.frunstarttime;
    frunendtime=BandCopy.frunendtime;

    fVectorEnergy=BandCopy.fVectorEnergy;
    fVectorArea=BandCopy.fVectorArea;
    fVectorResolution=BandCopy.fVectorResolution;
    fVectorBiais=BandCopy.fVectorBiais;
    fVectorInterEnergy=BandCopy.fVectorInterEnergy;
    fVectorInterArea=BandCopy.fVectorInterArea;
    fVectorInterResolution=BandCopy.fVectorInterResolution;
    fVectorInterBiais=BandCopy.fVectorInterBiais;

    fEthMC=BandCopy.fEthMC;
    fFirstEmcBinNumber=BandCopy.fFirstEmcBinNumber;
    fFirstEmcValFit=BandCopy.fFirstEmcValFit; // VIM
    fFirstEmcValDistrib=BandCopy.fFirstEmcValDistrib; // VIM
    fUseInstrumentEnergyDistribution=BandCopy.fUseInstrumentEnergyDistribution;
    fInterpolTypeForBand = BandCopy.fInterpolTypeForBand;

    fInterpolArea=0; // JLK
    fInterpolBiais=0; // JLK
    fInterpolResol=0; // JLK

    ebin.clear();
    for(unsigned int ibin(0); ibin<BandCopy.ebin.size(); ibin++) ebin.push_back(BandCopy.ebin[ibin]);

    if (fVectorArea.size()>4 && fVectorEnergy.size()>4) {
      fInterpolArea  = new ROOT::Math::Interpolator(fVectorEnergy,fVectorArea,fInterpolTypeForBand);
    }
    else {
      fInterpolArea  = 0;
    }
    
    if (fVectorResolution.size()>4 && fVectorEnergy.size()>4) {
      fInterpolResol  = new ROOT::Math::Interpolator(fVectorEnergy,fVectorResolution,fInterpolTypeForBand);
    }
    else {
      fInterpolResol  = 0;
    }
    
    if (fVectorBiais.size()>4 && fVectorEnergy.size()>4) {
      fInterpolBiais  = new ROOT::Math::Interpolator(fVectorEnergy,fVectorBiais,fInterpolTypeForBand);
    }
    else {
      fInterpolBiais  = 0;
    }
  
    fDistributionInstrumentMapTable=BandCopy.fDistributionInstrumentMapTable;
    InitDistributionInterpTable();

  }
  return (*this);
}

/**
 * \brief Function to Init all the interpolators at one. Can be usefull to rebuild them in case you get a Band from a file (because the interpolator are not streamed)
 **/
void START::Band::InitInterpolator() {
  SetGSLInterpolatorForArea(fVectorEnergy,fVectorArea);
  SetGSLInterpolatorForResolution(fVectorEnergy,fVectorResolution);
  SetGSLInterpolatorForBiais(fVectorEnergy,fVectorBiais);
  InitDistributionInterpTable();
}


void START::Band::Print(Option_t *) const {

  PrintBand();

}

/**
 * \brief Print caracteristics of the band
 *
 * Print : 
 * <ul> 
 *    <li> mean offset </li>
 *    <li> mean zenithal angle for ON events</li>
 *    <li> mean zenithal angle for OFF events</li>
 *    <li> Efficiency </li>
 *    <li> ON Time of observation </li>
 *    <li> OFF Time of observation </li>
 *    <li> \f$\alpha\f$ defined as the ratio of \f$T_{ON}\f$ and \f$T_{OFF}\f$ </li>
 *    <li> Telcode defined as the number of telescopes participating in the run/band </li>
 *    <li> Energy threshold of the band (first bin participating) </li>
 * </ul>
 */
void START::Band::PrintBand(std::ostream &os) const {

  //std::ostringstream ossbandparam;

  int integerprecision(0);
  int doubleprecision(2);
  int tripleprecision(3);
  std::ostringstream bandheader;
  bandheader.precision(integerprecision);
  bandheader << fkeepband << " | Run " << fnbrun;
  bandheader.precision(doubleprecision);
  bandheader.setf(std::ios::fixed);
  bandheader.precision(integerprecision);
  if(fTelCode>0) bandheader << " TelCode " << fTelCode;
  bandheader.precision(doubleprecision);
  bandheader << " RunStart(MJD) " << frunstarttime;
  bandheader << " RunEnd(MJD) " << frunendtime;
  bandheader << " LiveTime(h) " << flivetime;
  bandheader.precision(tripleprecision);
  bandheader << " Alpha " << falpharun;
  bandheader.precision(doubleprecision);
  bandheader << " Efficiency(%) " << fefficiency;
  bandheader << " Offset(deg) " << foffset;
  bandheader << " Zenith(deg) " << fzenon;
  bandheader << std::endl;
  bandheader << "   ";
  bandheader << " TotalOn " << GetNOnTot();
  bandheader << " TotalOff " << GetNOffTot();
  bandheader << " TotalExcess " << GetExcessTot();
  bandheader << " Significance " << GetSignificanceTot();
  bandheader << " EthMC(TeV) " << fEthMC;
  bandheader << " LiveTimeFraction(%) " << flivetimefraction;
  os << bandheader.str().c_str() << std::endl;

  std::ostringstream binheader;
  binheader.width(1);
  binheader << "";
  binheader.width(10);
  binheader << "Emin";
  binheader.width(10);
  binheader << "Emean";
  binheader.width(10);
  binheader << "Emax";
  binheader.width(9);
  binheader << "ON";
  binheader.width(11);
  binheader << "OFF";
  binheader.width(10);
  binheader << "Alpha";
  binheader.width(11);
  binheader << "Excess";
  binheader.width(9);
  binheader << "Sigma";
  binheader.width(12);
  binheader << "LiveTime";
  binheader.width(9);
  binheader << "EffArea";
  binheader.width(9);
  binheader << "ExpON";
  binheader.width(11);
  binheader << "ExpOFF";
  binheader.width(12);
  binheader << "ExpExcess";

  os << binheader.str().c_str() << std::endl;

  for(std::vector<EnergyBin>::const_iterator it = ebin.begin(); it != ebin.end(); ++it) {
    it->PrintEnergyBin(os);
  }

}

/**
 * \brief Return interpolated area for energy E
 */
double START::Band::GetInterpolatedArea(double E)
{
  if(fInterpolArea==0) {
    std::cout << "You have to initialize Interpolators before use them!!!! ==> EXIT" << std::endl;
    exit(EXIT_FAILURE);
  }
  Double_t area = fInterpolArea->Eval(E);
  return ( (area>=0) ? area : 0.);
}

/**
 * \brief Return interpolated biais for energy E
 */
double START::Band::GetInterpolatedBiais(double E)
{
  if(fInterpolBiais==0) {
    std::cout << "You have to initialize Interpolators before use them!!!! ==> EXIT" << std::endl;
    exit(EXIT_FAILURE);
  }
  return (fInterpolBiais->Eval(E));
}

/**
 * \brief Return interpolated resolution for energy E
 */
double START::Band::GetInterpolatedResolution(double E)
{
  if(fInterpolResol==0) {
    std::cout << "You have to initialize Interpolators before use them!!!! ==> EXIT" << std::endl;
    exit(EXIT_FAILURE);
  }
  return (fInterpolResol->Eval(E));
}

/**
 * \brief Set area interpolator in band
 */
void START::Band::SetGSLInterpolatorForArea(std::vector<double> &x, std::vector<double> &y)
{

  if (fInterpolArea!=0) delete fInterpolArea;
  fInterpolArea=0;

  if(fkeepband==1) {
    
    if(x.size()>4 && y.size()>4) {
      fInterpolArea  = new ROOT::Math::Interpolator(x,y,fInterpolTypeForBand);
    }
    else {
      //std::cout << "VIM : Je ne devrais pas etre la, c'est pas bon" << std::endl;
      delete fInterpolArea;
      fInterpolArea = 0;
      std::cout << "Failed to set GSL interpolator for area!!!" << std::endl;
    }

  }

}

/**
 * \brief Set biais interpolators in band
 */
void START::Band::SetGSLInterpolatorForBiais(std::vector<double> &x, std::vector<double> &y)
{

  if (fInterpolBiais!=0) delete fInterpolBiais;
  fInterpolBiais=0;

  if(fkeepband==1) {

    if(x.size()>4 && y.size()>4) {
      fInterpolBiais  = new ROOT::Math::Interpolator(x,y,fInterpolTypeForBand);
    }
    else {
      ///std::cout << "VIM : Je ne devrais pas etre la, c'est pas bon" << std::endl;
      delete fInterpolBiais;
      fInterpolBiais = 0;
      //std::cout << "Failed to set GSL interpolator for biais!!!" << std::endl;
    }

  }

}

/**
 * \brief Set resolution interpolators in band
 */
void START::Band::SetGSLInterpolatorForResolution(std::vector<double> &x, std::vector<double> &y)
{

  if (fInterpolResol!=0) delete fInterpolResol;
  fInterpolResol=0;

  if(fkeepband==1) {

    if(x.size()>4 && y.size()>4) {
      fInterpolResol  = new ROOT::Math::Interpolator(x,y,fInterpolTypeForBand);
    }
    else {
      //std::cout << "VIM : Je ne devrais pas etre la, c'est pas bon" << std::endl;
      delete fInterpolResol;
      fInterpolResol = 0;
      //std::cout << "Failed to set GSL interpolator for Resol!!! " << std::endl;
    }

  }

}

/**
 * \brief Function to init the Interpolator for the EnergyResolution taken from the Distributions.
 * This should be called after SetDistributionVectorTable
 */
void START::Band::InitDistributionInterpTable() {
  
  if (fDistributionInstrumentInterpTable.size()!=0) {
    std::cout << "VIM DEBUG> Aha... Strange, it should be empty at this point... We're going to reset the map" << std::endl;
    DeleteInstrumentInterpTable();
  }
  
  for (std::map<double, std::pair< std::vector<double>, std::vector<double> > >::iterator it = fDistributionInstrumentMapTable.begin(); it!= fDistributionInstrumentMapTable.end();++it) {
    ROOT::Math::Interpolator *interp = new ROOT::Math::Interpolator( (it->second).first,(it->second).second, fInterpolTypeForBand);
    //ROOT::Math::Interpolator *interp = new ROOT::Math::Interpolator( (it->second).first,(it->second).second, ROOT::Math::Interpolation::kLINEAR);
    fDistributionInstrumentInterpTable[it->first] = interp;
  }
}

/**
 * \brief Function to delete all the interpolators created in fDistributionInstrumentInterpTable
 */
void START::Band::DeleteInstrumentInterpTable() {
  
  for (std::map<double, ROOT::Math::Interpolator*>::iterator it = fDistributionInstrumentInterpTable.begin();
       it!=fDistributionInstrumentInterpTable.end();++it) {
    if (it->second) delete it->second;
    it->second = 0;
  }
  fDistributionInstrumentInterpTable.clear();
}

/**
 * \brief Return PDF from distribution. The interpolation between the distrbution around the wanted Etrue is linear.
 *
 * The Element that has created the fDistributionInstrumentInterpTable and fDistributionInstrumentMapTable are supposed to be normalized
 *
 * \param Ereco Reconstructed Energy
 * \param Etrue Reconstructed Energy
 * \return PDF of the resolution (i.e : pobability to reconstruct Ereco knowing Etrue)
 *
 * \todo Check the procedure with fresh eyes
 */
double START::Band::GetEResolProbabilityFromDistribution(double Ereco, double Etrue) {
  
  double e_logratio = TMath::Log(Ereco/Etrue);
  
  if (Etrue< (fDistributionInstrumentInterpTable.begin())->first) {
    //std::cout << "Etrue (" << Etrue << ") is less than the energy of the first bin inside the region : " << (fDistributionInstrumentInterpTable.begin())->first << " PROBLEM!!!!!!!!!!!" << std::endl;
    return 0.0;
  }
  
  if (fDistributionInstrumentInterpTable.size()==0) {
    if (fDistributionInstrumentMapTable.size()!=0) {
      InitDistributionInterpTable();
    }
    else {
      std::cout << "Can't get fDistributionInstrumentInterpTable, It does not exist nor is initialized --> Please Check!!!" << std::endl;
      return 0.0;
    }
  }

  std::vector<double> &first_vec_map_table = ((fDistributionInstrumentMapTable.begin())->second).first;
  if (e_logratio < first_vec_map_table[0] || e_logratio> (first_vec_map_table[first_vec_map_table.size()-1])) {
    //VIM : If we are here, this mean that the e_logratio request is less or above the definition used 
    //for the interpolations. So, the value is 0.0, no need to go further.
    return 0.;
  }
  
  std::map<double,ROOT::Math::Interpolator*>::const_iterator it_distrib;
  std::map<double,ROOT::Math::Interpolator*>::const_iterator it_distrib_next;
  std::map<double,ROOT::Math::Interpolator*>::iterator it_distrib_end = fDistributionInstrumentInterpTable.end();
  
  it_distrib      = fDistributionInstrumentInterpTable.begin();
  it_distrib_next = fDistributionInstrumentInterpTable.begin();
  ++it_distrib_next;
  double e0 = it_distrib->first;
  double e1 = it_distrib_next->first;
  //std::cout << "AVANT while : e0 = " << e0 << " e1 = " << e1 << " Etrue = " << Etrue << " Ereco = " << Ereco << " e_logratio = " << e_logratio << std::endl;
  Bool_t BinFound = false;
  while (!BinFound && it_distrib_next!=it_distrib_end) {
    if (e0<=Etrue && Etrue<e1) {
      //break;
      BinFound = true;
    }
    else {
      it_distrib = it_distrib_next;
      ++it_distrib_next;
      e0 = e1;
      e1 = it_distrib_next->first;
    }
    //std::cout << "Dans While : e0 = " << e0 << " e1 = " << e1 << std::endl;
  }

  // std::cout << "TOTOR Etrue = " << Etrue << " e0 = " << e0 << " e1 = " << e1 << " it_distrib->first = " << it_distrib->first << " it_distrib_next->first = " << it_distrib_next->first << std::endl;
  double prob=0.0;
  
  if (it_distrib_next==fDistributionInstrumentInterpTable.end()) {
    prob = (it_distrib->second)->Eval(e_logratio);
  }
  else {
    double prob1             = (it_distrib->second)->Eval(e_logratio);
    double prob2             = (it_distrib_next->second)->Eval(e_logratio);
    double weight1           = it_distrib_next->first - Etrue;
    double weight2           = Etrue - it_distrib->first;
    double norm_for_weight   = it_distrib_next->first - it_distrib->first;
    // Original Line :
    // double prob1             = (it_distrib->second)->Eval(e_logratio);
    // double prob2             = (it_distrib_next->second)->Eval(e_logratio);
    // double norm_for_weight   = it_distrib_next->first - it_distrib->first;
    // double weight1           = Ereco - it_distrib->first;
    // double weight2           = it_distrib_next->first - Ereco;
    //double prob = (weight1*prob1 + weight2*prob2)/norm_for_weight;
    prob = (weight1*prob1 + weight2*prob2)/norm_for_weight; // JLK
    //std::cout << "VIM ProbEncadre> prob = " << prob << " prob1 = " << prob1 << " prob2 = " << prob2 << " norm_for_weight = " << norm_for_weight << " weight1 = " << weight1 << " weight2 = " << weight2 << " Ereco = " << Ereco << " it->first = " << it_distrib->first << " it_next->first = " << it_distrib_next->first << " e1 = " << e_logratio1 << " e2 = " << e_logratio2 << std::endl;
  }

  return ( (prob>=0.) ? prob : 0. );   
}
  
/**
 * \brief Return PDF from the fit of the resolution distribution
 * \todo Improve Doc
 */
double START::Band::GetEResolProbabilityFromFit(double Ereco, double Etrue) {
  double e_logratio = TMath::Log(Ereco/Etrue);
  
  double mean;
  double sigma;
  double xresol = 0.;
  if (Etrue<GetVectorEnergy()[GetVectorEnergy().size()-1]) {
    mean = GetInterpolatedBiais(Etrue);
    sigma = GetInterpolatedResolution(Etrue);
  }
  else { //SP: avoid a time consuming warning from Interpolator:
    mean  = GetVectorBiais()[GetVectorBiais().size()-1];
    sigma = GetVectorResolution()[GetVectorResolution().size()-1];
    
    if(xresol<0.) {
      std::cout << "ComputeResult::FunctionResolution : interpolated sigma < 0! ==> EXIT!" << std::endl;
      std::cout << "Bande : " << GetNbRun() << std::endl;
      std::cout << "Etrue = " << Etrue << " & Ereco = " << Ereco << std::endl;
      std::cout << "interpolated resolution = " << sigma << " < 0" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  
  xresol = (1. / (TMath::Sqrt(2.*TMath::Pi())*sigma))*TMath::Exp(-0.5* TMath::Power((e_logratio-mean)/sigma,2));
  //return xresol;
  return ( (std::isnan(xresol)) ? 0. : xresol );
}

/**
 * \brief Function that return the probability to reconstruct the energy at Ereco from a true energy Etrue
 * (Select between the result of the fit or the distribution in function of the class member fUseInstrumentEnergyDistribution)
 *
 * @param Ereco Reconstruct Energy
 * @param Etrue True Energy
 * @return Probability to have Ereco knowing Etrue
 */
double START::Band::GetEResolProbability(double Ereco, double Etrue) {  
  
  return ( fUseInstrumentEnergyDistribution ? GetEResolProbabilityFromDistribution(Ereco,Etrue) : GetEResolProbabilityFromFit(Ereco,Etrue) );
  
  // if (fUseInstrumentEnergyDistribution) {
  //   return GetEResolProbabilityFromDistribution(Ereco,Etrue);
  // }
  // else {
  //   return GetEResolProbabilityFromFit(Ereco,Etrue);
  // }
}



/**
 * \brief Function to remove all the informations from a Band (usefull if you want to make a summarized band)
 *
 */
void START::Band::ClearBandInfo() {

  for (std::vector<EnergyBin>::iterator ibin = ebin.begin();ibin!=ebin.end();++ibin) {
    ibin->ClearEnergyBinInfo();
  }

  fnbrun    = 0;
  fTelCode  = 0;
  fzenon    = 0;
  fzenoff   = 0;
  foffset   = 0;
  fazimuth  = 0.;
  falpharun = 0;
  fefficiency    = 0;
  flivetime     = 0;
  flivetimefraction = 0.;
  fkeepband  = 1;
  frunstarttime = 0.;
  frunendtime = 0.;
  /*
  //VIM I don't really know what to do with them
  std::vector<double> fVectorEnergy;
  std::vector<double> fVectorArea;
  std::vector<double> fVectorResolution;
  std::vector<double> fVectorBiais;
  std::vector<double> fVectorInterEnergy;
  std::vector<double> fVectorInterArea;
  std::vector<double> fVectorInterResolution;
  std::vector<double> fVectorInterBiais;
  std::map< double, std::pair< std::vector<double>, std::vector<double> > > fDistributionInstrumentMapTable;
  std::map< double, ROOT::Math::Interpolator* > fDistributionInstrumentInterpTable;
  */

  fEthMC=0; // energy threshold (to be defined)
  fFirstEmcBinNumber=0; /* SP: bin number in vector MonteCarlo::GetEnergy() for which 
			   area, resol and bias are available.
			   In the estimation of the expected number of gammas, the
			   integral on true energy must avoid going below this energy.
			*/
  fFirstEmcValFit=0; //corresponding MC energy value     //VIM
  fFirstEmcValDistrib=0; //corresponding MC energy value //VIM

  /*
  // VIM : I don't know what to do with those variables : 
  //bool fUseInstrumentEnergyDistribution; // VIM : In order to tell to the program to use the distribution
  // Interpolator of area, resolution and biais
  ROOT::Math::Interpolator *fInterpolArea; 
  ROOT::Math::Interpolator *fInterpolBiais;
  ROOT::Math::Interpolator *fInterpolResol;
  */

  
}



/**
 * \brief Routine to add information from a vector of Band (basically if you want to summarize the information)
 * A sanity check is performed in order to know if the Energybin size is the same
 * The information from the current band is also added (but safer to reset it first). They are weighted by obeservation time.
 *
 * The average on Alpha is made by averaging over noff counts, and use the averaging over livetime as a backup solution when no off.
 * This last case is biased. So another solution should be found !
 *
 * \param std::vector<Band> list of band to add
 * \param IsUseEthresh boolean if true, then information will be add only if the bin have the GetKeepBin()=1 (above Ethreshold), default = false
 *
 * \todo Add frunstarttime and frunendtime 
 *
 */
void START::Band::AddInfoFromBands(const std::vector<Band> &VecBand, bool IsUseEthresh) {
  
  // VIM : Some sanity checks : 
  
  if(VecBand.size()==0) {
    WARN_OUT << "You have to give a vector that actually contains something" << std::endl;
    return;
  }

  // VIM : Check on the number of energy bin
  bool is_same_number_of_energybin = true;
  int badbandnumber=-1;

  for (unsigned int iband=0;iband<VecBand.size();++iband) {
    if (VecBand[iband].GetKeepBand()==0) continue;  // We are interested only in the selected bands
    if (VecBand[iband].ebin.size()!=ebin.size()) {
      is_same_number_of_energybin = false;
      badbandnumber = iband;
    }
  }

  // VIM : Check on the energyrange of energy bin
  if (!is_same_number_of_energybin) {
    WARN_OUT << "The vector entry " << badbandnumber << " have " << VecBand[badbandnumber].ebin.size() 
	     << " energy bin, compared to " << ebin.size() << " in the local (I only signalled the first bad entry)" << std::endl;
    return;
  }
  
  bool is_same_binsize = true;
  int badbinnumber=-1;

  // It could be good to rewrite this condition in a more cleaner way
  for (unsigned int iband=0;iband<VecBand.size();++iband) {
    if (VecBand[iband].GetKeepBand()==0) continue;  // We are interested only in the selected bands
    for (unsigned int ibin=0;ibin<VecBand[iband].ebin.size();++ibin) {
      if ((ebin[ibin].GetEmin() != VecBand[iband].ebin[ibin].GetEmin()) || (ebin[ibin].GetEmax() != VecBand[iband].ebin[ibin].GetEmax())) {
	is_same_binsize = false;
	badbandnumber = iband;
	badbinnumber = ibin;
      }
    }
  }
  
  if (!is_same_binsize) {
    WARN_OUT << "The Energy bin number : " << badbinnumber << " of the band number : " << badbandnumber 
	     << " does not have the same size than the local one " << std::endl;
    return;
  }
  
  // ********* VIM : Ok, everything is correct, let's add !


  // Computation of the global parameter for the summary :

  int nrun_add = 0;
  double alpha_total_mean =0.;
  double alpha_total_mean_backup = 0.;
  double total_time=0.;
  double total_noff=0.;
  double mean_zen_on=0.;
  double mean_zen_off=0;
  double meanefficiency=0;
  double meanoffset=0;
  double meanehtresh=0.;

  alpha_total_mean = GetAlphaRun()*GetNOffTot(IsUseEthresh);
  alpha_total_mean_backup = GetAlphaRun()*GetLiveTime();

  total_noff+=GetNOffTot(IsUseEthresh);
  total_time+=GetLiveTime();

  mean_zen_on+=(GetZenON()*GetLiveTime());
  mean_zen_off+=(GetZenOFF()*GetLiveTime());
  meanefficiency+=(GetEff()*GetLiveTime());
  meanoffset+=(GetOffset()*GetLiveTime());
  meanehtresh+=(GetEthMC()*GetLiveTime());
  

  if (GetNOnTot(IsUseEthresh)!=0 || GetNOffTot(IsUseEthresh)!=0) { // In this case, we assume that this band is filled, and must be count
    ++nrun_add;
  }
  
  
  for (std::vector<Band>::const_iterator itband = VecBand.begin();itband!=VecBand.end();++itband) {
    if (itband->GetKeepBand()==0) {
      continue;
    }

    ++nrun_add;
    double noff = itband->GetNOffTot(IsUseEthresh);
    total_noff+=noff;
    alpha_total_mean += itband->GetAlphaRun()*noff;
    alpha_total_mean_backup += itband->GetAlphaRun()*itband->GetLiveTime();
    total_time+=itband->GetLiveTime();
    
    mean_zen_on+=(itband->GetZenON()*itband->GetLiveTime());
    mean_zen_off+=(itband->GetZenOFF()*itband->GetLiveTime());
    meanefficiency+=(itband->GetEff()*itband->GetLiveTime());
    meanoffset+=(itband->GetOffset()*itband->GetLiveTime());
    meanehtresh+=(itband->GetEthMC()*itband->GetLiveTime());
  }
  
  if (total_noff) {alpha_total_mean=alpha_total_mean/total_noff;}
  if (total_time) {alpha_total_mean_backup=alpha_total_mean_backup/total_time;}
  if (alpha_total_mean==0.) {alpha_total_mean = alpha_total_mean_backup;}
  
  if (total_time) {
    mean_zen_on=mean_zen_on/total_time;
    mean_zen_off=mean_zen_off/total_time;
    meanefficiency=meanefficiency/total_time;
    meanoffset=meanoffset/total_time;
    meanehtresh=meanehtresh/total_time;
  }
  
  // We fill the band : BE  CAREFULL, AFTER THOSE LINE, THE MEMBER OF THIS CLASS WILL BE CHANGE !!!!
  SetKeepBand(0); // To be certain that this band won't be used
  SetTelCode(0); // To be certain that this band won't be used
  SetZenON(mean_zen_on);
  SetZenOFF(mean_zen_off);
  SetOffset(meanoffset);
  SetEff(meanefficiency);
  SetLiveTime(total_time);
  SetNbRun(nrun_add);
  SetAlphaRun(alpha_total_mean); // VIM Be carefull to not call bellow the value GetAlphaRun(), we have modify it with the average value...
  if (IsUseEthresh) {
    SetEthMC(meanehtresh);
  }
  else {
    SetEthMC(0.);
  }

  // VIM Now we're going to stack the different energy bin
  
  for (unsigned int ibin = 0;ibin<ebin.size();++ibin) {

    EnergyBin &ebin_local = ebin[ibin]; // Reference to avoid calling to this each time (should be faster in time, and I am lazy)
    
    // VIM : In this case we need to handle the treatment of alpha outside the AddInfoFromEBin, 
    // because it doesn't work properly when adding results from 2 bin with low statistics, which is the case here.
    double noffbin_tot = 0.;
    double alphabin_total_mean = 0.;
    
    for (std::vector<Band>::const_iterator itband = VecBand.begin();itband!=VecBand.end();++itband) {
      
      if (itband->GetKeepBand()==0) continue;  // We are interested only in the selected bands
      
      const EnergyBin &ebin_vector = itband->ebin[ibin]; 

      if (!IsUseEthresh || ebin_vector.GetKeepBin()==1) {
	noffbin_tot += ebin_vector.GetOff();
	alphabin_total_mean += ebin_vector.GetAlpha()*ebin_vector.GetOff();
	ebin_local.AddInfoFromEBin(ebin_vector);

	if (ebin_local.GetKeepBin()==1 || ebin_vector.GetKeepBin()==1) {
	  ebin_local.SetKeepBin(1);
	}	
      }
    }

    if (noffbin_tot!=0.) {
      alphabin_total_mean = alphabin_total_mean/noffbin_tot;
      ebin_local.SetAlpha(alphabin_total_mean);
    }
    // No Need for backup solution because it is already implemented inside AddInfoFromEBin
  }
  
  // JLK : add start and end of the run

  std::vector<double> startsrun;
  std::vector<double> endsrun;
  for(std::vector<Band>::const_iterator band=VecBand.begin(); band!=VecBand.end(); ++band) {
    startsrun.push_back(band->GetRunStartTime());
    endsrun.push_back(band->GetRunEndTime());
  }
  
  std::vector<double>::iterator min,max;
  min=min_element(startsrun.begin(),startsrun.end());
  max=max_element(endsrun.begin(),endsrun.end());

  SetRunStartTime(*min);
  SetRunEndTime(*max);

}


/**
 * \brief Function that will retrieve the Total number of On event inside the band
 *
 * @param IsUseEthresh boolean to say if you want the total number of event above Ethreshold
 *
 * \todo Maybe change also those function in order to interate in a energy range (from x to y TeV ? )
 */
double START::Band::GetNOnTot(bool IsUseEthresh) const {
  double non_tot = 0.0;
  for (std::vector<EnergyBin>::const_iterator itbin=ebin.begin();itbin!=ebin.end();++itbin) {
    if (!IsUseEthresh || itbin->GetKeepBin()==1) {
      non_tot+=itbin->GetOn();
    }
  }
  return non_tot;
}


/**
 * \brief Function that will retrieve the Total number of Off event inside the band
 *
 * @param IsUseEthresh boolean to say if you want the total number of event above Ethreshold
 *
 * \todo Maybe change also those function in order to interate in a energy range (from x to y TeV ? )
 */
double START::Band::GetNOffTot(bool IsUseEthresh) const {
  double noff_tot = 0.0;
  for (std::vector<EnergyBin>::const_iterator itbin=ebin.begin();itbin!=ebin.end();++itbin) {
    if (!IsUseEthresh || itbin->GetKeepBin()==1) {
      noff_tot+=itbin->GetOff();
    }
  }  
  return noff_tot;
}

/**
 * \brief Function that will retrieve the Total Excess inside the band
 *
 * @param IsUseEthresh boolean to say if you want the total number of event above Ethreshold
 *
 * \todo Maybe change also those function in order to interate in a energy range (from x to y TeV ? )
 */
double START::Band::GetExcessTot(bool IsUseEthresh) const {
  double nexcess_tot = GetNOnTot(IsUseEthresh)-(GetAlphaRun()*GetNOffTot(IsUseEthresh));
  return nexcess_tot;
}

/**
 * \brief Function that will retrieve the Significance of the band
 *
 * @param IsUseEthresh boolean to say if you want the total number of event above Ethreshold
 *
 * \todo Maybe change also those function in order to interate in a energy range (from x to y TeV ? )
 */
double START::Band::GetSignificanceTot(bool IsUseEthresh) const {
  double significance = STARTUtils::LiMaSignificance(TMath::Nint(GetNOnTot(IsUseEthresh)),TMath::Nint(GetNOffTot(IsUseEthresh)),GetAlphaRun());
  return significance;
}

/**
 * \brief Function that will retrieve the Total number of expected signal event inside the band
 *
 * @param IsUseEthresh boolean to say if you want the total number of event above Ethreshold
 *
 * \todo Actually it is only filled above te Ethreshold (in MinimizeFactory.C), when te Fit Converged. Maybe change that with an extrapolation
 * \todo Maybe change also those function in order to interate in a energy range (from x to y TeV ? )
 */
double START::Band::GetNSthTot(bool IsUseEthresh) const {
  double nsth_tot = 0.0;
  for (std::vector<EnergyBin>::const_iterator itbin=ebin.begin();itbin!=ebin.end();++itbin) {
    if (!IsUseEthresh || itbin->GetKeepBin()==1) {
      nsth_tot+=itbin->GetSth();
    }
  }
  return nsth_tot;
}

/**
 * \brief Function that will retrieve the Total number of expected backround event inside the band
 *
 * @param IsUseEthresh boolean to say if you want the total number of event above Ethreshold
 *
 * \todo Actually it is only filled above te Ethreshold (in MinimizeFactory.C), when te Fit Converged. Maybe change that with an extrapolation
 * \todo Maybe change also those function in order to interate in a energy range (from x to y TeV ? )
 */
double START::Band::GetNOffFitted(bool IsUseEthresh) const {
  double nofff_tot = 0.0;
  for (std::vector<EnergyBin>::const_iterator itbin=ebin.begin();itbin!=ebin.end();++itbin) {
    if (!IsUseEthresh || itbin->GetKeepBin()==1) {
      nofff_tot+=itbin->GetOffFitted();
    }
  }
  return nofff_tot;
}



/**
 * \brief Procedure to initialise the alpha value inside each bin
 * The interest of having the alpha value in the bin is when you want to stack the data, then alpha
 * might be different from bin to bin (depending on the threshold energy
 * Do not forget this function when you change the alpha value
 *
 * @param alpharun The alpha value of the run.
 *
 **/
void START::Band::InitEbinAlpha(double alpharun) {
  for (std::vector<EnergyBin>::iterator itbin = ebin.begin();itbin!=ebin.end();++itbin) {
    itbin->SetAlpha(alpharun);
  }
}



/**
 * \brief Procedure to initialise the livetime value inside each bin
 * The interest of having the livetime value in the bin is when you want to stack the data, then livetime
 * might be different from bin to bin (depending on the threshold energy)
 * Do not forget this function when you change the alpha value
 *
 * @param livetime The time value of the run.
 *
 **/
void START::Band::InitEbinLiveTime(double livetime) {
  for (std::vector<EnergyBin>::iterator itbin = ebin.begin();itbin!=ebin.end();++itbin) {
    itbin->SetLiveTime(livetime);
  }
}
  

  
/**
 * \brief Procedure to initialise the energy bins
 * This function is usefull in case you want to quickly initialise the energy bins of this band;
 * The bin size is done in logarithmic scale
 *
 * @param emin Minimum energy of the first bin
 * @param emax Maximum energy of the last bin
 * @param Nbins Number of bins (logarithmic scale)
 *
 **/

void START::Band::InitEnergyBins(Double_t emin, Double_t emax, Int_t Nbins)
{
  if (ebin.size()!=0) {
    std::cout << "EnergyBins already initiated !" << std::endl;
    return;
  }

  Double_t logemin  = TMath::Log10(emin);
  Double_t logemax  = TMath::Log10(emax);
  Double_t binwidth = (logemax - logemin)/Double_t(Nbins);
  
  for(Int_t enbin=0;enbin<Nbins;enbin++) {
    Double_t ebinmin = TMath::Power(10.,logemin+Double_t(enbin)*binwidth);
    Double_t ebinmax = TMath::Power(10.,logemin+Double_t(enbin+1)*binwidth);
    ebin.push_back( EnergyBin(ebinmin,ebinmax) );
  }

}


/**
 * \brief Function that return a pointer to the EnergyBin that contains a given energy
 *
 * @param ener Energy for which we search the EnergyBin object that contains it.
 *
 * A EnergyBin is considered as containing the use energy when GetEmin() <= ener < GetEmax()
 *
 */

START::EnergyBin* START::Band::GetEnergyBin(Double_t ener) {

  std::vector<EnergyBin>::iterator it = ebin.begin();
  std::vector<EnergyBin>::iterator itend = ebin.end();
  Bool_t foundbin = false;

  while(!foundbin && it!=itend) {
    if ( it->GetEmin()<=ener && ener<it->GetEmax() ) {
      foundbin = true;
    }
    else {
      ++it;
    }
  }

  if (!foundbin) {
    std::cout << "\033[31m" << "Can't find an energybin that contains " << ener << " TeV !" << "\033[0m" << std::endl;
  }
  
  return (foundbin ?  &(*it) : 0);
}



/**
 * \brief Construct and return the name of the Band
 *
 * \return Return the name of the band
 *
 **/
std::string START::Band::ConstructBandName() const {
  return ConstructBandNameOrTitle(true);
}

/**
 * \brief Construct and return the title of the Band
 *
 * \return Return the title of the band
 *
 **/
std::string START::Band::ConstructBandTitle()  const {
  return ConstructBandNameOrTitle(false);
}

/**
 * \brief Update the current name and title of the Band
 * This function actually touch the class, and REALLY modify the TNamed parameters !
 **/
void START::Band::UpdateNameAndTitle() {
  SetNameTitle( ConstructBandName().c_str(), ConstructBandTitle().c_str() );
}

/**
 * \brief Construct the name of title for the current band (If you want that name/title to be stored call UpdateNameAndTitle() function
 *
 * \param constructname Boolean to tell if the name has to be construct
 *
 * \return Return the Name if constructname is true and the Title otherwise
 *
 **/
std::string START::Band::ConstructBandNameOrTitle(Bool_t constructname) const {
  
  std::ostringstream bandname;
  bandname << "Run";
  bandname.precision(0);
  bandname << fnbrun;
  if (constructname) { bandname << "_eff";}
  else {bandname << ", eff=";}
  bandname.precision(2);
  bandname << fefficiency;
  if (constructname) { bandname << "_zen";}
  else {bandname << ", zen=";}
  bandname.precision(2);
  bandname << fzenon;
  if (constructname) { bandname << "_off";}
  else {bandname << ", off=";}
  bandname.precision(2);
  bandname << foffset;
  if (constructname) { bandname << "_TelCode";}
  else {bandname << ", TelCode=";}
  bandname.precision(0);
  bandname << fTelCode;
  // bandname.precision(2);
  // bandname << falpharun;
  return bandname.str();
}








/**
 * \brief Procedure that display the Area for this band
 *
 **/

void START::Band::DisplayArea(Bool_t save) {

  std::string can_areaname = "canArea_"+ConstructBandName();
  
  TCanvas *can_area = new TCanvas(can_areaname.c_str(),can_areaname.c_str());
  can_area->SetLogy();
  can_area->SetLogx();
  TGraph *gr_area=0;
  TGraph *gr_interparea=0;

  if (fVectorEnergy.size()!=0 && fVectorArea.size()!=0 && fVectorArea.size()==fVectorEnergy.size()) {
    std::ostringstream oss_grarea_title;
    oss_grarea_title << "Area_Zen" << GetZenON() << "_Off" << GetOffset() << "_Eff"  << GetEff();
    gr_area = new TGraph(fVectorArea.size());
    gr_area->SetNameTitle("FixedEnergy",oss_grarea_title.str().c_str());
    for (unsigned int i = 0 ; i<fVectorArea.size() ; ++i) {
      gr_area->SetPoint(i,fVectorEnergy[i],fVectorArea[i]);
    }
  }
  else {
    std::cout << "fVectorEnergy.size() = " << fVectorEnergy.size() << " fVectorArea.size() = " << fVectorArea.size() << std::endl;
  }

  if (fVectorInterEnergy.size()!=0 && fVectorInterArea.size()!=0 && fVectorInterArea.size()==fVectorInterEnergy.size()) {
    std::ostringstream oss_interpgrarea_title;
    oss_interpgrarea_title << "InterpArea_Zen" << GetZenON() << "_Off" << GetOffset() << "_Eff"  << GetEff();
    gr_interparea = new TGraph(fVectorInterArea.size());
    gr_interparea->SetNameTitle("Interpolated",oss_interpgrarea_title.str().c_str());
    for (unsigned int i = 0 ; i<fVectorInterArea.size() ; ++i) {
      gr_interparea->SetPoint(i,fVectorInterEnergy[i],fVectorInterArea[i]);
    }
  }
  else {
    std::cout << "fVectorInterEnergy.size() = " << fVectorInterEnergy.size() << " fVectorInterArea.size() = " << fVectorInterArea.size() << std::endl;
  }


  // Make root delete those when the canvas is deleted (if I understood correctly) :
  gr_area->SetBit(kCanDelete);
  gr_interparea->SetBit(kCanDelete); 
  
  gr_area->SetMarkerColor(4);
  gr_area->SetMarkerStyle(8);
  gr_area->Draw("AP");
  gr_area->GetHistogram()->GetXaxis()->SetTitle("Energy (TeV)");
  gr_area->GetHistogram()->GetYaxis()->SetTitle("Area (m^{2})");
    
  if (gr_interparea) {
    gr_interparea->SetLineWidth(2);
    gr_interparea->Draw("L SAME");
  }

  TLegend *leg_area = new TLegend(0.57,0.20,0.97,0.47);
  leg_area->SetFillColor(0);
  leg_area->AddEntry(gr_area,gr_area->GetTitle(),"p");
  if (gr_interparea!=0) {
    leg_area->AddEntry(gr_interparea,gr_interparea->GetTitle(),"l");
  }


  if (fEthMC>0.) {
    TLine *ThresholdLine = new TLine( fEthMC, gr_area->GetHistogram()->GetYaxis()->GetBinLowEdge(1), fEthMC, gr_area->GetHistogram()->GetYaxis()->GetBinUpEdge( gr_area->GetHistogram()->GetNbinsY() ) );
    ThresholdLine->SetLineWidth(4);
    ThresholdLine->SetLineColor(2);
    ThresholdLine->Draw("Same");
    ThresholdLine->SetBit(kCanDelete);
    leg_area->AddEntry(ThresholdLine,"Safe Threshold","l");
  }
  leg_area->Draw();

  if (save) {
    can_area->SaveAs(".png");
  }
}

/**
 * \brief Display the energy resolution of this band for a given True Energy
 * \param Energy_True Input true energy at which you want to display the resolution
 **/
void START::Band::DisplayEnergyResolutionForTrueEnergy(Double_t Energy_True) {
  DisplayEnergyResolutionForTrueOrRecoEnergy(Energy_True,true);
}

/**
 * \brief Display the energy resolution of this band for a given Reco Energy
 * \param Energy_Reco Input reco energy at which you want to display the resolution
 **/
void START::Band::DisplayEnergyResolutionForRecoEnergy(Double_t Energy_Reco) {
  DisplayEnergyResolutionForTrueOrRecoEnergy(Energy_Reco,false);
}

/**
 * \brief Display the energy resolution for a given True or Reco Energy
 * \param Energy Input True or Reco energy at which you want to display the resolution
 * \param IsTrueEnergy Boolean True if you input True Energy in the Energy Parameter
 **/
void START::Band::DisplayEnergyResolutionForTrueOrRecoEnergy(Double_t Energy, Bool_t IsTrueEnergy) {
  
  std::ostringstream oss_canname;
  oss_canname << "canResolution_" << (IsTrueEnergy ? "ETrue" : "EReco" ) << Energy << "TeV_" << ConstructBandName();
  std::string can_resolname = oss_canname.str();
  oss_canname.str("");
  
  TCanvas *can_resol = new TCanvas(can_resolname.c_str(),can_resolname.c_str());
  can_resol->cd();
  
  std::ostringstream oss_histname;
  oss_histname << "Resolution_" << (IsTrueEnergy ? "ETrue" : "EReco") << Energy << "TeV_" << ConstructBandName();
  std::ostringstream oss_histtitle;
  oss_histtitle << "Resolution " << (IsTrueEnergy ? "ETrue" : "EReco") << " = " << Energy << " TeV, " << ConstructBandTitle();
  
  
  TH1F *hreso = new TH1F(oss_histname.str().c_str(),oss_histtitle.str().c_str(),300,-3.,3.);
  hreso->GetXaxis()->SetTitle("ln( #frac{E_{Reco}}{E_{True}} )");
  // Make root delete those when the canvas is deleted (if I understood correctly) :
  hreso->SetBit(kCanDelete);
  
  // Creation d'un interpolateur
  std::vector<Double_t> vec_ln_er_o_et;
  std::vector<Double_t> vec_proba;

  for (Int_t ibin=1;ibin<=hreso->GetNbinsX();++ibin) {
    Double_t ln_er_o_et_a = hreso->GetXaxis()->GetBinLowEdge(ibin);
    Double_t ereco_a;
    Double_t prob_a;
    if (IsTrueEnergy) {
      ereco_a = Energy*TMath::Exp(ln_er_o_et_a); // RECONSTRUCTED ENERGY
      prob_a = GetEResolProbability(ereco_a,Energy);
    }
    else {
      ereco_a = Energy*TMath::Exp(-1.*ln_er_o_et_a); // TRUE ENERGY
      prob_a = GetEResolProbability(Energy,ereco_a);
    }
    vec_ln_er_o_et.push_back(ln_er_o_et_a);
    vec_proba.push_back(prob_a);

    if (ibin==hreso->GetNbinsX()) {
      Double_t ln_er_o_et_b = hreso->GetXaxis()->GetBinUpEdge(ibin);
      Double_t ereco_b;
      Double_t prob_b;
      if (IsTrueEnergy) {
	ereco_b = Energy*TMath::Exp(ln_er_o_et_b); // RECONSTRUCTED ENERGY
	prob_b = GetEResolProbability(ereco_b,Energy);
      }
      else {
	ereco_b = Energy*TMath::Exp(-1.*ln_er_o_et_b); // TRUE ENERGY
	prob_b = GetEResolProbability(Energy,ereco_b);
      }
      vec_ln_er_o_et.push_back(ln_er_o_et_b);
      vec_proba.push_back(prob_b);
    }
  }

  //ROOT::Math::Interpolator interp(vec_ln_er_o_et,vec_proba,ROOT::Math::Interpolation::kAKIMA);
  ROOT::Math::Interpolator interp(vec_ln_er_o_et,vec_proba,ROOT::Math::Interpolation::kLINEAR);
  
  for (Int_t ibin=1;ibin<=hreso->GetNbinsX();++ibin) {
    Double_t integ_bin = interp.Integ(hreso->GetXaxis()->GetBinLowEdge(ibin),hreso->GetXaxis()->GetBinUpEdge(ibin));
    hreso->SetBinContent(ibin,integ_bin);
  }

  hreso->Draw();
  
}


/**
 * \brief Display the resolution in a TH2D with ETrue in abscisse and EReco in ordonnee
 *
 **/
void START::Band::DisplayEnergyResolution(Bool_t save) {
  
  std::ostringstream oss_canname;
  oss_canname << "canResolution_" << ConstructBandName();
  std::string can_resolname = oss_canname.str();
  oss_canname.str("");
  
  TCanvas *can_resol = new TCanvas(can_resolname.c_str(),can_resolname.c_str());
  can_resol->SetLogx();
  can_resol->SetLogy();
  can_resol->cd();

  std::ostringstream oss_histname;
  oss_histname << "Resolution_" << ConstructBandName();
  std::ostringstream oss_histtitle;
  oss_histtitle << "Resolution " << ConstructBandTitle();

  std::vector<Double_t> EnergyTrue_bins;
  std::vector<Double_t> EnergyReco_bins;

  Int_t NbinsReco = 300;
  Double_t LogEnergyRecoMin = TMath::Log10(0.03);
  Double_t LogEnergyRecoMax = TMath::Log10(120.);
  Double_t LogERecoStep = (LogEnergyRecoMax-LogEnergyRecoMin)/Double_t(NbinsReco);
  for (Int_t i=0;i<=NbinsReco;++i) {
    Double_t ereco = TMath::Power(10.,LogEnergyRecoMin+LogERecoStep*Double_t(i));
    EnergyReco_bins.push_back(ereco);
  }

  Int_t NbinsTrue = 300;
  Double_t LogEnergyTrueMin = TMath::Log10(0.03);
  Double_t LogEnergyTrueMax = TMath::Log10(120.);
  Double_t LogETrueStep = (LogEnergyTrueMax-LogEnergyTrueMin)/Double_t(NbinsTrue);
  for (Int_t i=0;i<=NbinsTrue;++i) {
    Double_t etrue = TMath::Power(10.,LogEnergyTrueMin+LogETrueStep*Double_t(i));
    EnergyTrue_bins.push_back(etrue);
  }
  
  
  TH2D *hreso = new TH2D(oss_histname.str().c_str(),oss_histtitle.str().c_str(),NbinsTrue,&EnergyTrue_bins[0],NbinsReco,&EnergyReco_bins[0]);
  hreso->SetStats(0);
  hreso->GetXaxis()->SetTitle("E_{True}");
  hreso->GetYaxis()->SetTitle("E_{Reco}");
  hreso->SetBit(kCanDelete);
    
  // On va boucler comme un sagouin sur X, prendre le centre, et ensuite faire la boucle !
  
  for (Int_t ibin = 1;ibin<hreso->GetNbinsX();++ibin) {
    Double_t E_true_cur = hreso->GetXaxis()->GetBinCenter(ibin);
    
    std::vector<Double_t> vec_ln_er_o_et;
    std::vector<Double_t> vec_proba;
    
    for (Int_t jbin = 1;jbin<=hreso->GetNbinsY();++jbin) {
      Double_t E_reco_cur = hreso->GetYaxis()->GetBinLowEdge(jbin);
      vec_ln_er_o_et.push_back( TMath::Log(E_reco_cur/E_true_cur) );
      vec_proba.push_back( GetEResolProbability(E_reco_cur,E_true_cur) );
      if (jbin==hreso->GetNbinsY()) {
	E_reco_cur = hreso->GetYaxis()->GetBinUpEdge(jbin);
	vec_ln_er_o_et.push_back( TMath::Log(E_reco_cur/E_true_cur) );
	vec_proba.push_back( GetEResolProbability(E_reco_cur,E_true_cur) );
      }
    }
    
    //ROOT::Math::Interpolator interp(vec_ln_er_o_et,vec_proba,ROOT::Math::Interpolation::kLINEAR);
    ROOT::Math::Interpolator interp(vec_ln_er_o_et,vec_proba,ROOT::Math::Interpolation::kAKIMA);
    
    for (Int_t jbin=1;jbin<=hreso->GetNbinsY();++jbin) {
      Double_t x_binmin = TMath::Log(hreso->GetYaxis()->GetBinLowEdge(jbin)/E_true_cur);
      Double_t x_binmax = TMath::Log(hreso->GetYaxis()->GetBinUpEdge(jbin)/E_true_cur);
      //std::cout << "x_binmin = " << x_binmin << " x_binmax = " << x_binmax << std::endl;
      Double_t integ_bin = interp.Integ(x_binmin,x_binmax);
      hreso->SetBinContent(ibin,jbin,integ_bin);
    }
  }
  
  hreso->Draw("COLZ");


  if (fEthMC>0.) {
    TLine *ThresholdLine = new TLine( fEthMC, hreso->GetYaxis()->GetBinLowEdge(1), fEthMC, hreso->GetYaxis()->GetBinUpEdge( hreso->GetNbinsY() ) );
    ThresholdLine->SetLineWidth(4);
    ThresholdLine->SetLineColor(2);
    ThresholdLine->Draw("Same");
    ThresholdLine->SetBit(kCanDelete);
  }

if (save) {
    can_resol->SaveAs(".png");
  }
}



// /**
//  * \brief Stack Bands into one single Band
//  *
//  * \brief No IRF available
//  **/
// void START::Band::StackBandFromVectorOfBands(std::vector<Band> vBands) {




// }
