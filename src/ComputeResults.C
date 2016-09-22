// STL
#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>
#include <cstdlib>

// ROOT
#include "Math/Interpolator.h"
#include "Math/Integrator.h"
#include "Math/GSLIntegrator.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/IntegrationTypes.h"
#include "Math/Functor.h"
#include "TF1.h"
#include "TMath.h"

// START
#include "ComputeResults.hh"
#include "GSLError.h"

// Utilities
#define DEBUG 0
#include "debugging.hh"

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::ComputeResults)
#endif

/**
 * \brief Constructor used in FCN
 */
START::ComputeResults::ComputeResults(const std::vector<Band> &SelectedBands, Hypothesis &hypo)
:fHypothesis(0),
  myband_tmp(0),
  mynrj_tmp(0),
  nrjreco_tmp(0)
{

  fBandArray = SelectedBands;
  
  fHypothesis = &hypo;
  
  fabsoluteprecision = 0.;
  frelativeprecision = 0.0001;

  fIntegrationType = ROOT::Math::IntegrationOneDim::kNONADAPTIVE;

  //ROOT::Math::GSLError::GSLError();

}

/**
 * \brief Constructor used to compute integrated resolution
 * on reco energy
 */
START::ComputeResults::ComputeResults(const std::vector<Band> &SelectedBands)
:fHypothesis(0),
  myband_tmp(0),
  mynrj_tmp(0),
  nrjreco_tmp(0)
{

  fBandArray = SelectedBands; 

  fabsoluteprecision = 0.;
  frelativeprecision = 0.000001;

  fIntegrationType = ROOT::Math::IntegrationOneDim::kADAPTIVE;

}

/**
 * \brief Destructor
 */
START::ComputeResults::~ComputeResults()
{

}


/**
 * \brief Return effective Area function of band for energy Etrue
 */
double START::ComputeResults::FunctionEffectiveArea(double Etrue, Band *band)
{

  double xarea(0);
  
  // Make interpolation to get the area
  unsigned int lastbin = band->GetVectorEnergy().size()-1;
  
  if (Etrue<band->GetVectorEnergy()[lastbin]) {
    xarea = band->GetInterpolatedArea(Etrue);
    if(xarea!=xarea) {
      std::cout << "Effective area is NAN!!!" << std::endl;
      std::cout << "Etrue = " << Etrue << " | area = " << xarea << std::endl;      
      exit(EXIT_FAILURE);
    }
  } 
  else {
    //linear interpolation because kAKIMA remains constant after last energy point
    xarea = band->GetVectorArea()[lastbin] + 
      ((Etrue-band->GetVectorEnergy()[lastbin])/(band->GetVectorEnergy()[lastbin]-band->GetVectorEnergy()[lastbin-1]))*(band->GetVectorArea()[lastbin]-band->GetVectorArea()[lastbin-1]);
    if (xarea<0) xarea=0;
  }

#ifdef mydebug
  std::cout << "START::ComputeResults::FunctionEffectiveArea" << std::endl;
  std::cout << "Etrue=" << Etrue << " " << std::endl;
  std::cout << "xarea=" << xarea << std::endl;
#endif

  return xarea*TMath::Power(10,4); // CGS units
}


/**
 * \brief Return Resolution function of band for energy Etrue and energy Ereco
 *
 *SP:
 * We have to consider that:
 * - here Etrue can take values between the limits fixed in FunctionExpectedExcess 
 *   (typiquelly between Ereco/Exp(3) and Ereco*Exp(3)).
 * - but the bias or the sigma of the resolution function are not known (set to 0...)
 *   if Etrue is too low.
 * Proposition:
 * - make the energy vector to be defined only for energies corresponding to
 *   well defined A,B,R. This has to be made at the level of the energy vector
 *   definition for each band.
 * - le the interpolator to make extrapolations if Etrue is below the first value
 *   of the energy vector, but check that the difference between these two values
 *   is not too important. If necessary it is possible to redefine Etrue_min in 
 *   FunctionExpectedExcess
 *
 * VIM : 
 *: I add the treatment of thte distribution and fitted case.
 * In order to improve the readability, I create 2 functions in Band : 
 * GetEResolProbabilityFromDistribution(Ereco,Etrue)
 * GetEResolProbabilityFromFit(Ereco,Etrue)
 *
 * SP: why it is necessary to divide xresol by Ereco before integration: 
 * Remember that the gaussian function used for the resolution has been
 * built as function of x=log(Ereco/Etrue) so as to be normalized:
 * => integral(-inf;+inf) R(x)dx = 1
 * For a fixed value of Etrue, dx = d(log(Ereco)) = dEreco/Ereco
 * So, for a bin in Ereco between E1 and E2 we get:
 * => integral(E1;E2) R(x)*(1/Ereco)*dEreco
 */
double START::ComputeResults::FunctionResolution(double Etrue, double Ereco, Band *band)
{
  double xresol(0);
  
  /*
    if (band->GetUseOfInstrumentEnergyDistribution()) {
    xresol = band->GetEResolProbabilityFromDistribution(Ereco,Etrue);
  }
  else {
    
    xresol = band->GetEResolProbabilityFromFit(Ereco,Etrue);
    
  }
  */
  
  
  xresol = band->GetEResolProbability(Ereco,Etrue);
  

  xresol/=Ereco;
  
  return xresol;
}

/**
 * Return resolution and used for the integral of the resolution
 */

double START::ComputeResults::IntegrantFunctionResolution(double x)
{
  double integrant(0.);
  integrant = FunctionResolution(*nrjreco_tmp,x,myband_tmp);

  return integrant;
}

/**
 * \brief Make and copy in energybin a vector which contains the integral of
 * the resolution on the measured energy for each MC energy
 * \todo Explore if the -1.5/1.5 used is enough
 */

int START::ComputeResults::MakeVectorPartialIntegral(std::vector<Band> &SelectedBands)
{

  for(std::vector<Band>::iterator band = SelectedBands.begin();band!=SelectedBands.end();++band) {

    if(band->GetKeepBand()==0) continue;

    for (std::vector<EnergyBin>::iterator bin = band->ebin.begin();bin!=band->ebin.end();++bin) {

      double ebin_min = bin->GetEmin();
      double ebin_max = bin->GetEmax();
      
      
      double Etrue_min = ebin_min*TMath::Exp(-1.5);
      double Etrue_max = ebin_max*TMath::Exp(1.5);

      //std::cout << "Etrue_min = " << Etrue_min << "Etrue_max = " << Etrue_max << std::endl;
      //std::cout << "ebin_min = " << ebin_min << "ebin_max = " << ebin_max << std::endl;

      /* 
	 Be careful: instrument functions are not defined below FirstEmcVal
	 So, overwrite Etrue_min if necessary:
      */

      double FirstEmcVal = band->GetFirstEmcVal();

      if (Etrue_min<FirstEmcVal) {    
	//std::cout << "Overwriting Etrue_min from " << Etrue_min << " to "<< FirstEmcVal << std::endl;
	Etrue_min=FirstEmcVal;
	
	if (ebin_min/FirstEmcVal < TMath::Exp(0.3)) {
	  //std::cout <<"warning warning warning : new Etrue_min very close to ebin_min !"<< std::endl;
	  //std::cout <<"be careful => this implies than Area below Etrue_min very close to 0"<< std::endl; 
	}
	
      }

      std::vector<double> partialintegral;
      std::vector<double> energy;

      //double n = 30.;
      double n = 100.;

      for(double ien(TMath::Log(Etrue_min)); ien<=TMath::Log(Etrue_max); 
	  ien+=(TMath::Log(Etrue_max)-TMath::Log(Etrue_min))/n) {

	if(TMath::Exp(ien)<band->GetFirstEmcVal()) {
	  partialintegral.push_back(0.);
	  energy.push_back(TMath::Exp(ien));
	  continue;
	}

	double recoenergy = TMath::Exp(ien);
	nrjreco_tmp = &recoenergy;
	ROOT::Math::Functor1D PartialIntegrant(this,&START::ComputeResults::IntegrantFunctionResolution); 
	ROOT::Math::GSLIntegrator ig(fIntegrationType,fabsoluteprecision,frelativeprecision); 
	ig.SetFunction(PartialIntegrant);

	myband_tmp = (&*band);

	DEBUG_OUT_L(2) << "Compute Partial Integral For Band[" << band->GetNbRun() << "] ebin_min = " << ebin_min << " ebin_max = " << ebin_max << std::endl;
	double integral = ig.Integral(ebin_min,ebin_max);
	DEBUG_OUT_L(2) << "\t integral = " << integral << std::endl;
	if(integral<0.) {
	  std::cout << "probleme!!!!" << std::endl;
	  std::cout << "partial integral = " << integral << std::endl;
	  std::cout << "Integrated resolution is positive by definition... ==> EXIT!!" << std::endl;
	  return -1;
	}
	
	partialintegral.push_back(integral);
	energy.push_back(TMath::Exp(ien));

      } // loop on true energy

      std::pair<std::vector<double >, std::vector<double> > pair_en_int;

      pair_en_int = std::make_pair(energy,partialintegral);

      bin->SetPartialIntegral(pair_en_int);

      partialintegral.clear();
      energy.clear();

    } //loop on bin

  } // loop on band

  return 1;

}


/*
 * \brief Compute the expected excess Sth in one EnergyBin
 *
 * \param iband Band
 * \param iem EnergyBin
 * \param paramfit
 * \param etrueminlc true min energy for lightcurve. -1 if not light curve
 * \param etruemaxlc true max energy for lightcurve. -1 if not light curve
 */
double START::ComputeResults::FunctionExpectedExcess(Band *iband, EnergyBin *iem, 
																										 const std::vector<double> &paramfit)
{ 

  double ebin_min = iem->GetEmin();
  double ebin_max = iem->GetEmax();

  double Etrue_min(ebin_min/TMath::Exp(1.5));
  double Etrue_max(ebin_max*TMath::Exp(1.5));

  fHypothesis->SetParameters(paramfit);

  /* 
     Be careful: instrument functions are not defined below FirstEmcVal
     So, overwrite Etrue_min if necessary:
  */

  double FirstEmcVal = iband->GetFirstEmcVal();

  DEBUG_OUT_L(2) << "ebin_min = " << ebin_min << " ebin_max = " << ebin_max << std::endl;
  DEBUG_OUT_L(2) << "Etrue_min = " << Etrue_min << " Etrue_max = " << Etrue_max << std::endl;
  //DEBUG_OUT_L(2) << "etrueminlc = " << etrueminlc << " etruemaxlc = " << etruemaxlc << std::endl;
  DEBUG_OUT_L(2) << "FirstEmcVal = " << FirstEmcVal << std::endl;

  if(DEBUG>1) {
    std::vector<double> partial_bin_energy = iem->GetPartialIntegral().first;
    std::vector<double> partial_bin_resolution = iem->GetPartialIntegral().second;
    std::cout << "energy   p_resolution :" << std::endl;
    for(unsigned int i(0); i<partial_bin_energy.size(); i++)
      std::cout << partial_bin_energy[i] << "   " << partial_bin_resolution[i] << std::endl;
  }

  if (Etrue_min<FirstEmcVal) {    
    //std::cout<<"Overwriting Etrue_min from "<<Etrue_min<<" to "<<FirstEmcVal<<std::endl;
    Etrue_min=FirstEmcVal;

    if (ebin_min/FirstEmcVal < TMath::Exp(0.3)) {
      //std::cout <<"warning warning warning : new Etrue_min very close to ebin_min !"<<std::endl;
      //std::cout <<"be careful => this implies than Area below Etrue_min very close to 0"<<std::endl; 
    }
  }

  /*
  // TODO : JLK remplacer firstemcval par etruemin quand il faut
  // LightCurve case, etrue min and etrue max are given by user
  if(etrueminlc!=-1. && etruemaxlc!=-1.) {
    if(etrueminlc>Etrue_min) { // LC energy min is greater than true min energy 
      Etrue_min = etrueminlc;
      if(Etrue_min<FirstEmcVal) Etrue_min=FirstEmcVal;
      if(etruemaxlc<Etrue_max) Etrue_max=etruemaxlc;
      DEBUG_OUT_L(2) << "*CASE 1 Integration for true energy between " << Etrue_min << " and " << Etrue_max << " Tev" << std::endl;
    }
    else if(etrueminlc<Etrue_min) {
      //Etrue_min = FirstEmcVal;
      if(etruemaxlc<Etrue_max) Etrue_max=etruemaxlc;
      DEBUG_OUT_L(2) << "*CASE 2 Integration for true energy between " << Etrue_min << " and " << Etrue_max << " Tev" << std::endl;
      if(Etrue_min>Etrue_max) return 0.;
    }
    else if(etrueminlc>Etrue_max && etruemaxlc>Etrue_max) {
      DEBUG_OUT_L(2) << "*CASE 3 We return 0!" << std::endl;
      return 0.; 
    }
    else return 0.;
  }
  */

  myband_tmp = iband;
  mynrj_tmp = iem;

  ROOT::Math::Functor1D ExcessIntegrant(this,&START::ComputeResults::IntegrantExpectedExcess); // JLK on pourrait le mettre dans le constructeur
 
  ROOT::Math::GSLIntegrator ig(fIntegrationType,fabsoluteprecision,frelativeprecision); 

  ig.SetFunction(ExcessIntegrant);

  double S = ig.Integral(Etrue_min,Etrue_max);
  double TON = iband->GetLiveTime(); // in hours
  S*=TON*3600.; // CGS units (a peu de choses pres...)

  DEBUG_OUT_L(2) << "Expected excess = " << S << std::endl;

  return S;

}


/*
 * \brief Return the integrated resolution for a given bin and energy
 */
double START::ComputeResults::PartialIntegral(EnergyBin *ebin, double Etrue)
{

  double part(0.);
  
  std::vector<double> energy = ebin->GetPartialIntegral().first;
  std::vector<double> partialintegral = ebin->GetPartialIntegral().second;

  if (energy.size()==0 || partialintegral.size()==0) {
    // VIM : This mean  that the interpolation couldn't have been done previously (i.e problems, or too low on energy)
    return 0.;
  }
  
  int lastbin = ebin->GetPartialIntegral().first.size()-1;
  
  if (Etrue<energy[energy.size()-1]) {

    part = ebin->GetInterpolatedPartialIntegral(Etrue);

    if(part!=part) {
      std::cout << "Integrated resolution is NAN!!!" << std::endl;
      for(unsigned int i(0); i<energy.size(); i++) {
      	std::cout << "energy = " << energy[i] << " | partialintegral = " << partialintegral[i] << std::endl;
      }
      std::cout << "Etrue = " << Etrue << " | part = " << part << std::endl;
      
      exit(EXIT_FAILURE);
    }

    if(part<0.) {
      for(unsigned int i(0); i<energy.size(); i++) {
      	std::cout << "energy = " << energy[i] << " | partialintegral = " << partialintegral[i] << std::endl;
      }
      std::cout << "Integrated resolution < 0 ==> EXIT!!" << std::endl;
      std::cout << "Etrue = " << Etrue << " | part = " << part << std::endl;
      
      exit(EXIT_FAILURE);
    }
    if(part<1.E-20) {      
      return (1.E-20);
    }
  } 
  else {
    part = partialintegral[lastbin] + 
      ((Etrue-energy[lastbin])/(energy[lastbin]-energy[lastbin-1]))*(partialintegral[lastbin]-partialintegral[lastbin-1]);
    if(part<1.E-20) { return (1.E-20);}
  }
  return part;
}

/**
 * \brief Used to compute the integrand of the expected excess
 *
 * \param x Etrue
 */
double START::ComputeResults::IntegrantExpectedExcess(double x)
{

  double integrant(0.);

  integrant = fHypothesis->GetFlux(x);
  integrant*=FunctionEffectiveArea(x,myband_tmp);
  integrant*=PartialIntegral(mynrj_tmp,x);

  return integrant;
}


