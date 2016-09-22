// STL
#include <iostream>
#include <sstream>

// ROOT
#include <TMath.h>

// START
#include <Residuals.hh>

#define DEBUG 0
#include "debugging.hh"

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::Residuals)
#endif

/**
 * \brief Constructor
 */
START::Residuals::Residuals(unsigned int size) {
  if(size==0) {
    fEnergy.clear();
    fEnergyMinus.clear();
    fEnergyPlus.clear();
    fResiduals.clear();
	fResidualsOn.clear();
	fResidualsOff.clear();
    fResidualsSigmaPlus.clear();
    fResidualsSigmaMinus.clear();
    fResiduals3SigmaPlus.clear();
    fResiduals3SigmaMinus.clear();
    fFlux.clear();
    fFluxSigmaPlus.clear();
    fFluxSigmaMinus.clear();
    fFlux3SigmaPlus.clear();
    fFlux3SigmaMinus.clear();
    fIsUpperLimits.clear();
  }
  else {
    fEnergy.resize(size,0);
    fEnergyMinus.resize(size,0);
    fEnergyPlus.resize(size,0);
    fResiduals.resize(size,0);
	fResidualsOn.resize(size,0);
	fResidualsOff.resize(size,0);
    fResidualsSigmaPlus.resize(size,0);
    fResidualsSigmaMinus.resize(size,0);
    fResiduals3SigmaPlus.resize(size,0);
    fResiduals3SigmaMinus.resize(size,0);
    fFlux.resize(size,0);
    fFluxSigmaPlus.resize(size,0);
    fFluxSigmaMinus.resize(size,0);
    fFlux3SigmaPlus.resize(size,0);
    fFlux3SigmaMinus.resize(size,0);
    fIsUpperLimits.resize(size,false);
  }

  fMarkerColor = kBlack;
  fMarkerSize = 1;
  fMarkerStyle = 21;
  fLineColor = kBlack;
  fLineWidth = 1;
}

/**
 * \brief Destructor
 */
START::Residuals::~Residuals() {

}

/**
 * \brief Print
 */
void START::Residuals::Print(Option_t *option) const {
  PrintResiduals();
}

/**
 * \brief Print
 */
void START::Residuals::PrintResiduals(std::ostream &os) const {

  std::ostringstream header;
  header.width(7);
  header << "        ";
  header.width(12);
  header << " e_mean ";
  header.width(14);
  header << " residuals- ";
  header.width(10);
  header << " residuals ";
  header.width(13);
  header << " residuals+ ";
  header.width(8);
  header << " flux- ";
  header.width(12);
  header << " flux ";
  header.width(12);
  header << " flux+ ";
  header.width(12);
  header << " res_on ";
  header.width(12);
  header << " res_off ";
  
  os << header.str() << std::endl;

  for(unsigned int ires(0); ires<fEnergy.size(); ires++) {
    
    std::ostringstream oss1;
    oss1.width(7);
    oss1 << "1 sigma";
    oss1.width(12);
    oss1.precision(4);
    oss1 << fEnergy[ires];
    oss1.width(12);
    oss1 << fResidualsSigmaMinus[ires]-1.;
    oss1.width(12);
    oss1 << fResiduals[ires]-1.;
    oss1.width(12);
    oss1 << fResidualsSigmaPlus[ires]-1.;
    oss1.width(12);
    oss1 << fFluxSigmaMinus[ires];
    oss1.width(12);
    oss1 << fFlux[ires];
    oss1.width(12);
    oss1 << fFluxSigmaPlus[ires];
    if(fIsUpperLimits[ires]) {
      oss1.width(17);
      oss1 << " ==> UpperLimit! ";
    }
    oss1.width(12);
    oss1 << fResidualsOn[ires];
    oss1.width(12);
	oss1 << fResidualsOff[ires];
    os << oss1.str() << std::endl;
    
    std::ostringstream oss2;
    oss2.width(7);
    oss2 << "3 sigma";
    oss2.width(12);
    oss2.precision(4);
    oss2 << fEnergy[ires];
    oss2.width(12);
    oss2 << fResiduals3SigmaMinus[ires]-1.;
    oss2.width(12);
    oss2 << fResiduals[ires]-1.;
    oss2.width(12);
    oss2 << fResiduals3SigmaPlus[ires]-1.;
    oss2.width(12);
    oss2 << fFlux3SigmaMinus[ires];
    oss2.width(12);
    oss2 << fFlux[ires];
    oss2.width(12);
    oss2 << fFlux3SigmaPlus[ires];
	oss1.width(12);
    oss1 << fResidualsOn[ires];
    oss1.width(12);
	oss1 << fResidualsOff[ires];
    os << oss2.str() << std::endl;
  }

  std::cout << "---------------------------------------------------Flux 99.7% (3 sigma)ULs in vFv MeV/cm2/s -------------------------------------------------"<< std::endl; 
  for(unsigned int ires(0); ires<fEnergy.size(); ires++) {
    
    //ADA For SED plotting : Energies in MeV, vFv in MeV/cm2/s 
    std::ostringstream oss3;
    oss3.width(10);
    oss3.precision(3);
    //    oss3 << " ";
    oss3 << std::scientific ;
    oss3 << fEnergy[ires]*1e6;
    oss3 << "    ";
    oss3 << fEnergyMinus[ires]*1e6; //in MeV
    oss3 << "      ";
    oss3 << fEnergyPlus[ires]*1e6;
    oss3 << "    ";
    oss3.width(11);
    oss3.precision(5);
    if(fIsUpperLimits[ires]) {
      oss3 << fFlux3SigmaPlus[ires]*(fEnergy[ires]*1e6*fEnergy[ires]*1e6)*1e-6;// /TeV/cm2:S -> MeV/cm2/s 
      oss3 << "             0    1";
    }
    else {
      oss3 << fFlux[ires]*(fEnergy[ires]*1e6*fEnergy[ires]*1e6)*1e-6;// /TeV/cm2:S -> MeV/cm2/s 
      oss3 << "    ";
      oss3 << (fFluxSigmaPlus[ires]-fFlux[ires])*(fEnergy[ires]*1e6*fEnergy[ires]*1e6)*1e-6;// /TeV/cm2:S -> MeV/cm2/s 
      oss3 << "    0";
    }
    os << oss3.str() << std::endl;
    
    //1.5244e+05     2.1101e+05   1.7543e+05   4.132594e-06             0    1

  }
  std::cout << "---------------------------------------------------Flux 68.3% (1 sigma)ULs in vFv MeV/cm2/s -------------------------------------------------"<< std::endl; 
  for(unsigned int ires(0); ires<fEnergy.size(); ires++) {
    
    //ADA For SED plotting : Energies in MeV, vFv in MeV/cm2/s 
    std::ostringstream oss3;
    oss3.width(10);
    oss3.precision(3);
    //    oss3 << " ";
    oss3 << std::scientific ;
    oss3 << fEnergy[ires]*1e6;
    oss3 << "    ";
    oss3 << fEnergyMinus[ires]*1e6; //in MeV
    oss3 << "      ";
    oss3 << fEnergyPlus[ires]*1e6;
    oss3 << "    ";
    oss3.width(11);
    oss3.precision(5);
    oss3 << fFluxSigmaPlus[ires]*(fEnergy[ires]*1e6*fEnergy[ires]*1e6)*1e-6;// /TeV/cm2:S -> MeV/cm2/s 
    oss3 << "             0    1";
    os << oss3.str() << std::endl;


  }
  
}
