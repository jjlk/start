// STL
#include <iostream>

// ROOT
#include <TMath.h>

// START
#include "MonteCarlo.hh"

// Utilities
#define DEBUG 0
#include "debugging.hh"

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::MonteCarlo)
#endif

/**
 * \brief Constructor
 */
START::MonteCarlo::MonteCarlo(START::STARTUtils::IrfOpt Type)
{
	fIrfOpt = Type;
	InitMCParameters();
}
/**
 * \brief Destructor
 */
START::MonteCarlo::~MonteCarlo()
{

}

/**
 * \brief Build the MonteCarlo parameters used in simulations
 */
void START::MonteCarlo::InitMCParameters()
{
	switch(fIrfOpt) {
	case STARTUtils::HapFr_Hess_Stereo: {
		double eff[] = {50,60,70,80,90,100};
		double zen[] = {0,18,26,32,37,41,46,50,53,57,60,63,67,70};
		double off[] = {0,0.5,1,1.5,2.,2.5};
		double en[] = {0.02,0.03,0.05,0.08,0.125,0.2,0.3,0.5,0.8,1.25,
					   2.,3.,5.,8.,12.5,20.,30.,50.,80.,125.};
		TString teltype[] = {"","3Tel"};
		fEfficiencyMC = std::vector<double>(eff,
											eff+sizeof(eff) / sizeof(double));
		fZenithMC = std::vector<double>(zen,
										zen+sizeof(zen) / sizeof(double));
		fOffsetMC = std::vector<double>(off,
										off+sizeof(off) / sizeof(double));
		fEnergyMC = std::vector<double>(en,
										en+sizeof(en) / sizeof(double));
		fTelTypeMC = std::vector<TString>(teltype,
										  teltype+sizeof(teltype) / sizeof(TString));
		break;
	}
	case STARTUtils::HapFr_Hess_Hybrid: {
		double eff[] = {201401};
		//		double zen[] = {0,18,26,32,37,41,46,50,53,57,60,63,67,70};
		//		double off[] = {0,0.5,1,1.5,2.,2.5};
		//ADA : SimTelArray phase2a80 
		double zen[] = {0,20,30,40,45,50};
		double off[] = {0.5};
		double en[] = {0.02,0.03,0.05,0.08,0.125,0.2,0.3,0.5,0.8,1.25,
					   2.,3.,5.,8.,12.5,20.,30.,50.,80.,125.};
		TString teltype[] = {""};
		fEfficiencyMC = std::vector<double>(eff,
											eff+sizeof(eff) / sizeof(double));
		fZenithMC = std::vector<double>(zen,
										zen+sizeof(zen) / sizeof(double));
		fOffsetMC = std::vector<double>(off,
										off+sizeof(off) / sizeof(double));
		fEnergyMC = std::vector<double>(en,
										en+sizeof(en) / sizeof(double));
		fTelTypeMC = std::vector<TString>(teltype,
										  teltype+sizeof(teltype) / sizeof(TString));
		break;
	}
	case STARTUtils::HapFr_Hess_Mono: {
		double eff[] = {201401};
		double zen[] = {0,10,20,30,40,45,50};
		//double zen[] = {0,20,30,40,45,50,55};
		double off[] = {0.5};
		double en[] = {0.005, 0.008,0.0125, 0.02,0.03,0.05,0.08,0.125,0.2,0.3,0.5,0.8,1.25,
			       2.,3.,5.,8.,12.5,20.};
		//		double en[] = {0.02,0.03,0.05,0.08,0.125,0.2,0.3,0.5,0.8,1.25,
		//			   2.,3.,5.,8.,12.5,20.,30.,50.,80.,125.};
		TString teltype[] = {""};
		fEfficiencyMC = std::vector<double>(eff,
											eff+sizeof(eff) / sizeof(double));
		fZenithMC = std::vector<double>(zen,
										zen+sizeof(zen) / sizeof(double));
		fOffsetMC = std::vector<double>(off,
										off+sizeof(off) / sizeof(double));
		fEnergyMC = std::vector<double>(en,
										en+sizeof(en) / sizeof(double));
		fTelTypeMC = std::vector<TString>(teltype,
										  teltype+sizeof(teltype) / sizeof(TString));
		break;
	}
	default:
		exit(EXIT_FAILURE);
	}
	
}

/**
 * \brief Transform telcode in string use by the code
 */
TString START::MonteCarlo::GetStringFromTelcode(unsigned int telcode) const
{
	if(telcode<30) return "3Tel";
	//else if(telcode==30) return "";
	else if(telcode>=30) return ""; // JLK HACK for HESS-II
  
	return "Unknown";
}

/**
 * \brief Check if parameter is MC
 * \param zen zenith
 */
bool START::MonteCarlo::IsItMCZenith(double zen) const {

	for(std::vector<double>::const_iterator it=fZenithMC.begin(); it!=fZenithMC.end(); ++it) {
		if(zen==*it) return true;
	}

	return false;

}

/**
 * \brief Check if parameter is MC
 * \param eff efficiency
 */
bool START::MonteCarlo::IsItMCEfficiency(double eff) const {

	for(std::vector<double>::const_iterator it=fEfficiencyMC.begin(); it!=fEfficiencyMC.end(); ++it) {
		if(eff==*it) return true;
	}

	return false;

}

/**
 * \brief Check if parameter is MC
 * \param off offset
 */
bool START::MonteCarlo::IsItMCOffset(double off) const {

	for(std::vector<double>::const_iterator it=fOffsetMC.begin(); it!=fOffsetMC.end(); ++it) {
		if(off==*it) return true;
	}

	return false;

}

/**
 * \brief Check if parameter is MC
 * \param telcode telcode
 */
bool START::MonteCarlo::IsItMCTelCode(unsigned telcode) const {

	TString code = GetStringFromTelcode(telcode);

	for(std::vector<TString>::const_iterator it=fTelTypeMC.begin(); it!=fTelTypeMC.end(); ++it) {
		if(code==*it) return true;
	}

	return false;

}

/**
 * \brief Check if parameter is MC
 * \param energy energy
 */
bool START::MonteCarlo::IsItMCEnergy(double energy) const {

	for(std::vector<double>::const_iterator it=fEnergyMC.begin(); it!=fEnergyMC.end(); ++it) {
		if(energy==*it) return true;
	}

	return false;

}

void START::MonteCarlo::PrintMCParameters() const {

	INFO_OUT << "Efficiency : " << std::endl;
	for(std::vector<double>::const_iterator it=fEfficiencyMC.begin(); it!=fEfficiencyMC.end(); ++it) std::cout << *it << " ";
	std::cout << std::endl;
	INFO_OUT << "Zenith angle : " << std::endl;
	for(std::vector<double>::const_iterator it=fZenithMC.begin(); it!=fZenithMC.end(); ++it) std::cout << *it << " ";
	std::cout << std::endl;
	INFO_OUT << "Offset : " << std::endl;
	for(std::vector<double>::const_iterator it=fOffsetMC.begin(); it!=fOffsetMC.end(); ++it) std::cout << *it << " ";
	std::cout << std::endl;
	INFO_OUT << "Energy : " << std::endl;
	for(std::vector<double>::const_iterator it=fEnergyMC.begin(); it!=fEnergyMC.end(); ++it) std::cout << *it << " ";
	std::cout << std::endl;
	INFO_OUT << "TelType : " << std::endl;
	for(std::vector<TString>::const_iterator it=fTelTypeMC.begin(); it!=fTelTypeMC.end(); ++it) std::cout << *it << " ";
	std::cout << std::endl;

}


std::pair<double,double> START::MonteCarlo::GetMCLimitingEfficiency(double eff) const {

	if(fEfficiencyMC.size()==1)
		return std::make_pair(fEfficiencyMC[0],fEfficiencyMC[0]);

	for(unsigned int i(0); i<fEfficiencyMC.size()-1; i++) {
		if(eff >= fEfficiencyMC[i] && eff <= fEfficiencyMC[i+1]) return std::make_pair(fEfficiencyMC[i],fEfficiencyMC[i+1]); 
	}

	return std::make_pair(-1.,-1.);

}

std::pair<double,double> START::MonteCarlo::GetMCLimitingZenith(double zen) const {

	if(fZenithMC.size()==1)
		return std::make_pair(fZenithMC[0],fZenithMC[0]);

	for(unsigned int i(0); i<fZenithMC.size()-1; i++) {
		if(zen >= fZenithMC[i] && zen <= fZenithMC[i+1]) return std::make_pair(fZenithMC[i],fZenithMC[i+1]); 
	}

	return std::make_pair(-1.,-1.);

}

std::pair<double,double> START::MonteCarlo::GetMCLimitingOffset(double off) const {

	if(fOffsetMC.size()==1)
		return std::make_pair(fOffsetMC[0],fOffsetMC[0]);

	for(unsigned int i(0); i<fOffsetMC.size()-1; i++) {
		if(off >= fOffsetMC[i] && off <= fOffsetMC[i+1]) return std::make_pair(fOffsetMC[i],fOffsetMC[i+1]); 
	}

	return std::make_pair(-1.,-1.);

}

std::pair<double,double> START::MonteCarlo::GetMCLimitingEnergy(double en) const {

	if(fEnergyMC.size()==1)
		return std::make_pair(fEnergyMC[0],fEnergyMC[0]);
	
	for(unsigned int i(0); i<fEnergyMC.size()-1; i++) {
		if(en >= fEnergyMC[i] && en <= fEnergyMC[i+1]) return std::make_pair(fEnergyMC[i],fEnergyMC[i+1]); 
	}

	return std::make_pair(-1.,-1.);

}



