//STL
#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>

// ROOT
#include "Minuit2/FCNBase.h"
#include "Math/InterpolationTypes.h"
#include "Math/Interpolator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/WrappedMultiTF1.h"
#include "TF2.h"
#include "TMath.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "TFile.h"
#include "TStopwatch.h"

// START
#include "FCNLikelihood.hh"
#include "STARTUtils.hh"

// Utilities
#define DEBUG 0
#include "debugging.hh"

#define INFO std::cout << INFOCOLOR << "FCNLikelihood> " << RESETCOLOR
#define WARNING std::cout << WARNINGCOLOR << "FCNLikelihood> " << RESETCOLOR

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::FCNLikelihood)
#endif



/**
 * \brief Constructor
 */
START::FCNLikelihood::FCNLikelihood(const std::vector<Band> &SelectedBands, Hypothesis &hypo, Config &Configuration, bool verbose)
:fErrorDef(0.5),
	fverbose(verbose),
	fFirstCallDone(false),
	fHypothesis(0),
	fCompRes(0),
	fBandsGroupingType(START::STARTUtils::NoGrouping)
{
	fConfig = &Configuration;

	fBandArray = SelectedBands;

	fHypothesis = &hypo;

	fCompRes = new ComputeResults(SelectedBands,*fHypothesis);
}


/**
 * \brief Destructor
 */
START::FCNLikelihood::~FCNLikelihood()
{
	if (fCompRes!=0) delete fCompRes;
	fCompRes = 0;

}

/**
 * \brief Likelihood function, takes parameters to fit in argument
 */
double START::FCNLikelihood::operator()(const std::vector<double> &par) const
{
	TStopwatch wtimeoneinter;
	wtimeoneinter.Start();

	/*
	  log L = sum Band(B) sum bin(b) { nON_Bb * log(S_Bb + alpha_B * p_Bb)
	  + nOFF_Bb * log(p_Bb)
	  - (alpha_B +1) * p_Bb - S_Bb }

	  with : S = expected number of events
	  p_Bb = 1/(alpha_B*(alpha_B+1)) * (a_Bb
	  + sqrt(a_Bb**2 + 4 * alpha_B*(alpha_B+1) * nOFF_Bb * S_Bb))
	  and    a_Bb = alpha_B * (nON_Bb + nOFF_Bb) - (alpha_B - 1) * S_Bb


	*/

	double alpha(0.), alphaplus(0.), alphaless(0.);
	double nON(0.), nOFF(0.), a(0.), p(0.);

	double calc1(0.), calc2(0.), calc3(0.); // first three terms in L

	double fcnvalue(0); // result of the likelihood
	double fcnvaluetmp(0); // tmp result used to check if fcnvalue is not a number

	const double normalization = STARTUtils::GetNormalizationConstante(); // normalization for fixed parameters

	bool isunphysical = false;
	bool isthereNaN = false;

	for(unsigned int i(0); i<par.size(); i++) {
		if(TMath::IsNaN(fHypothesis->GetParameters()[i])) 
			isthereNaN=true;
	}

	// JLK Ajout le switch ici la
	const std::vector<START::Band>& vBandsToUse = (fBandsGroupingType==START::STARTUtils::NoGrouping) ? fBandArray : 
		(fBandsGroupingType==START::STARTUtils::StackInOneBand) ? GetStackBandFromAllBands(par, 0., 90.) : 
		(fBandsGroupingType==START::STARTUtils::StackInZenBands) ?  GetStackZenBandsFromAllBands(par) : fBandArray;

	for(std::vector<Band>::const_iterator band = vBandsToUse.begin(); band!=vBandsToUse.end(); ++band) { // VIM
    
		// We are interested only in the selected bands (and may the force be with you :) )
		if(0==band->GetKeepBand()) continue; 

		alpha = band->GetAlphaRun(); // alpha run
		alphaplus = alpha+1.;
		alphaless = alpha-1.;

		for (std::vector<EnergyBin>::const_iterator bin = band->ebin.begin();bin!=band->ebin.end();++bin) {

			if (bin->GetKeepBin() == 0) continue; // filtre pour bins en energie en dessous du seuil
			double S(0);
			if(fBandsGroupingType==START::STARTUtils::StackInOneBand || fBandsGroupingType==START::STARTUtils::StackInZenBands )
				S+=bin->GetSth();
			else{
				S = fCompRes->FunctionExpectedExcess(const_cast<Band*>(&*band),const_cast<EnergyBin*>(&*bin),par); 
				if(fHypothesis->GetVectorNormalizedParameters().size()>0) S*=normalization;
			}
			nON = bin->GetOn();
			nOFF = bin->GetOff();

			// JLK : S can't be < 0 by definition, see classes MinimizeFactory and Hypothesis
			// so if it's happening, there is a (nice) bug or you're doing contours/scans and
			// in these cases, it's allowed

			if (S<0 && fverbose)
				isunphysical = true;

			p = STARTUtils::GetExpectedOff(nON,nOFF,alpha,S);

			if (nON==0)
				calc1=0;
			else {
				calc1 = S + alpha*p;
				if (calc1>0.) {
					calc1=TMath::Log(calc1);
					calc1*= nON;
				}
				else {
					calc1=0.;
				}
			}

			if (p>0) {
				calc2 = TMath::Log(p);
				calc2*= nOFF;
			} 
			else calc2=0;
      
			calc3 = p;
			calc3*= alphaplus;
      
			fcnvaluetmp=0.;
			fcnvaluetmp+= calc1;
			fcnvaluetmp+= calc2;
			fcnvaluetmp-= calc3;
			fcnvaluetmp-= S;
      
			fcnvalue+=fcnvaluetmp;

		} // loop on bin
    
	} // loop on band
  
	fcnvalue*= -1.; // MINUIT2 does minimization

	wtimeoneinter.Stop();

	if(fverbose) {
		INFO << "fcnvalue = " << -fcnvalue << ", time to do it : " << wtimeoneinter.RealTime()
			 << "s" << std::endl;
		if(isunphysical || isthereNaN) WARNING << "Unphysical region touched : " << std::endl;
		for(unsigned int i(0); i<par.size(); i++) {
			std::cout << fHypothesis->GetParametersNames()[i] << " = " << fHypothesis->GetParameters()[i] << std::endl;
		}
	}
  
	return fcnvalue;
}


/**
 * \brief VIM passera faire le menage dans son propre code
 */
const std::vector<START::Band> START::FCNLikelihood::GetStackBandFromAllBands(const std::vector<double> &par, float zmin, float zmax) const {

	//ADA 
	int NBandsAdded =0;
	
	// Do an iteration where I fill the expected value for one band, and cumulate the ON and OFF !
	START::Band SumBand;
	SumBand.SetKeepBand(1);

	//while (true) {
	for(std::vector<Band>::const_iterator band = fBandArray.begin();band!=fBandArray.end();++band) {
		if (band->GetKeepBand()==0) { continue; }

		SumBand.ebin = band->ebin;
		break;
	}
  
	//SumBand.ebin = fBandArray.begin()->ebin;
	std::vector<Double_t> LiveTimePerBin;
	for (std::vector<START::EnergyBin>::iterator it = SumBand.ebin.begin(); it!=SumBand.ebin.end(); ++it) {
		it->ClearEnergyBinInfo();
		LiveTimePerBin.push_back(0.);
	}
  
	// Fill the SumBand with the sum of the contribution
	Double_t TotalOff = 0.;
	Double_t TotalAlphaOff = 0.;
	Double_t TotalLiveTimeAlpha = 0.;
	Double_t TotalLiveTime = 0.;

	const double normalization = STARTUtils::GetNormalizationConstante(); // normalization for fixed parameters

	for(std::vector<Band>::const_iterator band = fBandArray.begin();band!=fBandArray.end();++band) {
   
		if(0==band->GetKeepBand()) continue; // We are interested only in the selected bands (and may the force be with you :) )
		// ADA add zenith angle condition : can be extended to add offset, efficiency conditions as well	
		double Zen = band->GetZenON();
		if(!fFirstCallDone){
			std::cout << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< In GetStackBandFromAllBands : Zen of band ()= " << Zen << " Zenmin = " << zmin << " Zenmax = " << zmax << std::endl;  
		}
		if((Zen < zmin) || (Zen> zmax))
			continue;
		
		NBandsAdded++;
		TotalOff += band->GetNOffTot(false);
		TotalAlphaOff += band->GetNOffTot(false) * band->GetAlphaRun();
		TotalLiveTimeAlpha += band->GetAlphaRun() * band->GetLiveTime();
		TotalLiveTime += band->GetLiveTime();
    
		std::vector<Double_t>::iterator it_LT = LiveTimePerBin.begin();
		std::vector<EnergyBin>::iterator sum_bin = SumBand.ebin.begin();
		std::vector<EnergyBin>::iterator sum_bin_end = SumBand.ebin.end();
    
		for (std::vector<EnergyBin>::const_iterator bin = band->ebin.begin(), bin_end = band->ebin.end(); 
			 bin!=bin_end && sum_bin!=sum_bin_end;
			 ++bin, ++sum_bin, ++it_LT) {
      
			if (bin->GetKeepBin() != 0) { // If At least there is one valid energy bin in the band vector, then we keep the bin !
				sum_bin->SetKeepBin(1);
			}
      
			if (bin->GetKeepBin() == 0) continue; // filtre pour bins en énergie en dessous du seuil
			DEBUG_OUT << "gonna integrate between [" << bin->GetEmin() << "," << bin->GetEmax() << "]" << std::endl;
			double S = fCompRes->FunctionExpectedExcess(const_cast<Band*>(&*band),const_cast<EnergyBin*>(&*bin),par); 
			if(fHypothesis->GetVectorNormalizedParameters().size()>0) S*=normalization;

			sum_bin->SetSth( sum_bin->GetSth() + S);
			sum_bin->SetOn( sum_bin->GetOn() + bin->GetOn() );
			sum_bin->SetOff( sum_bin->GetOff() + bin->GetOff() );
			// For the moment, Alpha = Alpha*NOff (in order to compute the average )
			sum_bin->SetAlpha( sum_bin->GetAlpha()  + bin->GetAlpha() * bin->GetOff() );
			// For the moment, livetime = Alpha *Livetime, for the average alpha
			sum_bin->SetLiveTime( sum_bin->GetLiveTime() + bin->GetLiveTime() * bin->GetAlpha() ); 
			*it_LT = bin->GetLiveTime();
		}
	}

	if (TotalOff>0.01) {
		SumBand.SetAlphaRun( TotalAlphaOff/TotalOff );
	}
	else {
		SumBand.SetAlphaRun( TotalLiveTimeAlpha/ TotalLiveTime );
	}
	SumBand.SetLiveTime( TotalLiveTime );
  
	// Now compute the average alpha per bin :
	std::vector<Double_t>::iterator it_LT = LiveTimePerBin.begin();
	for (std::vector<EnergyBin>::iterator sum_bin = SumBand.ebin.begin(), sum_bin_end = SumBand.ebin.end();
		 sum_bin != sum_bin_end; ++sum_bin, ++it_LT) {
		if (sum_bin->GetKeepBin()==0) { continue; }
    
		if (sum_bin->GetOff()>0.01) {
			sum_bin->SetAlpha( sum_bin->GetAlpha() / sum_bin->GetOff() );
		}
		else {
			sum_bin->SetAlpha( sum_bin->GetLiveTime() / *it_LT ); // Backup Strat !
		}
		sum_bin->SetLiveTime( *it_LT ); // Put a correct time in the sum_bin !
	}
  
	LiveTimePerBin.clear(); // No Need Anymore

	std::vector<Band> BandArray;
	if(NBandsAdded>0){
		BandArray.push_back(SumBand);
	}
	return BandArray;
}

/**
 * \brief VIM passera faire le menage dans le code de ADA !
 */
const std::vector<START::Band> START::FCNLikelihood::GetStackZenBandsFromAllBands(const std::vector<double> &par) const {
  
	// Do an iteration where I fill the expected value for N bands, according to zenith angle intervals, and cumulate the ON and OFF !

	std::vector<Band> BandArray;
  
	//ADA Hard coded For Vela :  To be displaced into the config file 
	// Vela Zen=[22.0, 40.] deg, CosZen step = 0.02 and 10 bins  , cos(20.9)=0.934, cos(20.9)+10*0.02=0.734 = cos(42.8) 
	// double CosZenmin=0.934;
	// double CosZenmax=0.734;
	// int nbinZen=4;
	// double Zenbinsize=(CosZenmin-CosZenmax)/nbinZen;

	double pi=3.14159265359;
	double CosZenmin=cos(fConfig->GetBandZenMax()*pi/180);
	double CosZenmax=cos(fConfig->GetBandZenMin()*pi/180);
	int nbinZen=fConfig->GetBandZenbin_UseReproj();
	double Zenbinsize=(CosZenmin-CosZenmax)/nbinZen;



	for(int iz=0;iz<nbinZen;iz++){
		double czmin=CosZenmin - iz*Zenbinsize;
		double czmax=CosZenmin - (iz+1)*Zenbinsize;

		double zmin=180.*TMath::ACos(czmin)/TMath::Pi();
		double zmax=180.*TMath::ACos(czmax)/TMath::Pi();

		const std::vector<START::Band>& vBandinz = GetStackBandFromAllBands(par,zmin,zmax);
		if(vBandinz.size())
			BandArray.push_back(vBandinz.front()); // There is onlony one band 
	}
	//Flag for debugging
	if(!fFirstCallDone){
		std::cout << " -------------------------------------------------- First call for  GetStackZenBandsFromAllBands : ----------------------------------" << std::endl;
		for(std::vector<Band>::const_iterator band = BandArray.begin();band!=BandArray.end();++band) 
			band->Print();
		fFirstCallDone=true;
	}

	return BandArray;

}
