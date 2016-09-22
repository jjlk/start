//STL
#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>
#include <limits>
#include <iomanip>
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

#define FITMETHOD 3
// 0 Normal, 1 Cumulate All Run with same alpha, 2 Regis-VIM Method. Everything different from 0 is experimental !!!

#define INFO std::cout << INFOCOLOR << "FCNLikelihood> " << RESETCOLOR
#define WARNING std::cout << WARNINGCOLOR << "FCNLikelihood> " << RESETCOLOR

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::FCNLikelihood)
#endif



/**
 * \brief Constructor
 */
START::FCNLikelihood::FCNLikelihood(const std::vector<Band> &SelectedBands, Hypothesis &hypo,bool verbose)
:fErrorDef(0.5),
																  fverbose(verbose),
																  fHypothesis(0),
																  fCompRes(0)
{
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

#if FITMETHOD == 0

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
		if(TMath::IsNaN(fHypothesis->GetParameters()[i])) isthereNaN=true;
	}

	for(std::vector<Band>::const_iterator band = fBandArray.begin();band!=fBandArray.end();++band) { // VIM
    
		if(0==band->GetKeepBand()) continue; // We are interested only in the selected bands (and may the force be with you :) )

		alpha = band->GetAlphaRun(); // alpha run
		alphaplus = alpha+1.;
		alphaless = alpha-1.;

		for (std::vector<EnergyBin>::const_iterator bin = band->ebin.begin();bin!=band->ebin.end();++bin) {

			if (bin->GetKeepBin() == 0) continue; // filtre pour bins en énergie en dessous du seuil
 
			double S = fCompRes->FunctionExpectedExcess(const_cast<Band*>(&*band),const_cast<EnergyBin*>(&*bin),par); 
			if(fHypothesis->GetVectorNormalizedParameters().size()>0) S*=normalization;

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
			if (nON>0) {
				fcnvaluetmp-=(nON*TMath::Log(nON)-nON);
			}
			if (nOFF>0) {
				fcnvaluetmp-=(nOFF*TMath::Log(nOFF)-nOFF);
			}
			fcnvalue+=fcnvaluetmp;

		} // loop on bin
    
	} // loop on band
  
	fcnvalue*= -1.; // MINUIT2 does minimization

	wtimeoneinter.Stop();

	if(fverbose) {
		std::streamsize ssp = std::cout.precision();
		INFO << "fcnvalue = " << std::setprecision(std::numeric_limits<double>::digits10)
			 << fcnvalue << std::setprecision(ssp) << ", time to do it : " << wtimeoneinter.RealTime()
			 << "s" << std::endl;
		if(isunphysical || isthereNaN) WARNING << "Unphysical region touched : " << std::endl;
		for(unsigned int i(0); i<par.size(); i++) {
			std::cout << fHypothesis->GetParametersNames()[i] << " = " << fHypothesis->GetParameters()[i] << std::endl;
		}
	}
  
	return fcnvalue;
}

#elif FITMETHOD == 1
/**
 * \brief Cumulative Fit Method, i.e, neglect the atmospheric variation from run to run, and bruteforcely cumulate in "one" band
 * Assume that the EnergyBin have the same definition for all START::Band
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
		if(TMath::IsNaN(fHypothesis->GetParameters()[i])) isthereNaN=true;
	}

	// Principle : Regroup the band with the same alpha, because it seems valid to me from a mathematical point of view
	// The group is only done on the number of events (ON and OFF). The expected excess is still computed bandwise, and then sum !
	//
	// Note : I don't know how to regroup correctly multiple band when alpha is different. 
	// I think there is a way to do it, because what is fundamentally different between two observation with two different alpha
	// or two observation with the same alpha ??????? ==> but I am a quiche in Math :(((
  
	std::map<Double_t, START::Band > MapAlphaBand;
	std::map<Double_t, Int_t > MapAlphaRun;
  
	// Loop on the band in order to regroup the bands:
	std::vector<START::Band>::const_iterator it_band_end = fBandArray.end();
	for (std::vector<Band>::const_iterator it_band = fBandArray.begin(); it_band!=it_band_end;++it_band) {
		if (0==it_band->GetKeepBand()) {
			continue; // Rejected is the band !
		}
   
		Double_t NReg = Double_t(TMath::Nint(1./it_band->GetAlphaRun())); 
		// VIM : Be carefull, not good for other method than reflected (quick and dirty test)
		// Originally, I did it with alpha, but sometimes two "same" alpha yield two entries in the map
		// Should rewrite the map condition, to declare equivalent two very close alpha.
		// For that, see : http://stackoverflow.com/questions/6684573/floating-point-keys-in-stdmap

		// Find the band corresponding to the number of region in the run,
		std::map<Double_t,START::Band>::iterator it_mab = MapAlphaBand.find(NReg);
		if (it_mab==MapAlphaBand.end()) {
			START::Band b; // Create a simple band
			b.SetNbRun(0);
			b.SetAlphaRun(1./NReg);
			for (std::vector<EnergyBin>::const_iterator it_ebin = it_band->ebin.begin();it_ebin!=it_band->ebin.end();++it_ebin) {
				START::EnergyBin EnerBin(it_ebin->GetEmin(),it_ebin->GetEmax());
				EnerBin.SetKeepBin(0);
				EnerBin.SetOn(0.);
				EnerBin.SetOff(0.);
				EnerBin.SetSth(0.);
				EnerBin.SetOffFitted(0.);
				EnerBin.SetAlpha(1./NReg);
				b.ebin.push_back(EnerBin);
			}
      
			std::pair< std::map<Double_t,START::Band>::iterator, Bool_t > ret = MapAlphaBand.insert( std::pair<Double_t,START::Band>(NReg,b) );
			if (!ret.second) {
				std::cout << "Can't Insert a Band, that's weird, because if we enter the if, it is precisly because it does not exist.... !" << std::endl;
				continue;
			}
			else {
				it_mab = ret.first;	
			}
		}
    
		// Now I should have a valid band to point at, we can now loop and add the information !
   
		// Add the Band Wise information : 
		it_mab->second.SetNbRun( it_mab->second.GetNbRun() + 1 );
    
		// Add the Energy Wise information :
		// This for loops on both the EnergyBin of the Current Band (it_ebin) and the EnergyBin of the SumBand (it_sumband_ebin)
		std::vector<EnergyBin>::iterator it_sumband_ebin = it_mab->second.ebin.begin();
		for (std::vector<EnergyBin>::const_iterator it_ebin = it_band->ebin.begin();
			 it_ebin!=it_band->ebin.end() && it_sumband_ebin!=it_mab->second.ebin.end();
			 ++it_ebin,++it_sumband_ebin) {
			//std::cout << "Current EBin of SumBand = [" << it_sumband_ebin->GetEmin() << "," << it_sumband_ebin->GetEmax() << "], EBin of CurBand = [" << it_ebin->GetEmin() << "," << it_ebin->GetEmax() << "]" << std::endl;
			if (it_ebin->GetKeepBin()==0) {
				continue;
			}

			it_sumband_ebin->AddInfoFromEBin( *it_ebin );

			// Now compute the expected excess per bin : 
			double S = fCompRes->FunctionExpectedExcess(const_cast<Band*>(&*it_band),const_cast<EnergyBin*>(&*it_ebin),par); 
			if(fHypothesis->GetVectorNormalizedParameters().size()>0) S*=normalization;
			double Sprevious = it_sumband_ebin->GetSth();
			it_sumband_ebin->SetSth( S + Sprevious );
			// Ok, that might be dangerous, because that assume that it_ebin->GetSth(0) is 0 !
			// but I would prefer to avoid this line (const_cast<EnergyBin*>(&*bin))->SetSth(S);
			// because it seems harmless regarding the code, but my level of C++ is not that good !
		}
	}

  

      
	//   if (it_mab==MapAlphaBand.end())  {
	//     // Entry does not exist ==> Create One !
	//     START::Band b( *it_band );
	//     b.SetNbRun(0);
	//     //      MapAlphaBand[it_band->GetAlphaRun()] = START::Band( *it_band ); // Could be monstruously speed up by not doing a copy, but create a basic band but it is for testing purpose, so it's faster to do so !!!
	//     MapAlphaBand[NReg] = b; // Could be monstruously speed up by not doing a copy, but create a basic band but it is for testing purpose, so it's faster to do so !!!
	//     MapAlphaRun[NReg] = 1;
	//   }
	//   else {
	//     // Already exist !
	//     std::vector<START::Band> vec_b;
	//     vec_b.push_back( *it_band );
	//     it_mab->second.AddInfoFromBands(vec_b,true);
	//     ++(MapAlphaRun[NReg]);
	//   }
	// }
    
	/*
	  std::cout << "MapAlphaBand.size() = " << MapAlphaBand.size();
	  std::cout << " Detail : ";
	  Int_t NRunTotal = 0;
	  for (std::map<Double_t,START::Band>::iterator it_mab = MapAlphaBand.begin(); it_mab!=MapAlphaBand.end(); ++it_mab) {
	  std::cout << "[" << it_mab->first << "," << 1./it_mab->first << "," << it_mab->second.GetNbRun() << "],";
	  //std::cout << "[" << it_mab->first << "," << TMath::Nint(1./it_mab->first) << "," << MapAlphaRun[it_mab->first] << "],";
	  NRunTotal+= it_mab->second.GetNbRun();
    
	  }
	  std::cout << " TotalBand = " << NRunTotal;
	  std::cout << std::endl;
	*/
    

	// Now let's compute the Likelihood and the expected OFF !
	fcnvalue=0;
	fcnvaluetmp=0;
  
	for(std::map<Double_t,START::Band>::iterator it_mab = MapAlphaBand.begin(); it_mab!=MapAlphaBand.end();++it_mab) {
		for (std::vector<EnergyBin>::iterator bin = it_mab->second.ebin.begin();bin!=it_mab->second.ebin.end();++bin) {
			// if (bin->GetKeepBin() == 0) continue; // No Need for KeepBin, because it is only filled by above thresold information
			nON = bin->GetOn();
			nOFF = bin->GetOff();
			alpha = 1./it_mab->first;
			alphaplus = alpha+1.;
			alphaless = alpha-1.;
			double S = bin->GetSth();
			p = STARTUtils::GetExpectedOff(nON,nOFF,alpha,S);
			bin->SetOffFitted(p);

			if (S<0 && fverbose)
				isunphysical = true;

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
      
			Double_t measured_non = nON;
			Double_t measured_noff = nOFF;
			Double_t LogFactMesNON = ( (measured_non==0.) ? 0. : (measured_non*TMath::Log(measured_non)-measured_non) );
			Double_t LogFactMesNOFF = ( (measured_noff==0.) ? 0. : (measured_noff*TMath::Log(measured_noff)-measured_noff) );
			fcnvaluetmp = fcnvaluetmp - (LogFactMesNON + LogFactMesNOFF); // TEST !!!
			fcnvalue+=fcnvaluetmp;
		}
	}


	/*
	  static int SAMERENSLIP = 0;
	  if (SAMERENSLIP==60) {
	  for (std::map<Double_t,START::Band>::iterator it_mab = MapAlphaBand.begin(); it_mab!=MapAlphaBand.end(); ++it_mab) {
	  it_mab->second.Print();
	  }
	  }
	  ++SAMERENSLIP;
	*/
  
  
	fcnvalue*= -1.; // MINUIT2 does minimization

	wtimeoneinter.Stop();

	if(fverbose) {
		std::streamsize ssp = std::cout.precision();
		INFO << "fcnvalue = " << std::setprecision(std::numeric_limits<double>::digits10) << fcnvalue << std::setprecision(ssp) << ", time to do it : " << wtimeoneinter.RealTime()
			 << "s" << std::endl;
		// INFO << "fcnvalue = " << fcnvalue << ", time to do it : " << wtimeoneinter.RealTime()
		// 	 << "s" << std::endl;
		if(isunphysical || isthereNaN) WARNING << "Unphysical region touched : " << std::endl;
		for(unsigned int i(0); i<par.size(); i++) {
			std::cout << fHypothesis->GetParametersNames()[i] << " = " << fHypothesis->GetParameters()[i] << std::endl;
		}
	}
  
	return fcnvalue;
}


#elif FITMETHOD == 2
/**
 * \brief Regis-VIM Method, i.e assume the background behave like a PWL above a certain energy
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
		if(TMath::IsNaN(fHypothesis->GetParameters()[i])) isthereNaN=true;
	}
  
  
	for(std::vector<Band>::const_iterator band = fBandArray.begin();band!=fBandArray.end();++band) { // VIM
    
		if(0==band->GetKeepBand()) continue; // We are interested only in the selected bands (and may the force be with you :) )

		alpha = band->GetAlphaRun(); // alpha run
		alphaplus = alpha+1.;
		alphaless = alpha-1.;

		// Here I create an  intermediate step, in order to regroup the bin above a certain energy, in order to increase the per-bin statistic !
		START::Band RegroupBand = *band;
		// First pass, I compute all the expected excess : 
		Double_t ELimit = 3.; // Above ELimit TeV we assume that this behave like a PWL !
		Double_t ExpectedIndex = 2.7; // Spectral index assumed, above ELimit;
		Double_t nON_AboveELimit=0.;
		Double_t nOFF_AboveELimit=0.;
		Double_t S_AboveELimit=0.;
		Double_t p_AboveELimit=0.;
		Double_t AccTot_AboveELimit=0.; // This is not the acceptance, but just the particle flux 
    
		std::vector<EnergyBin>::const_iterator orig_bin = band->ebin.begin();
		for (std::vector<EnergyBin>::iterator bin = RegroupBand.ebin.begin();bin!=RegroupBand.ebin.end();++bin,++orig_bin) {
      
			if (bin->GetKeepBin() == 0) continue; // filtre pour bins en énergie en dessous du seuil
			double S = fCompRes->FunctionExpectedExcess(const_cast<Band*>(&*band),const_cast<EnergyBin*>(&*orig_bin),par); 
			if(fHypothesis->GetVectorNormalizedParameters().size()>0) { S*=normalization; }
			nON = bin->GetOn();
			nOFF = bin->GetOff();
			p = STARTUtils::GetExpectedOff(nON,nOFF,alpha,S);
			bin->SetSth(S);
			bin->SetOffFitted(p);
			double  acceff = TMath::Power( bin->GetEmin(), 1.-ExpectedIndex) - TMath::Power( bin->GetEmax(), 1.-ExpectedIndex) ;
			bin->SetAcceff( acceff );
      
			if (bin->GetEmin()>ELimit) {
				nON_AboveELimit+=nON;
				nOFF_AboveELimit+=nOFF;
				S_AboveELimit+=S;
				AccTot_AboveELimit+=acceff;
			}
		}
    
		p_AboveELimit = STARTUtils::GetExpectedOff(nON_AboveELimit,nOFF_AboveELimit,alpha,S_AboveELimit);
		// We loop again :
		for (std::vector<EnergyBin>::iterator bin = RegroupBand.ebin.begin();bin!=RegroupBand.ebin.end();++bin) {
			if (bin->GetKeepBin() == 0) continue; // filtre pour bins en énergie en dessous du seuil
			if (bin->GetEmin()>ELimit) {
				double p_est = p_AboveELimit * bin->GetAcceff() / AccTot_AboveELimit;
				bin->SetOffFitted(p_est);
			}
		}
    
		std::vector<EnergyBin>::const_iterator regroupband_bin = RegroupBand.ebin.begin();
		for (std::vector<EnergyBin>::const_iterator bin = band->ebin.begin();bin!=band->ebin.end();++bin,++regroupband_bin) {
      
			if (bin->GetKeepBin() == 0) continue; // filtre pour bins en énergie en dessous du seuil
      
			/* double S = fCompRes->FunctionExpectedExcess(const_cast<Band*>(&*band),const_cast<EnergyBin*>(&*bin),par); 
			   if(fHypothesis->GetVectorNormalizedParameters().size()>0) S*=normalization;
	 
			   nON = bin->GetOn();
			   nOFF = bin->GetOff();
	 
			   // JLK : S can't be < 0 by definition, see classes MinimizeFactory and Hypothesis
			   // so if it's happening, there is a (nice) bug or you're doing contours/scans and
			   // in these cases, it's allowed
	 
			   if (S<0 && fverbose)
			   isunphysical = true;
	 
			   p = STARTUtils::GetExpectedOff(nON,nOFF,alpha,S);
			*/
      
			double S = regroupband_bin->GetSth();
			nON = bin->GetOn();
			nOFF = bin->GetOff();
			p = regroupband_bin->GetOffFitted();
      
			if (S<0 && fverbose) {
				isunphysical = true;
			}
      
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
		std::streamsize ssp = std::cout.precision();
		INFO << "fcnvalue = " << std::setprecision(std::numeric_limits<double>::digits10)
			 << fcnvalue << std::setprecision(ssp) << ", time to do it : " << wtimeoneinter.RealTime()
			 << "s" << std::endl;
		// INFO << "fcnvalue = " << -fcnvalue << ", time to do it : " << wtimeoneinter.RealTime()
		// 	 << "s" << std::endl;
		if(isunphysical || isthereNaN) WARNING << "Unphysical region touched : " << std::endl;
		for(unsigned int i(0); i<par.size(); i++) {
			std::cout << fHypothesis->GetParametersNames()[i] << " = " << fHypothesis->GetParameters()[i] << std::endl;
		}
	}
  
	return fcnvalue;
}

#elif FITMETHOD == 3

/**
 * \brief Make a single band out of all bands. Method ULTRA BOURRIN !!!!!!!!
 * EXTREMELY EXPERIMENTAL
 * WHAT ABOUT ALPHA !
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
		if(TMath::IsNaN(fHypothesis->GetParameters()[i])) isthereNaN=true;
	}

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

	for(std::vector<Band>::const_iterator band = fBandArray.begin();band!=fBandArray.end();++band) {
   
		if(0==band->GetKeepBand()) continue; // We are interested only in the selected bands (and may the force be with you :) )
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
			sum_bin->SetAlpha( sum_bin->GetAlpha()  + bin->GetAlpha() * bin->GetOff() ); // For the moment, Alpha = Alpha*NOff (in order to compute the average )
			sum_bin->SetLiveTime( sum_bin->GetLiveTime() + bin->GetLiveTime() * bin->GetAlpha() ); // For the moment, livetime = Alpha *Livetime, for the average alpha
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

	std::vector<Band> VecBandDeSagouin;
	VecBandDeSagouin.push_back(SumBand);
  
	//for(std::vector<Band>::const_iterator band = fBandArray.begin();band!=fBandArray.end();++band) { // VIM
	for(std::vector<Band>::const_iterator band = VecBandDeSagouin.begin();band!=VecBandDeSagouin.end();++band) { // VIM
    
		if(0==band->GetKeepBand()) continue; // We are interested only in the selected bands (and may the force be with you :) )
    
		alpha = band->GetAlphaRun(); // alpha run
		alphaplus = alpha+1.;
		alphaless = alpha-1.;

    
		for (std::vector<EnergyBin>::const_iterator bin = band->ebin.begin();bin!=band->ebin.end();++bin) {
      
			if (bin->GetKeepBin() == 0) continue; // filtre pour bins en énergie en dessous du seuil
      
			// double S = fCompRes->FunctionExpectedExcess(const_cast<Band*>(&*band),const_cast<EnergyBin*>(&*bin),par); 
			// if(fHypothesis->GetVectorNormalizedParameters().size()>0) S*=normalization;
			double S = bin->GetSth();
			nON = bin->GetOn();
			nOFF = bin->GetOff();
      
			// JLK : S can't be < 0 by definition, see classes MinimizeFactory and Hypothesis
			// so if it's happening, there is a (nice) bug or you're doing contours/scans and
			// in these cases, it's allowed
      
			if (S<0 && fverbose)
				isunphysical = true;
      
			p = STARTUtils::GetExpectedOff(nON,nOFF,alpha,S);
      
			if (nON==0.)
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
      
			if (p>0.) {
				calc2 = TMath::Log(p);
				calc2*= nOFF;
			} 
			else calc2=0.;
      
			calc3 = p;
			calc3*= alphaplus;
      
			fcnvaluetmp=0.;
			fcnvaluetmp+= calc1;
			fcnvaluetmp+= calc2;
			fcnvaluetmp-= calc3;
			fcnvaluetmp-= S;
			if (nON>0.) {
				fcnvaluetmp-=(nON*TMath::Log(nON)-nON);
			}
			if (nOFF>0.) {
				fcnvaluetmp-=(nOFF*TMath::Log(nOFF)-nOFF);
			}
			fcnvalue+=fcnvaluetmp;
      
		} // loop on bin
    
	} // loop on band
  
	fcnvalue*= -1.; // MINUIT2 does minimization
  
	wtimeoneinter.Stop();
  
	if(fverbose) {
		std::streamsize ssp = std::cout.precision();
		INFO << "fcnvalue = " << std::setprecision(std::numeric_limits<double>::digits10) << fcnvalue << std::setprecision(ssp) << ", time to do it : " << wtimeoneinter.RealTime()
			 << "s" << std::endl;
		// INFO << "fcnvalue = " << -fcnvalue << ", time to do it : " << wtimeoneinter.RealTime()
		// 	 << "s" << std::endl;
		if(isunphysical || isthereNaN) WARNING << "Unphysical region touched : " << std::endl;
		for(unsigned int i(0); i<par.size(); i++) {
			std::cout << fHypothesis->GetParametersNames()[i] << " = " << fHypothesis->GetParameters()[i] << std::endl;
		}
	}
  
	return fcnvalue;
}


#else

double START::FCNLikelihood::operator()(const std::vector<double> &par) const
{
	return 0.;
}

#endif
