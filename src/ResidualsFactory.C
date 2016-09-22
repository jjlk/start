/**
 * \brief Compute residuals, flux, associated errors and mean energy bin for the graphical representation.
 * 
 * Residuals are computed with the Rolke method and it gives errors at 1 and 3 sigmas.
 *
 * \author HAP-Fr team
 */


// STL
#include <iostream>
#include <cstdlib>
#include <sstream>
// ROOT
#include <TMath.h>
#include <TRolke.h>

// START
#include "ResidualsFactory.hh"
#include "ComputeResults.hh"
#include "Hypothesis.hh"
#include "Band.hh"
#include "Residuals.hh"
#include "STARTUtils.hh"

// Utilities
#define DEBUG 0
#include "debugging.hh"

#define INFO std::cout << INFOCOLOR << "ResidualsFactory> " << RESETCOLOR
#define WARNING std::cout << WARNINGCOLOR << "ResidualsFactory> " << RESETCOLOR

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::ResidualsFactory)
#endif 


//#define DEBUG

/**
 * \brief Constructor needed by ROOT
 */
START::ResidualsFactory::ResidualsFactory()
{

	fHypothesisArray.clear();

	fresiduals.clear();
	fresiduals_on.clear();
	fresiduals_off.clear();
	fresidualssigmaplus.clear();
	fresidualssigmaminus.clear();
	fresiduals3sigmaplus.clear();
	fresiduals3sigmaminus.clear();
	fflux.clear();
	ffluxsigmaplus.clear();
	ffluxsigmaminus.clear();
	fflux3sigmaplus.clear();
	fflux3sigmaminus.clear();
	femean.clear();

}

/**
 * \brief Constructor for one hypothesis
 */
START::ResidualsFactory::ResidualsFactory(Hypothesis &hypo, bool verbose)
	:fverbose(verbose)
{

	fHypothesisArray.clear();

	fHypothesisArray.push_back(&hypo);
 
	fresiduals.clear();
	fresiduals_on.clear();
	fresiduals_off.clear();
	fresidualssigmaplus.clear();
	fresidualssigmaminus.clear();
	fresiduals3sigmaplus.clear();
	fresiduals3sigmaminus.clear();
	fflux.clear();
	ffluxsigmaplus.clear();
	ffluxsigmaminus.clear();
	fflux3sigmaplus.clear();
	fflux3sigmaminus.clear();
	femean.clear();

}

/**
 * \brief Constructor for more than one hypothesis
 */
START::ResidualsFactory::ResidualsFactory(std::vector<Hypothesis*> &hypo, bool verbose)
	:fverbose(verbose)
{

	fHypothesisArray.clear();

	for(std::vector<Hypothesis*>::iterator it=hypo.begin(); it!=hypo.end(); ++it) {
		fHypothesisArray.push_back(&(**it));
	}
 
	fresiduals.clear();
	fresiduals_on.clear();
	fresiduals_off.clear();
	fresidualssigmaplus.clear();
	fresidualssigmaminus.clear();
	fresiduals3sigmaplus.clear();
	fresiduals3sigmaminus.clear();
	fflux.clear();
	ffluxsigmaplus.clear();
	ffluxsigmaminus.clear();
	fflux3sigmaplus.clear();
	fflux3sigmaminus.clear();
	femean.clear();

}

/**
 * \brief Destructor
 */
START::ResidualsFactory::~ResidualsFactory()
{
	// No pointers
}


/**
 * \brief Compute Residuals with the Rolke method.
 *
 * \todo Study the alternative possibility to compute alpha when there is no off data (see comments inside the code)
 */

void START::ResidualsFactory::ComputeResidualsRolke()
{

	for(std::vector<Hypothesis*>::iterator hypo=fHypothesisArray.begin(); hypo!=fHypothesisArray.end(); hypo++) {

		DEBUG_OUT << "for hypothesis " << (*hypo)->GetName() << std::endl;

		if((*hypo)->GetConvergence()) {

			std::vector<Band> BandArray;
			BandArray = (*hypo)->GetBandArray();

			fresiduals.clear();
			fresiduals_on.clear();
			fresiduals_off.clear();
			fresidualssigmaplus.clear();
			fresidualssigmaminus.clear();
			fresiduals3sigmaplus.clear();
			fresiduals3sigmaminus.clear();
			fflux.clear();
			ffluxsigmaplus.clear();
			ffluxsigmaminus.clear();
			fflux3sigmaplus.clear();
			fflux3sigmaminus.clear();
			femean.clear();

			std::vector<bool> isupperlimits;
			int NbBin = BandArray[0].ebin.size();
			std::vector<double> energyminus;
			std::vector<double> energyplus;

			// JLK : ajout d'une mini-securite
			if(NbBin==0) {
				WARNING << "You have to initialize Bands before using them!!!" << std::endl;
				continue;
			}
      
			TRolke ConfInter;
			ConfInter.SetCL(0.6827); // 1 sigma
			ConfInter.SetBounding(true); // VIM : Set Bounding method to true (Rolke 2005) --> yield a bit more overcoverage when significant negative excess (but give results more according to our feeling)
      
			TRolke ConfInter3;
			ConfInter3.SetCL(0.9973); // 3 sigma
			ConfInter3.SetBounding(true); // VIM : Set Bounding method to true (Rolke 2005) --> yield a bit more overcoverage when significant negative excess (but give results more according to our feeling)

			for(int iem=0; iem<NbBin; iem++) {

				bool Eavailable=false;
				double sum_on(0);
				double sum_off(0);
				double sum_excess(0);
				double sum_excess_th(0);
				double sum_time(0);
				double alpha_mean(0);
				double alpha_mean_backup(0);
	
				double e_mean(-1);    
				double e_min(-1);    
				double e_max(-1);    
				int nband=0;
	
				for(unsigned int iband=0; iband<BandArray.size(); iband++) {

					if(BandArray[iband].GetKeepBand()==0) continue; // We are interested only in the selected bands
	  
					//For a given energy, use only bands for which this energy is above the threshold
					if (BandArray[iband].ebin[iem].GetKeepBin()==0) continue;
	  
					nband++;
	  
					double emean = BandArray[iband].ebin[iem].GetEmean();
					double emin = BandArray[iband].ebin[iem].GetEmin();
					double emax = BandArray[iband].ebin[iem].GetEmax();

					if (e_mean<0) { // for all bands
						e_mean=emean;
						e_min=emin;
						e_max=emax;
					}
					if (e_mean<0 && e_mean!=emean) {
						std::cout<<"Something wrong with emean... Skip hypothesis (TBC)"<<std::endl;
						exit(EXIT_FAILURE);
					}
	  
					double non = BandArray[iband].ebin[iem].GetOn();
					double noff = BandArray[iband].ebin[iem].GetOff();
					//double alpha = BandArray[iband].GetAlphaRun();
					double alpha = BandArray[iband].ebin[iem].GetAlpha();
					double sth = BandArray[iband].ebin[iem].GetSth();
					//double time = BandArray[iband].GetLiveTime();
					double time = BandArray[iband].ebin[iem].GetLiveTime();

					sum_on += non;
					sum_off += noff;
					sum_excess += (non - alpha*noff);	  
					sum_excess_th += sth;
					sum_time+=time;
					alpha_mean += (alpha*noff);
					/*
					  VIM : What worries me in this way (above) to compute the alpha_mean is that
					  if you have one band with noff=0, then you can biais you result.
					  i.e : If you have 2 runs 1 with an offset of 0.5, but 0 off, and the other 2.0 with N off,
					  then the result will be alpha_mean = (alpha05deg*0+alpha2deg*N)/(N+0) = alpha2deg
					*/
					alpha_mean_backup += (alpha*time); //for background killers or HE
					// VIM & ADA : Do we have to make the average using the off exposure ? But how we do it? Use the gamma  tables ? (assuming gamma=electron=hadrons ? 
				} // end loop on band
	
				if (e_mean>0) Eavailable=true; //AFTER the first loop over bands for a given energy
	
				if (Eavailable) {
	  
					femean.push_back(e_mean);
					energyminus.push_back(e_min);
					energyplus.push_back(e_max);

					DEBUG_OUT_L(2) << "alpha_mean=" << alpha_mean << std::endl;
					DEBUG_OUT_L(2) << "sum_off=" <<  sum_off << std::endl;
					DEBUG_OUT_L(2) << "alpha_mean_backup=" << alpha_mean_backup << std::endl;
					DEBUG_OUT_L(2) << "sum_time=" <<  sum_time << std::endl;
	  
					if (sum_off) alpha_mean /= sum_off;
					if (sum_time) alpha_mean_backup /= sum_time;
	  
					if (alpha_mean==0) alpha_mean=alpha_mean_backup;

					//residuals = excess_exp/excess_th

					double residual_exp(0);
					residual_exp=sum_excess/sum_excess_th;
					//fresiduals.push_back(residual_exp);
					if(residual_exp>0.) fresiduals.push_back(residual_exp);
					else fresiduals.push_back(0.); // JLK a reflechir (transformation en upperlimit ou inversion on/off?)

					// JLK ADD FOR ADA
					double sum_exp_on = STARTUtils::GetExpectedOn(sum_on, sum_off, alpha_mean, sum_excess_th);
					fresiduals_on.push_back(sum_on/sum_exp_on);

					double sum_exp_off = STARTUtils::GetExpectedOff(sum_on, sum_off, alpha_mean, sum_excess_th);
					fresiduals_off.push_back(sum_off/sum_exp_off);
						
					double minns68;
					double maxns68;
					ConfInter.SetPoissonBkgKnownEff((int)sum_on,(int)sum_off,1.0/alpha_mean,1.0);    
					ConfInter.GetLimits(minns68,maxns68);
	  
					double minns99;
					double maxns99;
					ConfInter3.SetPoissonBkgKnownEff((int)sum_on,(int)sum_off,1.0/alpha_mean,1.0);    
					ConfInter3.GetLimits(minns99,maxns99);	

					double residual_sigma_minus(0.);
					double residual_sigma_plus(0.);
					residual_sigma_minus=minns68/sum_excess_th;
					residual_sigma_plus=maxns68/sum_excess_th;
					fresidualssigmaplus.push_back(residual_sigma_plus);
					fresidualssigmaminus.push_back(residual_sigma_minus);
	  
					double residual_3sigma_minus(0.);
					double residual_3sigma_plus(0.); 
					residual_3sigma_minus=minns99/sum_excess_th;
					residual_3sigma_plus=maxns99/sum_excess_th;
					fresiduals3sigmaplus.push_back(residual_3sigma_plus);
					fresiduals3sigmaminus.push_back(residual_3sigma_minus);	
	  
					double flux_th = (*hypo)->GetFluxFitParams(e_mean);
					fflux.push_back(flux_th*residual_exp);
	  
					ffluxsigmaplus.push_back(residual_sigma_plus*flux_th);
					ffluxsigmaminus.push_back(residual_sigma_minus*flux_th);
	  
					fflux3sigmaplus.push_back(residual_3sigma_plus*flux_th);
					fflux3sigmaminus.push_back(residual_3sigma_minus*flux_th);

					if(ffluxsigmaminus.back()==0. || fflux.back()<0.) {
						isupperlimits.push_back(true);
					}
					else isupperlimits.push_back(false);

					DEBUG_OUT << "residual_exp=" << residual_exp << " sum_excess=" << sum_excess 
							  << " sum_excess_th=" << sum_excess_th 
							  << "residual_3sigma_plus=" << residual_3sigma_plus << std::endl;
	  
					DEBUG_OUT << "e_mean" << " " << "residual_exp" << " " << "residual_sigma_minus" << " " << "residual_sigma_plus + flux" << std::endl;
					DEBUG_OUT << "1sigma : " << e_mean << " " << residual_exp -1 << " " << residual_sigma_minus -1 << " " << residual_sigma_plus -1 << " " 
							  << flux_th*residual_exp << " " << residual_sigma_minus*flux_th << " " << residual_sigma_plus*flux_th << std::endl;
					DEBUG_OUT << "3sigma : " << e_mean << " " << residual_exp -1 << " " << residual_3sigma_minus -1 << " " << residual_3sigma_plus -1 << " "
							  << flux_th*residual_exp << " " << residual_3sigma_minus*flux_th << " " << residual_3sigma_plus*flux_th <<std::endl;
					//ADA For SED plotting : Energies in MeV, vFv in MeV/cm2/s 
					std::cout << "---------------------------------------------------Flux ULs in vFv --------------------------------------------------------------"<< std::endl; 
					std::cout <<  e_min*1e6 << " " << e_mean*1e6 << " " << e_max*1e6 << " " << residual_3sigma_plus*flux_th*(e_mean*1e6*e_mean*1e6) << std::endl;
	  
	  
	  
				} //end if	
	
			} //end loop over energies
      
			(*hypo)->SetResiduals(fresiduals);
			// JLK ADD FOR ADA
			(*hypo)->SetResidualsOn(fresiduals_on);
			(*hypo)->SetResidualsOff(fresiduals_off);
			(*hypo)->SetResidualsSigmaPlus(fresidualssigmaplus);
			(*hypo)->SetResidualsSigmaMinus(fresidualssigmaminus);
			(*hypo)->SetResiduals3SigmaPlus(fresiduals3sigmaplus);
			(*hypo)->SetResiduals3SigmaMinus(fresiduals3sigmaminus);
			(*hypo)->SetResidualsFlux(fflux);
			(*hypo)->SetResidualsFluxSigmaPlus(ffluxsigmaplus);
			(*hypo)->SetResidualsFluxSigmaMinus(ffluxsigmaminus);
			(*hypo)->SetResidualsFlux3SigmaPlus(fflux3sigmaplus);
			(*hypo)->SetResidualsFlux3SigmaMinus(fflux3sigmaminus);
			(*hypo)->SetMeanEnergy(femean);
      
			Residuals *Res = new Residuals();
			Res->SetEnergy(femean);
			Res->SetEnergyPlus(energyplus);
			Res->SetEnergyMinus(energyminus);
			Res->SetResiduals(fresiduals);
			Res->SetResidualsSigmaMinus(fresidualssigmaminus);
			Res->SetResidualsSigmaPlus(fresidualssigmaplus);
			Res->SetResiduals3SigmaMinus(fresiduals3sigmaminus);
			Res->SetResiduals3SigmaPlus(fresiduals3sigmaplus);
			Res->SetFlux(fflux);
			Res->SetFluxSigmaMinus(ffluxsigmaminus);
			Res->SetFluxSigmaPlus(ffluxsigmaplus);
			Res->SetFlux3SigmaMinus(fflux3sigmaminus);
			Res->SetFlux3SigmaPlus(fflux3sigmaplus);
			Res->SetIsUpperLimits(isupperlimits);
			Res->SetResidualsOn(fresiduals_on);
			Res->SetResidualsOff(fresiduals_off);
			(*hypo)->AddResiduals(*Res);
      
			INFO << "Residuals for "<< (*hypo)->GetName() << std::endl;
			Res->Print();

			delete Res; Res = 0;

		}
		else {
			WARNING << "Skip hypothesis " << (*hypo)->GetName() << " because it didn't converge!" << std::endl;
			continue;
		}
    
    
	} // end loop on hypothesis
  
	INFO << "Computing residuals... ok" << std::endl;

}



/*
 * Determine flux points, residuals and both errors using the
 * same method as in spectrumIR: the mesured flux is the mean
 * value of the measured fluxes for couples (band i, energy bin j)
 * with weights w_ij=1/sigma_ij**2, where sigma_ij are the gaussian 
 * assumed errors on the ij measured fluxes.
 * If FittedOnOff=false, sigma_ij is determined with measured numbers
 * If FittedOnOff=true, sigma_ij is determined with fitted numbers
 *
 * \warning don't use it
 */

/*
  void START::ResidualsFactory::ComputeResidualsGaus(bool FittedOnOff) 
  {

  for(std::vector<Hypothesis*>::iterator hypo=fHypothesisArray.begin(); hypo!=fHypothesisArray.end(); hypo++) {

  if((*hypo)->GetIsMinimumValid()==true && (*hypo)->GetAreFittedParametersValid()==true && (*hypo)->GetIsCovarianceValid()==true) {

  fBandArray = (*hypo)->GetBandArray();
      
  int NbBin = fBandArray[0].ebin.size();
      
  // JLK : ajout d'une mini-securite
      
  if(NbBin==0) {
  std::cout << "Residuals.C : You have to initialize the Bands Before used them!!! ==> EXIT" << std::endl;
  exit(EXIT_FAILURE);
  }
      
  for(unsigned int iem=0; iem<NbBin; iem++) {
	
  double flux(0);
  double residual(0);
  double dflux(0);
  double dresidual(0);
  double sum_weights(0);
  double sum_excess(0);
  double alpha_mean(0);
  double tot_time(0);
  double e_mean(-1);
	
  int sum_non(0);
  int sum_noff(0);
	
  for(unsigned int iband=0; iband<fBandArray.size(); iband++) {
	  
  if(fBandArray[iband].GetKeepBand()==0) continue; // We are interested only in the selected bands
	  
  //For a given energy, use only bands for which this energy is above the threshold
  if (fBandArray[iband].ebin[iem].GetKeepBin() == 0) continue;
	  
  double emin = fBandArray[iband].ebin[iem].GetEmin();
  double emean = fBandArray[iband].ebin[iem].GetEmean();
	  
  if (e_mean<0) e_mean=emean; // for all bands
  if (e_mean<0 && e_mean!=emean) {
  std::cout<<"Something wrong with emean... Exit"<<std::endl;
  exit(EXIT_FAILURE);
  }
	  
  double non = fBandArray[iband].ebin[iem].GetOn();
  double noff = fBandArray[iband].ebin[iem].GetOff();
  double non_fitted = fBandArray[iband].ebin[iem].GetOnFitted();
  double noff_fitted = fBandArray[iband].ebin[iem].GetOffFitted();
  double alpha = fBandArray[iband].GetAlphaRun();
  double time = fBandArray[iband].GetLiveTime();
	  
  alpha_mean+=(alpha*time);
  tot_time+=time;
  double excess = non - alpha*noff;
	  
  double sth = fBandArray[iband].ebin[iem].GetSth();
	  
  double fluxband_th = (*hypo)->GetFluxFitParams(emean); 
  double fluxband_exp = (excess/sth)*fluxband_th; 
	  
  double sigma_signal = 0;
  if (FittedOnOff) 
  sigma_signal = TMath::Sqrt(non_fitted+alpha*alpha*noff_fitted);
  else
  sigma_signal = TMath::Sqrt(non+alpha*alpha*noff);
	  
  double sigma=(fluxband_th/sth)*sigma_signal;
	  
  double weight = 0;
	  
  //std::cout<<"sigma="<<sigma<<"   non="<<non<<"  noff="<<noff<<std::endl;
	  
  if (sigma>0) {
  weight = 1./(sigma*sigma);
  flux += (weight * fluxband_exp);
  sum_weights += weight;
  sum_non+=(int)non;  //SP: pquoi EnergyBin::fOn et fOff sont double et non pas int ?
  sum_noff+=(int)noff;
  sum_excess += excess;	
  } else {
  //std::cout<<"out because sigma null"<<std::endl;
  continue;
  }
	  
  }
	
  double flux_th = (*hypo)->GetFluxFitParams(e_mean);
	
  //std::cout<<"sum_noff="<<sum_noff<<" sum_weights="<<sum_weights<<std::endl;
	
  if (tot_time>0)
  alpha_mean/=tot_time;
  else
  alpha_mean=0;
	
  DEBUG_OUT << "sum_weights = " << sum_weights << std::endl;
	
  if (sum_weights) {
  flux/=sum_weights;
  dflux=1./(TMath::Sqrt(sum_weights));
	  
  residual= (flux - flux_th)/flux_th;
  dresidual=dflux/flux_th;
	  
  std::cout<<"Gaus: "<<e_mean<<" "<<sum_non<<" "<<sum_noff<<" "
  <<alpha_mean<<" "<<sum_excess<<" "<<flux<<" "<<dflux<<" "<<residual<<" "<<dresidual<<std::endl;
	  
  } // end loop band
	
  } // end loop ebin
      
  } 
  
  } // en loop hypothesis
  
  }
*/

/*
 * For each combination of i,j (i for band, j for energy) we consider the 
 * noff and alpha*noff values (measured or fitted?) as the mean values of two
 * different poisson distributions, respectively for the background in the
 * ON and OFF regions. Similarly, the expected signal in the bin will be used 
 * as the mean value of a third poisson distribution.
 * 
 * Then, N random experiments will give values of the excesses S'ij
 * => we get the distribution of the residual Rmc_j=(somme_i S'ij)/Sth_j
 *    ==> this distribution should be centered around 1 (otherwiser there 
 *    would be a bias to undestand and eventually to take into account).
 *    ==> we get the asymetrical errors on the residua d+ and d-
 *
 * The residual is calculated as R_j = S_j/Sth_j where S_j is the measured
 * excess in the energy bin (S_j = (somme_i)S_ij), the errors come from the MC
 *
 * The points and errors are derived from residuals using R_j=S_j/Sth_j=phi_j/phith_j
 *
 * The method to generate random number is func->GetRandom() which uses TRandom
 * => this could be bad for Poisson distributions (TRandom3 seems to be better)
 * To be improved if necessary.
 *
 * \warning don't use it
 */

/*
  void START::ResidualsFactory::ComputeResidualsPoissonMC() 
  {

  for(std::vector<Hypothesis*>::iterator hypo=fHypothesisArray.begin(); hypo!=fHypothesisArray.end(); hypo++) {

  fBandArray = (*hypo)->GetBandArray();

  int NbBin = fBandArray[0].ebin.size();

  // JLK : ajout d'une mini-securite
  if(NbBin==0) {
  std::cout << "Residuals.C : You have to initialise Bands before using them!!! ==> EXIT" << std::endl;
  exit(EXIT_FAILURE);
  }
    
  TRandom3 rg(0); // mersenne twister
    
  for(unsigned int iem=0; iem<NbBin; iem++) {
  //for(unsigned int iem=9 ; iem<10; iem++) {
  //for(unsigned int iem=8 ; iem<9; iem++) {
      
  bool Eavailable=false;
  double sum_excess(0);
  double sum_excess_th(0);
  double e_mean(-1);
      
  TString hname="hresMCforE_";
  hname+=iem;
  TH1F *hresMCforE = new TH1F(hname.Data(),hname.Data(),10000,-10.,10.);
      
  int NMC=10000; //number loops over bands generating MC numbers
  for (int imc(0); imc<NMC ; imc++) {
	
  double sum_rmd_excess(0);
	
  for(unsigned int iband=0; iband<fBandArray.size(); iband++) {
	  
  if(fBandArray[iband].GetKeepBand()==0) continue; // We are interested only in the selected bands
	  
  //For a given energy, use only bands for which this energy is above the threshold
  if (fBandArray[iband].ebin[iem].GetKeepBin() == 0) continue;
	  
  double emean = fBandArray[iband].ebin[iem].GetEmean();
	  
  if (e_mean<0) e_mean=emean; // for all bands
  if (e_mean<0 && e_mean!=emean) {
  std::cout<<"Something wrong with emean... Exit"<<std::endl;
  exit(EXIT_FAILURE);
  }
	  
  double non = fBandArray[iband].ebin[iem].GetOn();
  double noff = fBandArray[iband].ebin[iem].GetOff();
  double non_fitted = fBandArray[iband].ebin[iem].GetOnFitted();
  double noff_fitted = fBandArray[iband].ebin[iem].GetOffFitted();
  double alpha = fBandArray[iband].GetAlphaRun();
  double sth = fBandArray[iband].ebin[iem].GetSth();
	  
  double rnd_non;
  double rnd_noff;
	  
  rnd_noff = rg.PoissonD(noff_fitted); // Noff "seen" 
  rnd_non  = rg.PoissonD(non_fitted); // Non "seen"
  double rnd_excess = rnd_non-alpha*rnd_noff;
	  
	  
  if (imc == 0) {
  double excess = non - alpha*noff;
  sum_excess += excess;	  
  sum_excess_th += sth;
  }
	  
  sum_rmd_excess += rnd_excess;
	  
  }
	
	
  hresMCforE->Fill(sum_rmd_excess/sum_excess_th);
	
	
  if (e_mean>0) Eavailable=true;//AFTER the first loop over bands for a given energy
	
  }
      
      
  if (Eavailable) {
	
  if (sum_excess>0) {
	  
  double residual_mc_mean = hresMCforE->GetMean();
  double residual_mc_sigma = hresMCforE->GetRMS();
	  
  TFile g("new.root","UPDATE");
  hresMCforE->Write();
  g.Close();
	  
  int nbins = hresMCforE->GetNbinsX();
  int binmean = hresMCforE->GetXaxis()->FindBin(residual_mc_mean);
	  
  int ibm=0;    
  int ibp=0;    
  double half68 = 0.5*0.6827;
	  
  for (int ib = 1 ;  ib<binmean ;  ib++) {
  double fraction = hresMCforE->Integral(ib,binmean)/hresMCforE->Integral() ;
  if (ib==1 && fraction < half68) {
  std::cout <<"Residuals: impossible to get sigma-"<<std::endl;	
  }
	    
  if (fraction < half68) {
  ibm=ib;
  break;
  }	
  }
	  
	  
  for (int ib = binmean ;  ib<=nbins ;  ib++) {
  double fraction = hresMCforE->Integral(binmean,ib)/hresMCforE->Integral() ;
	    
  if (fraction > half68) {
  ibp=ib;
  break;
  }
  }
	  

  double x1=hresMCforE->GetBinCenter(ibm-1);
  double y1=hresMCforE->Integral(ibm-1,binmean)/hresMCforE->Integral();
  double x2=hresMCforE->GetBinCenter(ibm);
  double y2=hresMCforE->Integral(ibm,binmean)/hresMCforE->Integral();
  double xm_interpol=x1+(x2-x1)*(half68-y1)/(y2-y1); //for y==(1/2)*68.27%
  #ifdef DEBUG 
  std::cout <<"-binA "<<ibm-1<<"  BinCenter="<<x1<<" histo="<<y1<<std::endl ;
  std::cout <<"-binB "<<ibm<<  "  BinCenter="<<x2<<" histo="<<y2<<std::endl ;
  std::cout <<"Interpolation: half68 reached at " <<xm_interpol<<std::endl;
  std::cout <<" "<<std::endl;
  #endif      

  x1=hresMCforE->GetBinCenter(ibp);
  y1=hresMCforE->Integral(binmean,ibp)/hresMCforE->Integral();
  x2=hresMCforE->GetBinCenter(ibp+1);
  y2=hresMCforE->Integral(binmean,ibp+1)/hresMCforE->Integral();
  double xp_interpol=x1+(x2-x1)*(half68-y1)/(y2-y1); //for y==(1/2)*68.27%
  #ifdef DEBUG 
  std::cout <<"+binA "<<ibp<<"  BinCenter="<<x1<<" histo="<<y1<<std::endl ;
  std::cout <<"+binB "<<ibp+1<<  "  BinCenter="<<x2<<" histo="<<y2<<std::endl ;
  std::cout <<"Interpolation: half68 reached at " <<xp_interpol<<std::endl;
  std::cout <<" "<<std::endl;
  #endif      
	  
  double residual_mc_sigma_minus=residual_mc_mean - xm_interpol;
  double residual_mc_sigma_plus=xp_interpol - residual_mc_mean;
	  
  hresMCforE->Delete();
	  
  double residual_exp = sum_excess/sum_excess_th;
	  
  double flux_th = (*hypo)->GetFluxFitParams(e_mean);
  double flux = flux_th * residual_exp;
	  
  std::cout<<"MC: "<<e_mean<<" "<<sum_excess<<" "<<sum_excess_th<<" "
  <<residual_exp<<" "<<residual_mc_mean<<" "<<residual_mc_sigma
  <<" "<<residual_mc_sigma_minus<<" "<<residual_mc_sigma_plus<<" "
  <<flux<<" "<<std::endl;
  } else {
  std::cout <<"MC: "<<e_mean<<" negative or zero excess"<<std::endl;
  } //end if over sum_excess
	
  } //end if over Eavailable
      
  } //end loop over energies
    
  } // end loop over hypothesis
  
  }


*/






/*

  TESTS:
  On recupere:   e_mean   residual_exp   residual_mc_mean   residual_mc_sigma
  more toto | grep -v "negative" | awk '{print $2"   "$5"   "$6"   "$7}' > res

  gnuplot> set logscale x
  gnuplot> plot 'resg' u 1:8:9 with yerrorbars,\
  gnuplot> 'res' u 1:3

*/


