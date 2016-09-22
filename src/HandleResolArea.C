/**
 * \brief Look for and stock in bands effective area, resolutions and biais stored
 * in MC DSTs. Additional infos are then stock in bands and bins.
 * 
 * Histograms are stored in maps and their keys are defined by efficacity, 
 * zenith angle and offset, then interpolated areas, biais and resolution
 * are stored in bands to speed up the computations of the expected excess.
 *
 * The energy threshold is computed and then used as a flag for each bin
 * in energy to determine if the bins can be used in the minimization.
 *
 * \author HAP-Fr team
 */

// STL
#include <string>
#include <algorithm>
#include <sstream>
#include <cstdlib>

// ROOT
#include <TSystem.h>
#include <TMath.h>
#include <Math/InterpolationTypes.h>
#include <Math/Interpolator.h>
#include <TROOT.h>
#include <TFile.h>
#include <TLatex.h>
#include <TLegend.h>

// START
#include "ComputeResults.hh"
#include "HandleResolArea.hh"

// Utilities
#define DEBUG  0
#include "debugging.hh"

#define INFO std::cout << INFOCOLOR << "HandleResolArea> " << RESETCOLOR
#define WARNING std::cout << WARNINGCOLOR << "HandleResolArea> " << RESETCOLOR

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::HandleResolArea)
#endif

/**
 * \brief Constructor to build all instrument functions
 *
 * \param configname name of the configuration
 * \param forceuseofdistribution false for point like source and true for extended source
 * \param readnow if true read all IRF
 */
	START::HandleResolArea::HandleResolArea(const START::Config &Conf,
											START::STARTUtils::IrfOpt Type,
											bool readnow)
:fMcProductionName(Conf.GetUserMcProductionName()),
	fBandArray(0),
	fForceUseDistribution(false),
	fConfigName(""),
	fAreaMin(4.),
	fResolFactor(2.),
	fEnergyBinThresholdCondition(Safe),
	fMc(0)
{
	// JLK add
	fConfigName = (Conf.GetUserForceAreaAnalysisConfig().CompareTo("")==0) ?
		Conf.GetUserAnalysisConfig() :
		Conf.GetUserForceAreaAnalysisConfig();
	
	fSafeThresholdMethod = AreaAndBiaisResolMethod;
	fFractionMaxArea = 10.;
  
	fIsInstrumentTableAlreadyLoad=readnow;

	// JLK add
	bool forceuseofdistribution = (fConfigName.Contains("thsq64"))  ? true : false;
	
	if (fConfigName.Contains("thsq64")) {
		fForceUseDistribution=true;
	}
	fForceUseDistribution = (fForceUseDistribution || forceuseofdistribution);
	if (fForceUseDistribution) {
		INFO_OUT << "We'll use the distribution" << std::endl;
	}
	// VIM : Je n'aime pas trop les quelques lignes au dessus.
  
	fMc = new START::MonteCarlo(Type); // JLK add
	fEffMin = fMc->GetEfficiency().front();
	fEffMax = fMc->GetEfficiency().back();
	fOffMin = fMc->GetOffset().front();
	fOffMax = fMc->GetOffset().back();
	fZenMin = fMc->GetZenith().front();
	fZenMax = fMc->GetZenith().back();

	//std::cout << "     effmin " << fEffMin << std::endl;	  
	//std::cout << "     effmax " << fEffMax << std::endl;	  
	//MPA: they are default values!

	if (readnow) {
		if(ReadAndStore()==0) {
			INFO << "Reading Instrument's functions... ok" << std::endl;
		}
		else {
			WARNING << "Reading Instrument's functions... FAILED" << std::endl;
			fIsInstrumentTableAlreadyLoad = false;
			exit(EXIT_FAILURE);
		}
	}
}


/**
 * \brief Destructor
 */

START::HandleResolArea::~HandleResolArea() 
{
	delete fMc; fMc = 0;
	
	ClearTablesVector();
}

/**
 * \brief Clear vectors filled by ReadAndStore
 */ 
void START::HandleResolArea::ClearTablesVector() {
	for(std::map<TString, TH1F *>::iterator it=fmapHistoResol.begin(); it!=fmapHistoResol.end(); ++it) {
		if (it->second!=0) delete it->second;
		it->second=0; // JLK
	}
	fmapHistoResol.clear();
  
	for(std::map<TString, TH1F *>::iterator it=fmapHistoBiais.begin(); it!=fmapHistoBiais.end(); ++it) {
		if (it->second!=0) delete it->second;
		it->second=0; // JLK
	}
	fmapHistoBiais.clear();
  
	for(std::map<TString, TH1F *>::iterator it=fmapHistoArea.begin(); it!=fmapHistoArea.end(); ++it) {
		if (it->second!=0) delete it->second;
		it->second=0; // JLK
	}
	fmapHistoArea.clear();
  
	for (std::map<TString, TH1F *>::iterator it=fmapHistoDistrib.begin(); it!=fmapHistoDistrib.end(); ++it) {
		if (it->second!=0) delete it->second;
		it->second=0; // JLK
	}
	fmapHistoDistrib.clear();
    
}
  
/**
 * \brief Set infos, vectors, interpolators in bands and bins
 */
void START::HandleResolArea::SetAndUpdateBands(std::vector<Band> &BandArray)
{

	fBandArray = &BandArray;

	//GetNeededMcRange(*fBandArray,fEffMin,fEffMax,fOffMin,fOffMax,fZenMin,fZenMax); // JLK ADD (no DataSummary anymore)
	GetNeededMcRange(*fBandArray,fEffMin,fEffMax,fOffMin,fOffMax,fZenMin,fZenMax); // JLK ADD (no DataSummary anymore)
	std::cout << "     effmin " << fEffMin << std::endl;	  
	std::cout << "     effmax " << fEffMax << std::endl;	  
	std::cout << "     effmin " << fOffMin << std::endl;	  
	std::cout << "     effmax " << fOffMax << std::endl;	  
	std::cout << "     effmin " << fZenMin << std::endl;	  
	std::cout << "     effmax " << fZenMax << std::endl;
	



	//MPA: values are set here!
	//MPA: problem if eff is exactly at the border of a band!!!

	
	std::cout << "allé" << std::endl;	
	FillBandsWithDistributionBoolean(); // VIM : To tell to each band to take the distribution or the mean and sigma. Should be call before the call of safethreshold
	std::cout << "allée" << std::endl;
  
	if (!fIsInstrumentTableAlreadyLoad) {
		// VIM : If fIsInstrumentTableAlreadyLoad is true, this mean all the table have been already loaded, 
		//no need to have a selective reading of those tables in this case !
		// VIM : Otherwise we reload the table at each time we call this function
		ClearTablesVector();

		if(ReadAndStore()==0) { 
			INFO << "Reading Instrument's functions... ok" << "\033[0m" << std::endl;
		}
		else {
			WARNING << "Reading Instrument's functions... FAILED" << std::endl;
			exit(EXIT_FAILURE);
		}
	}
  
	if (CopyVectorsInBands()==-1) {
		WARNING << "Copying Vectors in Band... FAILED" << std::endl;
		exit(EXIT_FAILURE);
	}
	else {
		INFO << "Copying Vectors in Band... ok" << std::endl;
	}
  
	if (CopyInterpolatorsInBands()==-1) {
		WARNING << "Assignment of GSLInterpolators in Band... FAILED" << std::endl;
		exit(EXIT_FAILURE);
	}
	else {
		INFO << "Assignment of GSLInterpolators in Band... ok" << std::endl;
	}

	SetInBandsFirstEmcBin();
	SetInBandsMCEnergyThreshold();
	SetInBinsKeepBin();
	SetInBinsEffectiveArea();


	ComputeResults *CompRes = new ComputeResults(*fBandArray);

	if(CompRes->MakeVectorPartialIntegral(*fBandArray)==-1) {
		WARNING << "Copying effective resolution in Bins... FAILED" << std::endl;
		exit(EXIT_FAILURE);
	}
	else {
		INFO << "Copying effective resolution in Bins... ok" << std::endl;
	}
	delete CompRes;

	if (CopyInterpolatorsInBins()==-1) {
		WARNING << "Assignment of GSLInterpolators in Bins... FAILED" << std::endl;
		exit(EXIT_FAILURE);
	}
	else {
		INFO << "Assignment of GSLInterpolators in EnergyBins... ok" << std::endl;
	}

}

/**
 * \brief Store in the band vector the boolean use to choose between distribution or mean and sigma instrument function
 */
void START::HandleResolArea::FillBandsWithDistributionBoolean() 
{

	for (std::vector<Band>::iterator it = (*fBandArray).begin(); it != (*fBandArray).end(); ++it) {
		it->SetUseOfInstrumentEnergyDistribution(fForceUseDistribution);
	}

}


/**
 * \brief Read and store histograms in maps
 */
int START::HandleResolArea::ReadAndStore()
{

#if DEBUG
	std::cout << "Debug begining START::HandleResolArea::ReadAndStore()" << std::endl;
#endif

	// Get MC values
	//MonteCarlo MC;  
	std::vector<double> effMC = fMc->GetEfficiency();
	std::vector<double> offMC = fMc->GetOffset();
	std::vector<double> zenMC = fMc->GetZenith();
	std::vector<TString> telMC = fMc->GetTelType();
	std::vector<double> enMC = fMc->GetEnergy();

	std::vector<TString>::iterator itel_begin = telMC.begin();
	std::vector<TString>::iterator itel_end = telMC.end();

  
	for(std::vector<TString>::iterator itel=itel_begin; itel!=itel_end; ++itel) {   

		for(std::vector<double>::iterator ieff=effMC.begin(); ieff!=effMC.end(); ++ieff) {
      
			// Keep only the interesting efficiencies
      
			if((*ieff) < fEffMin*0.99) continue; // SP A voir
			if((*ieff) > fEffMax*1.01) continue;
			//MPA: these 2 lines above can skipp everything without even a warning!!!
      
			TString areapath  = gSystem->Getenv("HESSCONFIG");
			TString resolpath = gSystem->Getenv("HESSCONFIG");
			TString distrpath = gSystem->Getenv("HESSCONFIG");

			areapath+="/";
			areapath+=fConfigName;
			if(*itel!="") {
				areapath+="_";
			}
			areapath+=*itel + "/";

			resolpath+="/";
			resolpath+=fConfigName;
			if(*itel!="") {
				resolpath+="_";
			}
			resolpath+=*itel + "/";

			distrpath+="/";
			distrpath+=fConfigName;
			if(*itel!="") {
				distrpath+="_";
			}
			distrpath+=*itel + "/";

			areapath+=fMcProductionName; areapath+="_eff";
			resolpath+=fMcProductionName; resolpath+="_eff";
			distrpath+=fMcProductionName; distrpath+="_eff";
     
			std::ostringstream oss_eff; // VIM : Better way, in order to use C++ library only !
			oss_eff.precision(0);
			oss_eff << std::fixed << *ieff;
			TString seff(oss_eff.str().c_str());
			oss_eff.str("");
      
			areapath+=seff;
			resolpath+=seff;
			distrpath+=seff;
			areapath+="/";
			areapath+="CollectionArea.root";
			resolpath+="/";
			resolpath+="EnergyResolution.root";
			distrpath+="/";
			distrpath+="EnergyDistribution.root";

			// Opening files containing resol and area

			TFile *areafile = new TFile(areapath, "READ");
			gROOT->cd();
			if(!areafile->IsOpen() || !areafile) {
				std::cout << "Can't find :" << areafile->GetName() << " ==> exit!" << std::endl;
				delete areafile; areafile=0; // VIM : It still opens a TFile that needs to be deleted
				return -1;
			}
      
			TFile *resolfile = new TFile(resolpath, "READ");
			gROOT->cd();
			if(!resolfile->IsOpen() || !resolfile) {
				std::cout << "Can't find :" << resolfile->GetName() << " ==> exit!" << std::endl;
				delete resolfile; resolfile = 0; // VIM : It still opens a TFile that needs to be deleted
				return -1;
			}
      
			TFile *distrfile = 0;
			if(fForceUseDistribution) {
				distrfile = new TFile(distrpath, "READ");
				gROOT->cd();    
				if(!distrfile->IsOpen()) {
					delete distrfile; distrfile = 0; // VIM : It still opens a TFile that needs to be deleted
					if (fForceUseDistribution) {
						std::cout << "Can't find :" << distrpath << " ==> exit !" << std::endl;
						return -1;
					}
				}
			}
      
      
			for(std::vector<double>::iterator izen=zenMC.begin(); izen!=zenMC.end(); ++izen) {
	
				// Keep only the interesting zenith angles
	
				if((*izen) < fZenMin*0.99) continue; // SP A voir
				if((*izen) > fZenMax*1.01) continue;
	
				std::ostringstream oss_zen; // VIM : The C++ way !
				oss_zen.fill('0');
				oss_zen.precision(0);
				oss_zen.width(2);
				oss_zen << std::fixed << *izen;
				TString szen(oss_zen.str().c_str());
				oss_zen.str("");
	
	
				for(std::vector<double>::iterator ioff=offMC.begin(); ioff!=offMC.end(); ++ioff) {
	  
					// Keep only the interesting offsets
	  
					if((*ioff) < fOffMin*0.99) continue; // SP A voir
					if((*ioff) > fOffMax*1.01) continue;
	  
					std::ostringstream oss_off;
					oss_off.precision(1);
					oss_off << std::fixed << *ioff;
					TString soff(oss_off.str().c_str());
					oss_off.str("");

					//Make histo's names
					TString areahistoname = "EffArea_";
					areahistoname+=szen;
					areahistoname+="deg_";
					areahistoname+=soff;
					areahistoname+="off_eff";
					areahistoname+=seff;
					areahistoname+="_FixedE";
	  
					TString resolhistoname = "Resol_Sigma_";
					resolhistoname+=szen;
					resolhistoname+="deg_";
					resolhistoname+=soff;
					resolhistoname+="off_eff";
					resolhistoname+=seff;
					resolhistoname+="_FixedE";
	  
					TString biaishistoname = "Resol_Biais_";
					biaishistoname+=szen;
					biaishistoname+="deg_";
					biaishistoname+=soff;
					biaishistoname+="off_eff";
					biaishistoname+=seff;
					biaishistoname+="_FixedE";
	  
#if DEBUG
					std::cout << "distrpath :" << distrpath << std::endl;	  
					std::cout << "resolpath :" << resolpath << std::endl;
					std::cout << "areapath :"<< areapath << std::endl;
					std::cout << "areahistoname : "<< areahistoname << std::endl;
					std::cout << "resolhistoname :" << resolhistoname << std::endl;
					std::cout << "biaishistoname :" << biaishistoname << std::endl;
					std::cout << "key : " << MakeKeyRun(*ieff, *izen, *ioff, *itel) << std::endl;
#endif


					// We build the map's keys and the histograms associated
					// And the map pairs vectors
	  
					TString KeyRun = MakeKeyRun(*ieff, *izen, *ioff, *itel); // VIM : Try to improve readability
	  
					TH1F *areahisto_tmp = 0;
					areahisto_tmp = (TH1F *)areafile->Get(areahistoname);
					if (areahisto_tmp) {
						//fmapHistoArea[KeyRun] = new TH1F(*areahisto_tmp);
						areahisto_tmp->SetDirectory(0);
						fmapHistoArea[KeyRun] = areahisto_tmp;
						if (areahisto_tmp->GetEntries()==0) {
							//WARN_OUT << "\033[31;40;1m" << "fmapHistoArea[" << KeyRun << "]->GetEntries() = " << areahisto_tmp->GetEntries() << "\033[0m" << std::endl;
							std::cout << "\033[31;40;1m" << "fmapHistoArea[" << KeyRun << "]->GetEntries() = " << areahisto_tmp->GetEntries() 
									  << " --> THIS IS A BIG PROBLEM CHECK THE TABLE AND REPROCESS IT!!" << "\033[0m" << std::endl;
						}
					}
					else {
						fmapHistoArea[KeyRun] = 0;
						DEBUG_OUT << "VIM DEBUG> I couldn't find the histo named = " << areahistoname << std::endl;
					}

					if(fmapHistoArea[KeyRun]) {
						std::vector<double> x,y; // will contain energy and resolution
						MakeVectorsFromMap(fmapHistoArea,KeyRun,x,y);
						fmapVectorsArea[KeyRun] = std::make_pair(x,y);
						x.clear(); y.clear();
					}
					else {
						std::cout << "Can't find TH1F : " << areahistoname << " ==> exit" << std::endl;
						return -1;
					}
	  
	  
					TH1F *biaishisto_tmp = 0;
					biaishisto_tmp = (TH1F *)resolfile->Get(biaishistoname);
					if (biaishisto_tmp) {
						//fmapHistoBiais[KeyRun] = new TH1F(*biaishisto_tmp);
						biaishisto_tmp->SetDirectory(0);
						fmapHistoBiais[KeyRun] = biaishisto_tmp;
						if (biaishisto_tmp->GetEntries()==0) {
							std::cout << "\033[31;40;1m" << "fmapHistoBiais[" << KeyRun << "]->GetEntries() = " << biaishisto_tmp->GetEntries() 
									  << " --> THIS IS A BIG PROBLEM CHECK THE TABLE AND REPROCESS IT!!" << "\033[0m" << std::endl;
						}

					}
					else {
						fmapHistoBiais[KeyRun] = 0;
						WARN_OUT << "VIM DEBUG> I couldn't find the histo named = " << biaishistoname << std::endl;
					}

					if(fmapHistoBiais[KeyRun]) {
						std::vector<double> x,y; // will contain energy and biais
						MakeVectorsFromMap(fmapHistoBiais,KeyRun,x,y);
						fmapVectorsBiais[KeyRun] = std::make_pair(x,y);
						x.clear(); y.clear();
					}
					else {
						std::cout << "Can't find TH1F : " << biaishistoname << " ==> exit" << std::endl;
						return -1;
					}

					TH1F *resolhisto_tmp = 0;
					resolhisto_tmp = (TH1F *)resolfile->Get(resolhistoname);
					if (resolhisto_tmp) {
						//fmapHistoResol[KeyRun] = new TH1F(*resolhisto_tmp);
						resolhisto_tmp->SetDirectory(0);
						fmapHistoResol[KeyRun] = resolhisto_tmp;
						if (resolhisto_tmp->GetEntries()==0) {
							std::cout << "\033[31;40;1m" << "fmapHistoResol[" << KeyRun << "]->GetEntries() = " << resolhisto_tmp->GetEntries() 
									  << " --> THIS IS A BIG PROBLEM CHECK THE TABLE AND REPROCESS IT!!" << "\033[0m" << std::endl;
						}
					}
					else {
						fmapHistoResol[KeyRun] = 0;
						DEBUG_OUT << "VIM DEBUG> I couldn't find the histo named = " << resolhistoname << std::endl;
					}

					if(fmapHistoResol[KeyRun]) {
						std::vector<double> x,y; // will contain energy and biais
						MakeVectorsFromMap(fmapHistoResol,KeyRun,x,y);
						fmapVectorsResol[KeyRun] = std::make_pair(x,y);
						x.clear(); y.clear();
					}
					else {
						WARN_OUT << "Can't find TH1F : " << resolhistoname << " ==> exit" << std::endl;
						return -1;
					}
	  
					// We build the map's keys and the pairs or vectors associated :
					// map pair : sigma and mean of resolution
					// map : area
					// map : biais
					// map : resolution

					for(std::vector<double>::iterator ien=enMC.begin(); ien!=enMC.end(); ++ien) {
	    
						// Build the name of the EnergyDistribution file
	    
						TString distrhistoname = "EnergyResol_";
						std::ostringstream oss_sen;
						oss_sen.precision(3);
						oss_sen << std::fixed << *ien;
						TString sen(oss_sen.str().c_str());
						oss_sen.str("");
	    
						distrhistoname+=sen;
						distrhistoname+="_";
						distrhistoname+=szen;
						distrhistoname+="deg_";
						distrhistoname+=soff;
						distrhistoname+="off_";	    
						distrhistoname+="eff";	    
						distrhistoname+=seff;
	    
						TString KeyRunEnergy = MakeKeyRunEnergy(*ieff, *izen, *ioff, *ien, *itel); // VIM : Try to improve "readability"
	    
						// Build the map's pairs
						//if (IsDistrFileExists) {
	    
						if (fForceUseDistribution) { // VIM : Fill the map is not usefull if we won't use them. Will spare time.
							TH1F *hdistrib = 0;
							TH1F *hdistrib_tmp = 0;
							hdistrib_tmp = (TH1F *)distrfile->Get(distrhistoname);

							// if (hdistrib_tmp) {
							////hdistrib = new TH1F(*hdistrib_tmp);
							//	hdistrib_tmp->SetDirectory(0);
							//	hdistrib = hdistrib_tmp;
							//}
							if (hdistrib_tmp) {
								//hdistrib = new TH1F(*hdistrib_tmp);
								hdistrib_tmp->SetDirectory(0);
								hdistrib = hdistrib_tmp;
							}
							else {
								//std::cout << REDCOLOR << "WARNING> The histo : " << distrhistoname << " in " << distrfile->GetName()
								//	  << " is not found ! --> I CREATE A DUMMY ONE (VERY RISKY!) ! YOU SHOULD SOLVE THE TABLE PROBLEM INSTEAD !" << RESETCOLOR << std::endl;
								hdistrib = new TH1F(distrhistoname,distrhistoname,100,-3.,3);
							}
	      
							if (!hdistrib) {
								//std::cout << "WARNING> The histo : " << distrhistoname << " in " << distrfile->GetName() 
								//	  << " is not found ! --> That's not normal ! " << std::endl;
								exit(EXIT_FAILURE);
							}
							else {
								// VIM : There is some Sanity Check to implement !!!!! Integral to 1. Put Zero when it is Nan.
								for (int iix=1;iix<=hdistrib->GetNbinsX();++iix) {
									double bincont = (double)hdistrib->GetBinContent(iix);
									//std::cout << "distrhistoname = " << distrhistoname << " hdistrib->GetBinContent(" << iix << ") = " << bincont << std::endl;
									if (TMath::IsNaN(bincont)==1) {hdistrib->SetBinContent(iix,0.);} // VIM : remove any Nan
								}
		
								double total_hdist = hdistrib->Integral(1,hdistrib->GetNbinsX());
								//if (TMath::IsNaN(total_hdist)==1) {std::cout << "STRANGE WE ARE NOT SUPPOSED TO BE HERE" << std::endl;}
		
								if (total_hdist!=0.0) {
									hdistrib->Scale(1./total_hdist); // VIM : We ensure that the total is 1.
								}
								//std::cout << "distrhistoname = " << distrhistoname << " hdistrib->Integral(1," 
								//<< hdistrib->GetNbinsX() << ") = " << hdistrib->Integral(1,hdistrib->GetNbinsX()) << std::endl;
							}
							fmapHistoDistrib[KeyRunEnergy] = hdistrib;
						}
	  
						// Build the map of Resolution
					} // energy loop
				} //offset loop
			} // zenith loop

			delete areafile;
			delete resolfile;
			if (distrfile!=0) delete distrfile; 
			distrfile=0;

		} // efficiency loop  
	} // teltype loop
  
#if DEBUG
	std::cout << "Debug ending START::HandleResolArea::ReadAndStore()" << std::endl;
#endif
  
	return 0;
}


/**
 * \brief Build two vectors x and y whith the values contained in the histogram
 * stored in the map "Map" with the key "Key"
 */
void START::HandleResolArea::MakeVectorsFromMap(std::map<TString, TH1F *> Map, TString Key,
												std::vector<double> &x, std::vector<double> &y) const
{

	x.clear();
	y.clear();
    
	//MonteCarlo MC;
	std::vector<double> enMC = fMc->GetEnergy();
    
	x.reserve(enMC.size()); // VIM : Avoid allocation/reallocation does at each push_back (gain a small amount of time)
	y.reserve(enMC.size());
    
	std::map<TString,TH1F*>::iterator it;
	it = Map.find(Key); 
    
	// VIM : In theory, there is one bin for MCenergy, so I should scan with MCenergy.
	for (unsigned int i=0;i<enMC.size();++i) {
		int hbin = it->second->FindBin(TMath::Log10(enMC[i]));
		if (hbin!=0 && hbin!=(it->second->GetNbinsX()+1)) {
			double hval = it->second->GetBinContent(hbin);
			x.push_back(enMC[i]);
			y.push_back(hval);
		}
		else {
			std::cout << "VIM : We should not be here !!! This mean the histogram does not contain the MC value (This imply a problem at the construction) !! enMC[" << i << "] = " << enMC[i] << " and hbin = " << hbin  << std::endl;
		}
	}
}

/**
 * \brief Make a key  as a function of MC efficiencys, zenith angle and offset
 * for the local storage of weights in a std::map
 */
TString START::HandleResolArea::MakeKeyWeight(double eff, double zen, double off)
{

	std::ostringstream oss_eff;
	std::ostringstream oss_off;
	std::ostringstream oss_zen;
  
	oss_eff.precision(0);
	oss_zen.precision(0);
	oss_off.precision(1);
  
	oss_eff << std::fixed <<  eff;
	oss_zen << std::fixed <<  zen;
	oss_off << std::fixed <<  off;
  
	std::ostringstream oss;
	oss << "weight_" << oss_eff.str() << "_zen" << oss_zen.str() << "_off" << oss_off.str();
  
	TString key2(oss.str().c_str());
	DEBUG_OUT << "key = " << key2 << std::endl;
	return key2;
}

/**
 * \brief Make a key  as a function of MC efficiencys, zenith angle and offset and energy
 * for the local storage of weight for the energy distribution
 */
TString START::HandleResolArea::MakeKeyWeightEnergy(double eff, double zen, double off, double ene)
{

	std::ostringstream oss_eff;
	std::ostringstream oss_off;
	std::ostringstream oss_zen;
	std::ostringstream oss_mce;
  
	oss_eff.precision(0);
	oss_zen.precision(0);
	oss_off.precision(1);
	oss_mce.precision(3);
  
	oss_eff << std::fixed <<  eff;
	oss_zen << std::fixed <<  zen;
	oss_off << std::fixed <<  off;
	oss_mce << std::fixed <<  ene;
  
	std::ostringstream oss;
	oss << "weight_" << oss_eff.str() << "_zen" << oss_zen.str() << "_off" << oss_off.str() << "_en" << oss_mce.str();
  
	TString key2(oss.str().c_str());
	DEBUG_OUT << "key = " << key2 << std::endl;
	return key2;
}


/**
 * \brief Make a key  as a function of efficiency, zenith angle and offset
 * for the std::maps which will contain histograms
 */
TString START::HandleResolArea::MakeKeyRun(double eff, double zen, double off, TString tel)
{
  
	std::ostringstream oss_eff;
	std::ostringstream oss_off;
	std::ostringstream oss_zen;
  
	oss_eff.precision(0);
	oss_zen.precision(0);
	oss_off.precision(1);
  
	oss_eff << std::fixed <<  eff;
	oss_zen << std::fixed <<  zen;
	oss_off << std::fixed <<  off;
  
	std::ostringstream oss;
	oss << "histo_eff" << oss_eff.str() << "_zen" << oss_zen.str() << "_off" << oss_off.str() << "_" << tel.Data();
	// VIM : I would prefer to have the _ at the end of the cline only if there is something in tel.data()
  
	TString key2(oss.str().c_str());
	//DEBUG_OUT << "key = " << key2 << std::endl;
	return key2;  
}

/**
 * \brief Make a string for the map with efficiency, zenith, offset, energy and teltcode
 */
TString START::HandleResolArea::MakeKeyRunEnergy(double eff, double zen, double off, double en, TString tel)
{	
  
	std::ostringstream oss_eff;
	std::ostringstream oss_off;
	std::ostringstream oss_zen;
	std::ostringstream oss_mce;
  
	oss_eff.precision(0);
	oss_zen.precision(0);
	oss_off.precision(1);
	oss_mce.precision(3);
  
	oss_eff << std::fixed <<  eff;
	oss_zen << std::fixed <<  zen;
	oss_off << std::fixed <<  off;
	oss_mce << std::fixed <<  en;
  
	std::ostringstream oss;
	oss << "Values_eff" << oss_eff.str() << "_zen" << oss_zen.str() << "_off" << oss_off.str() << "_en" << oss_mce.str() << "_" << tel.Data();
	TString key2(oss.str().c_str());
  
	DEBUG_OUT << "key = " << key2 << std::endl;
	return key2;  
}

/**
 * \brief Make a string for the telescope with MC values of teltype
 */
TString START::HandleResolArea::MakeStringTelCode(int telcode) {
	TString Teltype = "";
	//MonteCarlo MC;
	std::vector<TString > telMC = fMc->GetTelType();
	if(telcode==30) {
		Teltype = telMC[0];
	}
	else if(telcode<30) {
		Teltype = telMC[1];
	}
	return Teltype;
}

/**
 * \brief Make a weighted energy distribution for a given zenith, offset and efficiency.
 * \warning VIM : Since I am completely stupid, I wrote variable named log_etrue_over_ereco but what is stored in the EnergyResolution Table IS in ln(Ereco/Etrue) !!!!!
 * \warning The algorithm is a bit complex. I think it works, but I would be happy to have a x-check (detail in the code itself)
 * \todo Check the alidity of the method used to create the distribution instrument map.
 */
void START::HandleResolArea::MakeWeightedInstrumentDistributionForBand(Band myBand,
																	   std::map<double,std::pair< std::vector<double>,std::vector<double> > > &DistributionInstrumentMapTable) {
  
	DistributionInstrumentMapTable.clear();
  
	// min and max MC values for the band
	double effminMC = FindEffMinMC(myBand);
	double effmaxMC = FindEffMaxMC(myBand);
	double offminMC = FindOffMinMC(myBand);
	double offmaxMC = FindOffMaxMC(myBand);
	double zenminMC = FindZenMinMC(myBand);
	double zenmaxMC = FindZenMaxMC(myBand);

	double zenminMCCos = TMath::Cos(TMath::DegToRad()*zenminMC);
	double zenmaxMCCos = TMath::Cos(TMath::DegToRad()*zenmaxMC);
  
#if DEBUG
	std::cout << "DEBUG VIM> INFO : " << std::endl;
	myBand.Print();
	std::cout << "effmin : " << effminMC << std::endl;
	std::cout << "eff    : " << myBand.GetEff() << std::endl;
	std::cout << "effmax : " << effmaxMC << std::endl;
	std::cout << "offmin : " << offminMC << std::endl;
	std::cout << "off    : " << myBand.GetOffset() << std::endl;
	std::cout << "offmax : " << offmaxMC << std::endl;
	std::cout << "zenmin : " << zenminMC << std::endl;
	std::cout << "zen    : " << myBand.GetZenON() << std::endl;
	std::cout << "zenmax : " << zenmaxMC << std::endl;
	std::cout << "DEBUG VIM> FIN INFO" << std::endl;
#endif
  
	// MC's Vectors used in the following loop
	//MonteCarlo MC;
	std::vector<double> effMC = fMc->GetEfficiency();
	std::vector<double> offMC = fMc->GetOffset();
	std::vector<double> zenMC = fMc->GetZenith();
	std::vector<double>  enMC = fMc->GetEnergy();
	//std::cout << "sizeMC" << enMC.size() << std::endl;
	//for(int itest=0; itest<enMC.size(); itest++){
	//std::cout << enMC[itest] << std::endl;
	//}
	// Parameters of the band
	double effBand = myBand.GetEff();
	double zenBand = myBand.GetZenON();
	double zenBandCos = TMath::Cos(TMath::DegToRad()*zenBand);
	double offBand = myBand.GetOffset();
	int telcodeBand = myBand.GetTelCode();
	TString teltype = MakeStringTelCode(telcodeBand);
  
	/* 
	   Loop to determine weights of area or resol determined for each
	   combination of (effMC, zenMC, offsetMC).
	   The 8=2x2x2 weights for this band are stored in the map mapWeights.
	   This can be 2x2 if there efficiency is >100% or less than 50% (the current smaller efficiency used in MC)
	*/
  
	std::map<TString,double> mapWeights;
  
	unsigned int n_ok_bands=8;
	if (effminMC == effmaxMC) {
		n_ok_bands=4;
	}
  
	DEBUG_OUT << "VIM DEBUG> n_ok_bands = " << n_ok_bands << std::endl;
  
	double sumweight=0.;
	for(std::vector<double>::iterator ieff=effMC.begin(); ieff!=effMC.end(); ++ieff) {
    
		if(*ieff!=effminMC && *ieff!=effmaxMC) continue;
    
		// determine the weight for efficiency
		double effweight(0);
		if(effMC.size()>1) { // JLK add
			if(*ieff==effminMC) {
				effweight=(effmaxMC-effBand)/(effmaxMC-effminMC);
			}
			else if(*ieff==effmaxMC) {
				effweight=(effBand-effminMC)/(effmaxMC-effminMC);
			}
		}
		else // no interpolation
			effweight = 1.;
    
		for(std::vector<double>::iterator izen=zenMC.begin(); izen!=zenMC.end(); ++izen) {
      
			if(*izen!=zenminMC && *izen!=zenmaxMC) continue;
      
			// determine the weight for zenith
			// Cosine function decreases with theta between 0 and 90 so order is reversed
			double zenweight(0);
			if(zenMC.size()>1) { // JLK add 
				if(*izen==zenminMC) {
					zenweight=(zenBandCos-zenmaxMCCos)/(zenminMCCos-zenmaxMCCos);
				}
				else if (*izen==zenmaxMC) {
					zenweight=(zenminMCCos-zenBandCos)/(zenminMCCos-zenmaxMCCos);
				}
			}
			else // no interpolation
				zenweight = 1.;
			for(std::vector<double>::iterator ioff=offMC.begin(); ioff!=offMC.end(); ++ioff) {

				if(*ioff!=offminMC && *ioff!=offmaxMC) continue;
	
				// determine the weight for offset
				double offweight(0);
				if(offMC.size()>1) { // JLK add
					if(*ioff==offminMC) {
						offweight=(offmaxMC-offBand)/(offmaxMC-offminMC);
					}
					else if(*ioff==offmaxMC) {
						offweight=(offBand-offminMC)/(offmaxMC-offminMC);
					}
				}
				else
					offweight = 1.;
				
				// Il faut que je fasse ici la boucle sur l'energie !!
				for (std::vector<double>::iterator ien = enMC.begin();ien!=enMC.end();++ien)  {
					//std::vector<double>  enMC = MC.GetEnergy();
					mapWeights[MakeKeyWeightEnergy(*ieff,*izen,*ioff,*ien)] = effweight*zenweight*offweight;
					sumweight+=(effweight*zenweight*offweight);
					DEBUG_OUT << "weight = " << effweight*zenweight*offweight  << " sumweight = " << sumweight << std::endl;
	  
				} // loop energy
			} // loop offset
		} // loop zenith angle
	} // loop efficacity

	DEBUG_OUT << "The sum of the weight is : " << sumweight << std::endl;
  
  
	TH1F *hresol_referencehistogram = new TH1F("ReferenceHistos","E resolution",150,-2.9,2.9); // VIM : histogram temporaire qui ne sert a rien, 
	//juste pour permettre d'avoir une reference pour les points a la fin pour calculer l'interpolateur finale.

	/*
	  We're going to apply here the weight on each instrument function.
	  For the moment, I only compute tables if there is the 8 value for the weight.
	  It will need to change a bit in the case where the efficiency is below/aboveabove the minimum/maximum 
	  simulated value respectively.
	  The principle is quite ugly for the moment, but I don't really see how I can do better.
	  In few words the principle of the code below is : 
    
	  1) For each good value of the simulation that contains the "true" observation, I make an interpolation of the histogram 
	  of the resolution distribution. 
	  2) I put in a new histogram the value of the wighted sum of each interpolated histogram at the bincenter position. 
	  This new historam could be then considered as the weighted distribution.
	  I can't add directly the histos because they are rebined if they don't have enough statistic in order to make the fit 
	  work for the gaussian approximation. 
	  Moreover, I add also a "safe" selection on which histogram bins I take into account to make the interpolation. 
	  For this, If there is a bin that is empty, but the next and previous bin is filled, then I don't take this one into account, 
	  because I assume that this is a statistical fluctuation. I think this won't change so much the results, but IT NEEDED TO BE TESTED.
    
	  Oui, je sais, c'est un peu une usine a gaz...
	  Peut etre que TH1F::Merge pourrait marcher, mais j'ai un doute....

	  \warning VIM : Since I am completely stupid, I wrote variable named log_etrue_over_ereco but what is stored in the EnergyResolution Table IS in ln(Ereco/Etrue) !!!!!
	*/
 
	//std::map<double, std::pair< std::vector<double> , std::vector<double> > > DistributionInstrumentMapTable;

	for (std::vector<double>::iterator ien = enMC.begin();ien!=enMC.end();++ien)  {
		int nhistoadded=0;
    
		// VIM : Will contain the value of the 
		//Distribution Histogram for a given energy, zenith angle, offset, efficiency, etc...
		std::map< TString , std::pair< std::vector<double> , std::vector<double> > > map_EDistribForATrueEnergy;     
		for(std::vector<double>::iterator ieff=effMC.begin(); ieff!=effMC.end(); ++ieff) {
			if(*ieff!=effminMC && *ieff!=effmaxMC) continue;
			for(std::vector<double>::iterator izen=zenMC.begin(); izen!=zenMC.end(); ++izen) {
				if(*izen!=zenminMC && *izen!=zenmaxMC) continue;
				for(std::vector<double>::iterator ioff=offMC.begin(); ioff!=offMC.end(); ++ioff) {
					if(*ioff!=offminMC && *ioff!=offmaxMC) continue;
	  
					TString KeyRunEnergy = MakeKeyRunEnergy(*ieff, *izen, *ioff, *ien, teltype); 
					std::map<TString,TH1F*>::iterator fmapHistoDistrib_elem;
					fmapHistoDistrib_elem = fmapHistoDistrib.find(KeyRunEnergy);
	  
					if (fmapHistoDistrib_elem==fmapHistoDistrib.end()) { 
						//std::cout << "VIM DEBUG> The map element is associated to the key " 
						//<< KeyRunEnergy << "has not been found... THIS IS TRULY PROBLEMATIC!!!" << std::endl;
						exit(EXIT_FAILURE);
					}
	  
					TH1F *histotmp = fmapHistoDistrib_elem->second; //fmapHistoDistrib[KeyRunEnergy];
	  
					if (!histotmp) {
						WARN_OUT << "VIM DEBUG> The histogram associated with the key : " << KeyRunEnergy << " Does not exist... THIS IS TRULY PROBLEMATIC !!!" << std::endl;
						exit(EXIT_FAILURE); // VIM : Oui c'est violent...!
						// return;
					}
	  
					std::map<TString,double>::iterator mapWeights_elem;
					TString KeyWeightEnergy = MakeKeyWeightEnergy(*ieff,*izen,*ioff,*ien);
					mapWeights_elem = mapWeights.find(KeyWeightEnergy);
	  
					if (mapWeights_elem==mapWeights.end()) {
						WARN_OUT << "VIM DEBUG> I can't find the mapweight element associated to this key : " << KeyWeightEnergy << " ... THIS IS TRULY PROBLEMATIC !!!" << std::endl;
						exit(EXIT_FAILURE);
						// return;
					}
	  
					std::vector<double> log_etrue_over_ereco;
					std::vector<double> value_of_log_etrue_over_ereco;
					log_etrue_over_ereco.reserve(histotmp->GetNbinsX());
					value_of_log_etrue_over_ereco.reserve(histotmp->GetNbinsX());
	  
					if (histotmp->Integral(1,histotmp->GetNbinsX())>0.0) { // There is something in the map
	    
						// VIM : I scan the histogram. If this current bin and the following and previous is zero, then I put it to zero, 
						//otherwise I skip it. This is in order to prevent hole, and to take into account the low statistics that can give this hole.

						for (int iih=1;iih<=histotmp->GetNbinsX();++iih) {
							double histo_bin_val = histotmp->GetBinContent(iih);
	      
							//std::cout << "VIMVIMVIMVIM> histotmp->GetBinCenter(" << iih << ") = " << histotmp->GetBinCenter(iih) << std::endl;
	      
							if (histo_bin_val==0.0) {
								if (iih==1 || iih==histotmp->GetNbinsX()) {
									log_etrue_over_ereco.push_back(histotmp->GetBinCenter(iih));
									value_of_log_etrue_over_ereco.push_back(histo_bin_val);
									//std::cout << "iih = " << iih << " log_etrue_over_ereco.push_back(histotmp->GetBinCenter(iih)) = " << (histotmp->GetBinCenter(iih)) << " histo_bin_val = " << histo_bin_val << std::endl;
								}
								else {
									if ( histotmp->GetBinContent(iih-1)>0.0 && histotmp->GetBinContent(iih+1)>0.0 ) {
										// VIM : Cela veut dire que l'on est dans un trou, surement par manque de stat, on ne rempli pas le vecteur..
										//std::cout << "iih = " << iih << std::endl;
									}
									else {
										log_etrue_over_ereco.push_back(histotmp->GetBinCenter(iih));
										value_of_log_etrue_over_ereco.push_back(histo_bin_val);
										//std::cout << "iih = " << iih << " log_etrue_over_ereco.push_back(histotmp->GetBinCenter(iih)) = " << (histotmp->GetBinCenter(iih)) << " histo_bin_val = " << histo_bin_val << std::endl;
									}
								}
							}
							else {
								// VIM : On a une entree, alors on rempli normalement.
								log_etrue_over_ereco.push_back(histotmp->GetBinCenter(iih));
								value_of_log_etrue_over_ereco.push_back(histo_bin_val);
								//std::cout << "iih = " << iih << " log_etrue_over_ereco.push_back(histotmp->GetBinCenter(iih)) = " << (histotmp->GetBinCenter(iih)) << " histo_bin_val = " << histo_bin_val << std::endl;
							}
						}
						++nhistoadded;
						map_EDistribForATrueEnergy[KeyWeightEnergy] = std::make_pair(log_etrue_over_ereco,value_of_log_etrue_over_ereco);
					}
					else {
						// There is no entries in the histogram, this is mean that the histogram can't be use....
	    
						// VIM : We still fill the vector, just in case...
						for (int iih=1;iih<=histotmp->GetNbinsX();++iih) {
							log_etrue_over_ereco.push_back(histotmp->GetBinCenter(iih));
							value_of_log_etrue_over_ereco.push_back(histotmp->GetBinContent(iih));
						}
						// VIM : But we don't fill the map_EDistribForATrueEnergy..
					}
	  
					//std::cout << "VIM GROSDEBUG> histotmp->GetNbinsX() = " << histotmp->GetNbinsX() << " log_etrue_over_ereco.size() = " 
					//<< log_etrue_over_ereco.size() << std::endl;
					//std::cout << "VIM DEBUG> nhistoadded = " << nhistoadded << std::endl;
	  

#if DEBUG
					std::cout << "VIM DEBUG> OUTPUT OF THE VECTOR OF THE DISTRIBUTION OF THE BAND FOR EFF = " << *ieff << " ZEN = " 
							  << *izen << " OFF= " <<  *ioff  << " ENE = " << *ien << " TEL = " << teltype << std::endl;
					for (int iiv=0;iiv<log_etrue_over_ereco.size();++iiv) {
						//std::cout << log_etrue_over_ereco[iiv] << " " << value_of_log_etrue_over_ereco[iiv] << std::endl;
					}
					std::cout << "VIM DEBUG> END OF THE OUTPUT" << std::endl;
#endif
	  
					if (log_etrue_over_ereco.size()==0) {
						std::cout << "VIM DEBUG> log_etrue_over_ereco.size() = " << log_etrue_over_ereco.size() 
								  << " C'EST ETRANGE QUE L'ON SOIT LA, C'EST MEME PAS NORMAL DU TOUT" << std::endl;
						exit(EXIT_FAILURE);
					}
	  
					//map_EDistribForATrueEnergy[KeyWeightEnergy] = std::make_pair(log_etrue_over_ereco,value_of_log_etrue_over_ereco);
	  
				} // loop offset
			} // loop zenith angle
		} // loop efficacity
    
    
		// VIM NOTES: At this point we have stored the vector which represent the bin content of each energy distribution 
		// for each MC value that surround the actual observation value.
		// VIM NOTES: In the next step, we're going to combine them in order to have the distribution expected at each MC energies.
    
		std::vector<double> log_etrue_over_ereco_resample;
		std::vector<double> value_of_log_etrue_over_ereco_resample;
		log_etrue_over_ereco_resample.reserve(hresol_referencehistogram->GetNbinsX());
		value_of_log_etrue_over_ereco_resample.reserve(hresol_referencehistogram->GetNbinsX());
    
		for (int iinv=1;iinv<=hresol_referencehistogram->GetNbinsX();++iinv) {
			log_etrue_over_ereco_resample.push_back(hresol_referencehistogram->GetBinCenter(iinv));
			value_of_log_etrue_over_ereco_resample.push_back(0.0);
		}
    
		double mysumweight = 0.0;
		DEBUG_OUT << "VIM DEBUG> map_EDistribForATrueEnergy.size() = " << map_EDistribForATrueEnergy.size() << " n_ok_bands = " << n_ok_bands << " GetNbRun() = " << myBand.GetNbRun() << std::endl;
    
		if (map_EDistribForATrueEnergy.size()==n_ok_bands) {
      
			for (std::map< TString , std::pair< std::vector<double> , std::vector<double> > >::iterator it = map_EDistribForATrueEnergy.begin(); it != map_EDistribForATrueEnergy.end();++it) {
				double weight_value = mapWeights[it->first];
				mysumweight+=weight_value;
				//std::cout << "it->first = " << it->first << " weight_value = " << weight_value << " sumweight = " << mysumweight << std::endl;
				//ROOT::Math::Interpolator InterpolDist((it->second).first,(it->second).second,ROOT::Math::Interpolation::kAKIMA);
				ROOT::Math::Interpolator InterpolDist((it->second).first,(it->second).second,ROOT::Math::Interpolation::kLINEAR);
				//std::cout << "((it->second).first).size() = " << ((it->second).first).size() << std::endl;
				double xmin = ((it->second).first)[0];
				int index_xmax = (int)((it->second).first).size()-1;
				//std::cout << "index_xmax = " << index_xmax << std::endl;
				double xmax = ((it->second).first)[index_xmax];
				//std::cout << "xmin = " << xmin << " xmax = " << xmax << std::endl;
				double interp_integral = InterpolDist.Integ(xmin,xmax);
	
				for (unsigned int iv=0;iv<log_etrue_over_ereco_resample.size();++iv) {
					double x_histo = log_etrue_over_ereco_resample[iv];
					if (x_histo>=xmin && x_histo<=xmax) {
						double interpol_value = weight_value * (InterpolDist.Eval(x_histo))/interp_integral;
						value_of_log_etrue_over_ereco_resample[iv] = value_of_log_etrue_over_ereco_resample[iv]+interpol_value;
					}
				}
			}
      
			ROOT::Math::Interpolator InterpolGen(log_etrue_over_ereco_resample,value_of_log_etrue_over_ereco_resample,ROOT::Math::Interpolation::kLINEAR);
			double integral_interpol_gen = InterpolGen.Integ(hresol_referencehistogram->GetBinCenter(1),hresol_referencehistogram->GetBinCenter(hresol_referencehistogram->GetNbinsX()));
			//std::cout << "VIM DEBUG> L'integral vaut " <<  integral_interpol_gen << std::endl;
      
      
#if DEBUG
			for (int iv=0;iv<log_etrue_over_ereco_resample.size();++iv) {
				hresol_referencehistogram->SetBinContent(hresol_referencehistogram->FindBin(log_etrue_over_ereco_resample[iv]),value_of_log_etrue_over_ereco_resample[iv]);
			}
      
			std::ostringstream oss_eff;
			oss_eff.precision(0);
			oss_eff << std::fixed << effBand;
			TString seff(oss_eff.str().c_str());
			oss_eff.str("");
			std::ostringstream oss_zen;
			oss_zen.fill('0');
			oss_zen.precision(0);
			oss_zen.width(2);
			oss_zen << std::fixed << zenBand;
			TString szen(oss_zen.str().c_str());
			oss_zen.str("");
			std::ostringstream oss_off;
			oss_off.precision(1);
			oss_off << std::fixed << offBand;
			TString soff(oss_off.str().c_str());
			oss_off.str("");
			std::ostringstream oss_sen;
			oss_sen.precision(3);
			oss_sen << std::fixed << *ien;
			TString sen(oss_sen.str().c_str());
			oss_sen.str("");
      
			TString distrhistoname = "EnergyResol_Interpolate_";
			distrhistoname+=sen;
			distrhistoname+="_";
			distrhistoname+=szen;
			distrhistoname+="deg_";
			distrhistoname+=soff;
			distrhistoname+="off_";           
			distrhistoname+="eff";            
			distrhistoname+=seff;
			distrhistoname+=".root";

			TFile *file = new TFile(distrhistoname.Data(),"RECREATE");
			file->cd();
			hresol_referencehistogram->Write();
			delete file;
			gROOT->cd();
			hresol_referencehistogram->Reset("M");

#endif

			for (unsigned int iv=0;iv<log_etrue_over_ereco_resample.size();++iv) {
				value_of_log_etrue_over_ereco_resample[iv] = value_of_log_etrue_over_ereco_resample[iv]/integral_interpol_gen;
			}
			//std::cout << "VIM DEBUG> L'integral vaut " 
			//<<  InterpolGen.Integ(hresol_referencehistogram->GetBinCenter(1),hresol_referencehistogram->GetBinCenter(hresol_referencehistogram->GetNbinsX())) 
			//<< std::endl;

			// JLK not used?
			/*   
				 ROOT::Math::Interpolator InterpolGen2(log_etrue_over_ereco_resample,value_of_log_etrue_over_ereco_resample,ROOT::Math::Interpolation::kLINEAR);
				 double integral_interpol_gen2 = 
				 InterpolGen2.Integ(hresol_referencehistogram->GetBinCenter(1),hresol_referencehistogram->GetBinCenter(hresol_referencehistogram->GetNbinsX()));
			*/
			//std::cout << "VIM DEBUG> L'integral vaut " <<  integral_interpol_gen2 << std::endl;

			std::pair< std::vector<double> , std::vector<double> > pair_for_energy_distrib;
			pair_for_energy_distrib.first = log_etrue_over_ereco_resample;
			pair_for_energy_distrib.second= value_of_log_etrue_over_ereco_resample;
			//std::cout << pair_for_energy_distrib.first.size() << std::endl;
			//std::cout << pair_for_energy_distrib.second.size() << std::endl;
			//std::cout <<" energy" << *ien << std::endl;
			DistributionInstrumentMapTable[*ien] = pair_for_energy_distrib;
			//std::cout <<" size map" << DistributionInstrumentMapTable.size() << std::endl;
			// VIM  : 
			/*
			  On a sauvearder les "histos" dans la structure Map_Distribution.
			  Il faut maintenant remplir la vrai dsitribution qui sera tocker dans Band : fDistributionInstrumentMapTable;
			  Une fois que cela sera gerer, il faudra alors creer la fonction pour calculer la resolution en energie, et son integral, 
			  comme ce qui est fait dans partialintegral.
			*/
		}
		else {
			DEBUG_OUT << "I can't handle the case where there is not 8 (or 4 depends of the value of muon efficiency) valid parameter for the weight" << std::endl; 
		}
    
    
	} // loop energy

	delete hresol_referencehistogram;  
	return;
}

/**
 * \brief Make a vector of resolution or effective area depending on the band
 *
 * The vector y is made with an interpolation in three dimension (eff,zen,offset).
 * Values are stored in vectors (x : energy and y : resol, biais or area)
 * Determine also the minimal MC energy for which the instrument values are complete 
 */
void START::HandleResolArea::MakeVectorForBand(Band myBand,
											   std::map<TString,std::pair<std::vector<double >,std::vector<double > > > mapVector,
											   std::vector<double> &x,std::vector<double> &y)
{
	x.clear();
	y.clear();

	// min and max MC values for the band
	double effminMC = FindEffMinMC(myBand);
	double effmaxMC = FindEffMaxMC(myBand);
	double offminMC = FindOffMinMC(myBand);
	double offmaxMC = FindOffMaxMC(myBand);
	double zenminMC = FindZenMinMC(myBand);
	double zenmaxMC = FindZenMaxMC(myBand);

	double zenminMCCos = TMath::Cos(TMath::DegToRad()*zenminMC);
	double zenmaxMCCos = TMath::Cos(TMath::DegToRad()*zenmaxMC);

#if DEBUG
	myBand.Print();
	std::cout << "effmin : " << effminMC << std::endl;
	std::cout << "effmax : " << effmaxMC << std::endl;
	std::cout << "offmin : " << offminMC << std::endl;
	std::cout << "offmax : " << offmaxMC << std::endl;
	std::cout << "zenmin : " << zenminMC << std::endl;
	std::cout << "zenmin : " << zenmaxMC << std::endl;
#endif

	// MC's Vectors used in the following loop
	//MonteCarlo MC;
	std::vector<double> effMC = fMc->GetEfficiency();
	std::vector<double> offMC = fMc->GetOffset();
	std::vector<double> zenMC = fMc->GetZenith();
	std::vector<double>  enMC = fMc->GetEnergy();

	// Parameters of the band
	double effBand = myBand.GetEff();
	double zenBand = myBand.GetZenON();
	double zenBandCos = TMath::Cos(TMath::DegToRad()*zenBand);
	double offBand = myBand.GetOffset();
	int telcodeBand = myBand.GetTelCode();
	TString teltype = MakeStringTelCode(telcodeBand);

	/* 
	   Loop to determine weights of area or resol determined for each
	   combination of (effMC, zenMC, offsetMC).
	   The 8=2x2x2 weights for this band are stored in the map mapWeights
	*/
	std::map<TString,double> mapWeights;

	for(std::vector<double>::iterator ieff=effMC.begin(); ieff!=effMC.end(); ++ieff) {

		if(*ieff!=effminMC && *ieff!=effmaxMC) continue;

		// determine the weight for efficiency
		double effweight(0);
		if(effMC.size()>1) { // JLK add
			if(*ieff==effminMC) {
				effweight=(effmaxMC-effBand)/(effmaxMC-effminMC);
			}
			else if(*ieff==effmaxMC) {
				effweight=(effBand-effminMC)/(effmaxMC-effminMC);
			}
		}
		else // no interpolation
			effweight = 1.;
		
		for(std::vector<double>::iterator izen=zenMC.begin(); izen!=zenMC.end(); ++izen) {

			if(*izen!=zenminMC && *izen!=zenmaxMC) continue;

			// determine the weight for zenith
			// Cosine function decreases with theta between 0 and 90 so order is reversed
			double zenweight(0);
			if(zenMC.size()>1) { // JLK add 
				if(*izen==zenminMC) {
					zenweight=(zenBandCos-zenmaxMCCos)/(zenminMCCos-zenmaxMCCos);
				}
				else if (*izen==zenmaxMC) {
					zenweight=(zenminMCCos-zenBandCos)/(zenminMCCos-zenmaxMCCos);
				}
			}
			else // no interpolation
				zenweight = 1.;

			for(std::vector<double>::iterator ioff=offMC.begin(); ioff!=offMC.end(); ++ioff) {

				if(*ioff!=offminMC && *ioff!=offmaxMC) continue;

				// determine the weight for offset
				double offweight(0);
				if(offMC.size()>1) { // JLK add
					if(*ioff==offminMC) {
						offweight=(offmaxMC-offBand)/(offmaxMC-offminMC);
					}
					else if(*ioff==offmaxMC) {
						offweight=(offBand-offminMC)/(offmaxMC-offminMC);
					}
				}
				else
					offweight = 1.;

				mapWeights[MakeKeyWeight(*ieff, *izen, *ioff)] = effweight*zenweight*offweight;
				// 	std::cout <<"key="<<MakeKeyWeight(*ieff, *izen, *ioff)
				// 		  <<" Wzen="<<zenweight
				// 		  <<" W="<<effweight*zenweight*offweight<<std::endl;
			} // loop offset

		} // loop zenith angle

	} // loop efficacity


	/* test */
	double sum(0);  
	for (std::map<TString,double>::const_iterator wit = mapWeights.begin() ; 
		 wit!=mapWeights.end() ; ++wit) {

		//std::cout <<"mapWeights["<<wit->first<<"]="<<wit->second<<std::endl;
		sum+=wit->second;
	}
	//std::cout<<"sum="<<sum<< std::endl;


	/*
	  First loop on MC energy bins: 
	  - find the first bin (# k) for which the 8=2x2x2 values (of area, resol or bias) are known
	  - for bin k-1 (num_ebin_tofix below): it could be useful to determine missing values using 
	  extrapolations (over energy).
	  Note that this loop is stopped before the end (see break)
	*/

	double FirstEmcVal(-1);
	int FirstEmcBinNumber(-1);//first MC energy bin with valid instrument response
	int num_ebin_tofix(-1);
	int n_ok_previous(0);

	if(DEBUG) {
		for(int ien(0); ien<(int)enMC.size(); ++ien) {
			std::cout <<"enMC["<<ien<<"]="<<enMC[ien]<<std::endl;
		}
	}

	for(int ien(0); ien<(int)enMC.size(); ++ien) {
		//std::cout <<"enMC["<<ien<<"]="<<enMC[ien]<<std::endl;
    
		int n_ok(0);
    
		for(std::vector<double>::iterator ieff=effMC.begin(); ieff!=effMC.end(); ++ieff) {
			if(*ieff!=effminMC && *ieff!=effmaxMC) continue;
      
			for(std::vector<double>::iterator izen=zenMC.begin(); izen!=zenMC.end(); ++izen) {
				if(*izen!=zenminMC && *izen!=zenmaxMC) continue;
	
				for(std::vector<double>::iterator ioff=offMC.begin(); ioff!=offMC.end(); ++ioff) {
					if(*ioff!=offminMC && *ioff!=offmaxMC) continue;
					if(DEBUG) {
						std::cout << "effminMC" << effminMC << " effmaxMC" << effmaxMC << std::endl;
						std::cout << "zenminMC" << zenminMC << " zenmaxMC" << zenmaxMC << std::endl;
						std::cout << "offminMC" << offminMC << " offmaxMC" << offmaxMC << std::endl;

						std::cout <<"key: "<<MakeKeyRun(*ieff,*izen,*ioff,teltype)<< std::endl;
					}
					std::vector<double> energymap = mapVector[MakeKeyRun(*ieff,*izen,*ioff,teltype)].first;
					std::vector<double> resultmap = mapVector[MakeKeyRun(*ieff,*izen,*ioff,teltype)].second;

					if(DEBUG) {
						for(int jen(0); jen<(int)energymap.size(); ++jen) {
							std::cout <<"energymap["<<jen<<"]="<<energymap[jen]<<std::endl;
						}
					}

					//consistency: check that always enMC[ien]==energymap[ien]
					if (enMC[ien]!=energymap[ien]) {
						/*
						  std::cout << "energy size=" << energymap.size() << std::endl;

						  for(int i(0); i<resultmap.size(); i++) {
						  std::cout << "energy=" << energymap[i] << std::endl;
						  std::cout << "result=" << resultmap[i] << std::endl;
						  }
						*/
						std::cout<<"FATAL ERROR: enMC["<<ien<<"]="<<enMC[ien]
								 <<"!=energymap["<<ien<<"]="<<energymap[ien]<<std::endl;
						std::cout << "key : " << MakeKeyRun(*ieff,*izen,*ioff,teltype) << std::endl;
						exit(EXIT_FAILURE);
					}
	  
					//std::cout<<"resultmap["<<ien<<"]="<<resultmap[ien]<<std::endl;
					if (resultmap[ien]!=0) n_ok++; 
				}
			}
		}

		//std::cout <<"For enMC="<<enMC[ien]<<"  n_ok="<<n_ok<<std::endl;
		// JLK change, depending on which parameter the interpolations are done
		// the number of weights will change (x2 for each interpolation)
		int nmaxvalues(8);
		if(effMC.size()==1)
			nmaxvalues/=2;
		if(offMC.size()==1)
			nmaxvalues/=2;
		if(zenMC.size()==1)
			nmaxvalues/=2;
		
		//if (n_ok==8) {
		if (n_ok==nmaxvalues) {
			FirstEmcBinNumber=ien; //eventually will be set to (ien-1) below
			FirstEmcVal=enMC[ien];
			num_ebin_tofix=ien-1;
			break;
		}

		n_ok_previous=n_ok;
	}

	//std::cout <<"num_ebin_tofix="<<num_ebin_tofix<<"  available values: "<<n_ok_previous<<"/8"<<std::endl;

	/*
	  For energy bin num_ebin_tofix:
	  fixing the missing instrument values (area, resol or bias) (TO BE DONE)
	*/
	//std::cout<<"FirstEmcBinNumber="<<FirstEmcBinNumber<<" => energy="<<FirstEmcVal<<std::endl;
	bool FixDone=false;
	if (n_ok_previous>4) {
		for(std::vector<double>::iterator ieff=effMC.begin(); ieff!=effMC.end(); ++ieff) {
			if(*ieff!=effminMC && *ieff!=effmaxMC) continue;
      
			for(std::vector<double>::iterator izen=zenMC.begin(); izen!=zenMC.end(); ++izen) {
				if(*izen!=zenminMC && *izen!=zenmaxMC) continue;
	
				for(std::vector<double>::iterator ioff=offMC.begin(); ioff!=offMC.end(); ++ioff) {
					if(*ioff!=offminMC && *ioff!=offmaxMC) continue;
	  
					std::vector<double> resultmap = mapVector[MakeKeyRun(*ieff,*izen,*ioff,teltype)].second;
					/*
					  if (resultmap[num_ebin_tofix]==0)
					  std::cout<<"TBD here: Fixing resultmap["<<num_ebin_tofix<<"]="<<resultmap[num_ebin_tofix]
					  <<" using resultmap["<<num_ebin_tofix+1<<"]="<<resultmap[num_ebin_tofix+1]
					  <<" and resultmap["<<num_ebin_tofix+2<<"]="<<resultmap[num_ebin_tofix+2]
					  <<std::endl;
					*/	 
				}
			}
		}
		//set FixDone to true if ok
		if (FixDone) {
			FirstEmcBinNumber=num_ebin_tofix;
			FirstEmcVal=enMC[num_ebin_tofix];
			//std::cout<<"FirstEmcBinNumber="<<FirstEmcBinNumber<<" => energy="<<FirstEmcVal<<std::endl;
		}
	}


	/*
	  Determine interpolated instrument functions at MC energies
	*/  
	std::vector<double> result; // result is area or resolution
	std::vector<double> sumweight;
	for(int ie(0); ie< (int)enMC.size(); ie++) {
		result.push_back(0); // will contain area or resol for each MC energies
		sumweight.push_back(0); // will contain sumweight for each MC energies
	}
	for(std::vector<double>::iterator ieff=effMC.begin(); ieff!=effMC.end(); ++ieff) {
		if(*ieff!=effminMC && *ieff!=effmaxMC) continue;
    
		for(std::vector<double>::iterator izen=zenMC.begin(); izen!=zenMC.end(); ++izen) {
			if(*izen!=zenminMC && *izen!=zenmaxMC) continue;
      
			for(std::vector<double>::iterator ioff=offMC.begin(); ioff!=offMC.end(); ++ioff) {
				if(*ioff!=offminMC && *ioff!=offmaxMC) continue;
	
				std::vector<double> energymap = mapVector[MakeKeyRun(*ieff,*izen,*ioff,teltype)].first;
				std::vector<double> resultmap = mapVector[MakeKeyRun(*ieff,*izen,*ioff,teltype)].second;
				double branchweight = mapWeights[MakeKeyWeight(*ieff, *izen, *ioff)];

				for(int ie(FirstEmcBinNumber); ie<(int)result.size(); ie++) {
					//	  if (resultmap[ie]==0) {
					//std::cout<<"FATAL ERROR: resultmap["<<ie<<"]=0... this should not happen"<<std::endl;
					//   std::cout << "key : " << MakeKeyRun(*ieff,*izen,*ioff,teltype) << std::endl;
					//   exit(EXIT_FAILURE);
					// }
					result[ie]+=resultmap[ie]*branchweight;
					sumweight[ie]+=branchweight;
					// 	  std::cout <<"Key="<<MakeKeyWeight(*ieff, *izen, *ioff)
					// 		    <<"  energymap["<<ie<<"]="<<energymap[ie]
					// 		    <<"  resultmap["<<ie<<"]="<<resultmap[ie]
					// 		    <<"  branchweight="<<branchweight<<std::endl;
				}
			}
		}
	}

	// for(int ie(FirstEmcBinNumber); ie<(int)enMC.size(); ie++) {
	//   if (fabs(sumweight[ie]-1)>1e-3)  {
	//     std::cout<<"FATAL ERROR: sumweight["<<ie<<"]="<<sumweight[ie]<<"!=1 ... this should not happen"<<std::endl;
	//     exit(EXIT_FAILURE);
	//   }
	// }
  
	//   for(int ie(FirstEmcBinNumber); ie<result.size(); ie++) {
	//     std::cout <<"  enMC["<<ie<<"]="<<enMC[ie]
	// 	      <<"  result["<<ie<<"]="<<result[ie]<<std::endl;
	//   }


	x = enMC;
	y = result;

}


/**
 * \brief Make Interpolated function for area, resolution and biais
 *
 * Vectors xinter and yinter will contain energy and {resolution or area or biais}.
 * Those vectors will be stocked in object band and will be used (or not!) in the LikeLihood class. 
 * The values are interpolated with GSL
 */
void START::HandleResolArea::MakeSmoothFunction(const std::vector<double> x, const std::vector<double> y,
												std::vector<double> &xinter, std::vector<double> &yinter)
{
	xinter.clear();
	yinter.clear();


	//interpol
	// CSPLINE, LINEAR, POLYNOMIAL,
	// CSPLINE_PERIODIC, AKIMA, AKIMA_PERIODIC
	ROOT::Math::Interpolator Interpol(x,y,ROOT::Math::Interpolation::kAKIMA);
	//MonteCarlo MC;
	double emin(fMc->GetEnergy().front()), emax(fMc->GetEnergy().back());
	double eminirf(x.front()), emaxirf(x.back());
	double log10emin = TMath::Log10(emin);
	double log10emax = TMath::Log10(emax);
	double npoint(50);
	double binsize = (log10emax-log10emin)/npoint;

	DEBUG_OUT << "Original points :" << std::endl;
	for(unsigned int ipoint(0); ipoint<x.size(); ipoint++)
		DEBUG_OUT << "en=" << x[ipoint] << " value=" << y[ipoint] << std::endl;

	DEBUG_OUT << "emin = " << emin << " emax = " << emax << std::endl; //JLK
	// values for energy in TeV
	for(double ien(log10emin); ien<=log10emax*1.001; ien+=binsize) {
		double energy = TMath::Power(10.,ien);
		double result_inter;
		if(energy>=eminirf && energy<=emaxirf) result_inter=Interpol.Eval(energy);
		else {
			unsigned int lastx(0);
			if(energy>emaxirf) lastx = x.size()-1;
			else lastx = 1;
			result_inter = y[lastx]+(energy-x[lastx])/(x[lastx]-x[lastx-1])*(y[lastx]-y[lastx-1]);
		}
		xinter.push_back(energy);
		yinter.push_back(result_inter);
	}
  
	DEBUG_OUT << "Interpolated/Extrapolated points :" << std::endl;
	for(unsigned int ipoint(0); ipoint<xinter.size(); ipoint++)
		DEBUG_OUT << "en=" << xinter[ipoint] << " value=" << yinter[ipoint] << std::endl;

}

/**
 * \brief Copy interpolators in bands
 */
int START::HandleResolArea::CopyInterpolatorsInBands()
{
  
	for (std::vector<Band>::iterator iband=fBandArray->begin(); iband!=fBandArray->end();++iband) {
		iband->SetGSLInterpolatorForArea(iband->GetVectorEnergy(),iband->GetVectorArea());
		iband->SetGSLInterpolatorForBiais(iband->GetVectorEnergy(),iband->GetVectorBiais());
		iband->SetGSLInterpolatorForResolution(iband->GetVectorEnergy(),iband->GetVectorResolution());
		if (fForceUseDistribution) {
			iband->InitDistributionInterpTable();
		}
	}

	return 1;
}

//surcharge de la fonction precedente ou en entree on donne l objetband
int START::HandleResolArea::CopyInterpolatorsInBands(std::vector<Band> &InputBandArray)
{
  
	for (std::vector<Band>::iterator iband=InputBandArray.begin(); iband!=InputBandArray.end();++iband) {
		iband->SetGSLInterpolatorForArea(iband->GetVectorEnergy(),iband->GetVectorArea());
		iband->SetGSLInterpolatorForBiais(iband->GetVectorEnergy(),iband->GetVectorBiais());
		iband->SetGSLInterpolatorForResolution(iband->GetVectorEnergy(),iband->GetVectorResolution());
		if (fForceUseDistribution) {
			iband->InitDistributionInterpTable();
		}
	}

	return 1;
}
/**
 * \brief Copy interpolators in bins
 */
int START::HandleResolArea::CopyInterpolatorsInBins()
{
    
	for (std::vector<Band>::iterator iband=fBandArray->begin(); iband!=fBandArray->end();++iband) {
		for (std::vector<EnergyBin>::iterator iener = (iband->ebin).begin(); iener!= (iband->ebin).end(); ++iener) {
			iener->SetGSLInterpolatorForPartialIntegral(iener->GetPartialIntegral());
		}
	}

	return 1;
}
//surcharge de la fonction precedente ou en entree on donne l objetband
int START::HandleResolArea::CopyInterpolatorsInBins(std::vector<Band> &InputBandArray)
{
    
	for (std::vector<Band>::iterator iband=InputBandArray.begin(); iband!=InputBandArray.end();++iband) {
		for (std::vector<EnergyBin>::iterator iener = (iband->ebin).begin(); iener!= (iband->ebin).end(); ++iener) {
			iener->SetGSLInterpolatorForPartialIntegral(iener->GetPartialIntegral());
		}
	}

	return 1;
}
/**
 * Copy vectors of resolution, energy and area in object bands
 */
int START::HandleResolArea::CopyVectorsInBands()
{
	// VIM : If I understand the call to the functiion MakeVectorForBand, this mean that the first argument is a copy of 
	// the BandArray... It seems that we could improve time by putting a reference to it !!
	for(int iband(0); iband<(int)(*fBandArray).size(); iband++) {
  
		std::vector<double> energy, smooth_energy;
		std::vector<double> area, smooth_area, resol ,smooth_resol, smooth_biais, biais;
    
		// TOTOR LE WARRIOR
   
		if((*fBandArray)[iband].GetKeepBand()==0) continue;

		// Copy area in bands
		MakeVectorForBand((*fBandArray)[iband],fmapVectorsArea,energy,area);
		(*fBandArray)[iband].SetVectorEnergy(energy);
		(*fBandArray)[iband].SetVectorArea(area);


		MakeSmoothFunction(energy,area,smooth_energy,smooth_area);
		(*fBandArray)[iband].SetVectorInterEnergy(smooth_energy);
		(*fBandArray)[iband].SetVectorInterArea(smooth_area);

		// Copy resolution in bands
		MakeVectorForBand((*fBandArray)[iband],fmapVectorsResol,energy,resol);
		(*fBandArray)[iband].SetVectorResolution(resol);
		MakeSmoothFunction(energy,resol,smooth_energy,smooth_resol);
		(*fBandArray)[iband].SetVectorInterResolution(smooth_resol);
		// Copy biais in bands
		//JFK : look in here :
		MakeVectorForBand((*fBandArray)[iband],fmapVectorsBiais,energy,biais);
		(*fBandArray)[iband].SetVectorBiais(biais);
		MakeSmoothFunction(energy,biais,smooth_energy,smooth_biais);
		(*fBandArray)[iband].SetVectorInterBiais(smooth_biais);

		area.clear();
		resol.clear();
		biais.clear();
		energy.clear();
		smooth_area.clear();
		smooth_resol.clear();
		smooth_biais.clear();
		smooth_energy.clear();

		// VIM : Pour les resultions : 
		std::map<double,std::pair< std::vector<double>,std::vector<double> > > DistributionInstrumentMapTable;
		if (fForceUseDistribution) {
			MakeWeightedInstrumentDistributionForBand((*fBandArray)[iband],DistributionInstrumentMapTable);
			(*fBandArray)[iband].SetDistributionVectorTable(DistributionInstrumentMapTable);
		}
	}
	return 0;
}

/**
 * Determine value of the first MC energy bin with valid area, resol AND bias
 */
void START::HandleResolArea::SetInBandsFirstEmcBin() {

	for(int iband(0); iband<(int)(*fBandArray).size(); iband++) {

		if((*fBandArray)[iband].GetKeepBand()==0) continue;

		std::vector<double> energy = (*fBandArray)[iband].GetVectorEnergy();
		std::vector<double> resol  = (*fBandArray)[iband].GetVectorResolution();
		std::vector<double> area   = (*fBandArray)[iband].GetVectorArea();
		std::vector<double> biais  = (*fBandArray)[iband].GetVectorBiais();

		/*
		  Note that Band::fVectorArea, Band::fVectorResolution and Band::fVectorBiais
		  have been filled in MakeVectorForBand so as to be 0 if no correct information 
		  is available.
		*/

		int binnum(-1);
		for (int ie(0); ie<(int)energy.size() ; ie++) {
			double prod = TMath::Abs(area[ie]*resol[ie]*biais[ie]);
			//std::cout<<ie<<" "<<area[ie]<<" "<<resol[ie]<<" "<<biais[ie]<<" "<<prod<<std::endl;
			if (prod>1e-3) {
				binnum=ie;
				break;
			}
		}

		(*fBandArray)[iband].SetFirstEmcBinNumber(binnum);
		//(*fBandArray)[iband].SetFirstEmcVal(energy[binnum]);
		(*fBandArray)[iband].SetFirstEmcValFit(energy[binnum]); // VIM
		/*std::cout <<"energy["<<binnum<<"]="<<energy[binnum]
		  <<"  (*fBandArray)["<<iband<<"].GetFirstEmcVal()="
		  <<(*fBandArray)[iband].GetFirstEmcVal() << std::endl;*/
		energy.clear(); resol.clear();  area.clear();  biais.clear();  
    
		if (fForceUseDistribution) {
			// std::cout << "(((*fBandArray)[iband].GetDistributionVectorTable()).begin())->first = " 
			//<< ((((*fBandArray)[iband]).GetDistributionVectorTable()).begin())->first << std::endl;
			((*fBandArray)[iband]).SetFirstEmcValDistrib(((((*fBandArray)[iband]).GetDistributionVectorTable()).begin())->first);
		}
		//std::cout << "VIM DEBUG> GetFirstEmcValFit() = " << (*fBandArray)[iband].GetFirstEmcValFit() 
		//<< " GetFirstEmcValDistrib() = " << (*fBandArray)[iband].GetFirstEmcValDistrib() << std::endl;
	}
}

void START::HandleResolArea::SetInBandsFirstEmcBin(std::vector<Band> &InputBandArray) {

	for(int iband(0); iband< InputBandArray.size(); iband++) {

		if(InputBandArray[iband].GetKeepBand()==0) continue;

		std::vector<double> energy = InputBandArray[iband].GetVectorEnergy();
		std::vector<double> resol  = InputBandArray[iband].GetVectorResolution();
		std::vector<double> area   = InputBandArray[iband].GetVectorArea();
		std::vector<double> biais  = InputBandArray[iband].GetVectorBiais();

		/*
		  Note that Band::fVectorArea, Band::fVectorResolution and Band::fVectorBiais
		  have been filled in MakeVectorForBand so as to be 0 if no correct information 
		  is available.
		*/

		int binnum(-1);
		for (int ie(0); ie<(int)energy.size() ; ie++) {
			double prod = TMath::Abs(area[ie]*resol[ie]*biais[ie]);
			//std::cout<<ie<<" "<<area[ie]<<" "<<resol[ie]<<" "<<biais[ie]<<" "<<prod<<std::endl;
			if (prod>1e-3) {
				binnum=ie;
				break;
			}
		}

		InputBandArray[iband].SetFirstEmcBinNumber(binnum);
		//InputBandArray[iband].SetFirstEmcVal(energy[binnum]);
		InputBandArray[iband].SetFirstEmcValFit(energy[binnum]); // VIM
		/*std::cout <<"energy["<<binnum<<"]="<<energy[binnum]
		  <<"  InputBandArray["<<iband<<"].GetFirstEmcVal()="
		  <<InputBandArray[iband].GetFirstEmcVal() << std::endl;*/
		energy.clear(); resol.clear();  area.clear();  biais.clear();  
    
		if (fForceUseDistribution) {
			// std::cout << "((InputBandArray[iband].GetDistributionVectorTable()).begin())->first = " 
			//<< (((InputBandArray[iband]).GetDistributionVectorTable()).begin())->first << std::endl;
			(InputBandArray[iband]).SetFirstEmcValDistrib((((InputBandArray[iband]).GetDistributionVectorTable()).begin())->first);
		}
		//std::cout << "VIM DEBUG> GetFirstEmcValFit() = " << InputBandArray[iband].GetFirstEmcValFit() 
		//<< " GetFirstEmcValDistrib() = " << InputBandArray[iband].GetFirstEmcValDistrib() << std::endl;
	}
}


/**
 * \brief Determine the value of 0 for a function by a dichotomie function.
 *
 * The value of x which correspond to f(x)=0 is returned.
 * The zeroexist variables return the number of time the function cross 0 BUT the value is compute for the first crossing region !!
 * The yatzero variables return the value for which we have
 */
double START::HandleResolArea::FindZeroInBandArray(std::vector<double> x, std::vector<double> y, int &zeroexist) { 
  
	zeroexist = 1;
  
	// VIM : We begin to scan if there is multiple solutions
	double ym = y[0];
	double yp;
	int ncrosszero = 0;
	int firstbincrosszero=0;
	for (unsigned int i=1;i<x.size();++i) {
		yp=y[i];
		if ( (yp>0. && ym<=0.) || (yp<=0. && ym>0.) ) { // VIM : Pour gerer les fonctions croissantes et decroissantes
			// VIM : I should have put <= somewhere, but in order to avoid problem with undefined values
			// in the resolution and biais (if there is no event, the value is put at 0! ). 
			// But what is the probability that we have EXACTELY 0 at one fixed energy ?
			ncrosszero++;
			if (ncrosszero==1) {
				firstbincrosszero = i;
			}
		}
	}
  
	if (ncrosszero==0) {
		zeroexist = 0;
		return x[0];
	}
  
	zeroexist = ncrosszero;
  
	ROOT::Math::Interpolator InterpolY(x,y,ROOT::Math::Interpolation::kAKIMA);
	// VIM : This mean that it cross zero between the bin firstbincrosszero and firstbincrosszero-1
	// We'll going to use a dichotomy to find the value. (not the most powerfull, but easy to implement)
  
	double prec = 1.e-3; // VIM : Too much ? 
	double ymin = y[firstbincrosszero-1];
	double ymax = y[firstbincrosszero-1];
	double ymid;
	double xmin = x[firstbincrosszero-1];
	double xmax = x[firstbincrosszero];
	double xmid = 0.5*(xmin+xmax);
	double deltax=1.e6; //enermax-enermin;
  
	int niter = 0;
	while (deltax>prec) {
		niter++;
		xmid = 0.5*(xmax+xmin);
		ymid = InterpolY.Eval(xmid);
		double sign = ymin*ymid;
    
		if (sign>0) {
			// VIM : This mean that areamid and areamin have the same sign, so mid become min !
			xmin = xmid;
			ymin = ymid;
		}
		else if (sign<0) {
			xmax = xmid;
			ymax = ymid;
		}
		else {
			std::cout << "VIM : You are very lucky, you have found the right value ... Strange because it's nearly impossible!!" << std::endl;
			ymin = ymid;
			xmin = xmid;
			ymax = ymid;
			xmax = xmid;
		}
		deltax = xmax-xmin;
	}
  
	DEBUG_OUT << "NUMBER OF ITERATION = " << niter << std::endl;
  
	return xmid;
}


/**
 * \brief Copy Effective area in each bins
 */
void START::HandleResolArea::SetInBinsEffectiveArea() {

	for(std::vector<Band>::iterator band=fBandArray->begin(); band!=fBandArray->end(); ++band) {

		if(band->GetKeepBand()==0) continue; // only intrested in selected bands

		for(std::vector<EnergyBin>::iterator bin=band->ebin.begin(); bin!=band->ebin.end(); ++bin) {
      
			bin->SetAcceff(band->GetInterpolatedArea(bin->GetEmid())*1.e-4);

		}

	}

}

//surcharge de la fonction precedente
void START::HandleResolArea::SetInBinsEffectiveArea(std::vector<Band> &InputBandArray) {

	for(std::vector<Band>::iterator band=InputBandArray.begin(); band!=InputBandArray.end(); ++band) {

		if(band->GetKeepBand()==0) continue; // only intrested in selected bands

		for(std::vector<EnergyBin>::iterator bin=band->ebin.begin(); bin!=band->ebin.end(); ++bin) {
      
			bin->SetAcceff(band->GetInterpolatedArea(bin->GetEmid())*1.e-4);

		}

	}

}
/**
 * \brief Set the flag in bins. If keepbin=1 the bin is kept and if not
 * keepbin=0
 *
 */
void START::HandleResolArea::SetInBinsKeepBin()
{

	for(std::vector<Band>::iterator band=fBandArray->begin(); band!=fBandArray->end(); ++band) {

		if(band->GetKeepBand()==0) continue; // only interested in selected bands

		for(std::vector<EnergyBin>::iterator bin=band->ebin.begin(); bin!=band->ebin.end(); ++bin) {
      
			switch(fEnergyBinThresholdCondition) {
			case Safe:
				if(bin->GetEmin()>=band->GetEthMC() && bin->GetKeepBin()==1) bin->SetKeepBin(1);
				else bin->SetKeepBin(0);
				break;
			case Explorer:
				if(bin->GetEmin()>=band->GetEthMC() && bin->GetKeepBin()==1) bin->SetKeepBin(1);
				else if(band->GetEthMC()>bin->GetEmin() && band->GetEthMC()<bin->GetEmax() && 
						(TMath::Abs(bin->GetEmin()-band->GetEthMC())<TMath::Abs(bin->GetEmax()-band->GetEthMC())) && 
						bin->GetKeepBin()==1) bin->SetKeepBin(1);
				else bin->SetKeepBin(0);
				/*
				  if(bin->GetEmax()>=band->GetEthMC() && bin->GetKeepBin()==1) bin->SetKeepBin(1);
				  else if((bin->GetEmin()>=band->GetEthMC() && bin->GetEmax()<=band->GetEthMC()) && 
				  (TMath::Abs(bin->GetEmin()-band->GetEthMC())<TMath::Abs(bin->GetEmax()-band->GetEthMC())) && 
				  bin->GetKeepBin()==1) bin->SetKeepBin(1);
				*/
				break;
			case Insane:
				bin->SetKeepBin(1);
			}
	
		}

	}

}

void START::HandleResolArea::SetInBinsKeepBin(std::vector<Band> &InputBandArray)
{

	for(std::vector<Band>::iterator band=InputBandArray.begin(); band!=InputBandArray.end(); ++band) {

		if(band->GetKeepBand()==0) continue; // only interested in selected bands

		for(std::vector<EnergyBin>::iterator bin=band->ebin.begin(); bin!=band->ebin.end(); ++bin) {
      
			switch(fEnergyBinThresholdCondition) {
			case Safe:
				if(bin->GetEmin()>=band->GetEthMC() && bin->GetKeepBin()==1) bin->SetKeepBin(1);
				else bin->SetKeepBin(0);
				break;
			case Explorer:
				if(bin->GetEmin()>=band->GetEthMC() && bin->GetKeepBin()==1) bin->SetKeepBin(1);
				else if(band->GetEthMC()>bin->GetEmin() && band->GetEthMC()<bin->GetEmax() && 
						(TMath::Abs(bin->GetEmin()-band->GetEthMC())<TMath::Abs(bin->GetEmax()-band->GetEthMC())) && 
						bin->GetKeepBin()==1) bin->SetKeepBin(1);
				else bin->SetKeepBin(0);
				/*
				  if(bin->GetEmax()>=band->GetEthMC() && bin->GetKeepBin()==1) bin->SetKeepBin(1);
				  else if((bin->GetEmin()>=band->GetEthMC() && bin->GetEmax()<=band->GetEthMC()) && 
				  (TMath::Abs(bin->GetEmin()-band->GetEthMC())<TMath::Abs(bin->GetEmax()-band->GetEthMC())) && 
				  bin->GetKeepBin()==1) bin->SetKeepBin(1);
				*/
				break;
			case Insane:
				bin->SetKeepBin(1);
			}
	
		}

	}

}


/**
 * \brief Function that setup the condition and method to compute the sasfe threshold for the analysis
 *
 * \param emeth enum of ESafeThresholdMethod type that will help to select the method use to compute 
 * the threshold in SetInBandsMCEnergyThreshold()
 * \param areamin_or_areacondition parameter to give the areamin, or the fraction of the maximal area 
 * (according to the method you choose)
 * \param resolfactor parameter to give the value used in the Biais-Resolution relation, in the case of the 
 * AreaAndBiaisResolMethod (it's not used otherwise)
 */
void START::HandleResolArea::SetESafeThresholdCondition(ESafeThresholdMethod emeth, double areamin_or_areacondition, double resolfactor)
{
	switch(emeth) {
	case AreaMaxMethod :
		fSafeThresholdMethod = AreaMaxMethod;
		fFractionMaxArea = areamin_or_areacondition;
		break;
	case AreaAndBiaisResolMethod :
		fSafeThresholdMethod = AreaAndBiaisResolMethod;
		fAreaMin = areamin_or_areacondition;
		fResolFactor = resolfactor;
		break;
	}
}

/**
 * \brief Function that will fill the ESafeThreshold in each bands.
 *
 * The principle is to loop over each band, and send this band to the fonction that trully computes the Safe Threshold.
 * The function use is according to the enum ESafeThresholdMethod in the member fSafeThresholdMethod.
 * Use of a switch in order to be easy to modify to add a new method if someone has a good idea.
 *
 */
void START::HandleResolArea::SetInBandsMCEnergyThreshold()
{
	DEBUG_OUT << "VIM : Tu dois modifier cette fonction !!!" << std::endl;
  
	//fAreaMin=0.0;
	//fResolFactor;
	for(int iband=0; iband<(int)(*fBandArray).size(); iband++) {
		if((*fBandArray)[iband].GetKeepBand()==0) continue;

		switch(fSafeThresholdMethod) {
		case AreaMaxMethod :
			ComputeSafeThresholdFromAreaMaxFraction((*fBandArray)[iband],fFractionMaxArea);
			break;
		case AreaAndBiaisResolMethod :
			ComputeSafeThresholdFromAreaAndResolFactor((*fBandArray)[iband],fAreaMin,fResolFactor);
			break;
		default :
			WARN_OUT << "I DON'T KNOW WHICH METHOD TO APPLY FOR THE SAFE THRESHOLD..." << std::endl;
			exit(EXIT_FAILURE);
		}
	}
}
//surchage de la fonction precedente avec un objet band
void START::HandleResolArea::SetInBandsMCEnergyThreshold(std::vector<Band> &InputBandArray)
{
	DEBUG_OUT << "VIM : Tu dois modifier cette fonction !!!" << std::endl;
  
	//fAreaMin=0.0;
	//fResolFactor;
	for(int iband=0; iband<InputBandArray.size(); iband++) {
		if(InputBandArray[iband].GetKeepBand()==0) continue;

		switch(fSafeThresholdMethod) {
		case AreaMaxMethod :
			ComputeSafeThresholdFromAreaMaxFraction(InputBandArray[iband],fFractionMaxArea);
			break;
		case AreaAndBiaisResolMethod :
			ComputeSafeThresholdFromAreaAndResolFactor(InputBandArray[iband],fAreaMin,fResolFactor);
			break;
		default :
			WARN_OUT << "I DON'T KNOW WHICH METHOD TO APPLY FOR THE SAFE THRESHOLD..." << std::endl;
			exit(EXIT_FAILURE);
		}
	}
}

/**
 * \brief Determine MC energy threshold at save it in bands.
 *
 * We find where the area is above AreaMin (Dichotomie Method). 
 * After we take a look at if the biais is less than fResolFactor*resolution  (Dichotomie Method).
 *
 *
 * \param InputBand Reference to a given band. This Band is update at the end of this function with the computed SafeThreshold
 * \param areamin The Minimal Area required (in hectare)
 * \param resolfactor Parameter used as a moderator for the biais (used in the condition biais < resolfactor * resolution )
 *
 */
void START::HandleResolArea::ComputeSafeThresholdFromAreaAndResolFactor(Band &InputBand, double areamin, double resolfactor)
{
	DEBUG_OUT << "Dans ComputeSafeThresholdFromAreaAndResolFactor !!!" << std::endl;
  
	// VIM : WE CARE FOR THE ENERGY FOR WHICH THE AREA MATCH THE CONDITION Area(E) = fAreaMin !
	std::vector<double> energy = InputBand.GetVectorEnergy();
	std::vector<double> area   = InputBand.GetVectorArea();
  
	if (DEBUG) {
		double areamax=0.;
		int imax=0;
		for (unsigned int i=0;i<area.size();++i) {
			if (areamax<=area[i]) {
				areamax = area[i];
				imax = i;
			}
		}
		std::cout << "VIM DEBUG> The maximum value of the area for this band is : " 
				  << areamax << " at the position imax = " 
				  << imax << " !!!" << std::endl;
	}
    
	if (DEBUG) {
		for (unsigned int i=0;i<energy.size();++i) {
			std::cout << "i = " << i << " energy = " << energy[i] 
					  << " area = " << area[i]/1.e4 << std::endl;
		}
	}
    
    if (energy.size()==0 || area.size()==0) {
		std::cout << "WARNING : This band InputBand (run " << InputBand.GetNbRun() 
				  << " ) has not initialized vector... --> Skip it (i.e ethres = 1.e6)" << std::endl; 
		InputBand.SetEthMC(1.e6);
		InputBand.SetKeepBand(1);
		return;
    }
    
    // VIM : We Transform the area vector in order to have the new function vector : 
    std::vector<double> area_minus_refarea;
    std::vector<double> energy_area_minus_refarea;
    
    area_minus_refarea.reserve(area.size());
    energy_area_minus_refarea.reserve(energy.size());
  
    for (unsigned int ii=0;ii<area.size();++ii) {
		double area_tmp = (area[ii]/1.e4)-areamin; // area is in m^2 areamin in hectares !
		if (area_tmp==0.) continue; // VIM : 0 is nearly impossible, or areamin = 0 and you are in a undefined range ! 
		area_minus_refarea.push_back( area_tmp ); // area is in m^2 areamin in hectares !
		energy_area_minus_refarea.push_back( energy[ii] );
    }
    
    // VIM : DEBUG LINE : 
    if (DEBUG) {
		for (unsigned int ii=0;ii<energy_area_minus_refarea.size();++ii) {
			std::cout << "energy_area = " << energy_area_minus_refarea[ii] << " area = " << area_minus_refarea[ii] << std::endl;
		}
    }
    
    double Ener_ForArea;
    int zeroexist_ForArea; // contain the number of times that the value cross 0. MAYBE WE HAVE TO PUT A TEST ON IT !
    Ener_ForArea = FindZeroInBandArray(energy_area_minus_refarea,area_minus_refarea,zeroexist_ForArea);
    
    if (zeroexist_ForArea==0 && area_minus_refarea[0]<=0.) {
		DEBUG_OUT << "VIM : This mean that we never go up to " << areamin << " Hectares in this band InputBand, with Zen = " 
				  << InputBand.GetZenON() << " Offset = " << InputBand.GetOffset() << " Efficiency = " << InputBand.GetEff() << std::endl;
		std::cout << "The surface is never above " << areamin << " Hectares for the band : InputBand,with Zen = " 
				  << InputBand.GetZenON() << " Offset = " << InputBand.GetOffset() << " Efficiency = " << InputBand.GetEff() 
				  << " --> We skip this run ! " << std::endl; 
		InputBand.SetEthMC(1.e6);
		InputBand.SetKeepBand(1); // VIM : I don't know if it is usefull
		return;
    }


    
    // VIM : NOW WE CARE ABOUT THE RESOLUTION (FOR A POINT LIKE SOURCE,  NOT DEFINE YET FOR A EXTENDED SOURCES) : 
    std::vector<double> resol  = InputBand.GetVectorResolution();
    std::vector<double> biais  = InputBand.GetVectorBiais();
    
    if (resol.size()==0 || biais.size()==0) {
		std::cout << "WARNING : This band InputBand (run " << InputBand.GetNbRun() 
				  << " ) has not initialized vector for biais and/or resol ... --> Skip Run (i.e ethres = 1.e6)" 
				  << std::endl;

		InputBand.SetEthMC(1.e6);
		InputBand.SetKeepBand(1); // VIM : I don't know if it is usefull
		return;
    }

    ROOT::Math::Interpolator InterpolResol(energy,resol,ROOT::Math::Interpolation::kAKIMA);
    ROOT::Math::Interpolator InterpolBiais(energy,biais,ROOT::Math::Interpolation::kAKIMA);
    double Resol_AtEnerForArea = InterpolResol.Eval(Ener_ForArea);
    double Biais_AtEnerForArea = TMath::Abs(InterpolBiais.Eval(Ener_ForArea));

    //    if (Biais_AtEnerForArea<resolfactor*Resol_AtEnerForArea || fForceUseDistribution ) {
    if (Biais_AtEnerForArea<resolfactor*Resol_AtEnerForArea || InputBand.GetUseOfInstrumentEnergyDistribution() ) {
		// VIM : Stricctly inferior, because when the value is not filled, the value is put at 0 ! 
		// VIM : In the case of the distribution, the first condition does not really mean something
		DEBUG_OUT << "AreaAndBiaisResolMethod : Ethres = " << Ener_ForArea << " For Zen = " 
				  << InputBand.GetZenON() << " Offset = " << InputBand.GetOffset() 
				  << " Efficiency = " << InputBand.GetEff() << std::endl;
		InputBand.SetEthMC(Ener_ForArea);
		return;
    }
    
    std::vector<double> biais_minus_resolfactor_resol;
    std::vector<double> energy_for_biais_minus_resolfactor_resol;
    biais_minus_resolfactor_resol.reserve(resol.size());
    energy_for_biais_minus_resolfactor_resol.reserve(energy.size());
    
    double bi_tmp = Biais_AtEnerForArea-resolfactor*Resol_AtEnerForArea;
    
    if (bi_tmp!=0.) {
		biais_minus_resolfactor_resol.push_back(bi_tmp);
		energy_for_biais_minus_resolfactor_resol.push_back(Ener_ForArea);
    }
    
    for (unsigned int ii=0;ii<resol.size();++ii) {
		if (energy[ii]<Ener_ForArea) continue;
		double bi_minus_resol = TMath::Abs(biais[ii])-resolfactor*resol[ii];
		if (bi_minus_resol==0.) continue; // VIM : The value 0 is only reach when biais and resolution is equal to 0. (i.e undefined).
		biais_minus_resolfactor_resol.push_back(bi_minus_resol);
		energy_for_biais_minus_resolfactor_resol.push_back(energy[ii]);
    }
    
    // VIM : DEBUG LINE : 
    if (DEBUG) {
		for (unsigned int ii=0;ii<energy_for_biais_minus_resolfactor_resol.size();++ii) {
			std::cout << "energy_resol = " << energy_for_biais_minus_resolfactor_resol[ii] 
					  << " biais_minus_resol = " << biais_minus_resolfactor_resol[ii] << std::endl;
		}
    }
    
    double Ener_ForResolTest;
    int zeroexist_ForResolTest;
    Ener_ForResolTest = FindZeroInBandArray(energy_for_biais_minus_resolfactor_resol,biais_minus_resolfactor_resol,zeroexist_ForResolTest);
    
    if (zeroexist_ForResolTest && biais_minus_resolfactor_resol[0]>0.) {
		std::cout << "WARNING : This band InputBand (run " << InputBand.GetNbRun() 
				  << " ) has not initialized vector for biais and/or resol ... --> Skip Run (i.e ethres = 1.e6)" 
				  << std::endl;
		InputBand.SetEthMC(1.e6);
		InputBand.SetKeepBand(1); // VIM : I don't know if it is usefull
		return;
    }
    
    DEBUG_OUT << "AreaAndBiaisResolMethod : Ethres = " << Ener_ForResolTest << " For Zen = " 
			  << InputBand.GetZenON() << " Offset = " << InputBand.GetOffset() << " Efficiency = " 
			  << InputBand.GetEff() << std::endl;
    InputBand.SetEthMC(Ener_ForResolTest);
    
}






/**
 * \brief Method that compute the safe threshold by taking a certain fraction of the maximum of the area of a given Band.
 * The minimum is found by using a dichotomie method
 * VIM Note : If I am correct this is swhat is done in ParisAnalysis.
 *
 * \param InputBand Reference to a given band. This Band is update at the end of this function with the computed SafeThreshold
 * \param areamaxfraction Fraction of the maximum value of the area taken to find the safe energy threshold (In percentage)
 *
 **/

void START::HandleResolArea::ComputeSafeThresholdFromAreaMaxFraction(Band &InputBand, double areamaxfraction)
{
	DEBUG_OUT << "Dans ComputeSafeThresholdFromAreaMaxFraction !!!" << std::endl;
  
	// VIM : WE CARE FOR THE ENERGY FOR WHICH THE AREA MATCH THE CONDITION areamin = areamaxfraction * Max(Area(E))  !
	std::vector<double> energy = InputBand.GetVectorEnergy();
	std::vector<double> area   = InputBand.GetVectorArea();
  
	if (energy.size()==0 || area.size()==0) {
		std::cout << "WARNING : This band InputBand (run " << InputBand.GetNbRun() 
				  << " ) has not initialized vector... --> Skip Run (i.e ethres = 1.e6)" << std::endl; 
		InputBand.SetEthMC(1.e6);
		InputBand.SetKeepBand(1);
		//continue;
		return;
	}

	// We evaluate here the maximum of the function
	double areamax=0.;
	int imax=0;
	for (unsigned int i=0;i<area.size();++i) {
		if (areamax<=area[i]) {
			areamax = area[i];
			imax = i;
		}
	}
	double areamin = areamax*(areamaxfraction/100.);
  
	DEBUG_OUT << "The maximum value of the area for this band is : " << areamax/1.e4 
			  << " ha at the position imax = " << imax << " So the minimal area you take is " 
			  << areamin/1.e4 << " ha !!!" << std::endl;
  
  
	// VIM : We Transform the area vector in order to have the new function vector : 
	std::vector<double> area_minus_refarea;
	std::vector<double> energy_area_minus_refarea;
	area_minus_refarea.reserve(area.size());
	energy_area_minus_refarea.reserve(energy.size());
  
	for (unsigned int ii=0;ii<area.size();++ii) {
		double area_tmp = area[ii]-areamin ;
		if (area_tmp==0.) continue; // VIM : 0 is nearly impossible, or fAreaMin = 0 and you are in a undefined range ! 
		area_minus_refarea.push_back( area_tmp ); 
		energy_area_minus_refarea.push_back( energy[ii] );
	}
  
	// VIM : DEBUG LINE : 
	if (DEBUG) {
		for (unsigned int ii=0;ii<energy_area_minus_refarea.size();++ii) {
			std::cout << "energy_area = " << energy_area_minus_refarea[ii] << " area = " << area_minus_refarea[ii] << std::endl;
		}
	}
  
	double Ener_ForArea;
	int zeroexist_ForArea; // contain the number of times that the value cross 0. MAYBE WE HAVE TO PUT A TEST ON IT !
	Ener_ForArea = FindZeroInBandArray(energy_area_minus_refarea,area_minus_refarea,zeroexist_ForArea);
  
	if (zeroexist_ForArea==0 && area_minus_refarea[0]<=0.) {
		DEBUG_OUT << "VIM : This mean that we never go up to " << areamin/1.e4 
				  << " Hectares in this band InputBand, with Zen = " << InputBand.GetZenON() << " Offset = " 
				  << InputBand.GetOffset() << " Efficiency = " << InputBand.GetEff() << std::endl;
		std::cout << "The surface is never above " << areamin/1.e4 << " Hectares for the band : InputBand,with Zen = " 
				  << InputBand.GetZenON() << " Offset = " << InputBand.GetOffset() << " Efficiency = " 
				  << InputBand.GetEff() << " --> We skip this run ! " << std::endl; 
		InputBand.SetEthMC(1.e6);
		InputBand.SetKeepBand(1); // VIM : I don't know if it is usefull
		//continue;
		return;
	}

	DEBUG_OUT << "AreaMaxMethod : Ethres = " << Ener_ForArea << " For Zen = " 
			  << InputBand.GetZenON() << " Offset = " << InputBand.GetOffset() 
			  << " Efficiency = " << InputBand.GetEff() << std::endl;
	InputBand.SetEthMC(Ener_ForArea);
}


/**
 * \brief Find min Mc value of band's efficiency
 */
double START::HandleResolArea::FindEffMinMC(Band myBand)
{
	double eff(0);
	//MonteCarlo MC;
	std::vector<double> effMC = fMc->GetEfficiency();
	bool foundit=false;

	// JLK add, in case there is only one value
	if(effMC.size()==1)
		return effMC[0];

	// VIM : To say that we assume the minimum for the MC thing (warning, we assume that the first element of MCis the lowest value
	for(int i(0); i<(int)effMC.size()-1; ++i) {
		if(foundit==true) continue;
		// VIM MERGE NOTES : Need to document, and put warnings
		if(myBand.GetEff()<effMC[0]) {
			eff = effMC[0];
			foundit = true;
		}
		else if (myBand.GetEff()>=effMC[effMC.size()-1]) {
			eff = effMC[effMC.size()-1];
			foundit = true;
		}
		else if (myBand.GetEff()>=effMC[i] && myBand.GetEff()< effMC[i+1]) {
			eff = effMC[i];
			foundit=true;
		}
		else {
		}
	}
	return eff;
}

/**
 * \brief Find max Mc value of band's efficacity
 */
double START::HandleResolArea::FindEffMaxMC(Band myBand)
{
	double eff(0);
	//MonteCarlo MC;
	std::vector<double> effMC = fMc->GetEfficiency();
	bool foundit=false;

		// JLK add, in case there is only one value
	if(effMC.size()==1)
		return effMC[0];
	
	for(int i(0); i<(int)effMC.size()-1; ++i) {
		if(foundit==true) continue;
		// VIM MERGE NOTES : Need documentation and warnings
		if(myBand.GetEff()<effMC[0]) {
			eff = effMC[0];
			foundit = true;
		}
		else if (myBand.GetEff()>=effMC[effMC.size()-1]) {
			eff = effMC[effMC.size()-1];
			foundit = true;
		}
		else if (myBand.GetEff()>=effMC[i] && myBand.GetEff()< effMC[i+1]) {
			eff = effMC[i+1];
			foundit=true;
		}
		else {
		}
	}
	return eff;
}

/**
 * \brief Find min Mc value of band's offset
 */
double START::HandleResolArea::FindOffMinMC(Band myBand)
{
	double off(0);
	//MonteCarlo MC;
	std::vector<double> offMC = fMc->GetOffset();
	bool foundit=false;

	// JLK add, in case there is only one value
	if(offMC.size()==1)
		return offMC[0];
	
	for(int i(0); i<(int)offMC.size()-1; ++i) {
		if(foundit==true) continue;
		if(myBand.GetOffset()>=offMC[i] && myBand.GetOffset()< offMC[i+1]) {
			off = offMC[i];
			foundit=true;
		}
	}
	return off;
}

/**
 * \brief Find max Mc value of band's offset
 */
double START::HandleResolArea::FindOffMaxMC(Band myBand)
{
	double off(0);
	//MonteCarlo MC;
	std::vector<double> offMC = fMc->GetOffset();
	//  std::cout << "FindOffMaxMC" << " size: " << offMC.size() <<   std::endl;
	bool foundit=false;

	// JLK add, in case there is only one value
	if(offMC.size()==1)
		return offMC[0];
	
	for(int i(0); i<(int)offMC.size()-1; ++i) {
		//    std::cout << "offMC " << i << " " << offMC[i] << std::endl;
		if(foundit==true) continue;
		if(myBand.GetOffset()>=offMC[i] && myBand.GetOffset()<offMC[i+1]) {
			off = offMC[i+1];
			foundit=true;
		}
	}

	return off;
}

/**
 * \brief Find min Mc value of band's zenith
 */
double START::HandleResolArea::FindZenMinMC(Band myBand)
{
	double zen(0);
	//MonteCarlo MC;
	std::vector<double> zenMC = fMc->GetZenith();
	bool foundit=false;

	// JLK add, in case there is only one value
	if(zenMC.size()==1)
		return zenMC[0];
	
	for(int i(0); i<(int)zenMC.size()-1; ++i) {
		if(foundit==true) continue;
		if(myBand.GetZenON()>=zenMC[i] && myBand.GetZenON()< zenMC[i+1]) {
			zen = zenMC[i];
			foundit=true;
		}
	}  
	return zen;
}

/**
 * \brief Find max Mc value of band's zenith
 */
double START::HandleResolArea::FindZenMaxMC(Band myBand)
{
	double zen(0);
	//MonteCarlo MC;
	std::vector<double> zenMC = fMc->GetZenith();
	bool foundit=false;

	// JLK add, in case there is only one value
	if(zenMC.size()==1)
		return zenMC[0];
	
	for(int i(0); i<(int)zenMC.size()-1; ++i) {
		if(foundit==true) continue;
		if(myBand.GetZenON()>=zenMC[i] && myBand.GetZenON()< zenMC[i+1]) {
			zen = zenMC[i+1];
			foundit=true;
		}
	}
	return zen;
}

/**
 * \brief Draw all IRF
 *
 * \warning It's a lot of plots...
 */
void START::HandleResolArea::DrawAllIRF()
{

	buildirfindividualscanvas = true;

	std::vector<double> veff = fMc->GetEfficiency();
	std::vector<double> vzen = fMc->GetZenith();
	std::vector<double> voff = fMc->GetOffset();
	std::vector<TString> vteltype = fMc->GetTelType();

	for(std::vector<double>::const_iterator ieff=veff.begin(); ieff!=veff.end(); ++ieff) {

		for(std::vector<double>::const_iterator izen=vzen.begin(); izen!=vzen.end(); ++izen) {

			for(std::vector<double>::const_iterator ioff=voff.begin(); ioff!=voff.end(); ++ioff) {

				for(std::vector<TString>::const_iterator itel=vteltype.begin(); itel!=vteltype.end(); ++itel) {

					// JLK hack : je vais faire un truc bien dans MonteCarlo, en attendant je fais ça
					unsigned int telcode;
					if(*itel=="") telcode=30;
					else if(*itel=="3Tel") telcode=1;
					DrawIRF(*ieff,*izen,*ioff,telcode);

				}

			}

		}

	}

}

/**
 * \brief Draw IRF for given parameters
 * \param eff efficiency
 * \param zen zenith angle
 * \param off offset
 * \param telcode telcode (ie 1..30)
 */
void START::HandleResolArea::DrawIRF(double eff, double zen, double off, unsigned int telcode)
{

	buildirfindividualscanvas = true;

	//MonteCarlo MC;

	if(!fMc->IsItMCEfficiency(eff)) {
		WARN_OUT << "efficiency "<< eff << " is not an MC value" << std::endl;
		fMc->PrintMCParameters();
		return;
	}
	if(!fMc->IsItMCZenith(zen)) {
		WARN_OUT << "zenith "<< zen << " is not an MC value" << std::endl;
		fMc->PrintMCParameters();
		return;
	}
	if(!fMc->IsItMCOffset(off)) {
		WARN_OUT << "offset "<< off << " is not an MC value" << std::endl;
		fMc->PrintMCParameters();
		return;
	}
	if(!fMc->IsItMCTelCode(telcode)) {
		WARN_OUT << "telcode "<< telcode << " is not an MC value" << std::endl;
		fMc->PrintMCParameters();
		return;
	}

	BuildIRFEffectiveAreaCanvas(eff,zen,off,telcode);
	BuildIRFResolutionCanvas(eff,zen,off,telcode);
	BuildIRFBiaisCanvas(eff,zen,off,telcode);

	fMapIRFEffectiveAreaCanvas[MakeKeyIRFEffectiveAreaCanvas(eff,zen,off,telcode)]->Draw();
	fMapIRFEffectiveAreaCanvas[MakeKeyIRFEffectiveAreaCanvas(eff,zen,off,telcode)]->cd();
	fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,telcode)]->Draw("AP");

	fMapIRFResolutionCanvas[MakeKeyIRFResolutionCanvas(eff,zen,off,telcode)]->Draw();
	fMapIRFResolutionCanvas[MakeKeyIRFResolutionCanvas(eff,zen,off,telcode)]->cd();
	fMapIRFResolutionGraph[MakeKeyIRFResolutionGraph(eff,zen,off,telcode)]->Draw("AP");

	fMapIRFBiaisCanvas[MakeKeyIRFBiaisCanvas(eff,zen,off,telcode)]->Draw();
	fMapIRFBiaisCanvas[MakeKeyIRFBiaisCanvas(eff,zen,off,telcode)]->cd();
	fMapIRFBiaisGraph[MakeKeyIRFBiaisGraph(eff,zen,off,telcode)]->Draw("AP");

}

void START::HandleResolArea::BuildIRFBiaisCanvas(double eff, double zen, double off, unsigned int telcode) {

	//MonteCarlo MC;

	TString canvasbiaisname = MakeKeyIRFBiaisCanvas(eff,zen,off,telcode);
	TString canvasbiaistitle = fConfigName;
	canvasbiaistitle+="_";
	canvasbiaistitle+="Biais_eff";
	canvasbiaistitle+=eff;
	canvasbiaistitle+="_zen";
	canvasbiaistitle+=zen;
	canvasbiaistitle+="_off";
	canvasbiaistitle+=off;
	canvasbiaistitle+=fMc->GetStringFromTelcode(telcode);

	if(buildirfindividualscanvas) {

		if(!fMapIRFBiaisCanvas[canvasbiaisname]) {
			fMapIRFBiaisCanvas[canvasbiaisname] = new TCanvas(canvasbiaisname,canvasbiaistitle);
		}
		else {
			delete fMapIRFBiaisCanvas[canvasbiaisname];
			fMapIRFBiaisCanvas[canvasbiaisname] = 0;
			fMapIRFBiaisCanvas[canvasbiaisname] = new TCanvas(canvasbiaisname,canvasbiaistitle);
		}
		fMapIRFBiaisCanvas[canvasbiaisname]->SetLogx();
		fMapIRFBiaisCanvas[canvasbiaisname]->SetGridx();
		fMapIRFBiaisCanvas[canvasbiaisname]->SetGridy();
    
	}

	std::vector<double> xbiais, ybiais;

	TString graphname = MakeKeyIRFBiaisGraph(eff,zen,off,telcode);
	if(!fMapIRFBiaisGraph[graphname]) {
		fMapIRFBiaisGraph[graphname] = new TGraph();
	}
	else {
		delete fMapIRFBiaisGraph[graphname];
		fMapIRFBiaisGraph[graphname] = 0;
		fMapIRFBiaisGraph[graphname] = new TGraph();
	}

	if(fmapVectorsBiais[MakeKeyRun(eff,zen,off,fMc->GetStringFromTelcode(telcode))].first.size()>0) {
    
		xbiais = fmapVectorsBiais[MakeKeyRun(eff,zen,off,fMc->GetStringFromTelcode(telcode))].first;
		ybiais = fmapVectorsBiais[MakeKeyRun(eff,zen,off,fMc->GetStringFromTelcode(telcode))].second;
    
		for(unsigned int ipoint(0); ipoint<xbiais.size(); ipoint++) {
			fMapIRFBiaisGraph[graphname]->SetPoint(ipoint,xbiais[ipoint],ybiais[ipoint]);
		}

		TString grbiaistitle=canvasbiaistitle;
		fMapIRFBiaisGraph[graphname]->SetName(MakeKeyRun(eff,zen,off,fMc->GetStringFromTelcode(telcode)).Data());
		fMapIRFBiaisGraph[graphname]->SetTitle(grbiaistitle.Data());
		fMapIRFBiaisGraph[graphname]->SetMarkerStyle(2);
		fMapIRFBiaisGraph[graphname]->SetMarkerColor(2);
		fMapIRFBiaisGraph[graphname]->SetFillStyle(0);
		fMapIRFBiaisGraph[graphname]->SetFillColor(0);
		fMapIRFBiaisGraph[graphname]->GetHistogram()->GetXaxis()->SetTitle("Log E (TeV)");
		fMapIRFBiaisGraph[graphname]->GetHistogram()->GetYaxis()->SetTitle("Biais");
		fMapIRFBiaisGraph[graphname]->GetHistogram()->GetXaxis()->CenterTitle();
		fMapIRFBiaisGraph[graphname]->GetHistogram()->GetYaxis()->CenterTitle();
		fMapIRFBiaisGraph[graphname]->GetHistogram()->SetMinimum(-1.);
		fMapIRFBiaisGraph[graphname]->GetHistogram()->SetMaximum(1.);

	}
	else {
		WARN_OUT << "Can't find map for biais :" << MakeKeyRun(eff,zen,off,fMc->GetStringFromTelcode(telcode)) << std::endl;
	}

}

void START::HandleResolArea::BuildIRFResolutionCanvas(double eff, double zen, double off, unsigned int telcode) {

	//MonteCarlo MC;

	TString canvasresolutionname = MakeKeyIRFResolutionCanvas(eff,zen,off,telcode);
	TString canvasresolutiontitle = fConfigName;
	canvasresolutiontitle+="_";
	canvasresolutiontitle+="Resolution_eff";
	canvasresolutiontitle+=eff;
	canvasresolutiontitle+="_zen";
	canvasresolutiontitle+=zen;
	canvasresolutiontitle+="_off";
	canvasresolutiontitle+=off;
	canvasresolutiontitle+=fMc->GetStringFromTelcode(telcode);

	if(buildirfindividualscanvas) {

		if(!fMapIRFResolutionCanvas[canvasresolutionname]) {
			fMapIRFResolutionCanvas[canvasresolutionname] = new TCanvas(canvasresolutionname,canvasresolutiontitle);
		}
		else {
			delete fMapIRFResolutionCanvas[canvasresolutionname];
			fMapIRFResolutionCanvas[canvasresolutionname] = 0;
			fMapIRFResolutionCanvas[canvasresolutionname] = new TCanvas(canvasresolutionname,canvasresolutiontitle);
		}
		fMapIRFResolutionCanvas[canvasresolutionname]->SetLogx();
		fMapIRFResolutionCanvas[canvasresolutionname]->SetGridx();
		fMapIRFResolutionCanvas[canvasresolutionname]->SetGridy();
    
	}

	TString graphname = MakeKeyIRFResolutionGraph(eff,zen,off,telcode);
	if(!fMapIRFResolutionGraph[graphname]) {
		fMapIRFResolutionGraph[graphname] = new TGraph();
	}
	else {
		delete fMapIRFResolutionGraph[graphname];
		fMapIRFResolutionGraph[graphname] = 0;
		fMapIRFResolutionGraph[graphname] = new TGraph();
	}

	std::vector<double> xresolution, yresolution;

	if(fmapVectorsResol[MakeKeyRun(eff,zen,off,fMc->GetStringFromTelcode(telcode))].first.size()>0) {
    
		xresolution = fmapVectorsResol[MakeKeyRun(eff,zen,off,fMc->GetStringFromTelcode(telcode))].first;
		yresolution = fmapVectorsResol[MakeKeyRun(eff,zen,off,fMc->GetStringFromTelcode(telcode))].second;
    
		for(unsigned int ipoint(0); ipoint<xresolution.size(); ipoint++) {
			fMapIRFResolutionGraph[graphname]->SetPoint(ipoint,xresolution[ipoint],yresolution[ipoint]);
		}

		TString grresolutiontitle=canvasresolutiontitle;
		fMapIRFResolutionGraph[graphname]->SetName(MakeKeyRun(eff,zen,off,fMc->GetStringFromTelcode(telcode)).Data());
		fMapIRFResolutionGraph[graphname]->SetTitle(grresolutiontitle.Data());
		fMapIRFResolutionGraph[graphname]->SetMarkerStyle(2);
		fMapIRFResolutionGraph[graphname]->SetMarkerColor(2);
		fMapIRFResolutionGraph[graphname]->SetFillStyle(0);
		fMapIRFResolutionGraph[graphname]->SetFillColor(0);
		fMapIRFResolutionGraph[graphname]->GetHistogram()->GetXaxis()->SetTitle("Log E (TeV)");
		fMapIRFResolutionGraph[graphname]->GetHistogram()->GetYaxis()->SetTitle("Resolution");
		fMapIRFResolutionGraph[graphname]->GetHistogram()->GetXaxis()->CenterTitle();
		fMapIRFResolutionGraph[graphname]->GetHistogram()->GetYaxis()->CenterTitle();
		fMapIRFResolutionGraph[graphname]->GetHistogram()->SetMinimum(0.);
		fMapIRFResolutionGraph[graphname]->GetHistogram()->SetMaximum(1.);  
    
	}
	else {
		WARN_OUT << "Can't find map for resolution :" << MakeKeyRun(eff,zen,off,fMc->GetStringFromTelcode(telcode)) << std::endl;
	}

}

void START::HandleResolArea::BuildIRFEffectiveAreaCanvas(double eff, double zen, double off, unsigned int telcode) {

	//MonteCarlo MC;

	TString canvasareaname = MakeKeyIRFEffectiveAreaCanvas(eff,zen,off,telcode);
	TString canvasareatitle = fConfigName;
	canvasareatitle+="_";
	canvasareatitle+="EffectiveArea_eff";
	canvasareatitle+=eff;
	canvasareatitle+="_zen";
	canvasareatitle+=zen;
	canvasareatitle+="_off";
	canvasareatitle+=off;
	canvasareatitle+=fMc->GetStringFromTelcode(telcode);

	if(buildirfindividualscanvas) {

		if(!fMapIRFEffectiveAreaCanvas[canvasareaname]) {
			fMapIRFEffectiveAreaCanvas[canvasareaname] = new TCanvas(canvasareaname,canvasareatitle);
		}
		else {
			delete fMapIRFEffectiveAreaCanvas[canvasareaname];
			fMapIRFEffectiveAreaCanvas[canvasareaname] = 0;
			fMapIRFEffectiveAreaCanvas[canvasareaname] = new TCanvas(canvasareaname,canvasareatitle);
		}
		fMapIRFEffectiveAreaCanvas[canvasareaname]->SetLogx();
		fMapIRFEffectiveAreaCanvas[canvasareaname]->SetLogy();
		fMapIRFEffectiveAreaCanvas[canvasareaname]->SetGridx();
		fMapIRFEffectiveAreaCanvas[canvasareaname]->SetGridy();
    
	}

	TString graphname = MakeKeyIRFEffectiveAreaGraph(eff,zen,off,telcode);
	if(!fMapIRFEffectiveAreaGraph[graphname]) {
		fMapIRFEffectiveAreaGraph[graphname] = new TGraph();
	}
	else {
		delete fMapIRFEffectiveAreaGraph[graphname];
		fMapIRFEffectiveAreaGraph[graphname] = 0;
		fMapIRFEffectiveAreaGraph[graphname] = new TGraph();
	}

	std::vector<double> xarea, yarea;

	if(fmapVectorsArea[MakeKeyRun(eff,zen,off,fMc->GetStringFromTelcode(telcode))].first.size()>0) {
    
		xarea = fmapVectorsArea[MakeKeyRun(eff,zen,off,fMc->GetStringFromTelcode(telcode))].first;
		yarea = fmapVectorsArea[MakeKeyRun(eff,zen,off,fMc->GetStringFromTelcode(telcode))].second;
    
		for(unsigned int ipoint(0); ipoint<xarea.size(); ipoint++) {
			fMapIRFEffectiveAreaGraph[graphname]->SetPoint(ipoint,xarea[ipoint],yarea[ipoint]);
		}

		TString grareatitle=canvasareatitle;
		fMapIRFEffectiveAreaGraph[graphname]->SetName(MakeKeyRun(eff,zen,off,fMc->GetStringFromTelcode(telcode)).Data());
		fMapIRFEffectiveAreaGraph[graphname]->SetTitle(grareatitle.Data());
		fMapIRFEffectiveAreaGraph[graphname]->SetMarkerStyle(2);
		fMapIRFEffectiveAreaGraph[graphname]->SetMarkerColor(2);
		fMapIRFEffectiveAreaGraph[graphname]->SetFillStyle(0);
		fMapIRFEffectiveAreaGraph[graphname]->SetFillColor(0);
		fMapIRFEffectiveAreaGraph[graphname]->GetHistogram()->GetXaxis()->SetTitle("Log E (TeV)");
		fMapIRFEffectiveAreaGraph[graphname]->GetHistogram()->GetYaxis()->SetTitle("Effective area (m^{2})");
		fMapIRFEffectiveAreaGraph[graphname]->GetHistogram()->GetXaxis()->CenterTitle();
		fMapIRFEffectiveAreaGraph[graphname]->GetHistogram()->GetYaxis()->CenterTitle();
		fMapIRFEffectiveAreaGraph[graphname]->GetHistogram()->SetMinimum(1.e0);
		fMapIRFEffectiveAreaGraph[graphname]->GetHistogram()->SetMaximum(1.e6); 
    
	}
	else {
		WARN_OUT << "Can't find map effective area :" << MakeKeyRun(eff,zen,off,fMc->GetStringFromTelcode(telcode)) << std::endl;
	}

}

/**
 * \brief Draw plot with interpolated IRF from band and IRF used for interpolation
 * \param BandToDraw Band
 */
void START::HandleResolArea::DrawIRFBandDiagnostic(Band &BandToDraw) {

	buildirfindividualscanvas = false;

	BuildBandEffectiveAreaCanvas(BandToDraw);

}

/**
 * \brief Draw plot with interpolated IRF from band and IRF used for interpolation
 * \param BandToDraw Band
 */
void START::HandleResolArea::BuildBandEffectiveAreaCanvas(Band &BandToDraw) {

	//MonteCarlo MC;

	std::pair<double,double> paireff = fMc->GetMCLimitingEfficiency(BandToDraw.GetEff());
	double effmin(paireff.first), effmax(paireff.second);

	std::pair<double,double> pairzen = fMc->GetMCLimitingZenith(BandToDraw.GetZenON());
	double zenmin(pairzen.first), zenmax(pairzen.second);

	std::pair<double,double> pairoff = fMc->GetMCLimitingOffset(BandToDraw.GetOffset());
	double offmin(pairoff.first), offmax(pairoff.second);

	//// Area

	// Canvas

	if(!fMapBandIRFEffectiveAreaCanvas[MakeKeyIRFBandEffectiveAreaCanvas(BandToDraw)]) {
		fMapBandIRFEffectiveAreaCanvas[MakeKeyIRFBandEffectiveAreaCanvas(BandToDraw)] = 
			new TCanvas(MakeKeyIRFBandEffectiveAreaCanvas(BandToDraw),MakeKeyIRFBandEffectiveAreaCanvas(BandToDraw));
	}
	else {
		delete fMapBandIRFEffectiveAreaCanvas[MakeKeyIRFBandEffectiveAreaCanvas(BandToDraw)];
		fMapBandIRFEffectiveAreaCanvas[MakeKeyIRFBandEffectiveAreaCanvas(BandToDraw)] = 0;
		fMapBandIRFEffectiveAreaCanvas[MakeKeyIRFBandEffectiveAreaCanvas(BandToDraw)] = 
			new TCanvas(MakeKeyIRFBandEffectiveAreaCanvas(BandToDraw),MakeKeyIRFBandEffectiveAreaCanvas(BandToDraw));
	}

	fMapBandIRFEffectiveAreaCanvas[MakeKeyIRFBandEffectiveAreaCanvas(BandToDraw)]->SetLogx();
	fMapBandIRFEffectiveAreaCanvas[MakeKeyIRFBandEffectiveAreaCanvas(BandToDraw)]->SetLogy();
	fMapBandIRFEffectiveAreaCanvas[MakeKeyIRFBandEffectiveAreaCanvas(BandToDraw)]->SetGridx();
	fMapBandIRFEffectiveAreaCanvas[MakeKeyIRFBandEffectiveAreaCanvas(BandToDraw)]->SetGridy();
	fMapBandIRFEffectiveAreaCanvas[MakeKeyIRFBandEffectiveAreaCanvas(BandToDraw)]->SetFillColor(0);
	fMapBandIRFEffectiveAreaCanvas[MakeKeyIRFBandEffectiveAreaCanvas(BandToDraw)]->SetFillStyle(0);

	// Multigraph

	if(!fMapBandIRFEffectiveAreaMultiGraph[MakeKeyIRFBandEffectiveAreaMultiGraph(BandToDraw)]) {
		fMapBandIRFEffectiveAreaMultiGraph[MakeKeyIRFBandEffectiveAreaMultiGraph(BandToDraw)] = 
			new TMultiGraph(MakeKeyIRFBandEffectiveAreaMultiGraph(BandToDraw),"Band Diagnostic - EffectiveArea");
	}
	else {
		delete fMapBandIRFEffectiveAreaMultiGraph[MakeKeyIRFBandEffectiveAreaMultiGraph(BandToDraw)];
		fMapBandIRFEffectiveAreaMultiGraph[MakeKeyIRFBandEffectiveAreaMultiGraph(BandToDraw)] = 0;
		fMapBandIRFEffectiveAreaMultiGraph[MakeKeyIRFBandEffectiveAreaMultiGraph(BandToDraw)] = 
			new TMultiGraph(MakeKeyIRFBandEffectiveAreaMultiGraph(BandToDraw),"Band Diagnostic - Area");
	}
  
	// Adding graph to multigraph

	double eff,zen,off;

	eff=effmin;
	zen=zenmin;
	off=offmin;
	if(!fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,BandToDraw.GetTelCode())])
		BuildIRFEffectiveAreaCanvas(eff,zen,off,BandToDraw.GetTelCode());
	fMapBandIRFEffectiveAreaMultiGraph[MakeKeyIRFBandEffectiveAreaMultiGraph(BandToDraw)]->Add(fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,
																																					  zen,
																																					  off,
																																					  BandToDraw.GetTelCode())]);

	fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,BandToDraw.GetTelCode())]->SetMarkerStyle(7);
	fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,BandToDraw.GetTelCode())]->SetMarkerColor(1);

	eff=effmin;
	zen=zenmax;
	off=offmin;
	if(!fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,BandToDraw.GetTelCode())])
		BuildIRFEffectiveAreaCanvas(eff,zen,off,BandToDraw.GetTelCode());
	fMapBandIRFEffectiveAreaMultiGraph[MakeKeyIRFBandEffectiveAreaMultiGraph(BandToDraw)]->Add(fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,
																																					  zen,
																																					  off,
																																					  BandToDraw.GetTelCode())]);

	fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,BandToDraw.GetTelCode())]->SetMarkerStyle(7);
	fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,BandToDraw.GetTelCode())]->SetMarkerColor(2);

	eff=effmin;
	zen=zenmin;
	off=offmax;
	if(!fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,BandToDraw.GetTelCode())])
		BuildIRFEffectiveAreaCanvas(eff,zen,off,BandToDraw.GetTelCode());
	fMapBandIRFEffectiveAreaMultiGraph[MakeKeyIRFBandEffectiveAreaMultiGraph(BandToDraw)]->Add(fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,
																																					  zen,
																																					  off,
																																					  BandToDraw.GetTelCode())]);

	fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,BandToDraw.GetTelCode())]->SetMarkerStyle(7);
	fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,BandToDraw.GetTelCode())]->SetMarkerColor(3);

	eff=effmin;
	zen=zenmax;
	off=offmax;
	if(!fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,BandToDraw.GetTelCode())])
		BuildIRFEffectiveAreaCanvas(eff,zen,off,BandToDraw.GetTelCode());
	fMapBandIRFEffectiveAreaMultiGraph[MakeKeyIRFBandEffectiveAreaMultiGraph(BandToDraw)]->Add(fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,
																																					  zen,
																																					  off,
																																					  BandToDraw.GetTelCode())]);

	fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,BandToDraw.GetTelCode())]->SetMarkerStyle(7);
	fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,BandToDraw.GetTelCode())]->SetMarkerColor(4);

	eff=effmax;
	zen=zenmin;
	off=offmin;
	if(!fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,BandToDraw.GetTelCode())])
		BuildIRFEffectiveAreaCanvas(eff,zen,off,BandToDraw.GetTelCode());
	fMapBandIRFEffectiveAreaMultiGraph[MakeKeyIRFBandEffectiveAreaMultiGraph(BandToDraw)]->Add(fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,
																																					  zen,
																																					  off,
																																					  BandToDraw.GetTelCode())]);

	fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,BandToDraw.GetTelCode())]->SetMarkerStyle(7);
	fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,BandToDraw.GetTelCode())]->SetMarkerColor(5);

	eff=effmax;
	zen=zenmax;
	off=offmin;
	if(!fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,BandToDraw.GetTelCode())])
		BuildIRFEffectiveAreaCanvas(eff,zen,off,BandToDraw.GetTelCode());
	fMapBandIRFEffectiveAreaMultiGraph[MakeKeyIRFBandEffectiveAreaMultiGraph(BandToDraw)]->Add(fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,
																																					  zen,
																																					  off,
																																					  BandToDraw.GetTelCode())]);

	fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,BandToDraw.GetTelCode())]->SetMarkerStyle(7);
	fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,BandToDraw.GetTelCode())]->SetMarkerColor(6);

	eff=effmax;
	zen=zenmin;
	off=offmax;
	if(!fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,BandToDraw.GetTelCode())])
		BuildIRFEffectiveAreaCanvas(eff,zen,off,BandToDraw.GetTelCode());
	fMapBandIRFEffectiveAreaMultiGraph[MakeKeyIRFBandEffectiveAreaMultiGraph(BandToDraw)]->Add(fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,
																																					  zen,
																																					  off,
																																					  BandToDraw.GetTelCode())]);

	fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,BandToDraw.GetTelCode())]->SetMarkerStyle(7);
	fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,BandToDraw.GetTelCode())]->SetMarkerColor(7);

	eff=effmax;
	zen=zenmax;
	off=offmax;
	if(!fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,BandToDraw.GetTelCode())])
		BuildIRFEffectiveAreaCanvas(eff,zen,off,BandToDraw.GetTelCode());
	fMapBandIRFEffectiveAreaMultiGraph[MakeKeyIRFBandEffectiveAreaMultiGraph(BandToDraw)]->Add(fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,
																																					  zen,
																																					  off,
																																					  BandToDraw.GetTelCode())]);

	fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,BandToDraw.GetTelCode())]->SetMarkerStyle(7);
	fMapIRFEffectiveAreaGraph[MakeKeyIRFEffectiveAreaGraph(eff,zen,off,BandToDraw.GetTelCode())]->SetMarkerColor(8);

	// adding graph from band
	TGraph *bandgraph = new TGraph();
	bandgraph->SetName("FixedIRFband");
	bandgraph->SetTitle("Fixed IRF band");
	bandgraph->SetMarkerStyle(21);
	bandgraph->SetMarkerSize(0.8);
	bandgraph->SetMarkerColor(kBlack);
	bandgraph->SetFillStyle(0);
	bandgraph->SetFillColor(0);
	bandgraph->GetHistogram()->GetXaxis()->SetTitle("Log E (TeV)");
	bandgraph->GetHistogram()->GetYaxis()->SetTitle("Effective area (m^{2})");
	bandgraph->GetHistogram()->GetXaxis()->CenterTitle();
	bandgraph->GetHistogram()->GetYaxis()->CenterTitle();
	bandgraph->GetHistogram()->SetMinimum(1.e1);
	bandgraph->GetHistogram()->SetMaximum(1.e6);
	for(unsigned int ipoint(0); ipoint<BandToDraw.GetVectorEnergy().size(); ipoint++) {
		bandgraph->SetPoint(ipoint,BandToDraw.GetVectorEnergy()[ipoint],BandToDraw.GetVectorArea()[ipoint]);
	}

	TGraph *bandintergraph = new TGraph();
	bandintergraph->SetName("InterpolatedIRFband");
	bandintergraph->SetTitle("Interpolated IRF band");
	bandintergraph->SetMarkerStyle(22);
	bandintergraph->SetMarkerSize(0.8);
	bandintergraph->SetMarkerColor(kRed);
	bandintergraph->SetFillStyle(0);
	bandintergraph->SetFillColor(0);
	bandintergraph->GetHistogram()->GetXaxis()->SetTitle("Log E (TeV)");
	bandintergraph->GetHistogram()->GetYaxis()->SetTitle("Effective area (m^{2})");
	bandintergraph->GetHistogram()->GetXaxis()->CenterTitle();
	bandintergraph->GetHistogram()->GetYaxis()->CenterTitle();
	bandintergraph->GetHistogram()->SetMinimum(1.e1);
	bandintergraph->GetHistogram()->SetMaximum(1.e6);
	for(unsigned int ipoint(0); ipoint<BandToDraw.GetVectorInterEnergy().size(); ipoint++) {
		bandintergraph->SetPoint(ipoint,BandToDraw.GetVectorInterEnergy()[ipoint],BandToDraw.GetVectorInterArea()[ipoint]);
	}
  
	fMapBandIRFEffectiveAreaMultiGraph[MakeKeyIRFBandEffectiveAreaMultiGraph(BandToDraw)]->Add(bandintergraph);
	fMapBandIRFEffectiveAreaMultiGraph[MakeKeyIRFBandEffectiveAreaMultiGraph(BandToDraw)]->Add(bandgraph);
  
	// drawing canvas
	fMapBandIRFEffectiveAreaCanvas[MakeKeyIRFBandEffectiveAreaCanvas(BandToDraw)]->cd();

	fMapBandIRFEffectiveAreaMultiGraph[MakeKeyIRFBandEffectiveAreaMultiGraph(BandToDraw)]->Draw("AP");
	fMapBandIRFEffectiveAreaMultiGraph[MakeKeyIRFBandEffectiveAreaMultiGraph(BandToDraw)]->GetHistogram()->SetMinimum(1.e0);
	fMapBandIRFEffectiveAreaMultiGraph[MakeKeyIRFBandEffectiveAreaMultiGraph(BandToDraw)]->GetHistogram()->SetMaximum(1.e6);
	TLegend *leg = 
		fMapBandIRFEffectiveAreaCanvas[MakeKeyIRFBandEffectiveAreaCanvas(BandToDraw)]->BuildLegend(0.39,0.18,0.89,0.61,"");
	leg->SetFillColor(0);

	TLatex *bandheader = new TLatex();
	bandheader->SetNDC();
	bandheader->SetTextSize(0.032);
	bandheader->SetTextColor(1);
	bandheader->SetTextFont(52);
	bandheader->DrawLatex(0.54,0.90,MakeBandHeader(BandToDraw));
	delete bandheader; bandheader=0;

}

TString START::HandleResolArea::MakeBandHeader(Band &BandToDraw) {

	std::ostringstream header;
	header << "Run ";
	header.precision(0);
	header << BandToDraw.GetNbRun() << ", eff=";
	header.precision(2);
	header << BandToDraw.GetEff() << ", zen=";
	header.precision(2);
	header << BandToDraw.GetZenON() << " and off=";
	header.precision(2);
	header << BandToDraw.GetOffset();

	return header.str().c_str();
}

void START::HandleResolArea::BuildBandResolutionCanvas(Band &BandToDraw) {

}

void START::HandleResolArea::BuildBandBiaisCanvas(Band &BandToDraw) {

}

TString START::HandleResolArea::MakeKeyIRFBandEffectiveAreaCanvas(Band &BandToDraw) {

	std::ostringstream oss_key;
	oss_key << "Canvas_Band_EffectiveArea_run" << BandToDraw.GetNbRun();
	oss_key << "_eff";
	oss_key.precision(2);
	oss_key << BandToDraw.GetEff();
	oss_key << "_zen";
	oss_key.precision(2);
	oss_key << BandToDraw.GetZenON();
	oss_key << "_off";
	oss_key.precision(2);
	oss_key << BandToDraw.GetOffset();

	TString key = oss_key.str().c_str();

	return key; 
}

TString START::HandleResolArea::MakeKeyIRFBandResolutionCanvas(Band &BandToDraw) {

	TString key="Canvas_Band_Resolution";
	key+="run";
	key+=BandToDraw.GetNbRun();
	key+="_eff";
	key+=BandToDraw.GetEff();
	key+="_zen";
	key+=BandToDraw.GetZenON();
	key+="_off";
	key+=BandToDraw.GetOffset();
	return key; 
}

TString START::HandleResolArea::MakeKeyIRFBandBiaisCanvas(Band &BandToDraw) {

	TString key="Canvas_Band_Biais";
	key+="run";
	key+=BandToDraw.GetNbRun();
	key+="_eff";
	key+=BandToDraw.GetEff();
	key+="_zen";
	key+=BandToDraw.GetZenON();
	key+="_off";
	key+=BandToDraw.GetOffset();
	return key; 
}

TString START::HandleResolArea::MakeKeyIRFBandBiaisMultiGraph(Band &BandToDraw) {

	TString key="Multigraph_Biais_Band_";
	key+="run";
	key+=BandToDraw.GetNbRun();
	key+="_eff";
	key+=BandToDraw.GetEff();
	key+="_zen";
	key+=BandToDraw.GetZenON();
	key+="_off";
	key+=BandToDraw.GetOffset();
	return key; 
}

TString START::HandleResolArea::MakeKeyIRFBandResolutionMultiGraph(Band &BandToDraw) {

	TString key="Multigraph_Resolution_Band_";
	key+="run";
	key+=BandToDraw.GetNbRun();
	key+="_eff";
	key+=BandToDraw.GetEff();
	key+="_zen";
	key+=BandToDraw.GetZenON();
	key+="_off";
	key+=BandToDraw.GetOffset();
	return key; 
}

TString START::HandleResolArea::MakeKeyIRFBandEffectiveAreaMultiGraph(Band &BandToDraw) {

	TString key="Multigraph_Area_Band_";
	key+="run";
	key+=BandToDraw.GetNbRun();
	key+="_eff";
	key+=BandToDraw.GetEff();
	key+="_zen";
	key+=BandToDraw.GetZenON();
	key+="_off";
	key+=BandToDraw.GetOffset();
	return key; 
}

/**
 * \brief make key for IRF biais canvas
 */
TString START::HandleResolArea::MakeKeyIRFBiaisCanvas(double eff, double zen, double off, unsigned int telcode) {

	TString key="Canvas_Biais_";
	key+=eff;
	key+="_";
	key+=zen;
	key+="_";
	key+=off;
	key+="_";
	key+=telcode;
  
	return key;  

}

/**
 * \brief make key for IRF resolution canvas
 */
TString START::HandleResolArea::MakeKeyIRFResolutionCanvas(double eff, double zen, double off, unsigned int telcode) {

	TString key="Canvas_Resolution_";
	key+=eff;
	key+="_";
	key+=zen;
	key+="_";
	key+=off;
	key+="_";
	key+=telcode;
  
	return key;  

}

/**
 * \brief make key for IRF effective area canvas
 */
TString START::HandleResolArea::MakeKeyIRFEffectiveAreaCanvas(double eff, double zen, double off, unsigned int telcode) {

	TString key="Canvas_Area_";
	key+=eff;
	key+="_";
	key+=zen;
	key+="_";
	key+=off;
	key+="_";
	key+=telcode;
  
	return key;  

}

/**
 * \brief make key for IRF biais canvas
 */
TString START::HandleResolArea::MakeKeyIRFBiaisGraph(double eff, double zen, double off, unsigned int telcode) {

	TString key="Graph_Biais_";
	key+=eff;
	key+="_";
	key+=zen;
	key+="_";
	key+=off;
	key+="_";
	key+=telcode;
  
	return key;  

}

/**
 * \brief make key for IRF resolution canvas
 */
TString START::HandleResolArea::MakeKeyIRFResolutionGraph(double eff, double zen, double off, unsigned int telcode) {

	TString key="Graph_Resolution_";
	key+=eff;
	key+="_";
	key+=zen;
	key+="_";
	key+=off;
	key+="_";
	key+=telcode;
  
	return key;  

}

/**
 * \brief make key for IRF effective area canvas
 */
TString START::HandleResolArea::MakeKeyIRFEffectiveAreaGraph(double eff, double zen, double off, unsigned int telcode) {

	TString key="Graph_Area_";
	key+=eff;
	key+="_";
	key+=zen;
	key+="_";
	key+=off;
	key+="_";
	key+=telcode;
  
	return key;  

}
/**
 * \brief Allows to determine which are the needed MC ranges to do the interpolations
 */
void START::HandleResolArea::GetNeededMcRange(std::vector<Band> const &BandArray,
											  double &mcmineff, double &mcmaxeff,
											  double &mcminoff, double &mcmaxoff,
											  double &mcminzen, double &mcmaxzen)
{

	std::vector<double> dist_off, dist_eff, dist_zen;
	for(std::vector<Band>::const_iterator band=BandArray.begin(); band!=BandArray.end(); ++band) {
		if(band->GetKeepBand()==0) continue;
		dist_off.push_back(band->GetOffset());
		dist_eff.push_back(band->GetEff());
		dist_zen.push_back(band->GetZenON());
	}

	// Find the min and max values of the bands parameters off, zen
	// and eff

	std::vector<double>::const_iterator it1, it2, it3, it4, it5, it6;
	it1 = min_element(dist_off.begin(), dist_off.end());
	it2 = max_element(dist_off.begin(), dist_off.end());
	double minoff = *it1;
	double maxoff = *it2;

	it3 = min_element(dist_eff.begin(), dist_eff.end());
	it4 = max_element(dist_eff.begin(), dist_eff.end());
	double mineff = *it3;
	double maxeff = *it4;

	it5 = min_element(dist_zen.begin(), dist_zen.end());
	it6 = max_element(dist_zen.begin(), dist_zen.end());
	double minzen = *it5;
	double maxzen = *it6;

	// Determine efficiencies pair
	if(fMc->GetEfficiency().size()>1) {
		for(unsigned int i(0); i<fMc->GetEfficiency().size()-1; i++) {
			if(mineff > fMc->GetEfficiency()[i] && mineff < fMc->GetEfficiency()[i+1])
				mcmineff = fMc->GetEfficiency()[i]; 
		}
		for(unsigned int i(0); i<fMc->GetEfficiency().size()-1; i++) {
			if(maxeff > fMc->GetEfficiency()[i] && maxeff < fMc->GetEfficiency()[i+1])
				mcmaxeff = fMc->GetEfficiency()[i+1]; 
		}
	}
	else {
		mcmineff = fMc->GetEfficiency()[0];
		mcmaxeff = fMc->GetEfficiency()[0];
	}

	// Determine offset pair
	if(fMc->GetOffset().size()>1) {
		for(unsigned int i(0); i<fMc->GetOffset().size()-1; i++) {
			if(minoff > fMc->GetOffset()[i] && minoff < fMc->GetOffset()[i+1])
				mcminoff = fMc->GetOffset()[i]; 
		}
		for(unsigned int i(0); i<fMc->GetOffset().size()-1; i++) {
			if(maxoff > fMc->GetOffset()[i] && maxoff < fMc->GetOffset()[i+1])
				mcmaxoff = fMc->GetOffset()[i+1]; 
		}
	}
	else {
		mcminoff = fMc->GetOffset()[0];
		mcmaxoff = fMc->GetOffset()[0];
	}

	// Determine zenith pair
	if(fMc->GetZenith().size()>1) {
		for(unsigned int i(0); i<fMc->GetZenith().size()-1; i++) {
			if(minzen > fMc->GetZenith()[i] && minzen > fMc->GetZenith()[i+1])
				mcminzen = fMc->GetZenith()[i]; 
		}
		for(unsigned int i(0); i<fMc->GetZenith().size()-1; i++) {
			if(maxzen > fMc->GetZenith()[i] && maxzen < fMc->GetZenith()[i+1])
				mcmaxzen = fMc->GetZenith()[i+1]; 
		}
	}
	else {
		mcminoff = fMc->GetZenith()[0];
		mcmaxoff = fMc->GetZenith()[0];
	}
  
}


/**
 * \brief Copy constructor
 *
 * JLK : keep? Ca va crasher si on fait une copy. (VIM : This should not crash)
 * The fBandArray is given by the user, and not stored in this function (just a pointer to the BandArray given by the user).
 * So this function modify directly what the user give, and does not store the value
 */
START::HandleResolArea::HandleResolArea(HandleResolArea const &StoreCopy)
	:TObject(StoreCopy),
	 fBandArray(0),
	 fForceUseDistribution(StoreCopy.fForceUseDistribution),
	 fConfigName(StoreCopy.fConfigName),
	 fEffMin(StoreCopy.fEffMin), fEffMax(StoreCopy.fEffMax),
	 fOffMin(StoreCopy.fOffMin), fOffMax(StoreCopy.fOffMax),
	 fZenMin(StoreCopy.fZenMin), fZenMax(StoreCopy.fZenMax),
	 fpairResolMeanSigma(StoreCopy.fpairResolMeanSigma),
	 fmapAreaValues(StoreCopy.fmapAreaValues),
	 fmapBiaisValues(StoreCopy.fmapBiaisValues),
	 fmapResolValues(StoreCopy.fmapResolValues),
	 fmapVectorsResol(StoreCopy.fmapVectorsResol),
	 fmapVectorsBiais(StoreCopy.fmapVectorsBiais),
	 fmapVectorsArea(StoreCopy.fmapVectorsArea),
	 fAreaMin(StoreCopy.fAreaMin),
	 fResolFactor(StoreCopy.fResolFactor)
{
	fSafeThresholdMethod = StoreCopy.fSafeThresholdMethod;
	fFractionMaxArea = StoreCopy.fFractionMaxArea;
  
	fBandArray=StoreCopy.fBandArray; // JLK : faux!!!!!!!!!!!!!!! si on fait une copy ca va cracher
	// VIM : No, the principle is to change the BandArray of the user, so this class does not create a new instance, just copy the pointer value.
	// JLK : si on detruit l'original le pointeur ne pointe plus sur rien donc si tu detruis la nouvelle ça plante. Beam.
	// VIM : I am not sure this is a good idea (maybe we should continue to distinguish, can be ambiguous to the user, and so on a bug generator...)
  
	for(std::map<TString, TH1F *>::const_iterator it=StoreCopy.fmapHistoResol.begin(); it!=StoreCopy.fmapHistoResol.end(); ++it) {
		fmapHistoResol[it->first] = new TH1F(*(it->second));
	}
  
	for(std::map<TString, TH1F *>::const_iterator it=StoreCopy.fmapHistoBiais.begin(); it!=StoreCopy.fmapHistoBiais.end(); ++it) {
		fmapHistoBiais[it->first] = new TH1F(*(it->second));
	}

	for(std::map<TString, TH1F *>::const_iterator it=StoreCopy.fmapHistoArea.begin(); it!=StoreCopy.fmapHistoArea.end(); ++it) {
		fmapHistoArea[it->first] = new TH1F(*(it->second));
	}
  
	for(std::map<TString, TH1F *>::const_iterator it=StoreCopy.fmapHistoDistrib.begin(); it!=StoreCopy.fmapHistoDistrib.end(); ++it) {
		if (it->second!=0) {
			fmapHistoDistrib[it->first] = new TH1F(*(it->second));
		}
		else {
			// VIM : Peut etre mettre un avertissement ici ?
		}
	}
}

/**
 * \brief Assignment operator
 *
 * JLK : keep?
 */
START::HandleResolArea &START::HandleResolArea::operator=(HandleResolArea const& StoreCopy)
{
	if(this != &StoreCopy) {
		fForceUseDistribution=StoreCopy.fForceUseDistribution;
		fConfigName=StoreCopy.fConfigName;
		fEffMin=StoreCopy.fEffMin;
		fEffMax=StoreCopy.fEffMax;
		fOffMin=StoreCopy.fOffMin;
		fOffMax=StoreCopy.fOffMax;
		fZenMin=StoreCopy.fZenMin;
		fZenMax=StoreCopy.fZenMax;
		fpairResolMeanSigma=StoreCopy.fpairResolMeanSigma;
		fmapAreaValues=StoreCopy.fmapAreaValues;
		fmapBiaisValues=StoreCopy.fmapBiaisValues;
		fmapResolValues=StoreCopy.fmapResolValues;
		fmapVectorsResol=StoreCopy.fmapVectorsResol;
		fmapVectorsBiais=StoreCopy.fmapVectorsBiais;
		fmapVectorsArea=StoreCopy.fmapVectorsArea;
		fBandArray=StoreCopy.fBandArray;
		fAreaMin=StoreCopy.fAreaMin;
		fResolFactor=StoreCopy.fResolFactor;
		fSafeThresholdMethod = StoreCopy.fSafeThresholdMethod;
		fFractionMaxArea = StoreCopy.fFractionMaxArea;
    
		for(std::map<TString, TH1F *>::iterator it=fmapHistoResol.begin(); it!=fmapHistoResol.end(); ++it) {
			if (it->second!=0) delete it->second;
		}
		fmapHistoResol.clear();
    
		for(std::map<TString, TH1F *>::iterator it=fmapHistoBiais.begin(); it!=fmapHistoBiais.end(); ++it) {
			if (it->second!=0) delete it->second;
		}
		fmapHistoBiais.clear();
    
		for(std::map<TString, TH1F *>::iterator it=fmapHistoArea.begin(); it!=fmapHistoArea.end(); ++it) {
			if (it->second!=0) delete it->second;
		}
		fmapHistoArea.clear();
    
		for(std::map<TString, TH1F *>::const_iterator it=fmapHistoDistrib.begin(); it!=fmapHistoDistrib.end(); ++it) {
			if (it->second!=0) delete it->second;
		}
		fmapHistoDistrib.clear();
    
		for(std::map<TString, TH1F *>::const_iterator it=StoreCopy.fmapHistoResol.begin(); it!=StoreCopy.fmapHistoResol.end(); ++it) {
			fmapHistoResol[it->first] = new TH1F(*(it->second));
		}
    
		for(std::map<TString, TH1F *>::const_iterator it=StoreCopy.fmapHistoBiais.begin(); it!=StoreCopy.fmapHistoBiais.end(); ++it) {
			fmapHistoBiais[it->first] = new TH1F(*(it->second));
		}
    
		for(std::map<TString, TH1F *>::const_iterator it=StoreCopy.fmapHistoArea.begin(); it!=StoreCopy.fmapHistoArea.end(); ++it) {
			fmapHistoArea[it->first] = new TH1F(*(it->second));
		}
    
		for(std::map<TString, TH1F *>::const_iterator it=fmapHistoDistrib.begin(); it!=fmapHistoDistrib.end(); ++it) {
			fmapHistoDistrib[it->first] = new TH1F(*(it->second));
		}
	}
  
	return (*this);
}
