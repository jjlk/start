/*
 *  Macro using a configuration file to analyse data
 */

// STL
#include <iostream>
#include <TMath.h>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>

// ROOT
#include <TString.h>
#include <TFile.h>
#include <TSystem.h>
// fitsio
//#include "cfitsio/fitsio.h"
//#include "fitsio.h"
#include <fitsio.h>

#ifndef __CINT__
// START
#include "Config.hh"
#include "BandsFactory.hh"
#include "Band.hh"
#include "DataSummary.hh"
#include "HandleResolArea.hh"
#include "debugging.hh"
#endif


// Beware: splitted runs are not properly dealt with need to check


void MakeHessRSPfiles(TString fileconfig) {
	std::cout << "fileconfig = " << fileconfig << std::endl;
	/// Creation of an object Config. It will contain all infos for the minimization & Co
	START::Config myConfig(fileconfig);

	/// Creation of a vector myData which will contain all the data in object Bands
	std::vector<START::Band> myData;

	// JLK Get IrfOpt
	std::string StrIrfAnalysis(myConfig.GetUserIrfOpt());
	START::STARTUtils::IrfOpt IrfType;
	if(StrIrfAnalysis=="HapFr_Hess_Stereo")
		IrfType = START::STARTUtils::HapFr_Hess_Stereo;
	else if(StrIrfAnalysis=="HapFr_Hess_Hybrid")
		IrfType = START::STARTUtils::HapFr_Hess_Hybrid;
	else if(StrIrfAnalysis=="HapFr_Hess_Mono")
		IrfType = START::STARTUtils::HapFr_Hess_Mono;
	else {
		std::cout << "WARNING - No known IrfOpt such as '" << StrIrfAnalysis <<  "' ==> EXIT !!" << std::endl;
		std::cout << "### Available:" << std::endl;
		std::cout << " - HapFr_Hess_Stereo" << std::endl;
		std::cout << " - HapFr_Hess_Hybrid" << std::endl;
		std::cout << " - HapFr_Hess_Mono" << std::endl;
		return;
	}


	/// Creation of an object BandsFactory which will copy the data in the object myData
	START::BandsFactory BandsFact(myConfig, IrfType);
	BandsFact.MakeBands(myData);

	/// Creation of an object HandleResolArea wich will store the effective areas and 
	/// resolutions useful to determine the energy threshold of the bands
	START::HandleResolArea *HandleIRF = 0;
	HandleIRF = new START::HandleResolArea(myConfig,IrfType);  
	HandleIRF->SetESafeThresholdCondition(START::HandleResolArea::AreaAndBiaisResolMethod,myConfig.GetUserAreaMin(),2.);
	HandleIRF->SetAndUpdateBands(myData); 

	// Here starts the new part 
	TString sourcename = myConfig.GetUserSourceName(); 
	double E0_mes= myConfig.GetUserERangeMin();
	double E1_mes= myConfig.GetUserERangeMax();
	int bins_m =  myConfig.GetUserERangeBinsNumber();
	int bins_v = 15*bins_m;   // Random value to start with

	double E0_vrai = E0_mes * exp(-1.5);
	double E1_vrai = E1_mes * exp(1.5);
	double d = (log10(E1_vrai) - log10(E0_vrai))/bins_v;  

	double Etrue=0;


	// 1.  Extract information on bands (emin, emax, on, off, alpha) and the interpolated data (aeff, mean, sigma) 
	for(std::vector<START::Band>::iterator binvec = myData.begin();binvec!=myData.end();binvec++){

		if(binvec->GetKeepBand()==0) continue;

		int n_rsp_bins = binvec->ebin.size();
		long* ev_on = new long[n_rsp_bins];
		long* ev_off = new long[n_rsp_bins];
		double* ratio = new double[n_rsp_bins];
		double* e_min = new double[n_rsp_bins];
		double* e_max = new double[n_rsp_bins];
		int* channel = new int[n_rsp_bins];
		int* group = new int[n_rsp_bins];
		// use finer (bins_v) binning for ARF and RMF
		double* ek_min = new double[bins_v];
		double* ek_max = new double[bins_v];
		double* arf = new double[bins_v];
		double rmf[n_rsp_bins];
		int index=0;

		for(std::vector<START::EnergyBin>::iterator it = binvec->ebin.begin();it!=binvec->ebin.end();it++){
     
			double emin = it->GetEmin(); 
			double emax = it->GetEmax();
			double on =   it->GetOn();
			double off =  it->GetOff();
			double alpha =it->GetAlpha();

			// fill the vectors
			ev_on[index]= (long) on; 
			ev_off[index]= (long) off; 
			ratio[index]=alpha; 
			e_min[index]=emin*1e9; // convert to keV !!
			e_max[index]=emax*1e9; // convert to keV !!
			channel[index]=index;  // removed +1 for sherpa compatibility
			group[index]=1;
			index++;
		}
  
		std::string spec_name("!run");
		std::stringstream runnum;
		runnum << binvec->GetNbRun();
		int status4=0;
		std::string rmf_name("!run_rmf");
		rmf_name += runnum.str();
		rmf_name += ".fits";
   
		fitsfile* fptr_rmf;
		std::string rmf_full_name=rmf_name;
		rmf_full_name += "(hess-rmf.tpl)";
		std::cout << rmf_full_name.c_str() << std::endl;
  
		if (fits_create_file(&fptr_rmf, rmf_full_name.c_str(), &status4))
			{
				std::cout << "Error code :" << status4 << std::endl;
				return;
			}
		if (fits_movnam_hdu(fptr_rmf, ANY_HDU, "MATRIX", 0, &status4))
			{
				std::cout << "Error code :" << status4 << std::endl;
				return;
			}
		/*
		  int tlmin4=1;
		  int tlmax4=n_rsp_bins;
		  int tlmin1=1;
		  int tlmax1=n_rsp_bins;
		*/
		int tlmin4=0;
		int tlmax4=n_rsp_bins-1;
		int tlmin1=0;
		int tlmax1=n_rsp_bins-1;

		double exposure = binvec->GetLiveTime()*3600.;//convert to seconds
		double backscale1 = 1.0;
		double backscale2 = 1/binvec->GetAlphaRun();
		double areascale =1.0;
		double Eth = binvec->GetEthMC()*1e9;

		double run_offset =  binvec->GetOffset();
		double run_efficiency =  binvec->GetEff();
		double run_Zenith_ON =  binvec->GetZenON();
		double run_Zenith_OFF =  binvec->GetZenOFF();
		int run_TelCode =  binvec->GetTelCode();
		double run_TSTART =  binvec->GetRunStartTime();
		double run_TSTOP =  binvec->GetRunEndTime();
    
		int Ntels;
		int Numtels;

		switch(run_TelCode)
			{
			case 30: Ntels=4;Numtels=1234;break;
			case 28: Ntels=3;Numtels=123;break;
			case 26: Ntels=3;Numtels=124;break;
			case 22: Ntels=3;Numtels=134;break;
			case 24: Ntels=2;Numtels=12;break;
			case 20: Ntels=2;Numtels=13;break;
			case 18: Ntels=2;Numtels=14;break;
			case 12: Ntels=2;Numtels=23;break;
			case 10: Ntels=2;Numtels=24;break;
			case  6: Ntels=2;Numtels=34;break;   
			default: Ntels=1;Numtels=0;break;
			}

		fits_write_key(fptr_rmf, TINT, "TLMIN4", &tlmin4 ,"min. channel value", &status4);
		fits_write_key(fptr_rmf, TINT, "TLMAX4", &tlmax4 ,"max. channel value", &status4);
		fits_update_key(fptr_rmf, TINT, "DETCHANS", &n_rsp_bins, "Number of spectral bins", &status4);    

		if (fits_movnam_hdu(fptr_rmf, ANY_HDU, "EBOUNDS", 0, &status4))
			{
				std::cout << "Error code :" << status4 << std::endl;
				return;
			}
		fits_update_key(fptr_rmf, TINT, "DETCHANS", &n_rsp_bins, "Number of spectral bins", &status4);
		fits_write_key(fptr_rmf, TINT, "TLMIN1", &tlmin1 ,"min. channel value", &status4);
		fits_write_key(fptr_rmf, TINT, "TLMAX1", &tlmax1 ,"max. channel value", &status4);
  
		// now loop over E_vrai for ARF and RMF
		double e0 = log10(E0_vrai);
    
		//    int test= binvec->fkeepband;
		for(int k=0;k<bins_v;k++){// 40 bins in E_vrai
			double ev_min = pow(10,e0+((k-1)*d));
			double ev_max = pow(10,e0+(k*d));
			Etrue = sqrt(ev_min*ev_max);
			double Aeff = binvec->GetInterpolatedArea(Etrue);
			arf[k]=Aeff*1e4; // convert to cm2 !!
			ek_min[k]=ev_min*1e9;
			ek_max[k]=ev_max*1e9;
			int detchan=0;
			int fchan=0;

			for(unsigned int m=0;m<binvec->ebin.size();m++){
	
				double emin = binvec->ebin.at(m).GetEmin(); 
				double emax = binvec->ebin.at(m).GetEmax();
				if((ev_min < emax* exp(1.5))&&(ev_max > emin* exp(-1.5)) && binvec->ebin.at(m).GetKeepBin()!=0 ){
					double integral=binvec->ebin.at(m).GetInterpolatedPartialIntegral(Etrue);
					rmf[detchan++]=integral;
					// Change!
					//	  if(detchan==1){fchan = m+1;}
					if(detchan==1){fchan = m;}
					//std::cout <<" rmf "<<integral<<" detchan "<<detchan<<std::endl;
	  
				}     
	     
			}// end of m-loop

			if (fits_movnam_hdu(fptr_rmf, ANY_HDU, "MATRIX", 0, &status4))
				{
					std::cout << "Error code :" << status4 << std::endl;
					return;
				}
			int n_grp=1;
			int f_chan[1];
			int n_chan[1];
			f_chan[0]=fchan;
			n_chan[0]=detchan;
			fits_write_col(fptr_rmf, TDOUBLE, 1, k+1, 1, 1, &ek_min[k], &status4);
			fits_write_col(fptr_rmf, TDOUBLE, 2, k+1, 1, 1, &ek_max[k], &status4);
			fits_write_col(fptr_rmf, TINT,    3, k+1, 1, 1, &n_grp, &status4);
			fits_write_col(fptr_rmf, TINT,    4, k+1, 1, 1, &f_chan, &status4);
			fits_write_col(fptr_rmf, TINT,    5, k+1, 1, 1, &n_chan, &status4);
			fits_write_col(fptr_rmf, TDOUBLE, 6, k+1, 1, detchan, &rmf, &status4);
		}// end loop over k

		if (fits_movnam_hdu(fptr_rmf, ANY_HDU, "EBOUNDS", 0, &status4))
			{
				std::cout << "Error code :" << status4 << std::endl;
				return;
			}
  
		fits_write_col(fptr_rmf, TINT,    1, 1, 1, n_rsp_bins, channel, &status4); 
		fits_write_col(fptr_rmf, TDOUBLE, 2, 1, 1, n_rsp_bins, e_min, &status4);
		fits_write_col(fptr_rmf, TDOUBLE, 3, 1, 1, n_rsp_bins, e_max, &status4);
		std::cout <<"RMF Matrix written  " <<status4<<std::endl;
           
		fits_close_file(fptr_rmf, &status4);
      
		// 2.1 Build the SPECTRUM file + fill it
		int status1=0;
		spec_name += runnum.str();
		spec_name += ".pha";

		int status2=0;
		std::string bkg_name("!run_bkg");
		bkg_name += runnum.str();
		bkg_name += ".pha";

		int status3=0;
		std::string arf_name("!run_arf");
		arf_name += runnum.str();
		arf_name += ".fits";
    
		fitsfile* fptr_spec;
		if (fits_create_file(&fptr_spec, spec_name.c_str(), &status1))
			{
				std::cout << "Error code :" << status1 << std::endl;
				return;
			}

		long nrows, pcount;
		char* ttype[]={"CHANNEL","COUNTS","GROUPING"};
		char* tform[]={"I","J","I"};
		char* tunits[]={"","count",""};
		nrows=n_rsp_bins;
		pcount=0;

		fits_create_tbl(fptr_spec, BINARY_TBL, nrows, 3, ttype, tform, tunits, "SPECTRUM", &status1);
		std::cout << "Table creation SPEC:" << status1 << std::endl;
		fits_write_col(fptr_spec, TINT, 1, 1, 1, nrows,channel,&status1);
		fits_write_col(fptr_spec, TLONG,2, 1, 1, nrows,ev_on,  &status1);
		fits_write_col(fptr_spec, TINT, 3, 1, 1, nrows,group,  &status1);

		fits_write_key_template(fptr_spec,"hess-spec.tpl",&status1);
		fits_update_key(fptr_spec, TINT, "TLMIN1", &tlmin1 ,"min. channel value", &status1);
		fits_update_key(fptr_spec, TINT, "TLMAX1", &tlmax1 ,"max. channel value", &status1);
		fits_update_key(fptr_spec, TINT, "DETCHANS", &n_rsp_bins, "Number of spectral bins", &status1);
		fits_write_key(fptr_spec, TDOUBLE, "EXPOSURE", &exposure ,"exposure", &status1);
		fits_update_key(fptr_spec, TDOUBLE, "BACKSCAL", &backscale1 ," ", &status1);
		fits_write_key(fptr_spec, TDOUBLE, "AREASCAL", &areascale," ", &status1);
    
		char* bkg_file = (char*) bkg_name.substr(1).c_str();
		fits_write_key(fptr_spec, TSTRING, "BACKFILE", bkg_file,"/background  ", &status1);
		char* rmf_file = (char*) rmf_name.substr(1).c_str();
		fits_write_key(fptr_spec, TSTRING, "RESPFILE", rmf_file,"/redistribution  ", &status1);
		char* arf_file = (char*) arf_name.substr(1).c_str();    
		fits_write_key(fptr_spec, TSTRING, "ANCRFILE", arf_file,"/ancilliary response ", &status1);
		fits_write_key(fptr_spec, TDOUBLE, "ETH", &Eth,"Energy threshold of the run", &status1);
		fits_write_key(fptr_spec, TSTRING, "OBJECT", (char*) sourcename.Data(),"", &status1);

		fits_write_key(fptr_spec, TDOUBLE, "TSTART",&run_TSTART ,"Start time of the run", &status1);
		fits_write_key(fptr_spec, TDOUBLE, "TSTOP",&run_TSTOP ,"Stop time of the run", &status1);
		fits_write_key(fptr_spec, TINT, "TELCODE",&run_TelCode ,"Run TelCode", &status1);
		fits_update_key(fptr_spec, TINT, "N_TELS",&Ntels ,"number of telescopes participating", &status1);
		fits_write_key(fptr_spec, TINT, "NUM_TELS",&Numtels ,"Participating telescopes numbers", &status1);
		fits_write_key(fptr_spec, TDOUBLE, "ZENITH",&run_Zenith_ON ,"Mean zenith angle", &status1);
		fits_write_key(fptr_spec, TDOUBLE, "OFFSET",&run_offset ,"Run offset", &status1);
		fits_write_key(fptr_spec, TDOUBLE, "EFFICIEN",&run_efficiency ,"Run efficiency", &status1);

		std::cout << "Table SPEC written:" <<status1 << std::endl;

		fits_close_file(fptr_spec,&status1);

		// 2.2 Build the BACKGROUND file + fill it
		nrows=n_rsp_bins;

		fitsfile* fptr_bkg;
		if (fits_create_file(&fptr_bkg, bkg_name.c_str(), &status2))
			{
				std::cout << "Error code :" << status2 << std::endl;
				return;
			}

		fits_create_tbl(fptr_bkg, BINARY_TBL, nrows, 3, ttype, tform, tunits, "SPECTRUM", &status2);
		std::cout << "Table creation BKG:" << status2 << std::endl;
		fits_write_col(fptr_bkg, TINT, 1, 1, 1, nrows,channel,&status2);
		fits_write_col(fptr_bkg, TLONG,2, 1, 1, nrows,ev_off ,&status2);
		fits_write_col(fptr_bkg, TINT, 3, 1, 1, nrows,group  ,&status2);

		fits_write_key_template(fptr_bkg,"hess-spec.tpl",&status2);
		fits_update_key(fptr_bkg, TINT, "TLMIN1", &tlmin1 ,"min. channel value", &status2);
		fits_update_key(fptr_bkg, TINT, "TLMAX1", &tlmax1 ,"max. channel value", &status2);
		fits_update_key(fptr_bkg, TINT, "DETCHANS", &n_rsp_bins, "Number of spectral bins", &status2);
		fits_write_key(fptr_bkg, TDOUBLE, "EXPOSURE", &exposure ,"exposure", &status2);
		fits_update_key(fptr_bkg, TDOUBLE, "BACKSCAL", &backscale2 ,"", &status2);
		fits_write_key(fptr_bkg, TDOUBLE, "AREASCAL", &areascale," ", &status2);
		std::cout << "Table BGK written"<<status2 << std::endl;

		fits_close_file(fptr_bkg,&status2);

		// 2.3 Build the ARF file + fill it
		char* ttype1[]={"ENERG_LO","ENERG_HI","SPECRESP"};
		char* tform1[]={"1J","1J","1J"};
		char* tunits1[]={"keV","keV","cm^2"};
		nrows=bins_v;

		fitsfile* fptr_arf;
		if (fits_create_file(&fptr_arf, arf_name.c_str(), &status3))
			{
				std::cout << "Error code :" << status3 << std::endl;
				return;
			}
		fits_create_tbl(fptr_arf, BINARY_TBL, nrows, 3, ttype1, tform1, tunits1, "SPECRESP", &status3);
		std::cout << "Table creation ARF:" << status3 << std::endl;

		fits_write_key_template(fptr_arf,"hess-arf.tpl",&status3);
		std::cout << status3 << std::endl;
		fits_write_col(fptr_arf, TDOUBLE, 1, 1, 1, nrows, ek_min, &status3); 	 
		fits_write_col(fptr_arf, TDOUBLE, 2, 1, 1, nrows, ek_max, &status3); 	 
		fits_write_col(fptr_arf, TDOUBLE, 3, 1, 1, nrows, arf,    &status3);

		fits_close_file(fptr_arf, &status3);

	} 
}

#ifndef __CINT__

int main( int argc, const char* argv[] )
{

	int code=0;

	std::string help_string("Usage: ");
	help_string += argv[0];
	help_string += "\t ConfigFileName\n";
	help_string += "Environment variable HESS_TPL must be initialized";
  
	if (argc<2)
		{
			std::cerr << "Too few arguments. Quit!" << std::endl;
			std::cout << help_string << std::endl;
			code=-1;
		}
	else if (argc > 2)
		{
			std::cerr << "Too many arguments. Quit!" << std::endl;
			std::cout << help_string << std::endl;
			code=-1;
		}
	else if (std::string(argv[1]) == "--help")  
		{
			std::cout << help_string << std::endl;
			code = -1;
		}
	else if (std::string(argv[1]) == "-h")  
		{
			std::cout << help_string << std::endl;
			code = -1;
		}
	else if (std::string(argv[1]) == "-?")  
		{
			std::cout << help_string << std::endl;
			code = -1;
		}
	else
		{
			system("cp $HESS_TPL/hess-arf.tpl .");
			system("cp $HESS_TPL/hess-rmf.tpl .");
			system("cp $HESS_TPL/hess-rmf-matrix.tpl .");     
			system("cp $HESS_TPL/hess-rmf-ebounds.tpl ."); 
			system("cp $HESS_TPL/hess-spec.tpl .");
			MakeHessRSPfiles(std::string(argv[1]));
			system("rm -f hess-arf.tpl hess-rmf.tpl hess-rmf-matrix.tpl hess-rmf-ebounds.tpl hess-spec.tpl");
		}

	return 0;
}


#endif
