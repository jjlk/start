// STL
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>

// ROOT
#include <TMath.h>

// CFITSIO
//#include <cfitsio/fitsio.h>
#include <fitsio.h>

// START
#include <STARTUtils.hh>
#include <STARTFITSUtils.hh>

#define DEBUG 0
#include <debugging.hh>

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::STARTFITSUtils);
#endif

std::string START::STARTFITSUtils::GetSourceInfoHeader(std::string objectname, Int_t TelPattern,
						       Double_t ra_obj, Double_t dec_obj, Double_t ra_point, 
						       Double_t dec_point, Double_t alt_point, Double_t az_point) {
  
  std::ostringstream oss_header;
  std::string endline="\n";

  std::vector<Int_t> tellist = START::STARTUtils::GetTelescopeListFromPattern(TelPattern);
  
  oss_header << "OBJECT  = '" << objectname     << "' /  observed object"                      << endline;
  oss_header << "RA_OBJ  =  " << ra_obj         << "  /  RA of target object (deg)"            << endline;
  oss_header << "DEC_OBJ =  " << dec_obj        << "  /  Dec of target object (deg)"           << endline;
  oss_header << "RA_PNT  =  " << ra_point       << "  /  nominal pointing position RA (deg)"   << endline;
  oss_header << "DEC_PNT =  " << dec_point      << "  /  nominal pointing position Dec (deg)"  << endline;
  oss_header << "ALT_PNT =  " << alt_point      << "  / mean altitude of run (deg)"            << endline;
  oss_header << "AZ_PNT  =  " << az_point       << "  /  mean azimuth of run (deg)"            << endline;
  oss_header << "OBS_MODE= 'wobble  '           /  observation mode (wobble,scan,on,off)"      << endline;
  oss_header << "N_TELS  =  " << tellist.size() << "  /  number of telescopes participating"   << endline;
  oss_header << "TEL_CODE=  " << TelPattern     << "  /  telescope pattern of the obseration"  << endline;

  std::string header = oss_header.str();
  return header;
}

std::string START::STARTFITSUtils::GetSpectralHeader() {

  std::string endline="\n";
  std::ostringstream oss_header;
  oss_header <<  "PCOUNT  =                    0 / number of group parameters" << endline;
  oss_header <<  "GCOUNT  =                    1 / number of groups" << endline;
  oss_header <<  "TLMIN1  =                    0 / Lowest legal channel number" << endline;
  oss_header <<  "TLMAX1  =                   59 / Highest legal channel number" << endline;
  oss_header <<  "POISSERR= T                  / Poissonian errors to be assumed" << endline;
  oss_header <<  "STAT_ERR=                    0 / no statistical error specified" << endline;
  oss_header <<  "SYS_ERR =                    0 / no systematic error specified" << endline;
  oss_header <<  "QUALITY =                    0 / no data quality information specified" << endline;
  oss_header <<  "DETCHANS=                   60 / Total No. of Detector Channel available" << endline;
  oss_header <<  "CORRSCAL=                  1.0 / correlation scale factor" << endline;
  oss_header <<  "BACKSCAL=                  1.0 / background scale factor" << endline;
  oss_header <<  "BACKFILE= 'none    '           / background FITS file for" << endline;
  oss_header <<  "CORRFILE= 'none    '           /  correlation FITS file for" << endline;
  oss_header <<  "CHANTYPE= 'PHA      '           /  Channels assigned by detector electronics" << endline;
  oss_header <<  "HDUCLASS= 'OGIP    '           /  format conforms to OGIP standard" << endline;
  oss_header <<  "HDUCLAS1= 'SPECTRUM'           / PHA dataset (OGIP memo OGIP-92-007)" << endline;
  oss_header <<  "HDUVERS = '1.2.1   '           / Version of format (OGIP memo OGIP-92-007)" << endline;
  oss_header <<  "TELESCOP= 'HESS     '           / Telescope (mission) name" << endline;
  oss_header <<  "INSTRUME= 'HESS     '           /  Instrument name" << endline;
  oss_header <<  "FILTER  = 'none    '           / Instrument filter in use" << endline;

  std::string header = oss_header.str();  
  return header;
}


std::string START::STARTFITSUtils::GetARFHeader() {

  std::string endline="\n";
  std::ostringstream oss_header;

  oss_header <<  "SIMPLE T" << endline;
  oss_header <<  "BITPIX = 16" << endline;
  oss_header <<  "NAXIS = 0" << endline;
  oss_header <<  "EXTEND = T" << endline;
  oss_header <<  "XTENSION= 'BINTABLE'           / binary table extension" << endline;
  oss_header <<  "EXTNAME = 'SPECRESP'           / The name of this table" << endline;
  oss_header <<  "BITPIX  =                    8 / 8-bit bytes" << endline;
  oss_header <<  "NAXIS   =                    2 / 2-dimensional binary table" << endline;
  oss_header <<  "PCOUNT  =                    0 / size of special data area" << endline;
  oss_header <<  "GCOUNT  =                    1 / one data group (required keyword)" << endline;
  oss_header <<  "TFIELDS =                    3 / number of fields in each row" << endline;
  oss_header <<  "TTYPE1  = 'ENERG_LO'           / No comment" << endline;
  oss_header <<  "TFORM1  = 'E       '           / data format of field: 4-byte REAL" << endline;
  oss_header <<  "TUNIT1  = 'keV     '           / physical unit of field" << endline;
  oss_header <<  "TTYPE2  = 'ENERG_HI'           / No comment" << endline;
  oss_header <<  "TFORM2  = 'E       '           / data format of field: 4-byte REAL" << endline;
  oss_header <<  "TUNIT2  = 'keV     '           / physical unit of field" << endline;
  oss_header <<  "TTYPE3  = 'SPECRESP'           / No comment" << endline;
  oss_header <<  "TFORM3  = 'E       '           / data format of field: 4-byte REAL" << endline;
  oss_header <<  "TUNIT3  = 'cm2     '           / physical unit of field" << endline;
  oss_header <<  "TELESCOP= 'HESS    '           / mission/satellite name" << endline;
  oss_header <<  "INSTRUME= 'HESS    '           / required by XSPEC" << endline;
  oss_header <<  "FILTER  = '        '           / filter in use (X-ray relic)" << endline;
  oss_header <<  "HDUCLASS= 'OGIP    '           / Format conforms to OGIP standard" << endline;
  oss_header <<  "HDUCLAS1= 'RESPONSE'           / dataset relates to spectral response" << endline;
  oss_header <<  "HDUCLAS2= 'SPECRESP'           / dataset is a spectral response matrix" << endline;
  oss_header <<  "HDUVERS1= '1.3.0   '           / Obsolete - included for backwards compatibility" << endline;
  oss_header <<  "HDUVERS2= '1.3.0   '           / Obsolete - included for backwards compatibility" << endline;
  oss_header <<  "HDUVERS = '1.3.0   '           / the version of the HDUCLAS2 format in use" << endline;
  
  std::string header = oss_header.str();
  return header;
}

/**
 * \brief Implementation the cfitsio function fits_write_key_template to read a string
 * This is needed when you want to construct your own template in C++
 * VIM : Maybe an expert in cfitsio would have said that this function already exist, but I haven't found it.
 */
Int_t START::STARTFITSUtils::UpdateHeaderWithString(fitsfile *fileptr, std::string inputheader, int *status) {
  
  std::istringstream iss_inputheader(inputheader);

  std::string line;
  while( std::getline(iss_inputheader,line) ) {
    char card[FLEN_CARD];
    int keytype;
    //int status=0;
    if ( fits_parse_template(const_cast<char *>(line.c_str()),card,&keytype,status) > 0 ) {
      std::cout << YELLOWCOLOR << "Problem in parsing the line : " << line << " --> Skip it (\"le kangourou\") !" << RESETCOLOR << std::endl;
      //break;
      continue;
    }

    std::string card_str = std::string(card);
    std::string keyname = card_str.substr(0,std::min(8,FLEN_CARD)); // Compare with flen_card for safety reason. The 8 is from cfistio

    if (keytype == -2) {  // rename the card 
      if (card_str.size()<48) {
	std::cout << REDCOLOR << "size of card = " << card_str.size() << " (CARD = " << card_str << ") is less than what expected from cfitsio --> I can't rename the card !" << RESETCOLOR << std::endl;
      }
      else {
	std::string newname = card_str.substr(40,40+8);
	//ffmnam(fileptr, keyname, newname, status); 
#if (defined CFITSIO_MAJOR)
	fits_modify_name(fileptr,keyname.c_str(),newname.c_str(),status);     
#else
	fits_modify_name(fileptr,const_cast<char*>(keyname.c_str()),const_cast<char*>(newname.c_str()),status);     
#endif
      }
    }
    else if (keytype == -1) {     // delete the card
      //ffdkey(fileptr, keyname, status);	 
#if (defined CFITSIO_MAJOR) 
      fits_delete_key(fileptr, keyname.c_str(), status);	 
#else 
      fits_delete_key(fileptr, const_cast<char*>(keyname.c_str()), status);	 
#endif
    }
    else if (keytype == 0) {  // update the card
      //ffucrd(fileptr, keyname, card, status);
#if (defined CFITSIO_MAJOR && CFITSIO_MAJOR>=3 && CFITSIO_MINOR>=33)
	fits_update_card(fileptr, keyname.c_str(), card_str.c_str(), status);
#else
	fits_update_card(fileptr, const_cast<char*>(keyname.c_str()), const_cast<char*>(card_str.c_str()), status);    
#endif
    }
    else if (keytype == 1)  {   // append the card
      //ffprec(fileptr, card, status);
      fits_write_record(fileptr, card_str.c_str(), status);
    }
    else {
      std::cout << YELLOWCOLOR << __FUNCTION__ << "> I don't know what to do with this keytype = " << keytype << "... " << RESETCOLOR << std::endl;
    }
    
  }
  //std::getline(
  
  return *status;
}

