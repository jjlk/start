#ifndef _START_FITSUtils_
#define _START_FITSUtils_

// STL
#include <string>

// CFITSIO
//#include <cfitsio/fitsio.h>
#include <fitsio.h>

// ROOT
#include <TObject.h>

namespace START {

  /**
   * \brief Class that contains function usefull for the FITS conversion of START data
   *
   **/

  class STARTFITSUtils : public TObject
  {
  public :
    virtual ~STARTFITSUtils(){};

    //static std::string SpectralHeader();
    static std::string GetSourceInfoHeader(std::string objectname="", Int_t TelPattern = 30, Double_t ra_obj=-999., Double_t dec_obj=-999., Double_t ra_point=-999., Double_t dec_point=-999., Double_t alt_point=-999., Double_t az_point=-999.);
    static std::string GetSpectralHeader();
    static std::string GetARFHeader();
    // static std::string GetRMFMatrixHeader();
    // static std::string GetRMFEBoundHeader();
    static Int_t UpdateHeaderWithString(fitsfile *fileptr, std::string inputheader, int *status);

#ifdef CONSTRUCT_LIBRARY
  public:
    ClassDef(START::STARTFITSUtils,1);
#endif


  };

}
#endif
