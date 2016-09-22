#ifndef _MULTIWAVELENGTHFACTORY_
#define _MULTIWAVELENGTHFACTORY_

// STL
#include <vector>
#include <string>
#include <map>
#include <utility>

// ROOT
#include <TNamed.h>
#include <Rtypes.h>
#include <TPolyLine.h>
#include <TF1.h>

// START
#include "Residuals.hh"

namespace START {
  /**
   * 
   * \brief Contains data for a multi-wavelength plot
   *
   * \author HAP-Fr team
   */
  class MultiWaveLengthFactory : public TNamed
  {

  public:

    /**
     * \brief Energy units
     */
    typedef enum {Hz, keV, MeV, GeV, TeV} EnergyUnits;

    /**
     * \brief Flux units
     * \param vFv erg.cm-2.s-1
     */
    typedef enum {vFv,ESquareFlux_TeV,Differential_TeV,Differential_GeV,Differential_MeV,Differential_keV} FluxUnits;

    MultiWaveLengthFactory(std::string name, FluxUnits FUnits=ESquareFlux_TeV, EnergyUnits EUnits=TeV);
    ~MultiWaveLengthFactory();
    MultiWaveLengthFactory(MultiWaveLengthFactory const &MwlCopy); // Copy constructor
    MultiWaveLengthFactory &operator=(MultiWaveLengthFactory const &MwlCopy); // assignment operator

    void AddPoints(std::string name, Residuals &Res,
		   Color_t MarkerColor=kBlack, Style_t MarkerStyle=21, Size_t MarkerSize=1,
		   Color_t LineColor=kBlack, Width_t LineWidth=1);
    void AddPoints(std::string name, std::vector<double> x, std::vector<double> y,
		   Color_t MarkerColor=kBlack, Style_t MarkerStyle=21, Size_t MarkerSize=1,
		   Color_t LineColor=kBlack, Width_t LineWidth=1);
    void AddPoints(std::string name, std::vector<double> x, std::vector<double> y, 
		   std::vector<double> xerr, std::vector<double> yerr,
		   Color_t MarkerColor=kBlack, Style_t MarkerStyle=21, Size_t MarkerSize=1,
		   Color_t LineColor=kBlack, Width_t LineWidth=1);
    void AddPoints(std::string name, std::vector<double> x, std::vector<double> y, 
		   std::vector<double> xminus, std::vector<double> xplus,
		   std::vector<double> y1minus, std::vector<double> y1plus,
		   Color_t MarkerColor=kBlack, Style_t MarkerStyle=21, Size_t MarkerSize=1,
		   Color_t LineColor=kBlack, Width_t LineWidth=1);
    void AddUpperLimits(std::string name, std::vector<double> x, std::vector<double> yhigh,
			Color_t LineColor=kBlack, Width_t LineWidth=1);
    void AddButterfly(std::string name, TPolyLine &Butterfly,
		      Color_t FillColor=0, Style_t FillStyle=0,
		      Color_t LineColor=kBlack, Width_t LineWidth=1);
    void AddButterfly(std::string name, std::vector<double> x, std::vector<double> y,
		      Color_t FillColor=0, Style_t FillStyle=0,
		      Color_t LineColor=kBlack, Width_t LineWidth=1);
    void AddButterflyJLK(std::string name, TPolyLine &Butterfly,
			 Color_t FillColor=0, Style_t FillStyle=0,
			 Color_t LineColor=kBlack, Width_t LineWidth=1);
    void AddTF1(std::string name, TF1 &f, Color_t LineColor=kBlack, Style_t LineStyle=1, Width_t LineWidth=1);

    std::map<std::string,TPolyLine*> GetMapButterfly() {return fMapButterfly;};
    std::map<std::string,Residuals> GetMapResiduals() {return fMapResiduals;};
    std::map<std::string,TF1*> GetMapTF1() {return fMapTF1;}; 

    FluxUnits GetFluxUnits() const {return fFluxUnits;};
    EnergyUnits GetEnergyUnits() const {return fEnergyUnits;};

    std::pair<double,double> GetFluxRange();
    std::pair<double,double> GetEnergyRange();
  
  private:

    EnergyUnits fEnergyUnits;
    FluxUnits fFluxUnits;

    std::map<std::string,TPolyLine*> fMapButterfly;
    std::map<std::string,Residuals> fMapResiduals;
    std::map<std::string,TF1*> fMapTF1;

#ifdef CONSTRUCT_LIBRARY
  public:
    ClassDef(START::MultiWaveLengthFactory,1);
#endif

  };
}
#endif
