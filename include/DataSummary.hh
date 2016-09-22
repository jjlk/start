#ifndef _DATASUMMARY_
#define _DATASUMMARY_

// STL
#include <iostream>
#include <utility>

// ROOT
#include <TH1.h>
#include <TGraph.h>
#include <TObject.h>
#include <TString.h>
#include <TCanvas.h>

// START
namespace START {
  
  class Band;


  /**
   * \brief This class contains global informations on the band and can plot and print infos about bands.
   *
   * It contains Distributions of zenith angle, offset, 
   * efficiency and energy of the bands.
   * It also contains members used to determine which MC DSTs
   * have to be used for the minimization
   *
   * \author HAP-Fr team
   */
  class DataSummary : public TObject 
  {
  public :
    DataSummary(); // default constructor
    DataSummary(std::vector<Band> const &BandArray);
    ~DataSummary();

    DataSummary(DataSummary const &DataSumCopy); // Copy constructor
    DataSummary &operator=(DataSummary const& DataSumCopy); // 

    // Build the MC min and max values of eff, zen and off 
    void MakePairs(std::vector<Band> const &BandArray);
    void DetermineTeltype(std::vector<Band> const &BandArray);

    void PrintDataInfo(const std::vector<Band> &BandArray) const;
    void PrintDataInfoInFile(const std::vector<Band> &BandArray,std::ostream &os) const;
    void PrintSummaryBand(const std::vector<Band> &BandArray, bool UseEThreshold=true) const;

    void InitAndFillSummaryHistograms(std::vector<Band> const &BandArray);
    void DrawSummaryHistograms(std::vector<Band> const &BandArray);

    std::pair<double,double> GetPairEfficiency() const {return fPairEfficiency;}
    std::pair<double,double> GetPairZenith() const {return fPairZenith;}
    std::pair<double,double> GetPairOffset() const {return fPairOffset;}

    std::vector<TString> GetWhichTelType() const {return fWhichTelType;}

  private :

    void CleanROOTObjetcs();
    void CleanTH1();
    void CleanGraph();
    void CleanCanvas();

    // Contain values min and max of zen, off, eff
    // computed from the BandArray
    std::pair<double,double> fPairEfficiency;
    std::pair<double,double> fPairZenith;
    std::pair<double,double> fPairOffset;

    // Contain values of zen, off, eff include in Bands
    std::vector<double> fEfficiency;
    std::vector<double> fZenith;
    std::vector<double> fOffset;

    // Tell if it exists run with three telecopes
    // if all bands are with four telescopes then fWichTelType = 4
    // else fWichTelType = 3
    std::vector<TString> fWhichTelType;

    TCanvas *fCanvasEfficiency;
    TCanvas *fCanvasOffset;
    TCanvas *fCanvasZenithON;
    TCanvas *fCanvasZenithOFF;
    TCanvas *fCanvasEnergyON;
    TCanvas *fCanvasEnergyOFF;
    TCanvas *fCanvasLiveTime;
    TCanvas *fCanvasAlpha;
    TCanvas *fCanvasSignificanceVsLivetime;
    TCanvas *fCanvasOffEventsVsLivetime;
    TCanvas *fCanvasOnEventsVsLivetime;
    TCanvas *fCanvasExcessVsLivetime;
    TCanvas *fCanvasOffsetVsLivetime;
    TCanvas *fCanvasTimeMJD;

    TH1D *fHistoEfficiency;
    TH1D *fHistoOffset;
    TH1D *fHistoZenithON;
    TH1D *fHistoZenithOFF;
    TH1D *fHistoEnergyON;
    TH1D *fHistoEnergyOFF;
    TH1D *fHistoLiveTime;
    TH1D *fHistoAlpha;
    TGraph *fGraphSignificanceVsLivetime;
    TGraph *fGraphOffEventsVsLivetime;
    TGraph *fGraphOnEventsVsLivetime;
    TGraph *fGraphExcessVsLivetime;
    TGraph *fGraphOffsetVsLivetime;
    TH1D *fHistoTimeMJD;

#ifdef CONSTRUCT_LIBRARY
  public:
    ClassDef(START::DataSummary,1);
#endif
  };
}
#endif
