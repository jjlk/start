#ifndef _COMPUTERES_
#define _COMPUTERES_

#include <iostream>
#include <vector>
#include <map>
#include <utility>

#include <TString.h>
#include <Math/AllIntegrationTypes.h>

#include "Band.hh"
#include "Hypothesis.hh"

namespace START {
  /**
   * \brief Internal class for intermediate computations such as integrals.
   * 
   * This class is used to compute the expected excesses which is a double integral on
   * true energy \f$E\f$ and reco energy \f$ \widetilde{E} \f$:
   * \f[
   * S = \int_{\widetilde{E}_{min}}^{\widetilde{E}_{max}} 
   * d\widetilde{E} \int_{0}^{\infty} R(\widetilde{E},E) A(E) \phi(E) dE
   * \f]
   * where \f$ R(\widetilde{E},E) \f$ is the resolution and \f$ A(E) \f$ is the effective
   * area
   *
   * Some of the method of this class seems weird or redundent with other part of the code, but they are workaround to use ROOT Integrator
   *
   * \author HAP-Fr team
   */

  class ComputeResults : public TObject
  {

  public :

    ComputeResults(const std::vector<Band> &SelectedBands, Hypothesis &hypo); // use in the FCN
    ComputeResults(const std::vector<Band> &SelectedBands); // Use for pre-calculation of integrated resolution

    ~ComputeResults();

    double IntegrantExpectedExcess(double x); // expected excess integrant 
 
    // Compute the expected excess for one band and one energy bin
    double FunctionExpectedExcess(Band *iband, EnergyBin *iem, 
																	const std::vector<double> &paramfit); 

    int MakeVectorPartialIntegral(std::vector<Band> &SelectedBands); // Make and copy in energybin a 
    //vector which contains the integral of the resolution on the measured energy for each MC energy

    double PartialIntegral(EnergyBin *ebin, double Etrue);  // find the partial integral for Etrue

    double IntegrantFunctionResolution(double x);  // 
  
  

  private :
    double FunctionEffectiveArea(double Etrue, Band *band);
    double FunctionResolution(double Etrue, double Ereco, Band *band);

    std::vector<Band> fBandArray; ///< Band's vector

    Hypothesis *fHypothesis; ///< Hypothesis

    double fabsoluteprecision; ///< absolute precision for integration
    double frelativeprecision; ///< relative precision for integration
  
    ROOT::Math::IntegrationOneDim::Type fIntegrationType;

    Band *myband_tmp; //! ///< current band (used for the computation of the expected excess)
    EnergyBin *mynrj_tmp; //! ///< current bin (used for the computation of the expected excess)
    double *nrjreco_tmp; //! ///< current reco energy use tom compute the integral over reco energy 

    std::vector<double> fparamfromminimizer; ///< used for the bidouille, contains paramaters from Minimizer

#ifdef CONSTRUCT_LIBRARY
  public:
    ClassDef(START::ComputeResults,1);
#endif
  };
}

#endif
