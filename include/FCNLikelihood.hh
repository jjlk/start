#ifndef _LIKELIHOOD_
#define _LIKELIHOOD_

// STL
#include <iostream>
#include <vector>
#include <map>
#include <utility>

// ROOT
#include <TString.h>

// MINUIT2
#include "Minuit2/FCNBase.h"

// START
#include "ComputeResults.hh"
#include "Hypothesis.hh"
#include "Band.hh"
#include "STARTUtils.hh"
#include "Config.hh"

namespace START {
	/**
	 * \brief Internal class used by MinimizeFactory for minimization.
	 *
	 * This class is used by MINUIT2 via MinimizeFactory to compute the likelihood.
	 *
	 * \author HAP-Fr team
	 */
	class FCNLikelihood : public ROOT::Minuit2::FCNBase
	{

	public :

		FCNLikelihood() {}; // Needed by ROOT
	        FCNLikelihood(const std::vector<Band> &SelectedBands, Hypothesis &hypo, Config &Configuration, bool verbose=false);
		~FCNLikelihood();
		virtual double Up() const {return fErrorDef;} ///< We define the errordef here (logL=logL_min+1/2)
		virtual double operator()(const std::vector<double> &) const; ///< We compute the fcn here
		virtual void SetErrorDef(const double errordef) {fErrorDef=errordef;};
		
		void SetBandsGroupingType(const START::STARTUtils::BandsGroupingTypeForMinimization Type) 
		{fBandsGroupingType = Type;};

	  const std::vector<START::Band> GetStackBandFromAllBands(const std::vector<double> &par, float zmin, float zmax) const;
	  const std::vector<START::Band> GetStackZenBandsFromAllBands(const std::vector<double> &par) const;

	private :

		double fErrorDef; ///< Error defintion of our likelihood (=0.5)

		bool fverbose;
	  
	        mutable bool fFirstCallDone;

		std::vector<Band> fBandArray; ///< Band's vector
  
		Hypothesis *fHypothesis; ///< Spectrum hypothesis
   
		ComputeResults *fCompRes; ///< ComputeResult to save time !

        // Type of grouping
		START::STARTUtils::BandsGroupingTypeForMinimization fBandsGroupingType;
	  //ADA Give access to the configuration file 
	        Config *fConfig;
	  
#ifdef CONSTRUCT_LIBRARY
	public:
		ClassDef(START::FCNLikelihood,1);
#endif  
	};

}
#endif
