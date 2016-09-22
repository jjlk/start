#ifndef _RESIDUALSFACTORY_
#define _RESIDUALSFACTORY_

// STL
#include <vector>

// ROOT
#include <TObject.h>
#include <TString.h>

// START
namespace START {
	class Hypothesis;

	/**
	 * \brief Compute residuals, flux, associated errors and mean energy bin for the graphical representation.
	 * 
	 * Residuals are computed with the Rolke method and it gives errors at 1 and 3 sigmas.
	 *
	 * \author HAP-Fr team
	 */
	class ResidualsFactory : public TObject
	{
	public :
  
		ResidualsFactory(); // VIM : Root a besoin d'un constructeur par defaut (qui n'alloue pas dynamiquement de la memoire)
		ResidualsFactory(Hypothesis &hypo, bool verbose=false);   // JLK : Changement du constructeur, prend en plus les bandes
		ResidualsFactory(std::vector<Hypothesis*> &hypo, bool verbose=false);   // JLK : Changement du constructeur, prend en plus les bandes

		virtual ~ResidualsFactory();

		//void CopyExcessInBandsHypothesis(std::vector<Band> const &BandArray); ///< Copy the expected excess in the hypothesis'bands

		//void ComputeResidualsGaus(bool FittedOnOff);
		//void ComputeResidualsPoissonMC();
		void ComputeResidualsRolke();

	private :
  
		bool fverbose;

		std::vector<Hypothesis*> fHypothesisArray; ///< array of pointers which contains the adress of the hypothesis

		std::vector<double> fresiduals; ///< residuals
		// JLK ADD FOR ADA
		std::vector<double> fresiduals_on; ///< residuals
		std::vector<double> fresiduals_off; ///< residuals
		std::vector<double> fresidualssigmaplus; ///< residuals error + at 1 sigma
		std::vector<double> fresidualssigmaminus; ///< residuals error - at 1 sigma
		std::vector<double> fresiduals3sigmaplus; ///< residuals error + at 3 sigma
		std::vector<double> fresiduals3sigmaminus; ///< residuals error - at 3 sigma
		std::vector<double> fflux; ///< experimental flux
		std::vector<double> ffluxsigmaplus; ///< error + on experimental flux at 1 sigma
		std::vector<double> ffluxsigmaminus; ///< error - on experimental flux at 1 sigma
		std::vector<double> fflux3sigmaplus; ///< error + on experimental flux at 3 sigma
		std::vector<double> fflux3sigmaminus; ///< error - on experimental flux at 3 sigma
		std::vector<double> femean; ///< mean energy bin

#ifdef CONSTRUCT_LIBRARY
	public:
		ClassDef(START::ResidualsFactory,1);
#endif
	};
}
#endif
