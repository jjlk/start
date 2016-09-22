#ifndef _MONTECARLO_
#define _MONTECARLO_

// STL
#include <vector>
#include <utility>

// ROOT
#include <TObject.h>
#include <TString.h>

// START
#include "STARTUtils.hh"

namespace START {
	/**
	 * \brief Class used to store Montecarlo's parameters.
	 *
	 * Short class where the parameters of the simulation
	 * (efficacity, offset, and zenith angle telescope type)
	 * are stored.
	 *
	 * \author HAP-Fr
	 */
	class MonteCarlo : public TObject
	{
	public :

		MonteCarlo(START::STARTUtils::IrfOpt Type);
		~MonteCarlo();

		std::vector<double> GetZenith() const {return fZenithMC;}
		std::vector<double> GetEfficiency() const {return fEfficiencyMC;}
		std::vector<double> GetOffset() const {return fOffsetMC;}
		std::vector<TString> GetTelType() const {return fTelTypeMC;}
		std::vector<double> GetEnergy() const {return fEnergyMC;}

		TString GetStringFromTelcode(unsigned int telcode) const;

		std::pair<double,double> GetMCLimitingEfficiency(double eff) const;
		std::pair<double,double> GetMCLimitingZenith(double zen) const;
		std::pair<double,double> GetMCLimitingOffset(double off) const;
		std::pair<double,double> GetMCLimitingEnergy(double en) const;

		bool IsItMCZenith(double zen) const;
		bool IsItMCEfficiency(double eff) const;
		bool IsItMCOffset(double off) const;
		bool IsItMCTelCode(unsigned telcode) const;
		bool IsItMCEnergy(double energy) const;

		void PrintMCParameters() const;

	private :

		void InitMCParameters();

		std::vector<double> fZenithMC; ///< mc zenith
		std::vector<double> fEfficiencyMC; ///< mc efficiency
		std::vector<double> fOffsetMC; ///< mc offset
		std::vector<TString> fTelTypeMC; ///< mc teltype
		std::vector<double> fEnergyMC; ///< mc energy
		std::vector<unsigned int> fTelCodeMC; ///< telcode

		START::STARTUtils::IrfOpt fIrfOpt;
		
#ifdef CONSTRUCT_LIBRARY
	public:
		ClassDef(START::MonteCarlo,1);
#endif
	};
}
#endif
