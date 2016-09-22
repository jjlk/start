#ifndef _RESIDUALS_
#define _RESIDUALS_

// STL
#include <vector>

// ROOT
#include <TObject.h>
#include <Rtypes.h>
// START

namespace START {
	/**
	 * \brief Contains residuals, flux, 1sigma and 3sigmas errors
	 *
	 * \author HAP-Fr team
	 */

	class Residuals : public TObject
	{

	public:

		Residuals(unsigned int size=0);
		~Residuals();

		void Print(Option_t *option="") const;
		void PrintResiduals(std::ostream &os=std::cout) const;

		std::vector<double> GetEnergy() const {return fEnergy;};
		void SetEnergy(std::vector<double> Energy) {fEnergy=Energy;};

		std::vector<double> GetEnergyMinus() const {return fEnergyMinus;};
		void SetEnergyMinus(std::vector<double> EnergyMinus) {fEnergyMinus=EnergyMinus;};

		std::vector<double> GetEnergyPlus() const {return fEnergyPlus;};
		void SetEnergyPlus(std::vector<double> EnergyPlus) {fEnergyPlus=EnergyPlus;};

		std::vector<double> GetResiduals() const {return fResiduals;};
		void SetResiduals(std::vector<double> Res) {fResiduals=Res;};

		std::vector<double> GetResidualsOn() const {return fResidualsOn;};
		void SetResidualsOn(std::vector<double> Res) {fResidualsOn=Res;};

		std::vector<double> GetResidualsOff() const {return fResidualsOff;};
		void SetResidualsOff(std::vector<double> Res) {fResidualsOff=Res;};

		
		std::vector<double> GetResidualsSigmaPlus() const {return fResidualsSigmaPlus;};
		void SetResidualsSigmaPlus(std::vector<double> ResidualsSigmaPlus) {fResidualsSigmaPlus=ResidualsSigmaPlus;};

		std::vector<double> GetResiduals3SigmaPlus() const {return fResiduals3SigmaPlus;};
		void SetResiduals3SigmaPlus(std::vector<double> Residuals3SigmaPlus) {fResiduals3SigmaPlus=Residuals3SigmaPlus;};

		std::vector<double> GetResidualsSigmaMinus() const {return fResidualsSigmaMinus;};
		void SetResidualsSigmaMinus(std::vector<double> ResidualsSigmaMinus) {fResidualsSigmaMinus=ResidualsSigmaMinus;};

		std::vector<double> GetResiduals3SigmaMinus() const {return fResiduals3SigmaMinus;};
		void SetResiduals3SigmaMinus(std::vector<double> Residuals3SigmaMinus) {fResiduals3SigmaMinus=Residuals3SigmaMinus;};

		std::vector<double> GetFlux() const {return fFlux;};
		void SetFlux(std::vector<double> Flux) {fFlux=Flux;};

		std::vector<double> GetFluxSigmaPlus() const {return fFluxSigmaPlus;};
		void SetFluxSigmaPlus(std::vector<double> FluxSigmaPlus) {fFluxSigmaPlus=FluxSigmaPlus;};

		std::vector<double> GetFlux3SigmaPlus() const {return fFlux3SigmaPlus;};
		void SetFlux3SigmaPlus(std::vector<double> Flux3SigmaPlus) {fFlux3SigmaPlus=Flux3SigmaPlus;};

		std::vector<double> GetFluxSigmaMinus() const {return fFluxSigmaMinus;};
		void SetFluxSigmaMinus(std::vector<double> FluxSigmaMinus) {fFluxSigmaMinus=FluxSigmaMinus;};

		std::vector<double> GetFlux3SigmaMinus() const {return fFlux3SigmaMinus;};
		void SetFlux3SigmaMinus(std::vector<double> Flux3SigmaMinus) {fFlux3SigmaMinus=Flux3SigmaMinus;};

		std::vector<bool> GetIsUpperLimits() const {return fIsUpperLimits;};
		void SetIsUpperLimits(std::vector<bool> IsUpperLimits) {fIsUpperLimits=IsUpperLimits;};

		Color_t GetMarkerColor() const {return fMarkerColor;};
		void SetMarkerColor(Color_t MarkerColor) {fMarkerColor=MarkerColor;};

		Size_t GetMarkerSize() const {return fMarkerSize;};
		void SetMarkerSize(Size_t MarkerSize) {fMarkerSize=MarkerSize;};

		Style_t GetMarkerStyle() const {return fMarkerStyle;};
		void SetMarkerStyle(Style_t MarkerStyle) {fMarkerStyle=MarkerStyle;};

		Color_t GetLineColor() const {return fLineColor;};
		void SetLineColor(Color_t LineColor) {fLineColor=LineColor;};

		Width_t GetLineWidth() const {return fLineWidth;};
		void SetLineWidth(Width_t LineWidth) {fLineWidth=LineWidth;};

	private:
		std::vector<double> fEnergy; ///< mean energy bin
		std::vector<double> fEnergyMinus; ///< mean energy bin
		std::vector<double> fEnergyPlus; ///< mean energy bin
		std::vector<double> fResiduals; ///< residuals
		std::vector<double> fResidualsSigmaPlus; ///< residuals error + at 1 sigma
		std::vector<double> fResidualsSigmaMinus; ///< residuals error - at 1 sigma
		std::vector<double> fResiduals3SigmaPlus; ///< residuals error + at 3 sigma
		std::vector<double> fResiduals3SigmaMinus; ///< residuals error - at 3 sigma
		std::vector<double> fFlux; ///< experimental flux
		std::vector<double> fFluxSigmaPlus; ///< error + on experimental flux at 1 sigma
		std::vector<double> fFluxSigmaMinus; ///< error - on experimental flux at 1 sigma
		std::vector<double> fFlux3SigmaPlus; ///< error + on experimental flux at 3 sigma
		std::vector<double> fFlux3SigmaMinus; ///< error - on experimental flux at 3 sigma
		std::vector<bool> fIsUpperLimits; // true if upper limit
		// JLK ADD FOR ADA
		std::vector<double> fResidualsOn; ///< residuals
		std::vector<double> fResidualsOff; ///< residuals
		
		Color_t fMarkerColor; ///< Marker color index
		Size_t fMarkerSize; ///< Marker size
		Style_t fMarkerStyle; //< Marker style
		Color_t fLineColor; ///< Line color index
		Width_t fLineWidth; ///< Line width

#ifdef CONSTRUCT_LIBRARY
	public:
		ClassDef(START::Residuals,1);
#endif
	};
}
#endif
