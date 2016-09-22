#ifndef _HYPOTHESIS_
#define _HYPOTHESIS_

// STL
#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <map>

// ROOT
#include <TNamed.h>
#include <TF1.h>
#include <TString.h>
#include <TPolyLine.h>
#include <Math/GSLIntegrator.h>
#include <Math/Functor.h>
#include <Math/AllIntegrationTypes.h>
#include <Math/RootFinder.h>

// START
#include "Band.hh"
#include "TimeBinVector.hh"
#include "Residuals.hh"

namespace START {
  
	/**
	 * \brief Abstract class used to facilitate implementation and the use of hypothesis.
	 * Can not be instantiate
	 *
	 * \author HAP-Fr team
	 */

	class Hypothesis : public TNamed
	{

	public:

		/**
		 * \brief Options for the spetral type of the hypothesis.
		 * Main canvas with spectrum and residuals is only done for type Differential.
		 *
		 * \param Differential \f$ cm^{-2}.s^{-1}.TeV^{-1} \f$
		 * \param Integrated \f$ cm^{-2}.s^{-1} \f$
		 * \param EnergyFlux \f$ erg.cm^{-2}.s^{-1} \f$
		 */
		typedef enum {Differential, Integrated, EnergyFlux} SpectralType;

		/**
		 * \brief Options for the butterfly's drawing on the spectrum canvas.
		 * \param LinearButterfly is drawn with the covaricance of the flux at 1st order
		 * \param LogarithmButterfly butterfly is drawn with the log10 covaricance of the flux at 1st order
		 * \param ContoursButterfly butterfly is drawn with the contours which have to be computed before
		 */
		typedef enum {LinearButterfly, LogarithmButterfly, ContoursButterfly} ButterflyType;

		Hypothesis(); // default constructor for ROOT
		Hypothesis(TString hypothesisname); // constructor used by Hypothesis
		virtual ~Hypothesis(); // destructor

		virtual Hypothesis* clone() const = 0; // VIM : Tricks to create a copy constructor
  
		/* Functions to be implemented by the user */

		virtual double GetFlux(double x) const = 0; // return the value of the flux

		virtual double GetFluxFitParams(double x) const; // return the value of the flux with fitted parameters

		virtual double GetMeanBinEnergy(double x1,double x2) const; // return the mean energy bin for energy x1 and x2.
		// Implement this function if you know the analytic form of the flux
  
		virtual void DerivativesFormulae(double x) = 0; // used for the butterfly, compute first derivatives of flux for param 1,2.. n at energy x
  
		virtual void InitParametersCaracteristics() = 0; // Initialization of parameters names and units

		/*Return the flux at fixed parameters times the energy */

		virtual double GetFluxFitParamsTimesEnergy(double x) const; // return the value of the flux times the energy x
		virtual double GetFluxTimesEnergy(double x) const; // return the value of the flux times the energy x

		/* Function allowing to put an hypothesis in a TF1 */

		virtual double operator() (double *x, double *p=0);

		/* Return the variance of the flux */

		virtual double GetSigmaFlux(double x); // Return the variance of the flux at first order

		/* Return the width of the butterfly at a fixed energy */

		virtual double GetLinearWidthSigma(double x); // return the value of log10(flux+sigma) - log10(flux-sigma) at a fixed energy
		virtual double GetLogarithmWidthSigma(double x); // return the value of Phi*exp(+sigma/Phi(E)) - Phi*exp(-sigma/Phi(E)) at a fixed energy
		virtual double GetContoursWidthSigma(double x); // return the value of width of the butterfly defined by the contours at a fixed energy

		/* Return the derivative variance of the flux */

		virtual double GetFirstDerivativeLinearWidthSigma(double x); // return the value of the derivative of Hypothesis::GetLinearWidthSigma
		virtual double GetFirstDerivativeLogarithmWidthSigma(double x); // return the value of the derivative of Hypothesis::GetLogarithmWidthSigma
		virtual double GetFirstDerivativeContoursWidthSigma(double x); // return the value of the derivative of Hypothesis::GetContoursWidthSigma

		/* Find decorrelation energy */

		virtual double FindLinearDecorrelationEnergy(); // find decorrelation energy with linear covariance
		virtual double FindLogarithmDecorrelationEnergy(); // find decorrelation energy with log flux covariance
		virtual double FindContoursDecorrelationEnergy(); // find decorrelation energy with contours

		/* Return decorrelation energy */

		/*
		  virtual double GetDecorrelationEnergyFromLinearCovariance(); // return the decorrelation energy from flux covariance
		  virtual double GetDecorrelationEnergyFromLogarithmCovariance(); // return the decorrelation energy from log flux covariance
		  virtual double GetDecorrelationEnergyFromContours(); // return the decorrelation energy from contours
		*/

		virtual double GetLinearDecorrelationEnergy(); // return the decorrelation energy from flux covariance
		virtual double GetLogarithmDecorrelationEnergy(); // return the decorrelation energy from log flux covariance
		virtual double GetContoursDecorrelationEnergy(); // return the decorrelation energy from contours

		/* Normalization has to be used in constructors to define the normalization parameters*/

		virtual void SetLowerLimitOnNormalizationParameter(unsigned int param, double value=1.e-18); // set a lowerlimit on the normalisation parameter "param" with value "value"

		/* Fix a parameter for the minimization and can redo the minimization with the parameter released */

		virtual void FixParameter(unsigned int param, double value, bool redominimization=false);
		virtual void FixParameter(std::string param, double value, bool redominimization=false);
		virtual void LimitParameter(unsigned int param, double valuemin , double valuemax);

		virtual std::vector<std::pair<unsigned int, std::pair<double, bool> > > GetFixedParameters() const {return ffixedparameter;};
		//virtual void SetFixedParameters(std::vector<std::pair<unsigned int, std::pair<double, bool> > > fixedparameter) {ffixedparameter=fixedparameter}
		virtual std::vector<std::pair<unsigned int, std::pair<double, double> > > GetLimitedParameters() const {return flimitedparameter;};

		/* Set and get energy range for of integrated fit */

		virtual void SetIntegratedFitEnergyRange(double emin,double emax) {fintegratedfitenergyrange=std::make_pair(emin,emax);}; 
		///< Set range for integration for case Integrated and EnergyFlux

		virtual std::pair<double,double> GetIntegratedFitEnergyRange() const {return fintegratedfitenergyrange;};
		///< Get range for integration for case Integrated and EnergyFlux

		/* Integrated flux */
		virtual double GetFluxIntegral(double x1, double x2) const;
		virtual double GetFluxTimesEnergyIntegral(double x1, double x2) const;
		virtual double GetFluxIntegralFitParams(double x1, double x2) const;
		virtual double GetFluxTimesEnergyIntegralFitParams(double x1, double x2) const;

		/* Set functions */

		virtual void SetEref(double const eref) {feref=eref;}; ///< set reference energy
		virtual void SetParameters(std::vector<double> const param) {fparam=param;}; ///< set the parameters
		virtual void SetParametersErrors(std::vector<double> const paramerr) {fparamerr=paramerr;}; ///< set the parameters' errors
		virtual void SetFittedParameters(std::vector<double> const paramfit) {fparamfit=paramfit;}; ///< set the fitted parameters
		virtual void SetFittedParametersErrors(std::vector<double> const paramerrfit) {fparamerrfit = paramerrfit;}; ///< set the parameters' errors
		virtual void SetFittedCovarianceMatrix(std::vector<std::vector<double> > covariancefit) {fcovariancefit=covariancefit;}; ///< set the fitted covariance matrix
		virtual void SetParametersNames(std::vector<std::string> const paramname) {fparamname=paramname;}; ///< set the parameters' name
		virtual void SetLatexParametersNames(std::vector<std::string> const latexparamname) {flatexparamname=latexparamname;}; ///< set the parameters' name
		virtual void SetParametersUnits(std::vector<std::string> const paramunits) {fparamunits=paramunits;}; ///< set the parameters'units
		virtual void SetFittedIntegratedFlux(double flux, double errorflux) {ffittedintegratedflux=std::make_pair(flux,errorflux);}; ///< set fitted integrated flux and error
		virtual void SetFittedEnergyFlux(double energyflux, double errorenergyflux) {ffittedenergyflux=std::make_pair(energyflux,errorenergyflux);}; ///< set fitted energy flux and error

		virtual void SetMinimizationEnergyRange(double emin, double emax) {fminimizationenergyrange=std::make_pair(emin,emax);}; ///< set minimization energy range

		virtual void SetModelName(std::string modelname) {fmodelname=modelname;}; ///< set model name

		virtual void SetMinosErrors(std::vector<std::pair<double,double> > const minoserrors) {fminoserrors=minoserrors;}; ///< set minos erros

		virtual void SetSpectralType(SpectralType type) {fSpectralType=type;};

		virtual void SetBandArray(std::vector<START::Band> const BandArray) {fBandArray=BandArray;}; ///< set fBandArray

		virtual void AddTimeBinVector(std::string name, const TimeBinVector TimeBins) {fMapTimeBinVector[name]=TimeBins;}; ///< set fTimeBinVector

		virtual void SetLightCurveIntegratedFluxEnergyRange(double emin, double emax) 
		{fLightCurveIntegratedFluxEnergyRange=std::make_pair(emin,emax);}; ///< set lc integrated flux energy range
		virtual std::pair<double,double> GetLightCurveIntegratedFluxEnergyRange() 
		{return fLightCurveIntegratedFluxEnergyRange;}; ///< get lc integrated flux energy range

		virtual void AddResiduals(const Residuals Res, double sigma=0.) {fMapResiduals[MakeKeyResiduals(sigma)]=Res;}; ///< set fTimeBinVector

		virtual void AddButterfly(ButterflyType TypeOfButterfly, const TPolyLine Butt) 
		{fMapButterfly[MakeKeyButterfly(TypeOfButterfly)]=new TPolyLine(Butt);};

		/* Get functions */

		virtual double GetEref() const {return feref;}; ///< return reference energy
		virtual unsigned int GetParametersNb() const {return fnumberofparameters;}; ///< return the number of parameters for the hypothesis
		virtual std::vector<double> GetParameters() const {return fparam;}; ///< get the parameters
		virtual std::vector<double> GetParametersErrors() const {return fparamerr;}; ///< get the parameters' errors
		virtual std::vector<double> GetFittedParameters() const {return fparamfit;}; ///< get the parameters
		virtual std::vector<double> GetFittedParametersErrors() const {return fparamerrfit;}; ///< get the parameters' errors
		virtual std::vector<std::vector<double> > GetFittedCovarianceMatrix() const {return fcovariancefit;}; ///< get the fitted covariance matrix
		virtual std::vector<std::string> GetParametersNames() const {return fparamname;}; ///< get the parameters' name
		virtual std::vector<std::string> GetLatexParametersNames() const {return flatexparamname;}; ///< get the parameters' name
		virtual std::vector<std::string> GetParametersUnits() const {return fparamunits;}; ///< get the parameters' units
		virtual std::pair<double,double> GetFittedIntegratedFlux() const {return ffittedintegratedflux;}; ///< get fitted integrated flux and error
		virtual std::pair<double,double> GetFittedEnergyFlux() const {return ffittedenergyflux;}; ///< get fitted energy flux and error
		virtual std::vector<double> GetFirstDerivatives() const {return ffirstderivatives;}; ///< get first derivatives at current energy
		virtual std::string GetModelName() const {return fmodelname;}; ///< get model name

		virtual std::pair<double,double> GetMinimizationEnergyRange() const {return fminimizationenergyrange;}; ///< get minimization energy range

		virtual std::string GetIntegratedFluxUnits() const; ///< get energy flux units
		virtual std::string GetEnergyFluxUnits() const; ///< get integrated flux units
		virtual std::string GetFluxUnits() const; ///< get flux units

		virtual std::string GetLatexFormula() const {return flatexformula;}; ///< get latex formula

		virtual std::vector<std::pair<double,double> > GetMinosErrors() {return fminoserrors;}; ///< get minos erros

		virtual SpectralType GetSpectralType() const {return fSpectralType;};

		virtual std::vector<std::pair<unsigned int,double> > GetVectorNormalizedParameters() const {return fnormalizedparam;}; ///< get the list of fixed parameters

		virtual std::vector<START::Band> GetBandArray() const {return fBandArray;}; ///< get fBandArray

		virtual std::map<std::string,START::TimeBinVector> GetMapTimeBinVector() const {return fMapTimeBinVector;}; ///< get fTimeBinVector

		virtual std::map<std::string,START::Residuals> GetMapResiduals() const {return fMapResiduals;}; ///< get fMapResiduals

		virtual std::map<std::string,TPolyLine*> GetMapButterfly() const {return fMapButterfly;}; ///< get fMapResiduals

		/* Print functions */
		virtual void Print(Option_t *option="") const;
		virtual void PrintHypothesis(std::ostream &os=std::cout) const;
		virtual void PrintParameters(std::ostream &os=std::cout) const; // print parameters
		virtual void PrintFittedParameters(std::ostream &os=std::cout) const; // print fitted parameters
		virtual void PrintFittedCovarianceMatrix(std::ostream &os=std::cout) const; // print the covariance matrix
		virtual void PrintMinimizationCarateristics(std::ostream &os=std::cout) const; // print minimization parameters
		virtual void PrintResiduals(std::ostream &os=std::cout) const; // print residuals
		virtual void PrintButterflies(std::ostream &os=std::cout) const; // print residuals
		virtual void PrintSummaryBand(std::ostream &os=std::cout,bool UseEThreshold=true) const; // print bands summary
		virtual void PrintContours(std::ostream &os=std::cout) const; // print contours
		/* Residuals */

		virtual std::vector<double> GetResiduals() const {return fresiduals;}; ///< get residuals
		// JLK ADD FOR ADA
		virtual std::vector<double> GetResidualsOn() const {return fresiduals_on;}; ///< get residuals
		virtual std::vector<double> GetResidualsOff() const {return fresiduals_off;}; ///< get residuals
		virtual std::vector<double> GetResidualsSigmaPlus() const {return fresidualssigmaplus;}; ///< get residuals error + at 1 sigma
		virtual std::vector<double> GetResidualsSigmaMinus() const {return fresidualssigmaminus;}; ///< get residuals error - at 1 sigma
		virtual std::vector<double> GetResiduals3SigmaPlus() const {return fresiduals3sigmaplus;}; ///< get residuals error + at 3 sigma
		virtual std::vector<double> GetResiduals3SigmaMinus() const {return fresiduals3sigmaminus;}; ///< get residuals error - at 3 sigma
		virtual std::vector<double> GetResidualsFlux() const {return fflux;}; ///< get experimental flux
		virtual std::vector<double> GetResidualsFluxSigmaPlus() const {return ffluxsigmaplus;}; ///< get error + on experimental flux at 1 sigma
		virtual std::vector<double> GetResidualsFluxSigmaMinus() const {return ffluxsigmaminus;}; ///< get error - on experimental flux at 1 sigma
		virtual std::vector<double> GetResidualsFlux3SigmaPlus() const {return fflux3sigmaplus;}; ///< get error + on experimental flux at 3 sigma
		virtual std::vector<double> GetResidualsFlux3SigmaMinus() const {return fflux3sigmaminus;}; ///< get error - on experimental flux at 3 sigma
		virtual std::vector<double> GetMeanEnergy() const {return femean;}; ///< get experimental flux

		virtual void SetResiduals(std::vector<double> residuals) {fresiduals = residuals;}; ///< set residuals
		// JLK ADD FOR ADA
		virtual void SetResidualsOn(std::vector<double> residuals) {fresiduals_on = residuals;}; ///< set residuals
		virtual void SetResidualsOff(std::vector<double> residuals) {fresiduals_off = residuals;}; ///< set residuals
		virtual void SetResidualsSigmaPlus(std::vector<double> residualssigmaplus) {fresidualssigmaplus = residualssigmaplus;}; ///< set residuals error + at 1 sigma
		virtual void SetResidualsSigmaMinus(std::vector<double> residualssigmaminus) {fresidualssigmaminus = residualssigmaminus;}; ///< set residuals error - at 1 sigma
		virtual void SetResiduals3SigmaPlus(std::vector<double> residuals3sigmaplus) {fresiduals3sigmaplus = residuals3sigmaplus;}; ///< set residuals error + at 3 sigma
		virtual void SetResiduals3SigmaMinus(std::vector<double> residuals3sigmaminus) {fresiduals3sigmaminus = residuals3sigmaminus;}; ///< set residuals error - at 3 sigma
		virtual void SetResidualsFlux(std::vector<double> flux) {fflux = flux;}; ///< set experimental flux
		virtual void SetResidualsFluxSigmaPlus(std::vector<double> fluxsigmaplus) {ffluxsigmaplus = fluxsigmaplus;}; ///< set error + on experimental flux at 1 sigma
		virtual void SetResidualsFluxSigmaMinus(std::vector<double> fluxsigmaminus) {ffluxsigmaminus = fluxsigmaminus;}; ///< set error - on experimental flux at 1 sigma
		virtual void SetResidualsFlux3SigmaPlus(std::vector<double> flux3sigmaplus) {fflux3sigmaplus = flux3sigmaplus;}; ///< set error + on experimental flux at 3 sigma
		virtual void SetResidualsFlux3SigmaMinus(std::vector<double> flux3sigmaminus) {fflux3sigmaminus = flux3sigmaminus;}; ///< set error - on experimental flux at 3 sigma
		virtual void SetMeanEnergy(std::vector<double> emean) {femean = emean;}; ///< set mean energy bin

		/* MINUIT2's parameters */
		virtual double GetMaximumLikelihood() const {return fmaximumlikelihood;}; ///< get the maximum of the likelihood after minimization
		virtual double GetEDM() const {return fedm;}; ///< Get associated variable fEDMValue
		virtual int GetIteration() const {return fiteration;}; ///< Get associated variable fIteration
		virtual bool GetIsMinimumValid() const {return fisminimumvalid;}; ///< Get associated variable fIsMinimumValid
		virtual bool GetIsCovarianceValid() const {return fiscovariancevalid;}; ///< Get associated variable fIsCovarianceValid
		virtual bool GetAreFittedParametersValid() const {return farefittedparametersvalid;}; ///< Get associated variable fAreFittedParametersValid
		virtual bool GetConvergence() const {return fconvergence;}; ///< Get associated variable fconvergence

		virtual std::vector<std::pair<std::vector<std::pair<double,double> >, std::pair<int,int> > > GetContourSigma1() const {return fcontoursigma1;}; ///< get contours at sigma1 
		virtual std::vector<std::pair<std::vector<std::pair<double,double> >, std::pair<int,int> > > GetContourSigma2() const {return fcontoursigma2;}; ///< get contours at sigma2 
		virtual std::vector<std::pair<std::vector<std::pair<double,double> >, std::pair<int,int> > > GetContourSigma3() const {return fcontoursigma3;}; ///< get contours at sigma3

		virtual std::vector<std::pair<unsigned int,std::vector<std::pair<double,double> > > > GetScansLikelihood() const {return fscanslikelihood;};

		virtual void SetMaximumLikelihood(double const maximumlikelihood) {fmaximumlikelihood=maximumlikelihood;}; ///< set the maximum of the likelihood after minimization
		virtual void SetEDM(double const edm) {fedm=edm;}; ///< set expected edm after minimization
		virtual void SetIteration(int const iteration) {fiteration=iteration;}; ///< set expected edm after minimization
		virtual void SetIsMinimumValid(bool const isminimumvalid) {fisminimumvalid=isminimumvalid;}; ///< set associated variable fisminimumcalid
		virtual void SetIsCovarianceValid(bool const iscovariancevalid) {fiscovariancevalid=iscovariancevalid;}; ///< set associated variable fiscovariancevalid
		virtual void SetAreFittedParametersValid(bool const arefittedparametersvalid) 
		{farefittedparametersvalid=arefittedparametersvalid;}; ///< set associated variable arefittedparametersvalid
		virtual void SetConvergence(bool const convergence) {fconvergence = convergence;}; ///< set associated variable fconvergence

		virtual void SetContourSigma1(std::vector<std::pair<std::vector<std::pair<double,double> >, std::pair<int,int> > > contoursigma1) 
		{fcontoursigma1=contoursigma1;}; ///< set contours at sigma1
		virtual void SetContourSigma2(std::vector<std::pair<std::vector<std::pair<double,double> >, std::pair<int,int> > > contoursigma2)
		{fcontoursigma2=contoursigma2;}; ///< set contours at sigma2
		virtual void SetContourSigma3(std::vector<std::pair<std::vector<std::pair<double,double> >, std::pair<int,int> > > contoursigma3)
		{fcontoursigma3=contoursigma3;}; ///< set contours at sigma3

		virtual void SetScansLikelihood(std::vector<std::pair<unsigned int,std::vector<std::pair<double,double> > > > scanslikelihood)
		{fscanslikelihood=scanslikelihood;}; ///< set region scans

		static std::string MakeKeyResiduals(double sigma=0.);
		static std::string MakeKeyButterfly(ButterflyType TypeOfButterfly);

	protected:

		/* Parameters */
		unsigned int fnumberofparameters; ///< number of parameters of the hypothesis
		double feref; ///< reference energy
		std::vector<std::string> fparamunits; ///< strings with parameters' units
		mutable std::vector<double> fparam; ///< parameters of the hypothesis 
		std::vector<double> fparamfit; ///< fitted parameters of the hypothesis 
		std::vector<double> fparamerr; ///< parameters' errors 
		std::vector<double> fparamerrfit; ///< fitted parameters' errors 
		std::vector<std::vector<double> > fcovariancefit; ///< covariance matrix
		std::vector<std::string> fparamname; ///< parameters' names
		std::vector<std::string> flatexparamname; ///< parameters' names written in latex
		double fmaximumlikelihood; ///< maximum of the likelihood sor the fitted parameters
		double fedm; ///< normalized distance to the minimum
		int fiteration; ///< number of iteration fro the minimization
		bool fisminimumvalid; ///< if true the minimum is valid
		bool fiscovariancevalid; ///< if true the covariance is valid
		bool farefittedparametersvalid; ///< if true fitted parameters are valid
		bool fconvergence; ///< if true , fisminimumvalid && fiscovariancevalid && farefittedparametersvalid are true
		double flineardecorrelationenergy; // linear decorrelation energy
		double flogarithmdecorrelationenergy; // logarithm decorrelation energy
		double fcontoursdecorrelationenergy; // contours decorrelation energy

		std::string flatexformula; ///< latex formula of the hypothesis (root's way)
		std::string fmodelname; ///< name of the model

		std::pair<double,double> ffittedintegratedflux; ///< fitted integrated flux + errors
		std::pair<double,double> ffittedenergyflux; ///< fitted energy flux + errors

		std::pair<double,double> fintegratedfitenergyrange; ///< energy range for the integrated flux and the energy flux

		std::pair<double,double> fminimizationenergyrange; ///< minimization range

		/* Internal attributs*/ 
		std::vector<std::pair<unsigned int,double> > fnormalizedparam; ///< number of the normalized parameter and the lower limit
		std::vector<std::pair<unsigned int, std::pair<double, bool> > > ffixedparameter; ///< parameter index, value and redominimization
		//ADA
		std::vector<std::pair<unsigned int, std::pair<double, double> > > flimitedparameter; ///< number of the normalized parameter and the lower limit

		std::vector<double> ffirstderivatives; ///< first derivatives of the flux for parameters 0, 1, 2 .. N

		/* Contours */
		std::vector<std::pair<std::vector<std::pair<double,double> >, std::pair<int,int> > > fcontoursigma1; ///< contours sigma1
		std::vector<std::pair<std::vector<std::pair<double,double> >, std::pair<int,int> > > fcontoursigma2; ///< contours sigma2
		std::vector<std::pair<std::vector<std::pair<double,double> >, std::pair<int,int> > > fcontoursigma3; ///< contours sigma3

		/* Scans */
		std::vector<std::pair<unsigned int,std::vector<std::pair<double,double> > > > fscanslikelihood; ///< scans (param,points)

		/* Minos errors */
		std::vector<std::pair<double,double> > fminoserrors; ///< minos errors (-,+)

		/* Data binned in energy */
		std::vector<START::Band> fBandArray; ///< Bands containing data, bin emean, sth for residuals and plotfactory

		/* Data binned in time */
		std::map<std::string,START::TimeBinVector> fMapTimeBinVector; ///< Map of vector of bins containing data binned in time
		std::pair<double,double> fLightCurveIntegratedFluxEnergyRange;

		/* Residuals */
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

		std::map<std::string,Residuals> fMapResiduals;

		/* Butterflies */
		std::map<std::string,TPolyLine*> fMapButterfly; 

		/* Integrator */
		ROOT::Math::GSLIntegrator *fIntegrantFlux; //!
		ROOT::Math::GSLIntegrator *fIntegrantFluxTimesE; //!

		ROOT::Math::Functor1D *fFunctorIntegrantFlux; //!
		ROOT::Math::Functor1D *fFunctorIntegrantFluxTimesE; //!

		ROOT::Math::GSLIntegrator *fIntegrantFitFlux; //!
		ROOT::Math::GSLIntegrator *fIntegrantFitFluxTimesE; //!

		ROOT::Math::Functor1D *fFunctorIntegrantFitFlux; //!
		ROOT::Math::Functor1D *fFunctorIntegrantFitFluxTimesE; //!

		double frelativeprecision;
		double fabsoluteprecision;
		ROOT::Math::IntegrationOneDim::Type fIntegrationType;

		/* Root finder */

		ROOT::Math::RootFinder::EType fRootFinderType;

		/* Spectral type */

		SpectralType fSpectralType;

#ifdef CONSTRUCT_LIBRARY
	public:
		ClassDef(START::Hypothesis,1);
#endif
	};
}
#endif
