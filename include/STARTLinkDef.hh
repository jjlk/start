// STL
#include <utility>
#include <map>
#include <vector>
#include <string>

// ROOT
#include <TString.h>
#include <TH1F.h>
//#include <TPolyLine.h>
//#include <Rtypes.h>
// START
//#include "include/Band.hh"
//#include "include/Hypothesis.hh"

#ifdef __CINT__

//#pragma link off all globals;
//#pragma link off all classes;
//#pragma link off all functions;

#pragma link C++ class START::MonteCarlo+;
#pragma link C++ class START::ComputeResults+;
#pragma link C++ class START::Config+;
#pragma link C++ class START::DataSummary+;
#pragma link C++ class START::FCNLikelihood+;
#pragma link C++ class START::HandleResolArea+;
#pragma link C++ class START::Hypothesis+;
#pragma link C++ class START::PowerLaw+;
#pragma link C++ class START::ExpoCutOffPowerLaw+;
#pragma link C++ class START::LogParabolic+;
#pragma link C++ class START::BrokenPowerLaw+;
#pragma link C++ class START::SmoothBrokenPowerLaw+;
#pragma link C++ class START::SuperExpoCutOffPowerLaw+;
#pragma link C++ class START::EnergyBin+;
#pragma link C++ class START::Band+;
#pragma link C++ class START::BandsFactory+;
#pragma link C++ class START::MinimizeFactory+;
#pragma link C++ class START::PlotFactory+;
#pragma link C++ class START::STARTUtils+;
#pragma link C++ class START::ResidualsFactory+;
#pragma link C++ class START::Residuals+;
#pragma link C++ class START::SumHypothesis+;
#pragma link C++ class START::LightCurveFactory+;
#pragma link C++ class START::TimeBin+;
#pragma link C++ class START::TimeBinVector+;
#pragma link C++ class START::MultiWaveLengthFactory+;
#pragma link C++ class START::STARTFITSUtils+;
#pragma link C++ class START::Event+;


// iterator
#pragma link C++ class vector<START::Band>;
#pragma link C++ class vector<START::Band>::iterator;
#pragma link C++ class vector<START::Band>::const_iterator;

#pragma link C++ class vector<START::EnergyBin>;
#pragma link C++ class vector<START::EnergyBin>::iterator;
#pragma link C++ class vector<START::EnergyBin>::const_iterator;

#pragma link C++ class vector<START::TimeBin>;
#pragma link C++ class vector<START::TimeBin>::iterator;
#pragma link C++ class vector<START::TimeBin>::const_iterator;

#pragma link C++ class vector<START::Hypothesis*>;
#pragma link C++ class vector<START::Hypothesis*>::iterator;
#pragma link C++ class vector<START::Hypothesis*>::const_iterator;

#pragma link C++ class vector<START::Event>;
#pragma link C++ class vector<START::Event>::iterator;
#pragma link C++ class vector<START::Event>::const_iterator;

// Hypothesis transforme des warnings en errors. Pas sur de comprendre pourquoi j'ai besoin de ceux la...

//#pragma link C++ class pair<string,TimeBinVector>;
//#pragma link C++ class pair<string,Residuals>;
//#pragma link C++ class pair<string,TPolyLine*>;

//#pragma link C++ class map<string,TimeBinVector>;
//#pragma link C++ class map<string,TimeBinVector>::iterator;
//#pragma link C++ class map<string,TimeBinVector>::const_iterator;

//#pragma link C++ class map<string,Residuals>;
//#pragma link C++ class map<string,Residuals>::iterator;
//#pragma link C++ class map<string,Residuals>::const_iterator;

//#pragma link C++ class map<string,TPolyLine*>;
//#pragma link C++ class map<string,TPolyLine*>::iterator;
//#pragma link C++ class map<string,TPolyLine*>::const_iterator;

// Hypothesis
//#pragma link C++ class pair<double,bool>;
//#pragma link C++ class pair<unsigned int,double>;
//#pragma link C++ class vector<pair<unsigned int,double> >;
//#pragma link C++ class pair<unsigned int,pair<double,bool> >;
//#pragma link C++ class pair<unsigned int,vector<pair<double,double> > >;
//#pragma link C++ class vector<pair<vector<pair<double,double> >,pair<int,int> > >;
//#pragma link C++ class vector<pair<unsigned int,vector<pair<double,double> > > >;
//#pragma link C++ class pair<vector<double>,vector<double> >;
//#pragma link C++ class pair<vector<pair<double,double> >,pair<int,int> >;

// Hypothesis
#pragma link C++ class pair<unsigned int,double>;
#pragma link C++ class pair<double,bool>;
#pragma link C++ class vector<pair<double,double> >;
#pragma link C++ class pair<unsigned int,pair<double,bool> >;
#pragma link C++ class pair<unsigned int,vector<pair<double,double> > >;
#pragma link C++ class pair<vector<double>,vector<double> >+;
#pragma link C++ class pair<vector<pair<double,double> >,pair<int,int> >;
#pragma link C++ class vector<pair<vector<pair<double,double> >,pair<int,int> > >;
#pragma link C++ class vector<pair<unsigned int,vector<pair<double,double> > > >;
#pragma link C++ class vector<pair<unsigned int,double> >;

// HandleResolArea
#pragma link C++ class map<TString, TH1F *>;
#pragma link C++ class map<TString, TH1F *>;
#pragma link C++ class map<TString, TH1F *>;
#pragma link C++ class map<TString, pair<vector<double >, vector<double > > >;
#pragma link C++ class map<TString, pair<vector<double >, vector<double > > >;
#pragma link C++ class map<TString, pair<vector<double >, vector<double > > >;



#endif
