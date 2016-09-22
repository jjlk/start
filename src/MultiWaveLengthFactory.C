// STL
#include <algorithm>
// ROOT

// START
#include "STARTUtils.hh"
#include "MultiWaveLengthFactory.hh"

#define DEBUG 0
#include "debugging.hh"

#ifdef CONSTRUCT_LIBRARY
ClassImp(START::MultiWaveLengthFactory)
#endif

/**
 * \brief Constructor
 */
START::MultiWaveLengthFactory::MultiWaveLengthFactory(std::string name, FluxUnits FUnits, EnergyUnits EUnits)
:TNamed(name.c_str(),name.c_str())
{

  fEnergyUnits = EUnits;
  fFluxUnits = FUnits;

}


/**
 * \brief Destructor
 */
START::MultiWaveLengthFactory::~MultiWaveLengthFactory()
{

  for(std::map<std::string,TPolyLine*>::iterator butt=fMapButterfly.begin(); butt!=fMapButterfly.end(); ++butt) {
    if(butt->second!=0) delete butt->second;
    butt->second=0;
  }
  fMapButterfly.clear();

  for(std::map<std::string,TF1*>::iterator tf1=fMapTF1.begin(); tf1!=fMapTF1.end(); ++tf1) {
    if(tf1->second!=0) delete tf1->second;
    tf1->second=0;
  }
  fMapTF1.clear();

  fMapResiduals.clear();

}


/**
 * \brief Copy constructor
 */
START::MultiWaveLengthFactory::MultiWaveLengthFactory(MultiWaveLengthFactory const &MwlCopy)
  :TNamed(MwlCopy)
{

  for(std::map<std::string,TPolyLine*>::iterator butt=fMapButterfly.begin(); butt!=fMapButterfly.end(); ++butt) {
    butt->second=0;
  }
  fMapButterfly.clear();
  for(std::map<std::string,TF1*>::iterator tf1=fMapTF1.begin(); tf1!=fMapTF1.end(); ++tf1) {
    if(tf1->second!=0) delete tf1->second;
    tf1->second=0;
  }
  fMapTF1.clear();
  fMapResiduals.clear();

  for(std::map<std::string,TPolyLine*>::const_iterator butt=MwlCopy.fMapButterfly.begin(); butt!=MwlCopy.fMapButterfly.end(); ++butt) {
    if(butt->second!=0) fMapButterfly[butt->first] = new TPolyLine(*(butt->second));
  }

  for(std::map<std::string,TF1*>::const_iterator tf1=MwlCopy.fMapTF1.begin(); tf1!=MwlCopy.fMapTF1.end(); ++tf1) {
    if(tf1->second!=0) fMapTF1[tf1->first] = new TF1(*(tf1->second));
  }

  fMapResiduals = MwlCopy.fMapResiduals;

}

/**
 * \brief Assignment operator
 */
START::MultiWaveLengthFactory &START::MultiWaveLengthFactory::operator=(MultiWaveLengthFactory const &MwlCopy)
{
  if(this != &MwlCopy) {

    for(std::map<std::string,TPolyLine*>::iterator butt=fMapButterfly.begin(); butt!=fMapButterfly.end(); ++butt) {
      butt->second=0;
    }
    fMapButterfly.clear();
    for(std::map<std::string,TF1*>::iterator tf1=fMapTF1.begin(); tf1!=fMapTF1.end(); ++tf1) {
      if(tf1->second!=0) delete tf1->second;
      tf1->second=0;
    }
    fMapTF1.clear();
    fMapResiduals.clear();
    
    for(std::map<std::string,TPolyLine*>::const_iterator butt=MwlCopy.fMapButterfly.begin(); butt!=MwlCopy.fMapButterfly.end(); ++butt) {
      if(butt->second!=0) fMapButterfly[butt->first] = new TPolyLine(*(butt->second));
    }
    
    for(std::map<std::string,TF1*>::const_iterator tf1=MwlCopy.fMapTF1.begin(); tf1!=MwlCopy.fMapTF1.end(); ++tf1) {
      if(tf1->second!=0) fMapTF1[tf1->first] = new TF1(*(tf1->second));
    }
    
    fMapResiduals = MwlCopy.fMapResiduals;
    
  }

  return (*this);
}


/**
 * \brief Add points
 * \warning Energy have to be in TeV and Flux in cm-2.s-1.TeV-1
 */
void START::MultiWaveLengthFactory::AddPoints(std::string objname, Residuals &Res,
				       Color_t MarkerColor, Style_t MarkerStyle, Size_t MarkerSize,
				       Color_t LineColor, Width_t LineWidth) {

  Res.SetMarkerColor(MarkerColor);
  Res.SetMarkerStyle(MarkerStyle);
  Res.SetMarkerSize(MarkerSize);
  Res.SetLineColor(LineColor);
  Res.SetLineWidth(LineWidth);

  if(DEBUG) {
    std::cout << "Before :" << std::endl;
    Res.Print();
  }
    

  std::vector<double> x = Res.GetEnergy();
  std::vector<double> xminus = Res.GetEnergyMinus();
  std::vector<double> xplus = Res.GetEnergyPlus();

  std::vector<double> y = Res.GetFlux();
  std::vector<double> y1minus = Res.GetFluxSigmaMinus();
  std::vector<double> y1plus = Res.GetFluxSigmaPlus();
  std::vector<double> y3minus = Res.GetFlux3SigmaMinus();
  std::vector<double> y3plus = Res.GetFlux3SigmaPlus();
  
  for(unsigned int i(0); i<y.size(); i++) {

    DEBUG_OUT_L(2) << "Before x=:" << x[i] << " y=" << y[i] << std::endl; 
    DEBUG_OUT_L(2) << "Before x=:" << x[i] << " y3+=" << y3plus[i] << std::endl; 

    switch(fFluxUnits) {
    case vFv:
      y[i] = y[i]*x[i]*x[i]*STARTUtils::GetTeVToErg();
      y1minus[i] = y1minus[i]*x[i]*x[i]*STARTUtils::GetTeVToErg();
      y1plus[i] = y1plus[i]*x[i]*x[i]*STARTUtils::GetTeVToErg();
      y3minus[i] = y3minus[i]*x[i]*x[i]*STARTUtils::GetTeVToErg();
      y3plus[i] = y3plus[i]*x[i]*x[i]*STARTUtils::GetTeVToErg();
      break;
    case ESquareFlux_TeV:
      y[i] = y[i]*x[i]*x[i];
      y1minus[i] = y1minus[i]*x[i]*x[i];
      y1plus[i] = y1plus[i]*x[i]*x[i];
      y3minus[i] = y3minus[i]*x[i]*x[i];
      y3plus[i] = y3plus[i]*x[i]*x[i];
      break;
    case Differential_TeV:
      y[i]/=1.;
      y1minus[i]/=1.;
      y1plus[i]/=1.;
      break;
    case Differential_GeV:
      y[i]/=STARTUtils::GetTeVToGeV();
      y1minus[i]/=STARTUtils::GetTeVToGeV();
      y1plus[i]/=STARTUtils::GetTeVToGeV();
      y3minus[i]/=STARTUtils::GetTeVToGeV();
      y3plus[i]/=STARTUtils::GetTeVToGeV();
      break;
    case Differential_MeV:
      y[i]/=STARTUtils::GetTeVToMeV();
      y1minus[i]/=STARTUtils::GetTeVToMeV();
      y1plus[i]/=STARTUtils::GetTeVToMeV();
      y3minus[i]/=STARTUtils::GetTeVToMeV();
      y3plus[i]/=STARTUtils::GetTeVToMeV();
      break;
    case Differential_keV:
      y[i]/=STARTUtils::GetTeVTokeV();
      y1minus[i]/=STARTUtils::GetTeVTokeV();
      y1plus[i]/=STARTUtils::GetTeVTokeV();
      y3minus[i]/=STARTUtils::GetTeVTokeV();
      y3plus[i]/=STARTUtils::GetTeVTokeV();
    }
  }

  for(unsigned int i(0); i<x.size(); i++) {
    switch(fEnergyUnits) {
    case TeV:
      break;
    case GeV:
      x[i]*=STARTUtils::GetTeVToGeV();
      xplus[i]*=STARTUtils::GetTeVToGeV();
      xminus[i]*=STARTUtils::GetTeVToGeV();
      break;
    case MeV:
      x[i]*=STARTUtils::GetTeVToMeV();
      xplus[i]*=STARTUtils::GetTeVToMeV();
      xminus[i]*=STARTUtils::GetTeVToMeV();
      break;
    case keV:
      x[i]*=STARTUtils::GetTeVTokeV();
      xminus[i]*=STARTUtils::GetTeVTokeV();
      xplus[i]*=STARTUtils::GetTeVTokeV();
      break;
    case Hz:
      x[i]*=STARTUtils::GetTeVToHz();
      xplus[i]*=STARTUtils::GetTeVToHz();
      xminus[i]*=STARTUtils::GetTeVToHz();
    }

    DEBUG_OUT_L(2) << "After x=:" << x[i] << " y=" << y[i] << std::endl; 

  }

  Res.SetEnergy(x);
  Res.SetEnergyMinus(xminus);
  Res.SetEnergyPlus(xplus);

  Res.SetFlux(y);
  Res.SetFluxSigmaMinus(y1minus);
  Res.SetFluxSigmaPlus(y1plus);
  Res.SetFlux3SigmaMinus(y3minus);
  Res.SetFlux3SigmaPlus(y3plus);


  fMapResiduals[objname] = Res;

  if(DEBUG) {
    std::cout << "After:" << std::endl;
    Res.Print();
  }

}


/**
 * \brief Add points
 * \warning Energy have to be in TeV and Flux in cm-2.s-1.TeV-1
 */
void START::MultiWaveLengthFactory::AddPoints(std::string objname, std::vector<double> x, std::vector<double> y,
				       Color_t MarkerColor, Style_t MarkerStyle, Size_t MarkerSize,
				       Color_t LineColor, Width_t LineWidth) {
  Residuals Res(x.size());
  Res.SetMarkerColor(MarkerColor);
  Res.SetMarkerStyle(MarkerStyle);
  Res.SetMarkerSize(MarkerSize);
  Res.SetLineColor(LineColor);
  Res.SetLineWidth(LineWidth);  

  for(unsigned int i(0); i<y.size(); i++) {

    DEBUG_OUT_L(2) << "Before x=:" << x[i] << " y=" << y[i] << std::endl; 

    switch(fFluxUnits) {
    case vFv:
      y[i] = y[i]*x[i]*x[i]*STARTUtils::GetTeVToErg();
      break;
    case ESquareFlux_TeV:
      y[i] = y[i]*x[i]*x[i];
      break;
    case Differential_TeV:
      y[i]/=1.;
      break;
    case Differential_GeV:
      y[i]/=STARTUtils::GetTeVToGeV();
      break;
    case Differential_MeV:
      y[i]/=STARTUtils::GetTeVToMeV();
      break;
    case Differential_keV:
      y[i]/=STARTUtils::GetTeVTokeV();
    }
  }

  for(unsigned int i(0); i<x.size(); i++) {
    switch(fEnergyUnits) {
    case TeV:
      break;
    case GeV:
      x[i]*=STARTUtils::GetTeVToGeV();
      break;
    case MeV:
      x[i]*=STARTUtils::GetTeVToMeV();
      break;
    case keV:
      x[i]*=STARTUtils::GetTeVTokeV();
      break;
    case Hz:
      x[i]*=STARTUtils::GetTeVToHz();
    }
    DEBUG_OUT_L(2) << "After x=:" << x[i] << " y=" << y[i] << std::endl; 
  }

  Res.SetEnergy(x);
  Res.SetFlux(y);

  fMapResiduals[objname] = Res;

  if(DEBUG) Res.Print();
}


/**
 * \brief Add points
 * \warning Energy have to be in TeV and Flux in cm-2.s-1.TeV-1
 */
void START::MultiWaveLengthFactory::AddPoints(std::string objname, std::vector<double> x, std::vector<double> y, 
				       std::vector<double> xerr, std::vector<double> yerr,
				       Color_t MarkerColor, Style_t MarkerStyle, Size_t MarkerSize,
				       Color_t LineColor, Width_t LineWidth) {
  Residuals Res(x.size());
  Res.SetMarkerColor(MarkerColor);
  Res.SetMarkerStyle(MarkerStyle);
  Res.SetMarkerSize(MarkerSize);
  Res.SetLineColor(LineColor);
  Res.SetLineWidth(LineWidth);

  std::vector<double> xminus;
  std::vector<double> xplus;
  std::vector<double> y1minus;
  std::vector<double> y1plus;

  for(unsigned int i(0); i<x.size(); i++) {
    xminus[i]=x[i]-xerr[i];
    xplus[i]=x[i]+xerr[i];
    y1minus[i]=y[i]-yerr[i];
    y1plus[i]=y[i]+yerr[i];
  }

  for(unsigned int i(0); i<y.size(); i++) {

    DEBUG_OUT_L(2) << "Before x=:" << x[i] << " y=" << y[i] << std::endl; 

    switch(fFluxUnits) {
    case vFv:
      y[i] = y[i]*x[i]*x[i]*STARTUtils::GetTeVToErg();
      y1minus[i] = y1minus[i]*x[i]*x[i]*STARTUtils::GetTeVToErg();
      y1plus[i] = y1plus[i]*x[i]*x[i]*STARTUtils::GetTeVToErg();
      break;
    case ESquareFlux_TeV:
      y[i] = y[i]*x[i]*x[i];
      y1minus[i] = y1minus[i]*x[i]*x[i];
      y1plus[i] = y1plus[i]*x[i]*x[i];
      break;
    case Differential_TeV:
      y[i]/=1.;
      y1minus[i]/=1.;
      y1plus[i]/=1.;
      break;
    case Differential_GeV:
      y[i]/=STARTUtils::GetTeVToGeV();
      y1minus[i]/=STARTUtils::GetTeVToGeV();
      y1plus[i]/=STARTUtils::GetTeVToGeV();
      break;
    case Differential_MeV:
      y[i]/=STARTUtils::GetTeVToMeV();
      y1minus[i]/=STARTUtils::GetTeVToMeV();
      y1plus[i]/=STARTUtils::GetTeVToMeV();
      break;
    case Differential_keV:
      y[i]/=STARTUtils::GetTeVTokeV();
      y1minus[i]/=STARTUtils::GetTeVTokeV();
      y1plus[i]/=STARTUtils::GetTeVTokeV();
    }
  }

  for(unsigned int i(0); i<x.size(); i++) {
    switch(fEnergyUnits) {
    case TeV:
      break;
    case GeV:
      x[i]*=STARTUtils::GetTeVToGeV();
      xplus[i]*=STARTUtils::GetTeVToGeV();
      xminus[i]*=STARTUtils::GetTeVToGeV();
      break;
    case MeV:
      x[i]*=STARTUtils::GetTeVToMeV();
      xplus[i]*=STARTUtils::GetTeVToMeV();
      xminus[i]*=STARTUtils::GetTeVToMeV();
      break;
    case keV:
      x[i]*=STARTUtils::GetTeVTokeV();
      xminus[i]*=STARTUtils::GetTeVTokeV();
      xplus[i]*=STARTUtils::GetTeVTokeV();
      break;
    case Hz:
      x[i]*=STARTUtils::GetTeVToHz();
      xplus[i]*=STARTUtils::GetTeVToHz();
      xminus[i]*=STARTUtils::GetTeVToHz();
    }

    DEBUG_OUT_L(2) << "After x=:" << x[i] << " y=" << y[i] << std::endl; 

  }

  Res.SetEnergy(x);
  Res.SetEnergyMinus(xminus);
  Res.SetEnergyPlus(xplus);
  Res.SetFlux(y);
  Res.SetFluxSigmaMinus(y1minus);
  Res.SetFluxSigmaPlus(y1plus);
  
  fMapResiduals[objname] = Res;

  if(DEBUG) Res.Print();
}


/**
 * \brief Add points
 * \warning Energy have to be in TeV and Flux in cm-2.s-1.TeV-1
 */
void START::MultiWaveLengthFactory::AddPoints(std::string objname, std::vector<double> x, std::vector<double> y, 
				       std::vector<double> xminus, std::vector<double> xplus,
				       std::vector<double> y1minus, std::vector<double> y1plus,
				       Color_t MarkerColor, Style_t MarkerStyle, Size_t MarkerSize,
				       Color_t LineColor, Width_t LineWidth) {
  Residuals Res(x.size());
  Res.SetMarkerColor(MarkerColor);
  Res.SetMarkerStyle(MarkerStyle);
  Res.SetMarkerSize(MarkerSize);
  Res.SetLineColor(LineColor);
  Res.SetLineWidth(LineWidth);

  for(unsigned int i(0); i<y.size(); i++) {

    DEBUG_OUT_L(2) << "Before x=:" << x[i] << " y=" << y[i] << std::endl; 

    switch(fFluxUnits) {
    case vFv:
      y[i] = y[i]*x[i]*x[i]*STARTUtils::GetTeVToErg();
      y1minus[i] = y1minus[i]*x[i]*x[i]*STARTUtils::GetTeVToErg();
      y1plus[i] = y1plus[i]*x[i]*x[i]*STARTUtils::GetTeVToErg();
      break;
    case ESquareFlux_TeV:
      y[i] = y[i]*x[i]*x[i];
      y1minus[i] = y1minus[i]*x[i]*x[i];
      y1plus[i] = y1plus[i]*x[i]*x[i];
      break;
    case Differential_TeV:
      y[i]/=1.;
      y1minus[i]/=1.;
      y1plus[i]/=1.;
      break;
    case Differential_GeV:
      y[i]/=STARTUtils::GetTeVToGeV();
      y1minus[i]/=STARTUtils::GetTeVToGeV();
      y1plus[i]/=STARTUtils::GetTeVToGeV();
      break;
    case Differential_MeV:
      y[i]/=STARTUtils::GetTeVToMeV();
      y1minus[i]/=STARTUtils::GetTeVToMeV();
      y1plus[i]/=STARTUtils::GetTeVToMeV();
      break;
    case Differential_keV:
      y[i]/=STARTUtils::GetTeVTokeV();
      y1minus[i]/=STARTUtils::GetTeVTokeV();
      y1plus[i]/=STARTUtils::GetTeVTokeV();
    }
  }

  for(unsigned int i(0); i<x.size(); i++) {
    switch(fEnergyUnits) {
    case TeV:
      break;
    case GeV:
      x[i]*=STARTUtils::GetTeVToGeV();
      xplus[i]*=STARTUtils::GetTeVToGeV();
      xminus[i]*=STARTUtils::GetTeVToGeV();
      break;
    case MeV:
      x[i]*=STARTUtils::GetTeVToMeV();
      xplus[i]*=STARTUtils::GetTeVToMeV();
      xminus[i]*=STARTUtils::GetTeVToMeV();
      break;
    case keV:
      x[i]*=STARTUtils::GetTeVTokeV();
      xminus[i]*=STARTUtils::GetTeVTokeV();
      xplus[i]*=STARTUtils::GetTeVTokeV();
      break;
    case Hz:
      x[i]*=STARTUtils::GetTeVToHz();
      xplus[i]*=STARTUtils::GetTeVToHz();
      xminus[i]*=STARTUtils::GetTeVToHz();
    }

    DEBUG_OUT_L(2) << "After x=:" << x[i] << " y=" << y[i] << std::endl; 

  }

  Res.SetEnergy(x);
  Res.SetEnergyMinus(xminus);
  Res.SetEnergyPlus(xplus);
  Res.SetFlux(y);
  Res.SetFluxSigmaMinus(y1minus);
  Res.SetFluxSigmaPlus(y1plus);

  fMapResiduals[objname] = Res;

  if(DEBUG) Res.Print();
}


/**
 * \brief Add UpperLimits
 * \warning Energy have to be in TeV and Flux in cm-2.s-1.TeV-1
 */
void START::MultiWaveLengthFactory::AddUpperLimits(std::string objname, std::vector<double> x, 
					    std::vector<double> yhigh,
					    Color_t LineColor, Width_t LineWidth) {
  Residuals Res(x.size());
  Res.SetLineColor(LineColor);
  Res.SetLineWidth(LineWidth);

  std::vector<bool> isupperlimit;
  for(unsigned int i(0); i<x.size(); i++) isupperlimit.push_back(true); // upperlimits

  for(unsigned int i(0); i<yhigh.size(); i++) {

    DEBUG_OUT_L(2) << "Before x=:" << x[i] << " y=" << yhigh[i] << std::endl; 

    switch(fFluxUnits) {
    case vFv:
      yhigh[i] = yhigh[i]*x[i]*x[i]*STARTUtils::GetTeVToErg();
      break;
    case ESquareFlux_TeV:
      yhigh[i] = yhigh[i]*x[i]*x[i];
      break;
    case Differential_TeV:
      yhigh[i]/=1.;
      break;
    case Differential_GeV:
      yhigh[i]/=STARTUtils::GetTeVToGeV();
      break;
    case Differential_MeV:
      yhigh[i]/=STARTUtils::GetTeVToMeV();
      break;
    case Differential_keV:
      yhigh[i]/=STARTUtils::GetTeVTokeV();
    }
  }

  for(unsigned int i(0); i<x.size(); i++) {
    switch(fEnergyUnits) {
    case TeV:
      break;
    case GeV:
      x[i]*=STARTUtils::GetTeVToGeV();
      break;
    case MeV:
      x[i]*=STARTUtils::GetTeVToMeV();
      break;
    case keV:
      x[i]*=STARTUtils::GetTeVTokeV();
      break;
    case Hz:
      x[i]*=STARTUtils::GetTeVToHz();
    }

    DEBUG_OUT_L(2) << "After x=:" << x[i] << " y=" << yhigh[i] << std::endl; 

  }

  Res.SetEnergy(x);
  Res.SetFlux3SigmaPlus(yhigh);
  Res.SetIsUpperLimits(isupperlimit);

  fMapResiduals[objname] = Res;

  if(DEBUG) Res.Print();
}


/**
 * \brief Add Butterfly
 * \warning Energy have to be in TeV and Flux in cm-2.s-1.TeV-1
 */
void START::MultiWaveLengthFactory::AddButterfly(std::string objname, TPolyLine &Butterfly,
					  Color_t FillColor, Style_t FillStyle,
					  Color_t LineColor, Width_t LineWidth) {

  std::vector<double> x,y;

  for(int i(0); i<Butterfly.GetN(); i++) {
    x.push_back(Butterfly.GetX()[i]);
    y.push_back(Butterfly.GetY()[i]);
  }

  fMapButterfly[objname] = new TPolyLine(Butterfly);
  fMapButterfly[objname]->SetFillColor(FillColor);
  fMapButterfly[objname]->SetFillStyle(FillStyle);
  fMapButterfly[objname]->SetLineColor(LineColor);
  fMapButterfly[objname]->SetLineWidth(LineWidth);

  for(unsigned int i(0); i<x.size(); i++) {

    switch(fFluxUnits) {
    case vFv:
      y[i] = y[i]*x[i]*x[i]*STARTUtils::GetTeVToErg();
      break;
    case ESquareFlux_TeV:
      y[i] = y[i]*x[i]*x[i];
      break;
    case Differential_TeV:
      y[i]/=1.;
      break;
    case Differential_GeV:
      y[i]/=STARTUtils::GetTeVToGeV();
      break;
    case Differential_MeV:
      y[i]/=STARTUtils::GetTeVToMeV();
      break;
    case Differential_keV:
      y[i]/=STARTUtils::GetTeVTokeV();
    }
    
    switch(fEnergyUnits) {
    case TeV:
      break;
    case GeV:
      x[i]*=STARTUtils::GetTeVToGeV();
      break;
    case MeV:
      x[i]*=STARTUtils::GetTeVToMeV();
      break;
    case keV:
      x[i]*=STARTUtils::GetTeVTokeV();
      break;
    case Hz:
      x[i]*=STARTUtils::GetTeVToHz();
    }

    DEBUG_OUT_L(2) << "After x=:" << x[i] << " y=" << y[i] << std::endl; 

    fMapButterfly[objname]->SetPoint(i,x[i],y[i]);

  }
  
}

void START::MultiWaveLengthFactory::AddButterflyJLK(std::string objname, TPolyLine &Butterfly,
					     Color_t FillColor, Style_t FillStyle,
					     Color_t LineColor, Width_t LineWidth) {
  
  std::vector<double> x,y;
  
  for(int i(0); i<Butterfly.GetN(); i++) {
    x.push_back(Butterfly.GetX()[i]);
    y.push_back(Butterfly.GetY()[i]);
  }
  
  fMapButterfly[objname] = new TPolyLine(Butterfly);
  fMapButterfly[objname]->SetFillColor(FillColor);
  fMapButterfly[objname]->SetFillStyle(FillStyle);
  fMapButterfly[objname]->SetLineColor(LineColor);
  fMapButterfly[objname]->SetLineWidth(LineWidth);

  for(unsigned int i(0); i<x.size(); i++) {
    
    switch(fEnergyUnits) {
    case TeV:
      break;
    case GeV:
      x[i]*=STARTUtils::GetTeVToGeV();
      break;
    case MeV:
      x[i]*=STARTUtils::GetTeVToMeV();
      break;
    case keV:
      x[i]*=STARTUtils::GetTeVTokeV();
      break;
    case Hz:
      x[i]*=STARTUtils::GetTeVToHz();
    }

    fMapButterfly[objname]->SetPoint(i,x[i],y[i]);
    
  }
  
}

/**
 * \brief Add Butterfly
 * \warning Energy have to be in TeV and Flux in cm-2.s-1.TeV-1
 */
void START::MultiWaveLengthFactory::AddButterfly(std::string objname, std::vector<double> x, std::vector<double> y,
					  Color_t FillColor, Style_t FillStyle, Color_t LineColor, Width_t LineWidth) {
  fMapButterfly[objname] = new TPolyLine();
  fMapButterfly[objname]->SetFillColor(FillColor);
  fMapButterfly[objname]->SetFillStyle(FillStyle);
  fMapButterfly[objname]->SetLineColor(LineColor);
  fMapButterfly[objname]->SetLineWidth(LineWidth);

  for(unsigned int i(0); i<x.size(); i++) {

    switch(fFluxUnits) {
    case vFv:
      y[i] = y[i]*x[i]*x[i]*STARTUtils::GetTeVToErg();
      break;
    case ESquareFlux_TeV:
      y[i] = y[i]*x[i]*x[i];
      break;
    case Differential_TeV:
      y[i]/=1.;
      break;
    case Differential_GeV:
      y[i]/=STARTUtils::GetTeVToGeV();
      break;
    case Differential_MeV:
      y[i]/=STARTUtils::GetTeVToMeV();
      break;
    case Differential_keV:
      y[i]/=STARTUtils::GetTeVTokeV();
    }
    
    switch(fEnergyUnits) {
    case TeV:
      break;
    case GeV:
      x[i]*=STARTUtils::GetTeVToGeV();
      break;
    case MeV:
      x[i]*=STARTUtils::GetTeVToMeV();
      break;
    case keV:
      x[i]*=STARTUtils::GetTeVTokeV();
      break;
    case Hz:
      x[i]*=STARTUtils::GetTeVToHz();
    }

    DEBUG_OUT_L(2) << "After x=:" << x[i] << " y=" << y[i] << std::endl; 

    fMapButterfly[objname]->SetPoint(i,x[i],y[i]);

  }

}


/**
 * \brief Add TF1
 * \warning Should be in appropriate units (energy anf flux given in constructor)
 */
void START::MultiWaveLengthFactory::AddTF1(std::string objname, TF1 &f, Color_t LineColor,Style_t LineStyle, Width_t LineWidth) {
  TF1 *tf1 = new TF1(f);
  tf1->SetLineColor(LineColor);
  tf1->SetLineStyle(LineStyle);
  tf1->SetLineWidth(LineWidth);
  tf1->SetNpx(1000);
  fMapTF1[objname] = tf1;
}


/**
 * \brief Return minimum and maximum of the flux
 */
std::pair<double,double> START::MultiWaveLengthFactory::GetFluxRange() {
  double fmin(0.),fmax(0.);

  std::vector<double> vfmin,vfmax;

  // Butterflies
  for(std::map<std::string,TPolyLine*>::iterator butt=fMapButterfly.begin(); butt!=fMapButterfly.end(); ++butt) {

    TPolyLine *pol = new TPolyLine(*(butt->second));
    unsigned int n(pol->GetN());
    for(unsigned int i(0); i<n; i++) {
      vfmin.push_back(pol->GetY()[i]);
      vfmax.push_back(pol->GetY()[i]);
    }
    delete pol; pol=0;
  }

  // TF1
  for(std::map<std::string,TF1*>::iterator tf1=fMapTF1.begin(); tf1!=fMapTF1.end(); ++tf1) {
    vfmin.push_back(tf1->second->GetMinimum());
    vfmax.push_back(tf1->second->GetMaximum());
    DEBUG_OUT << "tf1_" << tf1->second->GetName() << " min=" << vfmin.back() << " max=" << vfmax.back() << std::endl;
  }

  // Points
  for(std::map<std::string,Residuals>::const_iterator CRes=fMapResiduals.begin(); 
      CRes!=fMapResiduals.end(); ++CRes) {
    Residuals Res = CRes->second;
    std::vector<double> flux = Res.GetFlux();
    std::vector<double> flux1minus = Res.GetFluxSigmaMinus();
    std::vector<double> flux1plus = Res.GetFluxSigmaPlus();
    std::vector<double> flux3plus = Res.GetFlux3SigmaPlus();
    std::vector<bool> isupperlimit = Res.GetIsUpperLimits();
    for(unsigned int i(0); i<flux.size(); i++) {
      if(flux[i]>0.) {
	vfmin.push_back(flux[i]);
	vfmax.push_back(flux[i]);
      }
      if(flux1minus[i]>0.) vfmin.push_back(flux1minus[i]);
      if(isupperlimit[i] && flux3plus[i]>0.) vfmax.push_back(flux3plus[i]);
      if(isupperlimit[i] && flux3plus[i]>0.) vfmin.push_back(flux3plus[i]*0.05);
      if(flux1plus[i]>0.) vfmax.push_back(flux1plus[i]);
    }
  }

  std::vector<double>::iterator min,max;
  min = min_element(vfmin.begin(), vfmin.end());
  max = max_element(vfmax.begin(), vfmax.end());

  fmin = *min;
  fmax = *max;

  for(unsigned int i(0); i< vfmin.size(); i++) {
    DEBUG_OUT_L(2) << "fmin=" << vfmin[i] << " fmax=" << vfmax[i] << std::endl; 
  }

  DEBUG_OUT << "fmin=" << fmin << " fmax=" << fmax << std::endl;

  return (std::make_pair(fmin,fmax));
}

/**
 * \brief Return minimum and maximum of the energy
 */
std::pair<double,double> START::MultiWaveLengthFactory::GetEnergyRange() {
  double emin(0.),emax(0.);

  std::vector<double> vemin,vemax;

  for(std::map<std::string,TPolyLine*>::iterator butt=fMapButterfly.begin(); butt!=fMapButterfly.end(); ++butt) {

    TPolyLine *pol = new TPolyLine(*(butt->second));
    unsigned int n(pol->GetN());
    for(unsigned int i(0); i<n; i++) {
      vemin.push_back(pol->GetX()[i]);
      vemax.push_back(pol->GetX()[i]);
    }

  }


  for(std::map<std::string,TF1*>::iterator tf1=fMapTF1.begin(); tf1!=fMapTF1.end(); ++tf1) {
    vemin.push_back(tf1->second->GetXmin());
    vemax.push_back(tf1->second->GetXmax());
  }

  for(std::map<std::string,Residuals>::const_iterator CRes=fMapResiduals.begin(); 
      CRes!=fMapResiduals.end(); ++CRes) {
    Residuals Res = CRes->second;
    std::vector<double> energy = Res.GetEnergy();
    std::vector<double> energyminus = Res.GetEnergyMinus();
    std::vector<double> energyplus = Res.GetEnergyPlus();
    for(unsigned int i(0); i<energy.size(); i++) {
      if(energy[i]>0.) {
	vemin.push_back(energy[i]);
	vemax.push_back(energy[i]);
      }
      if(energyminus[i]>0.) vemin.push_back(energyminus[i]);
      if(energyplus[i]>0.) vemax.push_back(energyplus[i]);
    }
  }

  std::vector<double>::iterator min,max;
  min = min_element(vemin.begin(), vemin.end());
  max = max_element(vemax.begin(), vemax.end());

  emin = *min;
  emax = *max;

  DEBUG_OUT << "emin=" << emin << " emax=" << emax << std::endl;

  return (std::make_pair(emin,emax));
}


