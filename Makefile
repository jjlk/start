## Load ROOT's Makefile 

include Makefile.arch

## Definitions

DOLIBRARY = yes

CURDIR = $(shell pwd)

INCDIR = $(shell pwd)/include
SRCDIR = $(shell pwd)/src
OUTDIR = $(shell pwd)/out
BINDIR = $(shell pwd)/bin
LIBDIR = $(shell pwd)/lib
DOCDIR = $(shell pwd)/doxy

DIRLIBROOT = $(shell root-config --libdir)

GSLCFLAGS = $(shell gsl-config --cflags)
GSLLIBS = $(shell gsl-config --libs)

ifneq "$(wildcard /usr/lib64/libcfitsio.so)" ""
CFITSIOFLAGS= -I/usr/include/cfitsio/
CFITSIOLIBS= -L/usr/lib64/ -lcfitsio
else
CFITSIOFLAGS=$(shell pkg-config cfitsio --cflags)
CFITSIOLIBS=$(shell pkg-config cfitsio --libs)
endif

CXXFLAGS += -I$(INCDIR) $(GSLCFLAGS) $(CFITSIOFLAGS)
EXPLLINKLIBS += $(GSLLIBS) -L$(DIRLIBROOT) -lMathMore -lMinuit2 $(CFITSIOLIBS)
LIBS += $(GSLLIBS) -L$(DIRLIBROOT) -lMathMore -lMinuit2 $(CFITSIOLIBS)

EXEC = $(BINDIR)/startfitcompile
EXFITS = $(BINDIR)/MakeHessRSPfiles
EXFITS2 = $(BINDIR)/MakeHessRSPfiles2

EXEC = 
EXFITS =
EXFITS2 =

LIBRARY = $(LIBDIR)/libstart.so

DIRECTORIES = $(LIBDIR) $(OUTDIR) $(BINDIR)

INCS = 	$(INCDIR)/EnergyBin.hh $(INCDIR)/Band.hh $(INCDIR)/ComputeResults.hh $(INCDIR)/Config.hh \
	$(INCDIR)/DataSummary.hh $(INCDIR)/FCNLikelihood.hh $(INCDIR)/HandleResolArea.hh \
	$(INCDIR)/Hypothesis.hh $(INCDIR)/BandsFactory.hh $(INCDIR)/MinimizeFactory.hh \
	$(INCDIR)/MonteCarlo.hh $(INCDIR)/ResidualsFactory.hh $(INCDIR)/STARTUtils.hh \
	$(INCDIR)/PlotFactory.hh $(INCDIR)/PowerLaw.hh $(INCDIR)/ExpoCutOffPowerLaw.hh \
	$(INCDIR)/LogParabolic.hh $(INCDIR)/BrokenPowerLaw.hh $(INCDIR)/SuperExpoCutOffPowerLaw.hh \
	$(INCDIR)/SmoothBrokenPowerLaw.hh $(INCDIR)/SumHypothesis.hh $(INCDIR)/LightCurveFactory.hh \
	$(INCDIR)/TimeBin.hh $(INCDIR)/Residuals.hh $(INCDIR)/TimeBinVector.hh $(INCDIR)/MultiWaveLengthFactory.hh \
	$(INCDIR)/STARTFITSUtils.hh $(INCDIR)/Event.hh 

SRCS =  $(SRCDIR)/EnergyBin.C $(SRCDIR)/Band.C $(SRCDIR)/ComputeResults.C $(SRCDIR)/Config.C \
	$(SRCDIR)/DataSummary.C $(SRCDIR)/FCNLikelihood.C $(SRCDIR)/HandleResolArea.C \
	$(SRCDIR)/Hypothesis.C $(SRCDIR)/BandsFactory.C $(SRCDIR)/MinimizeFactory.C \
	$(SRCDIR)/MonteCarlo.C $(SRCDIR)/ResidualsFactory.C $(SRCDIR)/STARTUtils.C \
	$(SRCDIR)/PlotFactory.C $(SRCDIR)/PowerLaw.C $(SRCDIR)/ExpoCutOffPowerLaw.C \
	$(SRCDIR)/LogParabolic.C $(SRCDIR)/BrokenPowerLaw.C $(SRCDIR)/SuperExpoCutOffPowerLaw.C \
	$(SRCDIR)/SmoothBrokenPowerLaw.C $(SRCDIR)/SumHypothesis.C $(SRCDIR)/LightCurveFactory.C \
	$(SRCDIR)/TimeBin.C $(SRCDIR)/TimeBinVector.C $(SRCDIR)/Residuals.C $(SRCDIR)/MultiWaveLengthFactory.C \
	$(SRCDIR)/STARTFITSUtils.C $(SRCDIR)/Event.C 

OBJS =  $(OUTDIR)/Config.o $(OUTDIR)/StringTools.o $(OUTDIR)/EnergyBin.o $(OUTDIR)/Band.o \
	$(OUTDIR)/BandsFactory.o $(OUTDIR)/FCNLikelihood.o $(OUTDIR)/ComputeResults.o \
	$(OUTDIR)/MinimizeFactory.o $(OUTDIR)/HandleResolArea.o $(OUTDIR)/DataSummary.o \
	$(OUTDIR)/MonteCarlo.o $(OUTDIR)/Hypothesis.o $(OUTDIR)/ResidualsFactory.o \
	$(OUTDIR)/STARTUtils.o $(OUTDIR)/PlotFactory.o $(OUTDIR)/PowerLaw.o \
	$(OUTDIR)/ExpoCutOffPowerLaw.o $(OUTDIR)/LogParabolic.o $(OUTDIR)/BrokenPowerLaw.o \
	$(OUTDIR)/SuperExpoCutOffPowerLaw.o $(OUTDIR)/SmoothBrokenPowerLaw.o $(OUTDIR)/SumHypothesis.o \
	$(OUTDIR)/LightCurveFactory.o $(OUTDIR)/TimeBin.o $(OUTDIR)/TimeBinVector.o $(OUTDIR)/Residuals.o \
	$(OUTDIR)/MultiWaveLengthFactory.o $(OUTDIR)/STARTFITSUtils.o $(OUTDIR)/Event.o 

ifeq ($(DOLIBRARY),yes)
CXXFLAGS += -DCONSTRUCT_LIBRARY
OBJS += $(OUTDIR)/STARTDict.o
endif

ifeq ($(DOLIBRARY),yes)
all: $(DIRECTORIES) $(LIBRARY) $(EXEC) $(EXFITS) $(EXFITS2) $(AREA2DFITS)
else 
all: $(DIRECTORIES) $(EXEC)
endif

## directories
$(DIRECTORIES):
	@mkdir -p $@

## Library
$(LIBRARY): $(OBJS)
	@echo "Generating library $@..."
ifeq ($(PLATFORM),macosx)
	$(LD) $(SOFLAGS)$@ $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
else
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
endif
	@echo "$@ done!"

## Dictionaries: 
$(OUTDIR)/STARTDict.C: $(INCS) $(INCDIR)/STARTLinkDef.hh
	@echo "Generating dictionary $@..."
	$(ROOTCINT) -f $@ -c -DCONSTRUCT_LIBRARY $(CXXFLAGS) -p $^
	@echo "$@ done!"

## executable
$(EXEC): $(OBJS) 
ifeq ($(DOLIBRARY),yes)
	@echo "Generating $@ with library..."
	$(LD) $(LDFLAGS) $^ $(LIBS) $(LIBRARY) $(OutPutOpt)$@
	@echo "$@ done!"
else
	@echo "Generating $@ without library..."
	$(LD) $(LDFLAGS) -o $(EXEC) $^ $(LIBS) $(OutPutOpt)$@
	@echo "$@ done!"
endif

# ## executable
# $(EXFITS): $(OBJS) $(OUTDIR)/MakeHessRSPfiles.o 
# 	$(LD) $(LDFLAGS) -o $(EXFITS) $^ $(LIBS) $(OutPutOpt)$@
# 	@echo "$@ done!"


# ## sources
# $(OUTDIR)/MakeHessRSPfiles.o: $(SRCS) $(INCS) $(SRCDIR)/MakeHessRSPfiles.C
# 	$(LD) -c $(SRCDIR)/MakeHessRSPfiles.C $(CXXFLAGS)
# 	@mv *.o $(OUTDIR)
# 	@echo "$@ done!"

# ## executable
# $(EXFITS2): $(OBJS) $(OUTDIR)/MakeHessRSPfiles2.o 
# 	$(LD) $(LDFLAGS) -o $(EXFITS2) $^ $(LIBS) $(OutPutOpt)$@
# 	@echo "$@ done!"


# ## sources
# $(OUTDIR)/MakeHessRSPfiles2.o: $(SRCS) $(INCS) $(SRCDIR)/MakeHessRSPfiles2.C
# 	$(LD) -c $(SRCDIR)/MakeHessRSPfiles2.C $(CXXFLAGS)
# 	@mv *.o $(OUTDIR)
# 	@echo "$@ done!"

# ## sources
# $(OUTDIR)/main.o: $(SRCS) $(INCS) $(SRCDIR)/main.C
# 	$(LD) -c $(SRCDIR)/main.C $(CXXFLAGS)
# 	@mv *.o $(OUTDIR)
# 	@echo "$@ done!"

$(OUTDIR)/Config.o: $(SRCDIR)/Config.C $(INCDIR)/Config.hh \
	$(SRCDIR)/StringTools.C $(INCDIR)/StringTools.hh
	$(LD) -c $(SRCDIR)/Config.C $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/STARTUtils.o : $(SRCDIR)/STARTUtils.C $(INCDIR)/STARTUtils.hh
	$(LD) -c $(SRCDIR)/STARTUtils.C $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/STARTFITSUtils.o : $(SRCDIR)/STARTFITSUtils.C $(INCDIR)/STARTFITSUtils.hh
	$(LD) -c $(SRCDIR)/STARTFITSUtils.C $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/StringTools.o: $(SRCDIR)/StringTools.C $(INCDIR)/StringTools.hh
	$(LD) -c $(SRCDIR)/StringTools.C $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/EnergyBin.o: $(SRCDIR)/EnergyBin.C $(INCDIR)/EnergyBin.hh \
	$(SRCDIR)/Event.C $(INCDIR)/Event.hh
	$(LD) -c $(SRCDIR)/EnergyBin.C $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/Band.o: $(SRCDIR)/Band.C $(INCDIR)/Band.hh \
	$(SRCDIR)/EnergyBin.C $(INCDIR)/EnergyBin.hh \
	$(SRCDIR)/STARTUtils.C $(INCDIR)/STARTUtils.hh
	$(LD) -c $(SRCDIR)/Band.C $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/BandsFactory.o: $(SRCDIR)/BandsFactory.C $(INCDIR)/BandsFactory.hh \
	$(SRCDIR)/Band.C $(INCDIR)/Band.hh $(SRCDIR)/Config.C $(INCDIR)/Config.hh \
	$(SRCDIR)/Event.C $(INCDIR)/Event.hh
	$(LD) -c $(SRCDIR)/BandsFactory.C $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/FCNLikelihood.o: $(SRCDIR)/FCNLikelihood.C $(INCDIR)/FCNLikelihood.hh \
	$(SRCDIR)/Band.C $(INCDIR)/Band.hh \
	$(SRCDIR)/Hypothesis.C $(INCDIR)/Hypothesis.hh \
	$(SRCDIR)/STARTUtils.C $(INCDIR)/STARTUtils.hh
	$(LD) -c $(SRCDIR)/FCNLikelihood.C $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/ComputeResults.o: $(SRCDIR)/ComputeResults.C $(INCDIR)/ComputeResults.hh \
	$(SRCDIR)/Band.C $(INCDIR)/Band.hh $(SRCDIR)/Hypothesis.C \
	$(INCDIR)/Hypothesis.hh $(INCDIR)/GSLError.h 
	$(LD) -c $(SRCDIR)/ComputeResults.C $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/Hypothesis.o: $(SRCDIR)/Hypothesis.C $(INCDIR)/Hypothesis.hh \
	$(SRCDIR)/Band.C $(INCDIR)/Band.hh $(SRCDIR)/TimeBinVector.C \
	$(INCDIR)/TimeBinVector.hh $(SRCDIR)/Residuals.C $(INCDIR)/Residuals.hh 
	$(LD) -c $(SRCDIR)/Hypothesis.C $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/MinimizeFactory.o: $(SRCDIR)/MinimizeFactory.C $(INCDIR)/MinimizeFactory.hh \
	$(SRCDIR)/FCNLikelihood.C $(INCDIR)/FCNLikelihood.hh $(SRCDIR)/ComputeResults.C \
	$(INCDIR)/ComputeResults.hh $(SRCDIR)/STARTUtils.C $(INCDIR)/STARTUtils.hh \
	$(SRCDIR)/Hypothesis.C $(INCDIR)/Hypothesis.hh $(INCDIR)/SumHypothesis.hh \
	$(SRCDIR)/SumHypothesis.C  $(SRCDIR)/SumHypothesis.C
	$(LD) -c $(SRCDIR)/MinimizeFactory.C $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/DataSummary.o: $(SRCDIR)/DataSummary.C $(INCDIR)/DataSummary.hh \
	$(SRCDIR)/Band.C $(INCDIR)/Band.hh $(SRCDIR)/Config.C $(INCDIR)/Config.hh \
	$(SRCDIR)/MonteCarlo.C $(INCDIR)/MonteCarlo.hh
	$(LD) -c $(SRCDIR)/DataSummary.C $(CXXFLAGS) 
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/HandleResolArea.o: $(SRCDIR)/HandleResolArea.C $(INCDIR)/HandleResolArea.hh \
	$(SRCDIR)/Band.C $(INCDIR)/Band.hh $(SRCDIR)/MonteCarlo.C $(INCDIR)/MonteCarlo.hh \
	$(SRCDIR)/DataSummary.C $(INCDIR)/DataSummary.hh $(SRCDIR)/ComputeResults.C \
	$(INCDIR)/ComputeResults.hh
	$(LD) -c $(SRCDIR)/HandleResolArea.C $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/MonteCarlo.o: $(SRCDIR)/MonteCarlo.C $(INCDIR)/MonteCarlo.hh \
	$(SRCDIR)/STARTUtils.C $(INCDIR)/STARTUtils.hh
	$(LD) -c $(SRCDIR)/MonteCarlo.C $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/ResidualsFactory.o: $(SRCDIR)/ResidualsFactory.C $(INCDIR)/ResidualsFactory.hh \
	$(SRCDIR)/Band.C $(INCDIR)/Band.hh \
	$(SRCDIR)/Hypothesis.C $(INCDIR)/Hypothesis.hh
	$(LD) -c $(SRCDIR)/ResidualsFactory.C $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/Residuals.o: $(SRCDIR)/Residuals.C $(INCDIR)/Residuals.hh
	$(LD) -c $(SRCDIR)/Residuals.C $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/TimeBin.o: $(SRCDIR)/TimeBin.C $(INCDIR)/TimeBin.hh \
	$(SRCDIR)/STARTUtils.C $(INCDIR)/STARTUtils.hh
	$(LD) -c $(SRCDIR)/TimeBin.C $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/TimeBinVector.o: $(SRCDIR)/TimeBinVector.C $(INCDIR)/TimeBinVector.hh \
	$(SRCDIR)/TimeBin.C $(INCDIR)/TimeBin.hh
	$(LD) -c $(SRCDIR)/TimeBinVector.C $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/LightCurveFactory.o: $(SRCDIR)/LightCurveFactory.C $(INCDIR)/LightCurveFactory.hh \
	$(SRCDIR)/TimeBinVector.C $(INCDIR)/TimeBinVector.hh $(SRCDIR)/Config.C \
	$(INCDIR)/Config.hh $(SRCDIR)/Hypothesis.C $(INCDIR)/Hypothesis.hh \
	$(SRCDIR)/Event.C $(INCDIR)/Event.hh
	$(LD) -c $(SRCDIR)/LightCurveFactory.C $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/PlotFactory.o: $(SRCDIR)/PlotFactory.C $(INCDIR)/PlotFactory.hh \
	$(SRCDIR)/Band.C $(INCDIR)/Band.hh $(SRCDIR)/ResidualsFactory.C \
	$(INCDIR)/ResidualsFactory.hh \
	$(SRCDIR)/Hypothesis.C $(INCDIR)/Hypothesis.hh \
	$(SRCDIR)/MultiWaveLengthFactory.C $(INCDIR)/MultiWaveLengthFactory.hh \
	$(SRCDIR)/Residuals.C $(INCDIR)/Residuals.hh \
	$(SRCDIR)/TimeBinVector.C $(INCDIR)/TimeBinVector.hh	
	$(LD) -c $(SRCDIR)/PlotFactory.C $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/STARTDict.o: $(OUTDIR)/STARTDict.C
	@echo "Generating object dictionary $@..."
	$(CXX) -c $< $(CXXFLAGS)
	@mv *.o $(OUTDIR)/
	@echo "$@ done!"

$(OUTDIR)/PowerLaw.o: $(SRCDIR)/PowerLaw.C  $(INCDIR)/PowerLaw.hh \
	$(SRCDIR)/Hypothesis.C $(INCDIR)/Hypothesis.hh 
	$(LD) -c $< $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/ExpoCutOffPowerLaw.o: $(SRCDIR)/ExpoCutOffPowerLaw.C  $(INCDIR)/ExpoCutOffPowerLaw.hh \
	$(SRCDIR)/Hypothesis.C $(INCDIR)/Hypothesis.hh 
	$(LD) -c $< $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/LogParabolic.o: $(SRCDIR)/LogParabolic.C  $(INCDIR)/LogParabolic.hh \
	$(SRCDIR)/Hypothesis.C $(INCDIR)/Hypothesis.hh 
	$(LD) -c $< $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/BrokenPowerLaw.o: $(SRCDIR)/BrokenPowerLaw.C  $(INCDIR)/BrokenPowerLaw.hh \
	$(SRCDIR)/Hypothesis.C $(INCDIR)/Hypothesis.hh 
	$(LD) -c $< $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/SmoothBrokenPowerLaw.o: $(SRCDIR)/SmoothBrokenPowerLaw.C  $(INCDIR)/SmoothBrokenPowerLaw.hh \
	$(SRCDIR)/Hypothesis.C $(INCDIR)/Hypothesis.hh 
	$(LD) -c $< $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/SuperExpoCutOffPowerLaw.o: $(SRCDIR)/SuperExpoCutOffPowerLaw.C \
	$(INCDIR)/SuperExpoCutOffPowerLaw.hh $(SRCDIR)/Hypothesis.C \
	$(INCDIR)/Hypothesis.hh 
	$(LD) -c $< $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/SumHypothesis.o: $(SRCDIR)/SumHypothesis.C $(INCDIR)/SumHypothesis.hh \
	$(SRCDIR)/Hypothesis.C $(INCDIR)/Hypothesis.hh 
	$(LD) -c $< $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/MultiWaveLengthFactory.o: $(SRCDIR)/MultiWaveLengthFactory.C \
	$(INCDIR)/MultiWaveLengthFactory.hh $(SRCDIR)/STARTUtils.C \
	$(INCDIR)/STARTUtils.hh $(SRCDIR)/Residuals.C $(INCDIR)/Residuals.hh
	$(LD) -c $< $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

$(OUTDIR)/Event.o: $(SRCDIR)/Event.C $(INCDIR)/Event.hh
	$(LD) -c $< $(CXXFLAGS)
	@mv *.o $(OUTDIR)
	@echo "$@ done!"

.PHONY: clean cleandoc distclean doc symlink hessinstall

doc: 
	@echo "Generating documentation with doxygen..."
	@rm -rf $(DOCDIR)
	@mkdir $(DOCDIR)
	doxygen Doxyfile
	@echo "Doc done! You can find it in ./doxy/html/index.html"

cleandoc:
	rm -rf $(DOCDIR)

clean: 
	rm -rf $(OUTDIR)/*.o
	rm -rf $(OUTDIR)/STARTDict.C
	rm -rf $(OUTDIR)/STARTDict.h
	rm -rf $(BINDIR)/*

distclean: clean cleandoc
	rm -rf $(EXEC)
	rm -rf $(LIBDIR)/*.so

symlinkdir:
	$(shell mkdir -p $(CURDIR)/../include/)
	$(shell mkdir -p $(CURDIR)/../lib/)

symlink: symlinkdir
	$(shell cd $(CURDIR)/../include/; ln -sf $(INCDIR) start)
	$(shell cd $(CURDIR)/../lib/; ln -sf $(LIBRARY) libstart.so)
	$(shell cd $(CURDIR))

hessinstall: all symlink

