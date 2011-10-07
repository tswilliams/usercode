import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


##                      _              _       
##   ___ ___  _ __  ___| |_ __ _ _ __ | |_ ___ 
##  / __/ _ \| '_ \/ __| __/ _` | '_ \| __/ __|
## | (_| (_) | | | \__ \ || (_| | | | | |_\__ \
##  \___\___/|_| |_|___/\__\__,_|_| |_|\__|___/
##                                              
################################################


isMC = True
## InputFileName = "testNewWrite.root"
InputFileName = "/home/ppd/nnd85574/Work/BstdZee/TagAndProbeEffi/CMSSW_4_2_4_patch1/src/PhysicsTools/TagAndProbe/test/testTandPEffiTree.root"
OutputFilePrefix = "efficiency-data-"


HLTDef = "probe_passingHLT"
PDFName = "pdfSignalPlusBackground"

if isMC:
    InputFileName = "/home/ppd/nnd85574/Work/BstdZee/TagAndProbeEffi/CMSSW_4_2_4_patch1/src/PhysicsTools/TagAndProbe/test/testTandPEffiTree.root"
    PDFName = "pdfSignalPlusBackground"
    OutputFilePrefix = "efficiency-mc-"
################################################

#Lists of parameter bins for efficiency calculations ...
Effi_EtEtaBins = cms.PSet(
    probe_sc_et = cms.vdouble( 25, 30, 35, 40, 45, 50, 200 ),
    probe_sc_eta = cms.vdouble( -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5 )
)


#### For data: except for HLT step
EfficiencyBinningSpecification = cms.PSet(
    #specifies what unbinned variables to include in the dataset, the mass is needed for the fit
    UnbinnedVariables = cms.vstring("mass"),
    #specifies the binning of parameters
    BinnedVariables = cms.PSet(Effi_EtEtaBins)#,
    #first string is the default followed by binRegExp - PDFname pairs
    #BinToPDFmap = cms.vstring(PDFName)
)

#### For MC truth: do truth matching
EfficiencyBinningSpecificationMC = cms.PSet(
    UnbinnedVariables = cms.vstring("mass"),
    BinnedVariables = cms.PSet(
    probe_gsfEle_et = cms.vdouble( 25, 30, 35, 40, 45, 50, 200 ),
    probe_gsfEle_eta = cms.vdouble( -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5 ),
    mcTrue = cms.vstring("true")
    )#,
    #BinToPDFmap = cms.vstring()  
)

#### For HLT step: just do cut & count
EfficiencyBinningSpecificationHLT = cms.PSet(
    UnbinnedVariables = cms.vstring("mass"),
    BinnedVariables = cms.PSet(Effi_EtEtaBins),
    BinToPDFmap = cms.vstring()  
)

##########################################################################################
############################################################################################
if isMC:
    mcTruthModules = cms.PSet(
##         MCtruth_WP95 = cms.PSet(
##         EfficiencyBinningSpecificationMC,
##         EfficiencyCategoryAndState = cms.vstring("probe_isWP95","pass"),
##         ),
        MCtruth_WP90 = cms.PSet( EfficiencyBinningSpecificationMC, EfficiencyCategoryAndState = cms.vstring("probe_isWP90","pass"), )
    )
else:
    mcTruthModules = cms.PSet()

############################################################################################
############################################################################################
####### GsfElectron->Id / selection efficiency 
############################################################################################
############################################################################################

process.GsfElectronToId = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    # IO parameters:
    InputFileNames = cms.vstring(InputFileName),
    InputDirectoryName = cms.string("GsfElectronToId"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string(OutputFilePrefix+"GsfElectronToId.root"),
    #numbrer of CPUs to use for fitting
    NumCPU = cms.uint32(1),
    # specifies wether to save the RooWorkspace containing the data for each bin and
    # the pdf object with the initial and final state snapshots
    SaveWorkspace = cms.bool(True),
    floatShapeParameters = cms.bool(True),
    #fixVars = cms.vstring("mean"),
                                                 
    # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "60.0", "120.0", "GeV/c^{2}"),
        probe_gsfEle_et  = cms.vstring("Probe E_{T}","0","1000","GeV"),
        probe_gsfEle_eta = cms.vstring("Probe #eta","-2.5","2.5",""),
        probe_sc_et = cms.vstring("Probe E_{T}", "0", "1000", "GeV/c"),
        probe_sc_eta = cms.vstring("Probe #eta", "-2.5", "2.5", ""),                
    ),

    # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
    Categories = cms.PSet(
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
        probe_isWP90 = cms.vstring("probe_isWP90", "dummy[pass=1,fail=0]")
    ),
    # defines all the PDFs that will be available for the efficiency calculations; uses RooFit's "factory" syntax;
    # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" and "signalFractionInPassing[0.9]" are used for initial values  
    PDFs = cms.PSet(
        pdfSignalPlusBackground = cms.vstring("CBExGaussShape::signalResPass(mass, meanP[0.], sigmaP[8.5695e-04, 0., 3.],alphaP[3.8296e-04], nP[6.7489e+00], sigmaP_2[2.5849e+00], fracP[6.5704e-01])",  ### signal resolution for "pass" sample
           "CBExGaussShape::signalResFail(mass, meanF[2.0946e-01, -5., 5.], sigmaF[8.5695e-04, 0., 5.],alphaF[3.8296e-04], nF[6.7489e+00], sigmaF_2[2.5849e+00], fracF[6.5704e-01])",  ### signal resolution for "fail" sample     
           "ZGeneratorLineShape::signalPhy(mass)", ### NLO line shape
           "RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], betaPass, peakPass[90.0])",
           "RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], betaFail, peakFail[90.0])",
           "FCONV::signalPass(mass, signalPhy, signalResPass)",
           "FCONV::signalFail(mass, signalPhy, signalResFail)",     
           "efficiency[0.9,0,1]",
           "signalFractionInPassing[1.0]"),
    ),

    # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
    # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
    Efficiencies = cms.PSet( #mcTruthModules,
       WP90 = cms.PSet( EfficiencyBinningSpecification,
                   #Specifying the 'category' which contains the probe pass/fail result...
                   EfficiencyCategoryAndState = cms.vstring("probe_isWP90","pass"), )
    )
)

############################################################################################
############################################################################################
####### HLT efficiency (As yet untested in my reduced version of config file)
############################################################################################
############################################################################################

if isMC:
    HLTmcTruthModules = cms.PSet(
        MCtruth_efficiency = cms.PSet(
        EfficiencyBinningSpecificationMC,
        EfficiencyCategoryAndState = cms.vstring( HLTDef, "pass" ),
        ),    
    )
else:
    HLTmcTruthModules = cms.PSet()


EfficienciesPset = cms.PSet(
    HLTmcTruthModules,
    efficiency = cms.PSet(
    EfficiencyBinningSpecificationHLT,
    EfficiencyCategoryAndState = cms.vstring( HLTDef, "pass" ),
    ),
)

########
process.WP95ToHLT = process.GsfElectronToId.clone()
process.WP95ToHLT.InputDirectoryName = cms.string("WP95ToHLT")
process.WP95ToHLT.OutputFileName = cms.string(OutputFilePrefix+"WP95ToHLT.root")
process.WP95ToHLT.Categories = cms.PSet(
    mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),           
    probe_passingHLT = cms.vstring("probe_passingHLT", "dummy[pass=1,fail=0]"), 
    )
process.WP95ToHLT.Efficiencies = EfficienciesPset
process.WP95ToHLT.Efficiencies.efficiency.BinToPDFmap = cms.vstring()
########  
process.WP90ToHLT = process.WP95ToHLT.clone()
process.WP90ToHLT.InputDirectoryName = cms.string("WP90ToHLT")
process.WP90ToHLT.OutputFileName = cms.string(OutputFilePrefix+"WP90ToHLT.root")


###########################################################################################
##---------------------------------------------------------------------------------------##
##                              GSF -> HEEP efficiency                                   ##
##---------------------------------------------------------------------------------------##
###########################################################################################

heepInputFName = "testTandPEffiTree.root"
heepInputFDir = "heepTagProbeTree"

process.Effi_GsfToHEEP = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    # *** I/O PARAMETERS
    InputFileNames = cms.vstring(heepInputFName),
    InputDirectoryName = cms.string(heepInputFDir),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string(OutputFilePrefix + heepInputFDir + ".root"),
                                        
    ## *** BRANCHES USED IN EFFICIENCY CALC'N *** 
    # ... real variables intended for use in the efficiencies
    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "60.0", "120.0", "GeV/c^{2}"),
        probe_heepEt  = cms.vstring("Probe E_{T}^{HEEP}","0.0","1000.0","GeV"),
        probe_eta = cms.vstring("Probe #eta","-2.5","2.5",""),
        probe_scEta = cms.vstring("Probe #eta_{SC}", "-2.5", "2.5", ""),                
    ),
    # ... discrete variables intended for use in the efficiency calculations
    Categories = cms.PSet(
        #mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
        probe_passesAllHEEP = cms.vstring("probe_passesAllHEEP", "dummy[pass=1,fail=0]")
    ),
                                        
    ## *** EFFICIENCY SPECIFICATIONS ***
    # Define the binning, and cut-flag, used for each efficiency calculation
    # # N.B: Fitting-based efficiency calc'n is disabled by 
    Efficiencies = cms.PSet( #mcTruthModules,
       AllHEEP_unbinned = cms.PSet( 
          #Specify the 'category' which contains the probe pass/fail result
          EfficiencyCategoryAndState = cms.vstring("probe_passesAllHEEP","pass"), 
          #specifies what unbinned variables to include in the dataset, the mass is needed for the fit
          UnbinnedVariables = cms.vstring("mass"),
          # Specify binning of parameters
          BinnedVariables = cms.PSet( probe_scEta = cms.vdouble(-2.5, -1.56, -1.442, 1.442, 1.56, 2.5) )
          # BinToPDFmap NOT SPECIFIED so that only counting efficiency calculated
          #BinToPDFmap = cms.vstring(PDFName),
       )
    ),
                                        
    ## *** FITING PARAMETERS***
    # Number of CPUs to use for fitting
    NumCPU = cms.uint32(1),
    # specifies wether to save the RooWorkspace containing the data for each bin and
    # the pdf object with the initial and final state snapshots
    SaveWorkspace = cms.bool(True),
    floatShapeParameters = cms.bool(True),
    #fixVars = cms.vstring("mean"),
    # defines all the PDFs that will be available for the efficiency calculations; uses RooFit's "factory" syntax;
    # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" and "signalFractionInPassing[0.9]" are used for initial values  
    PDFs = cms.PSet(
        pdfSignalPlusBackground = cms.vstring("CBExGaussShape::signalResPass(mass, meanP[0.], sigmaP[8.5695e-04, 0., 3.],alphaP[3.8296e-04], nP[6.7489e+00], sigmaP_2[2.5849e+00], fracP[6.5704e-01])",  ### signal resolution for "pass" sample
           "CBExGaussShape::signalResFail(mass, meanF[2.0946e-01, -5., 5.], sigmaF[8.5695e-04, 0., 5.],alphaF[3.8296e-04], nF[6.7489e+00], sigmaF_2[2.5849e+00], fracF[6.5704e-01])",  ### signal resolution for "fail" sample     
           "ZGeneratorLineShape::signalPhy(mass)", ### NLO line shape
           "RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], betaPass, peakPass[90.0])",
           "RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], betaFail, peakFail[90.0])",
           "FCONV::signalPass(mass, signalPhy, signalResPass)",
           "FCONV::signalFail(mass, signalPhy, signalResFail)",     
           "efficiency[0.9,0,1]",
           "signalFractionInPassing[1.0]"),
    )
)



###################################################################
###################################################################
###### CMS Path
###################################################################
###################################################################

process.fit = cms.Path( process.Effi_GsfToHEEP )
