import FWCore.ParameterSet.Config as cms

process = cms.Process("UFHZZ4LAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.categories.append('UFHZZ4LAna')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag='74X_mcRun2_asymptotic_v2'

process.Timing = cms.Service("Timing",
                             summaryOnly = cms.untracked.bool(True)
                             )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

myfilelist = cms.untracked.vstring()
myfilelist.extend( [
    'root://cmsxrootd.fnal.gov//store/mc/RunIISpring15MiniAODv2/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/30000/3E964C5D-1D6E-E511-8B9A-0050560207C5.root',
    'root://cmsxrootd.fnal.gov//store/mc/RunIISpring15MiniAODv2/WminusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/30000/D8CA6B54-056F-E511-BB1A-02163E014CE3.root',
    'root://cmsxrootd.fnal.gov//store/mc/RunIISpring15MiniAODv2/WplusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/D22BEE88-C26D-E511-B330-002590A81EF0.root',
    'root://cmsxrootd.fnal.gov//store/mc/RunIISpring15MiniAODv2/ttH_HToZZ_4LFilter_M125_13TeV_powheg_JHUgen_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/80000/84F62DD7-1475-E511-9F59-009C02AB98A6.root'
]
)

process.source = cms.Source("PoolSource",fileNames = myfilelist,
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
#                            eventsToProcess = cms.untracked.VEventRange('1:264711-1:264711'),
                            inputCommands = cms.untracked.vstring('keep *',
                                                                  'drop LHEEventProduct_*_*_*',
                                                                  'drop LHERunInfoProduct_*_*_*')
                            )

process.TFileService = cms.Service("TFileService",
#                                   fileName = cms.string("Sync_v2_74X.root")
                                   fileName = cms.string("test.root")
)

# clean muons by segments 
process.boostedMuons = cms.EDProducer("PATMuonCleanerBySegments",
				     src = cms.InputTag("slimmedMuons"),
				     preselection = cms.string("track.isNonnull"),
				     passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
				     fractionOfSharedSegments = cms.double(0.499),
				     )

# Electron Calibrations
process.calibratedElectrons = cms.EDProducer("CalibratedPatElectronProducerRun2",
    electrons = cms.InputTag("slimmedElectrons","","PAT"),
    grbForestName = cms.string("gedelectron_p4combination_25ns"),
    isMC = cms.bool(True),
    isSynchronization = cms.bool(True)
)

# Electron MVA ID producers
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff']
# add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag("calibratedElectrons")

process.mvaSpring15NonTrig25nsV1 = cms.EDProducer("SlimmedElectronMvaIDProducer",
                                     mvaValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
                                     electronsCollection = cms.InputTag("calibratedElectrons"),
                                     Trig = cms.bool(False),
                                     )

# FSR Photons
process.load('UFHZZAnalysisRun2.FSRPhotons.fsrPhotons_cff')

# Jet Energy Corrections
process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")

process.jetCorrFactors = process.patJetCorrFactorsUpdated.clone(
    src = cms.InputTag("slimmedJets"),
    levels = ['L1FastJet', 
              'L2Relative', 
              'L3Absolute'],
    payload = 'AK4PFchs' ) 

process.slimmedJetsJEC = process.patJetsUpdated.clone(
    jetSource = cms.InputTag("slimmedJets"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("jetCorrFactors"))
    )

process.Ana = cms.EDAnalyzer('UFHZZ4LAna',
                              photonSrc    = cms.untracked.InputTag("slimmedPhotons"),
                              electronSrc  = cms.untracked.InputTag("mvaSpring15NonTrig25nsV1","NonTrig"),
                              muonSrc      = cms.untracked.InputTag("boostedMuons"),
                              jetSrc       = cms.untracked.InputTag("slimmedJets"),
#                              jetSrc       = cms.untracked.InputTag("slimmedJetsJEC"),
                              metSrc       = cms.untracked.InputTag("slimmedMETs"),
                              vertexSrc    = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),
                              isMC         = cms.untracked.bool(True),
                              isSignal     = cms.untracked.bool(True),
                              mH           = cms.untracked.double(125.0),
                              CrossSection = cms.untracked.double(1.0),
                              FilterEff    = cms.untracked.double(1),
                              weightEvents = cms.untracked.bool(True),
                              elRhoSrc     = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                              muRhoSrc     = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                              rhoSrcSUS    = cms.untracked.InputTag("fixedGridRhoFastjetCentralNeutral"),
                              pileupSrc     = cms.untracked.InputTag("slimmedAddPileupInfo"),
                              reweightForPU = cms.untracked.bool(True),
                              triggerSrc = cms.untracked.InputTag("TriggerResults","","HLT"),
                              triggerObjects = cms.InputTag("selectedPatTrigger"),
                              triggerList = cms.untracked.vstring(
                                            'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
                                            'HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v',
                                            'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v',
                                            'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v',
                                            'HLT_TripleMu_12_10_5_v',
                                            'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v',
                                            'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
                                            'HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v',
                                            'HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v',
                                            'HLT_Ele27_WP85_Gsf_v',
                              ),
                              verbose = cms.untracked.bool(False)              
#                              verbose = cms.untracked.bool(True)              
                             )


process.p = cms.Path(process.fsrPhotonSequence*
                     process.boostedMuons*
                     process.calibratedElectrons*
                     process.electronMVAValueMapProducer*
                     process.mvaSpring15NonTrig25nsV1*
#                     process.jetCorrFactors*
#                     process.slimmedJetsJEC*
                     process.Ana
                     )
