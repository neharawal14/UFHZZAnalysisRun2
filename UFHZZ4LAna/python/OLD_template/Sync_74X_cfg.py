import FWCore.ParameterSet.Config as cms

process = cms.Process("UFHZZ4LAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.categories.append('UFHZZ4LAna')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag='MCRUN2_74_V9::All'

process.Timing = cms.Service("Timing",
                             summaryOnly = cms.untracked.bool(True)
                             )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

myfilelist = cms.untracked.vstring()
myfilelist.extend( [
       '/store/mc/RunIISpring15DR74/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/68791C0A-3013-E511-88FD-D4AE5269F5FF.root',
       '/store/mc/RunIISpring15DR74/WplusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/60000/04BD6860-9F08-E511-8A80-842B2B1858FB.root',
       '/store/mc/RunIISpring15DR74/WminusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/70000/4A9FED55-DF0C-E511-A4B2-3417EBE6471D.root',
       '/store/mc/RunIISpring15DR74/ZH_HToZZ_4LFilter_M125_13TeV_powheg-minlo-HZJ_JHUgen_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/104B7067-0C02-E511-8FFB-0030487D07BA.root'
]
)

process.source = cms.Source("PoolSource",fileNames = myfilelist,
#                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
#                             eventsToProcess = cms.untracked.VEventRange('1:1-1:1000')
#                             eventsToProcess = cms.untracked.VEventRange('1:66718-1:66718')
#                             eventsToProcess = cms.untracked.VEventRange('1:86646-1:86646','1:312169-1:312169','1:423353-1:423353')
                            )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("Sync_74x.root")
                                   )

# clean muons by segments 
process.boostedMuons = cms.EDProducer("PATMuonCleanerBySegments",
				     src = cms.InputTag("slimmedMuons"),
				     preselection = cms.string("track.isNonnull"),
				     passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
				     fractionOfSharedSegments = cms.double(0.499),
				     )

# Electron MVA ID producer
process.mvaNonTrigV025nsPHYS14 = cms.EDProducer("SlimmedElectronMvaIDProducer",
                                     electronsCollection = cms.InputTag("slimmedElectrons","","PAT"),
                                     method = cms.string("BDTSimpleCat"),
                                     mvaWeightFile = cms.vstring(
                                     "EgammaAnalysis/ElectronTools/data/PHYS14FIX/EIDmva_EB1_5_oldscenario2phys14FIX_BDT.weights.xml",
                                     "EgammaAnalysis/ElectronTools/data/PHYS14FIX/EIDmva_EB2_5_oldscenario2phys14FIX_BDT.weights.xml",
                                     "EgammaAnalysis/ElectronTools/data/PHYS14FIX/EIDmva_EE_5_oldscenario2phys14FIX_BDT.weights.xml",
                                     "EgammaAnalysis/ElectronTools/data/PHYS14FIX/EIDmva_EB1_10_oldscenario2phys14FIX_BDT.weights.xml",
                                     "EgammaAnalysis/ElectronTools/data/PHYS14FIX/EIDmva_EB2_10_oldscenario2phys14FIX_BDT.weights.xml",
                                     "EgammaAnalysis/ElectronTools/data/PHYS14FIX/EIDmva_EE_10_oldscenario2phys14FIX_BDT.weights.xml"
                                     ),
                                     Trig = cms.bool(False),
                                     )
     
process.Ana = cms.EDAnalyzer('UFHZZ4LAna',
                              photonSrc    = cms.untracked.InputTag("slimmedPhotons"),
                              electronSrc  = cms.untracked.InputTag("mvaNonTrigV025nsPHYS14","NonTrig"),
                              muonSrc      = cms.untracked.InputTag("boostedMuons"),
                              correctedJetSrc = cms.untracked.InputTag("slimmedJets"),
                              jetSrc       = cms.untracked.InputTag("slimmedJets"),
                              metSrc       = cms.untracked.InputTag("slimmedMETs"),
                              vertexSrc    = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"), #or selectedVertices 
                              isMC         = cms.untracked.bool(True),
                              isSignal     = cms.untracked.bool(True),
                              mH           = cms.untracked.double(125.0),
                              CrossSection = cms.untracked.double(0.00544046 ),
                              FilterEff    = cms.untracked.double(1),
                              weightEvents = cms.untracked.bool(True),
                              elRhoSrc     = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                              muRhoSrc     = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                              reweightForPU = cms.untracked.bool(True),
                              triggerSrc = cms.untracked.InputTag("TriggerResults","","HLT"),
                              triggerList = cms.untracked.vstring(
                                            'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1',
                                            'HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v1',
                                            'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1',
                                            'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1',
                                            'HLT_TripleMu_12_10_5_v1',
                                            'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v1',
                                            'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1',
                                            'HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v1',
                                            'HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v1',
                                            'HLT_Ele32_eta2p1_WP75_Gsf_v1' 
                              ),
                              verbose = cms.untracked.bool(False)                          
#                              verbose = cms.untracked.bool(True)
                             )

process.load('UFHZZAnalysisRun2.FSRPhotons.fsrPhotons_cff')

process.p = cms.Path(process.fsrPhotonSequence*
		     process.boostedMuons*
                     process.mvaNonTrigV025nsPHYS14*
                     process.Ana
                     )
