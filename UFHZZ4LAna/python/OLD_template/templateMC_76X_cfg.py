import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing


process = cms.Process("UFHZZ4LAnalysis")


## Options and Output Report
process.options   = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True)
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.categories.append('UFHZZ4LAna')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.GlobalTag.globaltag='76X_mcRun2_asymptotic_v12'

process.Timing = cms.Service("Timing",
                             summaryOnly = cms.untracked.bool(True)
                             )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

myfilelist = cms.untracked.vstring(DUMMYFILELIST)

process.source = cms.Source("PoolSource",fileNames = myfilelist,
                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("DUMMYFILENAME.root")
)

# clean muons by segments 
process.boostedMuons = cms.EDProducer("PATMuonCleanerBySegments",
				     src = cms.InputTag("slimmedMuons"),
				     preselection = cms.string("track.isNonnull"),
				     passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
				     fractionOfSharedSegments = cms.double(0.499),
				     )


# Kalman Muon Calibrations
process.calibratedMuons = cms.EDProducer("KalmanMuonCalibrationsProducer",
                                         muonsCollection = cms.InputTag("boostedMuons"),
                                         isMC = cms.bool(True),
                                         isSync = cms.bool(False)
                                         )

process.selectedElectrons = cms.EDFilter("PATElectronSelector",
                                         src = cms.InputTag("slimmedElectrons"),
                                         cut = cms.string("pt > 5 && abs(eta)<2.5")
                                         )

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedPatElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(SEED),
        engineName = cms.untracked.string('TRandom3')
    )
)

process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')
process.calibratedPatElectrons = cms.EDProducer("CalibratedPatElectronProducerRun2",
                                        # input collections
                                        electrons = cms.InputTag('selectedElectrons'),
                                        gbrForestName = cms.string("gedelectron_p4combination_25ns"),
                                        isMC = cms.bool(True),
                                        isSynchronization = cms.bool(False),
                                        correctionFile = cms.string("EgammaAnalysis/ElectronTools/data/76X_16DecRereco_2015")
                                        )

# FSR Photons
process.load('UFHZZAnalysisRun2.FSRPhotons.fsrPhotons_cff')

# Jet Energy Corrections
import os
from CondCore.DBCommon.CondDBSetup_cfi import *
era = "Fall15_25nsV1_MC"
dBFile = os.environ.get('CMSSW_BASE')+"/src/UFHZZAnalysisRun2/UFHZZ4LAna/data/"+era+".db"
process.jec = cms.ESSource("PoolDBESSource",
                           CondDBSetup,
                           connect = cms.string("sqlite_file:"+dBFile),
                           toGet =  cms.VPSet(
        cms.PSet(
            record = cms.string("JetCorrectionsRecord"),
            tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PF"),
            label= cms.untracked.string("AK4PF")
            ),
        cms.PSet(
            record = cms.string("JetCorrectionsRecord"),
            tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PFchs"),
            label= cms.untracked.string("AK4PFchs")
            ),

        cms.PSet(
            record = cms.string("JetCorrectionsRecord"),
            tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK8PFchs"),
            label= cms.untracked.string("AK8PFchs")
            ),
        )
)

process.es_prefer_jec = cms.ESPrefer("PoolDBESSource",'jec')

process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")

process.jetCorrFactors = process.patJetCorrFactorsUpdated.clone(
    src = cms.InputTag("slimmedJets"),
    levels = ['L1FastJet', 
              'L2Relative', 
              'L3Absolute'],
    payload = 'AK4PFchs' ) 

process.AK8PFJetCorrFactors = process.patJetCorrFactorsUpdated.clone(
    src = cms.InputTag("packedPatJetsAK8PFCHS"),#slimmedJetsAK8"),
    levels = ['L1FastJet',
              'L2Relative',
              'L3Absolute'],
    payload = 'AK8PFchs' )

process.slimmedJetsJEC = process.patJetsUpdated.clone(
    jetSource = cms.InputTag("slimmedJets"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("jetCorrFactors"))
    )
process.slimmedJetsAK8JEC = process.patJetsUpdated.clone(
    jetSource = cms.InputTag("packedPatJetsAK8PFCHS"),#slimmedJetsAK8"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("AK8PFJetCorrFactors"))
    )

# JER
process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
dBJERFile = os.environ.get('CMSSW_BASE')+"/src/UFHZZAnalysisRun2/UFHZZ4LAna/data/Summer15_25nsV6.db"
process.jer = cms.ESSource("PoolDBESSource",
        CondDBSetup,
        connect = cms.string("sqlite_file:"+dBJERFile),
        toGet = cms.VPSet(
            cms.PSet(
                record = cms.string('JetResolutionRcd'),
                tag    = cms.string('JR_Summer15_25nsV6_MC_PtResolution_AK4PFchs'),
                label  = cms.untracked.string('AK4PFchs_pt')
                ),
            cms.PSet(
                record = cms.string('JetResolutionRcd'),
                tag    = cms.string('JR_Summer15_25nsV6_MC_PhiResolution_AK4PFchs'),
                label  = cms.untracked.string('AK4PFchs_phi')
                )
            )
        )

process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')


#QGTag
qgDatabaseVersion = 'v2b'
QGPoolDBESSource = cms.ESSource("PoolDBESSource",
      CondDBSetup,
      toGet = cms.VPSet(),
      connect = cms.string('frontier://FrontierProd/CMS_COND_PAT_000'),
)

for type in ['AK4PFchs']:
  QGPoolDBESSource.toGet.extend(cms.VPSet(cms.PSet(
    record = cms.string('QGLikelihoodRcd'),
    tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_'+type),
    label  = cms.untracked.string('QGL_'+type)
  )))

process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets=cms.InputTag("slimmedJetsJEC")    
# Could be reco::PFJetCollection or pat::JetCollection (both AOD and miniAOD)
process.QGTagger.jetsLabel=cms.string('QGL_AK4PFchs')        
# Other options: see https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion
process.QGTagger.srcVertexCollection=cms.InputTag("offlinePrimaryVertices")
#Additional parameters:
#process.QGTagger.jec=cms.string("")
# Provide the jet correction service if your jets are uncorrected, otherwise keep empty
#process.QGTagger.systematicsLabel = cms.string('')     
# Produce systematic smearings (not yet available, keep empty)



# Analyzer
process.Ana = cms.EDAnalyzer('UFHZZ4LAna',
                              photonSrc    = cms.untracked.InputTag("slimmedPhotons"),
                              electronSrc  = cms.untracked.InputTag("calibratedPatElectrons"),
                              muonSrc      = cms.untracked.InputTag("calibratedMuons"),
                              jetSrc       = cms.untracked.InputTag("slimmedJetsJEC"),
                              mergedjetSrc = cms.untracked.InputTag("slimmedJetsAK8JEC"),
                              metSrc       = cms.untracked.InputTag("slimmedMETs"),
                              vertexSrc    = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),
                              beamSpotSrc  = cms.untracked.InputTag("offlineBeamSpot"),
                              conversionSrc  = cms.untracked.InputTag("reducedEgamma","reducedConversions"),
                              isMC         = cms.untracked.bool(True),
                              isSignal     = cms.untracked.bool(True),
                              mH           = cms.untracked.double(125.0),
                              CrossSection = cms.untracked.double(DUMMYCROSSSECTION),
                              FilterEff    = cms.untracked.double(1),
                              weightEvents = cms.untracked.bool(True),
                              elRhoSrc     = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                              muRhoSrc     = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                              rhoSrcSUS    = cms.untracked.InputTag("fixedGridRhoFastjetCentralNeutral"),
                              pileupSrc     = cms.untracked.InputTag("slimmedAddPileupInfo"),
                              pfCandsSrc   = cms.untracked.InputTag("packedPFCandidates"),
                              fsrPhotonsSrc = cms.untracked.InputTag("boostedFsrPhotons"),
                              prunedgenParticlesSrc = cms.untracked.InputTag("prunedGenParticles"),
                              packedgenParticlesSrc = cms.untracked.InputTag("packedGenParticles"),
                              genJetsSrc = cms.untracked.InputTag("slimmedGenJets"),
                              generatorSrc = cms.untracked.InputTag("generator"),
                              lheInfoSrc = cms.untracked.InputTag("externalLHEProducer"),
                              reweightForPU = cms.untracked.bool(True),
                              triggerSrc = cms.InputTag("TriggerResults","","HLT"),
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
                                            'HLT_Ele23_WPLoose_Gsf_v',
                              ),
                              verbose = cms.untracked.bool(False)              
#                              verbose = cms.untracked.bool(True)              
                             )


process.p = cms.Path(process.fsrPhotonSequence*
                     process.boostedMuons*
                     process.calibratedMuons*
                     process.selectedElectrons*
                     process.calibratedPatElectrons*
                     process.jetCorrFactors*
                     process.slimmedJetsJEC*
                     process.QGTagger*
                     process.AK8PFJetCorrFactors*
                     process.slimmedJetsAK8JEC*
                     process.Ana
                     )
