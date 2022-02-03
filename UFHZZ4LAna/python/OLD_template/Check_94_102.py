import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("UFHZZ4LAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.categories.append('UFHZZ4LAna')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.Services_cff')

process.GlobalTag.globaltag='94X_mc2017_realistic_v17'

process.Timing = cms.Service("Timing",
                             summaryOnly = cms.untracked.bool(True)
                             )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

myfilelist = cms.untracked.vstring(
      '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext2-v1/60000/FC213BBA-4138-E911-9D6E-549F3525C380.root'
      #'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext2-v1/60000/F2F57D81-F73F-E911-8E34-0CC47A13CD56.root',
    #'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext2-v1/60000/E8DFFED2-B038-E911-A8E8-842B2B68820B.root',
#'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext2-v1/60000/E47BDA90-0B39-E911-932B-008CFA197B24.root',
#'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext2-v1/60000/D0965A69-363B-E911-9159-842B2B689031.root'
#'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext2-v1/60000/D064D4E8-4555-E911-B2FF-A0369FD0B150.root'
#'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext2-v1/60000/B6C96279-3D3A-E911-9399-0CC47AA53D8A.root',


      # '/store/mc/RunIIFall17MiniAODv2/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/60000/3450B123-E8BF-E811-B895-FA163E9604CF.root',
        )

process.source = cms.Source("PoolSource",fileNames = myfilelist,
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("DUMMYFILENAME_newSC_FullRunII_v3.root")
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
                                         isSync = cms.bool(False),
                                         useRochester = cms.untracked.bool(True),
                                         year = cms.untracked.int32(2017)
                                         )

process.selectedElectrons = cms.EDFilter("PATElectronSelector",
                                         src = cms.InputTag("slimmedElectrons"),
                                         cut = cms.string("pt > 5 && abs(eta)<2.5")
                                         )

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedPatElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(123456), # for crab
        engineName = cms.untracked.string('TRandom3')
    )
)

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=False,
                       runVID=True,
                       phoIDModules=['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff'],
                       era='2017-Nov17ReReco')

process.load("RecoEgamma.EgammaTools.calibratedEgammas_cff")
process.calibratedPatElectrons.correctionFile = "EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2017_17Nov2017_v1_ele_unc"
process.calibratedPatElectrons.src = cms.InputTag("slimmedElectrons")

# FSR Photons
process.load('UFHZZAnalysisRun2.FSRPhotons.fsrPhotons_cff')

import os
# Jet Energy Corrections
from CondCore.DBCommon.CondDBSetup_cfi import *
era = "Fall17_17Nov2017_V32_94X_MC"
# for HPC
dBFile = os.environ.get('CMSSW_BASE')+"/src/UFHZZAnalysisRun2/UFHZZ4LAna/data/"+era+".db"
# for crab
#dBFile = "src/UFHZZAnalysisRun2/UFHZZ4LAna/data/"+era+".db"
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

process.jetCorrFactors = process.updatedPatJetCorrFactors.clone(
    src = cms.InputTag("slimmedJets"),
    levels = ['L1FastJet', 
              'L2Relative', 
              'L3Absolute'],
    payload = 'AK4PFchs' ) 

process.AK8PFJetCorrFactors = process.updatedPatJetCorrFactors.clone(
    src = cms.InputTag("slimmedJetsAK8"),
    levels = ['L1FastJet',
              'L2Relative',
              'L3Absolute'],
    payload = 'AK8PFchs' )

process.slimmedJetsJEC = process.updatedPatJets.clone(
    jetSource = cms.InputTag("slimmedJets"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("jetCorrFactors"))
    )

process.slimmedJetsAK8JEC = process.updatedPatJets.clone(
    jetSource = cms.InputTag("slimmedJetsAK8"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("AK8PFJetCorrFactors"))
    )

### add pileup id and discriminant to patJetsReapplyJEC
process.load("RecoJets.JetProducers.PileupJetID_cfi")
process.pileupJetIdUpdated = process.pileupJetId.clone(
    jets=cms.InputTag("slimmedJets"),
    inputIsCorrected=False,
    applyJec=True,
    vertexes=cms.InputTag("offlineSlimmedPrimaryVertices")
)
process.slimmedJetsJEC.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']
process.slimmedJetsJEC.userData.userInts.src += ['pileupJetIdUpdated:fullId']

# JER
process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
# for hpc
dBJERFile = os.environ.get('CMSSW_BASE')+"/src/UFHZZAnalysisRun2/UFHZZ4LAna/data/Fall17_V3_94X_MC.db"
# for crab
#dBJERFile = "src/UFHZZAnalysisRun2/UFHZZ4LAna/data/Fall17_V3_94X_MC.db"
process.jer = cms.ESSource("PoolDBESSource",
        CondDBSetup,
        connect = cms.string("sqlite_file:"+dBJERFile),
        toGet = cms.VPSet(
            cms.PSet(
                record = cms.string('JetResolutionRcd'),
                tag    = cms.string('JR_Fall17_V3_94X_MC_PtResolution_AK4PFchs'),
                label  = cms.untracked.string('AK4PFchs_pt')
                ),
            cms.PSet(
                record = cms.string('JetResolutionRcd'),
                tag    = cms.string('JR_Fall17_V3_94X_MC_PhiResolution_AK4PFchs'),
                label  = cms.untracked.string('AK4PFchs_phi')
                ),
            cms.PSet(
                record = cms.string('JetResolutionScaleFactorRcd'),
                tag    = cms.string('JR_Fall17_V3_94X_MC_SF_AK4PFchs'),
                label  = cms.untracked.string('AK4PFchs')
                )
            )
        )
process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')


#QGTag
process.load("CondCore.CondDB.CondDB_cfi")
qgDatabaseVersion = 'cmssw8020_v2'
# for hpc
QGdBFile = os.environ.get('CMSSW_BASE')+"/src/UFHZZAnalysisRun2/UFHZZ4LAna/data/QGL_"+qgDatabaseVersion+".db"
# for crab
#QGdBFile = "src/UFHZZAnalysisRun2/UFHZZ4LAna/data/QGL_"+qgDatabaseVersion+".db"
process.QGPoolDBESSource = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(messageLevel = cms.untracked.int32(1)),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
        cms.PSet(
            record = cms.string('QGLikelihoodRcd'),
            tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_AK4PFchs'),
            label  = cms.untracked.string('QGL_AK4PFchs')
        ),
      ),
      connect = cms.string('sqlite_file:'+QGdBFile)
)
process.es_prefer_qg = cms.ESPrefer('PoolDBESSource','QGPoolDBESSource')
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets = cms.InputTag( 'slimmedJetsJEC' )
process.QGTagger.jetsLabel = cms.string('QGL_AK4PFchs')
process.QGTagger.srcVertexCollection=cms.InputTag("offlinePrimaryVertices")

# compute corrected pruned jet mass
process.corrJets = cms.EDProducer ( "CorrJetsProducer",
                                    jets    = cms.InputTag( "slimmedJetsAK8JEC" ),
                                    vertex  = cms.InputTag( "offlineSlimmedPrimaryVertices" ), 
                                    rho     = cms.InputTag( "fixedGridRhoFastjetAll"   ),
                                    payload = cms.string  ( "AK8PFchs" ),
                                    isData  = cms.bool    (  False ),
                                    year    = cms.untracked.int32(2017))


# Recompute MET
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

runMetCorAndUncFromMiniAOD(process,
            isData=False,
            )

# STXS
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
    inputPruned = cms.InputTag("prunedGenParticles"),
    inputPacked = cms.InputTag("packedGenParticles"),
)
process.myGenerator = cms.EDProducer("GenParticles2HepMCConverter",
    genParticles = cms.InputTag("mergedGenParticles"),
    genEventInfo = cms.InputTag("generator"),
    signalParticlePdgIds = cms.vint32(25)
)
process.rivetProducerHTXS = cms.EDProducer('HTXSRivetProducer',
  HepMCCollection = cms.InputTag('myGenerator','unsmeared'),
  LHERunInfo = cms.InputTag('externalLHEProducer'),
  ProductionMode = cms.string('AUTO'),
)
# HZZ Fiducial from RIVET
process.rivetProducerHZZFid = cms.EDProducer('HZZRivetProducer',
  HepMCCollection = cms.InputTag('myGenerator','unsmeared'),
)



# Analyzer
process.Ana = cms.EDAnalyzer('UFHZZ4LAna',
                              photonSrc    = cms.untracked.InputTag("slimmedPhotons"),
                              electronUnSSrc  = cms.untracked.InputTag("slimmedElectrons"),
                              #electronSrc  = cms.untracked.InputTag("slimmedElectrons"),
                              electronSrc  = cms.untracked.InputTag("calibratedPatElectrons"),
                              muonSrc      = cms.untracked.InputTag("calibratedMuons"),
                              tauSrc      = cms.untracked.InputTag("slimmedTaus"),
                              jetSrc       = cms.untracked.InputTag("slimmedJetsJEC"),
                              mergedjetSrc = cms.untracked.InputTag("corrJets"),
                              metSrc       = cms.untracked.InputTag("slimmedMETs","","UFHZZ4LAnalysis"),
                              vertexSrc    = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),
                              beamSpotSrc  = cms.untracked.InputTag("offlineBeamSpot"),
                              conversionSrc  = cms.untracked.InputTag("reducedEgamma","reducedConversions"),
                              isMC         = cms.untracked.bool(True),
                              isSignal     = cms.untracked.bool(False),
                              mH           = cms.untracked.double(125.0),
                              CrossSection = cms.untracked.double(1),#DUMMYCROSSSECTION),
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
                              triggerObjects = cms.InputTag("slimmedPatTrigger"),
                              doJER = cms.untracked.bool(True),
                              doJEC = cms.untracked.bool(True),
                              doTriggerMatching = cms.untracked.bool(False),
                              triggerList = cms.untracked.vstring(
                                  # Single Lepton:
                                  'HLT_Ele35_WPTight_Gsf_v',
                                  'HLT_Ele38_WPTight_Gsf_v',
                                  'HLT_Ele40_WPTight_Gsf_v',
                                  'HLT_IsoMu27_v',
                                  # Dilepton
                                  'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v',
                                  'HLT_DoubleEle33_CaloIdL_MW_v',
                                  'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v',
                                  'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v',
                                  'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
                                  'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
                                  'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
                                  'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
                                  'HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v',
                                  'HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v',
                                  'HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v',
                                  # TriLepton
                                  'HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v',
                                  'HLT_TripleMu_10_5_5_DZ_v',
                                  'HLT_TripleMu_12_10_5_v',
                              ),
                              verbose = cms.untracked.bool(False),              
                              skimLooseLeptons = cms.untracked.int32(2),              
                              skimTightLeptons = cms.untracked.int32(2),              
                              #bestCandMela = cms.untracked.bool(False),
#                              verbose = cms.untracked.bool(True),
                              year = cms.untracked.int32(2017)              
                             )

process.p = cms.Path(process.fsrPhotonSequence*
                     process.boostedMuons*
                     process.calibratedMuons*
                     process.egmGsfElectronIDSequence*
                     process.egmPhotonIDSequence*
                     process.egammaPostRecoSeq*
                     process.calibratedPatElectrons*
                     process.jetCorrFactors*
                     process.pileupJetIdUpdated*
                     process.slimmedJetsJEC*
                     process.QGTagger*
                     process.AK8PFJetCorrFactors*
                     process.slimmedJetsAK8JEC*
                     process.fullPatMetSequence*
                     process.corrJets*
                     process.mergedGenParticles*process.myGenerator*process.rivetProducerHTXS*#process.rivetProducerHZZFid*
                     process.Ana
                     )
