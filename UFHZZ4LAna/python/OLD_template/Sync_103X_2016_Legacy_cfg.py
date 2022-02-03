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
#process.GlobalTag.globaltag='102X_upgrade2018_realistic_v15'
process.GlobalTag.globaltag='94X_mcRun2_asymptotic_v3'

process.Timing = cms.Service("Timing",
                             summaryOnly = cms.untracked.bool(True)
                             )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

myfilelist = cms.untracked.vstring(
'/store/mc/RunIISummer16MiniAODv2/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV709_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/20000/1A5C54BE-BED3-E711-B0A4-44A84224053C.root'
)

process.source = cms.Source("PoolSource",fileNames = myfilelist,
#                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
#                            eventsToProcess = cms.untracked.VEventRange('1:679835-1:679835','1:902899-1:902899','1:905104-1:905104','1:972090-1:972090','1:908104-1:908104','1:949166-1:949166')
                            )

process.TFileService = cms.Service("TFileService",
#                                   fileName = cms.string("Sync_102X_test.root")
#                                   fileName = cms.string("Sync_102X.root")
#                                   fileName = cms.string("Sync_1031_2016_check_6event.root")
                                   fileName = cms.string("Sync_1031_2016_v2.root")
)

# clean muons by segments 
process.boostedMuons = cms.EDProducer("PATMuonCleanerBySegments",
#				     #src = cms.InputTag("calibratedMuons"),#### was "slimmedMuons"
                                     src = cms.InputTag("slimmedMuons"),
				     preselection = cms.string("track.isNonnull"),
				     passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
				     fractionOfSharedSegments = cms.double(0.499),
				     )


# Kalman Muon Calibrations
process.calibratedMuons = cms.EDProducer("KalmanMuonCalibrationsProducer",
                                         muonsCollection = cms.InputTag("boostedMuons"),
                                         isMC = cms.bool(True),
                                         isSync = cms.bool(True),
                                         useRochester = cms.untracked.bool(True),
                                         year = cms.untracked.int32(2016)
                                         )

# Ghost cleaning
#process.cleanedMu = cms.EDProducer("PATMuonCleanerBySegments",
#                                    src = cms.InputTag("calibratedMuons"),
#                                    preselection = cms.string("track.isNonnull"),
#                                    passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
#                                    fractionOfSharedSegments = cms.double(0.499),
#                                    )


#from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
#process = regressionWeights(process)
#process.load('EgammaAnalysis.ElectronTools.regressionApplication_cff')

process.selectedElectrons = cms.EDFilter("PATElectronSelector",
                                         src = cms.InputTag("slimmedElectrons"),
#                                         cut = cms.string("pt > 5 && abs(eta)<2.5 && abs(-log(tan(superClusterPosition.theta/2)))<2.5")
                                         cut = cms.string("pt > 5 && abs(eta)<2.5")
                                         )

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedPatElectrons = cms.PSet(
        #initialSeed = cms.untracked.uint32(SEED), # for HPC
        initialSeed = cms.untracked.uint32(123456), # for crab
        engineName = cms.untracked.string('TRandom3')
    )
)
#process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')
#process.calibratedPatElectrons = cms.EDProducer("CalibratedPatElectronProducerRun2",
#                                        # input collections
#                                        electrons = cms.InputTag('selectedElectrons'),
#
#                                        gbrForestName = cms.vstring('electron_eb_ECALTRK_lowpt', 'electron_eb_ECALTRK',
#                                                                    'electron_ee_ECALTRK_lowpt', 'electron_ee_ECALTRK',
#                                                                    'electron_eb_ECALTRK_lowpt_var', 'electron_eb_ECALTRK_var',
#                                                                    'electron_ee_ECALTRK_lowpt_var', 'electron_ee_ECALTRK_var'),
#
#                                        isMC = cms.bool(True),
#                                        autoDataType = cms.bool(True),
#                                        isSynchronization = cms.bool(True),
#                                        correctionFile = cms.string("EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2017_17Nov2017_v1_ele_unc"),
#
#                                        recHitCollectionEB = cms.InputTag('reducedEgamma:reducedEBRecHits'),
#                                        recHitCollectionEE = cms.InputTag('reducedEgamma:reducedEERecHits')
#                                        )


#### new added
#from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
#setupEgammaPostRecoSeq(process,
#                          runEnergyCorrections=True,
#                          runVID=True,
#                          eleIDModules=['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Summer16_ID_ISO_cff','RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff'],
#                          era='2016-Legacy')
## and don't forget to add the producer 'process.egmGsfElectronIDSequence' to the path, i.e. process.electrons

###############################################
#####   mva calcution before calibrated   #####
###############################################

process.load("RecoEgamma.EgammaTools.calibratedEgammas_cff")
process.calibratedPatElectrons.correctionFile = "EgammaAnalysis/ElectronTools/data/ScalesSmearings/Legacy2016_07Aug2017_FineEtaR9_v3_ele_unc"
#process.calibratedPatElectrons.src = cms.InputTag("selectedElectrons")
process.calibratedPatElectrons.src = cms.InputTag("electronsMVA")


from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format to be DataFormat.MiniAOD, as appropriate
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
# define which IDs we want to produce
#my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff']
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Summer16_ID_ISO_cff','RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff']
# add them to the VID producer
for idmod in my_id_modules:
   setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

#process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag("calibratedPatElectrons")
#process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag("slimmedElectrons")

#process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('calibratedPatElectrons')
#process.electronMVAVariableHelper.srcMiniAOD = cms.InputTag('calibratedPatElectrons')
#process.electronMVAValueMapProducer.srcMiniAOD= cms.InputTag('calibratedPatElectrons')
process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('selectedElectrons')
process.electronMVAVariableHelper.srcMiniAOD = cms.InputTag('selectedElectrons')
process.electronMVAValueMapProducer.srcMiniAOD= cms.InputTag('selectedElectrons')

from RecoEgamma.EgammaTools.egammaObjectModificationsInMiniAOD_cff import egamma_modifications,egamma8XLegacyEtScaleSysModifier,egamma8XObjectUpdateModifier
egamma_modifications.append(egamma8XObjectUpdateModifier)

process.electronsMVA = cms.EDProducer("SlimmedElectronMvaIDProducer",
                                      #mvaValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV2RawValues"),
                                      mvaValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Summer16IdIsoValues"),
                                      #electronsCollection = cms.InputTag("calibratedPatElectrons"),
                                      #electronsCollection = cms.InputTag("slimmedElectrons"),
                                      electronsCollection = cms.InputTag("selectedElectrons"),
                                      idname = cms.string("ElectronMVAEstimatorRun2Summer16IdIsoValues"),
)

# FSR Photons
process.load('UFHZZAnalysisRun2.FSRPhotons.fsrPhotons_cff')

import os
# Jet Energy Corrections
from CondCore.DBCommon.CondDBSetup_cfi import *
#era = "Autumn18_V3_MC"
era = "Summer16_07Aug2017_V11_MC"
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
dBJERFile = os.environ.get('CMSSW_BASE')+"/src/UFHZZAnalysisRun2/UFHZZ4LAna/data/Summer16_25nsV1_MC.db"
## for crab
#dBFile = "src/UFHZZAnalysisRun2/UFHZZ4LAna/data/Autumn18_V1_MC.db"
process.jer = cms.ESSource("PoolDBESSource",
        CondDBSetup,
        connect = cms.string("sqlite_file:"+dBJERFile),
        toGet = cms.VPSet(
            cms.PSet(
                record = cms.string('JetResolutionRcd'),
                tag    = cms.string('JR_Summer16_25nsV1_MC_PtResolution_AK4PFchs'),
                label  = cms.untracked.string('AK4PFchs_pt')
                ),
            cms.PSet(
                record = cms.string('JetResolutionRcd'),
                tag    = cms.string('JR_Summer16_25nsV1_MC_PhiResolution_AK4PFchs'),
                label  = cms.untracked.string('AK4PFchs_phi')
                ),
            cms.PSet(
                record = cms.string('JetResolutionScaleFactorRcd'),
                tag    = cms.string('JR_Summer16_25nsV1_MC_SF_AK4PFchs'),
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
#process.QGTagger.srcJets = cms.InputTag( 'slimmedJets' )
process.QGTagger.jetsLabel = cms.string('QGL_AK4PFchs')
process.QGTagger.srcVertexCollection=cms.InputTag("offlinePrimaryVertices")

# compute corrected pruned jet mass
process.corrJets = cms.EDProducer ( "CorrJetsProducer",
                                    jets    = cms.InputTag( "slimmedJetsAK8JEC" ),
                                    vertex  = cms.InputTag( "offlineSlimmedPrimaryVertices" ), 
                                    rho     = cms.InputTag( "fixedGridRhoFastjetAll"   ),
                                    payload = cms.string  ( "AK8PFchs" ),
                                    isData  = cms.bool    (  False ),
                                    year    = cms.untracked.int32(2016))


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
                              electronUnSSrc  = cms.untracked.InputTag("electronsMVA"),
                              electronSrc  = cms.untracked.InputTag("calibratedPatElectrons"),
                              muonSrc      = cms.untracked.InputTag("calibratedMuons"),
#                              muonSrc      = cms.untracked.InputTag("boostedMuons"),
                              tauSrc      = cms.untracked.InputTag("slimmedTaus"),
                              jetSrc       = cms.untracked.InputTag("slimmedJetsJEC"),
#                              jetSrc       = cms.untracked.InputTag("slimmedJets"),
                              mergedjetSrc = cms.untracked.InputTag("corrJets"),
                              metSrc       = cms.untracked.InputTag("slimmedMETs","","UFHZZ4LAnalysis"),
#                              metSrc       = cms.untracked.InputTag("slimmedMETs"),
#                              vertexSrc    = cms.untracked.InputTag("goodPrimaryVertices"),####"offlineSlimmedPrimaryVertices"
                              vertexSrc    = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),
                              beamSpotSrc  = cms.untracked.InputTag("offlineBeamSpot"),
                              conversionSrc  = cms.untracked.InputTag("reducedEgamma","reducedConversions"),
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
                              doJER = cms.untracked.bool(True),
                              doJEC = cms.untracked.bool(True),
                              doTriggerMatching = cms.untracked.bool(False),
                              triggerList = cms.untracked.vstring(     
                                   'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
                                   'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
                                   'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v',
                                   'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v',
                                   'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v',
                                   'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v',
                                   'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v',
                                   'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v',
                                   'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v',
                                   'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
                                   'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
                                   'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
                                   'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
                                   'HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v',
                                   'HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v',
                                   'HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v',
                                   'HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v',
                                   'HLT_TripleMu_12_10_5_v',
                                   'HLT_Ele25_eta2p1_WPTight_Gsf_v',
                                   'HLT_Ele27_WPTight_Gsf_v',
                                   'HLT_Ele27_eta2p1_WPLoose_Gsf_v',
                                   'HLT_Ele32_eta2p1_WPTight_Gsf_v',
                                   'HLT_IsoMu20_v',
                                   'HLT_IsoTkMu20_v',
                                   'HLT_IsoMu22_v',
                                   'HLT_IsoTkMu22_v',
                                   'HLT_IsoMu24_v',
                                   'HLT_IsoTkMu24_v',
                              ),
                              verbose = cms.untracked.bool(False),              
                              skimLooseLeptons = cms.untracked.int32(4),              
                              skimTightLeptons = cms.untracked.int32(4),              
#                              verbose = cms.untracked.bool(True),              
                              year = cms.untracked.int32(2016),####for year put 2016,2017, or 2018 to select correct Muon training and electron MVA
                              #BTagCut = cms.untracked.double(0.4184)####2016: 0.6321; 2017: 0.4941; 2018: 0.4184
                             )

process.p = cms.Path(process.fsrPhotonSequence*
                     process.boostedMuons*
                     process.calibratedMuons*
#                     process.regressionApplication*
                     process.selectedElectrons*
#                     process.calibratedPatElectrons*
                     process.egmGsfElectronIDSequence*
                     process.electronMVAValueMapProducer*
                     process.electronsMVA*                     
                     #process.egmGsfElectronIDSequence*
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
