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
process.GlobalTag.globaltag='102X_upgrade2018_realistic_v18'

process.Timing = cms.Service("Timing",
                             summaryOnly = cms.untracked.bool(True)
                             )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

myfilelist = cms.untracked.vstring(
        '/store/mc/RunIIAutumn18MiniAOD/GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/260000/F82D86C0-17CB-164F-A1A5-2BF309636E74.root'
        
        
        #'/store/mc/RunIIAutumn18MiniAOD/GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/130000/8A25152B-D130-7141-8090-06294F2560D9.root',
        #'/store/mc/RunIIAutumn18MiniAOD/GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/130000/82DE55FB-7A5A-6B46-B9C5-60CD92DDC3C9.root',
        #'/store/mc/RunIIAutumn18MiniAOD/GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/130000/7F87014D-6536-5A44-A258-325987436451.root',
        #'/store/mc/RunIIAutumn18MiniAOD/GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/130000/3F8C3110-01ED-7144-BBBE-46780BEC8D95.root',
        #'/store/mc/RunIIAutumn18MiniAOD/GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/130000/783D58CF-E335-E44E-AFEA-48E794867E95.root',
        #'/store/mc/RunIIAutumn18MiniAOD/GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/130000/4BCB8F25-7CCA-5D4C-979D-CD5F3546E8AB.root',
        #'/store/mc/RunIIAutumn18MiniAOD/GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/130000/35295C52-4A47-CF4B-AA02-5329753FCE3B.root'
 
 
 
 
 #       '/store/mc/RunIIAutumn18MiniAOD/ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/60000/FBC570F1-EF01-674F-8A44-081A9481EF84.root',
  #      '/store/mc/RunIIAutumn18MiniAOD/ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/60000/F36E6483-C42E-9645-A0FB-E88F5707DB2F.root',
   #     '/store/mc/RunIIAutumn18MiniAOD/ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/60000/135CA2F1-2F60-364B-8E32-495982753E59.root',
    #    '/store/mc/RunIIAutumn18MiniAOD/ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/40483F05-74DE-F143-A833-64C1AD8016E3.root',
     #   '/store/mc/RunIIAutumn18MiniAOD/ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/3E16C789-AD09-A84C-BCCB-77B7B7C2390B.root',
      #  '/store/mc/RunIIAutumn18MiniAOD/ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/60000/DCB7927B-269F-3B4B-9DA3-EFE07A37FC9E.root',
       
        #DUMMYFILELIST
        )

process.source = cms.Source("PoolSource",fileNames = myfilelist,
           #                 duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
           #eventsToProcess = cms.untracked.VEventRange('1:12:66017')
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
                                         isSync = cms.bool(False),
                                         useRochester = cms.untracked.bool(True),
                                         year = cms.untracked.int32(2018)
                                         )

#from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
#process = regressionWeights(process)
#process.load('EgammaAnalysis.ElectronTools.regressionApplication_cff')

process.selectedElectrons = cms.EDFilter("PATElectronSelector",
                                         src = cms.InputTag("slimmedElectrons"),
                                         #cut = cms.string("pt > 5 && abs(eta)<2.5 && abs(-log(tan(superClusterPosition.theta/2)))<2.5")
                                         cut = cms.string("pt > 5 && abs(eta)<2.5 ")
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
#                                        #electrons = cms.InputTag('selectedElectrons'),
#                                        electrons = cms.InputTag('electronsMVA'),
#
#                                        gbrForestName = cms.vstring('electron_eb_ECALTRK_lowpt', 'electron_eb_ECALTRK',
#                                                                    'electron_ee_ECALTRK_lowpt', 'electron_ee_ECALTRK',
#                                                                    'electron_eb_ECALTRK_lowpt_var', 'electron_eb_ECALTRK_var',
#                                                                    'electron_ee_ECALTRK_lowpt_var', 'electron_ee_ECALTRK_var'),
#
#                                        isMC = cms.bool(True),
#                                        autoDataType = cms.bool(True),
#                                        isSynchronization = cms.bool(False),
#                                        #correctionFile = cms.string("EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2017_17Nov2017_v1_ele_unc"),
#                                        correctionFile = cms.string("EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2018_Step2Closure_CoarseEtaR9Gain_v2"),
#
#                                        recHitCollectionEB = cms.InputTag('reducedEgamma:reducedEBRecHits'),
#                                        recHitCollectionEE = cms.InputTag('reducedEgamma:reducedEERecHits')
#
#
#                                        )

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=False,
                       runVID=True,
                       eleIDModules=['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Autumn18_ID_ISO_cff','RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff'],
                       phoIDModules=['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff'],
                       era='2018-Prompt')

###############################################
#####   mva calcution before calibrated   #####
###############################################

process.load("RecoEgamma.EgammaTools.calibratedEgammas_cff")
#process.calibratedPatElectrons.correctionFile = "EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2018_Step2Closure_CoarseEtaR9Gain"
process.calibratedPatElectrons.correctionFile = "EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2018_Step2Closure_CoarseEtaR9Gain_v2"
#process.calibratedPatElectrons.src = cms.InputTag("selectedElectrons")
#process.calibratedPatElectrons.src = cms.InputTag("electronsMVA")
process.calibratedPatElectrons.src = cms.InputTag("slimmedElectrons")

##  from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
##  dataFormat = DataFormat.MiniAOD
##  switchOnVIDElectronIdProducer(process, dataFormat)
##  # define which IDs we want to produce
##  my_id_modules = [ 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Autumn18_ID_ISO_cff','RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff' ]
##  # add them to the VID producer
##  for idmod in my_id_modules:
##      setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
##  #process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag("calibratedPatElectrons")
##  #process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons')
##  #process.electronMVAVariableHelper.srcMiniAOD = cms.InputTag('slimmedElectrons')
##  #process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag("slimmedElectrons")
##  process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('selectedElectrons')
##  process.electronMVAVariableHelper.srcMiniAOD = cms.InputTag('selectedElectrons')
##  process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag('selectedElectrons')
##  
##  process.electronsMVA = cms.EDProducer("SlimmedElectronMvaIDProducer",
##                                        mvaValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Autumn18IdIsoValues"),
##  #                                      electronsCollection = cms.InputTag("calibratedPatElectrons"),
##                                        #electronsCollection = cms.InputTag("slimmedElectrons"),
##                                        electronsCollection = cms.InputTag("selectedElectrons"),
##                                        idname = cms.string("ElectronMVAEstimatorRun2Autumn18IdIsoValues"),
##  )

# FSR Photons
process.load('UFHZZAnalysisRun2.FSRPhotons.fsrPhotons_cff')

import os
# Jet Energy Corrections
from CondCore.DBCommon.CondDBSetup_cfi import *
era = "Autumn18_V19_MC"
# for HPC
dBFile = os.environ.get('CMSSW_BASE')+"/src/UFHZZAnalysisRun2/UFHZZ4LAna/data/"+era+".db"
# for crab
dBFile = "src/UFHZZAnalysisRun2/UFHZZ4LAna/data/"+era+".db"
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
dBJERFile = os.environ.get('CMSSW_BASE')+"/src/UFHZZAnalysisRun2/UFHZZ4LAna/data/Autumn18_V7_MC.db"   
# for crab
dBJERFile = "src/UFHZZAnalysisRun2/UFHZZ4LAna/data/Autumn18_V7_MC.db"
process.jer = cms.ESSource("PoolDBESSource",
        CondDBSetup,
        connect = cms.string("sqlite_file:"+dBJERFile),
        toGet = cms.VPSet(
            cms.PSet(
                record = cms.string('JetResolutionRcd'),
                tag    = cms.string('JR_Autumn18_V7_MC_PtResolution_AK4PFchs'),
                label  = cms.untracked.string('AK4PFchs_pt')
                ),
            cms.PSet(
                record = cms.string('JetResolutionRcd'),
                tag    = cms.string('JR_Autumn18_V7_MC_PhiResolution_AK4PFchs'),
                label  = cms.untracked.string('AK4PFchs_phi')
                ),
            cms.PSet(
                record = cms.string('JetResolutionScaleFactorRcd'),
                tag    = cms.string('JR_Autumn18_V7_MC_SF_AK4PFchs'),
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
QGdBFile = "src/UFHZZAnalysisRun2/UFHZZ4LAna/data/QGL_"+qgDatabaseVersion+".db"
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
                                    year = cms.untracked.int32(2018))


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
                              #electronSrc  = cms.untracked.InputTag("electronsMVA"),
                              #electronUnSSrc  = cms.untracked.InputTag("electronsMVA"),
                              electronUnSSrc  = cms.untracked.InputTag("slimmedElectrons"),
                              electronSrc  = cms.untracked.InputTag("calibratedPatElectrons"),
                              muonSrc      = cms.untracked.InputTag("calibratedMuons"),
#                              muonSrc      = cms.untracked.InputTag("boostedMuons"),
                              tauSrc      = cms.untracked.InputTag("slimmedTaus"),
                              jetSrc       = cms.untracked.InputTag("slimmedJetsJEC"),
#                              jetSrc       = cms.untracked.InputTag("slimmedJets"),
                              mergedjetSrc = cms.untracked.InputTag("corrJets"),
                              metSrc       = cms.untracked.InputTag("slimmedMETs","","UFHZZ4LAnalysis"),
                              #metSrc       = cms.untracked.InputTag("slimmedMETs"),
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
                              triggerObjects = cms.InputTag("selectedPatTrigger"),
                              doJER = cms.untracked.bool(True),
                              doJEC = cms.untracked.bool(True),
                              doTriggerMatching = cms.untracked.bool(False),
                              triggerList = cms.untracked.vstring(
                                  # Toni
                                  'HLT_Ele32_WPTight_Gsf_v', 
                                  'HLT_IsoMu24_v',
                                  'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v',
                                  'HLT_DoubleEle25_CaloIdL_MW_v',
                                  'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v',
                                  'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
                                  'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
                                  'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
                                  'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
                                  'HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v',
                                  'HLT_TripleMu_10_5_5_DZ_v',             
                                  'HLT_TripleMu_12_10_5_v',               
                                  'HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v',   
                                  'HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v',
                                  # OLD
#                                  'HLT_Ele32_WPTight_Gsf_v',
#                                  'HLT_IsoMu24_v',
#                                  'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v',
#                                  'HLT_DoubleEle25_CaloIdL_MW_v',
#                                  'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v',
#                                  'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v',
#                                  'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
#                                  'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
#                                  'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
#                                  'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
#                                  'HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v',
##                                  'HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v',
#                                  'HLT_TripleMu_10_5_5_DZ_v',
#                                  'HLT_TripleMu_12_10_5_v',
                              ),
                              verbose = cms.untracked.bool(False),              
                              skimLooseLeptons = cms.untracked.int32(4),              
                              skimTightLeptons = cms.untracked.int32(4),              
                              #bestCandMela = cms.untracked.bool(False),
#                              verbose = cms.untracked.bool(True),              
                              year = cms.untracked.int32(2018),####for year put 2016,2017, or 2018 to select correct setting
                             )


process.p = cms.Path(process.fsrPhotonSequence*
                     process.boostedMuons*
                     process.calibratedMuons*
#                     process.regressionApplication*
        #             process.selectedElectrons*
#                     process.calibratedPatElectrons*
                     process.egmGsfElectronIDSequence*
        #             process.electronMVAValueMapProducer*
        #             process.electronsMVA*
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

