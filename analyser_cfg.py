 
import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("estarTuple")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
        )

##given at the bottom
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

#process.source = cms.Source("PoolSource",
#fileNames = cms.untracked.vstring(
    #"dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/jshilpi/"
    #"dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/161/312/D45FD145-7959-E011-B9DF-003048F110BE.root"
    #"dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/161/312/6AC2F8D8-D657-E011-B171-0030487C6062.root"
    #"/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/008/8468EAFF-C855-E011-BB51-003048F117B4.root"
    #"/store/data/Run2011A/Photon/AOD/DoublePhoton-May10ReReco-v1/0000/AC418CF0-6C7D-E011-9415-001A64789D40.root"
    #"/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/071/5A2512AD-C37F-E011-AB7D-001D09F2932B.root"
    #"/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/166/947/8E055754-0F98-E011-8230-0030487CD704.root"
    #"/store/mc/Summer11/ZZ_TuneZ2_7TeV_pythia6_tauola/AODSIM/PU_S4_START42_V11-v1/0000/BC914F76-BCAC-E011-B71C-0026189438BC.root"
 #   )
#)

process.load("Configuration.StandardSequences.Geometry_cff")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = cms.string('GR_R_52_V7::All')
#process.GlobalTag.globaltag = cms.string('GR_R_52_V9::All')
###process.GlobalTag.globaltag = cms.string('GR_P_V40::All')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#from Configuration.AlCa.autoCond import autoCond
#process.GlobalTag.globaltag = autoCond['com10']
#process.GlobalTag.globaltag = cms.string('GR_P_V40_AN1::All')
#process.GlobalTag.globaltag = cms.string('GR_P_V41_AN2::All')
process.GlobalTag.globaltag = cms.string('FT_53_V6C_AN3::All')
#from Configuration.PyReleaseValidation.autoCond import autoCond
#process.GlobalTag.globaltag = autoCond['startup']

#process.GlobalTag.globaltag = cms.string('START311_V1G1::All')
process.load('Configuration/StandardSequences/MagneticField_cff')

#process.load("UserCode.HEEPSkims.HEEPSkim2Ele_cfi")

# PAT Layer 0+1
process.load('PhysicsTools/PatAlgos/patSequences_cff')
#process.load('PhysicsTools.PatAlgos.patTemplate_cfg')
#from PhysicsTools.PatAlgos import patTemplate_cfg 



#########################################SWITCHES TO CONFIGURE FILE################################################
##################################################################################################################
##Data Type
flagData = 'ON'
flagMC = 'OFF'

if(flagData == 'ON'):
    #from PhysicsTools.PatAlgos.patTemplate_cfg import *
    from PhysicsTools.PatAlgos.tools.coreTools import *
    removeMCMatching(process, ['All'],outputModules = [])

    
    

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    #"dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/jshilpi/"
    #"dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/161/312/D45FD145-7959-E011-B9DF-003048F110BE.root"
    #"dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/161/312/6AC2F8D8-D657-E011-B171-0030487C6062.root"
    #"/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/008/8468EAFF-C855-E011-BB51-003048F117B4.root"
    #"/store/data/Run2011A/Photon/AOD/DoublePhoton-May10ReReco-v1/0000/AC418CF0-6C7D-E011-9415-001A64789D40.root"
    #"/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/071/5A2512AD-C37F-E011-AB7D-001D09F2932B.root"
    #   "/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/166/947/8E055754-0F98-E011-8230-0030487CD704.root"
    #    "file:pickevents.root"
    #"/store/mc/Summer12/ZprimeToJJ_M-500_TuneD6T_8TeV-pythia6/AODSIM/PU_S7_START50_V15-v1/0000/FC649E5D-5B6E-E111-8C2A-E41F1318174C.root"
    #"/store/mc/Summer12/QCD_Pt_20_30_EMEnriched_TuneZ2star_8TeV_pythia6/AODSIM/PU_S7_START50_V15-v1/0000/F6922FFC-C975-E111-ABAC-00304867918E.root"
    #"/store/data/Run2012A/Photon/AOD/08Jun2012-v2/0000/FEB856EA-E9B2-E111-A38B-001E67398C87.root"
    #"/store/data/Run2012B/SinglePhoton/AOD/PromptReco-v1/000/193/752/D09AC7A3-789B-E111-B093-003048D375AA.root"
    #"/store/data/Run2012A/Photon/AOD/PromptReco-v1/000/190/450/3CC2524E-ED80-E111-BCFB-BCAEC518FF52.root"
    #"/store/data/Run2012A/Photon/AOD/PromptReco-v1/000/193/112/BA253BD8-1A96-E111-BCF9-001D09F27003.root"
    #"/store/data/Run2012A/Photon/AOD/PromptReco-v1/000/191/003/06AB5D3E-1C86-E111-8CD1-003048F1C832.root"
    #"/store/data/Run2012A/Photon/AOD/PromptReco-v1/000/193/556/D05474CD-A399-E111-BE72-001D09F24EE3.root"
    #"file:pickevents.root"
    #"/store/data/Run2012C/SinglePhoton/AOD/PromptReco-v2/000/199/008/7C77A6CC-8FD0-E111-90B2-003048D2BED6.root"
    "/store/data/Run2012C/SinglePhoton/AOD/PromptReco-v2/000/198/941/6A7AF04D-70CF-E111-801A-BCAEC5329719.root"
    )
                            )

    
#check flags
flags = [flagData, flagMC]
if (flags.count("ON") > 1):
    print "You are running too many things at the same time."
    sys.exit(-1)

#SKIMS
flagNoSkim = 'OFF'
flagSkim2El = 'ON'
flagSkimMultiSC = 'OFF'

applyLaserCorr = 'ON'

#Stores CMSSW level file as well as globe file
flagStoreCMSSW = 'OFF'
################################################################################################################
###############################################################################################################
if (flagData == 'ON'):
    hltLabel = 'HLT'
elif (flagMC == 'ON'):
    hltLabel = 'HLT'



#############################################for good collsions################################                                                                         
### HLT Filter##########################################                                                                                                                
#process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskAlgoTrigConfig_cff')
#process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
#process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')

#process.load("HLTrigger.special.HLTTriggerTypeFilter_cfi")
# 0=random, 1=physics, 2=calibration, 3=technical                                                                                                                       
#process.hltTriggerTypeFilter.SelectedTriggerType = 1

# L1 Trigger Bit Selection (bit 40 and 41 for BSC trigger)                                                                                                              
#process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')                                                                                  
#process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
#process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)                                                                                                          
#process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41)')                                                                                      

##################################good collisions############################################                                                                           
#process.L1T1coll=process.hltLevel1GTSeed.clone()
#process.L1T1coll.L1TechTriggerSeeding = cms.bool(True)
#process.L1T1coll.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39) AND NOT ((42 AND NOT 43) OR (43 AND NOT 42))')
##process.L1T1coll.L1SeedsLogicalExpression = cms.string('(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39) AND NOT ((42 AND NOT 43) OR (43 AND NOT 42))')                       
#process.l1tcollpath = cms.Path(process.L1T1coll)

#####instead of above bits : see here :                                                                                                                                 
##https://twiki.cern.ch/twiki/bin/viewauth/CMS/Collisions2010Recipes#Physics_Declared_bit_selection                                                                     

#process.load('HLTrigger.special.hltPhysicsDeclared_cfi')
#process.hltPhysicsDeclared.L1GtReadoutRecordTag = 'gtDigis'


#scraping event filter                                                                                                                                                  
process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(False),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
                                  #thresh = cms.untracked.double(0.20)                                                                                                  
                                  )
#good vertex                                                                                                                                                            
process.primaryVertexFilter = cms.EDFilter("VertexSelector",
                                           src = cms.InputTag("offlinePrimaryVertices"),
                                           cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"), # tracksSize() > 3 for the older cut           
                                           filter = cms.bool(True)   # otherwise it won't filter the events, just produce an empty vertex collection.                   
                                           )



###################################################################################################################################333                                  
#from PhysicsTools.PatAlgos.tools.metTools import *
#addPfMET(process, 'PF')
#addTcMET(process, 'TC')

process.demo = cms.EDAnalyzer('Analyser',
                              #rhoLabel              = cms.InputTag("kt6PFJets", "rho"),
                              ####monophoton code
                              #rhoLabel         = cms.InputTag("kt6PFJets25", "rho"),
                              #sigmaLabel       = cms.InputTag("kt6PFJets25", "sigma"),
                              #rhoLabel44       = cms.InputTag("kt6PFJets44", "rho"),
                              #sigmaLabel44     = cms.InputTag("kt6PFJets44", "sigma"),

                              rhoLabel         = cms.InputTag("kt6PFJets", "rho"),
                              sigmaLabel       = cms.InputTag("kt6PFJets", "sigma"),
                              rhoLabel44       = cms.InputTag("kt6PFJets", "rho"),
                              sigmaLabel44     = cms.InputTag("kt6PFJets", "sigma"),
                              PileupSrc                = cms.untracked.InputTag("addPileupInfo"),
                              #electronTag              = cms.untracked.InputTag("cleanPatElectrons"),
                              #photonTag                = cms.untracked.InputTag("cleanPatPhotons"),
                              electronTag              = cms.untracked.InputTag("selectedPatElectrons"),
                              photonTag                = cms.untracked.InputTag("selectedPatPhotons"),
                              PFmetTag                 = cms.untracked.InputTag("patMETsPF"),
                              TCmetTag                 = cms.untracked.InputTag("patMETsTC"),
                              metTag                   = cms.untracked.InputTag("patMETs"),
                              jetTag                   = cms.untracked.InputTag("cleanPatJets"),
                              pfjetTag         = cms.untracked.InputTag("selectedPatJetsAK5PF"),
                              HLTriggerResults         = cms.untracked.InputTag("TriggerResults","",hltLabel),
                              processName              = cms.untracked.string(hltLabel),
                              triggerEventTag          = cms.untracked.InputTag("hltTriggerSummaryAOD","",hltLabel),
                              hltTrigToMatch           = cms.PSet(
    nhlt                   = cms.untracked.int32(5),
    hltMatch               = cms.untracked.vstring("HLT_Jet15U","HLT_Jet30U","HLT_Jet50U","HLT_Jet70U","HLT_Jet100U"),
    hltTrigModule          = cms.untracked.vstring("HLT1CaloJet","HLT1CaloJet","HLT1CaloJet","HLT1CaloJet","HLT1CaloJet")
    ),
                              rechitBTag               = cms.untracked.InputTag("reducedEcalRecHitsEB"),
                              rechitETag               = cms.untracked.InputTag("reducedEcalRecHitsEE"),
                              HybridSuperClusterColl   = cms.untracked.InputTag("correctedHybridSuperClusters"),
                              EndcapSuperClusterColl   = cms.untracked.InputTag("correctedMulti5x5SuperClustersWithPreshower"),
                              genJetLabel              = cms.untracked.InputTag("sisCone5GenJets"),
                              ########adopted for UCSD modules                                                                                                          
                              useSecondaryTrigger = cms.untracked.bool(False),
                              secondaryTriggerON       = cms.bool(False),
                              ElectronColl_std = cms.InputTag("gsfElectrons"),
                              PhotonCollStd = cms.InputTag("photons"),
                              MuonColl = cms.InputTag("muons"),
                              #JetColl_it5 = cms.InputTag("iterativeCone5CaloJets"),
                              JetColl_it5 = cms.InputTag("ak5CaloJets"),
                              Debug_Level = cms.int32(0),
                              #VertexColl_std = cms.InputTag("offlinePrimaryVerticesWithBS"),
                              VertexColl_std = cms.InputTag("offlinePrimaryVertices"),
                              TrackColl = cms.InputTag("generalTracks"),

                              #############################################PF isolation#########################################
                              IsoDepPhoton = cms.VInputTag(cms.InputTag('phPFIsoDepositChargedPFIso'),
                                                           cms.InputTag('phPFIsoDepositGammaPFIso'),
                                                           cms.InputTag('phPFIsoDepositNeutralPFIso')),
                              IsoValPhoton = cms.VInputTag(cms.InputTag('phoPFIso:chIsoForGsfEle'),
                                                           cms.InputTag('phoPFIso:phIsoForGsfEle'),
                                                           cms.InputTag('phoPFIso:nhIsoForGsfEle'),
                                                           ),
                              Photons = cms.InputTag('photons'),
                              ########correction on 11th April
                              PFCandLabel = cms.InputTag("particleFlow"),
                              #############################

                              #######################################END of PF isolation#######################################
                              
                              # CUTS
                              GeneratorCuts = cms.PSet(EtCut = cms.double(10.0)),
                              GenJetCuts = cms.PSet(EtCut = cms.double(-1.0)),

                              SimHitCuts = cms.PSet(EnergyCut = cms.double(0.0)),
                              SimTrackCuts = cms.PSet(EnergyCut = cms.double(0.0)),

                              EcalHitCuts = cms.PSet(BarrelEnergyCut = cms.double(0.08),
                                                         EndcapEnergyCut = cms.double(0.24),
                                                         PreEnergyCut = cms.double(-999.0),
                                                         EcalMaxDR = cms.double(0.5),
                                                         KeepOutsideCone = cms.bool(True)),
                              HcalHitsCuts = cms.PSet(HBHEEnergyCut = cms.double(0.35),
                                                          HFEnergyCut = cms.double(1.0),
                                                          HOEnergyCut = cms.double(0.7),
                                                          HcalMaxDR = cms.double(0.6),
                                                          KeepOutsideCone = cms.bool(True)),

                              BasicClusterCuts = cms.PSet(EnergyCut = cms.double(-1.0)),
                              SuperClusterCuts = cms.PSet(EnergyCut = cms.double(-1.0)),
                              CaloTowerCuts = cms.PSet(EtCut = cms.double(1.0)),
                              TrackCuts = cms.PSet(PtCut = cms.double(-1.0)),
                              TPCuts = cms.PSet(tpLvpCut = cms.double(99999.0),
                                                    tpEtaCut = cms.double(2.5),
                                                    tpTvpCut = cms.double(99999.0),
                                                    tpPtCut = cms.double(1.0),
                                                    tpPdgidCut = cms.vint32(11, 22)),
                              IsoCuts = cms.PSet(InnerCone = cms.double(0.015),
                                                     OuterCone = cms.double(0.3)),
                              ElectronCuts = cms.PSet(EtCut = cms.double(0.0)),
MuonCuts = cms.PSet(PtCut = cms.double(-1.0)),
                              PhotonCuts = cms.PSet(EtCut = cms.double(0.0)),
                              ConvertedPhotonCuts = cms.PSet(EtCut = cms.double(0.0)),
                              JetCuts = cms.PSet(EnergyCut = cms.double(0.0)),
                              #####end of UCSD's variables
                              ######very useful variables for ele. pho, vertex and track are automatically stored
                              #######if u want to store teh other ones, make below 1
                              usefulpateleTreeVar      = cms.untracked.int32(0),
                              usefulpatphoTreeVar      = cms.untracked.int32(0),
                              usefulvertexTreeVar      = cms.untracked.int32(0),
                              usefultrackTreeVar       = cms.untracked.int32(0),
                              
                              hOverEConeSizeSC         = cms.untracked.double(0.15),
                              runGsftracks             = cms.untracked.int32(0),
                              runGsfelectrons          = cms.untracked.int32(0),
                              runGeneraltracks         = cms.untracked.int32(1),
                              runPreshower             = cms.untracked.int32(0),
                              runPatphotons            = cms.untracked.int32(1),
                              runPhoRechitInfo            = cms.untracked.int32(1), ####29th paril, 2012
                              runPatelectrons          = cms.untracked.int32(1),
                              runHLT                   = cms.untracked.int32(0),
                              runUCSDHLT               = cms.untracked.int32(1),
                              runSC                    = cms.untracked.int32(1),
                              rungenParticle           = cms.untracked.int32(0),
                              rungenJets               = cms.untracked.int32(0),
                              runpatJets               = cms.untracked.int32(0),
                              runPFmet                 = cms.untracked.int32(0),
                              runpfJets                = cms.untracked.int32(0),
                              runTCmet                 = cms.untracked.int32(0),
                              runmet                   = cms.untracked.int32(0),
                              runVertex                = cms.untracked.int32(1),
                              runPileupinfo            = cms.untracked.int32(0),
                              runSCremoval             = cms.untracked.int32(0),
                              debugHLT                 = cms.untracked.int32(0),
                              debugPho                 = cms.untracked.int32(0),
                              debugPhoPFIso            = cms.untracked.int32(0),
                              debugRho                 = cms.untracked.int32(0),
                              debugEle                 = cms.untracked.int32(0),
                              debugEventinfo           = cms.untracked.int32(0),
                              wantLocalFile            = cms.untracked.int32(1),
                              wantRFIOFile             = cms.untracked.int32(0),
                              loutputFile              = cms.untracked.string("data.root"),       
                              #loutputFile              = cms.untracked.string("estar_heep200.root"),
                              rfoutputFile             = cms.untracked.string("rfio:/castor/cern.ch/user/s/shilpi/realdata/ztoee/370p4/7tev/mcfilter/ztoeemc.root")

                              )
                              

#process.load("UserCode.HEEPSkims.HEEPSkim2Ele_cfi")
#from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *

#process.demo.barrelCuts = cms.PSet(heepBarrelCuts)
#process.demo.endcapCuts = cms.PSet(heepEndcapCuts)

#process.p = cms.Path(process.hltTriggerTypeFilter + process.noscraping + process.L1T1coll + process.primaryVertexFilter)
#process.p = cms.Path(process.noscraping + process.L1T1coll + process.primaryVertexFilter)
#process.p = cms.Path(process.noscraping + process.hltPhysicsDeclared + process.primaryVertexFilter)                                                                    

########################2 electrons filter > 20 GeV#######################################
process.goodElectronsCounter2 = cms.EDFilter("CandViewCountFilter",
                                             src = cms.InputTag("gsfElectrons"),
                                             minNumber = cms.uint32(2)
                                             )
###############################################################################################

#######################for SC filter############################################################
process.superClusterMerger =  cms.EDProducer("EgammaSuperClusterMerger",
                                             src = cms.VInputTag(cms.InputTag('correctedHybridSuperClusters'),
                                                                 cms.InputTag('correctedMulti5x5SuperClustersWithPreshower'))
                                             )

process.multiSCfilter = cms.EDFilter("multiSCfilter",
                                     superClusterCollection = cms.InputTag("superClusterMerger"),
                                     #nSCrequire = cms.int32(2),
                                     nSCrequire = cms.int32(1),
                                     minEt = cms.double(20.),
                                     maxEta = cms.double(2.5)
                                     )
##########################################################################################################3

##################NO SKIM################################################################

process.dummySelector = cms.EDFilter("CandViewCountFilter",
                                     src = cms.InputTag("gsfElectrons"),
                                     minNumber = cms.uint32(0)
                                     )
##################################################################################
process.load("UserCode.HEEPSkims.HEEPSkim2Ele_cfi")

if (flagSkimMultiSC == 'ON'):
    process.eventFilter = cms.Sequence(process.superClusterMerger*process.multiSCfilter)
elif (flagSkim2El == 'ON'):
    process.eventFilter = cms.Sequence(process.goodElectronsCounter2*process.HEEPSkim2Ele)
elif (flagNoSkim == 'ON'):
    process.eventFilter = cms.Sequence(process.dummySelector)
    


    


from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

if(flagStoreCMSSW=='ON'):
    process.out = cms.OutputModule("PoolOutputModule",
                                   fileName = cms.untracked.string('patTuple.root'),
                                   # save only events passing the full path
                                   SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('pathFilter') ),
                                   # save PAT Layer 1 output; you need a '*' to
                                   # unpack the list of commands 'patEventContent'
                                   #outputCommands = cms.untracked.vstring('drop *','keep *_kt6PFJets*_*_*' ,*patEventContent)

                                   #outputCommands = cms.untracked.vstring('drop *','keep *_kt6PFJets*_*_*' ,'keep L1GlobalTrigger*_*_*_*',*patEventContent)
                                   outputCommands = cms.untracked.vstring('drop *','keep *_kt6PFJets*_*_*',*patEventContent)
                                   )



###switch off MC matching
#if(flagData == 'ON'):
#    from PhysicsTools.PatAlgos.tools.coreTools import *                             
    #removeMCMatching(process, 'All')                                                
#    removeMCMatching(process, ["All"])
###############


####from 2011
# Add PF jets
from PhysicsTools.PatAlgos.tools.jetTools import *
addJetCollection(process,cms.InputTag('ak5PFJets'),
                                  'AK5', 'PF',
                                  doJTA        = True,
                                  doBTagging   = True,
                                  jetCorrLabel = ('AK5PF', cms.vstring(['L1FastJet','L2Relative', 'L3Absolute','L2L3Residual'])),
                                  doType1MET    = True,
                                  doL1Cleaning  = True,
                                  doL1Counters  = False,
                                  genJetCollection=cms.InputTag("ak5GenJets"),
                                  doJetID       = False,
                                  jetIdLabel    = "ak5"
                                 )

process.selectedPatJetsAK5PF.cut = cms.string('pt > 10')

    
    
from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, 'PF')
addTcMET(process, 'TC')

######used in 2011
# FastJet Subtraction
#process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#process.load("RecoJets.JetProducers.kt4PFJets_cfi")
#process.kt6PFJets = process.kt4PFJets.clone()
#process.kt6PFJets.rParam = 0.6
#process.kt6PFJets.doRhoFastjet = True
#process.kt6PFJets.Rho_EtaMax = cms.double(2.5)
#process.kt6PFJets.doAreaFastjet = True

#process.load('RecoJets.JetProducers.ak5PFJets_cfi')
#process.ak5PFJets.doAreaFastjet = True
#process.ak5PFJets.Ghost_EtaMax = cms.double(5.0)
#process.ak5PFJets.Rho_EtaMax = cms.double(5.0)

###monophoton code
#---------Fast Rho calculation-------------------------
#Rho for eta= 2.5
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load("RecoJets.JetProducers.kt4PFJets_cfi")
process.kt6PFJets25 = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJets25.Rho_EtaMax = cms.double(2.5)
process.kt6PFJets25.Ghost_EtaMax = cms.double(2.5)
process.fastjetSequence25 = cms.Sequence( process.kt6PFJets25 )

#Rho for eta=4.4
process.kt6PFJets44 = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJets44.Rho_EtaMax = cms.double(4.4)
process.kt6PFJets44.Ghost_EtaMax = cms.double(5.0)
process.fastjetSequence44 = cms.Sequence( process.kt6PFJets44 )

##################Finally using this not the above ones######################
###############################After mail from Poter#########################
#---------Fast Rho calculation-------------------------
#process.kt6PFJets.doRhoFastjet = True
#process.kt6PFJets.Rho_EtaMax = cms.double(4.4)

#################################################

#Fast jet correction
process.load("RecoJets.JetProducers.ak5PFJets_cfi")
process.ak5PFJets.doAreaFastjet = True


process.load('EGamma.EGammaAnalysisTools.photonIsoProducer_cfi')
#process.phoPFIso.verbose = True
process.phoPFIso.verbose = False


if(flagData == 'ON'):
    process.demo.rungenParticle           = cms.untracked.int32(0)
    process.demo.runPileupinfo            = cms.untracked.int32(0)
    process.demo.loutputFile              = cms.untracked.string("data.root")
    
if(flagMC == 'ON'):
    process.demo.rungenParticle           = cms.untracked.int32(1)
    process.demo.runPileupinfo            = cms.untracked.int32(1)
    process.demo.loutputFile              = cms.untracked.string("bkg.root")

if(flagStoreCMSSW=='ON'):
    process.pathFilter = cms.Path(process.eventFilter)
    
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.p0 = cms.Path(process.PFTau)


###luminosity to proceas
#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('193112:25-193112:35')
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

###process.analyserPath = cms.Sequence(process.patDefaultSequence+process.kt6PFJets+process.demo)
#process.analyserPath = cms.Sequence(process.kt6PFJets+process.ak5PFJets+process.patDefaultSequence+process.demo)
#process.analyserPath = cms.Sequence(process.kt6PFJets+process.ak5PFJets+process.patDefaultSequence+process.demo)
#process.analyserPath = cms.Sequence(process.fastjetSequence25+process.fastjetSequence44+process.ak5PFJets+process.phoPFIso+process.patDefaultSequence+process.demo)

#process.analyserPath = cms.Sequence(process.ak5PFJets+process.kt6PFJets+process.phoPFIso+process.patDefaultSequence+process.demo)

if(applyLaserCorr=='ON'):
    process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')
    process.analyserPath = cms.Sequence(process.ecalLaserCorrFilter*process.ak5PFJets+process.phoPFIso+process.patDefaultSequence+process.demo)
elif(applyLaserCorr=='OFF'):
    process.analyserPath = cms.Sequence(process.ak5PFJets+process.phoPFIso+process.patDefaultSequence+process.demo)

#process.analyserPath = cms.Sequence(process.ak5PFJets+process.phoPFIso+process.patDefaultSequence+process.demo)

#process.analyserPath = cms.Sequence(process.fastjetSequence25+process.fastjetSequence44+process.ak5PFJets+process.phoPFIso+process.demo)
process.p21 = cms.Path(process.eventFilter*process.analyserPath)

if(flagStoreCMSSW == 'ON'):
    process.endpath = cms.EndPath(process.out)
#print process.dumpPython()
