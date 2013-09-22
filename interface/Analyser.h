#include <iostream>
#include "TRFIOFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "Analyzer/Analyser/interface/myPatElectron.h"
#include "Analyzer/Analyser/interface/GlobeHLT.h"
#include "Analyzer/Analyser/interface/GlobeVertex.h"
#include "Analyzer/Analyser/interface/Limits.h"
//#include "Analyzer/Analyser/interface/HeepCutbits.h"
//#include "Analyzer/Analyser/interface/TightIDphocutbits.h"
//#include "SHarper/HEEPAnalyzer/interface/HEEPAnalyzerBarePAT.h"
//#include "SHarper/HEEPAnalyzer/interface/HEEPCutCodes.h"

///pile up reweight  
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

class Analyser : public edm::EDAnalyzer {
 public:
  explicit Analyser(const edm::ParameterSet&);
  ~Analyser();
  GlobeHLT* hlt;
  GlobeVertex* vertex_std;

 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup);

  //my functions
  double mye2overe9(const DetId id, const EcalRecHitCollection & recHits);

  // ----------member data ---------------------------
  // now using hltConfigProvider for prescales etc that can change run-by-run 
  // (particularly because this for sure happened to HLT_L1_ETT100)             
  HLTConfigProvider hltConfig_;

  edm::ESHandle<CaloTopology> theCaloTopo_;
  //edm::InputTag fastjetLabel_;
  edm::InputTag rhoLabel_;
  edm::InputTag sigmaLabel_;
  edm::InputTag rhoLabel44_;
  edm::InputTag sigmaLabel44_;
  edm::InputTag PileupSrc_;
  edm::InputTag rechitBLabel_;
  edm::InputTag rechitELabel_;
  edm::InputTag eleLabel_;
  edm::InputTag phoLabel_;
  edm::InputTag PFmetLabel_;
  edm::InputTag TCmetLabel_;
  edm::InputTag metLabel_;
  edm::InputTag jetLabel_;
  edm::InputTag pfjetLabel_;
  edm::InputTag hlTriggerResults_;  // Input tag for TriggerResults 
  edm::InputTag hybridSuperClusterColl_;
  edm::InputTag endcapSuperClusterColl_;
  edm::InputTag genJetLabel_;
  edm::InputTag vertexLabel_;
  

  /////PF photons - 20th July
  edm::InputTag recophoLabel_;
  std::vector<edm::InputTag> inputTagIsoDepPhotons_;
  std::vector<edm::InputTag> inputTagIsoValPhotonsPFId_;   

  /////PF photons - 11th April, 2013
  edm::InputTag PFCandLabel_;

  /////25th Dec, 2012 - PF isolation SC footprint removal
  double isolation_cone_size_forSCremoval_;
  edm::InputTag tag_pfCandidates_forSCremoval_;
  edm::InputTag tag_Vertices_forSCremoval_;
  double rechit_link_enlargement_forSCremoval_;
  

  std::string processName_;

  edm::InputTag triggerEventTag_;

  edm::LumiReWeighting LumiWeights_; //pileup weighting


  int nhlt_;
  edm::ParameterSet hltTrigToMatch_;
  std::vector<std::string> hltMatch_;
  std::vector<std::string> hltTrigModule_;
  std::string loutputFile_;
  std::string rfoutputFile_;

  //heep::EleSelector cuts_;

  int nevents;
  int ntracks5;

  bool init_;
  
  //const edm::TriggerNames &triggerNames_;  // TriggerNames class
  std::vector<std::string>  hlNames_;  // name of each HLT algorithm    

  std::map<std::string,int> HLT_chosen;
  std::map<std::string,int> L1_chosen;


  TFile *lf;
  TRFIOFile *rfiof;
  TTree *tree;

  /////25th Dec
  ///SC footprint removal
  edm::ParameterSet myiConfig;

  //rho correction factor
  double rho;
  double sigma;
  double rho44;
  double sigma44;
  
  //pile up
  int pvi;
  int pileup_bunchXing[100];
  int pileup_nvtx[100];
  double pileup_Avgnvtx;
  double pileupWeight;
  double pileupEventWeight;

  //fall11
  double pileuptrue;

  double hOverEConeSizeSC_;
  int runGsftracks_;
  int runGsfelectrons_;
  int runGeneraltracks_;
  int runVertex_;
  int runPreshower_;
  int runPatelectrons_;
  int runPatphotons_;
  int runPhoRechitInfo_; //29th paril, 2012
  int runPFmet_;
  int runTCmet_;
  int runmet_;
  int wantLocalFile_;
  int wantRFIOFile_;
  int runHLT_;
  int runUCSDHLT_;
  int runSC_;
  int rungenParticle_;
  int rungenJets_;
  int runpatJets_;
  int runpfJets_;
  int runPileupinfo_;

  int runSCremoval_;

  int debugHLT_;
  int debugPho_;
  int debugphoPFIso_;
  int debugEle_;
  int debugEventinfo_;
  int debugRho_;
  
  ////24th Jul, 2012
  int usefulpateleTreeVar_;
  int usefulpatphoTreeVar_;
  //int usefulvertexTreeVar_;
  int usefultrackTreeVar_;

  //tree variables
  
  //event properties
  int run;
  int event;
  int orbit;
  int bx;
  int lumis;
  int isData;

  int    gsfsize; 
  double gsfpt[100];
  int    gsfcharge[100];
  double gsfchi2[100];
  double gsfeta[100];
  int    gsfnumberOfLostHits[100];
  int    gsfnumberOfValidHits[100];
  int    gsflost[100];
  double gsfd0[100];
  double gsfdxy[100];
  double gsfdz[100];
  double gsfptin[100];
  double gsfptout[100];
  double gsffbrem[100];
  double gsfqoverp[100];
  double gsfvx[100];
  double gsfvy[100];
  double gsfvz[100];
  double gsfphi[100];
  double gsfndof[100];
  int    gsfrecHitsSize[100];
  double gsftheta[100];
  int    gsfqualityMask[100];
  double gsfouterX[100];
  double gsfouterY[100];
  double gsfouterZ[100];
  double gsfouterRadius[100];
  double gsfinnerX[100];
  double gsfinnerY[100];
  double gsfinnerZ[100];

  //elec variables
  int    elesize; 
  double elept[25];
  double elepx[25];
  double elepy[25];
  double elepz[25];
  double elephi[25];
  double eletheta[25];
  double elefbrem[25];
  int    eletrackerDrivenSeed[25];
  int    eleecalDrivenSeed[25];
  int    elecharge[25];

  double elesceta[25];
  double elescphi[25];
  double elescE[25];
  double elescpt[25];

  double elepfeta[25];
  double elepfphi[25];
  double elepfE[25];
  double elepfpt[25];


  double eledr03EcalRecHitSumEt[25];
  double eledr03HcalDepth1TowerSumEt[25];
  double eledr03HcalDepth2TowerSumEt[25];
  double eledr03HcalTowerSumEt[25];
  double eledr04EcalRecHitSumEt[25];
  double eledr04HcalDepth1TowerSumEt[25];
  double eledr04HcalDepth2TowerSumEt[25];
  double eledr04HcalTowerSumEt[25];
  double elee1x5[25];
  double elee2x5Max[25];
  double elee5x5[25];
  double eleeEleClusterOverPout[25];
  double eleeSeedClusterOverP[25];
  double eleeSeedClusterOverPout[25];
  double eleeSuperClusterOverP[25];
  double eleeta[25];
  double elehadronicOverEm[25];
  double elesigmaIetaIeta[25]; 
  double eledeltaEtaEleClusterTrackAtCalo[25];
  double eledeltaEtaSeedClusterTrackAtCalo[25];
  double eledeltaEtaSuperClusterTrackAtVtx[25];
  double eledeltaPhiEleClusterTrackAtCalo[25];
  double eledeltaPhiSeedClusterTrackAtCalo[25];
  double eledeltaPhiSuperClusterTrackAtVtx[25];
  double eleenergy[25];
  double elemva[25];
  int    elenumberOfTracks[25];
  
  double eledr03TkSumPt[25];
  double eledr04TkSumPt[25];

  //swiss cross
  double elemaxEnergyXtal[25];
  double eleswissCross[25];
  double eleswissBasedspikevar[25];

  //timing
  double eleseedtime[25];
  int  elerecoFlag[25];
  double eleseverityLevel[25];

  //gsf
  double eletrkpt[25];
  int    eletrkcharge[25];
  double eletrkchi2[25];
  double eletrketa[25];
  int    eletrknumberOfLostHits[25];
  int    eletrknumberOfValidHits[25];
  int    eletrklost[25];
  double eletrkd0[25];
  double eletrkdxy[25];
  double eletrkdz[25];
  double eletrkptin[25];
  double eletrkptout[25];
  double eletrkfbrem[25];
  double eletrkqoverp[25];
  double eletrkvx[25];
  double eletrkvy[25];
  double eletrkvz[25];
  double eletrkphi[25];
  double eletrkndof[25];
  int    eletrkrecHitsSize[25];
  double eletrktheta[25];
  int    eletrkqualityMask[25];
  double eletrkouterX[25];
  double eletrkouterY[25];
  double eletrkouterZ[25];
  double eletrkouterRadius[25];
  double eletrkinnerX[25];
  double eletrkinnerY[25];
  double eletrkinnerZ[25];

  //ctf
  double electfpt[25];
  int    electfcharge[25];
  double electfchi2[25];
  double electfeta[25];
  int    electfnumberOfLostHits[25];
  int    electfnumberOfValidHits[25];
  int    electflost[25];
  double electfd0[25];
  double electfdxy[25];
  double electfdz[25];
  double electfptin[25];
  double electfptout[25];
  double electffbrem[25];
  double electfqoverp[25];
  double electfvx[25];
  double electfvy[25];
  double electfvz[25];
  double electfphi[25];
  double electfndof[25];
  int    electfrecHitsSize[25];
  double electftheta[25];
  int    electfqualityMask[25];
  double electfouterX[25];
  double electfouterY[25];
  double electfouterZ[25];
  double electfouterRadius[25];
  double electfinnerX[25];
  double electfinnerY[25];
  double electfinnerZ[25];
  
  //general trk
  int gentrksize;
  TClonesArray *gentrkp4;
  TClonesArray *gentrk_vtxpos;

  int    gentrkcharge[MAX_TRACKS];
  double gentrkchi2[MAX_TRACKS];
  int    gentrknumberOfLostHits[MAX_TRACKS];
  int    gentrknumberOfValidHits[MAX_TRACKS];
  int    gentrklost[MAX_TRACKS];
  double gentrkd0[MAX_TRACKS];
  double gentrkdxy[MAX_TRACKS];
  double gentrkdz[MAX_TRACKS];
  double gentrkptin[MAX_TRACKS];
  double gentrkptout[MAX_TRACKS];
  double gentrkfbrem[MAX_TRACKS];
  double gentrkqoverp[MAX_TRACKS];
  double gentrkndof[MAX_TRACKS];
  int    gentrkrecHitsSize[MAX_TRACKS];
  int    gentrkqualityMask[MAX_TRACKS];
  double gentrkouterX[MAX_TRACKS];
  double gentrkouterY[MAX_TRACKS];
  double gentrkouterZ[MAX_TRACKS];
  double gentrkouterRadius[MAX_TRACKS];
  double gentrkinnerX[MAX_TRACKS];
  double gentrkinnerY[MAX_TRACKS];
  double gentrkinnerZ[MAX_TRACKS];
  Int_t gentrk_hpnvalid[MAX_TRACKS];
  Int_t gentrk_hpnlost[MAX_TRACKS];
  Int_t gentrk_hpnvalidpix[MAX_TRACKS];

  //ES info
  //PS clusters info                                                                                                                                              
  int ESclussizeX[50];
  double ESclusEtX[50][1000];
  double ESclusEX[50][1000];
  double ESclusetaX[50][1000];
  double ESclusphiX[50][1000];
  double ESclusnhitsX[50][1000];

  int ESclussizeY[50];
  double ESclusEtY[50][1000];
  double ESclusEY[50][1000];
  double ESclusetaY[50][1000];
  double ESclusphiY[50][1000];
  double ESclusnhitsY[50][1000];

  double ESclusEtot[50];
  
  //ES rechit info
  int EEscsize;
  double EEsceta[1000];
  double EEscet[1000];
  int ESrechitsize[50];
  double ESrechitE[50][10000];
  double ESrechitplane[50][10000];
  double ESrechitzside[50][10000];
  double ESrechitsix[50][10000];
  double ESrechitsiy[50][10000];

  //pat elec variables
  int    patelesize; 
  
  TClonesArray *patelsc;
  TClonesArray *patelp4;
  TClonesArray *patelmomvtx;
  TClonesArray *patelmomvtxconst;
  TClonesArray *patelmomcalo;
  TClonesArray *patelmomout;
  TClonesArray *patelposvtx;
  TClonesArray *patelposcalo;

  double patelept[MAX_ELECTRONS];
  double patelepx[MAX_ELECTRONS];
  double patelepy[MAX_ELECTRONS];
  double patelepz[MAX_ELECTRONS];
  double patelephi[MAX_ELECTRONS];
  double pateletheta[MAX_ELECTRONS];
  double patelefbrem[MAX_ELECTRONS];
  int    pateletrackerDrivenSeed[MAX_ELECTRONS];
  int    pateleecalDrivenSeed[MAX_ELECTRONS];
  int    patelecharge[MAX_ELECTRONS];

  double patelesceta[MAX_ELECTRONS];
  double patelescphi[MAX_ELECTRONS];
  double patelescE[MAX_ELECTRONS];
  double patelescpt[MAX_ELECTRONS];

  int patelescind[MAX_ELECTRONS];

  double patelepfeta[MAX_ELECTRONS];
  double patelepfphi[MAX_ELECTRONS];
  double patelepfE[MAX_ELECTRONS];
  double patelepfpt[MAX_ELECTRONS];


  double pateledr03EcalRecHitSumEt[MAX_ELECTRONS];
  double pateledr03HcalDepth1TowerSumEt[MAX_ELECTRONS];
  double pateledr03HcalDepth2TowerSumEt[MAX_ELECTRONS];
  double pateledr03HcalTowerSumEt[MAX_ELECTRONS];
  double pateledr04EcalRecHitSumEt[MAX_ELECTRONS];
  double pateledr04HcalDepth1TowerSumEt[MAX_ELECTRONS];
  double pateledr04HcalDepth2TowerSumEt[MAX_ELECTRONS];
  double pateledr04HcalTowerSumEt[MAX_ELECTRONS];
  double patelee1x5[MAX_ELECTRONS];
  double patelee2x5Max[MAX_ELECTRONS];
  double patelee5x5[MAX_ELECTRONS];
  double pateleeEleClusterOverPout[MAX_ELECTRONS];
  double pateleeSeedClusterOverP[MAX_ELECTRONS];
  double pateleeSeedClusterOverPout[MAX_ELECTRONS];
  double pateleeSuperClusterOverP[MAX_ELECTRONS];
  double pateleeta[MAX_ELECTRONS];
  double patelehadronicOverEm[MAX_ELECTRONS];
  double patelesigmaIetaIeta[MAX_ELECTRONS]; 
  double pateledeltaEtaEleClusterTrackAtCalo[MAX_ELECTRONS];
  double pateledeltaEtaSeedClusterTrackAtCalo[MAX_ELECTRONS];
  double pateledeltaEtaSuperClusterTrackAtVtx[MAX_ELECTRONS];
  double pateledeltaPhiEleClusterTrackAtCalo[MAX_ELECTRONS];
  double pateledeltaPhiSeedClusterTrackAtCalo[MAX_ELECTRONS];
  double pateledeltaPhiSuperClusterTrackAtVtx[MAX_ELECTRONS];
  double pateleenergy[MAX_ELECTRONS];
  double patelemva[MAX_ELECTRONS];
  int    patelenumberOfTracks[MAX_ELECTRONS];
  
  double pateledr03TkSumPt[MAX_ELECTRONS];
  double pateledr04TkSumPt[MAX_ELECTRONS];
  
  ///added on 20th sept
  double pateleecalEnergyError[MAX_ELECTRONS];
  double pateletrackMomentumError[MAX_ELECTRONS];
  //added on 26th Sept, 2011
  double pateleecalEnergy[MAX_ELECTRONS];
  double patelecaloEnergy[MAX_ELECTRONS];
  double patelecaloCorrectedEnergy[MAX_ELECTRONS];

  double patelecorrectedEcalEnergy[MAX_ELECTRONS];

  //gaps info
  int pateleisEB[MAX_ELECTRONS];
  int pateleisEBEEGap[MAX_ELECTRONS];
  int pateleisEBEtaGap[MAX_ELECTRONS];
  int pateleisEBGap[MAX_ELECTRONS];
  int pateleisEBPhiGap[MAX_ELECTRONS];
  int pateleisEE[MAX_ELECTRONS];
  int pateleisEEDeeGap[MAX_ELECTRONS];
  int pateleisEEGap[MAX_ELECTRONS];
  int pateleisEERingGap[MAX_ELECTRONS];

  //heep cutbit 
  int    pateleheepcutword[MAX_ELECTRONS];
  int    pateleheepid[MAX_ELECTRONS];
  int    samspateleheepid[MAX_ELECTRONS];

  //swiss cross
  double patelemaxEnergyXtal[MAX_ELECTRONS];
  double pateleswissCross[MAX_ELECTRONS];
  double pateleswissBasedspikevar[MAX_ELECTRONS];

  //timing
  double pateleseedtime[MAX_ELECTRONS];
  int patelerecoFlag[MAX_ELECTRONS];
  int patelecheckFlag[MAX_ELECTRONS];
  double patelekOutOfTime[MAX_ELECTRONS];
  double pateleseverityLevel[MAX_ELECTRONS];
  double patelee2e9[MAX_ELECTRONS];

  //gsf
  double pateletrkpt[MAX_ELECTRONS];
  int    pateletrkcharge[MAX_ELECTRONS];
  double pateletrkchi2[MAX_ELECTRONS];
  double pateletrketa[MAX_ELECTRONS];
  int    pateletrknumberOfLostHits[MAX_ELECTRONS];
  int    pateletrknumberOfValidHits[MAX_ELECTRONS];
  int    pateletrklost[MAX_ELECTRONS];
  double pateletrkd0[MAX_ELECTRONS];
  double pateletrkdxy[MAX_ELECTRONS];
  double pateletrkdz[MAX_ELECTRONS];
  double pateletrkptin[MAX_ELECTRONS];
  double pateletrkptout[MAX_ELECTRONS];
  double pateletrkfbrem[MAX_ELECTRONS];
  double pateletrkqoverp[MAX_ELECTRONS];
  double pateletrkvx[MAX_ELECTRONS];
  double pateletrkvy[MAX_ELECTRONS];
  double pateletrkvz[MAX_ELECTRONS];
  double pateletrkphi[MAX_ELECTRONS];
  double pateletrkndof[MAX_ELECTRONS];
  int    pateletrkrecHitsSize[MAX_ELECTRONS];
  double pateletrktheta[MAX_ELECTRONS];
  int    pateletrkqualityMask[MAX_ELECTRONS];
  double pateletrkouterX[MAX_ELECTRONS];
  double pateletrkouterY[MAX_ELECTRONS];
  double pateletrkouterZ[MAX_ELECTRONS];
  double pateletrkouterRadius[MAX_ELECTRONS];
  double pateletrkinnerX[MAX_ELECTRONS];
  double pateletrkinnerY[MAX_ELECTRONS];
  double pateletrkinnerZ[MAX_ELECTRONS];
  double patelepout[MAX_ELECTRONS];
  double patelepin[MAX_ELECTRONS];
  
  ////conversion
  int pateleExpectednumberOfHits[MAX_ELECTRONS];
  double pateleconvDist[MAX_ELECTRONS];
  double pateleconvDcot[MAX_ELECTRONS];
  double pateleconvRadius[MAX_ELECTRONS];

  ////1st Aug, 2013
  int patelenumberOfLostHits[MAX_ELECTRONS];
  //pat photons
  int    patphosize;
  TClonesArray *patphop4;
  TClonesArray *patphocalopos;
  TClonesArray *patphoscp4;

  double patphoecalRecHitSumEtConeDR03[MAX_PHOTONS];
  double patphohcalDepth1TowerSumEtConeDR03[MAX_PHOTONS];
  double patphohcalDepth2TowerSumEtConeDR03[MAX_PHOTONS];
  double patphohcalTowerSumEtConeDR03[MAX_PHOTONS];
  double patphotrkSumPtHollowConeDR03[MAX_PHOTONS];
  double patphotrkSumPtSolidConeDR03[MAX_PHOTONS] ;
  double patphonTrkHollowConeDR03[MAX_PHOTONS];
  double patphonTrkSolidConeDR03[MAX_PHOTONS];

  double patphoecalRecHitSumEtConeDR04[MAX_PHOTONS];
  double patphohcalDepth1TowerSumEtConeDR04[MAX_PHOTONS];
  double patphohcalDepth2TowerSumEtConeDR04[MAX_PHOTONS];
  double patphohcalTowerSumEtConeDR04[MAX_PHOTONS];
  double patphotrkSumPtHollowConeDR04[MAX_PHOTONS];
  double patphotrkSumPtSolidConeDR04[MAX_PHOTONS];
  double patphonTrkHollowConeDR04[MAX_PHOTONS];
  double patphonTrkSolidConeDR04[MAX_PHOTONS];

  double patphoe1x5[MAX_PHOTONS];
  double patphoe2x5[MAX_PHOTONS];
  double patphoe3x3[MAX_PHOTONS];
  double patphoe5x5[MAX_PHOTONS];
  double patphoeta[MAX_PHOTONS];
  double patphohadronicOverEm[MAX_PHOTONS];
  double patphosigmaIetaIeta[MAX_PHOTONS];
  double patphor1x5[MAX_PHOTONS];
  double patphor2x5[MAX_PHOTONS];
  double patphor9[MAX_PHOTONS];

  ///store the sigmaiphi
  double patphosigmaIetaIphi[MAX_PHOTONS];
  double patphosigmaIphiIphi[MAX_PHOTONS];

  double patphonumberOfTracks[MAX_PHOTONS];
  int    patphohasPixelSeed[MAX_PHOTONS];
  int    patphoisConvertedPhoton[MAX_PHOTONS];
  double patphomaxEnergyXtal[MAX_PHOTONS];
  
  //store index of the SC
  int patphoscind[MAX_PHOTONS];

  //gaps info
  int patphoisEB[MAX_PHOTONS];
  int patphoisEBEEGap[MAX_PHOTONS];
  int patphoisEBEtaGap[MAX_PHOTONS];
  int patphoisEBGap[MAX_PHOTONS];
  int patphoisEBPhiGap[MAX_PHOTONS];
  int patphoisEE[MAX_PHOTONS];
  int patphoisEEDeeGap[MAX_PHOTONS];
  int patphoisEEGap[MAX_PHOTONS];
  int patphoisEERingGap[MAX_PHOTONS];

  double patphomipChi2[MAX_PHOTONS];
  double patphomipIntercept[MAX_PHOTONS];
  int    patphomipIsHalo[MAX_PHOTONS];
  int patphomipNhitCone[MAX_PHOTONS];
  double patphomipSlope[MAX_PHOTONS];
  double patphomipTotEnergy[MAX_PHOTONS];

  //tight id
  int patphoTightIDcutword[MAX_PHOTONS];
  int patphotightid[MAX_PHOTONS];
  
  //swiss cross           
  double patphoswissCross[MAX_PHOTONS];
  double patphoswissBasedspikevar[MAX_PHOTONS];

  //timing
  double patphoseedtime[MAX_PHOTONS];
  int patphorecoFlag[MAX_PHOTONS];
  int patphocheckFlag[MAX_PHOTONS];
  double patphokOutOfTime[MAX_PHOTONS];
  double patphoseverityLevel[MAX_PHOTONS];
  double patphoe2e9[MAX_PHOTONS];
  
  //conversion
  int    patphoconvsize[MAX_PHOTONS];
  int    patphohasConversionTracks[MAX_PHOTONS];
  double patphoconvtxX[MAX_PHOTONS][10];
  double patphoconvtxY[MAX_PHOTONS][10];
  double patphoconvtxZ[MAX_PHOTONS][10];
  double patphoconvtxR[MAX_PHOTONS][10];

  

  //gsf
  double patphotrkpt[MAX_PHOTONS];
  int    patphotrkcharge[MAX_PHOTONS];
  double patphotrkchi2[MAX_PHOTONS];
  double patphotrketa[MAX_PHOTONS];
  int    patphotrknumberOfLostHits[MAX_PHOTONS];
  int    patphotrknumberOfValidHits[MAX_PHOTONS];
  int    patphotrklost[MAX_PHOTONS];
  double patphotrkd0[MAX_PHOTONS];
  double patphotrkdxy[MAX_PHOTONS];
  double patphotrkdz[MAX_PHOTONS];
  double patphotrkptin[MAX_PHOTONS];
  double patphotrkptout[MAX_PHOTONS];
  double patphotrkfbrem[MAX_PHOTONS];
  double patphotrkqoverp[MAX_PHOTONS];
  double patphotrkvx[MAX_PHOTONS];
  double patphotrkvy[MAX_PHOTONS];
  double patphotrkvz[MAX_PHOTONS];
  double patphotrkphi[MAX_PHOTONS];
  double patphotrkndof[MAX_PHOTONS];
  int    patphotrkrecHitsSize[MAX_PHOTONS];
  double patphotrktheta[MAX_PHOTONS];
  int    patphotrkqualityMask[MAX_PHOTONS];
  double patphotrkouterX[MAX_PHOTONS];
  double patphotrkouterY[MAX_PHOTONS];
  double patphotrkouterZ[MAX_PHOTONS];
  double patphotrkouterRadius[MAX_PHOTONS];
  double patphotrkinnerX[MAX_PHOTONS];
  double patphotrkinnerY[MAX_PHOTONS];
  double patphotrkinnerZ[MAX_PHOTONS];

  //gen pho info
  int patphomGenisJet[MAX_PHOTONS];
  int patphomGenisPhoton[MAX_PHOTONS];
  int patphomGenisElectron[MAX_PHOTONS];
  int patphomGenpdgId[MAX_PHOTONS];
  int patphomGenstatus[MAX_PHOTONS];
  int patphonummoth[MAX_PHOTONS];
  int patphomGenmompdgId[MAX_PHOTONS][10];
  int patphomGengranpdgId[MAX_PHOTONS];
  double patphomGentheta[MAX_PHOTONS];
  double patphomGeneta[MAX_PHOTONS];
  double patphomGenphi[MAX_PHOTONS];
  double patphomGenpt[MAX_PHOTONS];
  double patphomGenpx[MAX_PHOTONS];
  double patphomGenpy[MAX_PHOTONS];
  double patphomGenpz[MAX_PHOTONS];
  double patphomGenenergy[MAX_PHOTONS];


  ///29th paril, 2012 - rechit info
  int patphoncrys[MAX_PHOTONS];
  int patphocrysrecoFlag[MAX_PHOTONS][MAX_PHOCRYS];
  int patphocryscheckFlag[MAX_PHOTONS][MAX_PHOCRYS];
  int patphocrysrawId[MAX_PHOTONS][MAX_PHOCRYS];
  int patphocrysieta[MAX_PHOTONS][MAX_PHOCRYS];
  int patphocrysiphi[MAX_PHOTONS][MAX_PHOCRYS];
  double patphocrysenergy[MAX_PHOTONS][MAX_PHOCRYS];
  double patphocrystime[MAX_PHOTONS][MAX_PHOCRYS];
  double patphocrystimeErr[MAX_PHOTONS][MAX_PHOCRYS];

  ////PF info - 19th July, 2012
  double patphopfchargediso[MAX_PHOTONS];
  double patphopfphotoniso[MAX_PHOTONS];
  double patphopfneutraliso[MAX_PHOTONS];

  double patphonewpfchargediso[MAX_PHOTONS];
  double patphonewpfphotoniso[MAX_PHOTONS];
  double patphonewpfneutraliso[MAX_PHOTONS];
  
  ////worst vertex - 21st june
  double patphoPFphotonWorstChargedHadronIso[MAX_PHOTONS];

  //converson safe veto
  int patphopasselectronveto[MAX_PHOTONS];
  ///new hOe
  double patphohadTowOverEm[MAX_PHOTONS];
  
  //PFmet
  int pfmetsize;
  double pfmetpt[50];
  double pfmetphi[50];
  double pfmetsumEt[50];
  double pfmetpx[50];
  double pfmetpy[50];

   //TCmet
  int tcmetsize;
  double tcmetpt[50];
  double tcmetphi[50];
  double tcmetsumEt[50];
  double tcmetpx[50];
  double tcmetpy[50];

  //met
  int metsize;
  double metpt[50];
  double metphi[50];
  double metsumEt[50];
  double metpx[50];
  double metpy[50];

  //HLT
  int is_HLT_Ele10_SW_L1R_event;
  int is_HLT_Ele15_SW_L1R_event;
  int is_HLT_Ele15_SW_EleId_L1R_event;
  int is_HLT_Ele15_SW_LooseTrackIso_L1R_event;
  int is_HLT_Ele15_SC15_SW_LooseTrackIso_L1R_event;
  int is_HLT_Ele15_SC15_SW_EleId_L1R_event;
  int is_HLT_Ele20_SW_L1R_event;
  int is_HLT_Ele20_SC15_SW_L1R_event;
  int is_HLT_Ele25_SW_L1R_event;
  int is_HLT_Ele25_SW_EleId_LooseTrackIso_L1R_event;
  int is_HLT_DoubleEle5_SW_Jpsi_L1R_event;
  int is_HLT_DoubleEle5_SW_Upsilon_L1R_event;
  int is_HLT_DoubleEle10_SW_L1R_event;

  int is_HLT_L1SingleEG5_event;
  int is_HLT_L1SingleEG8_event;
  int is_HLT_Ele10_LW_L1R_event;
  int is_HLT_Ele10_LW_EleId_L1R_event;
  int is_HLT_Ele15_LW_L1R_event;
  int is_HLT_Ele15_SiStrip_L1R_event;
  int is_HLT_Ele15_SC10_LW_L1R_event;
  int is_HLT_Ele20_LW_L1R_event;
  int is_HLT_L1DoubleEG5_event;
  int is_HLT_DoubleEle5_SW_L1R_event;
  
  int is_HLT_Photon10_L1R_event;
  int is_HLT_Photon10_LooseEcalIso_TrackIso_L1R_event;
  int is_HLT_Photon15_L1R_event;
  int is_HLT_Photon20_LooseEcalIso_TrackIso_L1R_event;
  int is_HLT_Photon25_L1R_event;
  int is_HLT_Photon25_LooseEcalIso_TrackIso_L1R_event;
  int is_HLT_Photon30_L1R_1E31_event;
  int is_HLT_Photon70_L1R_event;
  int is_HLT_DoublePhoton10_L1R_event;
  int is_HLT_DoublePhoton15_L1R_event;
  int is_HLT_DoublePhoton15_VeryLooseEcalIso_L1R_event;
  
  int is_HLT_Photon20_Cleaned_L1R_event;
  int is_HLT_Photon30_Cleaned_L1R_event;
  int is_HLT_Photon50_Cleaned_L1R_event;
  int is_HLT_Photon70_Cleaned_L1R_event;
  int is_HLT_Photon20_L1R_event;
  int is_HLT_Photon25_Cleaned_L1R_event;
  int is_HLT_Photon70_Cleaned_L1R_v1_event;

  //jet hlt
  int is_HLT_Jet15U;
  int is_HLT_Jet30U;
  int is_HLT_Jet50U;
  int is_HLT_Jet70U;
  int is_HLT_Jet100U;

  int is_HLT_Jet70U_v2;
  int is_HLT_Jet100U_v2;


  //hlt prescale 
  int ntrigToMatch; 
  int hltprescale[100];

  //trig match
  int hltobjsize[100];
  double ohtrigpt[100][100];
  double ohtrigeta[100][100];
  double ohtrigphi[100][100];
  double ohltphodR[100][100];
  double ohltele1dR[100][100];
  double ohltele2dR[100][100];
  double ohltjetdR[100][100];
  
  double patphohltmindRmatch[100] [50];
  
  //SC info
  TClonesArray *scp4;
  TClonesArray *scxyz;
  
  double scetaWidth[MAX_SUPERCLUSTERS];
  double scphiWidth[MAX_SUPERCLUSTERS];
  double scrawEnergy[MAX_SUPERCLUSTERS];
  double scpreshowerEnergy[MAX_SUPERCLUSTERS];
  double scclustersSize[MAX_SUPERCLUSTERS];
  int scisEB[MAX_SUPERCLUSTERS];
  int scisEE[MAX_SUPERCLUSTERS];
  int scsize;
  int scEBsize;
  int scEEsize;
  double sceMax[MAX_SUPERCLUSTERS];
  double sce2x2[MAX_SUPERCLUSTERS];
  double sce3x3[MAX_SUPERCLUSTERS];
  double sce5x5[MAX_SUPERCLUSTERS];
  double scr4[MAX_SUPERCLUSTERS];
  double scr9[MAX_SUPERCLUSTERS];
  double scr25[MAX_SUPERCLUSTERS];
  double sce1bye4[MAX_SUPERCLUSTERS];
  double sce1bye9[MAX_SUPERCLUSTERS];
  double sce1bye25[MAX_SUPERCLUSTERS];
  double sce2e9[MAX_SUPERCLUSTERS];
  double scHoE1[MAX_SUPERCLUSTERS];
  double scHoE2[MAX_SUPERCLUSTERS];
  double scHoE[MAX_SUPERCLUSTERS];

  //geninfo
  int gensize;
  int genpdgid[MAX_GENERATOR];
  int genstatus[MAX_GENERATOR];
  int gencharge[MAX_GENERATOR];
  int genmother[MAX_GENERATOR];
  int genndau[MAX_GENERATOR];
  int gennmoth[MAX_GENERATOR];
  TClonesArray *genp4;
  TClonesArray *genvtx;


  //gen jets
  int genjetsize;
  double genjetpt[50];
  double genjetpx[50];
  double genjetpy[50];
  double genjetpz[50];
  double genjetenergy[50];
  double genjetphi[50];
  double genjeteta[50];
  double genjetemEnergy[50];
  double genjethadEnergy[50];
  double genjetinvisibleEnergy[50];
  
  //pat jets
  int patjetsize;
  double patjetpt[200];
  double patjetpx[200];
  double patjetpy[200];
  double patjetpz[200];
  double patjetet[200];
  double patjetenergy[200];
  double patjetphi[200];
  double patjeteta[200];
  int    patjethasOverlapsmu[200];
  int    patjethasOverlapsele[200];
  int    patjethasOverlapspho[200];
  int    patjethasOverlapstau[200];
  int    patjethasOverlapstkIsoele[200];
  double patjetchargedEmEnergy[200];
  double patjetchargedEmEnergyFraction[200];
  double patjetchargedHadronEnergy[200];
  double patjetchargedHadronEnergyFraction[200];
  double patjetchargedMultiplicity[200];
  double patjetemEnergyFraction[200];
  double patjetemEnergyInEB[200];
  double patjetemEnergyInEE[200];
  double patjetemEnergyInHF[200];
  double patjetenergyFractionHadronic[200];
  double patjethadEnergyInHB[200];
  double patjethadEnergyInHE[200];
  double patjethadEnergyInHF[200];
  double patjethadEnergyInHO[200]; 


  //pf jets
  int pfjetsize;
  double pfjetpt[MAX_JETS];
  double pfjetpx[MAX_JETS];
  double pfjetpy[MAX_JETS];
  double pfjetpz[MAX_JETS];
  double pfjetet[MAX_JETS];
  double pfjetenergy[MAX_JETS];
  double pfjetphi[MAX_JETS];
  double pfjeteta[MAX_JETS];
  int    pfjethasOverlapsmu[MAX_JETS];
  int    pfjethasOverlapsele[MAX_JETS];
  int    pfjethasOverlapspho[MAX_JETS];
  int    pfjethasOverlapstau[MAX_JETS];
  int    pfjethasOverlapstkIsoele[MAX_JETS];
  double pfjetchargedEmEnergy[MAX_JETS];
  double pfjetchargedEmEnergyFraction[MAX_JETS];
  double pfjetchargedHadronEnergy[MAX_JETS];
  double pfjetchargedHadronEnergyFraction[MAX_JETS];
  double pfjetchargedMultiplicity[MAX_JETS];
  
  double pfjetemEnergyFraction[MAX_JETS];
  double pfjetemEnergyInEB[MAX_JETS];
  double pfjetemEnergyInEE[MAX_JETS];
  double pfjetemEnergyInHF[MAX_JETS];
  double pfjetenergyFractionHadronic[MAX_JETS];
  double pfjethadEnergyInHB[MAX_JETS];
  double pfjethadEnergyInHE[MAX_JETS];
  double pfjethadEnergyInHF[MAX_JETS];
  double pfjethadEnergyInHO[MAX_JETS]; 
  double pfjet_jecCorr[MAX_JETS];
  
  //added on 31st October, 2011 for fake rate
  double pfjetchargedMuEnergyFraction[MAX_JETS];
  //double pfjethadEnergy[MAX_JETS];
  

  
};
