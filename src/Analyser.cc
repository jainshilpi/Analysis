// -*- C++ -*-
//
// Package:    Analyser
// Class:      Analyser
// 
/**\class Analyser Analyser.cc Analyzer/Analyser/src/Analyser.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Shilpi Jain,40 4-B28,+41227671589,
//         Created:  Fri May 21 19:03:17 CEST 2010
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"

//user files
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h" 
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h" 
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

//ES cluster collection
#include "DataFormats/EgammaReco/interface/PreshowerClusterFwd.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
//ES rechit
#include "DataFormats/EcalDetId/interface/ESDetId.h"

//geometry
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
//HLT
#include <DataFormats/Common/interface/TriggerResults.h>
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"

//sc
#include <DataFormats/EgammaReco/interface/SuperCluster.h>
//gen info
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
//genjets
#include <DataFormats/JetReco/interface/GenJet.h>

//used in timing and spike info
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatusCode.h"
//#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

///PF jets
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

////Jet corrections
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"



//for H/E calculation of a SC
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"
#include "TRFIOFile.h" 

#include "Analyzer/Analyser/interface/Analyser.h"
//#include "/uscms_data/d2/jshilpi/CMSSW_3_11_1/src/Analyzer/Analyser/interface/TightIDphocutbits.h"
//#include "/uscms_data/d2/jshilpi/CMSSW_3_11_1/src/Analyzer/Analyser/interface/HeepCutbits.h"
//#include "/uscms_data/d2/jshilpi/CMSSW_3_11_1/src/Analyzer/Analyser/src/ID.cc"


////SHarpers HEEP id bit
//#include "SHarper/HEEPAnalyzer/interface/HEEPAnalyzerBarePAT.h"
//#include "SHarper/HEEPAnalyzer/interface/HEEPCutCodes.h"

//cmssw utility functions
#include "DataFormats/Math/interface/deltaR.h"

//pile up
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

///pile up reweight
 #include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

/////single crystal info
 #include "DataFormats/DetId/interface/DetId.h"

///PF iso info
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"

//conversion safe veto
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

//////removal of SC footprint from PF
#include "PFIsolation/SuperClusterFootprintRemoval/interface/SuperClusterFootprintRemoval.h"
//utility function prototypes


/////April 11, 2013
// PF isolation from alternate code
#include "EGamma/EGammaAnalysisTools/src/PFIsolationEstimator.cc"


double deltaphi(double phi1, double phi2);
double correct_phi(double phi);
double delta_R(double phi,double eta);
double Theta(double eta);
double Pl(double P,double Pt);

//ID function prototypes
int heepeleID( std::vector<pat::Electron> &patelectron_container,int ipatel );
int tightphoID( std::vector<pat::Photon> &patphoton_container,int ipatpho );

//myEcalSeverityLevelAlgo.cc function prototypes ---- for calculating e2e9
float recHitE( const DetId id, const EcalRecHitCollection &recHits );
float recHitE( const DetId id, const EcalRecHitCollection & recHits, int dEta, int dPhi );
float recHitApproxEt( const DetId id, const EcalRecHitCollection &recHits );


///for PF iso
#include "Analyzer/Analyser/interface/myrecophoclass.h"
/////////sort criterium
class PtSortCriterium{
public:
  bool operator() (myrecophoclass p1,myrecophoclass p2){
    return p1.pt >= p2.pt;
  }
};


//
// class declaration
//

/////moved to interface

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Analyser::Analyser(const edm::ParameterSet& iConfig):
  rhoLabel_(iConfig.getParameter<edm::InputTag>("rhoLabel")),
  sigmaLabel_(iConfig.getParameter<edm::InputTag>("sigmaLabel")),
  rhoLabel44_(iConfig.getParameter<edm::InputTag>("rhoLabel44")),
  sigmaLabel44_(iConfig.getParameter<edm::InputTag>("sigmaLabel44")),
  PileupSrc_(iConfig.getUntrackedParameter<edm::InputTag>("PileupSrc")),
  eleLabel_(iConfig.getUntrackedParameter<edm::InputTag>("electronTag")),
  phoLabel_(iConfig.getUntrackedParameter<edm::InputTag>("photonTag")),
  vertexLabel_(iConfig.getParameter<edm::InputTag>("VertexColl_std")),
  ////PF photons - 20th July
  recophoLabel_(iConfig.getParameter<edm::InputTag>("Photons")),
  inputTagIsoDepPhotons_(iConfig.getParameter< std::vector<edm::InputTag> >("IsoDepPhoton")),
  inputTagIsoValPhotonsPFId_(iConfig.getParameter< std::vector<edm::InputTag> >("IsoValPhoton")),
  /////////new PF photons - 11th April
  PFCandLabel_(iConfig.getParameter<edm::InputTag>("PFCandLabel")),
  PFmetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("PFmetTag")),
  TCmetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("TCmetTag")),
  metLabel_(iConfig.getUntrackedParameter<edm::InputTag>("metTag")),
  jetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("jetTag")),
  pfjetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("pfjetTag")),
  hlTriggerResults_(iConfig.getUntrackedParameter<edm::InputTag>("HLTriggerResults")),
  processName_(iConfig.getUntrackedParameter<std::string>("processName")),
  triggerEventTag_(iConfig.getUntrackedParameter<edm::InputTag>("triggerEventTag")),
  hltTrigToMatch_(iConfig.getParameter<edm::ParameterSet>("hltTrigToMatch")),  //which trigger you want to consider for matching
  //hltTrigModule_(iConfig.getUntrackedParameter<std::string>("hltTrigModule")),  //name of last filter module of that trigger
  rechitBLabel_(iConfig.getUntrackedParameter<edm::InputTag>("rechitBTag")),
  rechitELabel_(iConfig.getUntrackedParameter<edm::InputTag>("rechitETag")),
  hybridSuperClusterColl_(iConfig.getUntrackedParameter<edm::InputTag>("HybridSuperClusterColl")),
  endcapSuperClusterColl_(iConfig.getUntrackedParameter<edm::InputTag>("EndcapSuperClusterColl")),
  genJetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("genJetLabel")),
  hOverEConeSizeSC_(iConfig.getUntrackedParameter<double>("hOverEConeSizeSC",0.15)),
  //cuts_(iConfig),
  /////for reducing the tree branch size
  usefulpateleTreeVar_(iConfig.getUntrackedParameter<int>("usefulpateleTreeVar",0)),
  usefulpatphoTreeVar_(iConfig.getUntrackedParameter<int>("usefulpatphoTreeVar",0)),
  //usefulvertexTreeVar_(iConfig.getUntrackedParameter<int>("usefulvertexTreeVar",0)),
  usefultrackTreeVar_(iConfig.getUntrackedParameter<int>("usefultrackTreeVar",0)),
  ////////25th dec - for SC footprint removal
  //isolation_cone_size_forSCremoval_(iConfig.getUntrackedParameter<double>("isolation_cone_size_forSCremoval",0.3)),
  isolation_cone_size_forSCremoval_(iConfig.getUntrackedParameter<double>("isolation_cone_size_forSCremoval",0.3)),
  tag_pfCandidates_forSCremoval_(iConfig.getUntrackedParameter<edm::InputTag>("tag_pfCandidates_forSCremoval",edm::InputTag("particleFlow"))),
  tag_Vertices_forSCremoval_(iConfig.getUntrackedParameter<edm::InputTag>("tag_Vertices_forSCremoval",edm::InputTag("offlinePrimaryVertices"))),
  rechit_link_enlargement_forSCremoval_(iConfig.getUntrackedParameter<double>("rechit_link_enlargement_forSCremoval",0.25)),
  runGsftracks_(iConfig.getUntrackedParameter<int>("runGsftracks",1)),
  runGsfelectrons_(iConfig.getUntrackedParameter<int>("runGsfelectrons",1)),
  runGeneraltracks_(iConfig.getUntrackedParameter<int>("runGeneraltracks",1)),
  runVertex_(iConfig.getUntrackedParameter<int>("runVertex",1)),
  runPreshower_(iConfig.getUntrackedParameter<int>("runPreshower",1)),
  runPatelectrons_(iConfig.getUntrackedParameter<int>("runPatelectrons",1)),
  runPatphotons_(iConfig.getUntrackedParameter<int>("runPatphotons",1)),
  runPhoRechitInfo_(iConfig.getUntrackedParameter<int>("runPhoRechitInfo",0)),
  runHLT_(iConfig.getUntrackedParameter<int>("runHLT",1)),
  runUCSDHLT_(iConfig.getUntrackedParameter<int>("runUCSDHLT",1)),
  runSC_(iConfig.getUntrackedParameter<int>("runSC",1)),
  rungenParticle_(iConfig.getUntrackedParameter<int>("rungenParticle",1)),
  rungenJets_(iConfig.getUntrackedParameter<int>("rungenJets",1)),
  runpatJets_(iConfig.getUntrackedParameter<int>("runpatJets",1)),
  runPFmet_(iConfig.getUntrackedParameter<int>("runPFmet",1)),
  runTCmet_(iConfig.getUntrackedParameter<int>("runTCmet",1)),
  runmet_(iConfig.getUntrackedParameter<int>("runmet",1)),
  runpfJets_(iConfig.getUntrackedParameter<int>("runpfJets")),
  runPileupinfo_(iConfig.getUntrackedParameter<int>("runPileupinfo",1)),
  runSCremoval_(iConfig.getUntrackedParameter<int>("runSCremoval",0)),
  debugHLT_(iConfig.getUntrackedParameter<int>("debugHLT",0)),
  debugPho_(iConfig.getUntrackedParameter<int>("debugPho",0)),
  debugphoPFIso_(iConfig.getUntrackedParameter<int>("debugphoPFIso",0)),
  debugEle_(iConfig.getUntrackedParameter<int>("debugEle",0)),
  debugEventinfo_(iConfig.getUntrackedParameter<int>("debugEventinfo",0)),
  debugRho_(iConfig.getUntrackedParameter<int>("debugRho",0)),
  wantLocalFile_(iConfig.getUntrackedParameter<int>("wantLocalFile",1)),
  wantRFIOFile_(iConfig.getUntrackedParameter<int>("wantRFIOFile",0)),
  loutputFile_(iConfig.getUntrackedParameter<std::string>("loutputFile", "gsftrack.root")),
  rfoutputFile_(iConfig.getUntrackedParameter<std::string>("rfoutputFile", "/uscms_data/d2/sushil/CMSSW/MonoPhoton/CMSSW_3_8_6/src/QCDFakeRate/Analyser/test/gsftrack.root")),
  init_(false)
{
  std::cout<<"inside the constructor"<<std::endl;
  //now do what ever initialization is needed
  nhlt_           = hltTrigToMatch_.getUntrackedParameter<int>("nhlt",1); 
  hltMatch_       = hltTrigToMatch_.getUntrackedParameter<std::vector<std::string> >("hltMatch");
  hltTrigModule_  = hltTrigToMatch_.getUntrackedParameter<std::vector<std::string> >("hltTrigModule");
  nevents  = 0;
  ntracks5 = 0;

  //call UCSD's HLT method
  if (runUCSDHLT_)
    hlt = new GlobeHLT(iConfig);
    
  //CALL UCSD's vertex method
  vertex_std   = new GlobeVertex(iConfig, "std");

  ////for SC footprint removal - 25th December, 2012
  myiConfig = iConfig;
  
  std::cout<<"going outside the constructor now"<<std::endl;
}


Analyser::~Analyser()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  std::cout<<"inside destructor"<<std::endl;
  if(wantLocalFile_)
    {
      lf->Close();
      delete lf;
    }

  if(wantRFIOFile_)
    {
      rfiof->Close();
      delete rfiof;
    }


}


void Analyser::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  if (runUCSDHLT_)
    hlt->beginRun(iRun, iSetup);
  
}
 
//
// member functions
//

// ------------ method called to for each event  ------------
void
Analyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;
   
   using namespace std;

   //std::cout<<"inside analyze"<<std::endl;
   //vectors to pat objects     
   run    = iEvent.id().run();
   event  = iEvent.id().event();
   orbit  = iEvent.orbitNumber();
   bx     = iEvent.bunchCrossing();
   lumis  = iEvent.luminosityBlock();
   isData = iEvent.isRealData();
  
   if(debugEventinfo_){
     cout<<"run:lumi:event"<<run<<":"<<lumis<<":"<<event<<endl;
   }


   /*///RHO CORRECTION FACTOR
   edm::Handle<double> fastjethandle ;
   iEvent.getByLabel(fastjetLabel_, fastjethandle);

   rho = *fastjethandle;
   */

   ///monophoton code
   // Rho correction with 2.5
   edm::Handle<double> rhoHandle;
   iEvent.getByLabel(rhoLabel_, rhoHandle);
   rho=0.;
   if(rhoHandle.isValid()) {
     rho= *(rhoHandle.product());
   }


   edm::Handle<double> sigmaHandle;
   iEvent.getByLabel(sigmaLabel_, sigmaHandle);
   sigma =0.;
   if(sigmaHandle.isValid()) {      
     sigma = *(sigmaHandle.product());
   }                 



   // Rho correction with max eta 4.4
   edm::Handle<double> rhoHandle44;
   iEvent.getByLabel(rhoLabel44_, rhoHandle44);
   rho44=0.;
   if(rhoHandle44.isValid()) {
     rho44= *(rhoHandle44.product());
   }


   edm::Handle<double> sigmaHandle44;
   iEvent.getByLabel(sigmaLabel44_, sigmaHandle44);
   sigma44 =0.;
   if(sigmaHandle44.isValid()) {      
     sigma44 = *(sigmaHandle44.product());
   }                 

   if(debugRho_)
     {
       cout<<"rho:"<<rho<<endl;
     }


   //pile up info
   //edm::InputTag PileupSrc_ = "addPileupInfo";
   
   if(runPileupinfo_)
     {
       Handle<std::vector< PileupSummaryInfo > >  PupInfo;
       iEvent.getByLabel(PileupSrc_, PupInfo);
       
       std::vector<PileupSummaryInfo>::const_iterator PVI;
       
       float sum_nvtx = 0;

       for(int ii=0; ii<100; ii++)
	 {
	   pileup_nvtx[ii] = -1;
	 }
       
       pvi = 0;
       // (then, for example, you can do)
       for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
	 
	 //std::cout << " Pileup Information: bunchXing, nvtx: " << PVI->getBunchCrossing() << " " << PVI->getPU_NumInteractions() << std::endl;
	 
	 pileup_bunchXing[pvi] = PVI->getBunchCrossing();

	 //if(pileup_bunchXing == 0) { 
	   pileup_nvtx[pvi] = PVI->getPU_NumInteractions();
	   sum_nvtx += double(pileup_nvtx[pvi]);
	   //continue;
	   //}
	   if(PVI->getBunchCrossing() == 0) { 
	     pileuptrue = PVI->getTrueNumInteractions(); ///valid for fall11 samples - same for all bunchxings
	   }
	   
	   
	   pvi++;
       }//for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI)
       

       float ave_nvtx = sum_nvtx/3.;
       pileup_Avgnvtx = ave_nvtx;
       
       //pileupWeight = LumiWeights_.weight( pileup_nvtx );
       //pileupWeight = LumiWeights_.weight3BX( ave_nvtx );
       //edm::EventBase* iEventB = dynamic_cast<edm::EventBase*>(&iEvent);
       //pileupEventWeight = LumiWeights_.weight( (*iEventB) );

     }//if(runPileupinfo_)
   



   //initialize hltConfigProvider ..... used in HLT-RECO object match
   /*bool changed(true);
   //if(hltConfig_.init(iEvent,processName_,changed))
   if(hltConfig_.init(iEvent.getRun(),iSetup,processName_,changed))
     {
       if(changed)
	 {
	   std::cout << "HLT Config has changed at run "<< run<<std::endl;
	   edm::LogVerbatim("hlt-RECOmatch") << "Successfully initialized HLTConfigProvider with process name: "<<processName_;
	 }
     }//if(hltConfig_.init(iEvent.getRun(),iSetup,processName_,changed))
   
   else
     {
       std::cout << "ERROR: Couldnt get HLT configuration with name " << processName_ <<std::endl;
     }
   */


     
   //common collections, eg:right now used for recoelectrons but can be used for pat electrons
   edm::Handle<EcalRecHitCollection> Brechit;//barrel
   edm::Handle<EcalRecHitCollection> Erechit;//endcap
   iEvent.getByLabel(rechitBLabel_,Brechit);
   iEvent.getByLabel(rechitELabel_,Erechit);
   EcalClusterLazyTools lazyTool( iEvent, iSetup, rechitBLabel_,rechitELabel_  );

   const EcalRecHitCollection* barrelRecHits= Brechit.product();  ///used in reco::electrons for swiss cross 
   const EcalRecHitCollection* endcapRecHits= Erechit.product();


   edm::ESHandle<CaloTopology> pTopology;
   iSetup.get<CaloTopologyRecord>().get(theCaloTopo_);
   const CaloTopology *topology = theCaloTopo_.product();

   // get the channel status from the DB ----- used in accessing severity level
   edm::ESHandle<EcalChannelStatus> chStatus;
   iSetup.get<EcalChannelStatusRcd>().get(chStatus);

   
   if(runGsftracks_)
     {
       edm::Handle<reco::GsfTrackCollection> gsfTracksH; 
       bool found=iEvent.getByLabel("electronGsfTracks",gsfTracksH); 
       if(!found ) { 
	 std::ostringstream  err; 
	 err<<" cannot get GsfTracks:GsfTracks " 
	    <<std::endl; 
	 edm::LogError("TrackExtra")<<err.str(); 
	 throw cms::Exception( "MissingProduct", err.str()); 
       }   
       
       std::vector<reco::GsfTrack>  gsftrack_container;
       
       unsigned ncandidates=gsfTracksH->size(); 
       
       int ntracks=0;
       /*for (unsigned igsf=0;igsf<ncandidates;++igsf)
	 {
	 if (std::sqrt((*gsfTracksH)[igsf].innerMomentum().Perp2())>5.0)
	 ++ntracks;
	 
	 }
       */
       
       for(int igsf = 0; igsf<100; igsf++)
	 {
	   gsfpt[igsf]                         = -999.0;
	   gsfcharge[igsf]                     = -999;
	   gsfchi2[igsf]                       = -999.0;
	   gsfeta[igsf]                        = -999.0;
	   gsfnumberOfLostHits[igsf]           = -999;
	   gsfnumberOfValidHits[igsf]          = -999;
	   gsflost[igsf]                       = -999;
	   gsfd0[igsf]                         = -999.0;
	   gsfdxy[igsf]                        = -999.0;
	   gsfdz[igsf]                         = -999.0;
	   gsfptin[igsf]                       = -999.0;
	   gsfptout[igsf]                      = -999.0;
	   gsffbrem[igsf]                      = -999.0;
	   gsfqoverp[igsf]                     = -999.0;
	   gsfvx[igsf]                         = -999.0;
	   gsfvy[igsf]                         = -999.0;
	   gsfvz[igsf]                         = -999.0;
	   gsfphi[igsf]                        = -999.0;
	   gsfndof[igsf]                       = -999.0;
	   gsfrecHitsSize[igsf]                = -999;
	   gsftheta[igsf]                      = -999.0;
	   gsfqualityMask[igsf]                = -999;
	   gsfouterX[igsf]                     = -999.0;
	   gsfouterY[igsf]                     = -999.0;
	   gsfouterZ[igsf]                     = -999.0;
	   gsfouterRadius[igsf]                = -999.0;
	   gsfinnerX[igsf]                     = -999.0;
	   gsfinnerY[igsf]                     = -999.0;
	   gsfinnerZ[igsf]                     = -999.0;
	 }
       gsfsize = -999;
       
       
       for (reco::GsfTrackCollection::const_iterator gsfit = gsfTracksH->begin(); gsfit != gsfTracksH->end(); gsfit++)
	 {
	   gsftrack_container.push_back(*gsfit);
	 }
       //  std::cout << " NTracks " << ntracks << std::endl;
       
       //if(ntracks>0) ntracks5++ ;
       //nevents++;
       
       if((int)gsftrack_container.size()!=0)
	 {
	   for(int igsf = 0; igsf<(int)gsftrack_container.size(); igsf++)
	     {
	       gsfpt[igsf]                         = gsftrack_container[igsf].pt();
	       gsfcharge[igsf]                     = gsftrack_container[igsf].charge();
	       gsfchi2[igsf]                       = gsftrack_container[igsf].chi2();
	       gsfeta[igsf]                        = gsftrack_container[igsf].eta();
	       gsfnumberOfLostHits[igsf]           = gsftrack_container[igsf].numberOfLostHits();
	       gsfnumberOfValidHits[igsf]          = gsftrack_container[igsf].numberOfValidHits();
	       gsflost[igsf]                       = gsftrack_container[igsf].lost();
	       gsfd0[igsf]                         = gsftrack_container[igsf].d0();
	       gsfdxy[igsf]                        = gsftrack_container[igsf].dxy();
	       gsfdz[igsf]                         = gsftrack_container[igsf].dz();
	       gsfptin[igsf]                       = sqrt(gsftrack_container[igsf].innerMomentum().Perp2());
	       gsfptout[igsf]                      = sqrt(gsftrack_container[igsf].outerMomentum().Perp2());
	       gsffbrem[igsf]                      = (gsfptin[igsf]-gsfptout[igsf])/gsfptin[igsf];
	       gsfqoverp[igsf]                     = gsftrack_container[igsf].qoverp();
	       gsfvx[igsf]                         = gsftrack_container[igsf].vx();
	       gsfvy[igsf]                         = gsftrack_container[igsf].vy();
	       gsfvz[igsf]                         = gsftrack_container[igsf].vz();
	       gsfphi[igsf]                        = gsftrack_container[igsf].phi();
	       gsfndof[igsf]                       = gsftrack_container[igsf].ndof();
	       gsfrecHitsSize[igsf]                = gsftrack_container[igsf].recHitsSize();
	       gsftheta[igsf]                      = gsftrack_container[igsf].theta();
	       gsfqualityMask[igsf]                = gsftrack_container[igsf].qualityMask();
	       gsfouterX[igsf]                     = gsftrack_container[igsf].outerX();
	       gsfouterY[igsf]                     = gsftrack_container[igsf].outerY();
	       gsfouterZ[igsf]                     = gsftrack_container[igsf].outerZ();
	       gsfouterRadius[igsf]                = gsftrack_container[igsf].outerRadius();
	       gsfinnerX[igsf]                     = gsftrack_container[igsf].innerPosition().X();
	       gsfinnerY[igsf]                     = gsftrack_container[igsf].innerPosition().Y();
	       gsfinnerZ[igsf]                     = gsftrack_container[igsf].innerPosition().Z();
	     }
	 }//if(gsftrack_container.size!=0)
       gsfsize = (int)gsftrack_container.size();
     }//if(runGsftracks_)       

   //gsfelectron
   if(runGsfelectrons_)
     {
       edm::Handle<reco::GsfElectronCollection> gsfelectronH;
       bool foundele=iEvent.getByLabel("gsfElectrons",gsfelectronH);
       if(!foundele ) {
	 std::ostringstream  err;
	 err<<" cannot get Gsfelectrons:gsfelectrons "
	    <<std::endl;
	 edm::LogError("gsfelectrons")<<err.str();
	 throw cms::Exception( "MissingProduct", err.str());
       }
       
       std::vector<reco::GsfElectron>  gsfelectron_container;
       
       for(int igsf = 0; igsf<25; igsf++)
	 {
	   elept[igsf]                            = -999.0;
	   elefbrem[igsf]                         = -999.0;
	   eletrackerDrivenSeed[igsf]             = -999; 
	   eleecalDrivenSeed[igsf]                = -999;
	   elecharge[igsf]                        = -999; 
	   eledr03EcalRecHitSumEt[igsf]           = -999.0;
	   eledr03HcalDepth1TowerSumEt[igsf]      = -999.0;
	   eledr03HcalDepth2TowerSumEt[igsf]      = -999.0;
	   eledr03HcalTowerSumEt[igsf]            = -999.0;
	   eledr04EcalRecHitSumEt[igsf]           = -999.0;
	   eledr04HcalDepth1TowerSumEt[igsf]      = -999.0;
	   eledr04HcalDepth2TowerSumEt[igsf]      = -999.0;
	   eledr04HcalTowerSumEt[igsf]            = -999.0;
	   elee1x5[igsf]                          = -999.0;
	   elee2x5Max[igsf]                       = -999.0;
	   elee5x5[igsf]                          = -999.0;
	   eleeEleClusterOverPout[igsf]           = -999.0;
	   eleeSeedClusterOverP[igsf]             = -999.0;
	   eleeSeedClusterOverPout[igsf]          = -999.0;
	   eleeSuperClusterOverP[igsf]            = -999.0;
	   eleeta[igsf]                           = -999.0;
	   elehadronicOverEm[igsf]                = -999.0;
	   elesigmaIetaIeta[igsf]                 = -999.0;
	   eledeltaEtaEleClusterTrackAtCalo[igsf] = -999.0;
	   eledeltaEtaSeedClusterTrackAtCalo[igsf]= -999.0;
	   eledeltaEtaSuperClusterTrackAtVtx[igsf]= -999.0;
	   eledeltaPhiEleClusterTrackAtCalo[igsf] = -999.0;
	   eledeltaPhiSeedClusterTrackAtCalo[igsf]= -999.0;
	   eledeltaPhiSuperClusterTrackAtVtx[igsf]= -999.0;
	   eleenergy[igsf]                        = -999.0;
	   elemva[igsf]                           = -999.0;
	   elenumberOfTracks[igsf]                = -999;
	   elemaxEnergyXtal[igsf]                 = -999.;
	   eleswissCross[igsf]                    = -999.;
	   eleswissBasedspikevar[igsf]            = -999.;
	   
	   
	   //gsf
	   eletrkpt[igsf]                         = -999.0;
	   eletrkcharge[igsf]                     = -999;
	   eletrkchi2[igsf]                       = -999.0;
	   eletrketa[igsf]                        = -999.0;
	   eletrknumberOfLostHits[igsf]           = -999;
	   eletrknumberOfValidHits[igsf]          = -999;
	   eletrklost[igsf]                       = -999;
	   eletrkd0[igsf]                         = -999.0;
	   eletrkdxy[igsf]                        = -999.0;
	   eletrkdz[igsf]                         = -999.0;
	   eletrkptin[igsf]                       = -999.0;
	   eletrkptout[igsf]                      = -999.0;
	   eletrkfbrem[igsf]                      = -999.0;
	   eletrkqoverp[igsf]                     = -999.0;
	   eletrkvx[igsf]                         = -999.0;
	   eletrkvy[igsf]                         = -999.0;
	   eletrkvz[igsf]                         = -999.0;
	   eletrkphi[igsf]                        = -999.0;
	   eletrkndof[igsf]                       = -999.0;
	   eletrkrecHitsSize[igsf]                = -999;
	   eletrktheta[igsf]                      = -999.0;
	   eletrkqualityMask[igsf]                = -999;
	   eletrkouterX[igsf]                     = -999.0;
	   eletrkouterY[igsf]                     = -999.0;
	   eletrkouterZ[igsf]                     = -999.0;
	   eletrkouterRadius[igsf]                = -999.0;
	   eletrkinnerX[igsf]                     = -999.0;
	   eletrkinnerY[igsf]                     = -999.0;
	   eletrkinnerZ[igsf]                     = -999.0;
	   
	   //ctf
	   electfpt[igsf]                         = -999.0;
	   electfcharge[igsf]                     = -999;
	   electfchi2[igsf]                       = -999.0;
	   electfeta[igsf]                        = -999.0;
	   electfnumberOfLostHits[igsf]           = -999;
	   electfnumberOfValidHits[igsf]          = -999;
	   electflost[igsf]                       = -999;
	   electfd0[igsf]                         = -999.0;
	   electfdxy[igsf]                        = -999.0;
	   electfdz[igsf]                         = -999.0;
	   electfptin[igsf]                       = -999.0;
	   electfptout[igsf]                      = -999.0;
	   electffbrem[igsf]                      = -999.0;
	   electfqoverp[igsf]                     = -999.0;
	   electfvx[igsf]                         = -999.0;
	   electfvy[igsf]                         = -999.0;
	   electfvz[igsf]                         = -999.0;
	   electfphi[igsf]                        = -999.0;
	   electfndof[igsf]                       = -999.0;
	   electfrecHitsSize[igsf]                = -999;
	   electftheta[igsf]                      = -999.0;
	   electfqualityMask[igsf]                = -999;
	   electfouterX[igsf]                     = -999.0;
	   electfouterY[igsf]                     = -999.0;
	   electfouterZ[igsf]                     = -999.0;
	   electfouterRadius[igsf]                = -999.0;
	   electfinnerX[igsf]                     = -999.0;
	   electfinnerY[igsf]                     = -999.0;
	   electfinnerZ[igsf]                     = -999.0;
	   
	 }
       
       
       
       
       for (reco::GsfElectronCollection::const_iterator gsfeit = gsfelectronH->begin(); gsfeit != gsfelectronH->end(); gsfeit++)
	 {
	   gsfelectron_container.push_back(*gsfeit);
	   //cout<<"H/E:"<<gsfeit->hadronicOverEm()<<endl;
	 }
       
       if(gsfelectron_container.size()!=0)
	 {
	   for(int igsfel=0; igsfel<(int)gsfelectron_container.size(); igsfel++)
	     {
	       
	       elept[igsfel]                            = gsfelectron_container[igsfel].pt();
	       elepx[igsfel]                            = gsfelectron_container[igsfel].px();
               elepy[igsfel]                            = gsfelectron_container[igsfel].py();
               elepz[igsfel]                            = gsfelectron_container[igsfel].pz();
               elephi[igsfel]                           = correct_phi(gsfelectron_container[igsfel].phi());
	       eletheta[igsfel]                         = gsfelectron_container[igsfel].theta();


	       elefbrem[igsfel]                         = gsfelectron_container[igsfel].fbrem();
	       eletrackerDrivenSeed[igsfel]             = gsfelectron_container[igsfel].trackerDrivenSeed();
	       eleecalDrivenSeed[igsfel]                = gsfelectron_container[igsfel].ecalDrivenSeed();
	       elecharge[igsfel]                        = gsfelectron_container[igsfel].charge();
	       eledr03EcalRecHitSumEt[igsfel]           = gsfelectron_container[igsfel].dr03EcalRecHitSumEt();
	       eledr03HcalDepth1TowerSumEt[igsfel]      = gsfelectron_container[igsfel].dr03HcalDepth1TowerSumEt();
	       eledr03HcalDepth2TowerSumEt[igsfel]      = gsfelectron_container[igsfel].dr03HcalDepth2TowerSumEt();
	       eledr03HcalTowerSumEt[igsfel]            = gsfelectron_container[igsfel].dr03HcalTowerSumEt();
	       eledr04EcalRecHitSumEt[igsfel]           = gsfelectron_container[igsfel].dr04EcalRecHitSumEt();
	       eledr04HcalDepth1TowerSumEt[igsfel]      = gsfelectron_container[igsfel].dr04HcalDepth1TowerSumEt();
	       eledr04HcalDepth2TowerSumEt[igsfel]      = gsfelectron_container[igsfel].dr04HcalDepth2TowerSumEt();
	       eledr04HcalTowerSumEt[igsfel]            = gsfelectron_container[igsfel].dr04HcalTowerSumEt();
	       elee1x5[igsfel]                          = gsfelectron_container[igsfel].e1x5();
	       elee2x5Max[igsfel]                       = gsfelectron_container[igsfel].e2x5Max();
	       elee5x5[igsfel]                          = gsfelectron_container[igsfel].e5x5();
	       eleeEleClusterOverPout[igsfel]           = gsfelectron_container[igsfel].eEleClusterOverPout();
	       eleeSeedClusterOverP[igsfel]             = gsfelectron_container[igsfel].eSeedClusterOverP();
	       eleeSeedClusterOverPout[igsfel]          = gsfelectron_container[igsfel].eSeedClusterOverPout();
	       eleeSuperClusterOverP[igsfel]            = gsfelectron_container[igsfel].eSuperClusterOverP();
	       eleeta[igsfel]                           = gsfelectron_container[igsfel].eta();
	       elehadronicOverEm[igsfel]                = gsfelectron_container[igsfel].hadronicOverEm();
	       //cout<<"inside gsf electrons:elehadronicOverEm[igsfel]:"<<elehadronicOverEm[igsfel]<<":gsfelectron_container[igsfel].hadronicOverEm():"<<gsfelectron_container[igsfel].hadronicOverEm()<<endl;
	       elesigmaIetaIeta[igsfel]                 = gsfelectron_container[igsfel].sigmaIetaIeta();
	       eledr03TkSumPt[igsfel]                   = gsfelectron_container[igsfel].dr03TkSumPt();
	       eledr04TkSumPt[igsfel]                   = gsfelectron_container[igsfel].dr04TkSumPt();
	       
	       eledeltaEtaEleClusterTrackAtCalo[igsfel] = gsfelectron_container[igsfel].deltaEtaEleClusterTrackAtCalo();
	       eledeltaEtaSeedClusterTrackAtCalo[igsfel]= gsfelectron_container[igsfel].deltaEtaSeedClusterTrackAtCalo();
	       eledeltaEtaSuperClusterTrackAtVtx[igsfel]= gsfelectron_container[igsfel].deltaEtaSuperClusterTrackAtVtx();
	       eledeltaPhiEleClusterTrackAtCalo[igsfel] = gsfelectron_container[igsfel].deltaPhiEleClusterTrackAtCalo();
	       eledeltaPhiSeedClusterTrackAtCalo[igsfel]= gsfelectron_container[igsfel].deltaPhiSeedClusterTrackAtCalo();
	       eledeltaPhiSuperClusterTrackAtVtx[igsfel]= gsfelectron_container[igsfel].deltaPhiSuperClusterTrackAtVtx();
	       eleenergy[igsfel]                        = gsfelectron_container[igsfel].energy();
	       elemva[igsfel]                           = gsfelectron_container[igsfel].mva();
	       elenumberOfTracks[igsfel]                = gsfelectron_container[igsfel].numberOfTracks();
	       //std::cout<<"eletrackerDrivenSeed[igsfel] = "<<eletrackerDrivenSeed[igsfel]<<std::endl;
	       
	       //sceta
	       reco::SuperClusterRef scref = gsfelectron_container[igsfel].superCluster();
	       
	       elesceta[igsfel]                         = scref->eta();
	       elescphi[igsfel]                         = scref->phi();
	       elescE[igsfel]                           = scref->energy();
	       //elescpt[igsfel]                          = scref->pt();
	       //maximum energy xtal
	       const reco::CaloClusterPtr seed = scref->seed();
	       elemaxEnergyXtal[igsfel]                 = lazyTool.eMax( *seed );
	       
	       
	       //PF cluster
	       reco::SuperClusterRef pfref = gsfelectron_container[igsfel].pflowSuperCluster(); 
	       if(gsfelectron_container[igsfel].pflowSuperCluster().isNonnull()) 
		 {
		   elepfeta[igsfel]                         = pfref->eta();
		   elepfphi[igsfel]                         = pfref->phi();
		   elepfE[igsfel]                           = pfref->energy();
		   //elepfpt[igsfel]                          = pfref->pt();
		 }
	       
	     
	       //trackinfo
	       /*reco::GsfTrackRef trackref = gsfelectron_container[igsfel].gsfTrack();
	       if(gsfelectron_container[igsfel].gsfTrack().isNonnull())
		 {
		   eletrkpt[igsfel]                         = trackref->pt();
		   eletrkcharge[igsfel]                     = trackref->charge();
		   eletrkchi2[igsfel]                       = trackref->chi2();
		   eletrketa[igsfel]                        = trackref->eta();
		   eletrknumberOfLostHits[igsfel]           = trackref->numberOfLostHits();
		   eletrknumberOfValidHits[igsfel]          = trackref->numberOfValidHits();
		   eletrklost[igsfel]                       = trackref->lost();
		   eletrkd0[igsfel]                         = trackref->d0();
		   eletrkdxy[igsfel]                        = trackref->dxy();
		   eletrkdz[igsfel]                         = trackref->dz();
		   eletrkptin[igsfel]                       = sqrt(trackref->innerMomentum().Perp2());
		   eletrkptout[igsfel]                      = sqrt(trackref->outerMomentum().Perp2());
		   eletrkfbrem[igsfel]                      = (eletrkptin[igsfel]-eletrkptout[igsfel])/eletrkptin[igsfel];
		   eletrkqoverp[igsfel]                     = trackref->qoverp();
		   eletrkvx[igsfel]                         = trackref->vx();
		   eletrkvy[igsfel]                         = trackref->vy();
		   eletrkvz[igsfel]                         = trackref->vz();
		   eletrkphi[igsfel]                        = trackref->phi();
		   eletrkndof[igsfel]                       = trackref->ndof();
		   eletrkrecHitsSize[igsfel]                = trackref->recHitsSize();
		   eletrktheta[igsfel]                      = trackref->theta();
		   eletrkqualityMask[igsfel]                = trackref->qualityMask();
		   eletrkouterX[igsfel]                     = trackref->outerX();
		   eletrkouterY[igsfel]                     = trackref->outerY();
		   eletrkouterZ[igsfel]                     = trackref->outerZ();
		   eletrkouterRadius[igsfel]                = trackref->outerRadius();
		   eletrkinnerX[igsfel]                     = trackref->innerPosition().X();
		   eletrkinnerY[igsfel]                     = trackref->innerPosition().Y();
		   eletrkinnerZ[igsfel]                     = trackref->innerPosition().Z();
		 }//if(gsfelectron_container[igsfel].gsfTrack().isNonnull())
	       
	       reco::TrackRef ctfref = gsfelectron_container[igsfel].closestCtfTrackRef();
	       if(gsfelectron_container[igsfel].closestCtfTrackRef().isNonnull())
		 {
		   electfpt[igsfel]                         = ctfref->pt();
		   electfcharge[igsfel]                     = ctfref->charge();
		   electfchi2[igsfel]                       = ctfref->chi2();
		   electfeta[igsfel]                        = ctfref->eta();
		   electfnumberOfLostHits[igsfel]           = ctfref->numberOfLostHits();
		   electfnumberOfValidHits[igsfel]          = ctfref->numberOfValidHits();
		   electflost[igsfel]                       = ctfref->lost();
		   electfd0[igsfel]                         = ctfref->d0();
		   electfdxy[igsfel]                        = ctfref->dxy();
		   electfdz[igsfel]                         = ctfref->dz();
		   electfptin[igsfel]                       = sqrt(ctfref->innerMomentum().Perp2());
		   electfptout[igsfel]                      = sqrt(ctfref->outerMomentum().Perp2());
		   electffbrem[igsfel]                      = (electfptin[igsfel]-electfptout[igsfel])/electfptin[igsfel];
		   electfqoverp[igsfel]                     = ctfref->qoverp();
		   electfvx[igsfel]                         = ctfref->vx();
		   electfvy[igsfel]                         = ctfref->vy();
		   electfvz[igsfel]                         = ctfref->vz();
		   electfphi[igsfel]                        = ctfref->phi();
		   electfndof[igsfel]                       = ctfref->ndof();
		   electfrecHitsSize[igsfel]                = ctfref->recHitsSize();
		   electftheta[igsfel]                      = ctfref->theta();
		   electfqualityMask[igsfel]                = ctfref->qualityMask();
		   electfouterX[igsfel]                     = ctfref->outerX();
		   electfouterY[igsfel]                     = ctfref->outerY();
		   electfouterZ[igsfel]                     = ctfref->outerZ();
		   electfouterRadius[igsfel]                = ctfref->outerRadius();
		   electfinnerX[igsfel]                     = ctfref->innerPosition().X();
		   electfinnerY[igsfel]                     = ctfref->innerPosition().Y();
		   electfinnerZ[igsfel]                     = ctfref->innerPosition().Z();
		 }//if(gsfelectron_container[igsfel].closestCtfTrackRef().isNonnull())
	       */	       
	       //spike info
	       if(gsfelectron_container[igsfel].isEB())
		 {
		   eleswissCross[igsfel]   = EcalClusterTools::eTop( *(gsfelectron_container[igsfel].superCluster()->seed()), &(*barrelRecHits), &(*topology))+ EcalClusterTools::eBottom( *(gsfelectron_container[igsfel].superCluster()->seed()), &(*barrelRecHits), &(*topology)) + EcalClusterTools::eLeft( *(gsfelectron_container[igsfel].superCluster()->seed()), &(*barrelRecHits), &(*topology)) + EcalClusterTools::eRight( *(gsfelectron_container[igsfel].superCluster()->seed()), &(*barrelRecHits), &(*topology));
		   //std::cout<<"etop: "<<EcalClusterTools::eTop( *(gsfelectron_container[igsfel].superCluster()->seed()), &(*barrelRecHits), &(*topology))<<endl;
		   //std::cout<<"etop: "<<EcalClusterTools::eTop( *(gsfelectron_container[igsfel].superCluster()->seed()), &(*barrelRecHits), &(*topology))<<endl;
		   //if(1-eleswissCross[igsfel]/elemaxEnergyXtal[igsfel] > 0.95) cout<<"This electron candidate is an ECAL spike identified by Swiss Cross algorithm."<<endl;
		   
		   eleswissBasedspikevar[igsfel]             = 1-eleswissCross[igsfel]/elemaxEnergyXtal[igsfel];
		 }//end of if(gsfelectron_container[igsfel].isEB())
	       else{
		 eleswissCross[igsfel]   = EcalClusterTools::eTop( *(gsfelectron_container[igsfel].superCluster()->seed()), &(*endcapRecHits), &(*topology))+ EcalClusterTools::eBottom( *(gsfelectron_container[igsfel].superCluster()->seed()), &(*endcapRecHits), &(*topology)) + EcalClusterTools::eLeft( *(gsfelectron_container[igsfel].superCluster()->seed()), &(*endcapRecHits), &(*topology)) + EcalClusterTools::eRight( *(gsfelectron_container[igsfel].superCluster()->seed()), &(*endcapRecHits), &(*topology));
		 //if(1-eleswissCross[igsfel]/elemaxEnergyXtal[igsfel] > 0.95) cout<<"This electron candidate is an ECAL spike identified by Swiss Cross algorithm."<<endl;
		 eleswissBasedspikevar[igsfel]             = 1-eleswissCross[igsfel]/elemaxEnergyXtal[igsfel];
	       }//end of else
	       
	       //timing info
               DetId eleSeedDetId = lazyTool.getMaximum(*seed).first;

               if ( gsfelectron_container[igsfel].isEB() ) {
		 EcalRecHitCollection::const_iterator eleebrhit = barrelRecHits->find(eleSeedDetId);
                 if ( eleebrhit != barrelRecHits->end() ) {
                   eleseedtime[igsfel]         = eleebrhit->time();
                   elerecoFlag[igsfel]         = eleebrhit->recoFlag();

                   //eleseverityLevel[igsfel]    = EcalSeverityLevelAlgo::severityLevel( eleSeedDetId, (*barrelRecHits), *chStatus );
                 }//if ( eleebrhit != barrelRecHits->end() )
               }//if ( electron_container[igsfel].isEB() && barrelRecHits.isValid() )

               else {
		 EcalRecHitCollection::const_iterator eleeerhit = endcapRecHits->find(eleSeedDetId);
                 if ( eleeerhit != endcapRecHits->end() ) {
                   eleseedtime[igsfel]         = eleeerhit->time();
                   elerecoFlag[igsfel]         = eleeerhit->recoFlag();
                   //eleseverityLevel[igsfel]    = EcalSeverityLevelAlgo::severityLevel( eleSeedDetId, (*endcapRecHits), *chStatus );
                 }//if ( eleeerhit != endcapRecHits->end() )
               }//else if ( endcapRecHits.isValid() )
	     	       
	       
	     }//for(int igsfel=0; igsfel<gsfelectron_container.size(); igsfel++)
	   
	 }//if(gsfelectron_container.size()!=0)
       
       elesize = (int)gsfelectron_container.size();
     }//if(runGsfelectrons_)

   //general tracks
   if(runGeneraltracks_)
     {
       edm::Handle<reco::TrackCollection> gentracksH;
       bool foundtrk=iEvent.getByLabel("generalTracks",gentracksH);
       if(!foundtrk) {
	 std::ostringstream  err;
	 err<<" cannot get GeneralTracks:generaltracks "
	    <<std::endl;
	 edm::LogError("gsfTracks")<<err.str();
	 throw cms::Exception( "MissingProduct", err.str());
       }
       
       for(int igentrk = 0; igentrk<1000; igentrk++)
	 {
	   //gen tracks
	   gentrkcharge[igentrk]                     = -999;
	   gentrkchi2[igentrk]                       = -999.0;
	   gentrknumberOfLostHits[igentrk]           = -999;
	   gentrknumberOfValidHits[igentrk]          = -999;
	   gentrklost[igentrk]                       = -999;
	   gentrkd0[igentrk]                         = -999.0;
	   gentrkdxy[igentrk]                        = -999.0;
	   gentrkdz[igentrk]                         = -999.0;
	   gentrkptin[igentrk]                       = -999.0;
	   gentrkptout[igentrk]                      = -999.0;
	   gentrkfbrem[igentrk]                      = -999.0;
	   gentrkqoverp[igentrk]                     = -999.0;
	   gentrkndof[igentrk]                       = -999.0;
	   gentrkrecHitsSize[igentrk]                = -999;
	   gentrkqualityMask[igentrk]                = -999;
	   gentrkouterX[igentrk]                     = -999.0;
	   gentrkouterY[igentrk]                     = -999.0;
	   gentrkouterZ[igentrk]                     = -999.0;
	   gentrkouterRadius[igentrk]                = -999.0;
	   gentrkinnerX[igentrk]                     = -999.0;
	   gentrkinnerY[igentrk]                     = -999.0;
	   gentrkinnerZ[igentrk]                     = -999.0;
	   
	 }
       
       
       
       gentrkp4->Clear();
       gentrk_vtxpos->Clear();

       int igentrk = 0;
       
       for (reco::TrackCollection::const_iterator trk = gentracksH->begin(); trk != gentracksH->end(); trk++)
	 {
	   //trackinfo
	   //cout<<"inside tracker"<<endl;
	   new ((*gentrkp4)[igentrk]) TLorentzVector();
	   ((TLorentzVector *)gentrkp4->At(igentrk))->SetXYZT(trk->px(), trk->py(), trk->pz(), trk->p());

	   new ((*gentrk_vtxpos)[igentrk]) TVector3();
	   ((TVector3 *)gentrk_vtxpos->At(igentrk))->SetXYZ(trk->vx(), trk->vy(), trk->vz());
	   
	   gentrkcharge[igentrk]                     = trk->charge();
	   //cout<<"trk->charge() = "<<trk->charge()<<endl;
	   gentrkchi2[igentrk]                       = trk->chi2();
	   //cout<<"trk->chi2() "<<trk->chi2()<<endl;
	   gentrknumberOfLostHits[igentrk]           = trk->numberOfLostHits();
	   //cout<<"trk->numberOfLostHits() "<<trk->numberOfLostHits()<<endl;
	   gentrknumberOfValidHits[igentrk]          = trk->numberOfValidHits();
	   //cout<<"trk->numberOfValidHits() "<<trk->numberOfValidHits()<<endl;
	   gentrklost[igentrk]                       = trk->lost();
	   //cout<<"trk->lost() "<<trk->lost()<<endl;
	   gentrkd0[igentrk]                         = trk->d0();
	   //cout<<"trk->d0() "<<trk->d0()<<endl;
	   gentrkdxy[igentrk]                        = trk->dxy();
	   //cout<<"trk->dxy() "<<trk->dxy()<<endl;
	   gentrkdz[igentrk]                         = trk->dz();
	   //cout<<"trk->dz() "<<trk->dz()<<endl;
	   
	   /*gentrkptin[igentrk]                       = sqrt(trk->innerMomentum().Perp2());
	   cout<<"sqrt(trk->innerMomentum().Perp2()) "<<sqrt(trk->innerMomentum().Perp2())<<endl;
	   gentrkptout[igentrk]                      = sqrt(trk->outerMomentum().Perp2());
	   cout<<"sqrt(trk->outerMomentum().Perp2()) "<<sqrt(trk->outerMomentum().Perp2())<<endl;
	   gentrkfbrem[igentrk]                      = (gentrkptin[igentrk]-gentrkptout[igentrk])/gentrkptin[igentrk];
	   cout<<"(gentrkptin[igentrk]-gentrkptout[igentrk])/gentrkptin[igentrk] "<<(gentrkptin[igentrk]-gentrkptout[igentrk])/gentrkptin[igentrk]<<endl;
	   */

	   gentrkqoverp[igentrk]                     = trk->qoverp();
	   //cout<<"trk->qoverp() "<<trk->qoverp()<<endl;
	   gentrkndof[igentrk]                       = trk->ndof();
	   //cout<<"trk->ndof() "<<trk->ndof()<<endl;
	   
	   //gentrkrecHitsSize[igentrk]                = trk->recHitsSize();
	   //cout<<"trk->recHitsSize() "<<trk->recHitsSize()<<endl;
	   
	   gentrkqualityMask[igentrk]                = trk->qualityMask();
	   //cout<<"trk->qualityMask() "<<trk->qualityMask()<<endl;
	   
	   /*gentrkouterX[igentrk]                     = trk->outerX();
	   cout<<"trk->outerX() "<<trk->outerX()<<endl;
	   gentrkouterY[igentrk]                     = trk->outerY();
	   cout<<"trk->outerY() "<<trk->outerY()<<endl;
	   gentrkouterZ[igentrk]                     = trk->outerZ();
	   cout<<"trk->outerZ() "<<trk->outerZ()<<endl;
	   gentrkouterRadius[igentrk]                = trk->outerRadius();
	   cout<<"trk->outerRadius() "<<trk->outerRadius()<<endl;
	   gentrkinnerX[igentrk]                     = trk->innerPosition().X();
	   cout<<"trk->innerPosition().X() "<<trk->innerPosition().X()<<endl;
	   gentrkinnerY[igentrk]                     = trk->innerPosition().Y();
	   cout<<"trk->innerPosition().X() "<<trk->innerPosition().X()<<endl;
	   gentrkinnerZ[igentrk]                     = trk->innerPosition().Z();
	   cout<<"trk->innerPosition().X() "<<trk->innerPosition().X()<<endl;
	   */
	   
	   gentrk_hpnvalid[igentrk]                  = trk->hitPattern().numberOfValidHits();
	   //cout<<"trk->hitPattern().numberOfValidHits() = "<<trk->hitPattern().numberOfValidHits()<<endl;
	   gentrk_hpnlost[igentrk]                   = trk->hitPattern().numberOfLostHits();
	   //cout<<"trk->hitPattern().numberOfLostHits() "<<trk->hitPattern().numberOfLostHits()<<endl;
	   gentrk_hpnvalidpix[igentrk]               = trk->hitPattern().numberOfValidPixelHits();
	   //cout<<"trk->hitPattern().numberOfValidPixelHits() "<<trk->hitPattern().numberOfValidPixelHits()<<endl;
	   
	   
	   igentrk++;
	 }//for (reco::TrackCollection::const_iterator trk = gentracksH->begin(); gentrkit != gentracksH->end(); gentrkit++)
       gentrksize = igentrk;
       //cout<<"outside trk loop"<<endl;
     }//if(runGeneraltracks_)
       
   
   ///vertices info
   /*if(runVertices_)
     {
       edm::Handle<reco::VertexCollection> vtxH;
       bool foundvtx=iEvent.getByLabel("offlinePrimaryVerticesWithBS",vtxH);
       if(!foundvtx) {
	 std::ostringstream  err;
         err<<" cannot get VertexCollection:offlinePrimaryVerticesWithBS "
            <<std::endl;
	 edm::LogError("VertexCollection")<<err.str();
         throw cms::Exception( "MissingProduct", err.str());
       }//if(!foundvtx)
       
       edm::Handle<reco::TrackCollection> tkH;
       iEvent.getByLabel("generalTracks", tkH);

       vtx_xyz->Clear();
       vtx_dxdydz->Clear();
       vtx_vectorp3->Clear();
       
       int ivtx = 0;
       for (reco::VertexCollection::const_iterator vtxk = vtxH->begin(); vtx != vtxH->end(); vtx++)
	 {
	   new((*vtx_xyz)[vtx_n]) TVector3();
	   ((TVector3 *) vtx_xyz->At(vtx_n))->SetXYZ(vtx->x(), vtx->y(), vtx->z());

	   new((*vtx_dxdydz)[vtx_n]) TVector3();
	   ((TVector3 *)vtx_dxdydz->At(vtx_n))->SetXYZ(vtx->xError(), vtx->yError(), vtx->zError());

	   math::XYZVector vecsumP3 = math::XYZVector(0.,0.,0.);
	   vtx_scalarpt[vtx_n]=0.;

	   for(std::vector<reco::TrackBaseRef>::const_iterator trkItr = vtx->tracks_begin();trkItr != vtx->tracks_end(); ++trkItr) {
	     vtx_scalarpt[vtx_n] += (*trkItr)->pt();
	     vecsumP3 += (*trkItr)->momentum();
	   }

	   new((*vtx_vectorp3)[vtx_n]) TVector3();
	   ((TVector3 *)vtx_vectorp3->At(vtx_n))->SetXYZ(vecsumP3.x(), vecsumP3.y(), vecsumP3.z());

	   

	 }//for (reco::VertexCollection::const_iterator vtxk = vtxH->begin(); vtx != vtxH->end(); vtx++)
     }//if(runVertices_)
   */



   //ES info
   if(runPreshower_)
     {
       //X ESclusters
       edm::Handle<reco::PreshowerClusterCollection> psclusCollX;
       iEvent.getByLabel("multi5x5SuperClustersWithPreshower","preshowerXClusters",psclusCollX);
       if (!psclusCollX.isValid()) {
	 edm::LogError("PSClustercollection") << "Error! can't get collection with label : multi5x5SuperClustersWithPreshower";
	 
       }
       
       
       //Y ESclusters      
       edm::Handle<reco::PreshowerClusterCollection> psclusCollY;
       iEvent.getByLabel("multi5x5SuperClustersWithPreshower","preshowerYClusters",psclusCollY);
       if (!psclusCollY.isValid()) {
	 edm::LogError("PSClustercollection") << "Error! can't get collection with label : multi5x5SuperClustersWithPreshower";
	 
       }
       
       //ES rechits
       // Get ES rechits. select ES rechits of those SC for which SC raw Et>2GeV and lie within the ES region
       edm::Handle<EcalRecHitCollection> PreshowerRecHits;
       iEvent.getByLabel(InputTag("ecalPreshowerRecHit","EcalRecHitsES"), PreshowerRecHits);
       if( PreshowerRecHits.isValid() ) EcalRecHitCollection preshowerHits(*PreshowerRecHits);
       
       //get Geometry
       ESHandle<CaloGeometry> caloGeometry;
       iSetup.get<CaloGeometryRecord>().get(caloGeometry);
       const CaloSubdetectorGeometry *geometry = caloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
       const CaloSubdetectorGeometry *& geometry_p = geometry;
       
       
       edm::Handle<reco::SuperClusterCollection> pEndcapRawSuperClusters;
       iEvent.getByLabel("multi5x5SuperClusters","multi5x5EndcapSuperClusters",pEndcapRawSuperClusters);
       if (!pEndcapRawSuperClusters.isValid()) {
	 std::cout<< "Error! can't get collection with label multi5x5SuperClusters"<<std::endl;
	 
       }
       
       const reco::SuperClusterCollection* endcapRawSuperClusters = pEndcapRawSuperClusters.product();
       
       
       int isc = 0;
       int esrechitsize;
       //ES clus size
       int clussizeX;
       int clussizeY;
       
       for(reco::SuperClusterCollection::const_iterator aClus = endcapRawSuperClusters->begin();aClus != endcapRawSuperClusters->end(); aClus++)
	 {
	   
	   double sceta = aClus->position().eta();
	   double scet  = aClus->energy()/std::cosh(aClus->position().eta());
	   
	   if( fabs(sceta)>1.65 && fabs(sceta)<2.5 && scet>2.0 )
	     {
	       esrechitsize  = 0;
	       clussizeX = 0;
	       clussizeY = 0;
	       
	       EEsceta[isc]        = aClus->position().eta();
	       EEscet[isc]         = aClus->energy()/std::cosh(aClus->position().eta());
	       
	       reco::CaloCluster_iterator bc_iter = aClus->clustersBegin();
	       for ( ; bc_iter !=aClus->clustersEnd(); ++bc_iter ) { 
		 if (geometry)
		   {
		     // Get strip position at intersection point of the line EE - Vertex:
		     double X = (*bc_iter)->x();
		     double Y = (*bc_iter)->y();
		     double Z = (*bc_iter)->z();        
		     const GlobalPoint point(X,Y,Z);    
		     
		     DetId tmp1 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(point, 1);
		     DetId tmp2 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(point, 2);
		     ESDetId esfid = (tmp1 == DetId(0)) ? ESDetId(0) : ESDetId(tmp1);
		     ESDetId esrid = (tmp2 == DetId(0)) ? ESDetId(0) : ESDetId(tmp2);     
		     
		     
		     const ESRecHitCollection *ESRH = PreshowerRecHits.product();
		     
		     //std::cout<<"ESrechit size = "<<ESrechitsize<<std::endl;
		     EcalRecHitCollection::const_iterator esrh_it;
		     for ( esrh_it = ESRH->begin(); esrh_it != ESRH->end(); esrh_it++) {
		       ESDetId esdetid = ESDetId(esrh_it->id());
		       if ( esdetid.plane()==1 ) 
			 {
			   if ( esdetid.zside() == esfid.zside() && esdetid.siy() == esfid.siy() )
			     {
			       ESrechitE[isc][esrechitsize]          = esrh_it->energy();
			       ESrechitplane[isc][esrechitsize]      = esdetid.plane();
			       ESrechitzside[isc][esrechitsize]      = esdetid.zside();
			       ESrechitsix[isc][esrechitsize]        = esdetid.six();
			       ESrechitsiy[isc][esrechitsize]        = esdetid.siy();
			       
			       esrechitsize++ ; 
			     }//if ( esdetid.zside() == esfid.zside() && esdetid.siy() == esfid.siy() )
			 }//if ( esdetid.plane()==1 )
		       
		       if ( esdetid.plane()==2 )
			 {
			   if ( esdetid.zside() == esfid.zside() && esdetid.six() == esfid.six() )
			     {
			       ESrechitE[isc][esrechitsize]          = esrh_it->energy();
			       ESrechitplane[isc][esrechitsize]      = esdetid.plane();
			       ESrechitzside[isc][esrechitsize]      = esdetid.zside();
			       ESrechitsix[isc][esrechitsize]        = esdetid.six();
			       ESrechitsiy[isc][esrechitsize]        = esdetid.siy();
			       
			       esrechitsize++ ;
			     }//if ( esdetid.zside() == esfid.zside() && esdetid.six() == esfid.six() )
			 }//if ( esdetid.plane()==2 ) 
		     }//for ( esrh_it = ESRH->begin(); esrh_it != ESRH->end(); esrh_it++)
		     
		     
		     
		     //ES cluster info
		     const reco::CaloClusterPtr ecalBasicClusterPtr = *(bc_iter);
		     
		     //std::cout<<"PS cluster size = "<<psclussize<<std::endl;
		     
		     for(reco::PreshowerClusterCollection::const_iterator pClus = psclusCollX->begin();pClus != psclusCollX->end(); pClus++)
		       {
			 const reco::CaloClusterPtr preshBasicCluster = pClus->basicCluster();
			 if (preshBasicCluster == ecalBasicClusterPtr) {
			   ESclusEtX[isc][clussizeX]              = pClus->et();
			   ESclusEX[isc][clussizeX]               = pClus->energy();
			   ESclusetaX[isc][clussizeX]             = pClus->eta();
			   ESclusphiX[isc][clussizeX]             = pClus->phi();
			   ESclusnhitsX[isc][clussizeX]           = pClus->nhits();
			   
			   clussizeX++;
			   // std::cout<<"PS et      = "<<et<<std::endl;
			   //std::cout<<"PS energy  = "<<energy<<std::endl;
			 }//if (preshBasicCluster == ecalBasicClusterPtr)
		       }//for(reco::PreshowerClusterCollection::const_iterator pClus = psclusCollX->begin();pClus != psclusCollX->end(); pClus++)
		     
		     
		     //std::cout<<"PS cluster size = "<<psclussize<<std::endl;
		     double eY = 0;
		     
		     for(reco::PreshowerClusterCollection::const_iterator pClus = psclusCollY->begin();pClus != psclusCollY->end(); pClus++)
		       {
			 const reco::CaloClusterPtr preshBasicCluster = pClus->basicCluster();
			 if (preshBasicCluster == ecalBasicClusterPtr) {
			   ESclusEtY[isc][clussizeY]              = pClus->et();
			   ESclusEY[isc][clussizeY]               = pClus->energy();
			   ESclusetaY[isc][clussizeY]             = pClus->eta();
			   ESclusphiY[isc][clussizeY]             = pClus->phi();
			   ESclusnhitsY[isc][clussizeY]           = pClus->nhits();
			   clussizeY++;
			 }//if (preshBasicCluster == ecalBasicClusterPtr)
		       }//for(reco::PreshowerClusterCollection::const_iterator pClus = psclusCollY->begin();pClus != psclusCollY->end(); pClus++)
		     
		   }//if(geomentry)
	       }//for ( ; bc_iter !=it_super->clustersEnd(); ++bc_iter )
	       ESrechitsize[isc] = esrechitsize;
	       ESclussizeX[isc]  = clussizeX;
	       ESclussizeY[isc]  = clussizeY;
	       //std::cout<<"ESrechitsize["<<isc<<"] = "<<ESrechitsize[isc]<<std::endl;
	       //std::cout<<"ESclussizeX["<<isc<<"] = "<<ESclussizeX[isc]<<std::endl;
	       //std::cout<<"ESclussizeY["<<isc<<"] = "<<ESclussizeY[isc]<<std::endl;
	       isc++;
	     }//if( fabs(sceta[isc])>1.65 && fabs(sceta[isc])<2.5 && scet[isc]>2.0 )
	   
	 }//for(reco::SuperClusterCollection::const_iterator aClus = endcapRawSuperClusters->begin();aClus != endcapRawSuperClusters->end(); aClus++)
       
       EEscsize = isc;
       //std::cout<<"EEscsize = "<<EEscsize<<std::endl;
     }//if(runPreshower_)

   //pat info
   
   //electrons
   if(runPatelectrons_)
     {
       edm::Handle<edm::View<pat::Electron> > electronHandle;
       try{iEvent.getByLabel(eleLabel_,electronHandle); 
	 if(debugEle_){
	   std::cout<<"Found the handle"<<std::endl;                                                                                                           
	 }
       }
       catch(...){std::cout<<"Didn't fine the pat::electron handle"<<std::endl;}
       const edm::View<pat::Electron> & electrons = *electronHandle;

       ///for dxy
       edm::Handle<reco::VertexCollection> vtx;
       iEvent.getByLabel(vertexLabel_, vtx);



       //barrel sc's   
       edm::Handle<reco::SuperClusterCollection> superClustersHybridH;
       iEvent.getByLabel(hybridSuperClusterColl_,superClustersHybridH);

       //EE sc's     
       edm::Handle<reco::SuperClusterCollection> superClustersEndcapH;
       iEvent.getByLabel(endcapSuperClusterColl_, superClustersEndcapH);

       for(int ipat = 0; ipat<MAX_ELECTRONS; ipat++)
	 {
	   patelept[ipat]                            = -999.0;
	   patelefbrem[ipat]                         = -999.0;
	   pateletrackerDrivenSeed[ipat]             = -999; 
	   pateleecalDrivenSeed[ipat]                = -999;
	   patelecharge[ipat]                        = -999; 
	   pateledr03EcalRecHitSumEt[ipat]           = -999.0;
	   pateledr03HcalDepth1TowerSumEt[ipat]      = -999.0;
	   pateledr03HcalDepth2TowerSumEt[ipat]      = -999.0;
	   pateledr03HcalTowerSumEt[ipat]            = -999.0;
	   pateledr04EcalRecHitSumEt[ipat]           = -999.0;
	   pateledr04HcalDepth1TowerSumEt[ipat]      = -999.0;
	   pateledr04HcalDepth2TowerSumEt[ipat]      = -999.0;
	   pateledr04HcalTowerSumEt[ipat]            = -999.0;
	   patelee1x5[ipat]                          = -999.0;
	   patelee2x5Max[ipat]                       = -999.0;
	   patelee5x5[ipat]                          = -999.0;
	   pateleeEleClusterOverPout[ipat]           = -999.0;
	   pateleeSeedClusterOverP[ipat]             = -999.0;
	   pateleeSeedClusterOverPout[ipat]          = -999.0;
	   pateleeSuperClusterOverP[ipat]            = -999.0;
	   pateleeta[ipat]                           = -999.0;
	   patelehadronicOverEm[ipat]                = -999.0;
	   patelesigmaIetaIeta[ipat]                 = -999.0;
	   pateledeltaEtaEleClusterTrackAtCalo[ipat] = -999.0;
	   pateledeltaEtaSeedClusterTrackAtCalo[ipat]= -999.0;
	   pateledeltaEtaSuperClusterTrackAtVtx[ipat]= -999.0;
	   pateledeltaPhiEleClusterTrackAtCalo[ipat] = -999.0;
	   pateledeltaPhiSeedClusterTrackAtCalo[ipat]= -999.0;
	   pateledeltaPhiSuperClusterTrackAtVtx[ipat]= -999.0;
	   pateleenergy[ipat]                        = -999.0;
	   patelemva[ipat]                           = -999.0;
	   patelenumberOfTracks[ipat]                = -999;
	   patelemaxEnergyXtal[ipat]                 = -999.;
	   pateleswissCross[ipat]                    = -999.;
	   pateleswissBasedspikevar[ipat]            = -999.;
	   
	   
	   //gsf
	   pateletrkpt[ipat]                         = -999.0;
	   pateletrkcharge[ipat]                     = -999;
	   pateletrkchi2[ipat]                       = -999.0;
	   pateletrketa[ipat]                        = -999.0;
	   pateletrknumberOfLostHits[ipat]           = -999;
	   pateletrknumberOfValidHits[ipat]          = -999;
	   pateletrklost[ipat]                       = -999;
	   pateletrkd0[ipat]                         = -999.0;
	   pateletrkdxy[ipat]                        = -999.0;
	   pateletrkdz[ipat]                         = -999.0;
	   pateletrkptin[ipat]                       = -999.0;
	   pateletrkptout[ipat]                      = -999.0;
	   pateletrkfbrem[ipat]                      = -999.0;
	   pateletrkqoverp[ipat]                     = -999.0;
	   pateletrkvx[ipat]                         = -999.0;
	   pateletrkvy[ipat]                         = -999.0;
	   pateletrkvz[ipat]                         = -999.0;
	   pateletrkphi[ipat]                        = -999.0;
	   pateletrkndof[ipat]                       = -999.0;
	   pateletrkrecHitsSize[ipat]                = -999;
	   pateletrktheta[ipat]                      = -999.0;
	   pateletrkqualityMask[ipat]                = -999;
	   pateletrkouterX[ipat]                     = -999.0;
	   pateletrkouterY[ipat]                     = -999.0;
	   pateletrkouterZ[ipat]                     = -999.0;
	   pateletrkouterRadius[ipat]                = -999.0;
	   pateletrkinnerX[ipat]                     = -999.0;
	   pateletrkinnerY[ipat]                     = -999.0;
	   pateletrkinnerZ[ipat]                     = -999.0;
	   pateleExpectednumberOfHits[ipat]          = -999;
	   pateleconvDist[ipat]                       = -999.;
	   pateleconvDcot[ipat]                       = -999.;
	   pateleconvRadius[ipat]                   = -999;

	   ///1st aug, 2013
	   patelenumberOfLostHits[ipat]          = -999;
	   
	   if(debugEle_)
	     {
	       cout<<"filled all the electron vars"<<endl;
	     }
	 }

       patelsc->Clear();
       patelp4->Clear();
       patelmomvtx->Clear();
       patelmomvtxconst->Clear();
       patelmomcalo->Clear();
       patelposvtx->Clear();
       patelposcalo->Clear();
       
       
       int ipatel = 0;
       for(edm::View<pat::Electron>::const_iterator patele = electrons.begin(); patele!=electrons.end(); ++patele)
	 {
	   if(ipatel >= MAX_ELECTRONS)
	     {
	       cout<<"WARNING: PatElectrons: event has "<<electrons.size()<<" electrons while array can store only "<<MAX_ELECTRONS<<endl;
	       break;
	     }
	   
	   if(debugEle_)
	     {
	       cout<<"Inside pat electrons loop"<<endl;
	     }

	   new ((*patelp4)[ipatel]) TLorentzVector();
	   ((TLorentzVector *)patelp4->At(ipatel))->SetXYZT(patele->px(), patele->py(), patele->pz(), patele->energy());

	   
	   patelefbrem[ipatel]                         = patele->fbrem();
	   pateletrackerDrivenSeed[ipatel]             = patele->trackerDrivenSeed();
	   pateleecalDrivenSeed[ipatel]                = patele->ecalDrivenSeed();
	   patelecharge[ipatel]                        = patele->charge();
	   pateledr03EcalRecHitSumEt[ipatel]           = patele->dr03EcalRecHitSumEt();
	   pateledr03HcalDepth1TowerSumEt[ipatel]      = patele->dr03HcalDepth1TowerSumEt();
	   pateledr03HcalDepth2TowerSumEt[ipatel]      = patele->dr03HcalDepth2TowerSumEt();
	   pateledr03HcalTowerSumEt[ipatel]            = patele->dr03HcalTowerSumEt();
	   pateledr04EcalRecHitSumEt[ipatel]           = patele->dr04EcalRecHitSumEt();
	   pateledr04HcalDepth1TowerSumEt[ipatel]      = patele->dr04HcalDepth1TowerSumEt();
	   pateledr04HcalDepth2TowerSumEt[ipatel]      = patele->dr04HcalDepth2TowerSumEt();
	   pateledr04HcalTowerSumEt[ipatel]            = patele->dr04HcalTowerSumEt();
	   patelee1x5[ipatel]                          = patele->e1x5();
	   patelee2x5Max[ipatel]                       = patele->e2x5Max();
	   patelee5x5[ipatel]                          = patele->e5x5();
	   pateleeEleClusterOverPout[ipatel]           = patele->eEleClusterOverPout();
	   pateleeSeedClusterOverP[ipatel]             = patele->eSeedClusterOverP();
	   pateleeSeedClusterOverPout[ipatel]          = patele->eSeedClusterOverPout();
	   pateleeSuperClusterOverP[ipatel]            = patele->eSuperClusterOverP();
	   patelehadronicOverEm[ipatel]                = patele->hadronicOverEm();

	   if(debugEle_){
	     cout<<"HadOe:"<<patelehadronicOverEm[ipatel]<<"patele->hadronicOverEm():"<<patele->hadronicOverEm()<<endl;
	   }

	   patelesigmaIetaIeta[ipatel]                 = patele->sigmaIetaIeta();
	   pateledr03TkSumPt[ipatel]                   = patele->dr03TkSumPt();
	   pateledr04TkSumPt[ipatel]                   = patele->dr04TkSumPt();
	   
	   pateledeltaEtaEleClusterTrackAtCalo[ipatel] = patele->deltaEtaEleClusterTrackAtCalo();
	   pateledeltaEtaSeedClusterTrackAtCalo[ipatel]= patele->deltaEtaSeedClusterTrackAtCalo();
	   pateledeltaEtaSuperClusterTrackAtVtx[ipatel]= patele->deltaEtaSuperClusterTrackAtVtx();
	   pateledeltaPhiEleClusterTrackAtCalo[ipatel] = patele->deltaPhiEleClusterTrackAtCalo();
	   pateledeltaPhiSeedClusterTrackAtCalo[ipatel]= patele->deltaPhiSeedClusterTrackAtCalo();
	   pateledeltaPhiSuperClusterTrackAtVtx[ipatel]= patele->deltaPhiSuperClusterTrackAtVtx();
	   patelemva[ipatel]                           = patele->mva();
	   patelenumberOfTracks[ipatel]                = patele->numberOfTracks();
	   
	   /////added on 20th sept, 2011
	   patelecaloEnergy[ipatel]                    = patele->caloEnergy();
	   pateleecalEnergy[ipatel]                    = patele->ecalEnergy();
	   pateleecalEnergyError[ipatel]               = patele->ecalEnergyError();
	   pateletrackMomentumError[ipatel]            = patele->trackMomentumError();
	   
	   ////added on 27th may, 2012 
	   //patelecaloCorrectedEnergy[ipatel]          = patele->caloCorrectedEnergy();
	   patelecorrectedEcalEnergy[ipatel]            = patele->correctedEcalEnergy();
	   
	   if(debugEle_){
	     cout<<"electron corrected energy"<<patelecorrectedEcalEnergy[ipatel]<<endl;
	   }

	   
	   new ((*patelmomvtxconst)[ipatel]) TVector3();
	   ((TVector3 *)patelmomvtxconst->At(ipatel))->SetXYZ(patele->trackMomentumAtVtxWithConstraint().x(),
							     patele->trackMomentumAtVtxWithConstraint().y(), patele->trackMomentumAtVtxWithConstraint().z());
	   
	   new ((*patelmomvtx)[ipatel]) TVector3();
	   ((TVector3 *)patelmomvtx->At(ipatel))->SetXYZ(patele->trackMomentumAtVtx().x(),
							patele->trackMomentumAtVtx().y(), patele->trackMomentumAtVtx().z());
	   
	   new ((*patelmomcalo)[ipatel]) TVector3();
	   ((TVector3 *)patelmomcalo->At(ipatel))->SetXYZ(patele->trackMomentumAtCalo().x(),
							 patele->trackMomentumAtCalo().y(), patele->trackMomentumAtCalo().z());
	   
	   new ((*patelmomout)[ipatel]) TVector3();
	   ((TVector3 *)patelmomout->At(ipatel))->SetXYZ(patele->trackMomentumOut().x(),
							patele->trackMomentumOut().y(), patele->trackMomentumOut().z());
	   
	   new ((*patelposvtx)[ipatel]) TVector3();
	   ((TVector3 *)patelposvtx->At(ipatel))->SetXYZ(patele->trackPositionAtVtx().x(),
							patele->trackPositionAtVtx().y(), patele->trackPositionAtVtx().z());
	   
	   new ((*patelposcalo)[ipatel]) TVector3();
	   ((TVector3 *)patelposcalo->At(ipatel))->SetXYZ(patele->trackPositionAtCalo().x(),
							  patele->trackPositionAtCalo().y(), patele->trackPositionAtCalo().z());
	   
	   
	   
	   patelepout[ipatel] = patele->trackMomentumOut().R();
	   patelepin[ipatel] = patele->trackMomentumAtVtx().R();
	   
	   
	   /*//my heep id
	   pateleheepcutword[ipatel]                   = heepeleID( patelectron_container, ipatel );       //for setting heepcutbit
	   int cutmask = ~0x0;
	   int cutword = (cutmask^pateleheepcutword[ipatel]);
	   //std::cout<<"cutword = "<<cutword<<std::endl;
	       if( cutword == 0x0 )
                 pateleheepid[ipatel] = 1; 
               else
                 pateleheepid[ipatel] = 0; 

	       //std::cout<<"heepcutbit  for this electron = "<<pateleheepcutword[ipatel]<<std::endl;
	       //std::cout<<"heep id for this ele = "<<pateleheepid[ipatel]<<std::endl;

	       //std::cout<<"pateletrackerDrivenSeed[ipatel] = "<<pateletrackerDrivenSeed[ipatel]<<std::endl;
	       */
	       //SHarpers heep id
	   /*if( cuts_.passCuts(*patele) )
	     //std::cout<<"passed SHarpers heep id"<<std::endl;
	     samspateleheepid[ipatel] = 1;
	   
	   if( !cuts_.passCuts(*patele) )
	     //std::cout<<"could not passed SHarpers heep id"<<std::endl;
	     samspateleheepid[ipatel] = 0;
	   */
	   
	       //gaps
	       pateleisEB[ipatel]                                  = patele->isEB();
	       pateleisEBEEGap[ipatel]                             = patele->isEBEEGap();
	       pateleisEBEtaGap[ipatel]                            = patele->isEBEtaGap();
	       pateleisEBGap[ipatel]                               = patele->isEBGap();
	       pateleisEBPhiGap[ipatel]                            = patele->isEBPhiGap();
	       pateleisEE[ipatel]                                  = patele->isEE();
	       pateleisEEDeeGap[ipatel]                            = patele->isEEDeeGap();
	       pateleisEEGap[ipatel]                               = patele->isEEGap();
	       pateleisEERingGap[ipatel]                           = patele->isEERingGap();

	       //sceta
	       reco::SuperClusterRef scref = patele->superCluster();

	       
	       if( !(patele->superCluster().isNonnull()) )
		 std::cout<<"scref of patelectrons is not valid"<<std::endl;
	       
	       double sceta                         = scref->eta();
	       double scphi                         = scref->phi();
	       double sce                           = scref->energy();

	       double sctheta = (2*atan(exp(-sceta)));
	       double scpx = sce*sin(sctheta)*cos(scphi);
	       double scpy = sce*sin(sctheta)*sin(scphi);
	       double scpz = sce*cos(sctheta);

	       new ((*patelsc)[ipatel]) TLorentzVector();
	       ((TLorentzVector *)patelsc->At(ipatel))->SetXYZT(scpx, scpy, scpz, sce);
	       

	       
	       //maximum energy xtal
	       const reco::CaloClusterPtr seed = scref->seed();
	       patelemaxEnergyXtal[ipatel]                 = lazyTool.eMax( *seed );

	       //spike info
	       if(patele->isEB())
		 {
		   pateleswissCross[ipatel]   = EcalClusterTools::eTop( *(patele->superCluster()->seed()), &(*barrelRecHits), &(*topology))+ EcalClusterTools::eBottom( *(patele->superCluster()->seed()), &(*barrelRecHits), &(*topology)) + EcalClusterTools::eLeft( *(patele->superCluster()->seed()), &(*barrelRecHits), &(*topology)) + EcalClusterTools::eRight( *(patele->superCluster()->seed()), &(*barrelRecHits), &(*topology));
		   //std::cout<<"etop: "<<EcalClusterTools::eTop( *(patele->superCluster()->seed()), &(*barrelRecHits), &(*topology))<<endl;
		   //std::cout<<"etop: "<<EcalClusterTools::eTop( *(patele->superCluster()->seed()), &(*barrelRecHits), &(*topology))<<endl;
		   //if(1-pateleswissCross[ipatel]/patelemaxEnergyXtal[ipatel] > 0.95) cout<<"This electron candidate is an ECAL spike identified by Swiss Cross algorithm."<<endl;
		   
		   pateleswissBasedspikevar[ipatel]             = 1-pateleswissCross[ipatel]/patelemaxEnergyXtal[ipatel];
		 }//end of if(patele->isEB())
	       else{
		 pateleswissCross[ipatel]   = EcalClusterTools::eTop( *(patele->superCluster()->seed()), &(*endcapRecHits), &(*topology))+ EcalClusterTools::eBottom( *(patele->superCluster()->seed()), &(*endcapRecHits), &(*topology)) + EcalClusterTools::eLeft( *(patele->superCluster()->seed()), &(*endcapRecHits), &(*topology)) + EcalClusterTools::eRight( *(patele->superCluster()->seed()), &(*endcapRecHits), &(*topology));
		 //if(1-pateleswissCross[ipatel]/patelemaxEnergyXtal[ipatel] > 0.95) cout<<"This electron candidate is an ECAL spike identified by Swiss Cross algorithm."<<endl;
		 pateleswissBasedspikevar[ipatel]             = 1-pateleswissCross[ipatel]/patelemaxEnergyXtal[ipatel];
	       }//end of else

	       
	       //timing info
	       DetId pateleSeedDetId = lazyTool.getMaximum(*seed).first;

	       if ( patele->isEB() ) {
		 EcalRecHitCollection::const_iterator eleebrhit = barrelRecHits->find(pateleSeedDetId);
		 if ( eleebrhit != barrelRecHits->end() ) { 
		   pateleseedtime[ipatel]         = eleebrhit->time(); 
		   patelerecoFlag[ipatel]         = eleebrhit->recoFlag();
		   //patelecheckFlag[ipatel]         = eleebrhit->checkFlag();
		   //patelekOutOfTime[ipatel]       = eleebrhit->kOutOfTime();
		   //pateleseverityLevel[ipatel]    = EcalSeverityLevelAlgo::severityLevel( pateleSeedDetId, (*barrelRecHits), *chStatus );
		   patelee2e9[ipatel]             = mye2overe9(pateleSeedDetId,*barrelRecHits);
		 }//if ( eleebrhit != barrelRecHits->end() ) 
	       }//if ( patele->isEB() && barrelRecHits.isValid() )
	       else {
		 EcalRecHitCollection::const_iterator eleeerhit = endcapRecHits->find(pateleSeedDetId);
		 if ( eleeerhit != endcapRecHits->end() ) { 
		   pateleseedtime[ipatel]         = eleeerhit->time(); 
		   patelerecoFlag[ipatel]         = eleeerhit->recoFlag();
		   //patelecheckFlag[ipatel]         = eleeerhit->checkFlag();
		   //patelekOutOfTime[ipatel]       = eleeerhit->kOutOfTime();
		   //pateleseverityLevel[ipatel]    = EcalSeverityLevelAlgo::severityLevel( pateleSeedDetId, (*endcapRecHits), *chStatus );
		   patelee2e9[ipatel]             = 0;
		 }//if ( eleeerhit != endcapRecHits->end() )
	       }//else if ( endcapRecHits.isValid() )

	       int index = 0;
	       patelescind[ipatel] = -1;
	       reco::GsfElectron localele = reco::GsfElectron(*patele);
	       // loop over the two SC collections
	       for(int z = 0; z<2; ++z) {
		 if (z == 0) {
		   for(reco::SuperClusterCollection::size_type j = 0; j<superClustersHybridH->size(); ++j) {
		     
		     reco::SuperClusterRef cluster(superClustersHybridH, j);
		     //if (&(*(patele->superCluster())) == &(*cluster)) {
		     if(&(*localele.superCluster()) == &(*cluster)){
			 
		       patelescind[ipatel] = index;
		     }
		     index++;
		   }//for(reco::SuperClusterCollection::size_type j = 0; j<superClustersHybridH->size(); ++j)
		 }//if (z == 0) 

		 
		 if (z == 1) {
		   for(reco::SuperClusterCollection::size_type j = 0; j<superClustersEndcapH->size(); ++j) {
		     
		     reco::SuperClusterRef cluster(superClustersEndcapH, j);
		     //if (&(*(patele->superCluster())) == &(*cluster)) {
		     if(&(*localele.superCluster()) == &(*cluster)){
		       
		       patelescind[ipatel] = index;
		     }//if (&(*(patele->superCluster())) == &(*cluster))
		     index++;
		   }//for(reco::SuperClusterCollection::size_type j = 0; j<superClustersEndcapH->size(); ++j)
		 }//if (z == 1)
	       }//for(int z = 0; z<2; ++z)


	       //PF cluster
	       reco::SuperClusterRef pfref = patele->pflowSuperCluster(); 
	       if(patele->pflowSuperCluster().isNonnull()) 
		 {
		   patelepfeta[ipatel]                         = pfref->eta();
		   patelepfphi[ipatel]                         = pfref->phi();
		   patelepfE[ipatel]                           = pfref->energy();
		   //patelepfpt[ipatel]                          = pfref->pt();
		 }
	       
	       
	       //trackinfo
	       reco::GsfTrackRef trackref = patele->gsfTrack();
	       
	       if(patele->gsfTrack().isNonnull())
		 {
		   pateletrkpt[ipatel]                         = trackref->pt();
		   pateletrkcharge[ipatel]                     = trackref->charge();
		   pateletrkchi2[ipatel]                       = trackref->chi2();
		   pateletrketa[ipatel]                        = trackref->eta();
		   pateletrknumberOfLostHits[ipatel]           = trackref->numberOfLostHits();
		   pateletrknumberOfValidHits[ipatel]          = trackref->numberOfValidHits();
		   pateletrklost[ipatel]                       = trackref->lost();

		   pateletrkd0[ipatel]                         = trackref->d0();
		   
		   if( vtx->size() > 0 ){
		     pateletrkdxy[ipatel]                        = trackref->dxy(vtx->front().position());
		     pateletrkdz[ipatel]                         = trackref->dz(vtx->front().position());
		   } else{
		     pateletrkdxy[ipatel]                        = trackref->dxy();
		     pateletrkdz[ipatel]                         = trackref->dz();
		   }

		   /*pateletrkptin[ipatel]                       = sqrt(trackref->innerMomentum().Perp2());
		   pateletrkptout[ipatel]                      = sqrt(trackref->outerMomentum().Perp2());
		   pateletrkfbrem[ipatel]                      = (pateletrkptin[ipatel]-pateletrkptout[ipatel])/pateletrkptin[ipatel];
		   */
		   pateletrkqoverp[ipatel]                     = trackref->qoverp();
		   pateletrkvx[ipatel]                         = trackref->vx();
		   pateletrkvy[ipatel]                         = trackref->vy();
		   pateletrkvz[ipatel]                         = trackref->vz();
		   pateletrkphi[ipatel]                        = trackref->phi();
		   pateletrkndof[ipatel]                       = trackref->ndof();
		   //pateletrkrecHitsSize[ipatel]                = trackref->recHitsSize();
		   pateletrktheta[ipatel]                      = trackref->theta();
		   pateletrkqualityMask[ipatel]                = trackref->qualityMask();
		   /*pateletrkouterX[ipatel]                     = trackref->outerX();
		   pateletrkouterY[ipatel]                     = trackref->outerY();
		   pateletrkouterZ[ipatel]                     = trackref->outerZ();
		   pateletrkouterRadius[ipatel]                = trackref->outerRadius();
		   pateletrkinnerX[ipatel]                     = trackref->innerPosition().X();
		   pateletrkinnerY[ipatel]                     = trackref->innerPosition().Y();
		   pateletrkinnerZ[ipatel]                     = trackref->innerPosition().Z();*/

		   /////////////////////////////////////////////////////////conversion rejection
		   //1. extpected hits
		   const reco::HitPattern& p_inner = trackref->trackerExpectedHitsInner(); 
		   pateleExpectednumberOfHits[ipatel] = p_inner.numberOfHits();

		   ///1st Aug
		   patelenumberOfLostHits[ipatel] = p_inner.numberOfLostHits();
		   
		   //std::cout<<"NumberOfExpectedInnerHits : "<<pateleExpectednumberOfHits[ipatel]<<std::endl;

		   ///2. dcot(theta) && Dist bet tracks
		   pateleconvDist[ipatel] = patele->convDist();
		   pateleconvDcot[ipatel] = patele->convDcot();
		   //pateleconvDcot[ipatel] = patele->convDcot();
		   pateleconvRadius[ipatel] = patele->convRadius();
		   
		   ///////////////////////////////////////////////////////////////////////////////////////////////////
		   
		 }//if(patele->gsfTrack().isNonnull())
	       
	       ipatel++;
	 }//for(edm::View<pat::Electron>::const_iterator patele = electrons.begin(); patele!=electrons.end(); ++patele)
       patelesize = ipatel;
     }//if(runPatelectrons_)

   //pat photons
   if(runPatphotons_)
     {
       edm::Handle<edm::View<pat::Photon> > photonHandle;
       try{iEvent.getByLabel(phoLabel_,photonHandle); 
         //std::cout<<"Found the handle"<<std::endl;                                                                                                           
       }
       catch(...){std::cout<<"Didn't fine the pat::photon handle"<<std::endl;}
       const edm::View<pat::Photon> & photons = *photonHandle;
       
       //////19th July
       ////////////////////////////////for PF photons ////////////////////////////////////////
       edm::Handle<reco::PhotonCollection> recophoH;
       bool found = iEvent.getByLabel(recophoLabel_,recophoH);
       if(!found ) {
	 std::ostringstream  err;
	 err<<" cannot get reco Photons: "
	    <<recophoLabel_<<std::endl;
	 edm::LogError("Analyzer")<<err.str();
	 throw cms::Exception( "MissingProduct", err.str());
       }

       // get the iso deposits. 3 (charged hadrons, photons, neutral hadrons)
       unsigned nTypes=3;
       
       typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > > IsoDepositMaps;
       typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;


       IsoDepositMaps photonIsoDep(nTypes);
       for (size_t j = 0; j<inputTagIsoDepPhotons_.size(); ++j) {
	 iEvent.getByLabel(inputTagIsoDepPhotons_[j], photonIsoDep[j]);
       }

       IsoDepositVals photonIsoValPFId(nTypes);
       
       for (size_t j = 0; j<inputTagIsoValPhotonsPFId_.size(); ++j) {
	 iEvent.getByLabel(inputTagIsoValPhotonsPFId_[j], photonIsoValPFId[j]);
       }
       
       const IsoDepositVals * photonIsoVals = &photonIsoValPFId;
       /////fill isolation below


       //////updated on 11th April, 2013
       // PF Candidates
       Handle<reco::PFCandidateCollection> pfCandidatesHandle;
       iEvent.getByLabel(PFCandLabel_, pfCandidatesHandle);
       const PFCandidateCollection thePfColl = *(pfCandidatesHandle.product());

       
       // Photon
       edm::Handle<reco::VertexCollection> recVtxs;
       iEvent.getByLabel(vertexLabel_, recVtxs);

       PFIsolationEstimator isolator;
       isolator.initializePhotonIsolation(kTRUE);
       isolator.setConeSize(0.3);

       ///for checks with Yong
       /*PFIsolationEstimator isolator04;
       isolator04.initializePhotonIsolation(kTRUE);
       isolator04.setConeSize(0.4);
       */

       
       VertexRef myVtxRef(recVtxs, 0);
       /////////updated on 11th april, 2013


       //////////////////////////////////////END of PF iso///////////////////////////////////////////////

       ////conversion safe electron veto for photon ID
       edm::Handle<reco::BeamSpot> bsHandle;
       iEvent.getByLabel("offlineBeamSpot", bsHandle);
       const reco::BeamSpot &beamspot = *bsHandle.product();

       edm::Handle<reco::ConversionCollection> hConversions;
       iEvent.getByLabel("allConversions", hConversions);

       edm::Handle<reco::GsfElectronCollection> hElectrons;
       iEvent.getByLabel("gsfElectrons", hElectrons);


       



       for(int ipat=0; ipat<MAX_PHOTONS; ipat++)
	 {
	   patphoecalRecHitSumEtConeDR03[ipat]           = -999.;
	   patphohcalDepth1TowerSumEtConeDR03[ipat]      = -999.;
	   patphohcalDepth2TowerSumEtConeDR03[ipat]      = -999.;
	   patphohcalTowerSumEtConeDR03[ipat]            = -999.;
	   patphotrkSumPtHollowConeDR03[ipat]            = -999.;
	   patphotrkSumPtSolidConeDR03[ipat]             = -999. ;
	   patphonTrkHollowConeDR03[ipat]                = -999.;
	   patphonTrkSolidConeDR03[ipat]                 = -999.;
	   
	   patphoecalRecHitSumEtConeDR04[ipat]           = -999.;
	   patphohcalDepth1TowerSumEtConeDR04[ipat]      = -999.;
	   patphohcalDepth2TowerSumEtConeDR04[ipat]      = -999.;
	   patphohcalTowerSumEtConeDR04[ipat]            = -999.;
	   patphotrkSumPtHollowConeDR04[ipat]            = -999.;
	   patphotrkSumPtSolidConeDR04[ipat]             = -999.;
	   patphonTrkHollowConeDR04[ipat]                = -999.;
	   patphonTrkSolidConeDR04[ipat]                 = -999.;
	   
	   patphoe1x5[ipat]                              = -999.;
	   patphoe2x5[ipat]                              = -999.;
	   patphoe3x3[ipat]                              = -999.;
	   patphoe5x5[ipat]                              = -999.;
	   patphoeta[ipat]                               = -999.;
	   patphohadronicOverEm[ipat]                    = -999.;
	   patphosigmaIetaIeta[ipat]                     = -999.;
	   patphor1x5[ipat]                              = -999.;
	   patphor2x5[ipat]                              = -999.;
	   patphor9[ipat]                                = -999.;
	   patphonumberOfTracks[ipat]                    = -999.;
	   patphohasPixelSeed[ipat]                      = -999;
	   patphoisConvertedPhoton[ipat]                 = -999.;
	   patphomaxEnergyXtal[ipat]                     = -999.;

	   patphosigmaIetaIphi[ipat]                     = -999.;
	   patphosigmaIphiIphi[ipat]                     = -999.;
	   
	   //swiss cross 
	   patphoswissCross[ipat]                        = -999.;
	   patphoswissBasedspikevar[ipat]                = -999.;
	   
	   ////mip tagger
	   patphomipChi2[ipat]                           = -999.;
	   patphomipIntercept[ipat]                      = -999.;
	   patphomipIsHalo[ipat]                         = -999.;
	   patphomipNhitCone[ipat]                       = -999.;
	   patphomipSlope[ipat]                          = -999.;
	   patphomipTotEnergy[ipat]                      = -999.;

	   //conversion  
	   patphoconvsize[ipat]                          = -999.;
	   for(int iconv=0; iconv<10; iconv++)
	     {
	       patphoconvtxX[ipat][iconv]                = -999.;
	       patphoconvtxY[ipat][iconv]                = -999.;
	       patphoconvtxZ[ipat][iconv]                = -999.;
	       patphoconvtxR[ipat][iconv]                = -999.;
	     }
	   patphotrkpt[ipat]                             = -999.;
	   patphotrkcharge[ipat]                         = -999.;
	   patphotrkchi2[ipat]                           = -999.;
	   patphotrketa[ipat]                            = -999.;
	   patphotrknumberOfLostHits[ipat]               = -999.;
	   patphotrknumberOfValidHits[ipat]              = -999.;
	   patphotrklost[ipat]                           = -999.;
	   patphotrkd0[ipat]                             = -999.;
	   patphotrkdxy[ipat]                            = -999.;
	   patphotrkdz[ipat]                             = -999.;
	   patphotrkptin[ipat]                           = -999.;
	   patphotrkptout[ipat]                          = -999.;
	   patphotrkfbrem[ipat]                          = -999.;
	   patphotrkqoverp[ipat]                         = -999.;
	   patphotrkvx[ipat]                             = -999.;
	   patphotrkvy[ipat]                             = -999.;
	   patphotrkvz[ipat]                             = -999.;
	   patphotrkphi[ipat]                            = -999.;
	   patphotrkndof[ipat]                           = -999.;
	   patphotrkrecHitsSize[ipat]                    = -999.;
	   patphotrktheta[ipat]                          = -999.;
	   patphotrkqualityMask[ipat]                    = -999.;
	   patphotrkouterX[ipat]                         = -999.;
	   patphotrkouterY[ipat]                         = -999.;
	   patphotrkouterZ[ipat]                         = -999.;
	   patphotrkouterRadius[ipat]                    = -999.;
	   patphotrkinnerX[ipat]                         = -999.;
	   patphotrkinnerY[ipat]                         = -999.;
	   patphotrkinnerZ[ipat]                         = -999.;
	  
	   patphomGenpdgId[ipat]                         = -999.;                     
	   patphomGentheta[ipat]                         = -999.;
	   patphomGenphi[ipat]                           = -999.; 
	   patphomGenpt[ipat]                            = -999.;
	   patphomGenpx[ipat]                            = -999.;
	   patphomGenpy[ipat]                            = -999.;
	   patphomGenpz[ipat]                            = -999.;
	   patphomGenenergy[ipat]                        = -999.;
	   patphomGenisJet[ipat]                         = -999.;
	   patphomGenisPhoton[ipat]                      = -999.; 
	   
	   for(int imoth=0; imoth<10; imoth++)
	     patphomGenmompdgId[ipat][imoth]             = -999.;

	   patphomGengranpdgId[ipat]                     = -999.; 
	   patphomGenisElectron[ipat]                    = -999.;

	   //PF iso - 19th July
	   patphopfchargediso[ipat]                      = -999.;
	   patphopfphotoniso[ipat]                       = -999.;
	   patphopfneutraliso[ipat]                      = -999.;

	   patphonewpfchargediso[ipat]                      = -999.;
	   patphonewpfphotoniso[ipat]                       = -999.;
	   patphonewpfneutraliso[ipat]                      = -999.;

	   
	 }
       
       patphosize = -999;
       
       

       
       std::vector<myrecophoclass> myrecopho;
       myrecophoclass myrecophoobj;

       
       unsigned npho=recophoH->size();
    
       for(unsigned ipho=0; ipho<npho;++ipho) {
	 reco::PhotonRef myPhotonRef(recophoH,ipho);
      
	 /*myrecophoobj.charged = ((*(*photonIsoVals)[0])[myPhotonRef]/myPhotonRef->pt());
	   myrecophoobj.photon  = ((*(*photonIsoVals)[1])[myPhotonRef]/myPhotonRef->pt());
	 myrecophoobj.neutral = ((*(*photonIsoVals)[2])[myPhotonRef]/myPhotonRef->pt());
	 */
	 
	 myrecophoobj.charged = ((*(*photonIsoVals)[0])[myPhotonRef]);
	 myrecophoobj.photon  = ((*(*photonIsoVals)[1])[myPhotonRef]);
	 myrecophoobj.neutral = ((*(*photonIsoVals)[2])[myPhotonRef]);
	 

	 myrecophoobj.pt = myPhotonRef->pt();
	 
	 ////////////11th april, 2013
	 // New PF isolation
	 //edm::Ptr<reco::Candidate> recoPhoRef = myPhotonRef.originalObjectRef();
	 //const reco::Photon *recoPhoton = dynamic_cast<const reco::Photon *>(recoPhoRef.get());

	 isolator.fGetIsolation(&*myPhotonRef, &thePfColl, myVtxRef, recVtxs);
	 myrecophoobj.recomencharged      = isolator.getIsolationCharged();
	 myrecophoobj.recomenphoton       = isolator.getIsolationPhoton();
	 myrecophoobj.recomenneutral      = isolator.getIsolationNeutral();


	 /////for checks with Yong
	 /*isolator04.fGetIsolation(&*myPhotonRef, &thePfColl, myVtxRef, recVtxs);
	 cout<<"CHARGED:  Cone of 0.3: CONE OF 0.4: " <<isolator.getIsolationCharged()<<":"<<isolator04.getIsolationCharged()<<endl;
	 cout<<"PHOTON :  Cone of 0.3: CONE OF 0.4: " <<isolator.getIsolationPhoton()<<":"<<isolator04.getIsolationPhoton()<<endl;
	 cout<<"NEUTRAL:  Cone of 0.3: CONE OF 0.4: " <<isolator.getIsolationNeutral()<<":"<<isolator04.getIsolationNeutral()<<endl;
	 */

	 ////////////11th april, 2013


	 /////////////21st june, 2013
	 ////////WORST VERTEX
	 //-------------------------------------------------------                                                                            
	 //Get VtxIso for the worst sum of charged hadron isolation                                                                           
	 //-------------------------------------------------------                                                                            

	 //loops over all vertices, takes maximum                                                                                              
	 myrecophoobj.PFphotonWorstChargedHadronIso=0;

	 for(int iv=0;iv<int(recVtxs->size());++iv)
	   {
	     reco::VertexRef thisVtxRef(recVtxs, iv);
	     //isolator03_.fGetIsolation( &*myPhotonRef , pfH.product(),thisVtxRef, VtxsH);
	     isolator.fGetIsolation(&*myPhotonRef, &thePfColl, thisVtxRef, recVtxs);
	     Float_t thisChargedHadronIso = isolator.getIsolationCharged();

	     if(thisChargedHadronIso > myrecophoobj.PFphotonWorstChargedHadronIso) 
	       myrecophoobj.PFphotonWorstChargedHadronIso = thisChargedHadronIso;
	   }//loop over vtx                                                                                                                


	 ///////WORST VERTEX
	 
	 if(debugPho_)
	   {
	     cout<<"=====PF info======="<<endl;
	     std::cout << "Photon: " << " run " << iEvent.id().run() << " lumi " << iEvent.id().luminosityBlock() << " event " << iEvent.id().event();
	     std::cout << " pt " <<  myPhotonRef->pt() << " eta " << myPhotonRef->eta() << " phi " << myPhotonRef->phi() << " charge " << myPhotonRef->charge()<< " : ";
	     std::cout << " ChargedIso " << myrecophoobj.charged ;
	     std::cout << " PhotonIso " <<  myrecophoobj.photon ;
	     std::cout << " NeutralHadron Iso " << myrecophoobj.neutral << std::endl;
	   }

	 if(debugphoPFIso_)
	   {
	     cout<<"=====PF info======="<<endl;
	     std::cout << "Photon: " << " run " << iEvent.id().run() << " lumi " << iEvent.id().luminosityBlock() << " event " << iEvent.id().event();
	     std::cout << " pt " <<  myPhotonRef->pt() << " eta " << myPhotonRef->eta() << " phi " << myPhotonRef->phi() << " charge " << myPhotonRef->charge()<< " : ";
	     std::cout << " ChargedIso : old : recommend: " << myrecophoobj.charged <<" : "<<myrecophoobj.recomencharged << std::endl;
	     std::cout << " PhotonIso : old : recommend: " <<  myrecophoobj.photon <<" : "<<myrecophoobj.recomenphoton << std::endl;
	     std::cout << " NeutralHadron Iso : old : recommend: " << myrecophoobj.neutral << " : "<<myrecophoobj.recomenneutral << std::endl;

	     if( myrecophoobj.charged != myrecophoobj.recomencharged ) {
	       
	       std::cout<<"WARNINF from debugphoPFIso_:old charged not the same as recom charged" << std::endl;
	     }
	     
	     if( myrecophoobj.photon != myrecophoobj.recomenphoton ) {
	       std::cout<<"WARNINF from debugphoPFIso_:old photon not the same as recom  photon "<< std::endl;
	     }

	     if( myrecophoobj.neutral != myrecophoobj.recomenneutral ) {
	       std::cout<<"WARNINF from debugphoPFIso_:old neutral not the same as recom neutral " << std::endl;
	     }

	   }
	 
	 ///////////conversion safe electron veto///////////////////////
	 bool issafeveto = !ConversionTools::hasMatchedPromptElectron(myPhotonRef->superCluster(), hElectrons, hConversions, beamspot.position());
	 
	 if(issafeveto)
	   myrecophoobj.passelectronveto = 0;

	 if(!issafeveto)
           myrecophoobj.passelectronveto = 1;
	 
	 

	 myrecopho.push_back(myrecophoobj);
       }//for(unsigned ipho=0; ipho<npho;++ipho)

       if(npho!=0)
	 {
	   std::sort(myrecopho.begin(),myrecopho.end(),PtSortCriterium());
	 }
       
       

       //barrel sc's
       edm::Handle<reco::SuperClusterCollection> superClustersHybridH;
       iEvent.getByLabel(hybridSuperClusterColl_,superClustersHybridH);
       
       //EE sc's
       edm::Handle<reco::SuperClusterCollection> superClustersEndcapH;
       iEvent.getByLabel(endcapSuperClusterColl_, superClustersEndcapH);
       
       int ipatpho = 0;
       
       patphop4->Clear();
       patphocalopos->Clear();
       patphoscp4->Clear();

       for(edm::View<pat::Photon>::const_iterator patpho = photons.begin(); patpho!=photons.end(); ++patpho)
	 {
	   if(ipatpho >= MAX_PHOTONS)
	     {
	       cout<<"WARNING: PatPhotons: event has "<<photons.size()<<" photons while array can store only "<<MAX_PHOTONS<<endl;
	       break;
	     }
	   
	   new ((*patphop4)[ipatpho]) TLorentzVector();
	   new ((*patphocalopos)[ipatpho]) TVector3();
	 
	   ((TLorentzVector *)patphop4->At(ipatpho))->SetXYZT(patpho->px(), patpho->py(), patpho->pz(), patpho->energy());
	   ((TVector3 *)patphocalopos->At(ipatpho))->SetXYZ(patpho->caloPosition().x(), patpho->caloPosition().y(), patpho->caloPosition().z());
	    
	    patphoecalRecHitSumEtConeDR03[ipatpho]       = patpho->ecalRecHitSumEtConeDR03();
	    patphohcalDepth1TowerSumEtConeDR03[ipatpho]  = patpho->hcalDepth1TowerSumEtConeDR03();
	    patphohcalDepth2TowerSumEtConeDR03[ipatpho]  = patpho->hcalDepth2TowerSumEtConeDR03();
	    patphohcalTowerSumEtConeDR03[ipatpho]        = patpho->hcalTowerSumEtConeDR03();
	    patphotrkSumPtHollowConeDR03[ipatpho]        = patpho->trkSumPtHollowConeDR03();
	    patphotrkSumPtSolidConeDR03[ipatpho]         = patpho->trkSumPtSolidConeDR03();
	    patphonTrkHollowConeDR03[ipatpho]            = patpho->nTrkHollowConeDR03();
	    patphonTrkSolidConeDR03[ipatpho]             = patpho->nTrkSolidConeDR03();
	    
	    patphoecalRecHitSumEtConeDR04[ipatpho]       = patpho->ecalRecHitSumEtConeDR04();
	    patphohcalDepth1TowerSumEtConeDR04[ipatpho]  = patpho->hcalDepth1TowerSumEtConeDR04();
	    patphohcalDepth2TowerSumEtConeDR04[ipatpho]  = patpho->hcalDepth2TowerSumEtConeDR04();
	    patphohcalTowerSumEtConeDR04[ipatpho]        = patpho->hcalTowerSumEtConeDR04();
	    patphotrkSumPtHollowConeDR04[ipatpho]        = patpho->trkSumPtHollowConeDR04();
	    patphotrkSumPtSolidConeDR04[ipatpho]         = patpho->trkSumPtSolidConeDR04();
	    patphonTrkHollowConeDR04[ipatpho]            = patpho->nTrkHollowConeDR04();
	    patphonTrkSolidConeDR04[ipatpho]             = patpho->nTrkSolidConeDR04();
	    
	    patphoe1x5[ipatpho]                          = patpho->e1x5();
	    patphoe2x5[ipatpho]                          = patpho->e2x5();
	    patphoe3x3[ipatpho]                          = patpho->e3x3();
	    patphoe5x5[ipatpho]                          = patpho->e5x5();
	    patphoeta[ipatpho]                           = patpho->eta();
	    patphohadronicOverEm[ipatpho]                = patpho->hadronicOverEm();
	    
	    //cout<<"patphohadronicOverEm["<<ipatpho<<"] "<<patphohadronicOverEm[ipatpho]<<endl;
	    
	    patphosigmaIetaIeta[ipatpho]                 = patpho->sigmaIetaIeta();
	    
	    patphor1x5[ipatpho]                          = patpho->r1x5();
	    patphor2x5[ipatpho]                          = patpho->r2x5();
	    patphor9[ipatpho]                            = patpho->r9();
	    
	    patphonumberOfTracks[ipatpho]                = patpho->numberOfTracks();
	    //std::cout<<"patphotrackerDrivenSeed[ipatpho] = "<<patphotrackerDrivenSeed[ipatpho]<<std::endl;
	    
	    patphohasPixelSeed[ipatpho]                  = patpho->hasPixelSeed();
	    patphoisConvertedPhoton[ipatpho]             = patpho->isConvertedPhoton();

	    //////////////////////////////////////20th July///////////////////////////////////////////////
	    /////PF isolation info
	    patphopfchargediso[ipatpho]                  = myrecopho[ipatpho].charged;
	    patphopfphotoniso[ipatpho]                   = myrecopho[ipatpho].photon;
	    patphopfneutraliso[ipatpho]                  = myrecopho[ipatpho].neutral;
	    
	    
	    /////worst vertex isolation - 21st june
	    patphoPFphotonWorstChargedHadronIso[ipatpho] = myrecopho[ipatpho].PFphotonWorstChargedHadronIso;
	    ///////11th april
	    /*patphopfchargediso[ipatpho]                  = myrecopho[ipatpho].recomencharged;
	    patphopfphotoniso[ipatpho]                   = myrecopho[ipatpho].recomenphoton;
	    patphopfneutraliso[ipatpho]                  = myrecopho[ipatpho].recomenneutral;
	    */

	    if(debugPho_)
	      {
		std::cout<<"======some chks for PF iso: if the reco::photon container is sorted right======"<<endl;
		std::cout<<"ipho:"<<ipatpho<<" :Pt from pat photon:"<<patpho->pt()<<" :and from the sorted reco container:"<<myrecopho[ipatpho].pt<<endl;
	      }
	    
	    ////conversion safe electron veto
	    patphopasselectronveto[ipatpho]                  = myrecopho[ipatpho].passelectronveto;

	    /////new HoE
	    patphohadTowOverEm[ipatpho]                 = patpho->hadTowOverEm();
	    /////////////////////////////////////////////////////////////////////////////////////////////////////

	    ////////sigmaiphiiphi - 23rd Aug
	    const reco::CaloClusterPtr phoSeed = (*patpho).superCluster()->seed();
	    vector<float> phoCov;
	    phoCov = lazyTool.localCovariances(*phoSeed);
	    patphosigmaIetaIphi[ipatpho] = sqrt(phoCov[1]);
	    patphosigmaIphiIphi[ipatpho] = sqrt(phoCov[2]);
	    

	    /*patphoTightIDcutword[ipatpho]                = tightphoID( patphoton_container, ipatpho );
	    int cutmask = ~0x0;
	    int cutword = (cutmask^patphoTightIDcutword[ipatpho]);
	    
	    if( cutword == 0x0 )
	      patphotightid[ipatpho] = 1;
	    else
	      patphotightid[ipatpho] = 0;
	    
	    if(debugPho_)
	      {
		std::cout<<"cutword = "<<cutword<<std::endl; 
		std::cout<<"patphoTightIDcutword  for this photon = "<<patphoTightIDcutword[ipatpho]<<std::endl;
		std::cout<<"tight id for this photon              = "<<patphotightid[ipatpho]<<std::endl;
	      }
	    
	    
	    if(debugPho_)
	      {
		if(patphotightid[ipatpho])
		  {
		    //some printouts for sanity checks        
		    if( patphoecalRecHitSumEtConeDR04[ipatpho]>4.2+0.003*patphopt[ipatpho] )
		      std::cout<<"WARNING!!!some problem with filling of patphotightid. ecalRecHitSumEt is supposed to pass this cut"<<std::endl;
		    if( patphohcalTowerSumEtConeDR04[ipatpho]>2.2+0.001*patphopt[ipatpho] )
		      std::cout<<"WARNING!!!some problem with filling of patphotightid. hcalTowerSumEt is supposed to pass this cut"<<std::endl;
		    if( patphotrkSumPtHollowConeDR04[ipatpho]>2.0+0.001*patphopt[ipatpho] )
		      std::cout<<"WARNING!!!some problem with filling of patphotightid. trkSumPtHollowCone is supposed to pass this cut"<<std::endl;
		    if(patpho->isEB() && patphosigmaIetaIeta[ipatpho]>0.013)
		      std::cout<<"WARNING!!!some problem with filling of patphotightid. sigmaietaieta in EB is supposed to pass this cut"<<std::endl;
		    if(patpho->isEE() && patphosigmaIetaIeta[ipatpho]>0.03)
		      std::cout<<"WARNING!!!some problem with filling of patphotightid. sigmaietaieta in EE is supposed to pass this cut"<<std::endl;
		    if(patphohadronicOverEm[ipatpho]>0.05)
		      std::cout<<"WARNING!!!some problem with filling of patphotightid. H/E is supposed to pass this cut"<<std::endl;

		    std::cout<<"haspixel after passing tight id = "<<patphohasPixelSeed[ipatpho]<<std::endl;
		  }//if(patphotightid[ipatpho])                
	      }//if(debugPho_)
	    
	    */
	    
	    //gaps
	    patphoisEB[ipatpho]                                  = patpho->isEB();
	    patphoisEBEEGap[ipatpho]                             = patpho->isEBEEGap();
	    patphoisEBEtaGap[ipatpho]                            = patpho->isEBEtaGap();
	    patphoisEBGap[ipatpho]                               = patpho->isEBGap();
	    patphoisEBPhiGap[ipatpho]                            = patpho->isEBPhiGap();
	    patphoisEE[ipatpho]                                  = patpho->isEE();
	    patphoisEEDeeGap[ipatpho]                            = patpho->isEEDeeGap();
	    patphoisEEGap[ipatpho]                               = patpho->isEEGap();
	    patphoisEERingGap[ipatpho]                           = patpho->isEERingGap();
	    
	    
	    /////mip tagger
	    patphomipChi2[ipatpho]                           = patpho->mipChi2();
	    patphomipIntercept[ipatpho]                      = patpho->mipIntercept();
	    
	    if(patpho->mipIsHalo())
	      patphomipIsHalo[ipatpho]                         = 0;
	    else
	      patphomipIsHalo[ipatpho]                         = 1;
	    
	    patphomipNhitCone[ipatpho]                       = patpho->mipNhitCone();
	    patphomipSlope[ipatpho]                          = patpho->mipSlope();
	    patphomipTotEnergy[ipatpho]                      = patpho->mipTotEnergy();

	    
	    //sceta
	    reco::SuperClusterRef scref = patpho->superCluster();
	    if( !patpho->superCluster().isNonnull() )
	      std::cout<<"scref of patphotons is not valid"<<std::endl;
	    
	    double sceta                         = scref->eta();
	    double scphi                         = scref->phi();
	    double sce                           = scref->energy();
	    
	    double sctheta = (2*atan(exp(-sceta)));
	    double scpx = sce*sin(sctheta)*cos(scphi);
	    double scpy = sce*sin(sctheta)*sin(scphi);
	    double scpz = sce*cos(sctheta);
	    
	    new ((*patphoscp4)[ipatpho]) TLorentzVector();
	    ((TLorentzVector *)patphoscp4->At(ipatpho))->SetXYZT(scpx, scpy, scpz, sce);
	   
	    patphomaxEnergyXtal[ipatpho]                 = patpho->maxEnergyXtal();
	    
	    //patphoscpt[ipatpho]                          = scref->pt();
	    //maximum energy xtal
	 
	    /////29th april, 2012 
	    ////rechit info
	    for(int icrys=0; icrys<MAX_PHOCRYS; icrys++)
	      {
		patphocrysrawId[ipatpho][icrys]   = -9999;
		patphocrysenergy[ipatpho][icrys]  = -9999;
		patphocrystime[ipatpho][icrys]    = -9999;
		patphocrystimeErr[ipatpho][icrys] = -9999;
		patphocrysrecoFlag[ipatpho][icrys]= -9999;
		//patphocryscheckFlag[ipatpho][icrys]= -9999;
		patphocrysieta[ipatpho][icrys]    = -9999;
		patphocrysiphi[ipatpho][icrys]    = -9999;
	      }

	    if(runPhoRechitInfo_){
	      std::vector<std::pair<DetId, float> > phoxtals = scref->hitsAndFractions();
	      vector< std::pair<DetId, float> >::const_iterator detitr;
	      int nphocrys = 0;
	      //int nEEphocrys = 0;
	      for(detitr = phoxtals.begin(); detitr != phoxtals.end(); ++detitr)
		{
		  if( ((*detitr).first).det() == DetId::Ecal && ((*detitr).first).subdetId() == EcalBarrel )
		    {
		      EcalRecHitCollection::const_iterator j= Brechit->find(((*detitr).first));
		      EcalRecHitCollection::const_iterator thishit;
		      if ( j!= Brechit->end())  thishit = j;
		      if ( j== Brechit->end()){
			continue;
		      }
		      EBDetId detId  = (EBDetId)((*detitr).first);
		      patphocrysrawId[ipatpho][nphocrys]  = thishit->id().rawId();
		      patphocrysenergy[ipatpho][nphocrys] = thishit->energy();
		      patphocrystime[ipatpho][nphocrys]   = thishit->time();
		      patphocrystimeErr[ipatpho][nphocrys] = thishit->timeError();
		      patphocrysrecoFlag[ipatpho][nphocrys] = thishit->recoFlag();
		      //patphocryscheckFlag[ipatpho][nphocrys] = thishit->checkFlag();
		      patphocrysieta[ipatpho][nphocrys]   = detId.ieta();
		      patphocrysiphi[ipatpho][nphocrys]   = detId.iphi();
		      nphocrys++;
		    }// if( ((*detitr).first).det() == DetID::Ecal && ((*detitr..
		  
		  
		  else if (((*detitr).first).det() == DetId::Ecal && ((*detitr).first).subdetId() == EcalEndcap){
		    EcalRecHitCollection::const_iterator j= Erechit->find(((*detitr).first));
		    EcalRecHitCollection::const_iterator thishit;
		    if ( j!= Erechit->end())  thishit = j;
		    if ( j== Erechit->end()){
		      continue;
		    }
		    EEDetId detId  = (EEDetId)((*detitr).first);
		    patphocrysenergy[ipatpho][nphocrys] = thishit->energy();
		    patphocrystime[ipatpho][nphocrys]   = thishit->time();
		    patphocrystimeErr[ipatpho][nphocrys] = thishit->timeError();
		    patphocrysrecoFlag[ipatpho][nphocrys] = thishit->recoFlag(); 
		    //patphocryscheckFlag[ipatpho][nphocrys] = thishit->checkFlag(); 
		    patphocrysrawId[ipatpho][nphocrys]  = -9999;
		    patphocrysieta[ipatpho][nphocrys]   = -9999;
		    patphocrysiphi[ipatpho][nphocrys]   = -9999;
		    nphocrys++;
		  }//else if (((*detitr).first).det() == DetId::Ecal && ((*d..EF)
		}//
	      
	      patphoncrys[ipatpho] = nphocrys;
	      if(nphocrys > MAX_PHOCRYS)
		{

		  cout<<"WARNING FROM ANALYZER.CC!!!!Number of crystals to be stored for a photon increased its maximum storage!!!PLEASE FIX:nphocrys:MAX_PHOCRYS:"<<nphocrys<<":"<<MAX_PHOCRYS<<endl;
		}
	    }//if(runPhoRechitInfo_)
	    


	    ///from UCSD
	    int index = 0;
	    
	    if(debugPho_)
              {
                cout<<"uperClustersHybridH->size() = "<<superClustersHybridH->size()<<endl;
		cout<<"superClustersEndcapH->size() = "<<superClustersEndcapH->size()<<endl;
              }
	    
	    patphoscind[ipatpho] = -1;
	    
	    reco::Photon localPho = reco::Photon(*patpho);

	    for(int isuperClusterType=0; isuperClusterType<2; ++isuperClusterType) { //in UCSD isuperClusterType<3
	      if (isuperClusterType == 0) {
		for(reco::SuperClusterCollection::size_type j = 0; j<superClustersHybridH->size(); ++j){
		  
		  reco::SuperClusterRef sc(superClustersHybridH, j);
		  
		  if(debugPho_)
		    {
		      cout<<"index EB= "<<endl;
		      //cout<<"(&(*(patpho->superCluster())) = "<<(&(*(patpho->superCluster())))<<endl;
		      cout<<"(&(*localPho.superCluster()) == &(*sc)) = "<<(&(*localPho.superCluster()) == &(*sc))<<endl;
		      cout<<"&(*sc) = "<<&(*sc)<<endl;
		      cout<<"(&(*localPho.superCluster()) = "<<(&(*localPho.superCluster()))<<endl;
		    }
		  
		  //if (&(*(patpho->superCluster())) == &(*sc)) { //address inside patpho->superCluster==add inside sc
		  if(&(*localPho.superCluster()) == &(*sc)){
		    

		    patphoscind[ipatpho] = index;
		    break;
		  }
		  index++;
		}//for(reco::SuperClusterCollection::size_type j = 0; j<superClustersHybridH->size(); ++j)
	      }
	      
	      if (isuperClusterType == 1) { //in UCSD isuperClusterType == 2
		for(reco::SuperClusterCollection::size_type j = 0; j<superClustersEndcapH->size(); ++j){
		  
		  reco::SuperClusterRef sc(superClustersEndcapH, j);
		  if(debugPho_)
		    {
		      cout<<"index EE= "<<endl;
		      //cout<<"(&(*(patpho->superCluster())) = "<<(&(*(patpho->superCluster())))<<endl;
		      cout<<"(&(*localPho.superCluster()) == &(*sc)) = "<<(&(*localPho.superCluster()) == &(*sc))<<endl;
		      cout<<"&(*sc) = "<<&(*sc)<<endl;
		      cout<<"(&(*localPho.superCluster()) = "<<(&(*localPho.superCluster()))<<endl;
		    }

		  
		  //if (&(*(patpho->superCluster())) == &(*sc)) {
		  if(&(*localPho.superCluster()) == &(*sc)){
		    patphoscind[ipatpho] = index;
		    break;
		  }
		  index++;
		}//for(reco::SuperClusterCollection::size_type j = 0; j<superClustersEndcapH->size(); ++j)
	      }//if (isuperClusterType == 2)
	    }//for(int isuperClusterType=0; isuperClusterType<3; ++isuperClusterType)

	    //int scindex = patphoscindex[ipatpho];

	    if(runSCremoval_){
	      ////done on 25th Dec
	      //SuperClusterFootprintRemoval remover(iEvent,myiConfig,iSetup);
	      SuperClusterFootprintRemoval remover(iEvent,iSetup, myiConfig);
	      //SuperClusterFootprintRemoval remover(iEvent,iSetup);
	      //SuperClusterRef scref(reco::SuperClusterRef(scHandle,scindex));
	      int vertexforchargediso = 0; // this is the index of the vertex selected in the event (see the following section for explanation)
	      // IF YOU WANT TO CALCULATE PF ISO
	      float chargediso = remover.PFIsolation("charged",scref,vertexforchargediso);
	      float neutraliso = remover.PFIsolation("neutral",scref);
	      float photoniso  = remover.PFIsolation("photon",scref);
	      
	      ////refilling the PF isolation
	      /////////////////////////////////////////////////earlier filled here:
	      //////////////////////////////////////20th July///////////////////////////////////////////////
	      /////PF isolation info
	      ////////////////////////////////////////////////end of earlier filled here
	      
	      patphonewpfchargediso[ipatpho]                  = chargediso;
	      patphonewpfphotoniso[ipatpho]                   = photoniso;
	      patphonewpfneutraliso[ipatpho]                  = neutraliso;
	    }

	    const reco::CaloClusterPtr seed = scref->seed();
	    
	    //spike info
	    if(patpho->isEB())
	      {
		patphoswissCross[ipatpho]   = EcalClusterTools::eTop( *(patpho->superCluster()->seed()), &(*barrelRecHits), &(*topology))+ EcalClusterTools::eBottom( *(patpho->superCluster()->seed()), &(*barrelRecHits), &(*topology)) + EcalClusterTools::eLeft( *(patpho->superCluster()->seed()), &(*barrelRecHits), &(*topology)) + EcalClusterTools::eRight( *(patpho->superCluster()->seed()), &(*barrelRecHits), &(*topology));
		//std::cout<<"etop: "<<EcalClusterTools::eTop( *(patpho->superCluster()->seed()), &(*barrelRecHits), &(*topology))<<endl;
		//std::cout<<"etop: "<<EcalClusterTools::eTop( *(patpho->superCluster()->seed()), &(*barrelRecHits), &(*topology))<<endl;
		//if(1-patphoswissCross[ipatpho]/patphomaxEnergyXtal[ipatpho] > 0.95) cout<<"This phoctron candidate is an ECAL spike identified by Swiss Cross algorithm."<<endl;
		   
		patphoswissBasedspikevar[ipatpho]             = 1-patphoswissCross[ipatpho]/patphomaxEnergyXtal[ipatpho];
	      }//end of if(patpho->isEB())
	    else{
	      patphoswissCross[ipatpho]   = EcalClusterTools::eTop( *(patpho->superCluster()->seed()), &(*endcapRecHits), &(*topology))+ EcalClusterTools::eBottom( *(patpho->superCluster()->seed()), &(*endcapRecHits), &(*topology)) + EcalClusterTools::eLeft( *(patpho->superCluster()->seed()), &(*endcapRecHits), &(*topology)) + EcalClusterTools::eRight( *(patpho->superCluster()->seed()), &(*endcapRecHits), &(*topology));
	      //if(1-patphoswissCross[ipatpho]/patphomaxEnergyXtal[ipatpho] > 0.95) cout<<"This phoctron candidate is an ECAL spike identified by Swiss Cross algorithm."<<endl;
	      patphoswissBasedspikevar[ipatpho]             = 1-patphoswissCross[ipatpho]/patphomaxEnergyXtal[ipatpho];
	       }//end of else
	    
	 
	        //timing info
	    DetId patphoSeedDetId = lazyTool.getMaximum(*seed).first;
	    
	    if ( patpho->isEB() ) {
	      EcalRecHitCollection::const_iterator phoebrhit = barrelRecHits->find(patphoSeedDetId);
	      if ( phoebrhit != barrelRecHits->end() ) {
		patphoseedtime[ipatpho]         = phoebrhit->time();
		patphorecoFlag[ipatpho]         = phoebrhit->recoFlag();
		//patphocheckFlag[ipatpho]         = phoebrhit->checkFlag();
		//patphokOutOfTime[ipatpho]       = phoebrhit->kOutOfTime();
		//patphoseverityLevel[ipatpho]    = EcalSeverityLevelAlgo::severityLevel( patphoSeedDetId, (*barrelRecHits), *chStatus );
		
		patphoe2e9[ipatpho] = mye2overe9(patphoSeedDetId,*barrelRecHits);
		
		if(debugPho_)
		  std::cout<<"patphoe2e9 = "<<patphoe2e9[ipatpho]<<std::endl;
		
	      }//if ( phoebrhit != barrelRecHits->end() )
	    }//if ( patpho->isEB() && barrelRecHits.isValid() )
            
	    else{
	      EcalRecHitCollection::const_iterator phoeerhit = endcapRecHits->find(patphoSeedDetId);
	      if ( phoeerhit != endcapRecHits->end() ) {
		patphoseedtime[ipatpho]         = phoeerhit->time();
		patphorecoFlag[ipatpho]         = phoeerhit->recoFlag();
		//patphocheckFlag[ipatpho]         = phoeerhit->checkFlag();
		//patphokOutOfTime[ipatpho]       = phoeerhit->kOutOfTime();
		//patphoseverityLevel[ipatpho]    = EcalSeverityLevelAlgo::severityLevel( patphoSeedDetId, (*endcapRecHits), *chStatus );

		patphoe2e9[ipatpho] = 0;
	      }//if ( phoeerhit != endcapRecHits->end() )
	    }//else if ( endcapRecHits.isValid() ) 
	    
	    //conversions
	    reco::ConversionRefVector conversions  = patpho->conversions();
	    patphohasConversionTracks[ipatpho] = patpho->hasConversionTracks();
	    if(conversions.size()>0)
	      {
		if(debugPho_)
		  {
		    std::cout<<"no of converted tracks of photons:"<<conversions.size()<<std::endl;
		  }
		
		for (unsigned int iConv=0; iConv<conversions.size(); iConv++) {
		  reco::ConversionRef aConv=conversions[iConv];
		  
		  if ( aConv->conversionVertex().isValid() )
		    {
		      if(debugPho_)
			{
			  double ntrks = aConv->nTracks();
			  std::cout<<"ntrks of conv:"<<ntrks<<std::endl;
			}

		      const reco::Vertex&  vtx                 = aConv->conversionVertex();
		      patphoconvtxX[ipatpho][iConv]      = vtx.x();
		      patphoconvtxY[ipatpho][iConv]      = vtx.y();
		      patphoconvtxZ[ipatpho][iConv]      = vtx.z();
		      patphoconvtxR[ipatpho][iConv]      = sqrt( pow(patphoconvtxX[ipatpho][iConv],2) + pow(patphoconvtxY[ipatpho][iConv],2) + pow(patphoconvtxZ[ipatpho][iConv],2) );
		    }//if ( aConv->conversionVertex().isValid() )
		}// for (unsigned int iConv=0; iConv<conversions.size(); iConv++)
		patphoconvsize[ipatpho] = (int)conversions.size();
	      }//f(patpho->conversions().isNonnull())
	    
	       //trackinfo
	    reco::GsfTrackRef trackref = patpho->gsfTrack();
	    if(patpho->gsfTrack().isNonnull())
	      {
		patphotrkpt[ipatpho]                         = trackref->pt();
		patphotrkcharge[ipatpho]                     = trackref->charge();
		patphotrkchi2[ipatpho]                       = trackref->chi2();
		patphotrketa[ipatpho]                        = trackref->eta();
		patphotrknumberOfLostHits[ipatpho]           = trackref->numberOfLostHits();
		patphotrknumberOfValidHits[ipatpho]          = trackref->numberOfValidHits();
		patphotrklost[ipatpho]                       = trackref->lost();
		patphotrkd0[ipatpho]                         = trackref->d0();
		patphotrkdxy[ipatpho]                        = trackref->dxy();
		patphotrkdz[ipatpho]                         = trackref->dz();
		patphotrkptin[ipatpho]                       = sqrt(trackref->innerMomentum().Perp2());
		patphotrkptout[ipatpho]                      = sqrt(trackref->outerMomentum().Perp2());
		patphotrkfbrem[ipatpho]                      = (patphotrkptin[ipatpho]-patphotrkptout[ipatpho])/patphotrkptin[ipatpho];
		patphotrkqoverp[ipatpho]                     = trackref->qoverp();
		patphotrkvx[ipatpho]                         = trackref->vx();
		patphotrkvy[ipatpho]                         = trackref->vy();
		patphotrkvz[ipatpho]                         = trackref->vz();
		patphotrkphi[ipatpho]                        = trackref->phi();
		patphotrkndof[ipatpho]                       = trackref->ndof();
		patphotrkrecHitsSize[ipatpho]                = trackref->recHitsSize();
		patphotrktheta[ipatpho]                      = trackref->theta();
		patphotrkqualityMask[ipatpho]                = trackref->qualityMask();
		patphotrkouterX[ipatpho]                     = trackref->outerX();
		patphotrkouterY[ipatpho]                     = trackref->outerY();
		patphotrkouterZ[ipatpho]                     = trackref->outerZ();
		patphotrkouterRadius[ipatpho]                = trackref->outerRadius();
		patphotrkinnerX[ipatpho]                     = trackref->innerPosition().X();
		patphotrkinnerY[ipatpho]                     = trackref->innerPosition().Y();
		patphotrkinnerZ[ipatpho]                     = trackref->innerPosition().Z();
	      }//if(patphoctron_container[ipatpho].gsfTrack().isNonnull())
	    
	       //gen info
	    if( !isData )
	      {
		if( patpho->genParticleRef().isNonnull() )
		  {
		    patphomGenstatus[ipatpho]                   = patpho->genParticleRef()->status();
		    patphomGenpdgId[ipatpho]                    = patpho->genParticleRef()->pdgId();
		    patphomGentheta[ipatpho]                    = patpho->genParticleRef()->theta();
		    patphomGeneta[ipatpho]                      = patpho->genParticleRef()->eta();
		    patphomGenphi[ipatpho]                      = patpho->genParticleRef()->phi();
		    patphomGenpt[ipatpho]                       = patpho->genParticleRef()->pt();
		    patphomGenpx[ipatpho]                       = patpho->genParticleRef()->px();
		    patphomGenpy[ipatpho]                       = patpho->genParticleRef()->py();
		    patphomGenpz[ipatpho]                       = patpho->genParticleRef()->pz();
		    patphomGenenergy[ipatpho]                   = patpho->genParticleRef()->energy();
		    patphomGenisJet[ipatpho]                    = patpho->genParticleRef()->isJet();
		    patphomGenisPhoton[ipatpho]                 = patpho->genParticleRef()->isPhoton();
		    
		    patphonummoth[ipatpho] = (int)patpho->genParticleRef()->numberOfMothers();
		    if(patphonummoth[ipatpho]!=0)
		      {
			for(int imoth=0; imoth<patphonummoth[ipatpho]; imoth++)
			  patphomGenmompdgId[ipatpho][imoth]                 = patpho->genParticleRef()->mother(imoth)->pdgId();
		      }
		    
		    int numgranmoth = (int)patpho->genParticleRef()->mother()->numberOfMothers();
		    if(numgranmoth!=0)
		      patphomGengranpdgId[ipatpho]                = patpho->genParticleRef()->mother()->mother()->pdgId();
		    
		    patphomGenisElectron[ipatpho]               = patpho->genParticleRef()->isElectron();
		    
		  }
		
		if(debugPho_)
		  {
		    std::cout<<"patphomGenpdgId["<<ipatpho<<"] = "<<patphomGenpdgId[ipatpho]<<std::endl;
		    std::cout<<"patphomGenmompdgId["<<ipatpho<<"] = "<<patphomGenmompdgId[ipatpho]<<std::endl;
		    std::cout<<"patphomGenpt["<<ipatpho<<"] = "<<patphomGenpt[ipatpho]<<std::endl;
		    std::cout<<"patphomGenpx["<<ipatpho<<"] = "<<patphomGenpx[ipatpho]<<std::endl;
		    std::cout<<"patphomGenpy["<<ipatpho<<"] = "<<patphomGenpy[ipatpho]<<std::endl;
		    std::cout<<"patphomGenpz["<<ipatpho<<"] = "<<patphomGenpz[ipatpho]<<std::endl;
		    std::cout<<"patphomGenenergy["<<ipatpho<<"] = "<<patphomGenenergy[ipatpho]<<std::endl;
		  }
	      }//if( !isData )
	    ipatpho++;
	 }//for(edm::View<pat::Photon>::const_iterator patpho = photons.begin(); patpho!=photons.end(); ++patpho)
       patphosize = ipatpho;
     }//if(runPatphotons_)
   




   //PFMET
  if(runPFmet_)
     {
       //met collection
       edm::Handle<edm::View<pat::MET> > pfmetH;
       bool foundpfmet = iEvent.getByLabel(PFmetLabel_,pfmetH);
       if(!foundpfmet ) {
	 std::ostringstream  err;
	 err<<" cannot get pfmet "
	    <<std::endl;
	 edm::LogError("met")<<err.str();
	 throw cms::Exception( "MissingProduct", err.str());
       }


       //PFmet
       for(int imet=0; imet<25; imet++)
	 {
	   pfmetpt[imet]                                 = -999.0;
	   pfmetphi[imet]                                = -999.0;
	   pfmetsumEt[imet]                              = -999.0;
	   pfmetpx[imet]                                 = -999.0;
	   pfmetpy[imet]                                 = -999.0;
	 }
       pfmetsize = -999;
       

       int pfmetsize = 0;
       int pfnmet    = 0;

       

       
       pfmetsize = pfmetH->size();
       if( pfmetsize > 1 ) std::cout<<"how is it possible! met size = "<<pfmetsize<<std::endl;
       
       int imet = 0;
       for( edm::View<pat::MET>::const_iterator pfmetit = pfmetH->begin(); pfmetit !=pfmetH->end(); pfmetit++ )
	 {
	   pfmetpt[imet]    = pfmetit->pt();
           pfmetphi[imet]    = correct_phi(pfmetit->phi());
           pfmetsumEt[imet] = pfmetit->sumEt();
           pfmetpx[imet]    = pfmetit->px();
           pfmetpy[imet]    = pfmetit->py();

	   imet++;
	 }

     }//if(runPFmet_)



   //TCMET
   if(runTCmet_)
     {
       //met collection
       edm::Handle<edm::View<pat::MET> > tcmetH;
       bool foundtcmet = iEvent.getByLabel(TCmetLabel_,tcmetH);
       if(!foundtcmet ) {
         std::ostringstream  err;
         err<<" cannot get tc met "
            <<std::endl;
         edm::LogError("met")<<err.str();
         throw cms::Exception( "MissingProduct", err.str());
       }


       //PFmet
       for(int imet=0; imet<25; imet++)
         {
           tcmetpt[imet]                                 = -999.0;
           tcmetphi[imet]                                = -999.0;
           tcmetsumEt[imet]                              = -999.0;
           tcmetpx[imet]                                 = -999.0;
           tcmetpy[imet]                                 = -999.0;
         }
       tcmetsize = -999;


       int tcmetsize = 0;
       int tcnmet    = 0;

  tcmetsize = tcmetH->size();
       if( tcmetsize > 1 ) std::cout<<"how is it possible! met size = "<<tcmetsize<<std::endl;

       int imet = 0;
       for( edm::View<pat::MET>::const_iterator tcmetit = tcmetH->begin(); tcmetit != tcmetH->end(); tcmetit++ )
         {
           tcmetpt[imet]    = tcmetit->pt();
           tcmetphi[imet]   = correct_phi(tcmetit->phi());
           tcmetsumEt[imet] = tcmetit->sumEt();
           tcmetpx[imet]    = tcmetit->px();
           tcmetpy[imet]    = tcmetit->py();
           imet++;
         }

     }//if(runtcmet_)



   //MET
   if(runmet_)
     {
       //met collection
       edm::Handle<edm::View<pat::MET> > met;
       bool foundmet = iEvent.getByLabel(metLabel_,met);
       if(!foundmet ) {
         std::ostringstream  err;
         err<<" cannot get met "
            <<std::endl;
         edm::LogError("met")<<err.str();
         throw cms::Exception( "MissingProduct", err.str());
       }


       //met
       for(int imet=0; imet<25; imet++)
         {
           metpt[imet]                                 = -999.0;
           metphi[imet]                                = -999.0;
           metsumEt[imet]                              = -999.0;
           metpx[imet]                                 = -999.0;
           metpy[imet]                                 = -999.0;
         }
       metsize = -999;


       int metsize = 0;
       int nmet    = 0;


       metsize = met->size();
       if( metsize > 1 ) std::cout<<"how is it possible! met size = "<<metsize<<std::endl;

       int imet = 0;
       for( edm::View<pat::MET>::const_iterator metit = met->begin(); metit !=met->end(); metit++ )
         {
           metpt[imet]    = metit->pt();
           metphi[imet]   = correct_phi(metit->phi());
           metsumEt[imet] = metit->sumEt();
           metpx[imet]    = metit->px();
           metpy[imet]    = metit->py();

           imet++;
         }

     }//if(runmet_)



   //pat jets
   if(runpatJets_)
     {
       
       edm::Handle<edm::View<pat::Jet> > jetHandle;
       iEvent.getByLabel(jetLabel_,jetHandle);
       const edm::View<pat::Jet> & jets = *jetHandle;
       std::vector<pat::Jet>  patJet_container;
       for(edm::View<pat::Jet>::const_iterator jet_iter = jets.begin(); jet_iter!=jets.end(); ++jet_iter){
	 patJet_container.push_back(*jet_iter);
       }
       
       if(patJet_container.size()!=0){
         for(unsigned int ijet=0;ijet < patJet_container.size();ijet++)
           {
	     patjetpt[ijet]                               = patJet_container[ijet].pt();
	     patjetpx[ijet]                               = patJet_container[ijet].px();
	     patjetpy[ijet]                               = patJet_container[ijet].py();
	     patjetpz[ijet]                               = patJet_container[ijet].pz();
	     patjetet[ijet]                               = patJet_container[ijet].et();
	     patjetenergy[ijet]                           = patJet_container[ijet].energy();
	     patjetphi[ijet]                              = correct_phi(patJet_container[ijet].phi());
	     patjeteta[ijet]                              = patJet_container[ijet].eta();
	     patjethasOverlapsmu[ijet]                    = patJet_container[ijet].hasOverlaps("muons");
	     patjethasOverlapsele[ijet]                   = patJet_container[ijet].hasOverlaps("electrons");
	     patjethasOverlapspho[ijet]                   = patJet_container[ijet].hasOverlaps("photons");
	     patjethasOverlapstau[ijet]                   = patJet_container[ijet].hasOverlaps("taus");
	     patjethasOverlapstkIsoele[ijet]              = patJet_container[ijet].hasOverlaps("tkIsoElectrons");

	     /*patjetchargedEmEnergy[ijet]                  = patJet_container[ijet].chargedEmEnergy();
	     patjetchargedEmEnergyFraction[ijet]          = patJet_container[ijet].chargedEmEnergyFraction();
	     patjetchargedHadronEnergy[ijet]              = patJet_container[ijet].chargedHadronEnergy();
	     patjetchargedHadronEnergyFraction[ijet]      = patJet_container[ijet].chargedHadronEnergyFraction();
	     patjetchargedMultiplicity[ijet]              = patJet_container[ijet].chargedMultiplicity();
	     patjetemEnergyFraction[ijet]                 = patJet_container[ijet].emEnergyFraction();
	     patjetemEnergyInEB[ijet]                     = patJet_container[ijet].emEnergyInEB();
	     patjetemEnergyInEE[ijet]                     = patJet_container[ijet].emEnergyInEE();
	     patjetemEnergyInHF[ijet]                     = patJet_container[ijet].emEnergyInHF();
	     patjetenergyFractionHadronic[ijet]           = patJet_container[ijet].energyFractionHadronic();
	     patjethadEnergyInHB[ijet]                    = patJet_container[ijet].hadEnergyInHB();
	     patjethadEnergyInHE[ijet]                    = patJet_container[ijet].hadEnergyInHE();
	     patjethadEnergyInHF[ijet]                    = patJet_container[ijet].hadEnergyInHF();
	     patjethadEnergyInHO[ijet]                    = patJet_container[ijet].hadEnergyInHO();
	     */
	     //patjethadEnergy[ijet]         = patJet_container[ijet].hadEnergy();
	   }//for(int ijet=0; ijet<patJet_container.size(); ijet++)
	 
	 
       }//if(patJet_container.size()!=0)
       patjetsize = patJet_container.size();
     }//if(runpatJets_)





   int ijet = 0;
   
   //pat jets
   if(runpfJets_)
     {
       //for jec uncert
       edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
       iSetup.get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl); 
       //JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
       //JetCorrectionUncertainty *pfjecUnc = new JetCorrectionUncertainty(JetCorPar);


       edm::Handle<edm::View<pat::Jet> > pfjetHandle;
       iEvent.getByLabel(pfjetLabel_,pfjetHandle);
       const edm::View<pat::Jet> & pfjets = *pfjetHandle;
  
       
       for(edm::View<pat::Jet>::const_iterator jet_iter = pfjets.begin(); jet_iter!=pfjets.end(); ++jet_iter){
	 pfjetpt[ijet]                               = jet_iter->pt();
	 pfjetpx[ijet]                               = jet_iter->px();
	 pfjetpy[ijet]                               = jet_iter->py();
	 pfjetpz[ijet]                               = jet_iter->pz();
	 pfjetet[ijet]                               = jet_iter->et();
	 pfjetenergy[ijet]                           = jet_iter->energy();
	 pfjetphi[ijet]                              = correct_phi(jet_iter->phi());
	 pfjeteta[ijet]                              = jet_iter->eta();
	 pfjethasOverlapsmu[ijet]                    = jet_iter->hasOverlaps("muons");
	 pfjethasOverlapsele[ijet]                   = jet_iter->hasOverlaps("electrons");
	 pfjethasOverlapspho[ijet]                   = jet_iter->hasOverlaps("photons");
	 pfjethasOverlapstau[ijet]                   = jet_iter->hasOverlaps("taus");
	 pfjethasOverlapstkIsoele[ijet]              = jet_iter->hasOverlaps("tkIsoElectrons");
	 
	 ///added on 31st October, 2011 for fake rate
	 pfjetchargedEmEnergyFraction[ijet]          = jet_iter->chargedEmEnergyFraction();
	 pfjetchargedHadronEnergyFraction[ijet]      = jet_iter->chargedHadronEnergyFraction();
	 pfjetchargedMuEnergyFraction[ijet]          = jet_iter->chargedMuEnergyFraction();
	 

	 if(jet_iter->jecFactor("Uncorrected")!= 0)
	   {pfjet_jecCorr[ijet]  = (1.0/jet_iter->jecFactor("Uncorrected")); 
	   }
         else{pfjet_jecCorr[ijet] =0.;}


	 
	 //patjethadEnergy[ijet]         = jet_iter->hadEnergy();
	 
	 ijet++;
       }//if(runpatJets_)
       pfjetsize = ijet;
   
     }//if(runpfJets_)
   


      //----------------------------------------------Trigger info--------------------------------------------------------     
   /////HLT
   int triggerIndex[100];  //for HLT-RECO object matching 
   if(runHLT_==1){
     Handle<TriggerResults> HLTR;
     iEvent.getByLabel(hlTriggerResults_,HLTR);
     
     int indexhltobj;  //for HLT-RECO object

     //get prescales here 
     /*int iprescale = 0;
     if(nhlt_>100)
       std::cout<<"WARNING!Analyzer.cc can accomodate only atmost 100 triggers to match. not more than that. change the array size of hltprescale to accomodate that much"<<std::endl;
     
     for(std::vector<std::string>::const_iterator itr=hltMatch_.begin();itr!=hltMatch_.end();itr++)
       {
	 hltprescale[iprescale] =  hltConfig_.prescaleValue(iEvent, iSetup, (*itr)) ;
	 if(debugHLT_==1)
	   {
	     std::cout<<"Name of trigger to match = "<<(*itr)<<std::endl;
	     std::cout<<"event = "<<event<<" hltprescale["<<iprescale<<"] = "<<hltprescale[iprescale]<<std::endl;
	   }
	      
	 iprescale++;
       }
     
     if(debugHLT_==1)
       {
       if(iprescale!=nhlt_)
	 std::cout<<"no. of triggers you WANT to match not equal to the ones PROVIDED(iprescale!=nhlt_)! correct it!"<<std::endl;
       }
     */

     ntrigToMatch = nhlt_;

     //if (!init_) {
       //init_=true;
       //init_=false;     
       //triggerNames_.init(*HLTR);
       const TriggerNames  &triggerNames_ = iEvent.triggerNames(*HLTR);
       hlNames_=triggerNames_.triggerNames();
       
       /*int itrigindex = 0;
       for(std::vector<std::string>::const_iterator itr=hltMatch_.begin();itr!=hltMatch_.end();itr++)
	 {
	   triggerIndex[itrigindex] = triggerNames_.triggerIndex((*itr));  //for HLT-RECO object matching
	   itrigindex++;
	 }//for(std::vector<std::string>::const
       */
       //}

     // decision for each HL algorithm
     const unsigned int n(hlNames_.size());
     //std::cout<<"n = "<<n<<std::endl;
     for(unsigned int i = 0; i<n ;i++)
       {
         //cout<<"inside loop"<<endl;
         //cout<<hlNames_[i]<<" :"<<HLTR->accept(i)<<endl;
	 HLT_chosen[ hlNames_[i]]=HLTR->accept(i);

	 if(hlNames_[i]=="HLT_Ele10_SW_L1R")
           is_HLT_Ele10_SW_L1R_event=HLTR->accept(i);
         //cout<<"is_HLT_Ele10_SW_L1R_event = "<<is_HLT_Ele10_SW_L1R_event<<endl;
	 
         if(hlNames_[i]=="HLT_Ele15_SW_L1R")
           is_HLT_Ele15_SW_L1R_event=HLTR->accept(i);

         if(hlNames_[i]=="HLT_Ele15_SW_EleId_L1R")
           is_HLT_Ele15_SW_EleId_L1R_event=HLTR->accept(i);

         if(hlNames_[i]=="HLT_Ele15_SW_LooseTrackIso_L1R")
           is_HLT_Ele15_SW_LooseTrackIso_L1R_event=HLTR->accept(i);

         if(hlNames_[i]=="HLT_Ele15_SC15_SW_LooseTrackIso_L1R")
           is_HLT_Ele15_SC15_SW_LooseTrackIso_L1R_event=HLTR->accept(i);
	 
	 if(hlNames_[i]=="HLT_Ele15_SC15_SW_EleId_L1R")
           is_HLT_Ele15_SC15_SW_EleId_L1R_event=HLTR->accept(i);

         if(hlNames_[i]=="HLT_Ele20_SW_L1R")
           is_HLT_Ele20_SW_L1R_event=HLTR->accept(i);

         if(hlNames_[i]=="HLT_Ele20_SC15_SW_L1R")
           is_HLT_Ele20_SC15_SW_L1R_event=HLTR->accept(i);

         if(hlNames_[i]=="HLT_Ele25_SW_L1R")
           is_HLT_Ele25_SW_L1R_event=HLTR->accept(i);

         if(hlNames_[i]=="HLT_Ele25_SW_EleId_LooseTrackIso_L1R")
           is_HLT_Ele25_SW_EleId_LooseTrackIso_L1R_event=HLTR->accept(i);

         if(hlNames_[i]=="HLT_DoubleEle5_SW_Jpsi_L1R")
           is_HLT_DoubleEle5_SW_Jpsi_L1R_event=HLTR->accept(i);
	 
	 if(hlNames_[i]=="HLT_DoubleEle5_SW_Upsilon_L1R")
           is_HLT_DoubleEle5_SW_Upsilon_L1R_event=HLTR->accept(i);

         if(hlNames_[i]=="HLT_DoubleEle10_SW_L1R")
           is_HLT_DoubleEle10_SW_L1R_event=HLTR->accept(i);

	 if(hlNames_[i]=="HLT_L1SingleEG5")
	   is_HLT_L1SingleEG5_event=HLTR->accept(i);

	 if(hlNames_[i]=="HLT_L1SingleEG8")
	   is_HLT_L1SingleEG8_event=HLTR->accept(i);

	 if(hlNames_[i]=="HLT_Ele10_LW_L1R")
	   is_HLT_Ele10_LW_L1R_event=HLTR->accept(i);
	 
	 if(hlNames_[i]=="HLT_Ele10_LW_EleId_L1R")
	   is_HLT_Ele10_LW_EleId_L1R_event=HLTR->accept(i);
	 
	 if(hlNames_[i]=="HLT_Ele15_LW_L1R")
	   is_HLT_Ele15_LW_L1R_event=HLTR->accept(i);

	 if(hlNames_[i]=="HLT_Ele15_SiStrip_L1R")
	   is_HLT_Ele15_SiStrip_L1R_event=HLTR->accept(i);

	 if(hlNames_[i]=="HLT_Ele15_SC10_LW_L1R")
	   is_HLT_Ele15_SC10_LW_L1R_event=HLTR->accept(i);

	 if(hlNames_[i]=="HLT_Ele20_LW_L1R")
	   is_HLT_Ele20_LW_L1R_event=HLTR->accept(i);
	 
	 if(hlNames_[i]=="HLT_L1DoubleEG5")
	   is_HLT_L1DoubleEG5_event=HLTR->accept(i);
	 
	 if(hlNames_[i]=="HLT_DoubleEle5_SW_L1R")
	   is_HLT_DoubleEle5_SW_L1R_event=HLTR->accept(i);

	 if(hlNames_[i]=="HLT_Photon10_L1R")
           is_HLT_Photon10_L1R_event=HLTR->accept(i);

         if(hlNames_[i]=="HLT_Photon10_LooseEcalIso_TrackIso_L1R")
           is_HLT_Photon10_LooseEcalIso_TrackIso_L1R_event=HLTR->accept(i);

         if(hlNames_[i]=="HLT_Photon15_L1R")
           is_HLT_Photon15_L1R_event=HLTR->accept(i);

         if(hlNames_[i]=="HLT_Photon20_LooseEcalIso_TrackIso_L1R")
           is_HLT_Photon20_LooseEcalIso_TrackIso_L1R_event=HLTR->accept(i);

         if(hlNames_[i]=="HLT_Photon25_L1R")
           is_HLT_Photon25_L1R_event=HLTR->accept(i);

         if(hlNames_[i]=="HLT_Photon25_LooseEcalIso_TrackIso_L1R")
           is_HLT_Photon25_LooseEcalIso_TrackIso_L1R_event=HLTR->accept(i);
	 
	 if(hlNames_[i]=="HLT_Photon30_L1R_1E31")
           is_HLT_Photon30_L1R_1E31_event=HLTR->accept(i);

         if(hlNames_[i]=="HLT_Photon70_L1R")
           is_HLT_Photon70_L1R_event=HLTR->accept(i);

         if(hlNames_[i]=="HLT_DoublePhoton10_L1R")
           is_HLT_DoublePhoton10_L1R_event=HLTR->accept(i);

         if(hlNames_[i]=="HLT_DoubleIsoPhoton15_L1R")
           is_HLT_DoublePhoton15_L1R_event=HLTR->accept(i);

         if(hlNames_[i]=="HLT_DoubleIsoPhoton15_VeryLooseEcalIso_L1R")
           is_HLT_DoublePhoton15_VeryLooseEcalIso_L1R_event=HLTR->accept(i);

	 if(hlNames_[i]=="HLT_Photon20_Cleaned_L1R"){
	   is_HLT_Photon20_Cleaned_L1R_event=HLTR->accept(i);
          // std::cout<<"Hlt 20 found"<<std::endl;
          }

	 if(hlNames_[i]=="HLT_Photon30_Cleaned_L1R") {
	   is_HLT_Photon30_Cleaned_L1R_event=HLTR->accept(i);
          // std::cout<<"HLT 30 found"<<std::endl; 
          }

	 if(hlNames_[i]=="HLT_Photon50_Cleaned_L1R"){
	   is_HLT_Photon50_Cleaned_L1R_event=HLTR->accept(i);
          //std::cout<<"HLT 50 found "<<std::endl;
          }


	 if(hlNames_[i]=="HLT_Photon70_Cleaned_L1R"){
	   is_HLT_Photon70_Cleaned_L1R_event=HLTR->accept(i);
          // std::cout<<"HLT 70 found"<<std::endl;
          }

	 //for e*
	 /*if(hlNames_[i]=="HLT_DoublePhoton20_L1R")
           is_HLT_DoublePhoton20_L1R_event=HLTR->accept(i);

	 if(hlNames_[i]=="HLT_DoublePhoton17_L1R")
           is_HLT_DoublePhoton17_L1R_event=HLTR->accept(i);

	 if(hlNames_[i]=="HLT_DoublePhoton17_SC17HE_L1R")
           is_HLT_DoublePhoton17_SC17HE_L1R_event=HLTR->accept(i);

	 if(hlNames_[i]=="HLT_DoublePhoton22_SC22HE_L1R")
           is_HLT_DoublePhoton10_L1R_event=HLTR->accept(i);
	 */

	 //jet
	 unsigned long int fhlt;
	 int found15 = 0, found30=0, found50=0, found70=0, found100=0;

	 fhlt = hlNames_[i].find("HLT_Jet15U_v",0);
	 if(fhlt!=std::string::npos)
	   found15 = 1;
	 if(hlNames_[i]=="HLT_Jet15U" || found15)
	   is_HLT_Jet15U=HLTR->accept(i);

	 
	 fhlt = hlNames_[i].find("HLT_Jet30U_v",0);
	 if(fhlt!=std::string::npos)
	   found30 = 1;
	 if(hlNames_[i]=="HLT_Jet30U" || found30)
	   is_HLT_Jet30U=HLTR->accept(i);

	 fhlt = hlNames_[i].find("HLT_Jet50U_v",0);
         if(fhlt!=std::string::npos)
           found50 = 1;
	 if(hlNames_[i]=="HLT_Jet50U" || found50)
	   is_HLT_Jet50U=HLTR->accept(i);
	 
	 fhlt = hlNames_[i].find("HLT_Jet70U_v",0);
         if(fhlt!=std::string::npos)
           found70 = 1;
	 if(hlNames_[i]=="HLT_Jet70U" || found70)
	   is_HLT_Jet70U=HLTR->accept(i);
	 
	 fhlt = hlNames_[i].find("HLT_Jet100U_v",0);
         if(fhlt!=std::string::npos)
           found100 = 1;
	 if(hlNames_[i]=="HLT_Jet100U" || found100)
	   is_HLT_Jet100U=HLTR->accept(i);

	 if(hlNames_[i]=="HLT_Jet70U_v2")
	   is_HLT_Jet70U_v2=HLTR->accept(i);
	 
	 if(hlNames_[i]=="HLT_Jet100U_v2")
	   is_HLT_Jet100U_v2=HLTR->accept(i);

	 if(found15 || found30 || found50 || found70 || found100)
	   {
	     std::cout<<"SOMETHING CHANGED IN HLT MENU and that is "<<hlNames_[i]<<std::endl;
	     for(int j=0;j <nhlt_; j++)
	       {
		 int comp = strncmp(hltMatch_[j].c_str(),hlNames_[i].c_str(),10);
		 if(comp==0)
		   {
		     std::cout<<"name of HLT earlier = "<<hltMatch_[j]<<std::endl;
		     hltMatch_[j] = hlNames_[i];
		     std::cout<<"HLT NAME CHANGED to "<<hltMatch_[j]<<std::endl;
		     
		   }
	       }
	   }

       }//for(int i = 0;i<n; i++)
   

     //get prescales here 
     int iprescale = 0;
     if(nhlt_>100)
       std::cout<<"WARNING!Analyzer.cc can accomodate only atmost 100 triggers to match. not more than that. change the array size of hltprescale to accomodate that much"<<std::endl;
     
     for(std::vector<std::string>::const_iterator itr=hltMatch_.begin();itr!=hltMatch_.end();itr++)
       {
	 hltprescale[iprescale] =  hltConfig_.prescaleValue(iEvent, iSetup, (*itr)) ;
	 if(debugHLT_==1)
	   {
	     std::cout<<"Name of trigger to match = "<<(*itr)<<std::endl;
	     std::cout<<"event = "<<event<<" hltprescale["<<iprescale<<"] = "<<hltprescale[iprescale]<<std::endl;
	   }
	 
	 iprescale++;
       }
     
     if(debugHLT_==1)
       {
	 if(iprescale!=nhlt_)
	   std::cout<<"no. of triggers you WANT to match not equal to the ones PROVIDED(iprescale!=nhlt_)! correct it!"<<std::endl;
       }
     
     int itrigindex = 0;
     for(std::vector<std::string>::const_iterator itr=hltMatch_.begin();itr!=hltMatch_.end();itr++)
       {
	 triggerIndex[itrigindex] = triggerNames_.triggerIndex((*itr));  //for HLT-RECO object matching                                                               
	 itrigindex++;
       }//for(std::vector<std::string>::const
     
     ///HLT-reco object match  : taken from :
     ///afs/cern.ch/user/c/chenders/public/forMara/EcalOnlySumEtAnalyser.cc
     // get trigger event     
     //const int triggerIndex = triggerNames_.triggerIndex("HLT_Jet15U");
     Handle<trigger::TriggerEvent> triggerEventHandle;
     iEvent.getByLabel(triggerEventTag_,triggerEventHandle);
     if (!triggerEventHandle.isValid()) {
       std::cout << "Error in getting TriggerEvent product from Event!" << std::endl;
       return;
     }
   
   // A HLT path consists of many different modules - producers and filters     
   // The event can be rejected at any filter stage along the path     
   //
   // As well as just the basic pass/fail info,       
   // the TriggerResults object stores the index of the module in the path        
   // which made the final decision
   // In the case where the event was accepted, this index is therefore just    
   // the index of the last module along the path  

     for(int imatch=0 ;imatch<nhlt_; imatch++)
       {

            //std::cout<<"HLTR->size(): "<<(HLTR->size())<<std::endl;
            //std::cout<<"triggerIndex[Index]"<<(triggerIndex[imatch])<<std::endl;

	 if (  !(fabs(triggerIndex[imatch])<HLTR->size()) ) 
	   {
             //std::cout<<" Index = "<<imatch<<std::endl;
	     std::cout<<"!!!no trigger info found!!!"<<std::endl;
	     return;
	   }
	 
	 int lastModuleIndex = HLTR->index(triggerIndex[imatch]);
	 
	 if(debugHLT_ == 1)
	   {
	     std::cout<<"triggerIndex corresponding to"<<hltMatch_[imatch]<<std::endl; //debug
	     std::cout<<"lastModuleIndex = "<<lastModuleIndex<<std::endl; //debug
	   }
	 //   cout << "Last module index = " << lastModuleIndex <<endl;
	 // the HLT ConfigProvider can provide the module labels and types       
	 // (eg type HLTLevel1GTSeed with label hltL1sETT100)
	 // for the modules in a given path 
	 // it seems that what is stored is the module labels as strings in a vector
	 // then to get the module type, the ConfigProvider finds the label 
	 // inside the full menu that it has parsed and stored (processPset) ...          
	 // so basically get the module index or label first, then ask the ConfigProvider  
	 // what the type for this module label is... 
	 //   hltConfig_.moduleLabel(triggerIndex,moduleIndex);          
	 // hltConfig_.moduleLabels(triggerIndex) returns the vector of module labels       
	 // note that seems one can also access these by triggerName, not just triggerIndex ...     
	 // hltConfig_.moduleType(moduleLabel)          
	 // takes a string for module label and returns the module type (also a string)     
	 // loop over the modules in the desired path 
	 // but only the ones which were run for this event!             
	 // ie stop at lastModuleIndex from above                        
	 
	 indexhltobj = 0;
	 hltobjsize[imatch]  = 0;
	 
	 for(int imodule=0;imodule<=lastModuleIndex;imodule++) {
	   
	   if(debugHLT_ == 1)
	     std::cout << "imodule = "<<imodule<<std::endl;
	   
	   const std::string moduleLabel = hltConfig_.moduleLabel(triggerIndex[imatch],imodule);
	   const std::string moduleType = hltConfig_.moduleType(moduleLabel);
	   if(debugHLT_ == 1)
	     std::cout << " modulelabel = " << moduleLabel << " (type is " << moduleType << ")" <<std::endl;        
	   // once you have a moduleLabel, you can check if this module   
	   // wrote any filter products into the TriggerEvent  
	   // technically, you need to access this by an InputTag 
	   // that you build from the label + processName                                          
	   const int filterIndex = triggerEventHandle->filterIndex(InputTag(moduleLabel,"",processName_));
	   if(debugHLT_ == 1)
	     {
	       std::cout << "Filter index = " << filterIndex <<std::endl; 
	       std::cout<<"size of filter = "<<triggerEventHandle->sizeFilters()<<std::endl;
	     }
	   // this will return the index of the filter within the TriggerEvent's      
	   // list of all event filters,      
	   // or if not found, then it will just return the size of that filter list                            
	   // so that's how you check if it correctly found this filter module:               
	   
	   if(filterIndex<triggerEventHandle->sizeFilters()) {
	     // then this module was indeed found in the list of event filters       
	     if(debugHLT_ == 1)
	       {
		 std::cout << "Found this module in TriggerEvent filters!" <<std::endl;              
		 std::cout << moduleLabel << " (type " << moduleType << ")" <<std::endl;       
	       }
	     // so therefore it has written some TriggerObjects  
	     // now let's find them                 
	     // TriggerObject is a simple class of 4-vector info       
	     // plus an Id to represent physics type,   
	     // eg TriggerPhoton, TriggerMET, TriggerTET, etc                                                   
	     // as defined in DataFormats/HLTReco/interface/TriggerTypeDefs.h                            
	     // Note that a generic HLT filter can make a decision based on more 
	     // than just one object               
	     // hence a generic filter has a vector of associated trigger objects                
	     // the TriggerEvent has a 'master list' of all the *unique* trigger objects                    
	     // in this event        
	     // ie, even if a certain HLT object passes more than one trigger           
	     // (eg a single HLT muon which passes triggers of different pt thresholds,                    
	     // and/or is one of a combination which passes some dimuon or muon+X trigggers)      
	     // then the object still only gets written once to the master list                     
	     // For an individual filter then, its associated objects are stored 
	     // in the form of a set of indices which refer to where its specific objects           
	     // are in the master list            
	     // So first we get this set of keys for the objects for this filter                      
	     // note that 'Keys' is just a vector<int> in the trigger typedefs                     
	     const trigger::Keys& keys = triggerEventHandle->filterKeys(filterIndex);
	     if(debugHLT_ == 1)
	       std::cout << "n keys = " << keys.size() <<std::endl;                             
	     // and get the master list of all Trigger Objects
	     // ( probably I can get this outside the module loop            
	     // because there is just one list per event                 
	     // but I dont think this really matters ...)           
	     const trigger::TriggerObjectCollection& allTriggerObjects = triggerEventHandle->getObjects();
	     // now loop over the filter list of keys,    
	     // and get the appropriate trigger object from the master list                      
	     
	     for(trigger::size_type ikey=0;ikey<keys.size();ikey++) {
	       const trigger::TriggerObject& trigObject = allTriggerObjects[keys[ikey]];
	       //if(moduleType=="hlt1jet15U")
	       if(moduleType==hltTrigModule_[imatch])
		 {
		   if(debugHLT_ == 1)
		     std::cout<<"moduleType "<<moduleType<<" MATCHED with hltTrigModule_ = "<<hltTrigModule_[imatch]<<std::endl;
		   
		   ohtrigpt[imatch][indexhltobj]  = trigObject.pt();
		   ohtrigeta[imatch][indexhltobj] = trigObject.eta();
		   ohtrigphi[imatch][indexhltobj] = correct_phi(trigObject.phi());
		   
		   if(debugHLT_ == 1)
		     {
		       std::cout << "Trig object: type " << trigObject.id()<< "; pt,eta,phi = " << trigObject.pt() << ", " << trigObject.eta() << ", " << trigObject.phi()  <<std::endl;
		       std::cout<<"ohtrigpt["<<imatch<<"]["<<indexhltobj<<"] = "<<ohtrigpt[imatch][indexhltobj]<<std::endl;
		       std::cout<<"ohtrigeta["<<imatch<<"]["<<indexhltobj<<"] = "<<ohtrigeta[imatch][indexhltobj]<<std::endl;
		       std::cout<<"ohtrigphi["<<imatch<<"]["<<indexhltobj<<"] = "<<ohtrigphi[imatch][indexhltobj]<<std::endl;
		     }
		   
		   //DeltaR between HLT object and highest pt photon
		   /*if( runPatphotons_ && patphosize>0 )
		     {
		       //ohltphodR[imatch][indexhltobj] = deltaR(patphoeta[0],patphophi[0],ohtrigeta[imatch][indexhltobj],ohtrigphi[imatch][indexhltobj]);   
		       if(debugHLT_ == 1)
			 {
			   std::cout<<"deltaR between highest pt photon & HLT trigger object = "<<ohltphodR[imatch][indexhltobj]<<std::endl;
			 }
		       //now do this for every photon and select the mimimum dR matched HLT object 
		     }
		   if( runPatelectrons_ && patelesize>0)
		     {
		       ohltele1dR[imatch][indexhltobj] = deltaR(pateleeta[0],patelephi[0],ohtrigeta[imatch][indexhltobj],ohtrigphi[imatch][indexhltobj]);
		       if(patelesize>1)
			 ohltele2dR[imatch][indexhltobj] = deltaR(pateleeta[1],patelephi[1],ohtrigeta[imatch][indexhltobj],ohtrigphi[imatch][indexhltobj]);
		       
		       if(debugHLT_ == 1)
			 {
			   std::cout<<"deltaR between highest pt electron & HLT trigger object = "<<ohltele1dR[imatch][indexhltobj]<<std::endl;
			   if(patelesize>1)
			     std::cout<<"deltaR between sec highest pt electron & HLT trigger object = "<<ohltele2dR[imatch][indexhltobj]<<std::endl;
			 }
		     }//if( runPatelectrons_ && patelesize>0)
		   */

		   if( runpatJets_&&patjetsize>0 )
		     {
		       ohltjetdR[imatch][indexhltobj] = deltaR(patjeteta[0],patjetphi[0],ohtrigeta[imatch][indexhltobj],ohtrigphi[imatch][indexhltobj]);
		       if(debugHLT_ == 1)
			 {
			   std::cout<<"deltaR between highest pt jet & HLT trigger object = "<<ohltjetdR[imatch][indexhltobj]<<std::endl;
			 }
		     }
		   
		   
		   indexhltobj++;
		 }//	 if(moduleType=="hlt1jet15U")
	     }//for(trigger::size_type ikey=0;ikey<keys.size();ikey++)
	   }//if(filterIndex<triggerEventHandle->sizeFilters())
	 }//for(int imodule=0;imodule<=lastModuleIndex;imodule++) 
       
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ////NOTE : now deltaphi is calculated(just here) bet highest pt photon and the trigger object passing HLT_Jet15U///// 
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ///////////end of HLT-reco object match  

	 hltobjsize[imatch] = indexhltobj;

	 /*//do this to reject the photons for datadriven study
	 if( runPatphotons_ && patphosize>0 )
	   {
	     for(int ii=0; ii<patphosize; ii++)
	       {
		 double mindR = 0;
		 double dR = 0;
		 for(int jj=0; jj<hltobjsize[imatch]; jj++)
		   {
		     dR = deltaR(patphoeta[ii],patphophi[ii],ohtrigeta[imatch][jj],ohtrigphi[imatch][jj]);
		     if(mindR < dR)
		       {
			 mindR = dR;
		       }//if(mindR < dR)
		     patphohltmindRmatch[imatch][ii] = mindR;
		   }//for(jj=0; jj<hltobjsize[imatch]; jj++)
	       }//for(int ii=0; ii<patphosize; ii++)
	   }//if( runPatphotons_ && patphosize>0 )
	 */

       }//for(int imatch=0 ;imatch<nhlt_; imatch++)
     //sort the array of dR
   //double index[100] = {0.0};
   //TMath::Sort(hltobjsize,ohltphodR,index);
   }//if(runHLT_ == 1)
   

   /////////////////////////////CALLING UCSD's METHODS
   
   ////////CALL UCSD's HLT method
   
   if (runUCSDHLT_)
     hlt->analyze(iEvent, iSetup);
   
   ////////CALL UCSD's Vertex method
   
   if(runVertex_)
     vertex_std->analyze(iEvent, iSetup);


   //SC's
   if(runSC_)
     {
       //barrel sc's
       edm::Handle<reco::SuperClusterCollection> superClustersHybridH;
       iEvent.getByLabel(hybridSuperClusterColl_,superClustersHybridH);
       std::vector<reco::SuperCluster>scEB_container;

       //EE sc's
       edm::Handle<reco::SuperClusterCollection> superClustersEndcapH;
       iEvent.getByLabel(endcapSuperClusterColl_, superClustersEndcapH);
       std::vector<reco::SuperCluster>scEE_container;
       
       //std::vector<reco::SuperCluster>sc_container;

       //also get handle to HCAL for calculating H/E ----> reference http://cmslxr.fnal.gov/lxr/source/RecoEgamma/EgammaPhotonProducers/src/PhotonProducer.cc
       edm::Handle<CaloTowerCollection>  hcalTowersHandle;
       try
	 {
	   iEvent.getByLabel("towerMaker",hcalTowersHandle);
	   //std::cout<<"Found the handle"<<std::endl; 
	 }
       catch(...){std::cout<<"Didn't fine the HCAL::tower  handle"<<std::endl;}

       if( !hcalTowersHandle.isValid() )
	 std::cout<<"hcaltower handle is not valid"<<std::endl;

       
       // calculate HoE
       const CaloTowerCollection* hcalTowersColl = hcalTowersHandle.product();
       EgammaTowerIsolation towerIso1(hOverEConeSizeSC_,0.,0.,1,hcalTowersColl);  
       EgammaTowerIsolation towerIso2(hOverEConeSizeSC_,0.,0.,2,hcalTowersColl);  

       
       
       
       //////here info from EB and EE sc's are going into the same array. 
       
       //EB
       scp4->Clear();
       scxyz->Clear();
       
       int isc = 0;
       const reco::SuperClusterCollection* scEB = superClustersHybridH.product();
       for(reco::SuperClusterCollection::const_iterator sc = scEB->begin(); sc!=scEB->end(); ++sc)
	 {
	   scisEB[isc]             = 1;
	   scisEE[isc]             = 0;
	   
	   double scenergy           = sc->energy();
	   double sceta              = sc->eta();
	   //debug 
	   //std::cout<<"isc : inside EB; sceta = "<<isc << " : "<<sceta<<std::endl;
	   
	   double sctheta            = (2*atan(exp(-sc->eta())));
	   double scphi              = correct_phi(sc->phi());
	   double scpx               = scenergy*sin(sctheta)*cos(scphi);
	   double scpy               = scenergy*sin(sctheta)*sin(scphi);
	   double scpz               = scenergy*cos(sctheta);
	   
	   new ((*scp4)[isc]) TLorentzVector();
	   ((TLorentzVector *)scp4->At(isc))->SetXYZT(scpx,scpy,scpz,scenergy);
	   
	   
	   /*scx[isc]                = sc->x();
	     scy[isc]                = sc->y();
	     scz[isc]                = sc->z();
	   */
	   
	   new ((*scxyz)[isc]) TVector3();
	   ((TVector3 *)scxyz->At(isc))->SetXYZ(sc->x(),sc->y(),sc->z());
	   
	   scrawEnergy[isc]        = sc->rawEnergy();
	   scpreshowerEnergy[isc]  = sc->preshowerEnergy();
	   scclustersSize[isc]     = sc->clustersSize();
	   
	   
	   // cluster shape                                                                                                                                    
	   const reco::CaloClusterPtr seed = sc->seed();
	   sceMax[isc]             = lazyTool.eMax( *seed );
	   sce2x2[isc]             = lazyTool.e2x2( *seed );
	   sce3x3[isc]             = lazyTool.e3x3( *seed );
	   sce5x5[isc]             = lazyTool.e5x5( *seed );
	   scr4[isc]               = sce2x2[isc]/scenergy;
	   scr9[isc]               = sce3x3[isc]/scenergy;
	   scr25[isc]              = sce5x5[isc]/scenergy;
	   sce1bye4[isc]           = sceMax[isc]/sce2x2[isc];
	   sce1bye9[isc]           = sceMax[isc]/sce3x3[isc];
	   sce1bye25[isc]          = sceMax[isc]/sce5x5[isc];
	   
	   std::vector <float> locCov = lazyTool.localCovariances( *seed );
	   scetaWidth[isc]         = sqrt(locCov[0]);
	   scphiWidth[isc]         = sqrt(locCov[2]);
	   
	   const DetId seedId = lazyTool.getMaximum(*seed).first;
	   sce2e9[isc]             = mye2overe9(seedId,*barrelRecHits);
	   
	   //isolation
	   if( hcalTowersHandle.isValid() )
	     {
	       scHoE1[isc]             = towerIso1.getTowerESum(&(*sc))/sc->energy();
	       scHoE2[isc]             = towerIso2.getTowerESum(&(*sc))/sc->energy(); 
	       scHoE[isc]              = scHoE1[isc] + scHoE2[isc];
	     }
	   
	   else 
	     {
	       scHoE1[isc]             = -999.0;
	       scHoE2[isc]             = -999.0;
	       scHoE[isc]              = -999.0;
	     }

	   isc++;
	 }//for(reco::SuperClusterCollection::const_iterator sc = sc->begin(); sc!=sc->end(); ++sc)

       scEBsize = isc;

       //EE
       const reco::SuperClusterCollection* scEE = superClustersEndcapH.product();
       for(reco::SuperClusterCollection::const_iterator sc = scEE->begin(); sc!=scEE->end(); ++sc)
         {
	   scisEB[isc]             = 0;
	   scisEE[isc]             = 1;
	   
	   double scenergy           = sc->energy();
	   double sceta              = sc->eta();
	   //debug 
	   //std::cout<<"isc : inside EE; sceta = "<<isc << " : "<<sceta<<std::endl;
	   
	   double sctheta            = (2*atan(exp(-sc->eta())));
	   double scphi              = correct_phi(sc->phi());
	   double scpx               = scenergy*sin(sctheta)*cos(scphi);
	   double scpy               = scenergy*sin(sctheta)*sin(scphi);
	   double scpz               = scenergy*cos(sctheta);
	   
	   new ((*scp4)[isc]) TLorentzVector();
	   ((TLorentzVector *)scp4->At(isc))->SetXYZT(scpx,scpy,scpz,scenergy);
	   
	   
	   new ((*scxyz)[isc]) TVector3();
	   ((TVector3 *)scxyz->At(isc))->SetXYZ(sc->x(),sc->y(),sc->z());
	   
	   scrawEnergy[isc]        = sc->rawEnergy();
	   scpreshowerEnergy[isc]  = sc->preshowerEnergy();
	   scclustersSize[isc]     = sc->clustersSize();
	   
	   
           
	   // cluster shape                                                                                                                                         
	   const reco::CaloClusterPtr seed = sc->seed();
	   sceMax[isc]             = lazyTool.eMax( *seed );
	   sce2x2[isc]             = lazyTool.e2x2( *seed );
	   sce3x3[isc]             = lazyTool.e3x3( *seed );
	   sce5x5[isc]             = lazyTool.e5x5( *seed );
	   scr4[isc]               = sce2x2[isc]/scenergy;
	   scr9[isc]               = sce3x3[isc]/scenergy;
	   scr25[isc]              = sce5x5[isc]/scenergy;
	   sce1bye4[isc]           = sceMax[isc]/sce2x2[isc];
	   sce1bye9[isc]           = sceMax[isc]/sce3x3[isc];
	   sce1bye25[isc]          = sceMax[isc]/sce5x5[isc];
	   
	   std::vector <float> locCov = lazyTool.localCovariances( *seed );
	   scetaWidth[isc]         = sqrt(locCov[0]);
	   scphiWidth[isc]         = sqrt(locCov[2]);
	   sce2e9[isc]             = 0;
	   
	   //isolation
	   if( hcalTowersHandle.isValid() )
	     {
	       scHoE1[isc]             = towerIso1.getTowerESum(&(*sc))/sc->energy();
	       scHoE2[isc]             = towerIso2.getTowerESum(&(*sc))/sc->energy();
	       scHoE[isc]             = scHoE1[isc] + scHoE2[isc];
	     }
	   
	   else
	     {
	       scHoE1[isc]             = -999.0;
	       scHoE2[isc]             = -999.0;
	       scHoE[isc]              = -999.0;
	     }
	   isc++;
	 }//for(reco::SuperClusterCollection::const_iterator sc = scEE->begin(); sc!=scEE->end(); ++sc)

       scEEsize = isc-scEBsize;
       
       scsize = scEBsize + scEEsize;
       
     }//if(runSC_)  



   //geninfo
   if( rungenParticle_ && !isData )
     {
       Handle<reco::GenParticleCollection> genParticles;
       try{ iEvent.getByLabel("genParticles", genParticles); } catch(...) {;}
       
       genp4->Clear();
       genvtx->Clear();
       
       int igen = 0;
       
       for (reco::GenParticleCollection::const_iterator gen = genParticles->begin(); gen!= genParticles->end(); gen++)
	 {
	   genpdgid[igen]               = gen->pdgId();
	   genstatus[igen]              = gen->status();
	   
	   new ((*genp4)[igen]) TLorentzVector();
	   ((TLorentzVector *)genp4->At(igen))->SetPtEtaPhiM(gen->pt(), gen->eta(), gen->phi(), gen->mass());
	   
	   new ((*genvtx)[igen]) TVector3();
	   ((TVector3 *)genvtx->At(igen))->SetXYZ(gen->vx(), gen->vy(), gen->vz());
	   
	   gencharge[igen]           = gen->charge();
	
	   if (gen->numberOfMothers() != 0)
	     genmother[igen] = gen->motherRef().key();
	   else
	     genmother[igen] = -1;
   
	   genndau[igen] = gen->numberOfDaughters(); 
	   gennmoth[igen]= gen->numberOfMothers();

	   igen++;
	 }//for (GenParticleCollection::const_iterator p = genParticles->begin(); p!= genParticles->end(); p++)
       gensize = igen;
     }//if( rungenParticle_ )

   //gen jets
   if(rungenJets_)
     {
       Handle<reco::GenJetCollection> genJets;
       try{ iEvent.getByLabel( genJetLabel_, genJets ); }
       catch(...) {;}

       std::vector<reco::GenJet>genJet_container;
	      
       for(GenJetCollection::const_iterator jet = genJets->begin(); jet != genJets->end() ; jet++ ) 
	 {
	   genJet_container.push_back(*jet);
	   
	 }//for(GenJetCollection::const_iterator gen = genJets->begin(); gen != genJets->end() ; ++ gen )

       if(genJet_container.size()!=0)
	 {
	   for(unsigned int ijet=0; ijet<genJet_container.size(); ijet++)
	     {
	       genjetpt[ijet]                = genJet_container[ijet].pt();
	       genjetpx[ijet]                = genJet_container[ijet].px();
	       genjetpy[ijet]                = genJet_container[ijet].py();
	       genjetpz[ijet]                = genJet_container[ijet].pz();
	       genjetenergy[ijet]            = genJet_container[ijet].energy();
	       genjetphi[ijet]               = genJet_container[ijet].phi();
	       genjeteta[ijet]               = genJet_container[ijet].eta();
	       genjetemEnergy[ijet]          = genJet_container[ijet].emEnergy();
	       genjethadEnergy[ijet]         = genJet_container[ijet].hadEnergy();
	       genjetinvisibleEnergy[ijet]   = genJet_container[ijet].invisibleEnergy();
	     }//for(int ijet=0; ijet<genJet_container.size(); ijet++)

	 }//if(genJet_container.size()!=0)
       genjetsize = genJet_container.size();
       
     }//if(rungenJets_)





   tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
Analyser::beginJob()
{
  std::cout<<"inside begin job"<<std::endl;

  if(wantLocalFile_)
    {
      //std::cout<<"inside rootFilename_"<<std::endl;                                                                                                                   
      lf     = new TFile(loutputFile_.c_str(), "RECREATE");
      tree   = new TTree("myEvent","a tree with histograms");
    }

  if(wantRFIOFile_)
    {
      //std::cout<<"inside rfiorootFilename_"<<std::endl;                                                                                                               
      rfiof = new TRFIOFile(rfoutputFile_.c_str(), "RECREATE");
      tree  = new TTree("myEvent","a tree with histograms");
    }

  ///pile up reweightinh
  std::vector< float > Wlumi ;
  //std::vector< float > PoissonIntDist ;
  
  Double_t PoissonIntDist[25] = {
    0.104109,
    0.0703573,
    0.0698445,
    0.0698254,
    0.0697054,
    0.0697907,
    0.0696751,
    0.0694486,
    0.0680332,
    0.0651044,
    0.0598036,
    0.0527395,
    0.0439513,
    0.0352202,
    0.0266714,
    0.019411,
    0.0133974,
    0.00898536,
    0.0057516,
    0.00351493,
    0.00212087,
    0.00122891,
    0.00070592,
    0.000384744,
    0.000219377
  };

  float WLumi_f[25] = {
    0.0857895,
    0.0627972,
    0.0666277,
    0.0673892,
    0.0681171,
    0.0675986,
    0.0679215,
    0.0687084,
    0.0690788,
    0.0671276,
    0.0632601,
    0.0564774,
    0.0488245,
    0.0394473,
    0.0311093,
    0.0231416,
    0.0164191,
    0.0113318,
    0.007602,
    0.00491717,
    0.00285379,
    0.00174977,
    0.000856368,
    0.000468688,
    0.000259225
  };


  std::vector< float > TrueDist2011;
  float TrueDist2011_f[25] = {
    0.0132558,
    0.0316993,
    0.0719455,
    0.115284,
    0.145239,
    0.152783,
    0.139182,
    0.112847,
    0.082904,
    0.055968,
    0.0351001,
    0.0206271,
    0.0114405,
    0.00602595,
    0.00303009,
    0.00146112,
    0.000678137,
    0.000303988,
    0.000132051,
    5.57086e-05,
    2.28897e-05,
    9.17508e-06,
    3.59522e-06,
    1.3797e-06,
    8.16915e-07
  };


  /*float TrueDist2011_f[25] = {
    0.019091,
    0.0293974,
    0.0667931,
    0.108859,
    0.139533,
    0.149342,
    0.138629,
    0.114582,
    0.0859364,
    0.059324,
    0.0381123,
    0.0229881,
    0.0131129,
    0.00711764,
    0.00369635,
    0.00184543,
    0.000889604,
    0.000415683,
    0.000188921,
    0.000146288,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0
  };
  */

  for( int i=0; i<25; ++i) {
    TrueDist2011.push_back(TrueDist2011_f[i]);
    //Wlumi.push_back(WLumi_f[i]);
    Wlumi.push_back(PoissonIntDist[i]);
  }

  LumiWeights_ = edm::LumiReWeighting(Wlumi, TrueDist2011);



  //event properties
  tree->Branch("run",&run,"run/I");
  tree->Branch("event",&event,"event/I");
  tree->Branch("orbit",&orbit,"orbit/I");
  tree->Branch("bx",&bx,"bx/I");
  tree->Branch("lumis",&lumis,"lumis/I");
  tree->Branch("isData",&isData,"isData/I");
  tree->Branch("rho",&rho,"rho/D");
  tree->Branch("sigma",&sigma,"sigma/D");

  /*tree->Branch("rho44",&rho44,"rho44/D");
  tree->Branch("sigma44",&sigma44,"sigma44/D");
  */
  
  //pile up info
  tree->Branch("pvi",&pvi,"pvi/I");
  tree->Branch("pileup_bunchXing",pileup_bunchXing,"pileup_bunchXing[pvi]/I");
  tree->Branch("pileup_nvtx",pileup_nvtx,"pileup_nvtx[pvi]/I");
  tree->Branch("pileup_Avgnvtx",&pileup_Avgnvtx,"pileup_Avgnvtx/D");
  tree->Branch("pileupWeight",&pileupWeight,"pileupWeight/D");
  tree->Branch("pileupEventWeight",&pileupEventWeight,"pileupEventWeight/D");
  //fall11
  tree->Branch("pileuptrue",&pileuptrue,"pileuptrue/D");

  if(runGsftracks_)
    {
      tree->Branch("gsfsize",&gsfsize,"gsfsize/I");
      tree->Branch("gsfpt",gsfpt,"gsfpt[gsfsize]/D");
      tree->Branch("gsfcharge",gsfcharge,"gsfcharge[gsfsize]/I");
      tree->Branch("gsfchi2",gsfchi2,"gsfchi2[gsfsize]/D");
      tree->Branch("gsfeta",gsfeta,"gsfeta[gsfsize]/D");
      tree->Branch("gsfnumberOfLostHits",gsfnumberOfLostHits,"gsfnumberOfLostHits[gsfsize]/I");
      tree->Branch("gsfnumberOfValidHits",gsfnumberOfValidHits,"gsfnumberOfValidHits[gsfsize]/I");
      tree->Branch("gsflost",gsflost,"gsflost[gsfsize]/I");
      tree->Branch("gsfd0",gsfd0,"gsfd0[gsfsize]/D");
      tree->Branch("gsfdxy",gsfdxy,"gsfdxy[gsfsize]/D");
      tree->Branch("gsfdz",gsfdz,"gsfdz[gsfsize]/D");
      tree->Branch("gsfptin",gsfptin,"gsfptin[gsfsize]/D");
      tree->Branch("gsfptout",gsfptout,"gsfptout[gsfsize]/D");
      tree->Branch("gsffbrem",gsffbrem,"gsffbrem[gsfsize]/D");
      tree->Branch("gsfqoverp",gsfqoverp,"gsfqoverp[gsfsize]/D");
      tree->Branch("gsfvx",gsfvx,"gsfvx[gsfsize]/D");
      tree->Branch("gsfvy",gsfvy,"gsfvy[gsfsize]/D");
      tree->Branch("gsfvz",gsfvz,"gsfvz[gsfsize]/D");
      tree->Branch("gsfphi",gsfphi,"gsfphi[gsfsize]/D");
      tree->Branch("gsfndof",gsfndof,"gsfndof[gsfsize]/D");
      tree->Branch("gsfrecHitsSize",gsfrecHitsSize,"gsfrecHitsSize[gsfsize]/I");
      tree->Branch("gsftheta",gsftheta,"gsftheta[gsfsize]/D");
      tree->Branch("gsfqualityMask",gsfqualityMask,"gsfqualityMask[gsfsize]/I");
      tree->Branch("gsfouterX",gsfouterX,"gsfouterX[gsfsize]/D");
      tree->Branch("gsfouterY",gsfouterY,"gsfouterY[gsfsize]/D");
      tree->Branch("gsfouterZ",gsfouterZ,"gsfouterZ[gsfsize]/D");
      tree->Branch("gsfouterRadius",gsfouterRadius,"gsfouterRadius[gsfsize]/D");
      tree->Branch("gsfinnerX",gsfinnerX,"gsfinnerX[gsfsize]/D");
      tree->Branch("gsfinnerY",gsfinnerY,"gsfinnerY[gsfsize]/D");
      tree->Branch("gsfinnerZ",gsfinnerZ,"gsfinnerZ[gsfsize]/D");
    }//if(runGsftracks_)

  //electron
  if(runGsfelectrons_)
    {
      tree->Branch("elesize",&elesize,"elesize/I");
      tree->Branch("elept",elept,"elept[elesize]/D");
      tree->Branch("elepx",elepx,"elepx[elesize]/D");
      tree->Branch("elepy",elepy,"elepy[elesize]/D");
      tree->Branch("elepz",elepz,"elepz[elesize]/D");
      tree->Branch("elephi",elephi,"elephi[elesize]/D");
      tree->Branch("eletheta",eletheta,"eletheta[elesize]/D");
      tree->Branch("elefbrem",elefbrem,"elefbrem[elesize]/D");
      tree->Branch("eletrackerDrivenSeed",eletrackerDrivenSeed,"eletrackerDrivenSeed[elesize]/I");
      tree->Branch("eleecalDrivenSeed",eleecalDrivenSeed,"eleecalDrivenSeed[elesize]/I");
      tree->Branch("elecharge",elecharge,"elecharge[elesize]/I");
      tree->Branch("elesceta",elesceta,"elesceta[elesize]/D");
      tree->Branch("eledr03EcalRecHitSumEt",eledr03EcalRecHitSumEt,"eledr03EcalRecHitSumEt[elesize]/D");
      tree->Branch("eledr03HcalDepth1TowerSumEt",eledr03HcalDepth1TowerSumEt,"eledr03HcalDepth1TowerSumEt[elesize]/D");
      tree->Branch("eledr03HcalDepth2TowerSumEt",eledr03HcalDepth2TowerSumEt,"eledr03HcalDepth2TowerSumEt[elesize]/D");
      tree->Branch("eledr03HcalTowerSumEt",eledr03HcalTowerSumEt,"eledr03HcalTowerSumEt[elesize]/D");
      tree->Branch("eledr04EcalRecHitSumEt",eledr04EcalRecHitSumEt,"eledr04EcalRecHitSumEt[elesize]/D");
      tree->Branch("eledr04HcalDepth1TowerSumEt",eledr04HcalDepth1TowerSumEt,"eledr04HcalDepth1TowerSumEt[elesize]/D");
      tree->Branch("eledr04HcalDepth2TowerSumEt",eledr04HcalDepth2TowerSumEt,"eledr04HcalDepth2TowerSumEt[elesize]/D");
      tree->Branch("eledr04HcalTowerSumEt",eledr04HcalTowerSumEt,"eledr04HcalTowerSumEt[elesize]/D");
      tree->Branch("elee1x5",elee1x5,"elee1x5[elesize]/D");
      tree->Branch("elee2x5Max",elee2x5Max,"elee2x5Max[elesize]/D");
      tree->Branch("elee5x5",elee5x5,"elee5x5[elesize]/D");
      tree->Branch("eleeEleClusterOverPout",eleeEleClusterOverPout,"eleeEleClusterOverPout[elesize]/D");
      tree->Branch("eleeSeedClusterOverP",eleeSeedClusterOverP,"eleeSeedClusterOverP[elesize]/D");
      tree->Branch("eleeSeedClusterOverPout",eleeSeedClusterOverPout,"eleeSeedClusterOverPout[elesize]/D");
      tree->Branch("eleeSuperClusterOverP",eleeSuperClusterOverP,"eleeSuperClusterOverP[elesize]/D");
      tree->Branch("eleeta",eleeta,"eleeta[elesize]/D");
      tree->Branch("elehadronicOverEm",elehadronicOverEm,"elehadronicOverEm[elesize]/D");
      tree->Branch("elesigmaIetaIeta",elesigmaIetaIeta,"elesigmaIetaIeta[elesize]/D");
      tree->Branch("eledeltaEtaEleClusterTrackAtCalo",eledeltaEtaEleClusterTrackAtCalo,"eledeltaEtaEleClusterTrackAtCalo[elesize]/D");
      tree->Branch("eledeltaEtaSeedClusterTrackAtCalo",eledeltaEtaSeedClusterTrackAtCalo,"eledeltaEtaSeedClusterTrackAtCalo[elesize]/D");
      tree->Branch("eledeltaEtaSuperClusterTrackAtVtx",eledeltaEtaSuperClusterTrackAtVtx,"eledeltaEtaSuperClusterTrackAtVtx[elesize]/D");
      tree->Branch("eledeltaPhiEleClusterTrackAtCalo",eledeltaPhiEleClusterTrackAtCalo,"eledeltaPhiEleClusterTrackAtCalo[elesize]/D");
      tree->Branch("eledeltaPhiSeedClusterTrackAtCalo",eledeltaPhiSeedClusterTrackAtCalo,"eledeltaPhiSeedClusterTrackAtCalo[elesize]/D");
      tree->Branch("eledeltaPhiSuperClusterTrackAtVtx",eledeltaPhiSuperClusterTrackAtVtx,"eledeltaPhiSuperClusterTrackAtVtx[elesize]/D");
      tree->Branch("eleenergy",eleenergy,"eleenergy[elesize]/D");
      tree->Branch("elemva",elemva,"elemva[elesize]/D");
      tree->Branch("elenumberOfTracks",elenumberOfTracks,"elenumberOfTracks[elesize]/D");
      tree->Branch("eledr03TkSumPt",eledr03TkSumPt,"eledr03TkSumPt[elesize]/D");
      tree->Branch("eledr04TkSumPt",eledr04TkSumPt,"eledr04TkSumPt[elesize]/D");
      
      tree->Branch("elescphi",elescphi,"elescphi[elesize]/D");
      tree->Branch("elescE",elescE,"elescE[elesize]/D");
      //tree->Branch("elescpt",elescpt,"elescpt[elesize]/D");
  
      //PF info
      tree->Branch("elepfeta",elepfeta,"elepfeta[elesize]/D");
      tree->Branch("elepfphi",elepfphi,"elepfphi[elesize]/D");
      tree->Branch("elepfE",elepfE,"elepfE[elesize]/D");
      //tree->Branch("elepfpt",elepfpt,"elepfpt[elesize]/D");
      


      //spikes
      tree->Branch("elemaxEnergyXtal",elemaxEnergyXtal,"elemaxEnergyXtal[elesize]/D");
      tree->Branch("eleswissCross",eleswissCross,"eleswissCross[elesize]/D");
      tree->Branch("eleswissBasedspikevar",eleswissBasedspikevar,"eleswissBasedspikevar[elesize]/D");

      //seed timing
      tree->Branch("eleseedtime",eleseedtime,"eleseedtime[elesize]/D");
      tree->Branch("elerecoFlag",elerecoFlag,"elerecoFlag[elesize]/I");
      //tree->Branch("elecheckFlag",elecheckFlag,"elecheckFlag[elesize]/I");
      //tree->Branch("eleseverityLevel",eleseverityLevel,"eleseverityLevel[elesize]/D");
  
      //gsf
      tree->Branch("eletrkpt",eletrkpt,"eletrkpt[elesize]/D");
      tree->Branch("eletrkcharge",eletrkcharge,"eletrkcharge[elesize]/I");
      tree->Branch("eletrkchi2",eletrkchi2,"eletrkchi2[elesize]/D");
      tree->Branch("eletrketa",eletrketa,"eletrketa[elesize]/D");
      tree->Branch("eletrknumberOfLostHits",eletrknumberOfLostHits,"eletrknumberOfLostHits[elesize]/I");
      tree->Branch("eletrknumberOfValidHits",eletrknumberOfValidHits,"eletrknumberOfValidHits[elesize]/I");
      tree->Branch("eletrklost",eletrklost,"eletrklost[elesize]/I");
      tree->Branch("eletrkd0",eletrkd0,"eletrkd0[elesize]/D");
      tree->Branch("eletrkdxy",eletrkdxy,"eletrkdxy[elesize]/D");
      tree->Branch("eletrkdz",eletrkdz,"eletrkdz[elesize]/D");
      tree->Branch("eletrkptin",eletrkptin,"eletrkptin[elesize]/D");
      tree->Branch("eletrkptout",eletrkptout,"eletrkptout[elesize]/D");
      tree->Branch("eletrkfbrem",eletrkfbrem,"eletrkfbrem[elesize]/D");
      tree->Branch("eletrkqoverp",eletrkqoverp,"eletrkqoverp[elesize]/D");
      tree->Branch("eletrkvx",eletrkvx,"eletrkvx[elesize]/D");
      tree->Branch("eletrkvy",eletrkvy,"eletrkvy[elesize]/D");
      tree->Branch("eletrkvz",eletrkvz,"eletrkvz[elesize]/D");
      tree->Branch("eletrkphi",eletrkphi,"eletrkphi[elesize]/D");
      tree->Branch("eletrkndof",eletrkndof,"eletrkndof[elesize]/D");
      tree->Branch("eletrkrecHitsSize",eletrkrecHitsSize,"eletrkrecHitsSize[elesize]/I");
      tree->Branch("eletrktheta",eletrktheta,"eletrktheta[elesize]/D");
      tree->Branch("eletrkqualityMask",eletrkqualityMask,"eletrkqualityMask[elesize]/I");
      tree->Branch("eletrkouterX",eletrkouterX,"eletrkouterX[elesize]/D");
      tree->Branch("eletrkouterY",eletrkouterY,"eletrkouterY[elesize]/D");
      tree->Branch("eletrkouterZ",eletrkouterZ,"eletrkouterZ[elesize]/D");
      tree->Branch("eletrkouterRadius",eletrkouterRadius,"eletrkouterRadius[elesize]/D");
      tree->Branch("eletrkinnerX",eletrkinnerX,"eletrkinnerX[elesize]/D");
      tree->Branch("eletrkinnerY",eletrkinnerY,"eletrkinnerY[elesize]/D");
      tree->Branch("eletrkinnerZ",eletrkinnerZ,"eletrkinnerZ[elesize]/D");
      
      //ctf
      tree->Branch("electfpt",electfpt,"electfpt[elesize]/D");
      tree->Branch("electfcharge",electfcharge,"electfcharge[elesize]/I");
      tree->Branch("electfchi2",electfchi2,"electfchi2[elesize]/D");
      tree->Branch("electfeta",electfeta,"electfeta[elesize]/D");
      tree->Branch("electfnumberOfLostHits",electfnumberOfLostHits,"electfnumberOfLostHits[elesize]/I");
      tree->Branch("electfnumberOfValidHits",electfnumberOfValidHits,"electfnumberOfValidHits[elesize]/I");
      tree->Branch("electflost",electflost,"electflost[elesize]/I");
      tree->Branch("electfd0",electfd0,"electfd0[elesize]/D");
      tree->Branch("electfdxy",electfdxy,"electfdxy[elesize]/D");
      tree->Branch("electfdz",electfdz,"electfdz[elesize]/D");
      tree->Branch("electfptin",electfptin,"electfptin[elesize]/D");
      tree->Branch("electfptout",electfptout,"electfptout[elesize]/D");
      tree->Branch("electffbrem",electffbrem,"electffbrem[elesize]/D");
      tree->Branch("electfqoverp",electfqoverp,"electfqoverp[elesize]/D");
      tree->Branch("electfvx",electfvx,"electfvx[elesize]/D");
      tree->Branch("electfvy",electfvy,"electfvy[elesize]/D");
      tree->Branch("electfvz",electfvz,"electfvz[elesize]/D");
      tree->Branch("electfphi",electfphi,"electfphi[elesize]/D");
      tree->Branch("electfndof",electfndof,"electfndof[elesize]/D");
      tree->Branch("electfrecHitsSize",electfrecHitsSize,"electfrecHitsSize[elesize]/I");
      tree->Branch("electftheta",electftheta,"electftheta[elesize]/D");
      tree->Branch("electfqualityMask",electfqualityMask,"electfqualityMask[elesize]/I");
      tree->Branch("electfouterX",electfouterX,"electfouterX[elesize]/D");
      tree->Branch("electfouterY",electfouterY,"electfouterY[elesize]/D");
      tree->Branch("electfouterZ",electfouterZ,"electfouterZ[elesize]/D");
      tree->Branch("electfouterRadius",electfouterRadius,"electfouterRadius[elesize]/D");
      tree->Branch("electfinnerX",electfinnerX,"electfinnerX[elesize]/D");
      tree->Branch("electfinnerY",electfinnerY,"electfinnerY[elesize]/D");
      tree->Branch("electfinnerZ",electfinnerZ,"electfinnerZ[elesize]/D");
      
    }//if(runGsfelectrons_)

  //general track
  if(runGeneraltracks_)
    {
      gentrkp4 = new TClonesArray("TLorentzVector", MAX_TRACKS);
      gentrk_vtxpos = new TClonesArray("TVector3", MAX_TRACKS);
      
      tree->Branch("gentrksize",&gentrksize,"gentrksize/I");
      tree->Branch("gentrkqualityMask",gentrkqualityMask,"gentrkqualityMask[gentrksize]/I");

      if(usefultrackTreeVar_){
	
	tree->Branch("gentrkp4", "TClonesArray", &gentrkp4, 32000, 0);
	tree->Branch("gentrk_vtxpos", "TClonesArray", &gentrk_vtxpos, 32000, 0);
	
	tree->Branch("gentrkcharge",gentrkcharge,"gentrkcharge[gentrksize]/I");
	tree->Branch("gentrkchi2",gentrkchi2,"gentrkchi2[gentrksize]/D");
	tree->Branch("gentrknumberOfLostHits",gentrknumberOfLostHits,"gentrknumberOfLostHits[gentrksize]/I");
	tree->Branch("gentrknumberOfValidHits",gentrknumberOfValidHits,"gentrknumberOfValidHits[gentrksize]/I");
	tree->Branch("gentrklost",gentrklost,"gentrklost[gentrksize]/I");
	tree->Branch("gentrkd0",gentrkd0,"gentrkd0[gentrksize]/D");
	tree->Branch("gentrkdxy",gentrkdxy,"gentrkdxy[gentrksize]/D");
	tree->Branch("gentrkdz",gentrkdz,"gentrkdz[gentrksize]/D");
	tree->Branch("gentrkptin",gentrkptin,"gentrkptin[gentrksize]/D");
	tree->Branch("gentrkptout",gentrkptout,"gentrkptout[gentrksize]/D");
	tree->Branch("gentrkfbrem",gentrkfbrem,"gentrkfbrem[gentrksize]/D");
	tree->Branch("gentrkqoverp",gentrkqoverp,"gentrkqoverp[gentrksize]/D");
	tree->Branch("gentrkndof",gentrkndof,"gentrkndof[gentrksize]/D");
	tree->Branch("gentrkrecHitsSize",gentrkrecHitsSize,"gentrkrecHitsSize[gentrksize]/I");
	
	tree->Branch("gentrkouterX",gentrkouterX,"gentrkouterX[gentrksize]/D");
	tree->Branch("gentrkouterY",gentrkouterY,"gentrkouterY[gentrksize]/D");
	tree->Branch("gentrkouterZ",gentrkouterZ,"gentrkouterZ[gentrksize]/D");
	tree->Branch("gentrkouterRadius",gentrkouterRadius,"gentrkouterRadius[gentrksize]/D");
	tree->Branch("gentrkinnerX",gentrkinnerX,"gentrkinnerX[gentrksize]/D");
	tree->Branch("gentrkinnerY",gentrkinnerY,"gentrkinnerY[gentrksize]/D");
	tree->Branch("gentrkinnerZ",gentrkinnerZ,"gentrkinnerZ[gentrksize]/D");
	tree->Branch("gentrk_hpnvalid", gentrk_hpnvalid, "gentrk_hpnvalid[gentrksize]/I");
	tree->Branch("gentrk_hpnlost", gentrk_hpnlost, "gentrk_hpnlost[gentrksize]/I");
	tree->Branch("gentrk_hpnvalidpix", gentrk_hpnvalidpix, "gentrk_hpnvalidpix[gentrksize]/I");
      }//if(usefultrackTreeVar_)

    }//if(runGeneraltracks_)

  //ES rechits
  if(runPreshower_)
    {
      tree->Branch("EEscsize",&EEscsize,"EEscsize/I");
      tree->Branch("EEsceta",EEsceta,"EEsceta[EEscsize]/D");
      tree->Branch("EEscet",EEscet,"EEscet[EEscsize]/D");
      
      //ES info
      tree->Branch("ESclussizeX",ESclussizeX,"ESclussizeX[EEscsize]/I");
      tree->Branch("ESclusEtX",ESclusEtX,"ESclusEtX[EEscsize][1000]/D");
      tree->Branch("ESclusEX",ESclusEX,"ESclusEX[EEscsize][1000]/D");
      tree->Branch("ESclusetaX",ESclusetaX,"ESclusetaX[EEscsize][1000]/D");
      tree->Branch("ESclusphiX",ESclusphiX,"ESclusphiX[EEscsize][1000]/D");
      tree->Branch("ESclusnhitsX",ESclusnhitsX,"ESclusnhitsX[EEscsize][1000]/D");
      
      tree->Branch("ESclussizeY",ESclussizeY,"ESclussizeY[EEscsize]/I");
      tree->Branch("ESclusEtY",ESclusEtY,"ESclusEtY[EEscsize][1000]/D");
      tree->Branch("ESclusEY",ESclusEY,"ESclusEY[EEscsize][1000]/D");
      tree->Branch("ESclusetaY",ESclusetaY,"ESclusetaY[EEscsize][1000]/D");
      tree->Branch("ESclusphiY",ESclusphiY,"ESclusphiY[EEscsize][1000]/D");
      tree->Branch("ESclusnhitsY",ESclusnhitsY,"ESclusnhitsY[EEscsize][1000]/D");
      
      
      
      tree->Branch("ESrechitsize",ESrechitsize,"ESrechitsize[EEscsize]/I");
      tree->Branch("ESrechitE",ESrechitE,"ESrechitE[EEscsize][10000]/D");
      tree->Branch("ESrechitplane",ESrechitplane,"ESrechitplane[EEscsize][10000]/D");
      tree->Branch("ESrechitzside",ESrechitzside,"ESrechitzside[EEscsize][10000]/D");
      tree->Branch("ESrechitsix",ESrechitsix,"ESrechitsix[EEscsize][10000]/D");
      tree->Branch("ESrechitsiy",ESrechitsiy,"ESrechitsiy[EEscsize][10000]/D");
    }//if(runPreshower_)

  //pat electrons
  if(runPatelectrons_)
    {
      patelsc = new TClonesArray("TLorentzVector", MAX_ELECTRONS);
      patelp4 = new TClonesArray("TLorentzVector", MAX_ELECTRONS);
      patelmomvtx = new TClonesArray("TVector3", MAX_ELECTRONS);
      patelmomvtxconst = new TClonesArray("TVector3", MAX_ELECTRONS);
      patelmomcalo = new TClonesArray("TVector3", MAX_ELECTRONS);
      patelmomout = new TClonesArray("TVector3", MAX_ELECTRONS);
      patelposvtx = new TClonesArray("TVector3", MAX_ELECTRONS);
      patelposcalo = new TClonesArray("TVector3", MAX_ELECTRONS);
      
      
      
      tree->Branch("patelesize",&patelesize,"patelesize/I");
      tree->Branch("patelsc", "TClonesArray", &patelsc, 32000, 0);
      tree->Branch("patelp4", "TClonesArray", &patelp4, 32000, 0);
      tree->Branch("patelescind",patelescind,"patelescind[patelesize]/I");
      tree->Branch("patelefbrem",patelefbrem,"patelefbrem[patelesize]/D");
      tree->Branch("pateletrackerDrivenSeed",pateletrackerDrivenSeed,"pateletrackerDrivenSeed[patelesize]/I");
      tree->Branch("pateleecalDrivenSeed",pateleecalDrivenSeed,"pateleecalDrivenSeed[patelesize]/I");
      tree->Branch("patelecharge",patelecharge,"patelecharge[patelesize]/I");
      tree->Branch("pateledr03EcalRecHitSumEt",pateledr03EcalRecHitSumEt,"pateledr03EcalRecHitSumEt[patelesize]/D");
      tree->Branch("pateledr03HcalDepth1TowerSumEt",pateledr03HcalDepth1TowerSumEt,"pateledr03HcalDepth1TowerSumEt[patelesize]/D");
      tree->Branch("pateledr03HcalDepth2TowerSumEt",pateledr03HcalDepth2TowerSumEt,"pateledr03HcalDepth2TowerSumEt[patelesize]/D");
      tree->Branch("pateledr03HcalTowerSumEt",pateledr03HcalTowerSumEt,"pateledr03HcalTowerSumEt[patelesize]/D");
      tree->Branch("pateledr04EcalRecHitSumEt",pateledr04EcalRecHitSumEt,"pateledr04EcalRecHitSumEt[patelesize]/D");
      tree->Branch("pateledr04HcalDepth1TowerSumEt",pateledr04HcalDepth1TowerSumEt,"pateledr04HcalDepth1TowerSumEt[patelesize]/D");
      tree->Branch("pateledr04HcalDepth2TowerSumEt",pateledr04HcalDepth2TowerSumEt,"pateledr04HcalDepth2TowerSumEt[patelesize]/D");
      tree->Branch("pateledr04HcalTowerSumEt",pateledr04HcalTowerSumEt,"pateledr04HcalTowerSumEt[patelesize]/D");
      tree->Branch("patelee1x5",patelee1x5,"patelee1x5[patelesize]/D");
      tree->Branch("patelee2x5Max",patelee2x5Max,"patelee2x5Max[patelesize]/D");
      tree->Branch("patelee5x5",patelee5x5,"patelee5x5[patelesize]/D");
      tree->Branch("patelehadronicOverEm",patelehadronicOverEm,"patelehadronicOverEm[patelesize]/D");
      tree->Branch("patelesigmaIetaIeta",patelesigmaIetaIeta,"patelesigmaIetaIeta[patelesize]/D");
      tree->Branch("pateledeltaEtaSuperClusterTrackAtVtx",pateledeltaEtaSuperClusterTrackAtVtx,"pateledeltaEtaSuperClusterTrackAtVtx[patelesize]/D");
      tree->Branch("pateledeltaPhiSuperClusterTrackAtVtx",pateledeltaPhiSuperClusterTrackAtVtx,"pateledeltaPhiSuperClusterTrackAtVtx[patelesize]/D");
      tree->Branch("patelenumberOfTracks",patelenumberOfTracks,"patelenumberOfTracks[patelesize]/I");
      tree->Branch("pateledr03TkSumPt",pateledr03TkSumPt,"pateledr03TkSumPt[patelesize]/D");
      tree->Branch("pateledr04TkSumPt",pateledr04TkSumPt,"pateledr04TkSumPt[patelesize]/D");
      //added on Sept 26th, 2011
      tree->Branch("pateleecalEnergy", pateleecalEnergy, "pateleecalEnergy[patelesize]/D");
      tree->Branch("patelecaloEnergy", patelecaloEnergy, "patelecaloEnergy[patelesize]/D");
      //tree->Branch("patelecaloCorrectedEnergy", patelecaloCorrectedEnergy, "patelecaloCorrectedEnergy[patelesize]/D");
      tree->Branch("patelecorrectedEcalEnergy", patelecorrectedEcalEnergy, "patelecorrectedEcalEnergy[patelesize]/D");
      //spikes
      tree->Branch("patelemaxEnergyXtal",patelemaxEnergyXtal,"patelemaxEnergyXtal[patelesize]/D");
      tree->Branch("pateleswissCross",pateleswissCross,"pateleswissCross[patelesize]/D");
      tree->Branch("pateleswissBasedspikevar",pateleswissBasedspikevar,"pateleswissBasedspikevar[patelesize]/D");
      
      //seed timing
      tree->Branch("pateleseedtime",pateleseedtime,"pateleseedtime[patelesize]/D");
      tree->Branch("patelerecoFlag",patelerecoFlag,"patelerecoFlag[patelesize]/I");
      //tree->Branch("patelecheckFlag",patelecheckFlag,"patelecheckFlag[patelesize]/I");
      //tree->Branch("patelekOutOfTime",patelekOutOfTime,"patelekOutOfTime[patelesize]/D");
      //tree->Branch("pateleseverityLevel",pateleseverityLevel,"pateleseverityLevel[patelesize]/D");
      tree->Branch("patelee2e9",patelee2e9,"patelee2e9[patelesize]/D");
      /////conversion rejection
      tree->Branch("pateleExpectednumberOfHits",pateleExpectednumberOfHits,"pateleExpectednumberOfHits[patelesize]/I");


      tree->Branch("pateleconvDist",pateleconvDist,"pateleconvDist[patelesize]/D");
      tree->Branch("pateleconvDcot",pateleconvDcot,"pateleconvDcot[patelesize]/D");
      tree->Branch("pateleconvRadius",pateleconvRadius,"pateleconvRadius[patelesize]/D");
      
      /////for dxy
      tree->Branch("pateletrkd0",pateletrkd0,"pateletrkd0[patelesize]/D");
      tree->Branch("pateletrkdxy",pateletrkdxy,"pateletrkdxy[patelesize]/D");
      tree->Branch("pateletrkdz",pateletrkdz,"pateletrkdz[patelesize]/D");
      
      tree->Branch("patelmomvtx", "TClonesArray", &patelmomvtx, 32000, 0);
      
      tree->Branch("patelposvtx", "TClonesArray", &patelposvtx, 32000, 0);
      
      ///1st Aug,2013
      tree->Branch("patelenumberOfLostHits",patelenumberOfLostHits,"patelenumberOfLostHits[patelesize]/I");
	

      if(usefulpateleTreeVar_){
	tree->Branch("patelmomvtxconst", "TClonesArray", &patelmomvtxconst, 32000, 0);
	tree->Branch("patelmomcalo", "TClonesArray", &patelmomcalo, 32000, 0);
	tree->Branch("patelmomout", "TClonesArray", &patelmomout, 32000, 0);
	
	tree->Branch("patelposcalo", "TClonesArray", &patelposcalo, 32000, 0);
	
	
	tree->Branch("pateleeEleClusterOverPout",pateleeEleClusterOverPout,"pateleeEleClusterOverPout[patelesize]/D");
	tree->Branch("pateleeSeedClusterOverP",pateleeSeedClusterOverP,"pateleeSeedClusterOverP[patelesize]/D");
	tree->Branch("pateleeSeedClusterOverPout",pateleeSeedClusterOverPout,"pateleeSeedClusterOverPout[patelesize]/D");
	tree->Branch("pateleeSuperClusterOverP",pateleeSuperClusterOverP,"pateleeSuperClusterOverP[patelesize]/D");
	
	tree->Branch("pateledeltaEtaEleClusterTrackAtCalo",pateledeltaEtaEleClusterTrackAtCalo,"pateledeltaEtaEleClusterTrackAtCalo[patelesize]/D");
	tree->Branch("pateledeltaEtaSeedClusterTrackAtCalo",pateledeltaEtaSeedClusterTrackAtCalo,"pateledeltaEtaSeedClusterTrackAtCalo[patelesize]/D");
	
	tree->Branch("pateledeltaPhiEleClusterTrackAtCalo",pateledeltaPhiEleClusterTrackAtCalo,"pateledeltaPhiEleClusterTrackAtCalo[patelesize]/D");
	tree->Branch("pateledeltaPhiSeedClusterTrackAtCalo",pateledeltaPhiSeedClusterTrackAtCalo,"pateledeltaPhiSeedClusterTrackAtCalo[patelesize]/D");
	
	tree->Branch("patelemva",patelemva,"patelemva[patelesize]/D");
	
	tree->Branch("pateletrknumberOfLostHits",pateletrknumberOfLostHits,"pateletrknumberOfLostHits[patelesize]/I");
	
	tree->Branch("patelepin", patelepin, "patelepin[patelesize]/D");
	tree->Branch("patelepout", patelepout, "patelepout[patelesize]/D");
	
	//added on 20th sept, 2011
	tree->Branch("pateleecalEnergyError", pateleecalEnergyError, "pateleecalEnergyError[patelesize]/D");
	tree->Branch("pateletrackMomentumError", pateletrackMomentumError, "pateletrackMomentumError[patelesize]/D");
	
	
	//PF info
	tree->Branch("patelepfeta",patelepfeta,"patelepfeta[patelesize]/D");
	tree->Branch("patelepfphi",patelepfphi,"patelepfphi[patelesize]/D");
	tree->Branch("patelepfE",patelepfE,"patelepfE[patelesize]/D");
	//tree->Branch("patelepfpt",patelepfpt,"patelepfpt[patelesize]/D");
	
	//gaps info
	tree->Branch("pateleisEB",pateleisEB,"pateleisEB[patelesize]/I");
	tree->Branch("pateleisEBEEGap",pateleisEBEEGap,"pateleisEBEEGap[patelesize]/I");
	tree->Branch("pateleisEBEtaGap",pateleisEBEtaGap,"pateleisEBEtaGap[patelesize]/I");
	tree->Branch("pateleisEBGap",pateleisEBGap,"pateleisEBGap[patelesize]/I");
	tree->Branch("pateleisEBPhiGap",pateleisEBPhiGap,"pateleisEBPhiGap[patelesize]/I");
	tree->Branch("pateleisEE",pateleisEE,"pateleisEE[patelesize]/I");
	tree->Branch("pateleisEEDeeGap",pateleisEEDeeGap,"pateleisEEDeeGap[patelesize]/I");
	tree->Branch("pateleisEEGap",pateleisEEGap,"pateleisEEGap[patelesize]/I");
	tree->Branch("pateleisEERingGap",pateleisEERingGap,"pateleisEERingGap[patelesize]/I");
	
	
	//gsf
	tree->Branch("pateletrkpt",pateletrkpt,"pateletrkpt[patelesize]/D");
	tree->Branch("pateletrkcharge",pateletrkcharge,"pateletrkcharge[patelesize]/I");
	tree->Branch("pateletrkchi2",pateletrkchi2,"pateletrkchi2[patelesize]/D");
	tree->Branch("pateletrketa",pateletrketa,"pateletrketa[patelesize]/D");
	tree->Branch("pateletrknumberOfValidHits",pateletrknumberOfValidHits,"pateletrknumberOfValidHits[patelesize]/I");
	tree->Branch("pateletrklost",pateletrklost,"pateletrklost[patelesize]/I");
	/*tree->Branch("pateletrkd0",pateletrkd0,"pateletrkd0[patelesize]/D");
	tree->Branch("pateletrkdxy",pateletrkdxy,"pateletrkdxy[patelesize]/D");
	tree->Branch("pateletrkdz",pateletrkdz,"pateletrkdz[patelesize]/D");
	*/
	tree->Branch("pateletrkptin",pateletrkptin,"pateletrkptin[patelesize]/D");
	tree->Branch("pateletrkptout",pateletrkptout,"pateletrkptout[patelesize]/D");
	tree->Branch("pateletrkfbrem",pateletrkfbrem,"pateletrkfbrem[patelesize]/D");
	tree->Branch("pateletrkqoverp",pateletrkqoverp,"pateletrkqoverp[patelesize]/D");
	tree->Branch("pateletrkvx",pateletrkvx,"pateletrkvx[patelesize]/D");
	tree->Branch("pateletrkvy",pateletrkvy,"pateletrkvy[patelesize]/D");
	tree->Branch("pateletrkvz",pateletrkvz,"pateletrkvz[patelesize]/D");
	tree->Branch("pateletrkphi",pateletrkphi,"pateletrkphi[patelesize]/D");
	tree->Branch("pateletrkndof",pateletrkndof,"pateletrkndof[patelesize]/D");
	tree->Branch("pateletrkrecHitsSize",pateletrkrecHitsSize,"pateletrkrecHitsSize[patelesize]/I");
	tree->Branch("pateletrktheta",pateletrktheta,"pateletrktheta[patelesize]/D");
	tree->Branch("pateletrkqualityMask",pateletrkqualityMask,"pateletrkqualityMask[patelesize]/I");
	tree->Branch("pateletrkouterX",pateletrkouterX,"pateletrkouterX[patelesize]/D");
	tree->Branch("pateletrkouterY",pateletrkouterY,"pateletrkouterY[patelesize]/D");
	tree->Branch("pateletrkouterZ",pateletrkouterZ,"pateletrkouterZ[patelesize]/D");
	tree->Branch("pateletrkouterRadius",pateletrkouterRadius,"pateletrkouterRadius[patelesize]/D");
	tree->Branch("pateletrkinnerX",pateletrkinnerX,"pateletrkinnerX[patelesize]/D");
	tree->Branch("pateletrkinnerY",pateletrkinnerY,"pateletrkinnerY[patelesize]/D");
	tree->Branch("pateletrkinnerZ",pateletrkinnerZ,"pateletrkinnerZ[patelesize]/D");
	
	
	//heepid
	tree->Branch("pateleheepcutword",pateleheepcutword,"pateleheepcutword[patelesize]/I");
	tree->Branch("pateleheepid",pateleheepid,"pateleheepid[patelesize]/I");
	tree->Branch("samspateleheepid",samspateleheepid,"samspateleheepid[patelesize]/I");
      }//if(usefulpateleTreeVar_)
    }//if(runPatelectrons_)
      
  //pat photons
  if(runPatphotons_)
    {
      
      patphop4 = new TClonesArray("TLorentzVector", MAX_PHOTONS);
      patphocalopos = new TClonesArray("TVector3", MAX_PHOTONS);
      patphoscp4 = new TClonesArray("TLorentzVector", MAX_PHOTONS);

      tree->Branch("patphosize",&patphosize,"patphosize/I");
      tree->Branch("patphop4", "TClonesArray", &patphop4, 32000, 0);
      tree->Branch("patphocalopos", "TClonesArray", &patphocalopos, 32000, 0);
      tree->Branch("patphoscp4", "TClonesArray", &patphoscp4, 32000, 0);
      tree->Branch("patphoscind",patphoscind,"patphoscind[patphosize]/I");

      tree->Branch("patphoecalRecHitSumEtConeDR03",patphoecalRecHitSumEtConeDR03,"patphoecalRecHitSumEtConeDR03[patphosize]/D");
      tree->Branch("patphohcalDepth1TowerSumEtConeDR03",patphohcalDepth1TowerSumEtConeDR03,"patphohcalDepth1TowerSumEtConeDR03[patphosize]/D");
      tree->Branch("patphohcalDepth2TowerSumEtConeDR03",patphohcalDepth2TowerSumEtConeDR03,"patphohcalDepth2TowerSumEtConeDR03[patphosize]/D");
      tree->Branch("patphohcalTowerSumEtConeDR03",patphohcalTowerSumEtConeDR03,"patphohcalTowerSumEtConeDR03[patphosize]/D");
      tree->Branch("patphotrkSumPtHollowConeDR03",patphotrkSumPtHollowConeDR03,"patphotrkSumPtHollowConeDR03[patphosize]/D");
      tree->Branch("patphotrkSumPtSolidConeDR03",patphotrkSumPtSolidConeDR03,"patphotrkSumPtSolidConeDR03[patphosize]/D");
      tree->Branch("patphonTrkHollowConeDR03",patphonTrkHollowConeDR03,"patphonTrkHollowConeDR03[patphosize]/D");
      tree->Branch("patphonTrkSolidConeDR03",patphonTrkSolidConeDR03,"patphonTrkSolidConeDR03[patphosize]/D");

      tree->Branch("patphoecalRecHitSumEtConeDR04",patphoecalRecHitSumEtConeDR04,"patphoecalRecHitSumEtConeDR04[patphosize]/D");
      tree->Branch("patphohcalDepth1TowerSumEtConeDR04",patphohcalDepth1TowerSumEtConeDR04,"patphohcalDepth1TowerSumEtConeDR04[patphosize]/D");
      tree->Branch("patphohcalDepth2TowerSumEtConeDR04",patphohcalDepth2TowerSumEtConeDR04,"patphohcalDepth2TowerSumEtConeDR04[patphosize]/D");
      tree->Branch("patphohcalTowerSumEtConeDR04",patphohcalTowerSumEtConeDR04,"patphohcalTowerSumEtConeDR04[patphosize]/D");
      tree->Branch("patphotrkSumPtHollowConeDR04",patphotrkSumPtHollowConeDR04,"patphotrkSumPtHollowConeDR04[patphosize]/D");
      tree->Branch("patphotrkSumPtSolidConeDR04",patphotrkSumPtSolidConeDR04,"patphotrkSumPtSolidConeDR04[patphosize]/D");
      tree->Branch("patphonTrkHollowConeDR04",patphonTrkHollowConeDR04,"patphonTrkHollowConeDR04[patphosize]/D");
      tree->Branch("patphonTrkSolidConeDR04",patphonTrkSolidConeDR04,"patphonTrkSolidConeDR04[patphosize]/D");
            
      tree->Branch("patphoe1x5",patphoe1x5,"patphoe1x5[patphosize]/D");
      tree->Branch("patphoe2x5",patphoe2x5,"patphoe2x5[patphosize]/D");
      tree->Branch("patphoe5x5",patphoe5x5,"patphoe5x5[patphosize]/D");
      tree->Branch("patphoeta",patphoeta,"patphoeta[patphosize]/D");
      tree->Branch("patphohadronicOverEm",patphohadronicOverEm,"patphohadronicOverEm[patphosize]/D");
      tree->Branch("patphosigmaIetaIeta",patphosigmaIetaIeta,"patphosigmaIetaIeta[patphosize]/D");
      tree->Branch("patphonumberOfTracks",patphonumberOfTracks,"patphonumberOfTracks[patphosize]/I");
      tree->Branch("patphor9",patphor9,"patphor9[patphosize]/D");
      
      tree->Branch("patphohasPixelSeed",patphohasPixelSeed,"patphohasPixelSeed[patphosize]/I");
      tree->Branch("patphoisConvertedPhoton",patphoisConvertedPhoton,"patphoisConvertedPhoton[patphosize]/I");
      tree->Branch("patphomaxEnergyXtal",patphomaxEnergyXtal,"patphomaxEnergyXtal[patphosize]/D");

      //gaps info
      tree->Branch("patphoisEB",patphoisEB,"patphoisEB[patphosize]/I");
      tree->Branch("patphoisEE",patphoisEE,"patphoisEE[patphosize]/I");
      
      //spikes
      tree->Branch("patphoswissCross",patphoswissCross,"patphoswissCross[patphosize]/D");
      tree->Branch("patphoswissBasedspikevar",patphoswissBasedspikevar,"patphoswissBasedspikevar[patphosize]/D");

      //seed timing
      tree->Branch("patphoseedtime",patphoseedtime,"patphoseedtime[patphosize]/D");
      tree->Branch("patphorecoFlag",patphorecoFlag,"patphorecoFlag[patphosize]/I");
      //tree->Branch("patphocheckFlag",patphocheckFlag,"patphocheckFlag[patphosize]/I");
      //tree->Branch("patphokOutOfTime",patphokOutOfTime,"patphokOutOfTime[patphosize]/D");
      //tree->Branch("patphoseverityLevel",patphoseverityLevel,"patphoseverityLevel[patphosize]/D");
      tree->Branch("patphoe2e9",patphoe2e9,"patphoe2e9[patphosize]/D");

      //conversion
      tree->Branch("patphoconvsize",patphoconvsize,"patphosconvsize[patphosize]/I");
      tree->Branch("patphohasConversionTracks",patphohasConversionTracks,"patphosconvsize[patphosize]/I");
      tree->Branch("patphoconvtxX",patphoconvtxX,"patphosconvtxX[patphosize][25]/D");
      tree->Branch("patphoconvtxY",patphoconvtxY,"patphosconvtxX[patphosize][25]/D");
      tree->Branch("patphoconvtxZ",patphoconvtxZ,"patphosconvtxX[patphosize][25]/D");
      tree->Branch("patphoconvtxR",patphoconvtxR,"patphosconvtxX[patphosize][25]/D");
      ////PF photon iso
      ///19th July
      tree->Branch("patphopfchargediso",patphopfchargediso,"patphopfchargediso[patphosize]/D");
      tree->Branch("patphopfphotoniso",patphopfphotoniso,"patphopfphotoniso[patphosize]/D");
      tree->Branch("patphopfneutraliso",patphopfneutraliso,"patphopfneutraliso[patphosize]/D");

      tree->Branch("patphoPFphotonWorstChargedHadronIso",patphoPFphotonWorstChargedHadronIso,"patphoPFphotonWorstChargedHadronIso[patphosize]/D");

      if(runSCremoval_){
	tree->Branch("patphonewpfchargediso",patphonewpfchargediso,"patphonewpfchargediso[patphosize]/D");
	tree->Branch("patphonewpfphotoniso",patphonewpfphotoniso,"patphonewpfphotoniso[patphosize]/D");
	tree->Branch("patphonewpfneutraliso",patphonewpfneutraliso,"patphonewpfneutraliso[patphosize]/D");
      }
      
      ///conversion safe electron veto
      tree->Branch("patphopasselectronveto",patphopasselectronveto,"patphopasselectronveto[patphosize]/I");

      ////new HoE for photon
      tree->Branch("patphohadTowOverEm",patphohadTowOverEm,"patphohadTowOverEm[patphosize]/D");

      /////20th Oct
      tree->Branch("patphomipChi2",patphomipChi2,"patphomipChi2[patphosize]/D");
      tree->Branch("patphomipIntercept",patphomipIntercept,"patphomipIntercept[patphosize]/D");
      tree->Branch("patphomipIsHalo",patphomipIsHalo,"patphomipIsHalo[patphosize]/I");
      
      tree->Branch("patphomipNhitCone",patphomipNhitCone,"patphomipNhitCone[patphosize]/I");

      tree->Branch("patphomipSlope",patphomipSlope,"patphomipSlope[patphosize]/D");
      tree->Branch("patphomipTotEnergy",patphomipTotEnergy,"patphomipTotEnergy[patphosize]/D");
      
      //gen match
      
      tree->Branch("patphomGenpdgId",patphomGenpdgId,"patphomGenpdgId[patphosize]/I");
      tree->Branch("patphomGenstatus",patphomGenstatus,"patphomGenstatus[patphosize]/I");
      tree->Branch("patphonummoth",patphonummoth,"patphonummoth[patphosize]/I");
      tree->Branch("patphomGenmompdgId",patphomGenmompdgId,"patphomGenmompdgId[patphosize][10]/I");
      tree->Branch("patphomGengranpdgId",patphomGengranpdgId,"patphomGengranpdgId[patphosize]/I");
      
      if(usefulpatphoTreeVar_){
	//tight id
	tree->Branch("patphoTightIDcutword",patphoTightIDcutword,"patphoTightIDcutword[patphosize]/I");
	tree->Branch("patphotightid",patphotightid,"patphotightid[patphosize]/I");
	
	//gaps info
	tree->Branch("patphoisEBEEGap",patphoisEBEEGap,"patphoisEBEEGap[patphosize]/I");
	tree->Branch("patphoisEBEtaGap",patphoisEBEtaGap,"patphoisEBEtaGap[patphosize]/I");
	tree->Branch("patphoisEBGap",patphoisEBGap,"patphoisEBGap[patphosize]/I");
	tree->Branch("patphoisEBPhiGap",patphoisEBPhiGap,"patphoisEBPhiGap[patphosize]/I");
	
	tree->Branch("patphoisEEDeeGap",patphoisEEDeeGap,"patphoisEEDeeGap[patphosize]/I");
	tree->Branch("patphoisEEGap",patphoisEEGap,"patphoisEEGap[patphosize]/I");
	tree->Branch("patphoisEERingGap",patphoisEERingGap,"patphoisEERingGap[patphosize]/I");
	
	
	
	//gsf
	tree->Branch("patphotrkpt",patphotrkpt,"patphotrkpt[patphosize]/D");
	tree->Branch("patphotrkcharge",patphotrkcharge,"patphotrkcharge[patphosize]/I");
	tree->Branch("patphotrkchi2",patphotrkchi2,"patphotrkchi2[patphosize]/D");
	tree->Branch("patphotrketa",patphotrketa,"patphotrketa[patphosize]/D");
	tree->Branch("patphotrknumberOfLostHits",patphotrknumberOfLostHits,"patphotrknumberOfLostHits[patphosize]/I");
	tree->Branch("patphotrknumberOfValidHits",patphotrknumberOfValidHits,"patphotrknumberOfValidHits[patphosize]/I");
	tree->Branch("patphotrklost",patphotrklost,"patphotrklost[patphosize]/I");
	tree->Branch("patphotrkd0",patphotrkd0,"patphotrkd0[patphosize]/D");
	tree->Branch("patphotrkdxy",patphotrkdxy,"patphotrkdxy[patphosize]/D");
	tree->Branch("patphotrkdz",patphotrkdz,"patphotrkdz[patphosize]/D");
	tree->Branch("patphotrkptin",patphotrkptin,"patphotrkptin[patphosize]/D");
	tree->Branch("patphotrkptout",patphotrkptout,"patphotrkptout[patphosize]/D");
	tree->Branch("patphotrkfbrem",patphotrkfbrem,"patphotrkfbrem[patphosize]/D");
	tree->Branch("patphotrkqoverp",patphotrkqoverp,"patphotrkqoverp[patphosize]/D");
	tree->Branch("patphotrkvx",patphotrkvx,"patphotrkvx[patphosize]/D");
	tree->Branch("patphotrkvy",patphotrkvy,"patphotrkvy[patphosize]/D");
	tree->Branch("patphotrkvz",patphotrkvz,"patphotrkvz[patphosize]/D");
	tree->Branch("patphotrkphi",patphotrkphi,"patphotrkphi[patphosize]/D");
	tree->Branch("patphotrkndof",patphotrkndof,"patphotrkndof[patphosize]/D");
	tree->Branch("patphotrkrecHitsSize",patphotrkrecHitsSize,"patphotrkrecHitsSize[patphosize]/I");
	tree->Branch("patphotrktheta",patphotrktheta,"patphotrktheta[patphosize]/D");
	tree->Branch("patphotrkqualityMask",patphotrkqualityMask,"patphotrkqualityMask[patphosize]/I");
	tree->Branch("patphotrkouterX",patphotrkouterX,"patphotrkouterX[patphosize]/D");
	tree->Branch("patphotrkouterY",patphotrkouterY,"patphotrkouterY[patphosize]/D");
	tree->Branch("patphotrkouterZ",patphotrkouterZ,"patphotrkouterZ[patphosize]/D");
	tree->Branch("patphotrkouterRadius",patphotrkouterRadius,"patphotrkouterRadius[patphosize]/D");
	tree->Branch("patphotrkinnerX",patphotrkinnerX,"patphotrkinnerX[patphosize]/D");
	tree->Branch("patphotrkinnerY",patphotrkinnerY,"patphotrkinnerY[patphosize]/D");
	tree->Branch("patphotrkinnerZ",patphotrkinnerZ,"patphotrkinnerZ[patphosize]/D");

	///gen match
	tree->Branch("patphomGenisJet",patphomGenisJet,"patphomGenisJet[patphosize]/I");
	tree->Branch("patphomGenisPhoton",patphomGenisPhoton,"patphomGenisPhoton[patphosize]/I");
	tree->Branch("patphomGenisElectron",patphomGenisElectron,"patphomGenisElectron[patphosize]/I");
	tree->Branch("patphomGentheta",patphomGentheta,"patphomGentheta[patphosize]/D");
	tree->Branch("patphomGeneta",patphomGeneta,"patphomGeneta[patphosize]/D");
	tree->Branch("patphomGenphi",patphomGenphi,"patphomGenphi[patphosize]/D"); 
	tree->Branch("patphomGenpt",patphomGenpt,"patphomGenpt[patphosize]/D");
	tree->Branch("patphomGenpx",patphomGenpx,"patphomGenpx[patphosize]/D");
	tree->Branch("patphomGenpy",patphomGenpy,"patphomGenpy[patphosize]/D");
	tree->Branch("patphomGenpz",patphomGenpz,"patphomGenpz[patphosize]/D");
	tree->Branch("patphomGenenergy",patphomGenenergy,"patphomGenenergy[patphosize]/D");
	
      }//if(usefulpatphotreeVar_)
      
      
      if(runPhoRechitInfo_)
	{
	  tree->Branch("patphoncrys",patphoncrys,"patphoncrys[patphosize]/I");
	  tree->Branch("patphocrysrawId",patphocrysrawId,"patphocrysrawId[patphosize][1000]/I");
	  tree->Branch("patphocrysieta",patphocrysieta,"patphocrysieta[patphosize][1000]/I");
	  tree->Branch("patphocrysiphi",patphocrysiphi,"patphocrysiphi[patphosize][1000]/I");
	  tree->Branch("patphocrysrecoFlag",patphocrysrecoFlag,"patphocrysrecoFlag[patphosize][1000]/I");
	  //tree->Branch("patphocryscheckFlag",patphocryscheckFlag,"patphocryscheckFlag[patphosize][1000]/I");
	  tree->Branch("patphocrysenergy",patphocrysenergy,"patphocrysenergy[patphosize][1000]/D");
	  tree->Branch("patphocrystime",patphocrystime,"patphocrystime[patphosize][1000]/D");
	  tree->Branch("patphocrystimeErr",patphocrystimeErr,"patphocrystimeErr[patphosize][1000]/D");

	  ////23rd Aug
	  tree->Branch("patphosigmaIetaIphi",patphosigmaIetaIphi,"patphosigmaIetaIphi[patphosize]/D");
	  tree->Branch("patphosigmaIphiIphi",patphosigmaIphiIphi,"patphosigmaIphiIphi[patphosize]/D");
	}

      
    }//if(runPatphotons_)

  //PFMET
  if(runPFmet_)
    {
      tree->Branch("pfmetsize",&pfmetsize,"pfmetsize/I");
      tree->Branch("pfmetpt",pfmetpt,"pfmetpt[pfmetsize]/D");
      tree->Branch("pfmetphi",pfmetphi,"pfmetphi[pfmetsize]/D");
      tree->Branch("pfmetsumEt",pfmetsumEt,"pfmetsumEt[pfmetsize]/D");
      tree->Branch("pfmetpx",pfmetpx,"pfmetpx[pfmetsize]/D");
      tree->Branch("pfmetpy",pfmetpy,"pfmetpy[pfmetsize]/D");
    }

   
    //TCMET
  if(runTCmet_)
    {
      tree->Branch("tcmetsize",&tcmetsize,"tcmetsize/I");
      tree->Branch("tcmetpt",tcmetpt,"tcmetpt[tcmetsize]/D");
      tree->Branch("tcmetphi",tcmetphi,"tcmetphi[tcmetsize]/D");
      tree->Branch("tcmetsumEt",tcmetsumEt,"tcmetsumEt[tcmetsize]/D");
      tree->Branch("tcmetpx",tcmetpx,"tcmetpx[tcmetsize]/D");
      tree->Branch("tcmetpy",tcmetpy,"tcmetpy[tcmetsize]/D");
    }


   //MET
  if(runmet_)
    {
      tree->Branch("metsize",&metsize,"metsize/I");
      tree->Branch("metpt",metpt,"metpt[metsize]/D");
      tree->Branch("metphi",metphi,"metphi[metsize]/D");
      tree->Branch("metsumEt",metsumEt,"metsumEt[metsize]/D");
      tree->Branch("metpx",metpx,"metpx[metsize]/D");
      tree->Branch("metpy",metpy,"metpy[metsize]/D");
    }
    




  //HLT
  if(runHLT_ == 1)
    {
      tree->Branch("is_HLT_Ele10_SW_L1R_event",&is_HLT_Ele10_SW_L1R_event,"is_HLT_Ele10_SW_L1R_event/I");
      tree->Branch("is_HLT_Ele15_SW_L1R_event",&is_HLT_Ele15_SW_L1R_event,"is_HLT_Ele15_SW_L1R_event/I");
      tree->Branch("is_HLT_Ele15_SW_EleId_L1R_event",&is_HLT_Ele15_SW_EleId_L1R_event,"is_HLT_Ele15_SW_EleId_L1R_event/I");
      tree->Branch("is_HLT_Ele15_SW_LooseTrackIso_L1R_event",&is_HLT_Ele15_SW_LooseTrackIso_L1R_event,"is_HLT_Ele15_SW_LooseTrackIso_L1R_event/I");
      tree->Branch("is_HLT_Ele15_SC15_SW_LooseTrackIso_L1R_event",&is_HLT_Ele15_SC15_SW_LooseTrackIso_L1R_event,"is_HLT_Ele15_SC15_SW_LooseTrackIso_L1R_event/I");
      tree->Branch("is_HLT_Ele15_SC15_SW_EleId_L1R_event",&is_HLT_Ele15_SC15_SW_EleId_L1R_event,"is_HLT_Ele15_SC15_SW_EleId_L1R_event/I");
      tree->Branch("is_HLT_Ele20_SW_L1R_event",&is_HLT_Ele20_SW_L1R_event,"is_HLT_Ele20_SW_L1R_event/I");
      tree->Branch("is_HLT_Ele20_SC15_SW_L1R_event",&is_HLT_Ele20_SC15_SW_L1R_event,"is_HLT_Ele20_SC15_SW_L1R_event/I");
      tree->Branch("is_HLT_Ele25_SW_L1R_event",&is_HLT_Ele25_SW_L1R_event,"is_HLT_Ele25_SW_L1R_event/I");
      tree->Branch("is_HLT_Ele25_SW_EleId_LooseTrackIso_L1R_event",&is_HLT_Ele25_SW_EleId_LooseTrackIso_L1R_event,"is_HLT_Ele25_SW_EleId_LooseTrackIso_L1R_event/I");
      tree->Branch("is_HLT_DoubleEle5_SW_Jpsi_L1R_event",&is_HLT_DoubleEle5_SW_Jpsi_L1R_event,"is_HLT_DoubleEle5_SW_Jpsi_L1R_event/I");
      tree->Branch("is_HLT_DoubleEle5_SW_Upsilon_L1R_event",&is_HLT_DoubleEle5_SW_Upsilon_L1R_event,"is_HLT_DoubleEle5_SW_Upsilon_L1R_event/I");
      tree->Branch("is_HLT_DoubleEle10_SW_L1R_event",&is_HLT_DoubleEle10_SW_L1R_event,"is_HLT_DoubleEle10_SW_L1R_event/I");

      tree->Branch("is_HLT_L1SingleEG5_event",&is_HLT_L1SingleEG5_event,"is_HLT_L1SingleEG5_event/I");
      tree->Branch("is_HLT_L1SingleEG8_event",&is_HLT_L1SingleEG8_event,"is_HLT_L1SingleEG8_event/I");
      tree->Branch("is_HLT_Ele10_LW_L1R_event",&is_HLT_Ele10_LW_L1R_event,"is_HLT_Ele10_LW_L1R_event/I");
      tree->Branch("is_HLT_Ele10_LW_EleId_L1R_event",&is_HLT_Ele10_LW_EleId_L1R_event,"is_HLT_Ele10_LW_EleId_L1R_event/I");
      tree->Branch("is_HLT_Ele15_LW_L1R_event",&is_HLT_Ele15_LW_L1R_event,"is_HLT_Ele15_LW_L1R_event/I");
      tree->Branch("is_HLT_Ele15_SiStrip_L1R_event",&is_HLT_Ele15_SiStrip_L1R_event,"is_HLT_Ele15_SiStrip_L1R_event/I");
      tree->Branch("is_HLT_Ele15_SC10_LW_L1R_event",&is_HLT_Ele15_SC10_LW_L1R_event,"is_HLT_Ele15_SC10_LW_L1R_event/I");
      tree->Branch("is_HLT_Ele20_LW_L1R_event",&is_HLT_Ele20_LW_L1R_event,"is_HLT_Ele20_LW_L1R_event/I");
      tree->Branch("is_HLT_L1DoubleEG5_event",&is_HLT_L1DoubleEG5_event,"is_HLT_L1DoubleEG5_event/I");
      tree->Branch("is_HLT_DoubleEle5_SW_L1R_event",&is_HLT_DoubleEle5_SW_L1R_event,"is_HLT_DoubleEle5_SW_L1R_event/I");

      tree->Branch("is_HLT_Photon10_L1R_event",&is_HLT_Photon10_L1R_event,"is_HLT_Photon10_L1R_event/I");
      tree->Branch("is_HLT_Photon10_LooseEcalIso_TrackIso_L1R_event",&is_HLT_Photon10_LooseEcalIso_TrackIso_L1R_event,"is_HLT_Photon10_LooseEcalIso_TrackIso_L1R_event/I");
      tree->Branch("is_HLT_Photon15_L1R_event",&is_HLT_Photon15_L1R_event,"is_HLT_Photon15_L1R_event/I");
      tree->Branch("is_HLT_Photon20_LooseEcalIso_TrackIso_L1R_event",&is_HLT_Photon20_LooseEcalIso_TrackIso_L1R_event,"is_HLT_Photon20_LooseEcalIso_TrackIso_L1R_event/I");
      tree->Branch("is_HLT_Photon25_L1R_event",&is_HLT_Photon25_L1R_event,"is_HLT_Photon25_L1R_event/I");
      tree->Branch("is_HLT_Photon25_LooseEcalIso_TrackIso_L1R_event",&is_HLT_Photon25_LooseEcalIso_TrackIso_L1R_event,"is_HLT_Photon25_LooseEcalIso_TrackIso_L1R_event/I");
      tree->Branch("is_HLT_Photon30_L1R_1E31_event",&is_HLT_Photon30_L1R_1E31_event,"is_HLT_Photon30_L1R_1E31_event/I");
      tree->Branch("is_HLT_Photon70_L1R_event",&is_HLT_Photon70_L1R_event,"is_HLT_Photon70_L1R_event/I");
      tree->Branch("is_HLT_DoublePhoton10_L1R_event",&is_HLT_DoublePhoton10_L1R_event,"is_HLT_DoublePhoton10_L1R_event/I");
      tree->Branch("is_HLT_DoublePhoton15_L1R_event",&is_HLT_DoublePhoton15_L1R_event,"is_HLT_DoublePhoton15_L1R_event/I");
      tree->Branch("is_HLT_DoublePhoton15_VeryLooseEcalIso_L1R_event",&is_HLT_DoublePhoton15_VeryLooseEcalIso_L1R_event,"is_HLT_DoublePhoton15_VeryLooseEcalIso_L1R_event/I");

      tree->Branch("is_HLT_Photon20_Cleaned_L1R_event",&is_HLT_Photon20_Cleaned_L1R_event,"is_HLT_Photon20_Cleaned_L1R_event/I");
      tree->Branch("is_HLT_Photon30_Cleaned_L1R_event",&is_HLT_Photon30_Cleaned_L1R_event,"is_HLT_Photon30_Cleaned_L1R_event/I");
      tree->Branch("is_HLT_Photon50_Cleaned_L1R_event",&is_HLT_Photon50_Cleaned_L1R_event,"is_HLT_Photon50_Cleaned_L1R_event/I");
      tree->Branch("is_HLT_Photon70_Cleaned_L1R_event",&is_HLT_Photon70_Cleaned_L1R_event,"is_HLT_Photon70_Cleaned_L1R_event/I");
      //e*
      //tree->Branch("is_HLT_DoublePhoton20_L1R_event",&is_HLT_DoublePhoton20_L1R_event,"is_HLT_DoublePhoton20_L1R_event/I");
      //tree->Branch("is_HLT_DoublePhoton17_L1R_event",&is_HLT_DoublePhoton17_L1R_event,"is_HLT_DoublePhoton10_L1R_event/I");

      tree->Branch("is_HLT_Jet15U",&is_HLT_Jet15U,"is_HLT_Jet15U/I");
      tree->Branch("is_HLT_Jet30U",&is_HLT_Jet30U,"is_HLT_Jet30U/I");
      tree->Branch("is_HLT_Jet50U",&is_HLT_Jet50U,"is_HLT_Jet50U/I");
      tree->Branch("is_HLT_Jet70U",&is_HLT_Jet70U,"is_HLT_Jet70U/I");
      tree->Branch("is_HLT_Jet100U",&is_HLT_Jet100U,"is_HLT_Jet100U/I");

      tree->Branch("is_HLT_Jet70U_v2",&is_HLT_Jet70U_v2,"is_HLT_Jet70U_v2/I");
      tree->Branch("is_HLT_Jet100U_v2",&is_HLT_Jet100U_v2,"is_HLT_Jet100U_v2/I");
	

      //HLT-RECO object match
      tree->Branch("ntrigToMatch",&ntrigToMatch,"ntrigToMatch/I");
      tree->Branch("hltprescale",hltprescale,"hltprescale[ntrigToMatch]/I");
	    
      tree->Branch("hltobjsize",hltobjsize,"hltobjsize[ntrigToMatch]/I");
      tree->Branch("ohtrigpt",ohtrigpt,"ohtrigpt[ntrigToMatch][100]/D");
      tree->Branch("ohtrigeta",ohtrigeta,"ohtrigeta[ntrigToMatch][100]/D");
      tree->Branch("ohtrigphi",ohtrigphi,"ohtrigphi[ntrigToMatch][100]/D");


      /*if( runPatphotons_ )
	{
	  tree->Branch("ohltphodR",ohltphodR,"ohltphodR[ntrigToMatch][100]/D");
	  tree->Branch("patphohltmindRmatch",patphohltmindRmatch,"patphohltmindRmatch[ntrigToMatch][patphosize]/D");
	}
      
      if( runPatelectrons_ )
	{
	  tree->Branch("ohltele1dR",ohltele1dR,"ohltele1dR[ntrigToMatch][100]/D");
	  tree->Branch("ohltele2dR",ohltele2dR,"ohltele2dR[ntrigToMatch][100]/D");  
	}
      */

      if( runpatJets_ )
	tree->Branch("ohltjetdR",ohltjetdR,"ohltjetdR[ntrigToMatch][100]/D");

    }//if(runHLT_ == 1)




  if(runSC_)
    {
      
        scp4 = new TClonesArray("TLorentzVector", MAX_SUPERCLUSTERS);
	scxyz = new TClonesArray("TVector3",MAX_SUPERCLUSTERS);
	tree->Branch("scsize",&scsize,"scsize/I");
	tree->Branch("scEBsize",&scEBsize,"scEBsize/I");
	tree->Branch("scEEsize",&scEEsize,"scEEsize/I");
	
	tree->Branch("scp4", "TClonesArray", &scp4, 32000, 0);
	tree->Branch("scxyz", "TClonesArray", &scxyz, 32000, 0);
	
	tree->Branch("scetaWidth",scetaWidth,"scetaWidth[scsize]/D");
	tree->Branch("scphiWidth",scphiWidth,"scphiWidth[scsize]/D");
	tree->Branch("scrawEnergy",scrawEnergy,"scrawEnergy[scsize]/D");
	tree->Branch("scpreshowerEnergy",scpreshowerEnergy,"scpreshowerEnergy[scsize]/D");
	tree->Branch("scclustersSize",scclustersSize,"scclustersSize[scsize]/D");
	tree->Branch("scisEB",scisEB,"scisEB[scsize]/I");
	tree->Branch("scisEE",scisEE,"scisEE[scsize]/I");
	//cluster shape info
	tree->Branch("sceMax",sceMax,"sceMax[scsize]/D");
	tree->Branch("sce2x2",sce2x2,"sce2x2[scsize]/D");
	tree->Branch("sce3x3",sce3x3,"sce3x3[scsize]/D");
	tree->Branch("sce5x5",sce5x5,"sce5x5[scsize]/D");
	tree->Branch("scr4",scr4,"scr4[scsize]/D");
	tree->Branch("scr9",scr9,"scr9[scsize]/D");
	tree->Branch("scr25",scr25,"scr25[scsize]/D");
	tree->Branch("sce1bye4",sce1bye4,"sce1bye4[scsize]/D");
	tree->Branch("sce1bye9",sce1bye9,"sce1bye9[scsize]/D");
	tree->Branch("sce1bye25",sce1bye25,"sce1bye25[scsize]/D");
	tree->Branch("sce2e9",sce2e9,"sce2e9[scsize]/D");
	tree->Branch("scHoE1",scHoE1,"scHoE1[scsize]/D");
	tree->Branch("scHoE2",scHoE2,"scHoE2[scsize]/D");
	tree->Branch("scHoE",scHoE,"scHoE[scsize]/D");
    }//if(runSC_)
  
  if( rungenParticle_ ) 
    {
      genp4 = new TClonesArray("TLorentzVector", MAX_GENERATOR);
      genvtx = new TClonesArray("TVector3", MAX_GENERATOR);

      tree->Branch("gensize",&gensize,"gensize/I");
      tree->Branch("genp4", "TClonesArray", &genp4, 32000, 0);
      tree->Branch("genvtx", "TClonesArray", &genvtx, 32000, 0);
      tree->Branch("genstatus", genstatus, "genstatus[gensize]/I");
      tree->Branch("gencharge", gencharge, "gencharge[gensize]/I");
      tree->Branch("genpdgid", genpdgid, "genpdgid[gensize]/I");
      tree->Branch("genmother", genmother, "genmother[gensize]/I");
      tree->Branch("genndau", genndau, "genndau[gensize]/I");
      tree->Branch("gennmoth", gennmoth, "gennmoth[gensize]/I");
      
      
    }//if(rungenParticle_)

  //genjets
  if(rungenJets_)
    {
      tree->Branch("genjetsize",&genjetsize,"genjetsize/I");
      tree->Branch("genjetpt",genjetpt,"genjetpt[genjetsize]/D");
      tree->Branch("genjetpx",genjetpx,"genjetpx[genjetsize]/D");
      tree->Branch("genjetpy",genjetpy,"genjetpy[genjetsize]/D");
      tree->Branch("genjetpz",genjetpz,"genjetpz[genjetsize]/D");
      tree->Branch("genjetenergy",genjetenergy,"genjetenergy[genjetsize]/D");
      tree->Branch("genjetphi",genjetphi,"genjetphi[genjetsize]/D");
      tree->Branch("genjeteta",genjeteta,"genjeteta[genjetsize]/D");
      tree->Branch("genjetemEnergy",genjetemEnergy,"genjetemEnergy[genjetsize]/D");
      tree->Branch("genjethadEnergy",genjethadEnergy,"genjethadEnergy[genjetsize]/D");
      tree->Branch("genjetinvisibleEnergy",genjetinvisibleEnergy,"genjetinvisibleEnergy[genjetsize]/D");

    }// if(rungenJets_)


  //patjets
  if(runpatJets_)
    {
      tree->Branch("patjetsize",&patjetsize,"patjetsize/I");
      tree->Branch("patjetpt",patjetpt,"patjetpt[patjetsize]/D");
      tree->Branch("patjetpx",patjetpx,"patjetpx[patjetsize]/D");
      tree->Branch("patjetpy",patjetpy,"patjetpy[patjetsize]/D");
      tree->Branch("patjetpz",patjetpz,"patjetpz[patjetsize]/D");
      tree->Branch("patjetet",patjetpt,"patjetet[patjetsize]/D");
      tree->Branch("patjetenergy",patjetenergy,"patjetenergy[patjetsize]/D");
      tree->Branch("patjetphi",patjetphi,"patjetphi[patjetsize]/D");
      tree->Branch("patjeteta",patjeteta,"patjeteta[patjetsize]/D");
      tree->Branch("patjethasOverlapsmu",patjethasOverlapsmu,"patjethasOverlapsmu[patjetsize]/D");
      tree->Branch("patjethasOverlapsele",patjethasOverlapsele,"patjethasOverlapsele[patjetsize]/I");
      tree->Branch("patjethasOverlapspho",patjethasOverlapspho,"patjethasOverlapspho[patjetsize]/I");
      tree->Branch("patjethasOverlapstau",patjethasOverlapstau,"patjethasOverlapstau[patjetsize]/I");
      tree->Branch("patjethasOverlapstkIsoele",patjethasOverlapstkIsoele,"patjethasOverlapstkIsoele[patjetsize]/I");
      tree->Branch("patjetchargedEmEnergy",patjetchargedEmEnergy,"patjetchargedEmEnergy[patjetsize]/D");
      tree->Branch("patjetchargedEmEnergyFraction",patjetchargedEmEnergyFraction,"patjetchargedEmEnergyFraction[patjetsize]/D");
      tree->Branch("patjetchargedHadronEnergy",patjetchargedHadronEnergy,"patjetchargedHadronEnergy[patjetsize]/D");
      tree->Branch("patjetchargedHadronEnergyFraction",patjetchargedHadronEnergyFraction,"patjetchargedHadronEnergyFraction[patjetsize]/D");
      tree->Branch("patjetchargedMultiplicity",patjetchargedMultiplicity,"patjetchargedMultiplicity[patjetsize]/D");
      tree->Branch("patjetemEnergyFraction",patjetemEnergyFraction,"patjetemEnergyFraction[patjetsize]/D");
      tree->Branch("patjetemEnergyInEB",patjetemEnergyInEB,"patjetemEnergyInEB[patjetsize]/D");
      tree->Branch("patjetemEnergyInEE",patjetemEnergyInEE,"patjetemEnergyInEE[patjetsize]/D");
      tree->Branch("patjetemEnergyInHF",patjetemEnergyInHF,"patjetemEnergyInHF[patjetsize]/D");
      tree->Branch("patjetenergyFractionHadronic",patjetenergyFractionHadronic,"patjetenergyFractionHadronic[patjetsize]/D");
      tree->Branch("patjethadEnergyInHB",patjethadEnergyInHB,"patjethadEnergyInHB[patjetsize]/D");
      tree->Branch("patjethadEnergyInHE",patjethadEnergyInHE,"patjethadEnergyInHE[patjetsize]/D");
      tree->Branch("patjethadEnergyInHF",patjethadEnergyInHF,"patjethadEnergyInHF[patjetsize]/D");
      tree->Branch("patjethadEnergyInHO",patjethadEnergyInHO,"patjethadEnergyInHO[patjetsize]/D");
      //tree->Branch("patjethadEnergy",patjethadEnergy,"patjethadEnergy[patjetsize]/D");

    }// if(runpatJets_)


  //pfjets
  if(runpfJets_)
    {
      tree->Branch("pfjetsize",&pfjetsize,"pfjetsize/I");
      tree->Branch("pfjetpt",pfjetpt,"pfjetpt[pfjetsize]/D");
      tree->Branch("pfjetpx",pfjetpx,"pfjetpx[pfjetsize]/D");
      tree->Branch("pfjetpy",pfjetpy,"pfjetpy[pfjetsize]/D");
      tree->Branch("pfjetpz",pfjetpz,"pfjetpz[pfjetsize]/D");
      tree->Branch("pfjetet",pfjetpt,"pfjetet[pfjetsize]/D");
      tree->Branch("pfjetenergy",pfjetenergy,"pfjetenergy[pfjetsize]/D");
      tree->Branch("pfjetphi",pfjetphi,"pfjetphi[pfjetsize]/D");
      tree->Branch("pfjeteta",pfjeteta,"pfjeteta[pfjetsize]/D");
      tree->Branch("pfjethasOverlapsmu",pfjethasOverlapsmu,"pfjethasOverlapsmu[pfjetsize]/D");
      tree->Branch("pfjethasOverlapsele",pfjethasOverlapsele,"pfjethasOverlapsele[pfjetsize]/I");
      tree->Branch("pfjethasOverlapspho",pfjethasOverlapspho,"pfjethasOverlapspho[pfjetsize]/I");
      tree->Branch("pfjethasOverlapstau",pfjethasOverlapstau,"pfjethasOverlapstau[pfjetsize]/I");
      tree->Branch("pfjethasOverlapstkIsoele",pfjethasOverlapstkIsoele,"pfjethasOverlapstkIsoele[pfjetsize]/I");
      tree->Branch("pfjetchargedEmEnergy",pfjetchargedEmEnergy,"pfjetchargedEmEnergy[pfjetsize]/D");
      tree->Branch("pfjetchargedEmEnergyFraction",pfjetchargedEmEnergyFraction,"pfjetchargedEmEnergyFraction[pfjetsize]/D");
      tree->Branch("pfjetchargedHadronEnergy",pfjetchargedHadronEnergy,"pfjetchargedHadronEnergy[pfjetsize]/D");
      tree->Branch("pfjetchargedHadronEnergyFraction",pfjetchargedHadronEnergyFraction,"pfjetchargedHadronEnergyFraction[pfjetsize]/D");
      tree->Branch("pfjetchargedMultiplicity",pfjetchargedMultiplicity,"pfjetchargedMultiplicity[pfjetsize]/D");
      tree->Branch("pfjetemEnergyFraction",pfjetemEnergyFraction,"pfjetemEnergyFraction[pfjetsize]/D");
      tree->Branch("pfjetemEnergyInEB",pfjetemEnergyInEB,"pfjetemEnergyInEB[pfjetsize]/D");
      tree->Branch("pfjetemEnergyInEE",pfjetemEnergyInEE,"pfjetemEnergyInEE[pfjetsize]/D");
      tree->Branch("pfjetemEnergyInHF",pfjetemEnergyInHF,"pfjetemEnergyInHF[pfjetsize]/D");
      tree->Branch("pfjetenergyFractionHadronic",pfjetenergyFractionHadronic,"pfjetenergyFractionHadronic[pfjetsize]/D");
      tree->Branch("pfjethadEnergyInHB",pfjethadEnergyInHB,"pfjethadEnergyInHB[pfjetsize]/D");
      tree->Branch("pfjethadEnergyInHE",pfjethadEnergyInHE,"pfjethadEnergyInHE[pfjetsize]/D");
      tree->Branch("pfjethadEnergyInHF",pfjethadEnergyInHF,"pfjethadEnergyInHF[pfjetsize]/D");
      tree->Branch("pfjethadEnergyInHO",pfjethadEnergyInHO,"pfjethadEnergyInHO[pfjetsize]/D");
      tree->Branch("pfjet_jecCorr",pfjet_jecCorr,"pfjet_jecCorr[pfjetsize]/D");
      
      //added on 31st October, 2011 for fake rate 
      tree->Branch("pfjetchargedEmEnergyFraction",pfjetchargedEmEnergyFraction,"pfjetchargedEmEnergyFraction[pfjetsize]/D");
      tree->Branch("pfjetchargedHadronEnergyFraction",pfjetchargedHadronEnergyFraction,"pfjetchargedHadronEnergyFraction[pfjetsize]/D");
      tree->Branch("pfjetchargedMuEnergyFraction",pfjetchargedMuEnergyFraction,"pfjetchargedMuEnergyFraction[pfjetsize]/D");
      
      //tree->Branch("pfjethadEnergy",pfjethadEnergy,"pfjethadEnergy[pfjetsize]/D");

    }// if(runpfJets_)


  //UCSD's Functions
  if(runUCSDHLT_)
    hlt->defineBranch(tree);
  
  if(runVertex_)
    vertex_std->defineBranch(tree);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Analyser::endJob() {
  //std::cout<<"No of events in this file = "<<nevents<<std::endl;
  //std::cout<<"No of events in which atleast one track > 5GeV = "<<ntracks5<<std::endl;
  if(wantLocalFile_)
    lf->WriteTObject(tree);

  if(wantRFIOFile_)
    rfiof->WriteTObject(tree);

  delete tree;
  std::cout<<"Written the root file"<<std::endl;

}

double Analyser::mye2overe9( const DetId id, const EcalRecHitCollection & recHits)
{
  ///////////start calculating e2/e9
  ////http://cmslxr.fnal.gov/lxr/source/RecoLocalCalo/EcalRecAlgos/src/EcalSeverityLevelAlgo.cc#240
  // compute e2overe9
  //  
  //   | | | |
  //   +-+-+-+
  //   | |1|2|
  //   +-+-+-+
  //   | | | |
  //
  //   1 - input hit,  2 - highest energy hit in a 3x3 around 1
// 
//   rechit 1 must have E_t > recHitEtThreshold
//   rechit 2 must have E_t > recHitEtThreshold2
//
//   function returns value of E2/E9 centered around 1 (E2=energy of hits 1+2) if energy of 1>2
//
//   if energy of 2>1 and KillSecondHit is set to true, function returns value of E2/E9 centered around 2
//   *provided* that 1 is the highest energy hit in a 3x3 centered around 2, otherwise, function returns 0
  
  
  float recHitEtThreshold = 10.0; 
float recHitEtThreshold2 = 1.0;
bool avoidIeta85=false;
bool KillSecondHit=true;



if ( id.subdetId() == EcalBarrel ) {

  EBDetId ebId( id );

  // avoid recHits at |eta|=85 where one side of the neighbours is missing
  if ( abs(ebId.ieta())==85 && avoidIeta85) return 0;

  // select recHits with Et above recHitEtThreshold
  float e1 = recHitE( id, recHits );
  float ete1=recHitApproxEt( id, recHits );


  // check that rechit E_t is above threshold
  if (ete1 < std::min(recHitEtThreshold,recHitEtThreshold2) ) return 0;

  if (ete1 < recHitEtThreshold && !KillSecondHit ) return 0;
  
  float e2=-1;
  float ete2=0;
  float s9 = 0;

  // coordinates of 2nd hit relative to central hit
  int e2eta=0;
  int e2phi=0;

  // LOOP OVER 3x3 ARRAY CENTERED AROUND HIT 1
  
  for ( int deta = -1; deta <= +1; ++deta ) {
    for ( int dphi = -1; dphi <= +1; ++dphi ) {

      // compute 3x3 energy 
      float etmp=recHitE( id, recHits, deta, dphi );
      s9 += etmp;

      EBDetId idtmp=EBDetId::offsetBy(id,deta,dphi);
      float eapproxet=recHitApproxEt( idtmp, recHits );
      
      // remember 2nd highest energy deposit (above threshold) in 3x3 array
      if (etmp>e2 && eapproxet>recHitEtThreshold2 && !(deta==0 && dphi==0)) {

	e2=etmp;
	ete2=eapproxet;
	e2eta=deta;
	e2phi=dphi;

      }

    }
  }

  if ( e1 == 0 )  return 0;

  // return 0 if 2nd hit is below threshold
  if ( e2 == -1 ) return 0;

  // compute e2/e9 centered around 1st hit
  float e2nd=e1+e2;
  float e2e9=0;

  if (s9!=0) e2e9=e2nd/s9;

  // if central hit has higher energy than 2nd hit
  //  return e2/e9 if 1st hit is above E_t threshold
  
  if (e1 > e2 && ete1>recHitEtThreshold) return e2e9;

  // if second hit has higher energy than 1st hit
  if ( e2 > e1 ) {


    // return 0 if user does not want to flag 2nd hit, or
    // hits are below E_t thresholds - note here we
    // now assume the 2nd hit to be the leading hit.
    
    if (!KillSecondHit || ete2<recHitEtThreshold || ete1<recHitEtThreshold2) {

      return 0;

    }


    else {

      // LOOP OVER 3x3 ARRAY CENTERED AROUND HIT 2 
      float s92nd=0;

      float e2nd_prime=0;
      int e2prime_eta=0;
      int e2prime_phi=0;

      EBDetId secondid=EBDetId::offsetBy(id,e2eta,e2phi);

      for ( int deta = -1; deta <= +1; ++deta ) {
	for ( int dphi = -1; dphi <= +1; ++dphi ) {

	  // compute 3x3 energy
	  float etmp=recHitE( secondid, recHits, deta, dphi );
	  s92nd += etmp;

	  if (etmp>e2nd_prime && !(deta==0 && dphi==0)) {
	    e2nd_prime=etmp;
	    e2prime_eta=deta;
	    e2prime_phi=dphi;
	  }

	}
      }

      // if highest energy hit around E2 is not the same as the input hit, return 0;
      if (!(e2prime_eta==-e2eta && e2prime_phi==-e2phi))
	{
	  return 0;
	}


      // compute E2/E9 around second hit 
      float e2e9_2=0;
      if (s92nd!=0) e2e9_2=e2nd/s92nd;

      //   return the value of E2/E9 calculated around 2nd hit
      return e2e9_2;


    }

  }
  
 } else if ( id.subdetId() == EcalEndcap ) {
  // only used for EB at the moment 
  return 0;
 }
return 0;
}//double Analyser::mye2overe9( const DetId id, const EcalRecHitCollection & recHits)


//define this as a plug-in
DEFINE_FWK_MODULE(Analyser);
