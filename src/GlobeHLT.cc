#include "Analyzer/Analyser/interface/Limits.h"

//#ifdef CMSSW_VERSION_210
#include "Analyzer/Analyser/interface/GlobeHLT.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
//Shilpi
#include "FWCore/Common/interface/TriggerNames.h"


#include "TLorentzVector.h"

#include <bitset>

/*void set(unsigned int& x, int bit) {
  x |= (1 << bit);
}

void reset(unsigned int& x, int bit) {
  x &= ~ (1 << bit);


bool check(unsigned int x, int bit) {
  return (x && (1 << bit));
}
*/

void set(ULong64_t& x, ULong64_t bit) {
  //Shilpi
  //std::cout<<"inside set; BEFORE SETTING x ; bit  "<<x<<" : "<<bit<<std::endl;
  ULong64_t a = 1;
  //x |= (1 << bit);
  x |= (a << bit);
  //std::cout<<("a << ")<<bit<<" : "<<(a << bit)<<std::endl;
  //std::cout<<"inside set; AFTER SETTING x = "<<x<<std::endl;
}

void reset(ULong64_t& x, ULong64_t bit) {
  ULong64_t a = 1;
  x &= ~ (a << bit);
}

bool check(ULong64_t x, ULong64_t bit) {
  ULong64_t a = 1;
  return (x && (a << bit));
}

GlobeHLT::GlobeHLT(const edm::ParameterSet& iConfig, const char* n): nome(n) {
  
  /*edm::ParameterSet psetHLT = iConfig.getParameter<edm::ParameterSet>("HLTParameters");
  inputTag_ = psetHLT.getParameter<edm::InputTag>("TriggerResultsTag");
  secondaryTriggerON  = psetHLT.getParameter<bool>("useSecondaryTrigger");
  hlt1Tag_  = psetHLT.getParameter<edm::InputTag>("PrimaryTriggerResultsTag");
  */ ///this is what UCSD uses. path:/uscms_data/d2/jshilpi/ShilpiGAPAT386/src/Analyzers/GlobeAnalyzer/python/globeanalyzer_37X_cfi.py
  
  /////for my work, I modify it
  inputTag_ = iConfig.getUntrackedParameter<edm::InputTag>("triggerEventTag");
  secondaryTriggerON  = iConfig.getUntrackedParameter<bool>("useSecondaryTrigger");
  hlt1Tag_  = iConfig.getUntrackedParameter<edm::InputTag>("HLTriggerResults");
  
  hlt_path_names_HLT1_1 = new std::vector<std::string>; hlt_path_names_HLT1_1->clear();
  hlt_path_names_HLT1_2 = new std::vector<std::string>; hlt_path_names_HLT1_2->clear();
  hlt_path_names_HLT1_3 = new std::vector<std::string>; hlt_path_names_HLT1_3->clear();
  hlt_path_names_HLT1_4 = new std::vector<std::string>; hlt_path_names_HLT1_4->clear();
  hlt_path_names_HLT1_5 = new std::vector<std::string>; hlt_path_names_HLT1_5->clear();
  hlt_path_names_HLT1_6 = new std::vector<std::string>; hlt_path_names_HLT1_6->clear();

  hlt_path_names_HLT1_7 = new std::vector<std::string>; hlt_path_names_HLT1_7->clear();
  hlt_path_names_HLT1_8 = new std::vector<std::string>; hlt_path_names_HLT1_8->clear();
  hlt_path_names_HLT1_9 = new std::vector<std::string>; hlt_path_names_HLT1_9->clear();
  hlt_path_names_HLT1_10 = new std::vector<std::string>; hlt_path_names_HLT1_10->clear();

 
  hlt_path_names_HLT1 = new std::vector<std::string>; hlt_path_names_HLT1->clear();

  ///prescales - Shilpi
  hlt_prescales_HLT1_1 = new std::vector<int>; hlt_prescales_HLT1_1->clear();
  hlt_prescales_HLT1_2 = new std::vector<int>; hlt_prescales_HLT1_2->clear();
  hlt_prescales_HLT1_3 = new std::vector<int>; hlt_prescales_HLT1_3->clear();
  hlt_prescales_HLT1_4 = new std::vector<int>; hlt_prescales_HLT1_4->clear();
  hlt_prescales_HLT1_5 = new std::vector<int>; hlt_prescales_HLT1_5->clear();
  hlt_prescales_HLT1_6 = new std::vector<int>; hlt_prescales_HLT1_6->clear();

  hlt_prescales_HLT1_7 = new std::vector<int>; hlt_prescales_HLT1_7->clear();
  hlt_prescales_HLT1_8 = new std::vector<int>; hlt_prescales_HLT1_8->clear();
  hlt_prescales_HLT1_9 = new std::vector<int>; hlt_prescales_HLT1_9->clear();
  hlt_prescales_HLT1_10 = new std::vector<int>; hlt_prescales_HLT1_10->clear();
  
  

  if(secondaryTriggerON){
    //hlt2Tag_  = psetHLT.getParameter<edm::InputTag>("SecondaryTriggerResultsTag"); ///UCSD
    hlt2Tag_  = iConfig.getParameter<edm::InputTag>("HLTriggerResults");
    hlt_path_names_HLT2_1 = new std::vector<std::string>; hlt_path_names_HLT2_1->clear();
    hlt_path_names_HLT2_2 = new std::vector<std::string>; hlt_path_names_HLT2_2->clear();
    hlt_path_names_HLT2_3 = new std::vector<std::string>; hlt_path_names_HLT2_3->clear();
    hlt_path_names_HLT2_4 = new std::vector<std::string>; hlt_path_names_HLT2_4->clear();
    hlt_path_names_HLT2_5 = new std::vector<std::string>; hlt_path_names_HLT2_5->clear();
    hlt_path_names_HLT2_6 = new std::vector<std::string>; hlt_path_names_HLT2_6->clear();
  }
  theHLTLabels.clear();
  hlt_label_names_1 = new std::vector<std::string>; hlt_label_names_1->clear();
  hlt_label_names_2 = new std::vector<std::string>; hlt_label_names_2->clear();
  hlt_label_names_3 = new std::vector<std::string>; hlt_label_names_3->clear();
  hlt_label_names_4 = new std::vector<std::string>; hlt_label_names_4->clear();
  hlt_label_names_5 = new std::vector<std::string>; hlt_label_names_5->clear();
  hlt_label_names_6 = new std::vector<std::string>; hlt_label_names_6->clear();

  hlt_label_names_7 = new std::vector<std::string>; hlt_label_names_7->clear();
  hlt_label_names_8 = new std::vector<std::string>; hlt_label_names_8->clear();
  hlt_label_names_9 = new std::vector<std::string>; hlt_label_names_9->clear();
  hlt_label_names_10 = new std::vector<std::string>; hlt_label_names_10->clear();
  
  electronColl   = iConfig.getParameter<edm::InputTag>("ElectronColl_std");
  photonCollStd  = iConfig.getParameter<edm::InputTag>("PhotonCollStd");
  muonColl       = iConfig.getParameter<edm::InputTag>("MuonColl");
  jetColl        = iConfig.getParameter<edm::InputTag>("JetColl_it5");
  
  debug_level    = iConfig.getParameter<int>("Debug_Level");
  
  gCUT           = new GlobeCuts(iConfig);
  
  debugHLT_ = iConfig.getUntrackedParameter<int>("debugHLT");
}

void GlobeHLT::defineBranch(TTree* tree) {
  
  hlt_p4  = new TClonesArray("TLorentzVector",   MAX_HLT);
  
  // Event Trigger
  tree->Branch("hlt1_bit_1", &hlt1_bit_1, "hlt1_bit_1/l");
  tree->Branch("hlt1_bit_2", &hlt1_bit_2, "hlt1_bit_2/l");
  tree->Branch("hlt1_bit_3", &hlt1_bit_3, "hlt1_bit_3/l");
  tree->Branch("hlt1_bit_4", &hlt1_bit_4, "hlt1_bit_4/l");
  tree->Branch("hlt1_bit_5", &hlt1_bit_5, "hlt1_bit_5/l");
  tree->Branch("hlt1_bit_6", &hlt1_bit_6, "hlt1_bit_6/l");
  
  tree->Branch("hlt1_bit_7", &hlt1_bit_7, "hlt1_bit_7/l");
  tree->Branch("hlt1_bit_8", &hlt1_bit_8, "hlt1_bit_8/l");
  tree->Branch("hlt1_bit_9", &hlt1_bit_9, "hlt1_bit_9/l");
  tree->Branch("hlt1_bit_10", &hlt1_bit_10, "hlt1_bit_10/l");

  //tree->Branch("hlt1_bit", &hlt1_bit, "hlt1_bit[6]/l");
  tree->Branch("hlt1_bit", &hlt1_bit, "hlt1_bit[10]/l");

  tree->Branch("hlt_path_names_HLT1_1", "std::vector<std::string>", &hlt_path_names_HLT1_1);
  tree->Branch("hlt_path_names_HLT1_2", "std::vector<std::string>", &hlt_path_names_HLT1_2);
  tree->Branch("hlt_path_names_HLT1_3", "std::vector<std::string>", &hlt_path_names_HLT1_3);
  tree->Branch("hlt_path_names_HLT1_4", "std::vector<std::string>", &hlt_path_names_HLT1_4);
  tree->Branch("hlt_path_names_HLT1_5", "std::vector<std::string>", &hlt_path_names_HLT1_5);
  tree->Branch("hlt_path_names_HLT1_6", "std::vector<std::string>", &hlt_path_names_HLT1_6);

  tree->Branch("hlt_path_names_HLT1_7", "std::vector<std::string>", &hlt_path_names_HLT1_7);
  tree->Branch("hlt_path_names_HLT1_8", "std::vector<std::string>", &hlt_path_names_HLT1_8);
  tree->Branch("hlt_path_names_HLT1_9", "std::vector<std::string>", &hlt_path_names_HLT1_9);
  tree->Branch("hlt_path_names_HLT1_10", "std::vector<std::string>", &hlt_path_names_HLT1_10);

  tree->Branch("hlt_path_names_HLT1", "std::vector<std::string>", &hlt_path_names_HLT1);

  ///prescales - 21 may
  tree->Branch("hlt_prescales_HLT1_1", "std::vector<int>", &hlt_prescales_HLT1_1);
  tree->Branch("hlt_prescales_HLT1_2", "std::vector<int>", &hlt_prescales_HLT1_2);
  tree->Branch("hlt_prescales_HLT1_3", "std::vector<int>", &hlt_prescales_HLT1_3);
  tree->Branch("hlt_prescales_HLT1_4", "std::vector<int>", &hlt_prescales_HLT1_4);
  tree->Branch("hlt_prescales_HLT1_5", "std::vector<int>", &hlt_prescales_HLT1_5);
  tree->Branch("hlt_prescales_HLT1_6", "std::vector<int>", &hlt_prescales_HLT1_6);

  tree->Branch("hlt_prescales_HLT1_7", "std::vector<int>", &hlt_prescales_HLT1_7);
  tree->Branch("hlt_prescales_HLT1_8", "std::vector<int>", &hlt_prescales_HLT1_8);
  tree->Branch("hlt_prescales_HLT1_9", "std::vector<int>", &hlt_prescales_HLT1_9);
  tree->Branch("hlt_prescales_HLT1_10", "std::vector<int>", &hlt_prescales_HLT1_10);

  //
  if(secondaryTriggerON){
    tree->Branch("hlt2_bit_1", &hlt2_bit_1, "hlt2_bit_1/l");
    tree->Branch("hlt2_bit_2", &hlt2_bit_2, "hlt2_bit_2/l");
    tree->Branch("hlt2_bit_3", &hlt2_bit_3, "hlt2_bit_3/l");
    tree->Branch("hlt2_bit_4", &hlt2_bit_4, "hlt2_bit_4/l");
    tree->Branch("hlt2_bit_5", &hlt2_bit_5, "hlt2_bit_5/l");
    tree->Branch("hlt2_bit_6", &hlt2_bit_6, "hlt2_bit_6/l");
    tree->Branch("hlt2_bit", &hlt2_bit, "hlt2_bit[6]/l");

    tree->Branch("hlt_path_names_HLT2_1", "std::vector<std::string>", &hlt_path_names_HLT2_1);
    tree->Branch("hlt_path_names_HLT2_2", "std::vector<std::string>", &hlt_path_names_HLT2_2);
    tree->Branch("hlt_path_names_HLT2_3", "std::vector<std::string>", &hlt_path_names_HLT2_3);
    tree->Branch("hlt_path_names_HLT2_4", "std::vector<std::string>", &hlt_path_names_HLT2_4);
    tree->Branch("hlt_path_names_HLT2_5", "std::vector<std::string>", &hlt_path_names_HLT2_5);
    tree->Branch("hlt_path_names_HLT2_6", "std::vector<std::string>", &hlt_path_names_HLT2_6);
    tree->Branch("hlt_path_names_HLT2", "std::vector<std::string>", &hlt_path_names_HLT2);
  }
  // Trigger Candidates
  tree->Branch("hlt_n", &hlt_n, "hlt_n/I");
  tree->Branch("hlt_p4", "TClonesArray", &hlt_p4, 32000, 0);
  tree->Branch("hlt_candpath_1"    , &hlt_candpath_1    , "hlt_candpath_1[hlt_n]/l"    );
  tree->Branch("hlt_candpath_2"    , &hlt_candpath_2    , "hlt_candpath_2[hlt_n]/l"    );
  tree->Branch("hlt_candpath_3"    , &hlt_candpath_3    , "hlt_candpath_3[hlt_n]/l"    );
  tree->Branch("hlt_candpath_4"    , &hlt_candpath_4    , "hlt_candpath_4[hlt_n]/l"    );
  tree->Branch("hlt_candpath_5"    , &hlt_candpath_5    , "hlt_candpath_5[hlt_n]/l"    );
  tree->Branch("hlt_candpath_6"    , &hlt_candpath_6    , "hlt_candpath_6[hlt_n]/l"    );
  
  tree->Branch("hlt_candpath_7"    , &hlt_candpath_7    , "hlt_candpath_7[hlt_n]/l"    );
  tree->Branch("hlt_candpath_8"    , &hlt_candpath_8    , "hlt_candpath_8[hlt_n]/l"    );
  tree->Branch("hlt_candpath_9"    , &hlt_candpath_9    , "hlt_candpath_9[hlt_n]/l"    );
  tree->Branch("hlt_candpath_10"    , &hlt_candpath_10    , "hlt_candpath_10[hlt_n]/l"    );
  
  //tree->Branch("hlt_candpath"    , &hlt_candpath    , "hlt_candpath[6][hlt_n]/l");
  tree->Branch("hlt_candpath"    , &hlt_candpath    , "hlt_candpath[10][hlt_n]/l");

  /*tree->Branch("hlt_id"            , &hlt_id            , "hlt_id[hlt_n]/I"            );
  tree->Branch("hlt_el_offlineind" , &hlt_el_offlineind , "hlt_el_offlineind[hlt_n]/I" );
  tree->Branch("hlt_ph_offlineind" , &hlt_ph_offlineind , "hlt_ph_offlineind[hlt_n]/I" );
  tree->Branch("hlt_mu_offlineind" , &hlt_mu_offlineind , "hlt_mu_offlineind[hlt_n]/I" );
  tree->Branch("hlt_jet_offlineind", &hlt_jet_offlineind, "hlt_jet_offlineind[hlt_n]/I");
  */  ///not working in myPlot.C - give seg violation
  //
  tree->Branch("hlt_label_names_1", "std::vector<std::string>", &hlt_label_names_1);
  tree->Branch("hlt_label_names_2", "std::vector<std::string>", &hlt_label_names_2);
  tree->Branch("hlt_label_names_3", "std::vector<std::string>", &hlt_label_names_3);
  tree->Branch("hlt_label_names_4", "std::vector<std::string>", &hlt_label_names_4);
  tree->Branch("hlt_label_names_5", "std::vector<std::string>", &hlt_label_names_5);
  tree->Branch("hlt_label_names_6", "std::vector<std::string>", &hlt_label_names_6);

  tree->Branch("hlt_label_names_7", "std::vector<std::string>", &hlt_label_names_7);
  tree->Branch("hlt_label_names_8", "std::vector<std::string>", &hlt_label_names_8);
  tree->Branch("hlt_label_names_9", "std::vector<std::string>", &hlt_label_names_9);
  tree->Branch("hlt_label_names_10", "std::vector<std::string>", &hlt_label_names_10);
  
}


bool GlobeHLT::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  bool changed = false;
  bool result = configProvider.init(iRun,iSetup,hlt1Tag_.process(),changed);
  if(debugHLT_){
    std::cout<<"result of configprovider:"<<result<<std::endl;
  }

  return 1;

}

bool GlobeHLT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  edm::Handle<reco::GsfElectronCollection> elH;
  iEvent.getByLabel(electronColl, elH);
  
  edm::Handle<reco::MuonCollection> muH;
  iEvent.getByLabel(muonColl, muH);
  
  edm::Handle<reco::PhotonCollection> phH;
  iEvent.getByLabel(photonCollStd, phH);
  
  edm::Handle<reco::CaloJetCollection> jetH;
  iEvent.getByLabel(jetColl, jetH);  
  
  edm::Handle<trigger::TriggerEvent> triggerObj;
  edm::Handle<trigger::TriggerEventWithRefs> triggerObjWithRef;
  
  iEvent.getByLabel(inputTag_, triggerObj);
  
  // HLT1
  hlt1_bit_1 = 0;
  hlt1_bit_2 = 0;
  hlt1_bit_3 = 0;
  hlt1_bit_4 = 0;
  hlt1_bit_5 = 0;
  hlt1_bit_6 = 0;

  hlt1_bit_7 = 0;
  hlt1_bit_8 = 0;
  hlt1_bit_9 = 0;
  hlt1_bit_10 = 0;

  ///moved to beginRun
  /*
  bool changed = false;
  bool result = configProvider.init(iEvent.getRun(),iSetup,hlt1Tag_.process(),changed);
  if(debugHLT_){
    std::cout<<"result of configprovider:"<<result<<std::endl;
  }
  */

  edm::Handle<edm::TriggerResults> h_triggerResults_HLT1;
  iEvent.getByLabel(hlt1Tag_, h_triggerResults_HLT1);
  if (h_triggerResults_HLT1.isValid()) {
    hlt_path_names_HLT1_1->clear();
    hlt_path_names_HLT1_2->clear();
    hlt_path_names_HLT1_3->clear();
    hlt_path_names_HLT1_4->clear();
    hlt_path_names_HLT1_5->clear();
    hlt_path_names_HLT1_6->clear();

    hlt_path_names_HLT1_7->clear();
    hlt_path_names_HLT1_8->clear();
    hlt_path_names_HLT1_9->clear();
    hlt_path_names_HLT1_10->clear();

    hlt_path_names_HLT1->clear();

    //////presclaes
    hlt_prescales_HLT1_1->clear();
    hlt_prescales_HLT1_2->clear();
    hlt_prescales_HLT1_3->clear();
    hlt_prescales_HLT1_4->clear();
    hlt_prescales_HLT1_5->clear();
    hlt_prescales_HLT1_6->clear();

    hlt_prescales_HLT1_7->clear();
    hlt_prescales_HLT1_8->clear();
    hlt_prescales_HLT1_9->clear();
    hlt_prescales_HLT1_10->clear();
    

    //Shilpi
    edm::TriggerNames const  &triggerNames_ = iEvent.triggerNames(*h_triggerResults_HLT1);

    if(debug_level > 9) std::cout << "Fill names HLT1" << std::endl;
    
    /////shilpi
    if(debugHLT_){
      std::cout<<"configProvider.size() = "<<configProvider.size()<<std::endl;
    }
    
    //hltprescale[iprescale] =  hltConfig_.prescaleValue(iEvent, iSetup, (*itr)) ; 


    for (size_t i = 0; i < configProvider.size(); ++i) {
      //int j=(int)(i/32);
      //Shilpi
      int j=(int)(i/64);
      hlt_path_names_HLT1->push_back(configProvider.triggerName(i));
      
      ///Shilpi//////////////////////////////////////////
      if(debugHLT_){
	std::cout<<"hltPath = "<<configProvider.triggerName(i)<<std::endl;
	//if(configProvider.triggerName(i) == "HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1" || configProvider.triggerName(i) =="HLT_Photon32_CaloIdL_Photon26_CaloIdL_v1")
	{
	  
	  int triggerIndex = triggerNames_.triggerIndex(configProvider.triggerName(i));
	  int lastModuleIndex = h_triggerResults_HLT1->index(triggerIndex);
	  std::cout<<"trigger name "<<configProvider.triggerName(i)<<" accept/ reject : "<<h_triggerResults_HLT1->accept(i)<<std::endl;
	  //from config provider
	  const std::vector<std::string>& mlabelConfig = configProvider.moduleLabels(triggerIndex);
	  for(int ii=0; ii<(int)mlabelConfig.size(); ii++)
	    {
	      std::cout<<"from COnfig PRovider...... mlable : "<<mlabelConfig[ii]<<std::endl;
	    }
	  
	  
	  for(int imodule=0;imodule<=lastModuleIndex;imodule++) {
	    const std::string moduleLabel = configProvider.moduleLabel(triggerIndex,imodule);
	    
	    std::cout<<"FOUND from trigger results : module label = "<<moduleLabel<<std::endl;
	    std::cout<<"index of this filter is "<<triggerObj->filterIndex(edm::InputTag(moduleLabel,"","HLT"))<<std::endl;
	  }
	
	}//if(hltPath == "HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1")
      }//if(debugHLT_) 
          

    ///////////END//////////////////////////
	
      if(debugHLT_)
	{
	
	  //if( i< h_triggerResults_HLT1->size())
	  //if( configProvider.triggerIndex(configProvider.triggerName(i))< h_triggerResults_HLT1->size() )
	  std::cout<<"prescale of trigger:"<<configProvider.triggerName(i)<<" is:"<<configProvider.prescaleValue(iEvent, iSetup, configProvider.triggerName(i))<<std::endl;
	  std::cout<<"i:"<<i<<" :index of trigger:"<<triggerNames_.triggerIndex(configProvider.triggerName(i))<<": size of hlt trigger results:"<<h_triggerResults_HLT1->size()<<std::endl;
	  std::cout<<"trigger index from config:"<<configProvider.triggerIndex(configProvider.triggerName(i))<<std::endl;
	  //std::cout<<"prescale set:"<<configProvider.prescaleSet(iEvent, iSetup)<<std::endl;
	}
    
      //if( i< h_triggerResults_HLT1->size()){ ////remove this - just a chk
      //if( configProvider.triggerIndex(configProvider.triggerName(i))< h_triggerResults_HLT1->size() ){ ////remove this - just a chk
      if(j<1) {     
        hlt_path_names_HLT1_1->push_back(configProvider.triggerName(i));
	hlt_prescales_HLT1_1->push_back(configProvider.prescaleValue(iEvent, iSetup, configProvider.triggerName(i)));
      } else if(j<2){
        hlt_path_names_HLT1_2->push_back(configProvider.triggerName(i));
	hlt_prescales_HLT1_2->push_back(configProvider.prescaleValue(iEvent, iSetup, configProvider.triggerName(i)));
      } else if(j<3){
        hlt_path_names_HLT1_3->push_back(configProvider.triggerName(i));
	hlt_prescales_HLT1_3->push_back(configProvider.prescaleValue(iEvent, iSetup, configProvider.triggerName(i)));
      } else if(j<4){
        hlt_path_names_HLT1_4->push_back(configProvider.triggerName(i));
	hlt_prescales_HLT1_4->push_back(configProvider.prescaleValue(iEvent, iSetup, configProvider.triggerName(i)));
      } else if(j<5){
        hlt_path_names_HLT1_5->push_back(configProvider.triggerName(i));
	hlt_prescales_HLT1_5->push_back(configProvider.prescaleValue(iEvent, iSetup, configProvider.triggerName(i)));
      } else if(j<6){
        hlt_path_names_HLT1_6->push_back(configProvider.triggerName(i));
	hlt_prescales_HLT1_6->push_back(configProvider.prescaleValue(iEvent, iSetup, configProvider.triggerName(i)));
      } else if(j<7){
        hlt_path_names_HLT1_7->push_back(configProvider.triggerName(i));
	hlt_prescales_HLT1_7->push_back(configProvider.prescaleValue(iEvent, iSetup, configProvider.triggerName(i)));
      } else if(j<8){
        hlt_path_names_HLT1_8->push_back(configProvider.triggerName(i));
	hlt_prescales_HLT1_8->push_back(configProvider.prescaleValue(iEvent, iSetup, configProvider.triggerName(i)));
      } else if(j<9){
        hlt_path_names_HLT1_9->push_back(configProvider.triggerName(i));
	hlt_prescales_HLT1_9->push_back(configProvider.prescaleValue(iEvent, iSetup, configProvider.triggerName(i)));
      } else if(j<10){
        hlt_path_names_HLT1_10->push_back(configProvider.triggerName(i));
	hlt_prescales_HLT1_10->push_back(configProvider.prescaleValue(iEvent, iSetup, configProvider.triggerName(i)));
      }
    }
    //}
  
    // Trigger Results
    if(debug_level > 99) std::cout << "### Trigger Results 1 :" << hlt1Tag_.process() << std::endl;
    for (size_t i = 0; i < configProvider.size(); ++i) {
      if(debug_level > 99) std::cout << i << "\t" << configProvider.triggerName(i) << " " << (h_triggerResults_HLT1->accept(i) ? "passed" : "failed") << std::endl;
      //int j=(int)(i/32);
      int j=(int)(i/64);

      if(h_triggerResults_HLT1->accept(i))
        set(hlt1_bit[j], i-64*j);

      if(j<1) {     
        if(h_triggerResults_HLT1->accept(i))
          set(hlt1_bit_1,i-64*j);
      } else if(j<2){
        if(h_triggerResults_HLT1->accept(i))
          set(hlt1_bit_2,i-64*j);
      } else if(j<3){
        if(h_triggerResults_HLT1->accept(i))
          set(hlt1_bit_3,i-64*j);
      } else if(j<4){
        if(h_triggerResults_HLT1->accept(i))
          set(hlt1_bit_4,i-64*j);
      } else if(j<5){
        if(h_triggerResults_HLT1->accept(i))
          set(hlt1_bit_5,i-64*j);
      }else if(j<6){
        if(h_triggerResults_HLT1->accept(i))
          set(hlt1_bit_6,i-64*j);
      } else if(j<7){
        if(h_triggerResults_HLT1->accept(i))
          set(hlt1_bit_7,i-64*j);
      } else if(j<8){
        if(h_triggerResults_HLT1->accept(i))
          set(hlt1_bit_8,i-64*j);
      } else if(j<9){
        if(h_triggerResults_HLT1->accept(i))
          set(hlt1_bit_9,i-64*j);
      } else if(j<10){
        if(h_triggerResults_HLT1->accept(i))
          set(hlt1_bit_10,i-64*j);
      }

    
      
      //
      if(debug_level > 999) {
        std::bitset<64> binary_hlt1(hlt1_bit_1);
        std::bitset<64> binary_hlt2(hlt1_bit_2);
        std::bitset<64> binary_hlt3(hlt1_bit_3);
        std::bitset<64> binary_hlt4(hlt1_bit_4);
        std::bitset<64> binary_hlt5(hlt1_bit_5);
        std::bitset<64> binary_hlt6(hlt1_bit_6);

        std::bitset<64> binary_hlt7(hlt1_bit_7);
        std::bitset<64> binary_hlt8(hlt1_bit_8);
        std::bitset<64> binary_hlt9(hlt1_bit_9);
        std::bitset<64> binary_hlt10(hlt1_bit_10);

        std::cout << "HLT Path 1: " << binary_hlt1  << " ? "<< check(hlt1_bit_1,i-64*j) << std::endl;
        std::cout << "HLT Path 2: " << binary_hlt2  << " ? "<< check(hlt1_bit_2,i-64*j) << std::endl;
        std::cout << "HLT Path 3: " << binary_hlt3  << " ? "<< check(hlt1_bit_3,i-64*j) << std::endl;
        std::cout << "HLT Path 4: " << binary_hlt4  << " ? "<< check(hlt1_bit_4,i-64*j) << std::endl;
        std::cout << "HLT Path 5: " << binary_hlt5  << " ? "<< check(hlt1_bit_5,i-64*j) << std::endl;
        std::cout << "HLT Path 6: " << binary_hlt6  << " ? "<< check(hlt1_bit_6,i-64*j) << std::endl;

	std::cout << "HLT Path 7: " << binary_hlt7  << " ? "<< check(hlt1_bit_7,i-64*j) << std::endl;
        std::cout << "HLT Path 8: " << binary_hlt8  << " ? "<< check(hlt1_bit_8,i-64*j) << std::endl;
        std::cout << "HLT Path 9: " << binary_hlt9  << " ? "<< check(hlt1_bit_9,i-64*j) << std::endl;
        std::cout << "HLT Path 10: " << binary_hlt10  << " ? "<< check(hlt1_bit_10,i-64*j) << std::endl;
      }
      //
    }
    if(debug_level > 99) std::cout << "\t Final result = " << hlt1_bit_1 << " " << hlt1_bit_2 << " " << hlt1_bit_3 << " " << hlt1_bit_4 << " " << hlt1_bit_5 << " " << hlt1_bit_6 << " " << hlt1_bit_7 << " " << hlt1_bit_8 << " " << hlt1_bit_9 << " " << hlt1_bit_10 <<std::endl;
  } else {
    if(debug_level > 9) std::cout << "TriggerResults not valid " << hlt1Tag_ << std::endl;
  }
  
  // HLT2
  if(secondaryTriggerON){
    hlt2_bit_1 = 0;
    hlt2_bit_2 = 0;
    hlt2_bit_3 = 0;
    hlt2_bit_4 = 0;
    hlt2_bit_5 = 0;
    hlt2_bit_6 = 0;
    bool changed = false;
    configProvider.init(iEvent.getRun(),iSetup,hlt2Tag_.process(),changed);
    edm::Handle<edm::TriggerResults> h_triggerResults_HLT2;
    iEvent.getByLabel(hlt2Tag_, h_triggerResults_HLT2);
    if (h_triggerResults_HLT2.isValid()) {
      hlt_path_names_HLT2_1->clear();
      hlt_path_names_HLT2_2->clear();
      hlt_path_names_HLT2_3->clear();
      hlt_path_names_HLT2_4->clear();
      hlt_path_names_HLT2_5->clear();
      hlt_path_names_HLT2_6->clear();

      hlt_path_names_HLT2->clear();
      
      if(debug_level > 9) std::cout << "Fill names HLT2" << std::endl;
      for (size_t i = 0; i < configProvider.size(); ++i) {
        //int j=(int)(i/32);
	int j=(int)(i/64);
        hlt_path_names_HLT2->push_back(configProvider.triggerName(i));
      
        if(j<1) {     
          hlt_path_names_HLT2_1->push_back(configProvider.triggerName(i));
        } else if(j<2){
          hlt_path_names_HLT2_2->push_back(configProvider.triggerName(i));
        } else if(j<3){
          hlt_path_names_HLT2_3->push_back(configProvider.triggerName(i));
        } else if(j<4){
          hlt_path_names_HLT2_4->push_back(configProvider.triggerName(i));
        }else if(j<5){
          hlt_path_names_HLT2_5->push_back(configProvider.triggerName(i));
        }else if(j<6){
          hlt_path_names_HLT2_6->push_back(configProvider.triggerName(i));
        }
      }
      // Trigger Results
      if(debug_level > 99) std::cout << "### Trigger Results 2: " << hlt2Tag_.process() << std::endl;
      for (size_t i = 0; i < configProvider.size(); ++i) {
        if(debug_level > 99) std::cout << i << "\t" << configProvider.triggerName(i) << " " << (h_triggerResults_HLT2->accept(i) ? "passed" : "failed") << std::endl;
        //int j=(int)(i/32);
	int j=(int)(i/64);

        if(h_triggerResults_HLT2->accept(i))
          set(hlt2_bit[j], i-64*j);

        if(j<1) {     
          if(h_triggerResults_HLT2->accept(i))
            set(hlt2_bit_1,i-64*j);
        } else if(j<2){
          if(h_triggerResults_HLT2->accept(i))
            set(hlt2_bit_2,i-64*j);
        } else if(j<3){
          if(h_triggerResults_HLT2->accept(i))
            set(hlt2_bit_3,i-64*j);
        } else if(j<4){
          if(h_triggerResults_HLT2->accept(i))
            set(hlt2_bit_4,i-64*j);
        } else if(j<5){
          if(h_triggerResults_HLT2->accept(i))
            set(hlt2_bit_5,i-64*j);
        } else if(j<6){
          if(h_triggerResults_HLT2->accept(i))
            set(hlt2_bit_6,i-64*j);
        }
        //
        if(debug_level > 999) {
          std::bitset<64> binary_hlt1(hlt2_bit_1);
          std::bitset<64> binary_hlt2(hlt2_bit_2);
          std::bitset<64> binary_hlt3(hlt2_bit_3);
          std::bitset<64> binary_hlt4(hlt2_bit_4);
          std::bitset<64> binary_hlt5(hlt2_bit_5);
          std::bitset<64> binary_hlt6(hlt2_bit_6);
          std::cout << "HLT Path 1: " << binary_hlt1  << " ? "<< check(hlt2_bit_1,i-32*j) << std::endl;
          std::cout << "HLT Path 2: " << binary_hlt2  << " ? "<< check(hlt2_bit_2,i-32*j) << std::endl;
          std::cout << "HLT Path 3: " << binary_hlt3  << " ? "<< check(hlt2_bit_3,i-32*j) << std::endl;
          std::cout << "HLT Path 4: " << binary_hlt4  << " ? "<< check(hlt2_bit_4,i-32*j) << std::endl;
          std::cout << "HLT Path 5: " << binary_hlt5  << " ? "<< check(hlt2_bit_5,i-64*j) << std::endl;
          std::cout << "HLT Path 6: " << binary_hlt6  << " ? "<< check(hlt2_bit_6,i-64*j) << std::endl;
        }
        //
      }
      if(debug_level > 99) std::cout << "\t Final result = " << hlt2_bit_1 << " " << hlt2_bit_2 << " " << hlt2_bit_3 << " " << hlt2_bit_4 << " " << hlt2_bit_5 << " " << hlt2_bit_6 << std::endl;
    } else {
      if(debug_level > 9) std::cout << "TriggerResults not valid " << hlt2Tag_ << std::endl;
    }
  }


  
  if(!triggerObj.isValid()) 
    throw(cms::Exception("Release Validation Error") << "RAW-type HLT results not found" );
  // Do it only the first event (when theHLTLabels is still empty)
  if(debug_level > 9) std::cout << "Trigger Paths" << std::endl;
  hlt_label_names_1->clear();
  hlt_label_names_2->clear();
  hlt_label_names_3->clear();
  hlt_label_names_4->clear();
  hlt_label_names_5->clear();
  hlt_label_names_6->clear();

  hlt_label_names_7->clear();
  hlt_label_names_8->clear();
  hlt_label_names_9->clear();
  hlt_label_names_10->clear();

  theHLTLabels.clear();
  
  //Shilpi
  if(debugHLT_){
    std::cout<<"triggerObj->sizeFilters() = "<<triggerObj->sizeFilters()<<std::endl;
  }


  for(int i=0; i<triggerObj->sizeFilters(); ++i) {
    if(debug_level > 9) std::cout << triggerObj->filterTag(i) << std::endl;
    theHLTLabels.push_back( triggerObj->filterTag(i));
    //Shilpi
    if(debugHLT_){
      std::cout<<"sizeFilters : hlt lable = "<<triggerObj->sizeFilters()<<" : "<<triggerObj->filterTag(i).label()<<std::endl;
    }

    //int j=(int)(i/32);
    int j=(int)(i/64);
    if(j<1) {
      hlt_label_names_1->push_back(triggerObj->filterTag(i).label());
    } else if(j<2) {
      hlt_label_names_2->push_back(triggerObj->filterTag(i).label());
    } else if(j<3) {
      hlt_label_names_3->push_back(triggerObj->filterTag(i).label());
    } else if(j<4) {
      hlt_label_names_4->push_back(triggerObj->filterTag(i).label());
    } else if(j<5) {
      hlt_label_names_5->push_back(triggerObj->filterTag(i).label());
    } else if(j<6) {
      hlt_label_names_6->push_back(triggerObj->filterTag(i).label());
    } else if(j<7) {
      hlt_label_names_7->push_back(triggerObj->filterTag(i).label());
    } else if(j<8) {
      hlt_label_names_8->push_back(triggerObj->filterTag(i).label());
    } else if(j<9) {
      hlt_label_names_9->push_back(triggerObj->filterTag(i).label());
    } else if(j<10) {
      hlt_label_names_10->push_back(triggerObj->filterTag(i).label());
    }
  }


  hlt_n = 0;
  hlt_p4->Clear();

  trigger::TriggerObjectCollection triggerObjs = triggerObj->getObjects();
  if(debug_level > 99) std::cout << "Trigger Objects found " << triggerObjs.size() << std::endl;
  
  for (unsigned int iCand=0; iCand<triggerObjs.size(); ++iCand ) {
    if(debug_level > 99) std::cout << iCand << "=" << hlt_n << std::endl;
    if (hlt_n >= MAX_HLT) {
      std::cout << "GlobeHLT: WARNING TOO MANY HLT CANDIDATES: (allowed " << MAX_HLT << ")" << std::endl;
      break;
    }
    
    trigger::TriggerObject object = triggerObjs[iCand];
    hlt_candpath_1[hlt_n] = 0;
    hlt_candpath_2[hlt_n] = 0;
    hlt_candpath_3[hlt_n] = 0;
    hlt_candpath_4[hlt_n] = 0;
    hlt_candpath_5[hlt_n] = 0;
    hlt_candpath_6[hlt_n] = 0;

    hlt_candpath_7[hlt_n] = 0;
    hlt_candpath_8[hlt_n] = 0;
    hlt_candpath_9[hlt_n] = 0;
    hlt_candpath_10[hlt_n] = 0;
    
    for (int i=0; i<10; i++) 
      hlt_candpath[i][hlt_n] = 0;
    
    for(unsigned int n=0; n<theHLTLabels.size(); n++) {
      
      trigger::size_type elIndex = triggerObj->filterIndex(theHLTLabels[n]);
      ///shilpi
      if(debugHLT_){
	std::cout<<"obj:::::iCand : hltlable : "<<iCand<<" :" <<theHLTLabels[n].label()<<std::endl;
      }
      
      // Check HLT
      bool firedHLT = false;
      if (!(elIndex >= triggerObj->sizeFilters())) {
        const trigger::Keys & k = triggerObj->filterKeys(elIndex);
        for (trigger::Keys::const_iterator ki = k.begin(); ki !=k.end(); ++ki ) {
          if(*ki == iCand)
            firedHLT = true;
        }
      }
    
      
      if(firedHLT) {
        if(debug_level > 99) std::cout << n << "\t" << theHLTLabels[n].label() << "\t fired" << std::endl;
	///Shilpi
	if(debugHLT_){
	  std::cout<<"PASSED HLT obj:::::iCand : hltlable"<<iCand<<" :" <<theHLTLabels[n].label()<<std::endl;
	}

        //int j=(int)(n/32);
	int j=(int)(n/64);
	
        for (unsigned int i=0; i<configProvider.size(); i++) {
	  
          unsigned int nModules = configProvider.moduleLabels(i).size();
          unsigned int moduleIndex = configProvider.moduleIndex(i, theHLTLabels[n].label());
          //for(unsigned int y=0; y<nModules; y++) {
          //if (configProvider.moduleLabel(i, y) == theHLTLabels[n].label()) 
          //    std::cout << nModules << " " << y << " " << theHLTLabels[n].label() << " " << configProvider.triggerName(i) << std::endl;
          //}
          
          //std::cout <<  theHLTLabels[n].label() << " ";
          //std::cout << configProvider.moduleLabel(i, moduleIndex-1) << " " << nModules << " " << moduleIndex << std::endl;
          if (moduleIndex == (nModules-2)) {
            //int j=(int)(i/32);
	    int j=(int)(i/64);
            //std::cout << iEvent.id().event() << " " << configProvider.moduleLabel(i, moduleIndex) << " ";
            //std::cout <<  theHLTLabels[n].label() << " ";
            //std::cout << j << " " << i-32*j << " " <<  configProvider.triggerName(i) << std::endl;
            //std::cout << i-32*j << std::endl;
            //std::cout << std::endl;
            //set(hlt_candpath[j][hlt_n], i-32*j);
	    set(hlt_candpath[j][hlt_n], i-64*j);
          }
        }
     
	//SHILPI
	//std::cout<<"n : j : n-64*j "<<n<<" : "<<j<<" : "<<n-64*j<<std::endl;
        if(j<1) {
	  //std::cout<<"setting 1 "<<std::endl;
          set(hlt_candpath_1[hlt_n], n-64*j);
        } else if(j<2) {
	  //std::cout<<"setting 2 "<<std::endl;
	  //std::cout<<"BEFORE CALLING SET, hlt_candpath_2["<<hlt_n<<"] "<<hlt_candpath_2[hlt_n]<<std::endl;
          set(hlt_candpath_2[hlt_n], n-64*j);
	  //std::cout<<"hlt_candpath_2["<<hlt_n<<"] "<<hlt_candpath_2[hlt_n]<<std::endl;
	  //std::cout<<"hlt_candpath_2["<<hlt_n<<"] >> (n-64*j) & 1 : "<<(hlt_candpath_2[hlt_n] >> (n-64*j) & 1)<<std::endl;
        } else if(j<3) {
	  //std::cout<<"setting 3 "<<std::endl;
          set(hlt_candpath_3[hlt_n], n-64*j);
        } else if(j<4) {
	  //std::cout<<"setting 4 "<<std::endl;
          set(hlt_candpath_4[hlt_n], n-64*j);
        } else if(j<5) {
	  //std::cout<<"setting 5 "<<std::endl;
          set(hlt_candpath_5[hlt_n], n-64*j);
        } else if(j<6) {
	  //std::cout<<"setting 6 "<<std::endl;
          set(hlt_candpath_6[hlt_n], n-64*j);
        } else if(j<7) {
	  set(hlt_candpath_7[hlt_n], n-64*j);
        } else if(j<8) {
	  set(hlt_candpath_8[hlt_n], n-64*j);
        } else if(j<9) {
	  set(hlt_candpath_9[hlt_n], n-64*j);
        } else if(j<10) {
	  set(hlt_candpath_10[hlt_n], n-64*j);
        } 
	
      }//if(firedHLT) 
      
    }//for(unsigned int n=0; n<theHLTLabels.size(); n++)
    
    
    // Skip if no triggers were fired
    if(hlt_candpath_1[hlt_n]!=0 || hlt_candpath_2[hlt_n]!=0 || hlt_candpath_3[hlt_n]!=0 || hlt_candpath_4[hlt_n]!=0 || hlt_candpath_5[hlt_n]!=0  || hlt_candpath_6[hlt_n]!=0 || hlt_candpath_7[hlt_n]!=0 || hlt_candpath_8[hlt_n]!=0 || hlt_candpath_9[hlt_n]!=0 || hlt_candpath_10[hlt_n]!=0) {
      
        // Set variables
      hlt_id[hlt_n] = object.id();
    
      // Set HLT candidate p4
      TLorentzVector lv(object.px(), object.py(), object.pz(), 0);
      new ((*hlt_p4)[hlt_n]) TLorentzVector();
      ((TLorentzVector *)hlt_p4->At(hlt_n))->SetXYZT(object.px(), object.py(), object.pz(), object.energy());
      
      // OFFLINE matchings           
      /*float pt = object.pt(); 
      float eta = object.eta(); 
      float phi = object.phi(); 
      TVector3 v1;
      v1.SetPtEtaPhi(pt, eta, phi);
      // ELECTRON Matching
      float dRmin = 1.0;
      int index = -1;
      for(unsigned int z=0; z<elH->size(); z++) {
        reco::GsfElectronRef e(elH, z);
        // apply the cuts
        if(gCUT->cut(*e))
          continue;
        // passed cuts
        TVector3 v2;
        v2.SetPtEtaPhi(e->pt(), e->eta(), e->phi());
        float dR = v1.DeltaR(v2); 
        if (dR < dRmin) {
          dRmin = dR;
          index = z;
        }
      }
      hlt_el_offlineind[hlt_n] = index;
      
      // PHOTON matching
      dRmin = 1.0;
      index = -1;
      for(unsigned int z=0; z<phH->size(); z++) {
        reco::PhotonRef e(phH, z);
        // apply the cuts
        if(gCUT->cut(*e))
          continue;
        // passed cuts
        TVector3 v2;
        v2.SetPtEtaPhi(e->pt(), e->eta(), e->phi());
        float dR = v1.DeltaR(v2); 
        if (dR < dRmin) {
          dRmin = dR;
          index = z;
        }
      }
      hlt_ph_offlineind[hlt_n] = index;
      
      // MUON matching
      dRmin = 1.0;
      index = -1;
      for(unsigned int z=0; z<muH->size(); z++) {
        reco::MuonRef e(muH, z);
        // apply the cuts
        if(gCUT->cut(*e))
          continue;
        // passed cuts
        TVector3 v2;
        v2.SetPtEtaPhi(e->pt(), e->eta(), e->phi());
        float dR = v1.DeltaR(v2); 
        if (dR < dRmin) {
          dRmin = dR;
          index = z;
        }
      }
      hlt_mu_offlineind[hlt_n] = index; 
      
      // JET matching
      dRmin = 1.0;
      index = -1;
      for(unsigned int z=0; z<jetH->size(); z++) {
        reco::CaloJetRef e(jetH, z);
        // apply the cuts
        if(gCUT->cut(*e))
          continue;
        // passed cuts
        TVector3 v2;
        v2.SetPtEtaPhi(e->pt(), e->eta(), e->phi());
        float dR = v1.DeltaR(v2); 
        if (dR < dRmin) {
          dRmin = dR;
          index = z;
        }
      }
      hlt_jet_offlineind[hlt_n] = index; 
      
      //
      if(debug_level > 99) {
        std::bitset<64> binary_hlt1(hlt_candpath_1[hlt_n]);
        std::bitset<64> binary_hlt2(hlt_candpath_2[hlt_n]);
        std::bitset<64> binary_hlt3(hlt_candpath_3[hlt_n]);
        std::bitset<64> binary_hlt4(hlt_candpath_4[hlt_n]);
        std::bitset<64> binary_hlt5(hlt_candpath_5[hlt_n]);
        std::bitset<64> binary_hlt6(hlt_candpath_6[hlt_n]);
        std::cout << iCand << "=" << hlt_n
                  <<"\t pT=" << object.pt() << "\t eta=" << object.eta() << "\t particle " << object.id() << std::endl
                  << "\t candpath1=" << binary_hlt1 << "\t candpath2=" << binary_hlt2
                  << "\t candpath3=" << binary_hlt3 << "\t candpath4=" << binary_hlt4
                  << "\t candpath5=" << binary_hlt5 << "\t candpath6=" << binary_hlt6
                  << "\t Matching"
                  << "\t electron " << hlt_el_offlineind[ hlt_n]
                  << "\t photon "   << hlt_ph_offlineind[ hlt_n]
                  << "\t muon "     << hlt_mu_offlineind[ hlt_n]
                  << "\t jet "      << hlt_jet_offlineind[hlt_n]
                  << std::endl;
		  }*/
    //
      hlt_n++; 
    
    } // Store candidate which fired at least 1 HLT
    

    
  } // TriggerCandidate's Loop
  

  ///Shilpi
  /*for(int ii=0; ii<hlt_n; ii++)
    {
      std::cout<<"Now Looking at the object "<<std::endl;
      std::cout<<"hlt_candpath_1["<<ii<<"] "<<hlt_candpath_1[ii]<<std::endl;
      std::cout<<"hlt_candpath_2["<<ii<<"] "<<hlt_candpath_2[ii]<<std::endl;
      std::cout<<"hlt_candpath_2["<<ii<<"] >>47 &1"<<(hlt_candpath_2[ii] >> 47 & 1)<<std::endl;
      std::cout<<"hlt_candpath_3["<<ii<<"] "<<hlt_candpath_3[ii]<<std::endl;
      std::cout<<"hlt_candpath_4["<<ii<<"] "<<hlt_candpath_4[ii]<<std::endl;
      std::cout<<"hlt_candpath_5["<<ii<<"] "<<hlt_candpath_5[ii]<<std::endl;
      std::cout<<"hlt_candpath_6["<<ii<<"] "<<hlt_candpath_6[ii]<<std::endl;
    }
  */

  if(debug_level > 99) std::cout << "Trigger Objects stored " << hlt_n << std::endl;
  
  return true;
}

// USEFUL LINKS
// https://twiki.cern.ch/twiki/bin/view/CMS/TSG_13_VIII_09_8E29
// http://cms-project-confdb-hltdev.web.cern.ch/cms-project-confdb-hltdev/browser/convert2Html.jsp?dbName=HLTDEV&configName=/CMSSW/CMSSW_3_1_2/8E29/V1
// https://twiki.cern.ch/twiki/bin/view/CMS/TSG_13_VIII_09_1E31
// http://cms-project-confdb-hltdev.web.cern.ch/cms-project-confdb-hltdev/browser/convert2Html.jsp?dbName=HLTDEV&configName=/CMSSW/CMSSW_3_1_2/1E31/V1
// http://cms-service-sdtweb.web.cern.ch/cms-service-sdtweb/doxygen/CMSSW_3_2_5/doc/html/da/d01/TriggerTypeDefs_8h-source.html
//  8E29 Trigger Tables
// https://twiki.cern.ch/twiki/bin/view/CMS/TriggerMenuDescription8E29Devel
//  1E31 Trigger Tables
// https://twiki.cern.ch/twiki/bin/view/CMS/TriggerMenuDescription1E31Devel
