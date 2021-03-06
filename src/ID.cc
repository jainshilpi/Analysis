#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include <iostream>
#include <vector> 

#include "Analyzer/Analyser/interface/HeepCutbits.h"
#include "Analyzer/Analyser/interface/TightIDphocutbits.h"
//#include "/afs/cern.ch/user/s/shilpi/scratch0/CMSSW_3_7_0_patch4/src/Analyzer/Analyser/interface/myPatElectron.h"


int heepeleID( std::vector<pat::Electron> &patelectron_container,int ipatel ){
  
  using namespace std;
  
  double pateleet                            = patelectron_container[ipatel].et();
  double patelept                            = patelectron_container[ipatel].pt();

  //sceta
  reco::SuperClusterRef scref = patelectron_container[ipatel].superCluster();
  double patelesceta                         = scref->eta();
  
  int pateleecalDrivenSeed                = patelectron_container[ipatel].ecalDrivenSeed();
  double pateledr03EcalRecHitSumEt           = patelectron_container[ipatel].dr03EcalRecHitSumEt();
  double pateledr03HcalDepth1TowerSumEt      = patelectron_container[ipatel].dr03HcalDepth1TowerSumEt();
  double pateleIsoEMHad1                     = pateledr03EcalRecHitSumEt + pateledr03HcalDepth1TowerSumEt;
 
  double pateledr03HcalDepth2TowerSumEt      = patelectron_container[ipatel].dr03HcalDepth2TowerSumEt();
  double pateledr03HcalTowerSumEt            = patelectron_container[ipatel].dr03HcalTowerSumEt();
  double patelee1x5                          = patelectron_container[ipatel].e1x5();
  double patelee2x5Max                       = patelectron_container[ipatel].e2x5Max();
  double patelee5x5                          = patelectron_container[ipatel].e5x5();
  double patelehadronicOverEm                = patelectron_container[ipatel].hadronicOverEm();
  double patelesigmaIetaIeta                 = patelectron_container[ipatel].sigmaIetaIeta();
  double pateledr03TkSumPt                   = patelectron_container[ipatel].dr03TkSumPt();
         
  double pateledeltaEtaSuperClusterTrackAtVtx= patelectron_container[ipatel].deltaEtaSuperClusterTrackAtVtx();
  double pateledeltaPhiSuperClusterTrackAtVtx= patelectron_container[ipatel].deltaPhiSuperClusterTrackAtVtx();
  
  //int heepcutbit;
  int heepcutbit = ~0x0;

  //std::cout<<""<<std::endl;
  //std::cout<<"ipatel = "<<ipatel<<std::endl;
  //std::cout<<"eta = "<<patelectron_container[ipatel].eta()<<std::endl;
  
  if( !(patelectron_container[ipatel].isEB()) && !(patelectron_container[ipatel].isEE()) )
    heepcutbit = 0x0;

  if( patelectron_container[ipatel].isEB() )
    {
      
      //std::cout<<"its barrel electron "<<std::endl;
      //et
      if( pateleet < 25.0 ) 
	{
	  //cout<<"failed et "<<endl;
	  heepcutbit ^= et;
	  //cout<<"heepcut bit in ID.cc = "<<heepcutbit<<endl;
	}
      //isEcalDriven
      if( !pateleecalDrivenSeed )
	{
	  //cout<<"failed ecaldriven "<<endl;
	  heepcutbit ^= ecalDrivenSeed;
	  //cout<<"heepcut bit in ID.cc = "<<heepcutbit<<endl;
	}
      //H/E
      if( patelehadronicOverEm > 0.05)
	{
	  //cout<<"failed H/E "<<endl;
	  heepcutbit ^= hadronicOverEm;
	  //cout<<"heepcut bit in ID.cc = "<<heepcutbit<<endl;
	}
      //E2x5/E5x5 || E1x5/E5x5
      double e2x5Bye5x5 = patelee2x5Max/patelee5x5;
      double e1x5Bye5x5 = patelee1x5/patelee5x5;
      
      if( !(e2x5Bye5x5 > 0.94 || e1x5Bye5x5 > 0.83) )
	{
	  //cout<<"failed ratiocut "<<endl;
	  heepcutbit ^= ratioto5x5;
	  //cout<<"heepcut bit in ID.cc = "<<heepcutbit<<endl;
	}
      //sigmaieta -----NOT in EB
      heepcutbit |= sigmaIetaIeta;
      //cout<<"heepcut bit in ID.cc = "<<heepcutbit<<endl;

      //deltaEtain
      if( !(fabs(pateledeltaEtaSuperClusterTrackAtVtx)<0.005) )
	{
	  //cout<<"failed delEtain"<<endl;
	  heepcutbit ^= eleppaseddeltaEtain;
	  //cout<<"heepcut bit in ID.cc = "<<heepcutbit<<endl;
	}
      //deltaPhiin
      if( !(fabs(pateledeltaPhiSuperClusterTrackAtVtx)<0.09) )
	{
	  //cout<<"failed delPhiin"<<endl;
	  heepcutbit ^= eleppaseddeltaPhiin;
	  //cout<<"heepcut bit in ID.cc = "<<heepcutbit<<endl;
	}
      //isoEMHad
      if( !(pateleIsoEMHad1<2+0.03*pateleet) )
	{
	  //cout<<"failed EM"<<endl;
	  heepcutbit ^= dr03IsolEMHadDepth1;
	  //cout<<"heepcut bit in ID.cc = "<<heepcutbit<<endl;
	}
      //haddepth2 ------- NOT applicable in EB
      //if( pateledr03HcalDepth2TowerSumEt )
      heepcutbit |= IsoHadDepth2;
      //cout<<"heepcut bit in ID.cc = "<<heepcutbit<<endl;

      //trkIso
      if( pateledr03TkSumPt>7.5 )
	{
	  //cout<<"failed trkIso"<<endl;
	  heepcutbit ^= dr03TkSumPt;
	  //cout<<"heepcut bit in ID.cc = "<<heepcutbit<<endl;
	}
    }//if( patelectron_container[ipatel].isEB() )
  
  if( patelectron_container[ipatel].isEE() )
    {
  
      //std::cout<<"its endcap electron"<<std::endl;
      //et
      if( pateleet < 25.0 ) 
	heepcutbit ^= et;

      //isEcalDriven
      if( !pateleecalDrivenSeed )
	heepcutbit ^= ecalDrivenSeed;
      
      //H/E
      if( patelehadronicOverEm > 0.05)
	heepcutbit ^= hadronicOverEm;
      
      //E2x5/E5x5 || E1x5/E5x5
      //double e2x5Bye5x5 = patelee2x5Max/patelee5x5;
      //double e1x5Bye5x5 = patelee1x5/patelee5x5;
      
      //if( e2x5Bye5x5 > 0.94 || e1x5Bye5x5 > 0.83 ) ------- NOT applicable in EE
      heepcutbit |= ratioto5x5;
      
      //sigmaieta -----NOT in EB
      if( patelesigmaIetaIeta>0.03 )
	heepcutbit ^= sigmaIetaIeta;
      
      //deltaEtain
      if( !(fabs(pateledeltaEtaSuperClusterTrackAtVtx)<0.007) )
	heepcutbit ^= eleppaseddeltaEtain;
      
      //deltaPhiin
      if( !(fabs(pateledeltaPhiSuperClusterTrackAtVtx)<0.09) )
	heepcutbit ^= eleppaseddeltaPhiin;
      
      //isoEMHad
      double cut;
      if( pateleet<50 )
	cut = 2.5;
      else
	cut = 2.5+0.03*(pateleet - 50.0 );
      
      if( pateleIsoEMHad1>cut )
	heepcutbit ^= dr03IsolEMHadDepth1;
      
      //haddepth2 ------- NOT applicable in EB
      if( pateledr03HcalDepth2TowerSumEt>0.5 )
	heepcutbit ^= IsoHadDepth2;
      
      //trkIso
      if( pateledr03TkSumPt>15.0 )
	heepcutbit ^= dr03TkSumPt;
    }//if( patelectron_container[ipatel].isEE() )

  //heepcutbit = heepcutbit | exbit1 | exbit2 | exbit3 | exbit4 | exbit5 | exbit6;
  
  //std::cout<<"heepcut bit in ID.cc = "<<heepcutbit<<std::endl;
  return heepcutbit;
  
}//int heepeleID( std::vector<pat::Electron> &patelectron_container,int ipatel )


int tightphoID( std::vector<pat::Photon> &patphoton_container,int ipatpho ){
  
  double patphopt                                = patphoton_container[ipatpho].pt();
  double patphoecalRecHitSumEtConeDR04           = patphoton_container[ipatpho].ecalRecHitSumEtConeDR04();
  double patphohcalTowerSumEtConeDR04            = patphoton_container[ipatpho].hcalTowerSumEtConeDR04();
  double patphohadronicOverEm                    = patphoton_container[ipatpho].hadronicOverEm();
  double patphotrkSumPtHollowConeDR04            = patphoton_container[ipatpho].trkSumPtHollowConeDR04();
  double patphosigmaIetaIeta                     = patphoton_container[ipatpho].sigmaIetaIeta();
  int patphohasPixelSeed                         = patphoton_container[ipatpho].hasPixelSeed();

  int phoidcutbit = ~0x0;
  
  if( !(patphoton_container[ipatpho].isEB()) && !(patphoton_container[ipatpho].isEE()) )
    phoidcutbit = 0x0;


  if( patphoton_container[ipatpho].isEB() )
    {
      if( patphoecalRecHitSumEtConeDR04>4.2+0.003*patphopt )
	phoidcutbit ^= ecalrechitSum04;
      
      if( patphohcalTowerSumEtConeDR04>2.2+0.001*patphopt )
	phoidcutbit ^= hcaltowerSum04;

      if( patphohadronicOverEm>0.05 )
	phoidcutbit ^= phohadronicOverEm;

      if( patphotrkSumPtHollowConeDR04>2.0+0.001*patphopt )
	phoidcutbit ^= trkSumptHollow04;

      if( patphosigmaIetaIeta>0.013 )
	phoidcutbit ^= phosigmaIetaIeta;

      if(patphohasPixelSeed)
	phoidcutbit ^= haspixelseed;
    }//if( patphoton_container[ipatpho].isEB() )

  if( patphoton_container[ipatpho].isEE() )
    {
      if( patphoecalRecHitSumEtConeDR04>4.2+0.006*patphopt )
        phoidcutbit ^= ecalrechitSum04;

      if( patphohcalTowerSumEtConeDR04>2.2+0.0025*patphopt )
        phoidcutbit ^= hcaltowerSum04;

      if( patphohadronicOverEm>0.05 )
        phoidcutbit ^= phohadronicOverEm;

      if( patphotrkSumPtHollowConeDR04>3.5+0.001*patphopt )
        phoidcutbit ^= trkSumptHollow04;

      if( patphosigmaIetaIeta>0.03 ) ///-------- now applied in EE
        phoidcutbit ^= phosigmaIetaIeta;

      if(patphohasPixelSeed)
        phoidcutbit ^= haspixelseed;

      /*//if( patphosigmaIetaIeta>0.03 ) ///-------- not applied in EE
        phoidcutbit |= phosigmaIetaIeta;
      */
    }//if( patphoton_container[ipatpho].isEB() )

  return phoidcutbit;
}
