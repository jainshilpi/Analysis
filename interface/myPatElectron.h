#include "TLorentzVector.h"


class myPatElectron: public pat::Electron{
 public:
  //methods
  myPatElectron();
  ~myPatElectron();

  //ID flags 1(true) ; 0(false)
  int elecutbits;
  int heepcutbit;
  /*  int elepassedet;
  int elepassedecalDrivenSeed;
  int elepassedhadronicOverEm;
  int elepassede1x5Bye5x5;
  int elepassedsigmaIetaIeta;
  int eleppaseddeltaEtain;
  int eleppaseddeltaPhiin;
  int elepasseddr03IsolEMHadDepth1;
  int elepassedIsoHadDepth2;
  int elepasseddr03TkSumPt;
  */
};

myPatElectron::myPatElectron()
{}
myPatElectron::~myPatElectron()
{}
