#include "TLorentzVector.h"

class myrecophoclass: public TObject{
 public:
  // methods 
  myrecophoclass();
  ~myrecophoclass();
  
  double charged;
  double photon;
  double neutral;
  /////teh ones for recommende ID
  double recomencharged;
  double recomenphoton;
  double recomenneutral;

  double pt;
  int passelectronveto;

  ////for worst vertex isolation
  double PFphotonWorstChargedHadronIso;
};

myrecophoclass::myrecophoclass()
{}
myrecophoclass::~myrecophoclass()
{}
