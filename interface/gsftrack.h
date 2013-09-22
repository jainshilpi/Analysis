#include "TLorentzVector.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

class gsftrack: public TObject{
 public:
  //methods
  gsftrack();
  ~gsftrack();

  //variables
  int charge;
  int numberOfLostHits;
  int numberOfValidHits;
  int lost;
  int recHitsSize;
  int qualityMask;


  double pt;
  double eta;
  double chi2;
  double d0;
  double dxy;
  double dz;
  double ptin;
  double ptout;
  double fbrem;
  double qoverp;
  double vx;
  double vy;
  double vz;
  double phi;
  double theta;
  double ndof;
  double outerX;
  double outerY;
  double outerZ;
  double outerRadius;
  double innerX;
  double innerY;
  double innerZ;

  
};

gsftrack::gsftrack()
{}
gsftrack::~gsftrack()
{}
