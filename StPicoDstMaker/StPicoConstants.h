#ifndef StPicoConstants_h
#define StPicoConstants_h
#include "Rtypes.h"

  enum {
     nTrk     = 50000,
  };


class Pico {
public:
  // event selectioin
  static Float_t mVzMax;// = 40;
  static Float_t mVrMax;// = 2;
  static Int_t mRefMultMin;// = 8;                    // >=

  // track selection
  static Float_t mPtMin;// = 0.1;
  static Int_t mNHitsFitMin;// = 15;                 // >=
  static Float_t mRatioMin;// = 0; //0.52;  -> removed for Run14 AuAu200GeV
  static Float_t mGDcaMax;// = 10.;

  static UShort_t USHORTMAX;// = 65535;
  static UChar_t UCHARMAX;  // = 255;
  ClassDef(Pico, 1)
};
#endif
