#ifndef StHadronTrack_hh
#define StHadronTrack_hh

#include <cmath>

class StPicoTrack;
class StPicoDst;
class StPicoBTofPidTraits;
class StPicoEmcPidTraits;
class StDcaGeometry;

#include "TObject.h"
#include "StThreeVectorF.hh"
#include "TVector3.h"
#include <stdio.h>
#include <math.h>
#include "StEvent/StDcaGeometry.h"

// Macro to control EMC variables
#define EMCON 1
class StPicoEvent;

class StHadronTrack : public TObject {
 public:
  StHadronTrack();
  ~StHadronTrack();
  StHadronTrack(StPicoDst *picoDst, StPicoTrack *t, Int_t idx);
  virtual void Print(const Char_t *option = "") const;  ///< Print track info
            
  Int_t   id() const             { return (Int_t)mId; }
  //StThreeVectorF gMom() const;
  StThreeVectorF pMom() const    { return mPMom; }
  Short_t charge() const         { return (mNHitsFit>0) ? +1 : -1; }
  Int_t   nHitsFit() const       { return (mNHitsFit>0) ? (Int_t)mNHitsFit : (Int_t)(-1*mNHitsFit); }
  Int_t   nHitsDedx() const      { return (Int_t)mNHitsDedx; }
  //Int_t   nHitsMax() const       { return (Int_t)mNHitsMax; }
  Float_t nSigmaPion() const     { return (Float_t)mNSigmaPion/1000.; }
  Float_t nSigmaKaon() const     { return (Float_t)mNSigmaKaon/1000.; }
  Float_t dca() const           { return (Float_t)abs(mDca)/10000.; }
  //Float_t dcaXY() const            { return (Float_t)mDcaXY/10000.; }
  //Float_t dcaZ() const            { return (Float_t)mDcaZ/10000.; }
  Float_t beta() const           { return (Float_t)mBeta/20000.; }
  Float_t localY() const          {return (Float_t)mLocalY/1000.;}
  Bool_t isHFTTrack() const { return mDca<0?true:false; }
  Int_t towerId() const          {return mBTOWId;}
  Short_t adc0() const          {return mBTOWADC0;}
  Int_t emcTriggerId() const          {return mEmcTrgId;}

 protected:
  Short_t mId;               // track index in picoDst
  //Short_t  mDcaXY;              // dcaXY * 10000
  //Short_t  mDcaZ;              // dcaZ * 10000
  StThreeVectorF mPMom;
  //StThreeVectorF mGMom;

 // UShort_t  mGPt;             //global pt * 1000
 // Short_t  mGEta;            //global eta * 10000
 // Short_t  mGPhi;            //global phi * 10000

  // pidTraits
  UShort_t mBeta;  // *20000 
  Short_t  mLocalY; // *1000

  Short_t  mDca;              // dca * 10000 * (isHFT?-1:1)
  Char_t   mNHitsFit;         // q*nHitsFit
  //Char_t   mNHitsMax;         // nHitsMax - TPC
  UChar_t  mNHitsDedx;        // nHitsDedx
  Short_t  mNSigmaPion;       // nsigmaPi * 1000
  Short_t  mNSigmaKaon;       // nsigmaKaon * 1000

  // From Stadard BEMC Clustering Algorithm or Single Tower ID
  Short_t  mBTOWADC0;         // adc0 higest adc in the cluster
  Short_t  mBTOWId;           // projected tower Id 1-4800
  Short_t  mEmcTrgId;
  
  friend class StPicoDst;

  ClassDef(StHadronTrack, 1)
};

/*inline StThreeVectorF StHadronTrack::gMom() const{
   TVector3 tmp(0,0,0);
   tmp.SetPtEtaPhi(mGPt/1000.,mGEta/10000.,mGPhi/10000.);
   return StThreeVectorF(tmp.X(),tmp.Y(),tmp.Z());
}*/
#endif
