#ifndef StElectronTrack_hh
#define StElectronTrack_hh

#include <cmath>

class StPicoTrack;
class StPicoDst;
class StPicoBTofPidTraits;
class StPicoEmcPidTraits;
class StDcaGeometry;

#include "TObject.h"
#include "StThreeVectorF.hh"
#include "TVector2.h"
#include <stdio.h>
#include <math.h>
#include "StEvent/StDcaGeometry.h"
#include "PhysicalConstants.h"

class StElectronTrack : public TObject {
 public:
  StElectronTrack();
  ~StElectronTrack();
  StElectronTrack(StPicoDst *picoDst, StPicoTrack *t, Int_t idx);
  virtual void Print(const Char_t *option = "") const;  ///< Print track info
            
  Int_t   id() const             { return (Int_t)mId; }
  //Float_t chi2() const           { return (Float_t)mChi2/1000.; }
  Float_t gPt() const;
  Float_t gEta() const;
  Float_t gPhi() const;
  StThreeVectorF gMom() const;
  StThreeVectorF const & pMom() const    { return mPMom; }
  Short_t charge() const         { return (mNHitsFit>0) ? +1 : -1; }
  Int_t   nHitsFit() const       { return (mNHitsFit>0) ? (Int_t)mNHitsFit : (Int_t)(-1*mNHitsFit); }
  //Int_t   nHitsMax() const       { return (Int_t)mNHitsMax; }
  Int_t   nHitsDedx() const      { return (Int_t)mNHitsDedx; }
  //Int_t   nHitsMapHFT() const    { return (Int_t)mNHitsMapHFT; }
  //Int_t   firstTpcHitRow() const { return (Int_t)mFirstTpcHitRow; }
  //Int_t   lastTpcHitRow() const  { return (Int_t)mLastTpcHitRow; }  
  //Int_t   nHitsMapHFT() const    { return (Int_t)(mMap0 >> 1 & 0x7F); }
  //UInt_t  map0() const { return (UInt_t)mMap0; }
  //UInt_t  map1() const { return (UInt_t)mMap1; }
  //Float_t dEdx() const           { return (Float_t)mDedx/1000.; }
  //Float_t tofMatchFlag() const           { return mTofMatchFlag; }
  Float_t beta() const           { return (Float_t)mBeta/20000.; }
  Float_t localY() const          {return (Float_t)mLocalY/1000.;}
  //Float_t localZ() const          {return (Float_t)mLocalZ/1000.;}
  Float_t nSigmaElectron() const { return (Float_t)mNSigmaElectron/1000.; }
  Float_t dca() const           { return (Float_t)abs(mDca)/10000.; }
  Float_t dcaXY() const           { return (Float_t)mDcaXY/10000.; }
  Float_t dcaZ() const           { return (Float_t)mDcaZ/10000.; }
  Float_t e0() const             { return (Float_t)mBTOWE0/1000.;} 
  Float_t e() const             { return (Float_t)mBTOWE/1000.;} 
  Float_t pve() const; 
  UChar_t nEta() const          {return mBSMDNEta;}
  UChar_t nPhi() const          {return mBSMDNPhi;}
  Float_t zDist() const          {return mBEMCDistZ/1000.;}
  Float_t phiDist() const          {return mBEMCDistPhi/10000.;}
  Float_t etaTowDist() const          {return mBTOWDistEta/10000.;}
  Float_t phiTowDist() const          {return mBTOWDistPhi/10000.;}
  Int_t towerId() const          {return mBTOWId;}
  Short_t adc0() const          {return mBTOWADC0;}
  Int_t emcTriggerId() const          {return mEmcTrgId;}
  void	setEmcTriggerId(Int_t id) 	{mEmcTrgId=id;}

  //StDcaGeometry dcaGeometry() const;
  StPhysicalHelixD helix(float bField) const;
  Bool_t isHFTTrack() const { return mDca<0?true:false; }

 protected:
  Short_t mId;               // track index in picoDst
  StThreeVectorF mPMom;  // primary momentum, (0.,0.,0.) if none
  StThreeVectorF mGMom;
  //UShort_t mDedx;             // dEdx*1000
  Short_t  mDca;              // dca * 10000 * (isHFT?-1:1)
  Short_t  mDcaXY;            // dcaXY * 10000
  Short_t  mDcaZ;             // dcaZ * 10000
  Char_t   mNHitsFit;         // q*nHitsFit
  //Char_t   mNHitsMax;         // nHitsMax
  UChar_t  mNHitsDedx;        // nHitsDedx
  Short_t  mNSigmaElectron;   // nsigmaE * 1000
  //UInt_t   mMap0;             // TopologyMap data0 HFT + TPC
  //UInt_t   mMap1;             // TopologyMap data1 TPC + Others
  StThreeVectorF mOrigin;
  
  // pidTraits
  //Char_t   mTofMatchFlag;
  UShort_t mBeta;   // *20000
  Short_t  mLocalY; // *1000

  // these variables are extracted from the standard BEMC cluster algorithm
  Short_t  mBTOWADC0;         // adc0 higest adc in the cluster
  Short_t  mBTOWE0;           // E0*1000 highest tower in the cluster
  Short_t  mBTOWE;            // EMC point E*1000 
  Short_t  mBEMCDistZ;        // z*1000
  Short_t  mBEMCDistPhi;      // phi*10000
  Short_t  mBTOWDistEta;      // eta*10000 distance between track and matched tower center
  Short_t  mBTOWDistPhi;      // phi*10000 distance between track and matched tower center
  UChar_t  mBSMDNEta;         // # of hits in eta
  UChar_t  mBSMDNPhi;         // # of hits in phi

  // these variables are purely from single tower or nearby towers  
  Short_t  mBTOWId;           // projected tower Id 1-4800
  Short_t  mEmcTrgId;

  friend class StPicoDst;

  ClassDef(StElectronTrack, 1)
};

inline Float_t StElectronTrack::gPt() const
{
  //return 1./fabs(mPar[3]);
  return mGMom.perp();
}
inline Float_t StElectronTrack::gEta() const
{

  //float ptt = gPt();
  //StThreeVectorF gmom(ptt*std::cos(mPar[2]),ptt*std::sin(mPar[2]),ptt*mPar[4]);
  //return gmom.pseudoRapidity();
  return mGMom.pseudoRapidity();
}
inline Float_t StElectronTrack::gPhi() const
{
  //float ptt = gPt();
  //StThreeVectorF gmom(ptt*std::cos(mPar[2]),ptt*std::sin(mPar[2]),ptt*mPar[4]);
  //return gmom.phi();
  return mGMom.phi();;
}

inline StThreeVectorF StElectronTrack::gMom() const
{
  //float ptt = gPt();
  //return StThreeVectorF(ptt*std::cos(mPar[2]),ptt*std::sin(mPar[2]),ptt*mPar[4]);
  return mGMom;
}

//inline StDcaGeometry StElectronTrack::dcaGeometry() const
//{
//  StDcaGeometry a;
//  a.set(mPar, mErrMatrix);
//  return a;
//}
//      
//inline StPhysicalHelixD StElectronTrack::helix() const
//{
//  return dcaGeometry().helix();
//}        

inline StPhysicalHelixD StElectronTrack::helix(float bField) const
{
  	//return StPhysicalHelixD(mCurv,mDip,mPhase,mOrigin,mH);
  	return StPhysicalHelixD(mGMom,mOrigin,bField*kilogauss,charge());
}        

inline Float_t StElectronTrack::pve() const
{
    float p = pMom().mag();
    float eTow = e();
    if(eTow!=0) return p/eTow;
    else return -999.;
}        
#endif
