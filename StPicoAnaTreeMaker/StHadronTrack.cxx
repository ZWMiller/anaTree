#include "StHadronTrack.h"
#include "StMessMgr.h"
#include "TMath.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StPicoDstMaker/StPicoEmcPidTraits.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StBTofUtil/tofPathLength.hh"
#include "PhysicalConstants.h"
#include <climits>

ClassImp(StHadronTrack)

//----------------------------------------------------------------------------------
StHadronTrack::StHadronTrack() : mId(-1),
  //mGPt(0),mGEta(32768),mGPhi(32768),
   mBeta(0), mDca(32768), 
   mNHitsFit(0), mNHitsDedx(0), mNSigmaPion(32768), mNSigmaKaon(32768)
{

}

/////////////////////////////////////////////////////////////////////////////////////////
// t - the global track.  p - the associated primary track from the first primary vertex
/////////////////////////////////////////////////////////////////////////////////////////
//----------------------------------------------------------------------------------
StHadronTrack::StHadronTrack(StPicoDst *picoDst, StPicoTrack* t, Int_t idx)
   : mId(idx)
     //,mGPt(0),mGEta(32768),mGPhi(32768)
{
      mPMom      = t->pMom();
      int q      = t->charge();
      mNHitsFit  = t->nHitsFit()*q;
      mNHitsDedx = (UChar_t)(t->nHitsDedx());
      mNSigmaPion = (fabs(t->nSigmaPion() * 1000.) > 32768) ? 32768: (Short_t)(TMath::Nint(t->nSigmaPion() * 1000.));
      mNSigmaKaon = (fabs(t->nSigmaKaon() * 1000.) > 32768) ? 32768 : (Short_t)(TMath::Nint(t->nSigmaKaon() * 1000.));
      //mNSigmaElectron = (fabs(t->nSigmaElectron() * 100.) > 32768) ? 32768 : (Short_t)(TMath::Nint(t->nSigmaElectron() * 100.));

	   StThreeVectorF vertexPos = picoDst->event()->primaryVertex();
      StPhysicalHelixD helix = t->helix();
      double thePath = helix.pathLength(vertexPos);
      StThreeVectorF dcaPos = helix.at(thePath);
      mDca = fabs((dcaPos-vertexPos).mag()*10000.)>32768? 32768: (Short_t)((dcaPos-vertexPos).mag()*10000.);
      bool isHft = t->isHFTTrack();
      if(isHft) mDca *= -1;
      
    //  StThreeVectorF gMom = t->gMom(vertexPos,picoDst->event()->bField());
    //  mGPt = gMom.perp()*1000.>65535?65535:(UShort_t)(gMom.perp()*1000.);
    //  mGEta = fabs(gMom.pseudoRapidity()*10000)>32768? 32768:(Short_t)(gMom.pseudoRapidity()*10000.);
    //  mGPhi = fabs(gMom.phi()*10000)>32768? 32768:(Short_t)(gMom.phi()*10000.);

      //StThreeVectorF dcaPoint = helix.at(helix.pathLength(vertexPos.x(), vertexPos.y()));
      //float dcaZ = (dcaPoint.z() - vertexPos.z())*10000.;
      //float dcaXY = (helix.geometricSignedDistance(vertexPos.x(),vertexPos.y()))*10000.;
      //mDcaZ = dcaZ>32768?32768:(Short_t)dcaZ;
      //mDcaXY = dcaXY>32768?32768:(Short_t)dcaXY;

      int index2TofPid = t->bTofPidTraitsIndex();
      Float_t localY = -999.;
      if (index2TofPid>=0){
        StPicoBTofPidTraits *tofPid = picoDst->btofPidTraits(index2TofPid);
        //mTofMatchFlag = tofPid->btofMatchFlag();
        //Float_t mom = mGMom.mag();
        localY = tofPid->btofYLocal();
        mLocalY = tofPid->btofYLocal()*1000;
        //mLocalZ = tofPid->btofZLocal();
        Float_t beta = tofPid->btofBeta();
        if(beta<1e-4||beta>=(USHRT_MAX-1)/20000){
           Float_t tof = tofPid->btof();
           StThreeVectorF btofHitPos = tofPid->btofHitPos();
           float L = tofPathLength(&vertexPos, &btofHitPos, helix.curvature()); 
           beta = L/(tof*(c_light/1.0e9));
        }
        mBeta = (UShort_t)(beta*20000);
        
      }

   int index2EmcPid = t->emcPidTraitsIndex();
   if (index2EmcPid>=0){
      StPicoEmcPidTraits *emcPid = picoDst->emcPidTraits(index2EmcPid);
      mBTOWADC0 = emcPid->adc0();
      mBTOWId = emcPid->btowId();
      mEmcTrgId = -1;
   }else{
      mBTOWADC0 = 0;
      mBTOWId = 0;
      mEmcTrgId = -1;
   }
}

//----------------------------------------------------------------------------------
StHadronTrack::~StHadronTrack()
{
   /* noop */
}
//----------------------------------------------------------------------------------
void StHadronTrack::Print(const Char_t *option) const
{
      LOG_INFO << "id=" << id() << endm;
      LOG_INFO << "pMom=" << pMom() << endm;
      LOG_INFO << " nHitsFit = " << nHitsFit() << endm;
}


