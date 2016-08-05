#ifndef StPicoQAMaker_h
#define StPicoQAMaker_h

#include "StMaker.h"
class StPicoDst;
class StPicoEvent;
class StPicoDstMaker;
class StPicoMtdHit;
class StPicoTrack;
class TString;
class TH1F;
class TH2F;
class TH3F;
class TFile;

class StPicoQAMaker : public StMaker {
   public:
      StPicoQAMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName);
      Int_t Init();
      Int_t Make();
      void  Clear(Option_t *opt="");
      Int_t Finish();
      void  setRunList(TString name) {mRunFileName = name;};
      void  setNumberOfRuns(int rns) {mTotalRun = rns;};
   private:
      Bool_t Isgoodtrack(StPicoTrack* ); 
      Bool_t Ismuontrack(StPicoTrack* );
      Bool_t Ispasseventcuts(StPicoEvent*);
      Bool_t isBHT1(StPicoEvent*);
      Bool_t isBHT2(StPicoEvent*);
      Bool_t isBHT3(StPicoEvent*);
      Double_t RotatePhi(Double_t phi) const;   
      void    DeclareHistograms();

      StPicoDstMaker *mPicoDstMaker;
      StPicoDst      *mPicoDst;
      TString    mOutName;
      Int_t   mNBadRuns;       
      TString mRunFileName;
      Int_t mTotalRun;

      map<Int_t,Int_t> mTotalRunId;
      TFile*	   fout;
      //-----event QA-----
      TH1F*      hNEvents;
      //TH1F*    	htriggerindex;
      //TH1F*      mVz_vpd;
      //TH1F*      mVz_tpc;
      TH2F*      mVz_vpdtpc;
      //TH1F*      mdVz;
      TH2F*      mdVz_tpcVz;
      TH2F*      mVxy;
      TH2F*      mVRvsVZ;

      
      TH2F*      mRanking_nPtrks;
      TH2F*      mnPtrks_nGtrks;
      TH2F*      mnRefMult_nGRefMult;
      TH2F*      mnRefMult_nPtrks;
      TH2F*      mnRefMult_nGtrks;
      TH2F*      mnGRefMult_nPtrks;
      TH2F*      mnGRefMult_nGtrks;

      TH2F*      mnPtrks_nTofHits;
      TH2F*      mnPtrks_nMtdHits;

      TH2F*      mnTofHits_nMtdHits;
      TH2F*      mnTofMatch_nTofHits;


      //---------Track QA----------
      //-------TPC information------
      TH1F*      mNptracks;
      TH1F*      mNgtracks;
      TH1F*      mtrkpt;
      TH1F*      mtrkphi;
      TH1F*      mtrketa;
      TH1F*      mtrkdca;
      TH1F*      mtrkdcawHFT;
      TH1F*      mtrkdcaXY;
      TH1F*      mtrkdcaXYwHFT;
      TH1F*      mtrkdcaZ;
      TH1F*      mtrkdcaZwHFT;
      TH1F*      mnhitsdedx;
      TH1F*      mnhitsfit;
      TH1F*      mnhitsfitRatio;
      TH1F*      mnhitsdedxwHFT;
      TH1F*      mnhitsfitwHFT;
      TH1F*      mnhitsfitRatiowHFT;

      TH1F*      mnsigmaPI;
      TH1F*      mnsigmaK;
      TH1F*      mnsigmaE;
      TH1F*      mnsigmaP;

      TH2F*      mtrketaphi;
      TH2F*      mtrketa_pt;
      TH2F*      mtrkphi_pt;
      TH2F*      mtrketa_phiPos;
      TH2F*      mtrketa_phiNeg;
      TH2F*      mtrkdca_pt;
      TH2F*      mtrkdcaXY_pt;
      TH2F*      mtrkdcaZ_pt;
      TH2F*      mtrkdcawHFT_pt;
      TH2F*      mtrkdcaXYwHFT_pt;
      TH2F*      mtrkdcaZwHFT_pt;
      TH2F*      mnhitsfit_pt;
      TH2F*      mnhitsdedx_pt;
      TH2F*      mnhitsRatio_pt;
      TH2F*      mnhitsRatio_ptwHFT;
      TH2F*      mdedx_P;

      TH2F*      mnsigmaPI_P;
      TH2F*      mnsigmaP_P;
      TH2F*      mnsigmaK_P;
      TH2F*      mnsigmaE_P;

      //----TPC information end-----

      //-----TOF INFORMATION start----
      TH2F*      minvsBeta_P;
      TH2F*      mtofM2_P;
      //        mp vs nsigmaE
      TH2F*      mtoftray_localY;
      TH2F*      mtoftray_localZ;
      //TH3F*      mtofhitPosXYZ;
      TH2F*      mtoftray_matchflag;
      TH2F*      mtoftray_module;

      //----TOF information end----

      //---------mtd information-------------
      TH2F*      mmtdbgcell_module;//backleg*12+cell vs module
      TH2F*      mmtddeltaT_pt;
      TH2F*      mmtddeltaZ_pt;
      TH2F*      mmtddeltaY_pt;
      TH2F*      mmtdBeta_P;
      TH2F*      mmtdmatchflag_channel;
      TH2F*      mmtddeltaT_channel;
      TH2F*      mmtddeltaY_channel;
      TH2F*      mmtddeltaZ_channel;

      TH2F*      mmtddeltaT_channel1;
      TH2F*      mmtddeltaY_channel1;
      TH2F*      mmtddeltaZ_channel1;

      TH2F*      mmtdBeta_channel;

      //-----------BEMC start---------
      TH2F*      mbTowId_P;
      TH2F*      mbTowId2_P;
      TH2F*      mbTowId3_P;
      //TH2F*      mbTowE1_P;
      TH2F*      mBEMCe0_P;
      TH2F*      mBEMCe1_P;
      TH2F*      mBEMCe2_P;
      TH2F*      mBEMCe3_P;
      TH2F*      mBEMCE_P;
      TH2F*      mBEMCadc0_P;
      TH2F*      mBEMCzdist_P;
      TH2F*      mBEMCphidist_P;
      TH2F*      mBEMCneta_P;
      TH2F*      mBEMCnphi_P;
      TH2F*      mBEMCneta_nphi;
      

      //======run by run QA=====

      TH2F *hTPCVzvsRunIndex;
      TH2F *hTPCVxvsRunIndex;
      TH2F *hTPCVyvsRunIndex;
      TH2F *hVPDVzvsRunIndex;
      TH2F *hDeltaZvsRunIndex;

      TH2F *hRefMultvsRunIndex;
      TH2F *hGRefMultvsRunIndex;
      TH2F *hRefMultCorrvsRunIndex;  
      TH2F *hnPrimaryvsRunIndex; 
      TH2F *hnGlobalvsRunIndex;

      TH2F *hBBCXvsRunIndex;   
      TH2F *hZDCXvsRunIndex;

      TH2F *hnhitsfitvsRunIndex; 
      TH2F *hDedxvsRunIndex;
      TH2F *hnhitsdedxvsRunIndex;
      TH2F *hPtvsRunIndex;
      TH2F *hEtavsRunIndex;
      TH2F *hPhivsRunIndex;

      TH2F *hNSigmaEvsRunIndex;
      TH2F *hNSigmaPivsRunIndex;
      TH2F *hNSigmaKvsRunIndex;
      TH2F *hNSigmaPvsRunIndex;

      TH2F *hDcavsRunIndex;
      TH2F *hDcawHFTvsRunIndex;
      TH2F *hDcaXYvsRunIndex;
      TH2F *hDcaXYwHFTvsRunIndex;
      TH2F *hDcaZvsRunIndex;
      TH2F *hDcaZwHFTvsRunIndex;

      //TH2F *hHTth0vsRunIndex;   
      //TH2F *hHTth1vsRunIndex;   
      //TH2F *hHTth2vsRunIndex;
      //TH2F *hHTth3vsRunIndex;   


      TH2F *hntofmatchvsRunIndex;
      TH2F *hnbemchitsvsRunIndex;
      TH2F *hnbemcmatchvsRunIndex;
      TH2F *hnmtdhitsvsRunIndex;
      TH2F *hnmtdmatchvsRunIndex;
      TH2F *hnhftmatchvsRunIndex;  
      TH2F *hNtofElectronvsRunIndex;
      TH2F *hNbemcElectronvsRunIndex;
      TH2F *hNmuonvsRunIndex;

      TH2F *hnmtdhitsOVnGRvsRunIndex;
      TH2F *hntofmatchOVnGRvsRunIndex;
      TH2F *hnbemcmatchOVnGRvsRunIndex;
      TH2F *hnmtdmatchOVnGRvsRunIndex;
      TH2F *hnhftmatchOVnGRvsRunIndex;  
      TH2F *hNtofElectronOVnGRvsRunIndex;
      TH2F *hNbemcElectronOVnGRvsRunIndex;
      TH2F *hNmuonOVnGRvsRunIndex;

      TH2F *hnmtdhitsOVnRvsRunIndex;
      TH2F *hntofmatchOVnRvsRunIndex;
      TH2F *hnbemcmatchOVnRvsRunIndex;
      TH2F *hnmtdmatchOVnRvsRunIndex;
      TH2F *hnhftmatchOVnRvsRunIndex;  
      TH2F *hNtofElectronOVnRvsRunIndex;
      TH2F *hNbemcElectronOVnRvsRunIndex;
      TH2F *hNmuonOVnRvsRunIndex;

      //-----------BEMC end-------------------

      //----------after pid------------
      /*
      TH2F *htofPiondcavsPT;
      TH2F *htofKaondcavsPT;
      TH2F *htofProtondcavsPT;
      TH2F *htofElectrondcavsPT;
      TH2F *hbemcElectrondcavsPT;
      TH2F *hMuondcavsPT;
      //after pid -----wHFT------
      TH2F *htofPiondcawHFTvsPT;
      TH2F *htofKaondcawHFTvsPT;
      TH2F *htofProtondcawHFTvsPT;
      TH2F *htofElectrondcawHFTvsPT;
      TH2F *hbemcElectrondcawHFTvsPT;
      TH2F *hMuondcawHFTvsPT;

      //-----only cuts on tpc |nsigma(pi/k/p)|--------------
      TH2F *htpcPiondcavsPT;
      TH2F *htpcKaondcavsPT;
      TH2F *htpcProtondcavsPT;
      //-----tpc cuts + tof cuts----------
      TH2F *htpctofPiondcavsPT;
      TH2F *htpctofKaondcavsPT;
      TH2F *htpctofProtondcavsPT;

      //-----only cuts on tpc |nsigma(pi/k/p)|-----wHFT---------
      TH2F *htpcPiondcawHFTvsPT;
      TH2F *htpcKaondcawHFTvsPT;
      TH2F *htpcProtondcawHFTvsPT;
      //-----tpc cuts + tof cuts---wHFT-------
      TH2F *htpctofPiondcawHFTvsPT;
      TH2F *htpctofKaondcawHFTvsPT;
      TH2F *htpctofProtondcawHFTvsPT;
      */
      ClassDef(StPicoQAMaker,1)
};

#endif
