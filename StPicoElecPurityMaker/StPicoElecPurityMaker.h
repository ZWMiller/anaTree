#ifndef StPicoElecPurityMaker_h
#define StPicoElecPurityMaker_h
#include "StRoot/StPicoDstMaker/StPicoTrack.h"
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
class THnSparse;
class TFile;
class StPicoElecPurityMaker : public StMaker {
 public:
     StPicoElecPurityMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName);
     virtual ~StPicoElecPurityMaker();
     Bool_t passGoodTrack(StPicoEvent*, StPicoTrack*, int ); // ZWM
     Bool_t passGoodTrack_NoEta(StPicoEvent*, StPicoTrack*, int ); // ZWM
     Bool_t passBEMCCuts(StPicoEvent*, StPicoTrack*, int );   // ZWM
     Bool_t passTOFCuts(StPicoEvent*, StPicoTrack*, int );   // ZWM
     Int_t passSMDCuts(StPicoEvent*, StPicoTrack*, int );   // ZWM, return INT becuase could fit multiple cases
     Bool_t Ismuontrack(StPicoEvent*, StPicoTrack* );
     Bool_t IspassTOFcuts(StPicoTrack*);
     Bool_t IspassBEMCcuts(StPicoTrack*);
     Bool_t passEventCuts(StPicoEvent*,int);  // ZWM
     Bool_t passEventCuts_NodVz(StPicoEvent*,int);  // ZWM
     Bool_t checkHotTower(int, int); //ZWM
     Bool_t checkTriggers(StPicoEvent*, int);
     Bool_t isMB(StPicoEvent*);
     Bool_t isBHT0(StPicoEvent*);
     Bool_t isBHT1(StPicoEvent*);
     Bool_t isBHT2(StPicoEvent*);
     Bool_t isBHT3(StPicoEvent*);
     Double_t RotatePhi(Double_t phi) const;   
     void setRunList(TString name) { mRunFileName = name;};
     //virtual Int_t Init();
     //virtual Int_t Make();
     //virtual void  Clear(Option_t *opt="");
     //virtual Int_t Finish();
     Int_t Init();
     Int_t Make();
     void  Clear(Option_t *opt="");
     Int_t Finish();

     void    DeclareHistograms();
     void    WriteHistograms();
     Int_t    FillHistograms(Int_t, StPicoEvent*);
     Int_t    dVzStudy(StPicoEvent*);

     void  SetDefaultCuts();               // ZWM

      // ZWM Functions
     void   setNSigECuts(float l, float h)  { nSigELow = l;  nSigEHigh = h; };
     void   setNSigPCuts(float l, float h)  { nSigPLow = l;  nSigPHigh = h; };
     void   setNSigKCuts(float l, float h)  { nSigKLow = l;  nSigKHigh = h; };
     void   setNSigPiCuts(float l, float h) { nSigPiLow = l; nSigEHigh = h; };
     void   setvZCuts(int tr, float vz, float dvz)  { vZcut[tr] = vz; dvZcut[tr] = dvz; };
     void   setvZCutsHFT(int tr, float vz, float dvz)  { vZcutHFT[tr] = vz; dvZcutHFT[tr] = dvz; };
     void   setPrimaryPtCut(float tpt, float bpt ) { bemcPtCut = bpt; tofPtCut = tpt; };
     void   setPrimaryEtaCut(float et)      { etaCut = et; };
     void   setPrimaryDCACut(float dca)     { dcaCut = dca; };
     void   setNhitsCuts(float dedx, float fit, float ratio) 
        { nhitsdEdxCut = dedx; nhitsFitCut = fit; nhitsRatioCut = ratio; };
     void   setPoECut(float pEl, float pEh) { poeCutLow = pEl; poeCutHigh = pEh; };
     void   setToFBetaCut(float iB)         { tofInvBetaCut = iB; };
     void   setToFLocalyCut(float lY)       { toflocalyCut = lY; };
     void   setKaonEnhCut(float kel, float keh) { kaonEnhCutLow = kel; kaonEnhCutHigh = keh; };
     void   setPionEnhCut(float pel, float peh) { pionEnhCutLow = pel; pionEnhCutHigh = peh; };
     void   setProtonEnhCut(float pel, float peh) { protonEnhCutLow = pel; protonEnhCutHigh = peh; };
     void   setDsmAdcCut(int trg, int val) { dsmAdcCut[trg] = val; };
     int    getDsmAdcCut(int trg)            { return dsmAdcCut[trg]; };
     void   setSMDCuts(int ne, int np, float zd, float pd) 
     {nEtaCut = ne; nPhiCut = np; zDistCut = zd; phiDistCut = pd; };
     void   setSMDCuts2(int ne, int np, float zd, float pd) 
     {nEtaCut2 = ne; nPhiCut2 = np; zDistCut2 = zd; phiDistCut2 = pd; };

     void addTrigger(int tr, int id) { triggers[id].push_back(tr); };
     void setRunMode(int md) { runMode = md; };

  private:
   StPicoDstMaker *mPicoDstMaker;
   StPicoDst      *mPicoDst;
   
   TString mRunFileName;
   map<Int_t,Int_t> mTotalRunId;
   const int pp = 0;
   const int AuAu = 1;
   int runMode;
   vector<int> triggers[5]; //0-HT0, 1-HT1 ... 4-MB
   // dsm adc
   int dsmAdcCut[4];
   // Trigger Tracking
   int numTrigs;
   int trig;
   // Event cut vars
   float vZcut[4];
   float dvZcut[4];
   float vZcutHFT[4];
   float dvZcutHFT[4];
   // Track cut vars
   float bemcPtCut, tofPtCut, etaCut, dcaCut;
   float nhitsdEdxCut, nhitsFitCut, nhitsRatioCut;
   float poeCutLow, poeCutHigh;
   float tofInvBetaCut,toflocalyCut;
   
   // SMD cuts (1 = Xiaozhi run 12 cuts, 2 = Daniel run 10 cuts [tighter])
   int nEtaCut, nPhiCut;
   float zDistCut, phiDistCut;
   int nEtaCut2, nPhiCut2;
   float zDistCut2, phiDistCut2;

   //TPC cuts (not used for purity, here in case)
   float nSigELow, nSigEHigh;
   float nSigPLow, nSigPHigh;
   float nSigKLow, nSigKHigh;
   float nSigPiLow, nSigPiHigh;
   float pionEnhCutLow,pionEnhCutHigh;
   float protonEnhCutLow,protonEnhCutHigh;
   float kaonEnhCutLow,kaonEnhCutHigh;
   int trigCounter;
   TString    mOutName;

     Int_t   mNBadRuns;       
     Int_t   mNHotTower1;
     Int_t   mNHotTower2;
     Int_t   mNHotTower3;
     Int_t  trkHFTflag;
    
   TFile*	   fout;
    //-----event QA-----

    TH1F*      trigType;
    TH1F*      hNEvents[4];
    TH1F*    	htriggerindex[4];
    TH1F*      mVz_vpd[4];
    TH1F*      mVz_tpc[4];
    TH2F*      mVz_vpdtpc[4];
    TH1F*      mdVz[4];
    TH2F*      mdVz_tpcVz[4];
    TH2F*      mVxy[4];
    TH2F*      mVRvsVZ[4];

   // TH2F*      mRanking_nPtrks[4];
  //  TH2F*      mnPtrks_nGtrks[4];
  //  TH2F*      mnRefMult_nGRefMult[4];
  //  TH2F*      mnRefMult_nPtrks[4];
  //  TH2F*      mnRefMult_nGtrks[4];
  //  TH2F*      mnGRefMult_nPtrks[4];
  //  TH2F*      mnGRefMult_nGtrks[4];

    //TH2F*      mnPtrks_nTofHits[4];
    //TH2F*      mnPtrks_nMtdHits[4];

   // TH2F*      mnTofHits_nMtdHits[4];
   // TH2F*      mnTofMatch_nTofHits[4];
    

    //---------Track QA----------
    //-------TPC information------
    TH1F*      mNptracks[4];
    TH1F*      mNgtracks[4];
    TH1F*      mtrkpt[4];
    TH1F*      mtrkphi[4];
    TH1F*      mtrketa[4];
    TH1F*      mnsigmaK[4];
    TH1F*      mnsigmaE[4];
    TH1F*      mnsigmaP[4];
    TH1F*      mnsigmaPI[4];

    TH2F*      mtrketaphi[4];
    TH2F*      mtrketa_pt[4];
    TH2F*      mtrkphi_pt[4];
    TH2F*      mdedx_Pt[4];

    TH2F*      mnsigmaPI_Pt[4][2]; // 0 = no HFT requirement, 1 = w/HFT Only, Only do this for hists that are used in fits 
    TH2F*      mnsigmaP_Pt[4][2];
    TH2F*      mnsigmaK_Pt[4][2];
    TH2F*      mnsigmaE_Pt[4][2]; 

  /*  TH2F*      mnSigmaE_Pt_SMD[4][2];
    TH2F*      mnSigmaP_Pt_SMD[4][2];
    TH2F*      mnSigmaK_Pt_SMD[4][2];
    TH2F*      mnSigmaPI_Pt_SMD[4][2];

    TH2F*      mnSigmaE_Pt_SMD2[4][2];
    TH2F*      mnSigmaP_Pt_SMD2[4][2];
    TH2F*      mnSigmaK_Pt_SMD2[4][2];
    TH2F*      mnSigmaPI_Pt_SMD2[4][2];

    TH2F*      mnSigmaE_Pt_TOF[4][2];
    TH2F*      mnSigmaP_Pt_TOF[4][2];
    TH2F*      mnSigmaK_Pt_TOF[4][2];
    TH2F*      mnSigmaPI_Pt_TOF[4][2];
    
    TH2F*      mnSigmaE_Pt_BEMC[4][2];
    TH2F*      mnSigmaP_Pt_BEMC[4][2];
    TH2F*      mnSigmaK_Pt_BEMC[4][2];
    TH2F*      mnSigmaPI_Pt_BEMC[4][2];
    */
    
    //----TPC information end-----

    //-----TOF INFORMATION start----
    /*TH2F*      minvsBeta_Pt[4];
    TH2F*      mtofM2_Pt[4];
    //        mp vs nsigmaE
    TH2F*      mtoftray_localY[4];
    TH2F*      mtoftray_localZ[4];
    //TH3F*      mtofhitPosXYZ[4];
    TH2F*      mtoftray_matchflag[4];
    TH2F*      mtoftray_module[4];
*/
    //----TOF information end----
    
    // For Purity
  /*  TH2F*      mnSigmaEvsBeta[4];
    TH2F*      mnSigmaPIvsBeta[4];
    TH2F*      mnSigmaKvsBeta[4];
    TH2F*      mnSigmaPvsBeta[4];
    TH2F*      mdedxvsBeta[4];
    TH2F*      mtofm2vsBeta[4];
	*/
    TH1F*      hNTracks[4];
    /*TH2F*      mnSigmaE_KEnh_Pt[4][2];
    TH2F*      mnSigmaE_PiEnh_Pt[4][2];
    TH2F*      mnSigmaE_PEnh_Pt[4][2];
    TH2F*      mnSigmaK_KEnh_Pt[4][2];
    TH2F*      mnSigmaPi_PiEnh_Pt[4][2];
    TH2F*      mnSigmaP_PEnh_Pt[4][2];
*/
    // For Eta Dependence Study
    THnSparse*      mnSigmaE_TOF[4][2];
    THnSparse*      mnSigmaE_SMD[4][2];
    THnSparse*      mnSigmaE_SMD2[4][2];
    THnSparse*      mnSigmaE_BEMC[4][2];

    // For Centrality Study
    TH1F* gRefMult[4];
    TH1F* gRefMultCor[4];
    TH1F* gRefMultCorWg[4];
    TH1F* centrality16[4];


    // For dVz cut study
    TH2F* mTPCvsVPD_Vz;
    TH2F* mTPCvsDVz;

    ClassDef(StPicoElecPurityMaker, 1)
};

#endif
