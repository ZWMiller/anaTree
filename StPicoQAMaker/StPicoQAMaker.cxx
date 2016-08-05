#include "StPicoQAMaker.h"
#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoTrack.h"
//#include "StRoot/StPicoDstMaker/StPicoV0.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoDstMaker/StPicoEvent.h"
#include "StRoot/StPicoDstMaker/StPicoMtdHit.h"
#include "StRoot/StPicoDstMaker/StPicoConstants.h"
#include "StRoot/StPicoDstMaker/StPicoMtdPidTraits.h"
#include "StRoot/StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StRoot/StPicoDstMaker/StPicoEmcPidTraits.h"
#include "StDcaGeometry.h"

#include "StThreeVectorF.hh"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include <map>
#include <iostream>
#include <fstream>
#include "mBadRunList.h"

Bool_t fillhistflag=1;
ofstream runidfiles;
Int_t runIndex;
Int_t randomId;


ClassImp(StPicoQAMaker)
  //-----------------------------------------------------------------------------
  StPicoQAMaker::StPicoQAMaker(const char* name, StPicoDstMaker *picoMaker, const char* outName)
: StMaker(name)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = NULL;
  TH1F:: SetDefaultSumw2();//zaochen add
  mOutName = outName;
  mTotalRun = 1000;

  mNBadRuns = sizeof(mBadRuns)/sizeof(int);
  mRunFileName = "StRoot/StPicoQAMaker/mTotalRunList.dat";

}

//----------------------------------------------------------------------------- 
Int_t StPicoQAMaker::Init() {
  if(mOutName!="") {
    fout = new TFile(mOutName.Data(),"RECREATE");
  }else{
    fout = new TFile("picoQA_test.root","RECREATE");
  }
  LOG_INFO<<"StPicoQAMaker::Init  mOutName = "<<mOutName.Data()<<endm;
  DeclareHistograms();

  //runidfiles.open("/star/u/zye20/zye20/zaochen/analysis/run14picoQA/runnumber.txt", ios::out|ios::app);


  if(fillhistflag){
    //read in the runlist.dat
    ifstream indata;
    indata.open(mRunFileName);
    mTotalRunId.clear();
    if(indata.is_open()){
      cout<<"read in total run number list and recode run number ...";
      Int_t oldId;
      Int_t newId=0;
      while(indata>>oldId){
        mTotalRunId[oldId] = newId;
        newId++;
      }
      cout<<" [OK]"<<endl;  

    }else{
      cout<<"Failed to load the total run number list !!!"<<endl;
      return kFALSE;
    }

    indata.close();

    //for(map<Int_t,Int_t>::iterator iter=mTotalRunId.begin();iter!=mTotalRunId.end();iter++)
    //cout<<iter->second<<" \t"<<iter->first<<endl;
    //cout<<endl;
    // 
  }//

  return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StPicoQAMaker::Finish() {
  fout->cd();
  fout->Write();
  fout->Close();
  LOG_INFO<<"Histogram saved"<<endm;
  return kStOK;
}

//-----------------------------------------------------------------------------
void StPicoQAMaker::DeclareHistograms() {

  fout->cd();
  hNEvents   = new TH1F("hNEvents","number of events: 0 for total and 2 for MBs events", 10, 0, 10 );
  //htriggerindex =new TH1F("htriggerindex","triggerindex", 64,0,64);
  //mVz_vpd     = new TH1F("vz_vpd", "VZ_VPD distribution (cm)",400,-200,200);
  //mVz_tpc     = new TH1F("vz_tpc","the Vz_TPC distribution (cm) ",400,-200,200);
  mVz_vpdtpc  = new TH2F("Vz_vpdtpc","VertexZ: VPD VS TPC;TPC;VPD;", 400,-200,200,400,-200,200);
  //mdVz        = new TH1F("dVz", "VertexZ 'TPC-VPD' (cm)", 400,-200,200);
  mdVz_tpcVz  = new TH2F("dVz_tpcVz","VertexZ 'TPC-VPD vs TPC' (cm);VZtpc; VZ(tpc-vpd)",400,-200,200,200,-100,100);
  mVxy        = new TH2F("vxy","Vertex: Vy Vs Vx distribution (cm);vx;vy; ",200,-10,10,200,-10,10);
  mVRvsVZ    = new TH2F("VRvsVZ","Vertex: VR vs VZ  (cm);VZ;VR; ", 400, -200, 200,200,0,20);


  mRanking_nPtrks=new TH2F("Ranking_nPtrks","Ranking Vs nPrimaryTracks;nprimarytrks; ranking;",200,0,2000,100,-10,10);
  mnPtrks_nGtrks=new TH2F("nPtrks_nGtrks","nPrimarytrks_nGlobaltrks;nglobaltrks;nprimarytrks;",500,0,5000,200,0,2000);

  mnRefMult_nGRefMult=new TH2F("nRefMult_nGRefMult","nRefMult_nGRefMult;nGRefMult;nRefmult;",200,0,1000,200,0,1000);
  mnRefMult_nPtrks=new TH2F("nRefMult_nPtrks","nRefMult_nPrimarytrks;nprimarytrks;nRefmult;",200,0,2000,200,0,1000);
  mnRefMult_nGtrks=new TH2F("nRefMult_nGtrks","nRefMult_nGlobaltrks;nglobaltrks;nRefmult;",250,0,5000,200,0,1000);

  mnGRefMult_nPtrks=new TH2F("nGRefMult_nPtrks","nGRefMult_nPrimarytrks;nprimarytrks;nGRefmult;",200,0,2000,200,0,1000);
  mnGRefMult_nGtrks=new TH2F("nGRefMult_nGtrks","nGRefMult_nGlobaltrks;nglobaltrks;nGRefmult;",250,0,5000,200,0,1000);

  mnPtrks_nTofHits=new TH2F("nPtrks_nTofHits","nPrimarytrks_nTofHits;ntofhits;nprimarytrks;",250,0,2000,200,0,2000);
  mnPtrks_nMtdHits=new TH2F("nPtrks_nMtdHits","nPrimarytrks_nMtdHits;nmtdhits;nprimarytrks;",50,0,50,500,0,2000);
  mnTofHits_nMtdHits=new TH2F("nTofHits_nMtdHits","nTofHits_nMtdHits;nmtdhits;ntofhits;",50,0,50,200,0,2000);
  mnTofMatch_nTofHits=new TH2F("nTofMatch_nTofHits","nTofMatch_nTofHits;ntofhits;ntofmatch;",300,0,1500,300,0,1500);


  mNptracks    = new TH1F("Nptracks","#primary tracks",200,0,2000);    
  mNgtracks    = new TH1F("Ngtracks","# global tracks ",200,0,5000);    


  mnhitsfit_pt   = new TH2F("nhitsfit_pt","nhitsfit_pt; pt; nhitsfit;",200,0,20,50,0,50);
  mnhitsRatio_pt = new TH2F("nhitsRatio_pt","nhitsRatio_pt; pt; nhitsRatio;",200,0,20,100,0,1.2);
  mnhitsRatio_ptwHFT = new TH2F("nhitsRatio_ptwHFT","nhitsRatio_ptwHFT; pt; nhitsRatio;",200,0,20,100,0,1.2);      
  mnhitsdedx_pt  = new TH2F("nhitsdedx_pt","nhitsdedx_pt; pt; nhitsdedx;",200,0,20,50,0,50);

  mnhitsfit   = new TH1F("nhitsfit", "nhitsfit ",50,0,50);
  mnhitsfitRatio   = new TH1F("nhitsfitRatio", "nhitsfitRatio ",100,0,1.0);
  mnhitsfitwHFT   = new TH1F("nhitsfitwHFT", "nhitsfitwHFT",50,0,50);
  mnhitsfitRatiowHFT   = new TH1F("nhitsfitRatiowHFT", "nhitsfitRatiowHFT ",100,0,1.0);
  mnhitsdedx  = new TH1F("nhitsdedx", "nhitsfitdedx",50,0,50);
  mnhitsdedxwHFT  = new TH1F("nhitsdedxwHFT", "nhitsfitdedxwHFT",50,0,50);

  mtrketaphi  = new TH2F("trketaphi","the eta vs phi; phi; eta;",200,-0.5,6.5,200,-2,2);     
  mtrketa_pt  = new TH2F("trketa_pt", "trketa_pt; pt; eta;", 200,0,20,200,-2,2);
  mtrkphi_pt  = new TH2F("trkphi_pt","trkphi_pt; pt; phi;",200,0,20,200,0,6.3);
  mtrketa_phiPos = new TH2F("trketa_phiPos","trketa_phi(Positive); phi; eta;",200,0,6.3,200,-2,2);
  mtrketa_phiNeg = new TH2F("trketa_phiNeg","trketa_phi(Negative); Phi; eta;",200,0,6.3,200,-2,2);

  mtrkpt      = new TH1F("trkpt","the pt distribution of all tracks",200,0,20);
  mtrketa     = new TH1F("trketa","eta ",200,-2,2);
  mtrkphi     = new TH1F("trkphi","the phi distribution of all tracks",200,0,6.3);

  mtrkdca_pt     = new TH2F("trkdca_pt","trkdca_pt",200,0,20,250,0,10);
  mtrkdcaXY_pt     = new TH2F("trkdcaXY_pt","the dcaXY_PT; pt; dcaXY;",200,0,20,250,-10,10);
  mtrkdcaZ_pt      = new TH2F("trkdcaZ_pt","the dcaZ_PT; pt; dcaZ;",200,0,20,250,-10,10);

  mtrkdcawHFT_pt     = new TH2F("trkdcawHFT_pt","trkdcawHFT_pt",200,0,20,250,0,10);
  mtrkdcaXYwHFT_pt     = new TH2F("trkdcaXYwHFT_pt","dcaXYwHFT_PT; pt; dcaXY;",200,0,20,250,-10,10);
  mtrkdcaZwHFT_pt      = new TH2F("trkdcaZwHFT_pt","dcaZwHFT_PT; pt; dcaZ;",200,0,20,250,-10,10);

  mtrkdca     = new TH1F("trkdca","the dca of all tracks",500,0,10);
  mtrkdcawHFT       = new TH1F("trkdcawHFT","the dca With HFT",500,0,10);     
  mtrkdcaXY     = new TH1F("trkdcaXY","the dcaXY of all tracks",500,-10,10);
  mtrkdcaXYwHFT     = new TH1F("trkdcaXYwHFT","the dcaXYwHFT of all tracks",500,-10,10);
  mtrkdcaZ     = new TH1F("trkdcaZ","the dcaZ of all tracks",500,-10,10);
  mtrkdcaZwHFT      = new TH1F("trkdcaZwHFT","the dcaZwHFT of all tracks",500,-10,10);


  mnsigmaPI_P   = new TH2F("nsigmaPI_P", "nsigmapion vs P of all tracks; P; nsigmaPI;",200,0,20,100,-20,20);
  mnsigmaP_P   = new TH2F("nsigmaP_P", "nsigmaproton vs P of all tracks; P; nsigmaP;",200,0,20,100,-20,20);
  mnsigmaK_P   = new TH2F("nsigmaK_P", "nsigmaK vs P of all tracks; P; nsigmaK;",200,0,20,100,-20,20);
  mnsigmaE_P   = new TH2F("nsigmaE_P", "nsigmaE vs P of all tracks; P; nsigmaE;",200,0,20,100,-20,20);

  mnsigmaPI   = new TH1F("nsigmaPI", "nsigmapion of all tracks",50,-15,15);
  mnsigmaK    = new TH1F("nsigmaK", "nsigmaKaon of all tracks",50,-15,15);
  mnsigmaE    = new TH1F("nsigmaE", "nsigmaElectron of all tracks",50,-15,15);
  mnsigmaP    = new TH1F("nsigmaP", "nsigmaProton of all tracks",50,-15.,15.);

  mdedx_P        = new TH2F("Dedx_P","dedx(keV/cm) vs P; P; Dedx;",200,-20,20,250,0,10);
  minvsBeta_P  = new TH2F("invsBeta_P","1/Beta VS momentum; P; 1/beta;",300,0,6,200,0.90,1.12);
  mtofM2_P     = new TH2F("tofM2_P", "tofM2 VS momentum; P; tofM2;",200,0,20,200,-0.5,4.5);


  mtoftray_localY = new TH2F("toftray_localY","localY VS toftray; tray; localY;",120,0,120,100,-4,4);
  mtoftray_localZ = new TH2F("toftray_localZ","localZ VS toftray; tray; localZ;",120,0,120,100,-4,4);
  //mtofhitPosXYZ   = new TH3F("tofhitPosXYZ","tofhitPosXYZ",400,-200,200,400,-200,200,400,-200,200);
  mtoftray_matchflag = new TH2F("toftray_matchflag","toftray_matchflag;tofTray; tofmatchflag;",120,0,120,4,0,4);
  mtoftray_module = new TH2F("toftray_module","toftray_module;tofTray; module",120,0,120,40,0,40);


  //--bemc-----
  mbTowId_P= new TH2F("bTowId_P","bTowId_P; p (GeV/c); BEMC btowId;",200,0,20,200,-100,5000);
  mbTowId2_P= new TH2F("bTowId2_P","bTowId2_P; p (GeV/c); BEMC btowId2;",200,0,20,10,-1.5,8.5);
  mbTowId3_P= new TH2F("bTowId3_P","bTowId3_P; p (GeV/c); BEMC btowId3;",200,0,20,10,-1.5,8.5);
  mBEMCe0_P= new TH2F("BEMCe0_P","BEMCe0_P; p (GeV/c); BEMC e0(GeV);",200,0,20,200,0,10.0);
  mBEMCe1_P= new TH2F("BEMCe1_P","BEMCe1_P; p (GeV/c); BEMC e1(GeV);",200,0,20,200,0,10.0);
  mBEMCe2_P= new TH2F("BEMCe2_P","BEMCe2_P; p (GeV/c); BEMC e2(GeV);",200,0,20,200,0,10.0);
  mBEMCe3_P= new TH2F("BEMCe3_P","BEMCe3_P; p (GeV/c); BEMC e3(GeV);",200,0,20,200,0,10.0);
  mBEMCE_P= new TH2F("BEMCE_P","BEMCE_P; p (GeV/c); BEMC E (GeV);",200,0,20,200,0,40.0);
  mBEMCadc0_P= new TH2F("BEMCadc0_P","BEMCadc0_P; p (GeV/c); BEMCadc0;",200,0,20,200,0,2000);
  mBEMCzdist_P= new TH2F("BEMCzdist_P","BEMCzdist_P; p (GeV/c); BEMCzdist (cm);",200,0,20,200,-30,30);
  mBEMCphidist_P= new TH2F("BEMCphidist_P","BEMCphidist_P;p (GeV/c); phidist;",200,0,20,200,-0.2,0.2);
  mBEMCneta_P= new TH2F("BEMCneta_P","BEMCneta_P;p (GeV/c); neta;",200,0,20,10,0,10);
  mBEMCnphi_P= new TH2F("BEMCnphi_P","BEMCnphi_P;p (GeV/c); nphi;",200,0,20,10,0,10);
  mBEMCneta_nphi= new TH2F("BEMCneta_nphi","BEMCneta_nphi;neta;nphi;",10,0,4,10,0,4);

  //---------------------mtd------------------------
  mmtdbgcell_module = new TH2F("mtdbgcell_module","backleg*12+cell vs module; module; backleg*12+cell;",7,0,7,400,0,400);
  mmtddeltaT_pt     = new TH2F("mtddeltaT_pt","mtddeltaT_pt; pt (GeV/c); mtdDT;",200,0,20,200,-10,20);
  mmtddeltaZ_pt = new TH2F("mtddeltaZ_pt","detalZ VS pt; pt (GeV/c); mtdDZ;",200,0,20,200,-200,200);
  mmtddeltaY_pt = new TH2F("mtddeltaY_pt","detalY VS pt; pt (GeV/c); mtdDY;",200,0,20,200,-100,100);
  mmtdBeta_P   = new TH2F("mtdBeta_P","mtdBeta VS P; p (GeV/c); mtdBeta;",200,0,20,200,0.95,1.12);
  mmtdmatchflag_channel = new TH2F("mtdmatchflag_channel","mtdmatchflag_channel; channel; mtdmatchflag;",300,0,300,5,0,5);
  mmtddeltaT_channel = new TH2F("mtddeltaT_channel","mtddeltaT_channel; channel; mtdDT;",300,0,1800,300,-5,15);
  mmtddeltaY_channel = new TH2F("mtddeltaY_channel","mtddeltaY_channel; channel; mtdDY;",300,0,1800,300,-100,100);
  mmtddeltaZ_channel = new TH2F("mtddeltaZ_channel","mtddeltaZ_channel; channel; mtdDZ;",300,0,1800,300,-200,200);

  mmtddeltaT_channel1 = new TH2F("mtddeltaT_channel1","mtddeltaT_channel_300-500; channel; mtdDT;",200,300,500,200,-5,15);
  mmtddeltaY_channel1 = new TH2F("mtddeltaY_channel1","mtddeltaY_channel_300-500; channel; mtdDY;",200,300,500,200,-100,100);
  mmtddeltaZ_channel1 = new TH2F("mtddeltaZ_channel1","mtddeltaZ_channel_300-500; channel; mtdDZ;",200,300,500,200,-200,200);
  mmtdBeta_channel   = new TH2F("mtdBeta_channel","mtdBeta_channel; channel; mtdbeta;",300,0,1800,300,0.95,1.12);





  

  //run by run QA
  hTPCVxvsRunIndex = new TH2F("hTPCVxvsRunIndex","hTPCVxvsRunIndex;runIndex;TPC Vx (cm)",mTotalRun,0,mTotalRun,100,-10,10);
  hTPCVyvsRunIndex = new TH2F("hTPCVYvsRunIndex","hTPCVYvsRunIndex;runIndex;TPC VY (cm)",mTotalRun,0,mTotalRun,100,-10,10);
  hTPCVzvsRunIndex = new TH2F("hTPCVzvsRunIndex","hTPCVzvsRunIndex;runIndex;TPC Vz (cm)",mTotalRun,0,mTotalRun,400,-10,10);
  hVPDVzvsRunIndex = new TH2F("hVPDVzvsRunIndex","hVPDVzvsRunIndex;runIndex;VPD Vz (cm)",mTotalRun,0,mTotalRun,160,-40,40);
  hDeltaZvsRunIndex = new TH2F("hDeltaZvsRunIndex","hDeltaZvsRunIndex;runIndex; Vz_{TPC} - Vz_{VPD} (cm)",mTotalRun,0,mTotalRun,200,-10,10);     


  hRefMultvsRunIndex = new TH2F("hRefMultvsRunIndex","hRefMultvsRunIndex;runIndex; refMult",mTotalRun,0,mTotalRun,200,0,500);
  hGRefMultvsRunIndex = new TH2F("hGRefMultvsRunIndex","hGRefMultvsRunIndex;runIndex; GrefMult",mTotalRun,0,mTotalRun,200,0,500);
  //  hRefMultCorrvsRunIndex = new TH2F("hRefMultCorrvsRunIndex","hRefMultCorrvsRunIndex;runIndex; refMultCorr",mTotalRun,0,mTotalRun,200,0,1000);
  hnPrimaryvsRunIndex = new TH2F("hnPrimaryvsRunIndex","hnPrimaryvsRunIndex;runIndex; nPrimary",mTotalRun,0,mTotalRun,300,0,1200);
  hnGlobalvsRunIndex = new TH2F("hnGlobalvsRunIndex","hnGlobalvsRunIndex;runIndex; nGlobal",mTotalRun,0,mTotalRun,400,0,4000);


  hZDCXvsRunIndex = new TH2F("hZDCXvsRunIndex","hZDCXvsRunIndex;runIndex;zdcRate (KHz)",mTotalRun,0,mTotalRun,200,0.,60);
  hBBCXvsRunIndex = new TH2F("hBBCXvsRunIndex","hBBCXvsRunIndex;runIndex;bbcRate (KHz)",mTotalRun,0,mTotalRun,200,0,60);



  hnhitsfitvsRunIndex = new TH2F("hnhitsfitvsRunIndex","hhitsfitvsRunIndex;runIndex;nhitsfit",mTotalRun,0,mTotalRun,50,0,50);
  hnhitsdedxvsRunIndex = new TH2F("hnhitsdedxvsRunIndex","hhitsdedxvsRunIndex;runIndex;nhitsdedx",mTotalRun,0,mTotalRun,50,0,50);

  hPtvsRunIndex = new TH2F("hPtvsRunIndex","hPtvsRunIndex;runIndex;p_{T} (GeV/c)",mTotalRun,0,mTotalRun,100,0,10);
  hEtavsRunIndex = new TH2F("hEtavsRunIndex","hEtavsRunIndex;runIndex;#eta",mTotalRun,0,mTotalRun,200,-2,2);
  hPhivsRunIndex = new TH2F("hPhivsRunIndex","hPhivsRunIndex;runIndex;#phi",mTotalRun,0,mTotalRun,360,0,2*TMath::Pi());     
  hDedxvsRunIndex = new TH2F("hDedxvsRunIndex","hDedxvsRunIndex;runIndex;dE/dx (KeV/cm)",mTotalRun,0,mTotalRun,200,0,10);  

  hNSigmaEvsRunIndex = new TH2F("hNSigmaEvsRunIndex","hNSigmaEvsRunIndex;runIndex;n#sigma_{e}",mTotalRun,0,mTotalRun,600,-30-1.e-6,30-1.e-6);
  hNSigmaPivsRunIndex = new TH2F("hNSigmaPivsRunIndex","hNSigmaPivsRunIndex;runIndex;n#sigma_{#pi}",mTotalRun,0,mTotalRun,600,-30-1.e-6,30-1.e-6);
  hNSigmaKvsRunIndex = new TH2F("hNSigmaKvsRunIndex","hNSigmaKvsRunIndex;runIndex;n#sigma_{k}",mTotalRun,0,mTotalRun,600,-30-1.e-6,30-1.e-6);
  hNSigmaPvsRunIndex = new TH2F("hNSigmaPvsRunIndex","hNSigmaPvsRunIndex;runIndex;n#sigma_{p}",mTotalRun,0,mTotalRun,600,-30-1.e-6,30-1.e-6);


  //hHTth0vsRunIndex = new TH2F("hHTTH0vsRunIndex","hHTth0vsRunIndex;runIndex; HT_TH[0]",mTotalRun,0,mTotalRun, 50,-0.5,49.5);
  //hHTth1vsRunIndex = new TH2F("hHTTH1vsRunIndex","hHTth1vsRunIndex;runIndex; HT_TH[1]",mTotalRun,0,mTotalRun,50,-0.5,49.5);
  //hHTth2vsRunIndex = new TH2F("hHTTH2vsRunIndex","hHTth2vsRunIndex;runIndex; HT_TH[2]",mTotalRun,0,mTotalRun,50,-0.5,49.5);
  //hHTth3vsRunIndex = new TH2F("hHTTH3vsRunIndex","hHTth3vsRunIndex;runIndex; HT_TH[3]",mTotalRun,0,mTotalRun,50,-0.5,49.5);



  hntofmatchvsRunIndex = new TH2F("hntofmatchvsRunIndex","hntofmatchvsRunIndex;runIndex; ntofmatch",mTotalRun,0,mTotalRun,200,0,1000);
  hnbemcmatchvsRunIndex = new TH2F("hnbemcmatchvsRunIndex","hnbemcmatchvsRunIndex;runIndex; nbemcmatch",mTotalRun,0,mTotalRun,200,0,1000);
  hnmtdmatchvsRunIndex = new TH2F("hnmtdmatchvsRunIndex","hnmtdmatchvsRunIndex;runIndex; nmtdmatch",mTotalRun,0,mTotalRun,20,0,20);
  hnmtdhitsvsRunIndex = new TH2F("hnmtdhitsvsRunIndex","hnmtdhitsvsRunIndex;runIndex; nmtdhits",mTotalRun,0,mTotalRun,50,0,50);
  hnhftmatchvsRunIndex = new TH2F("hnhftmatchvsRunIndex","hnhftmatchvsRunIndex;runIndex; nhftmatch",mTotalRun,0,mTotalRun,200,0,500);


  hntofmatchOVnGRvsRunIndex = new TH2F("hntofmatchOVnGRvsRunIndex","hntofmatchOVnGRvsRunIndex;runIndex; ntofmatch/nGreff",mTotalRun,0,mTotalRun,200,0,1);
  hnbemcmatchOVnGRvsRunIndex = new TH2F("hnbemcmatchOVnGRvsRunIndex","hnbemcmatchOVnGRvsRunIndex;runIndex; nbemcmatch/nGreff",mTotalRun,0,mTotalRun,200,0,1);
  hnmtdhitsOVnGRvsRunIndex = new TH2F("hnmtdhitsOVnGRvsRunIndex","hnmtdhits/NGtrks vsRunIndex;runIndex; nmtdhits/nGreff",mTotalRun,0,mTotalRun,100,0,1);
  hnmtdmatchOVnGRvsRunIndex = new TH2F("hnmtdmatchOVnGRvsRunIndex","hnmtdmatchOVnGRvsRunIndex;runIndex; nmtdmatch/nGreff",mTotalRun,0,mTotalRun,100,0,1);
  hnhftmatchOVnGRvsRunIndex = new TH2F("hnhftmatchOVnGRvsRunIndex","hnhftmatchOVnGRvsRunIndex;runIndex; nhftmatch/nGreff",mTotalRun,0,mTotalRun,200,0,1);

  hntofmatchOVnRvsRunIndex = new TH2F("hntofmatchOVnRvsRunIndex","hntofmatchOVnRvsRunIndex;runIndex; ntofmatch/nReff",mTotalRun,0,mTotalRun,200,0,1);
  hnbemcmatchOVnRvsRunIndex = new TH2F("hnbemcmatchOVnRvsRunIndex","hnbemcmatchOVnRvsRunIndex;runIndex; nbemcmatch/nReff",mTotalRun,0,mTotalRun,200,0,1);
  hnmtdhitsOVnRvsRunIndex = new TH2F("hnmtdhitsOVnRvsRunIndex","hnmtdhits/NPtrks vsRunIndex;runIndex; nmtdhits/nReffs",mTotalRun,0,mTotalRun,100,0,1);
  hnmtdmatchOVnRvsRunIndex = new TH2F("hnmtdmatchOVnRvsRunIndex","hnmtdmatchOVnRvsRunIndex;runIndex; nmtdmatch/nReff",mTotalRun,0,mTotalRun,100,0,1);
  hnhftmatchOVnRvsRunIndex = new TH2F("hnhftmatchOVnRvsRunIndex","hnhftmatchOVnRvsRunIndex;runIndex; nhftmatch/nReff",mTotalRun,0,mTotalRun,100,0,1);


  // run by run after pid cuts       
  hNtofElectronvsRunIndex = new TH2F("hNtofElectronvsRunIndex","hNtofElectronvsRunIndex;runIndex; Ntofelectron",mTotalRun,0,mTotalRun,100,0,30);
  hNbemcElectronvsRunIndex = new TH2F("hNbemcElectronvsRunIndex","hNbemcElectronvsRunIndex;runIndex; Nbemcelectron",mTotalRun,0,mTotalRun,100,0,30);
  hNmuonvsRunIndex = new TH2F("hNmuonvsRunIndex","hNmuonvsRunIndex;runIndex; n#muon",mTotalRun,0,mTotalRun,10,0,10);
  // run by run after pid cuts       
  hNtofElectronOVnGRvsRunIndex = new TH2F("hNtofElectronOVnGRvsRunIndex","hNtofElectronOVnGRvsRunIndex;runIndex; Ntofelectron/nGreff",mTotalRun,0,mTotalRun,100,0,1);
  hNbemcElectronOVnGRvsRunIndex = new TH2F("hNbemcElectronOVnGRvsRunIndex","hNbemcElectronOVnGRvsRunIndex;runIndex; Nbemcelectron/nGreff",mTotalRun,0,mTotalRun,100,0,1);
  hNmuonOVnGRvsRunIndex = new TH2F("hNmuonOVnGRvsRunIndex","hNmuonOVnGRvsRunIndex;runIndex; n#muon/nGreff",mTotalRun,0,mTotalRun,10,0,1);


  hNtofElectronOVnRvsRunIndex = new TH2F("hNtofElectronOVnRvsRunIndex","hNtofElectronOVnRvsRunIndex;runIndex; Ntofelectron/nReff",mTotalRun,0,mTotalRun,100,0,1);
  hNbemcElectronOVnRvsRunIndex = new TH2F("hNbemcElectronOVnRvsRunIndex","hNbemcElectronOVnRvsRunIndex;runIndex; Nbemcelectron/nReff",mTotalRun,0,mTotalRun,100,0,1);
  hNmuonOVnRvsRunIndex = new TH2F("hNmuonOVnRvsRunIndex","hNmuonOVnRvsRunIndex;runIndex; n#muon/nReff",mTotalRun,0,mTotalRun,10,0,1);


  hDcavsRunIndex = new TH2F("hDcavsRunIndex","hDcavsRunIndex;runIndex;dca (cm)",mTotalRun,0,mTotalRun,300,0,10);
  hDcawHFTvsRunIndex = new TH2F("hDcawHFTvsRunIndex","hDcawHFTvsRunIndex;runIndex;dca (cm)",mTotalRun,0,mTotalRun,300,0,10);
  hDcaXYvsRunIndex = new TH2F("hDcaXYvsRunIndex","hDcaXYvsRunIndex;runIndex;dcaXY (cm)",mTotalRun,0,mTotalRun,300,-10,10);
  hDcaXYwHFTvsRunIndex = new TH2F("hDcaXYwHFTvsRunIndex","hDcaXYwHFTvsRunIndex;runIndex;dcaXY (cm)",mTotalRun,0,mTotalRun,300,-10,10);
  hDcaZvsRunIndex = new TH2F("hDcaZvsRunIndex","hDcaZvsRunIndex;runIndex;dcaZ (cm)",mTotalRun,0,mTotalRun,300,-10,10);
  hDcaZwHFTvsRunIndex = new TH2F("hDcaZwHFTvsRunIndex","hDcaZwHFTvsRunIndex;runIndex;dcaZ (cm)",mTotalRun,0,mTotalRun,300,-10,10);
  

    /*
    //after pid QA pi,k,p,e,mu
    htofPiondcavsPT = new TH2F("htofPiondcavsPT","tofPiondca vs PT;pt(GeV/c); Piondca(cm);",200,0,20,250,0,10);
    htofKaondcavsPT = new TH2F("htofKaondcavsPT","tofKaondca vs PT;pt(GeV/c); Kaondca(cm);",200,0,20,250,0,10);
    htofProtondcavsPT  = new TH2F("htofProtondcavsPT","tofProtondca vs PT;pt(GeV/c); Protondca(cm);",200,0,20,250,0,10);
    htofElectrondcavsPT = new TH2F("htofElectrondcavsPT","tofElectron vs PT;pt(GeV/c); Electrondca(cm);",200,0,20,250,0,10);

    hbemcElectrondcavsPT = new TH2F("hbemcElectrondcavsPT","bemcElectron vs PT;pt(GeV/c); Electrondca(cm);",200,0,20,250,0,10);
    hMuondcavsPT = new TH2F("hMuondcavsPT","Muon vs PT;pt(GeV/c); Muondca(cm);",200,0,20,250,0,10);


    //after pid QA pi,k,p,e,mu-----wHFT----
    htofPiondcawHFTvsPT = new TH2F("htofPiondcawHFTvsPT","tofPiondcawHFT vs PT;pt(GeV/c); Piondca(cm);",200,0,20,250,0,3);
    htofKaondcawHFTvsPT = new TH2F("htofKaondcawHFTvsPT","tofKaondcawHFT vs PT;pt(GeV/c); Kaondca(cm);",200,0,20,250,0,3);
    htofProtondcawHFTvsPT  = new TH2F("htofProtondcawHFTvsPT","tofProtondcawHFTvsPT;pt(GeV/c); Protondca(cm);",200,0,20,250,0,3);
    htofElectrondcawHFTvsPT = new TH2F("htofElectrondcawHFTvsPT","tofElectrondcawHFTvsPT;pt(GeV/c); Electrondca(cm);",200,0,20,250,0,3);
    hbemcElectrondcawHFTvsPT = new TH2F("hbemcElectrondcawHFTvsPT","bemcElectrondcawHFTvsPT;pt(GeV/c); Electrondca(cm);",200,0,20,250,0,3);
    hMuondcawHFTvsPT = new TH2F("hMuondcawHFTvsPT","MuondcawHFTvsPT;pt(GeV/c); Muondca(cm);",200,0,20,250,0,3);

    //---only tpc pid----pi/k/p----------
    htpcPiondcavsPT = new TH2F("htpcPiondcavsPT","Piondca vs PT (|n#sigma_{#pi}|<2);pt(GeV/c); Piondca(cm);",200,0,20,250,0,3);
    htpcKaondcavsPT = new TH2F("htpcKaondcavsPT","Kaondca vs PT (|n#sigma_{k}|<2);pt(GeV/c); Kaondca(cm);",200,0,20,250,0,3);
    htpcProtondcavsPT = new TH2F("htpcProtondcavsPT","Protondca vs PT (|n#sigma_{p}|<2);pt(GeV/c); Protondca(cm);",200,0,20,250,0,3);


    //--- tpc pid + tof----pi/k/p------pt<2----
    htpctofPiondcavsPT = new TH2F("htpctofPiondcavsPT","PiondcavsPT (|n#sigma_{#pi}|<2&tofm2);pt(GeV/c); Piondca(cm);",200,0,20,250,0,10);
    htpctofKaondcavsPT = new TH2F("htpctofKaondcavsPT","KdcavsPT (|n#sigma_{k}|<2&tofm2);pt(GeV/c); Kaondca(cm);",200,0,20,250,0,10);
    htpctofProtondcavsPT = new TH2F("htpctofProtondcavsPT","PdcavsPT (|n#sigma_{p}|<2&tofm2);pt(GeV/c); Protondca(cm);",200,0,20,250,0,10);

    //---only tpc pid----pi/k/p---wHFT-------
    htpcPiondcawHFTvsPT = new TH2F("htpcPiondcawHFTvsPT","PiondcawHFT vs PT (|n#sigma_{#pi}|<2);pt(GeV/c); PiondcawHFT(cm);",200,0,20,250,0,10);
    htpcKaondcawHFTvsPT = new TH2F("htpcKaondcawHFTvsPT","KaondcawHFT vs PT (|n#sigma_{k}|<2);pt(GeV/c); KaondcawHFT(cm);",200,0,20,250,0,10);
    htpcProtondcawHFTvsPT = new TH2F("htpcProtondcawHFTvsPT","ProtondcawHFT vs PT (|n#sigma_{p}|<2);pt(GeV/c); Protondca(cm);",200,0,20,250,0,10);


    //--- tpc pid + tof----pi/k/p----wHFT--pt<2----
    htpctofPiondcawHFTvsPT = new TH2F("htpctofPiondcawHFTvsPT","PiondcawHFTvsPT (|n#sigma_{#pi}|<2&tofm2);pt(GeV/c); Piondca(cm);",200,0,20,250,0,3);
    htpctofKaondcawHFTvsPT = new TH2F("htpctofKaondcawHFTvsPT","KdcawHFTvsPT (|n#sigma_{k}|<2&tofm2);pt(GeV/c); Kaondca(cm);",200,0,20,250,0,3);
    htpctofProtondcawHFTvsPT = new TH2F("htpctofProtondcawHFTvsPT","PdcawHFTvsPT (|n#sigma_{p}|<2&tofm2);pt(GeV/c); Protondca(cm);",200,0,20,250,0,3);
    */

}// er chen

//----------------------------------------------------------------------------- 
void StPicoQAMaker::Clear(Option_t *opt) {

  return StMaker::Clear(opt);
}

//----------------------------------------------------------------------------- 
Int_t StPicoQAMaker::Make() {

  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  mPicoDst = mPicoDstMaker->picoDst();

  if(!mPicoDst) {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
  }


  StPicoEvent* event=mPicoDst->event();
  //=================zaochen add====================
  hNEvents->Fill(0);
  if(!event) return kStOK;
  //=================event selection================
  int runId = mPicoDst->event()->runId();
  for(int i=0;i<mNBadRuns;i++){
    if(runId==mBadRuns[i]) return kStOK;
  }


  //if( ! mPicoDst->event()->isMinBias()) return kStOK;
  //if( ! isBHT1( mPicoDst->event() ) ) return kStOK;
  //=================event selection================

  hNEvents->Fill(2);

  if(fillhistflag){	

    Int_t RUNID = mPicoDst->event()->runId();
    map<Int_t, Int_t>::iterator iter = mTotalRunId.find(RUNID);
    if(iter != mTotalRunId.end())
      runIndex = iter->second;
    else{
      runIndex = -1;
      cout<<"sorry, no runNumber in runNumber list"<<endl;
      cout<<"the RUNID:="<<RUNID<<endl;
    }

    if(runIndex<0)return kStOK;

  }//

  //int triggerWord = event->triggerWord();
  //	if(triggerWord>>19 & 0x1
  Int_t triggerWORD=mPicoDst->event()->triggerWord();
  //for(Int_t i=0; i<25; i++){
  //   if( triggerWORD>>i & 0x1 ) htriggerindex->Fill(i);

  //}

  //	runidfiles<<mPicoDst->event()->runId()<<endl;

  //cout<<"run number is:"<<mPicoDst->event()->runId()<<endl;

  //======================ALL QA PLOTS ===============




  //---------------event information----------------------------

  Double_t vzvpd=mPicoDst->event()->vzVpd();
  Double_t vztpc=mPicoDst->event()->primaryVertex().z();
  Double_t vxtpc=mPicoDst->event()->primaryVertex().x();
  Double_t vytpc=mPicoDst->event()->primaryVertex().y();
  Double_t dvz=vztpc-vzvpd;
  Double_t vr=sqrt(vxtpc*vxtpc+vytpc*vytpc);
  if(fabs(vxtpc)<1.0e-5)return kStOK;
  if(fabs(vytpc)<1.0e-5)return kStOK;
  if(fabs(vztpc)<1.0e-5)return kStOK;

  if(fillhistflag){
    //mVz_tpc->Fill(vztpc);
    //mVz_vpd->Fill(vzvpd);
    //mdVz->Fill(dvz);  
    mVz_vpdtpc->Fill(vztpc,vzvpd);
    mdVz_tpcVz->Fill(vztpc,dvz);
    mVxy->Fill( mPicoDst->event()->primaryVertex().x(), mPicoDst->event()->primaryVertex().y() );
    mVRvsVZ->Fill(vztpc, vr);
  }//


  //Int_t Nptrks = mPicoDst->numberOfTracks();
  Float_t Ranking = mPicoDst->event()->ranking();
  Float_t zdcx = mPicoDst->event()->ZDCx();
  Float_t bbcx = mPicoDst->event()->BBCx();
  zdcx=zdcx/1000.;
  bbcx=bbcx/1000.;
  Int_t Nmtdhits = mPicoDst->numberOfMtdHits();
  Int_t Ntofhits = mPicoDst->numberOfBTofHits();
  Int_t NRefmultPos=mPicoDst->event()->refMultPos();
  Int_t NRefmultNeg=mPicoDst->event()->refMultNeg();
  Int_t NGnremult=mPicoDst->event()->grefMult();
  Int_t NRefmult=mPicoDst->event()->refMult();
  Int_t NGtrks = mPicoDst->event()->numberOfGlobalTracks();
  Int_t Ntofmatch = mPicoDst->event()->nBTOFMatch();
  //	Int_t Nbemchits = mPicoDst->event()->
  Int_t Nbemcmatch = mPicoDst->event()->nBEMCMatch();
  if(fillhistflag){


    //run by run
    if( Ispasseventcuts(event) ){

      //hHTth0vsRunIndex->Fill(runIndex, mPicoDst->event()->ht_th(0));
      //hHTth1vsRunIndex->Fill(runIndex, mPicoDst->event()->ht_th(1));
      //hHTth2vsRunIndex->Fill(runIndex, mPicoDst->event()->ht_th(2));
      //hHTth3vsRunIndex->Fill(runIndex, mPicoDst->event()->ht_th(3));
      
         hTPCVzvsRunIndex->Fill(runIndex, vztpc);
         hTPCVxvsRunIndex->Fill(runIndex, vxtpc);
         hTPCVyvsRunIndex->Fill(runIndex, vytpc);
         hVPDVzvsRunIndex->Fill(runIndex, vzvpd);
         hDeltaZvsRunIndex->Fill(runIndex, dvz);
         hZDCXvsRunIndex->Fill(runIndex,  zdcx);
         hBBCXvsRunIndex->Fill(runIndex,  bbcx);
         hRefMultvsRunIndex->Fill(runIndex, NRefmult);
         hGRefMultvsRunIndex->Fill(runIndex, NGnremult);

         hnGlobalvsRunIndex->Fill(runIndex, NGtrks);
         hnbemcmatchvsRunIndex->Fill(runIndex,Nbemcmatch);
         
    }//


  }
  //----------track information------------------------  
  Int_t numberoftracks = mPicoDst->numberOfTracks();
  //if(fillhistflag){	
  //   mNptracks->Fill(numberoftracks);
  //   mNgtracks->Fill(NGtrks);
  //}//

  //tVzTPC=mPicoDst->event()->primaryVertex().z()

  StThreeVectorF vertexPos;
  vertexPos = mPicoDst->event()->primaryVertex();

  Int_t Nprimarytracks=0;
  Int_t ntofmatchcount=0;
  Int_t nmtdmatchcount=0;
  Int_t nbemcmatchcount=0;
  Int_t nhftmatchcount=0;
  Int_t nmuons=0;
  Int_t ntofelecton=0;
  Int_t nbemcelectron=0;
  Float_t particleM[3]={0.938,0.140,0.494};
  for(int i=0; i<numberoftracks; i++){
    StPicoTrack* track=(StPicoTrack*) mPicoDst->track(i);
    StPhysicalHelixD helix = track->helix();

    StThreeVectorF dcaPoint = helix.at( helix.pathLength(vertexPos.x(), vertexPos.y())  );
    double dcamag= (dcaPoint-vertexPos).mag();
    StThreeVectorF dcaP = helix.momentumAt( vertexPos.x(),vertexPos.y() );
    double dcaXY= ( (dcaPoint-vertexPos).x()*dcaP.y()-(dcaPoint-vertexPos).y()*dcaP.x() )/dcaP.perp();
    double dcaZ= dcaPoint.z() - vertexPos.z();

    Double_t meta,mpt,mphi,mcharge,mdedx;
    //change to global track
    meta=track->pMom().pseudoRapidity();
    if(track->pMom().mag()!=0) Nprimarytracks++;
    mphi=RotatePhi(track->pMom().phi());
    mpt=track->pMom().perp();
    mcharge=track->charge();
    mdedx=track->dEdx();
    if(mcharge==0||meta==0||mphi==0||mdedx==0)continue;
    if(track->isHFTTrack()){
      nhftmatchcount++;
    }

    Float_t mmomentum=track->pMom().mag();
    Double_t rationhits = 0.00;
    rationhits=(double)track->nHitsFit()/track->nHitsMax();
    // cout<<"track->nHitsMax()="<<track->nHitsMax()<<"track->nHitsFit()"<<track->nHitsFit()<<endl;
    // cout<<"ratiohits"<<rationhits<<endl;
    Double_t nsigpi=track->nSigmaPion();
    Double_t nsigk=track->nSigmaKaon();
    Double_t nsigp=track->nSigmaProton();	     
    if(fillhistflag){

      mtrkpt->Fill( mpt );
      mtrketa->Fill( meta );
      mtrkphi->Fill( mphi );
      mtrketaphi->Fill( mphi,meta );
      mtrkdca->Fill( dcamag );
      mtrkdcaXY->Fill( dcaXY );
      mtrkdcaZ->Fill( dcaZ );

      mnhitsfit->Fill( track->nHitsFit() );   
      mnhitsdedx->Fill( track->nHitsDedx() );
      mnhitsfitRatio->Fill(rationhits);
      if( track->isHFTTrack() ){
        mtrkdcawHFT->Fill( dcamag );
        mtrkdcaXYwHFT->Fill( dcaXY );
        mtrkdcaZwHFT->Fill( dcaZ );
        mnhitsfitwHFT->Fill( track->nHitsFit() );
        mnhitsdedxwHFT->Fill( track->nHitsDedx() );
        mnhitsfitRatiowHFT->Fill(rationhits);
      }
      mnsigmaPI->Fill( track->nSigmaPion() );
      mnsigmaP->Fill( track->nSigmaProton() );
      mnsigmaK->Fill( track->nSigmaKaon() );
      mnsigmaE->Fill( track->nSigmaElectron() );

      //tpc pid on pi, k, p
      /*if(fabs(nsigpi)<2.0) htpcPiondcavsPT->Fill(mpt,dcamag);
        if(fabs(nsigk)<2.0)  htpcKaondcavsPT->Fill(mpt,dcamag);
        if(fabs(nsigp)<2.0)  htpcProtondcavsPT->Fill(mpt,dcamag);
      //tpc pid on pi, k, p---wHFT
      if(track->isHFTTrack()){
      if(fabs(nsigpi)<2.0) htpcPiondcawHFTvsPT->Fill(mpt,dcamag);
      if(fabs(nsigk)<2.0)  htpcKaondcawHFTvsPT->Fill(mpt,dcamag);
      if(fabs(nsigp)<2.0)  htpcProtondcawHFTvsPT->Fill(mpt,dcamag);
      }*/

      mtrketa_pt->Fill(mpt*mcharge,meta);
      mtrkphi_pt->Fill(mpt*mcharge,mphi);
      if(mcharge>0)mtrketa_phiPos->Fill(mphi,meta);
      if(mcharge<0)mtrketa_phiNeg->Fill(mphi,meta);

      mtrkdca_pt->Fill(mpt, dcamag);
      mtrkdcaZ_pt->Fill(mpt,dcaZ);
      mtrkdcaXY_pt->Fill(mpt,dcaXY);

      if(track->isHFTTrack()){
        mtrkdcawHFT_pt->Fill(mpt, dcamag);
        mtrkdcaZwHFT_pt->Fill(mpt,dcaZ);
        mtrkdcaXYwHFT_pt->Fill(mpt,dcaXY);
      }

      mnhitsfit_pt->Fill(mpt,track->nHitsFit());
      mnhitsdedx_pt->Fill(mpt,track->nHitsDedx());

      //if(track->nHitsMax()>0) rationhits=track->nHitsFit()/track->nHitsMax()*1.;

      mnhitsRatio_pt->Fill(mpt,rationhits);
      if(track->isHFTTrack()) mnhitsRatio_ptwHFT->Fill(mpt,rationhits);

      mdedx_P->Fill(mmomentum*mcharge,track->dEdx());
      mnsigmaPI_P->Fill(mmomentum,track->nSigmaPion());
      mnsigmaP_P->Fill(mmomentum,track->nSigmaProton());
      mnsigmaE_P->Fill(mmomentum,track->nSigmaElectron());
      mnsigmaK_P->Fill(mmomentum,track->nSigmaKaon());



      //------tpc information end-----

      

      //run by run QA
      if( Ispasseventcuts(event) ){
      hPtvsRunIndex->Fill(runIndex,mpt);
      hEtavsRunIndex->Fill(runIndex,meta);
      hPhivsRunIndex ->Fill(runIndex,mphi);
      hnhitsfitvsRunIndex ->Fill(runIndex,track->nHitsFit());
      hnhitsdedxvsRunIndex ->Fill(runIndex,track->nHitsDedx());
      hDedxvsRunIndex->Fill(runIndex,track->dEdx());

      hDcavsRunIndex->Fill(runIndex,dcamag);
      hDcaZvsRunIndex->Fill(runIndex,dcaZ);
      hDcaXYvsRunIndex->Fill(runIndex,dcaXY);
      if( track->isHFTTrack() ){
      hDcawHFTvsRunIndex->Fill(runIndex,dcamag);
      hDcaZwHFTvsRunIndex->Fill(runIndex,dcaZ);
      hDcaXYwHFTvsRunIndex->Fill(runIndex,dcaXY);
      }

      hNSigmaEvsRunIndex->Fill(runIndex,track->nSigmaElectron());
      hNSigmaPivsRunIndex->Fill(runIndex,track->nSigmaPion());
      hNSigmaKvsRunIndex->Fill(runIndex,track->nSigmaKaon());
      hNSigmaPvsRunIndex->Fill(runIndex,track->nSigmaProton());
      }
      
      }//
      


      Int_t tofpidid=track->bTofPidTraitsIndex();
      if(tofpidid>0){
        ntofmatchcount++;
        StPicoBTofPidTraits* btofpidtrait=(StPicoBTofPidTraits*) mPicoDst->btofPidTraits(tofpidid);

        //------tof information start----------
        Float_t tofbeta=btofpidtrait->btofBeta();
        /*
           if(fillhistflag){		   
           minvsBeta_P->Fill(mmomentum,1/tofbeta);
           if(tofbeta>0){
           Double_t tofm2=mmomentum*mmomentum*( 1.0/(tofbeta*tofbeta)-1.0);
           mtofM2_P->Fill(mmomentum,tofm2);
           }
           }//
           */


        Int_t tofcellid=   btofpidtrait->btofCellId();
        Int_t toftray= (int)tofcellid/192 + 1;
        Int_t tofmodule= (int)((tofcellid%192)/6.)+1;

        Float_t toflocaly = btofpidtrait->btofYLocal();
        Float_t toflocalz = btofpidtrait->btofZLocal();
        // Float_t tofhitPosx = btofpidtrait->btofHitPos().x();
        // Float_t tofhitPosy = btofpidtrait->btofHitPos().y();
        // Float_t tofhitPosz = btofpidtrait->btofHitPos().z();
        /*
           if(fillhistflag){
           mtoftray_localY->Fill(toftray,toflocaly);
           mtoftray_localZ->Fill(toftray,toflocalz);
           mtoftray_matchflag->Fill(toftray,btofpidtrait->btofMatchFlag());
           mtoftray_module->Fill(toftray,tofmodule);

           }//
           */

        double mnsigE=track->nSigmaElectron();
        /*
           if(fillhistflag){ 
        //mpt>0.2&&mpt<2
        if( mpt>0.2){

        if(mnsigE>-1.5&&mnsigE<3&&fabs(toflocaly)<1.8&&fabs(tofbeta)>0.97&&fabs(tofbeta)<1.03){
        ntofelecton++;

        htofElectrondcavsPT->Fill(mpt,dcamag);
        if(track->isHFTTrack()) {
        htofElectrondcawHFTvsPT->Fill(mpt,dcamag);
        }

        }//electron			      



        //pi----k----p----
        if(tofbeta>0){
        double tofm2=mmomentum*mmomentum*( 1.0/(tofbeta*tofbeta)-1.0);

        if( fabs(tofm2-particleM[0]*particleM[0])<0.05 ){
        htofProtondcavsPT->Fill(mpt,dcamag);

        if(track->isHFTTrack()){
        htofProtondcawHFTvsPT->Fill(mpt,dcamag);
        }

        if(abs(nsigp)<2.0){ 
        htpctofProtondcavsPT->Fill(mpt,dcamag);
        if(track->isHFTTrack()) {
        htpctofProtondcawHFTvsPT->Fill(mpt,dcamag);
        }
        }
        }//p

        if( fabs(tofm2-particleM[1]*particleM[1])<0.01 ){
        htofPiondcavsPT->Fill(mpt, dcamag);

        if(track->isHFTTrack()){
        htofPiondcawHFTvsPT->Fill(mpt,dcamag);
        }

        if(fabs(nsigpi)<2.0){
        htpctofPiondcavsPT->Fill(mpt,dcamag);
        if(track->isHFTTrack()){
        htpctofPiondcawHFTvsPT->Fill(mpt,dcamag);
        }

        }
        }//pion

        if( fabs(tofm2-particleM[2]*particleM[2])<0.02 ){
        htofKaondcavsPT->Fill(mpt, dcamag);

        if(track->isHFTTrack()){
        htofKaondcawHFTvsPT->Fill(mpt,dcamag);
        }

        if( fabs(nsigk)<2.0 ) {
        htpctofKaondcavsPT->Fill(mpt, dcamag);
        if(track->isHFTTrack()){
        htpctofKaondcawHFTvsPT->Fill(mpt,dcamag);
        }

        }
        }//kaon

        }

        }

      }//if fill hist
      */


      }//-------tof information end 


      //-------mtd information start-----------

      Int_t mtdpidtraitsid=track->mtdPidTraitsIndex();
      if(mtdpidtraitsid>=0){
        nmtdmatchcount++;
        StPicoMtdPidTraits* mtdpidtraits=(StPicoMtdPidTraits*) mPicoDst->mtdPidTraits(mtdpidtraitsid);
        Int_t mtdbackleg = mtdpidtraits->backleg();
        Int_t mtdmodule = mtdpidtraits->module();
        Int_t mtdcell   = mtdpidtraits->cell();
        Float_t mtddT = mtdpidtraits->deltaTimeOfFlight();
        Float_t mtddz = mtdpidtraits->deltaZ();
        Float_t mtddy = mtdpidtraits->deltaY();
        Float_t mtdbeta = mtdpidtraits->beta();
        Int_t mtdbgcell=mtdbackleg*12+mtdcell;
        Int_t mtdchannel=(mtdbackleg-1)*60+(mtdmodule-1)*12+mtdcell;// (backleg-1) * 60 + (module-1) * 12 + cell
        /*
           if(fillhistflag){
           mmtdbgcell_module->Fill(mtdmodule,mtdbgcell);
           mmtddeltaT_pt->Fill(mpt,mtddT);
           mmtddeltaZ_pt->Fill(mpt,mtddz);
           mmtddeltaY_pt->Fill(mpt,mtddy);

           mmtdBeta_P->Fill(mmomentum,mtdbeta);
           mmtdmatchflag_channel->Fill(mtdchannel,mtdpidtraits->matchFlag());
           mmtddeltaT_channel->Fill(mtdchannel,mtddT);
           mmtddeltaY_channel->Fill(mtdchannel,mtddy);
           mmtddeltaZ_channel->Fill(mtdchannel,mtddz);
           if(mtdchannel>300&&mtdchannel<500){
           mmtddeltaT_channel1 ->Fill(mtdchannel,mtddT);
           mmtddeltaY_channel1 ->Fill(mtdchannel,mtddy);
           mmtddeltaZ_channel1 ->Fill(mtdchannel,mtddz);
           }

           mmtdBeta_channel->Fill(mtdchannel,mtdbeta);
           }//
           */
      }

      /*
         if(Ismuontrack(track)) {
         nmuons++; 

         if(fillhistflag){	 
         hMuondcavsPT->Fill(mpt,dcamag);

         if(track->isHFTTrack()){
         hMuondcawHFTvsPT->Fill(mpt,dcamag);
         }

         }//		
         }   
         */

      //----------mtd information end------------

      //----------BEMC----------------
      Int_t emcpidtraitsid=track->emcPidTraitsIndex();
      /*
         if(emcpidtraitsid>=0){
         nbemcmatchcount++;
         StPicoEmcPidTraits* emcpidtraits=(StPicoEmcPidTraits*) mPicoDst->emcPidTraits(emcpidtraitsid);
         if(fillhistflag){  
         mbTowId_P->Fill(mmomentum,emcpidtraits->btowId());
         mbTowId2_P->Fill(mmomentum,emcpidtraits->btowId2());
         mbTowId3_P->Fill(mmomentum,emcpidtraits->btowId3());
         mBEMCe0_P->Fill(mmomentum,emcpidtraits->e0());
         mBEMCe1_P->Fill(mmomentum,emcpidtraits->e1());
         mBEMCe2_P->Fill(mmomentum,emcpidtraits->e2());
         mBEMCe3_P->Fill(mmomentum,emcpidtraits->e3());
         mBEMCE_P->Fill(mmomentum,emcpidtraits->e());
         mBEMCadc0_P->Fill(mmomentum,emcpidtraits->adc0());
         mBEMCzdist_P->Fill(mmomentum,emcpidtraits->zDist());
         mBEMCphidist_P->Fill(mmomentum,emcpidtraits->phiDist());
         mBEMCneta_P->Fill(mmomentum,emcpidtraits->nEta());
         mBEMCnphi_P->Fill(mmomentum,emcpidtraits->nPhi());
         mBEMCneta_nphi->Fill(emcpidtraits->nPhi(),emcpidtraits->nEta());
         }//



         float PovE=mmomentum/emcpidtraits->e();

         if( Ispasseventcuts(event) ){
         if(mpt>1.50&&PovE>0.3&&PovE<1.5){
         nbemcelectron++;		      
         if(fillhistflag){  
         hbemcElectrondcavsPT->Fill(mpt,dcamag);
         if(track->isHFTTrack()){
         hbemcElectrondcawHFTvsPT->Fill(mpt,dcamag);
         }
         }//
         }
         }
         }//end of bemc 
         */

  }//loop of all tracks

  Int_t Nptrks=Nprimarytracks;

  if(fillhistflag){
    //event information QA

    mRanking_nPtrks->Fill(Nptrks,Ranking);
    mnPtrks_nGtrks->Fill(NGtrks,Nptrks);
    mnRefMult_nPtrks->Fill(Nptrks,NRefmultPos+NRefmultNeg);
    mnRefMult_nGtrks->Fill(NGtrks,NRefmultPos+NRefmultNeg);
    mnGRefMult_nPtrks->Fill(Nptrks, NGnremult);
    mnGRefMult_nGtrks->Fill(NGtrks, NGnremult);	   

    mnRefMult_nGRefMult->Fill(NGnremult,NRefmultPos+NRefmultNeg);
    mnPtrks_nMtdHits->Fill(Nmtdhits,Nptrks);
    mnPtrks_nTofHits->Fill(ntofmatchcount,Nptrks);
    mnTofHits_nMtdHits->Fill(Nmtdhits,ntofmatchcount);                
    mnTofMatch_nTofHits->Fill(ntofmatchcount,Ntofmatch);  	


    //------------------------------------------------------//
    //---------- make run by run QA of ---------------------//
    //------------------------------------------------------//
    
       if( Ispasseventcuts(event) ){
       hnPrimaryvsRunIndex->Fill(runIndex,Nptrks);

       hntofmatchvsRunIndex->Fill(runIndex,ntofmatchcount);
       hnmtdmatchvsRunIndex->Fill(runIndex,nmtdmatchcount);
       hnbemcmatchvsRunIndex->Fill(runIndex,nbemcmatchcount);
       hnhftmatchvsRunIndex->Fill(runIndex,nhftmatchcount);
       hnmtdhitsvsRunIndex->Fill(runIndex,Nmtdhits);
       hNmuonvsRunIndex->Fill(runIndex,nmuons);
       hNtofElectronvsRunIndex->Fill(runIndex,ntofelecton);
       hNbemcElectronvsRunIndex->Fill(runIndex,nbemcelectron);

       float tofmatchratio= (float)ntofmatchcount/NGnremult;
       float mtdmatchratio=(float)nmtdmatchcount/NGnremult;
       float bemcmatchratio=(float)nbemcmatchcount/NGnremult;
       float hftmatchratio=(float)nhftmatchcount/NGnremult;
       float mtdratio=(float)Nmtdhits/NGnremult;		
       float muonratio=(float)nmuons/NGnremult;
       float toferatio=(float)ntofelecton/NGnremult;
       float bemceratio=(float)nbemcelectron/NGnremult;
       hntofmatchOVnGRvsRunIndex->Fill(runIndex,tofmatchratio);
       hnmtdmatchOVnGRvsRunIndex->Fill(runIndex, mtdmatchratio);
       hnbemcmatchOVnGRvsRunIndex->Fill(runIndex,bemcmatchratio);
       hnhftmatchOVnGRvsRunIndex->Fill(runIndex,hftmatchratio);
       hnmtdhitsOVnGRvsRunIndex->Fill(runIndex,mtdratio);
       hNmuonOVnGRvsRunIndex->Fill(runIndex, muonratio);
       hNtofElectronOVnGRvsRunIndex->Fill(runIndex,toferatio);
       hNbemcElectronOVnGRvsRunIndex->Fill(runIndex,bemceratio);


       float tofmatchratio2= (float)ntofmatchcount/NRefmult;
       float mtdmatchratio2=(float)nmtdmatchcount/NRefmult;
       float bemcmatchratio2=(float)nbemcmatchcount/NRefmult;
       float hftmatchratio2=(float)nhftmatchcount/NRefmult;
       float mtdratio2=(float)Nmtdhits/NRefmult;		
       float muonratio2=(float)nmuons/NRefmult;
       float toferatio2=(float)ntofelecton/NRefmult;
       float bemceratio2=(float)nbemcelectron/NRefmult;
       hntofmatchOVnRvsRunIndex->Fill(runIndex,tofmatchratio2);
       hnmtdmatchOVnRvsRunIndex->Fill(runIndex, mtdmatchratio2);
       hnbemcmatchOVnRvsRunIndex->Fill(runIndex,bemcmatchratio2);
       hnhftmatchOVnRvsRunIndex->Fill(runIndex,hftmatchratio2);
       hnmtdhitsOVnRvsRunIndex->Fill(runIndex,mtdratio2);
       hNmuonOVnRvsRunIndex->Fill(runIndex, muonratio2);
       hNtofElectronOVnRvsRunIndex->Fill(runIndex,toferatio2);
       hNbemcElectronOVnRvsRunIndex->Fill(runIndex,bemceratio2);
       cout << "runIndex: " << runIndex << endl;
       }
  }//
  //  mPicoDst->Print();
  //  mPicoDst->printTracks();

  return kStOK;
}//end of main fucntion

Bool_t StPicoQAMaker::isBHT1(StPicoEvent *event)
{
  int triggerWord = event->triggerWord();
  if(triggerWord>>19 & 0x1 || triggerWord>>20 & 0x1) return true;
  else return false;
}

//-----------------------------------------                                              
Bool_t StPicoQAMaker::isBHT2(StPicoEvent *event)
{
  int triggerWord = event->triggerWord();
  if(triggerWord>>21 & 0x1 || triggerWord>>22 & 0x1) return true;
  else return false;
}

//---------------------------------------------------  
Bool_t StPicoQAMaker::isBHT3(StPicoEvent *event)
{
  int triggerWord = event->triggerWord();
  if(triggerWord>>23 & 0x1 || triggerWord>>24 & 0x1) return true;
  else return false;
}


//------------------------------------------------------------- 

Bool_t StPicoQAMaker::Isgoodtrack(StPicoTrack* track)
{
  double pt,eta,fithitsfrac=0.,chargeq; 
  pt=track->pMom().perp();
  eta=track->pMom().pseudoRapidity();
  fithitsfrac=(track->nHitsFit())/45.0;
  chargeq=track->charge();
  if(pt>0.2&&fabs(eta)<1.0&&track->nHitsFit()>=15&&track->nHitsDedx()>=10&&fithitsfrac>0.52&&fabs(chargeq)>0 ) return true;
  else return false;

}
//-------------------------------------------------------------
Bool_t StPicoQAMaker::Ismuontrack(StPicoTrack* track)
{
  double pt;
  pt=track->pMom().perp();
  Int_t mtdpid=track->mtdPidTraitsIndex();
  if(mtdpid<=0)return false;
  StPicoMtdPidTraits* mtdpidtrait=(StPicoMtdPidTraits*) mPicoDst->mtdPidTraits(mtdpid);      
  double mtddz=mtdpidtrait->deltaZ();
  double mtddt=mtdpidtrait->deltaTimeOfFlight();
  if(track->nSigmaPion()<3.0&&track->nSigmaPion()>-1.0&&pt>1.0&&fabs(mtddz)<20.0&&mtddt>-1.0&&mtddt<1.0) return true;
  else return false;

}

//----------------------------------------------------------------
Bool_t StPicoQAMaker::Ispasseventcuts(StPicoEvent* event)
{
  double vzTPC=event->primaryVertex().z();
  if(fabs(vzTPC)<40)return true;
  else return false;
}

//---------------------------------------------------------------
Double_t StPicoQAMaker::RotatePhi(Double_t phi) const
{
  Double_t outPhi = phi;
  Double_t pi=TMath::Pi();
  while(outPhi<0) outPhi += 2*pi;
  while(outPhi>2*pi) outPhi -= 2*pi;
  return outPhi;
}




//===================================================================



