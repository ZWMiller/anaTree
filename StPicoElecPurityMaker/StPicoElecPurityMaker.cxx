#include "StPicoElecPurityMaker.h"
#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoTrack.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoDstMaker/StPicoEvent.h"
#include "StRoot/StPicoDstMaker/StPicoMtdHit.h"
#include "StRoot/StPicoDstMaker/StPicoConstants.h"
#include "StRoot/StPicoDstMaker/StPicoMtdPidTraits.h"
#include "StRoot/StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StRoot/StPicoDstMaker/StPicoEmcPidTraits.h"
#include "StRoot/StPicoDstMaker/StPicoEmcTrigger.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StDcaGeometry.h"

#include "StThreeVectorF.hh"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TFile.h"
#include <map>
#include<iostream>
#include<fstream>
#include "mBadRunList.h"
#include "mHotTowerList.h"
#include "StBTofUtil/tofPathLength.hh"
#define USHRT_MAX 65535

Bool_t fillhistflag=1;
ofstream runidfiles;
Int_t runIndex;
Int_t randomId;
Int_t mTotalRun = 599;//603
Bool_t DEBUG = kFALSE;

ClassImp(StPicoElecPurityMaker)

  //-----------------------------------------------------------------------------
  StPicoElecPurityMaker::StPicoElecPurityMaker(const char* name, StPicoDstMaker *picoMaker, const char* outName) 
: StMaker(name)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = 0;
  runMode = 0;
  TH1F:: SetDefaultSumw2();//zaochen add
  mOutName = outName;

  mNBadRuns = sizeof(mBadRuns)/sizeof(int);
  mNHotTower1 = sizeof(mHotTower1)/sizeof(int);
  mNHotTower2 = sizeof(mHotTower2)/sizeof(int);
  mNHotTower3 = sizeof(mHotTower3)/sizeof(int);

  // ZWM Add (declare default cut values, can change in Ana with same function)
  SetDefaultCuts();
}

//----------------------------------------------------------------------------- 
StPicoElecPurityMaker::~StPicoElecPurityMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t StPicoElecPurityMaker::Init() {
  cout << " --------------------------------------------------- in purity init ------------------------------------------------------- " << endl;
  cout << " mOutName = " << mOutName.Data() << endl;
  if(mOutName!="") {
    fout = new TFile(mOutName.Data(),"RECREATE");
  }else{
    fout = new TFile("picoQA_test.root","RECREATE");
  }
  DeclareHistograms();

  //runidfiles.open("/star/u/zye20/zye20/zaochen/analysis/run14picoQA/runnumber.txt", ios::out|ios::app);
  TString modeDescription = "p+p (no Centrality)";
  if(runMode == AuAu) modeDescription = "Au+Au (Centrality On)";
  cout << "StElecPurityMaker running in " << modeDescription << " mode." << endl;

  if(fillhistflag){
    //read in the runlist.dat
    ifstream indata;
    indata.open("StRoot/StPicoElecPurityMaker/mTotalRunList15.dat");
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

    for(map<Int_t,Int_t>::iterator iter=mTotalRunId.begin();iter!=mTotalRunId.end();iter++)
      cout<<iter->second<<" \t"<<iter->first<<endl;
    cout<<endl;
    // 
    //read in the runidlist


  }//


  return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StPicoElecPurityMaker::Finish() {
  fout->cd();
  for(int t=0; t<4; t++) // THnSparse must be manually written to file
  {
    cout << "Writing out purity Hists trig: " << t << endl;
    mnSigmaE_SMD[t][0]->Write();
    mnSigmaE_SMD[t][1]->Write();
    mnSigmaE_SMD2[t][0]->Write();
    mnSigmaE_SMD2[t][1]->Write();
    mnSigmaE_BEMC[t][0]->Write();
    mnSigmaE_BEMC[t][1]->Write();
    mnSigmaE_TOF[t][0]->Write();
    mnSigmaE_TOF[t][1]->Write();
  }
  fout->Write();
  fout->Close();
  return kStOK;
}

//-----------------------------------------------------------------------------
void StPicoElecPurityMaker::DeclareHistograms() {

  // Define parameters for histo definitions that change often
  Float_t nSigLim = 20+1e-6; // To fix bin undulation
  Float_t betaLow = 0.8;
  Float_t betaHigh = 1.8;
  fout->cd();

  trigType = new TH1F("trigType","1:MB,2:HT1,4:HT2,8:HT3. Sum all triggers in hist",40,0,20);
  for(int tr=0; tr<numTrigs; tr++)
  {
    hNEvents[tr] = new TH1F(Form("hNEvents_%i",tr),"number of events: 0 for total and 2 for MBs events and 4 for Vz<6cm", 10, 0, 10 );
    htriggerindex[tr] =new TH1F(Form("htriggerindex_%i",tr),"triggerindex", 25,0,25);
    mVz_vpd[tr] = new TH1F(Form("vz_vpd_%i",tr),"VZ_VPD distribution (cm)",400,-200,200);
    mVz_tpc[tr] = new TH1F(Form("vz_tpc_%i",tr),"the Vz_TPC distribution (cm) ",400,-200,200);
    mVz_vpdtpc[tr] = new TH2F(Form("Vz_vpdtpc_%i",tr),"VertexZ: VPD VS TPC;TPC;VPD;", 400,-200,200,400,-200,200);
    mdVz[tr] = new TH1F(Form("dVz_%i",tr),"VertexZ 'TPC-VPD' (cm)", 400,-200,200);
    mdVz_tpcVz[tr] = new TH2F(Form("dVz_tpcVz_%i",tr),"VertexZ 'TPC-VPD vs TPC' (cm);VZtpc; VZ(tpc-vpd)",400,-200,200,200,-100,100);
    mVxy[tr] = new TH2F(Form("vxy_%i",tr),"Vertex: Vy Vs Vx distribution (cm);vx;vy; ",200,-10,10,200,-10,10);
    mVRvsVZ[tr] = new TH2F(Form("VRvsVZ_%i",tr),"Vertex: VR vs VZ  (cm);VZ;VR; ", 400, -200, 200,200,0,20);

    mNptracks[tr] = new TH1F(Form("Nptracks_%i",tr),"#primary tracks",200,0,2000);    
    mNgtracks[tr] = new TH1F(Form("Ngtracks_%i",tr),"# global tracks ",200,0,5000);    

    mtrketa_pt[tr] = new TH2F(Form("trketa_pt_%i",tr),"trketa_pt; pt; eta;", 200,0,20,200,-2,2);
    mtrkphi_pt[tr] = new TH2F(Form("trkphi_pt_%i",tr),"trkphi_pt; pt; phi;",200,0,20,200,0,6.3);

    mtrkpt[tr] = new TH1F(Form("trkpt_%i",tr),"the pt distribution of all tracks",200,0,20);
    mtrketa[tr] = new TH1F(Form("trketa_%i",tr),"eta ",200,-2,2);
    mtrkphi[tr] = new TH1F(Form("trkphi_%i",tr),"the phi distribution of all tracks",200,0,6.3);


    mnsigmaPI[tr] = new TH1F(Form("nsigmaPI_%i",tr),"nsigmapion of all tracks",200,-nSigLim,nSigLim);
    mnsigmaK[tr] = new TH1F(Form("nsigmaK_%i",tr),"nsigmaKaon of all tracks",200,-nSigLim,nSigLim);
    mnsigmaE[tr] = new TH1F(Form("nsigmaE_%i",tr),"nsigmaElectron of all tracks",200,-nSigLim,nSigLim);
    mnsigmaP[tr] = new TH1F(Form("nsigmaP_%i",tr),"nsigmaProton of all tracks",200,-nSigLim,nSigLim);

    mdedx_Pt[tr] = new TH2F(Form("Dedx_Pt_%i",tr),"dedx(keV/cm) vs Pt; Pt; Dedx;",200,-20,20,250,0,10);
    hNTracks[tr] = new TH1F(Form("hNTracks_%i",tr),"Number of Tracks before (0) and after cuts (2);",10,0,10);

    //--------------------Eta Dependence Study----------------
    int ndims = 5; //pT, nSigmaE, Eta, Centrality, Vz
    int nbins[5] = {400,200,40,20,400};
    double xmin[5] = {0,-nSigLim,-1,0,-100};
    double xmax[5] = {20,nSigLim,1,20,100};
    mnSigmaE_SMD[tr][0]  = new THnSparseF(Form("nSigmaE_SMD_%i",tr),"nSigmaE nDim Hist",ndims,nbins,xmin,xmax);
    mnSigmaE_SMD2[tr][0] = new THnSparseF(Form("nSigmaE_SMD2_%i",tr),"nSigmaE nDim Hist",ndims,nbins,xmin,xmax);
    mnSigmaE_BEMC[tr][0] = new THnSparseF(Form("nSigmaE_BEMC_%i",tr),"nSigmaE nDim Hist",ndims,nbins,xmin,xmax);
    mnSigmaE_TOF[tr][0]  = new THnSparseF(Form("nSigmaE_TOF_%i",tr),"nSigmaE nDim Hist",ndims,nbins,xmin,xmax);
    mnSigmaE_TOF[tr][1]  = new THnSparseF(Form("nSigmaE_TOF_HFT_%i",tr),"nSigmaE nDim Hist",ndims,nbins,xmin,xmax);
    mnSigmaE_SMD[tr][1]  = new THnSparseF(Form("nSigmaE_SMD_HFT_%i",tr),"nSigmaE nDim Hist",ndims,nbins,xmin,xmax);
    mnSigmaE_SMD2[tr][1] = new THnSparseF(Form("nSigmaE_SMD2_HFT_%i",tr),"nSigmaE nDim Hist",ndims,nbins,xmin,xmax);
    mnSigmaE_BEMC[tr][1] = new THnSparseF(Form("nSigmaE_BEMC_HFT_%i",tr),"nSigmaE nDim Hist",ndims,nbins,xmin,xmax);
    mnSigmaE_SMD[tr][0]->Sumw2();
    mnSigmaE_SMD[tr][1]->Sumw2();
    mnSigmaE_SMD2[tr][0]->Sumw2();
    mnSigmaE_SMD2[tr][1]->Sumw2();
    mnSigmaE_BEMC[tr][0]->Sumw2();
    mnSigmaE_BEMC[tr][1]->Sumw2();
    mnSigmaE_TOF[tr][0]->Sumw2();
    mnSigmaE_TOF[tr][1]->Sumw2();

    // -------- For Centrality Study --------
    gRefMult[tr] = new TH1F(Form("gRefMult_%i",tr),"gRefMult; gRefMult; Counts",100,0,1000);
    gRefMultCor[tr] = new TH1F(Form("gRefMultCor_%i",tr),"gRefMultCor; gRefMult; Counts",100,0,1000);
    gRefMultCorWg[tr] = new TH1F(Form("gRefMultCorWg_%i",tr),"gRefMultCor Reweight; gRefMult; Counts",100,0,1000);
    centrality16[tr] = new TH1F(Form("centrality16_%i",tr),"Centrality; Centrality; Counts",20,0,20);
  }

  // -------- dVz Study -------
  mTPCvsVPD_Vz = new TH2F("TPCvsVPD_Vz","TPC vs VPD Vz (no dVz Cut)",400,-200,200,400,-200,200);
  mTPCvsDVz = new TH2F("TPCvsDVz","TPC vs DVz (no dVz Cut)",400,-200,200,400,-200,200);
}// er chen

//----------------------------------------------------------------------------- 
void StPicoElecPurityMaker::Clear(Option_t *opt) {

}

//----------------------------------------------------------------------------- 
Int_t StPicoElecPurityMaker::Make() {
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
  if(DEBUG)std::cout << "Zach Out: Before Bad Run" << endl;
  //=================zaochen add====================
  if(!event) return kStOK;
  //=================event selection================
  int runId = event->runId();
  for(int i=0;i<mNBadRuns;i++){
    if(runId==mBadRuns[i]) return kStOK;
  }
  if(DEBUG)std::cout << "Zach Out: After Bad Run" << endl;
  trig = 99;
  trigCounter=0.; // Used to find overlaps in samples. Modified after checking event cuts. Re-zero each event.
  if( isMB  ( event ) )   {FillHistograms(0,event);} // Decide what type of trigger you have, use to select what histos to fill
  if( isBHT1( event ) )   {FillHistograms(1,event);}
  if( isBHT2( event ) )   {FillHistograms(2,event);}
  if( isBHT3( event ) )   {FillHistograms(3,event); dVzStudy(event);}
  if(trigCounter == 0) 
    trigCounter = -99;
  trigType->Fill(trigCounter); // Use to detect sample overlaps
  if( trig == 99 ) return kStOK; // if no trigger match, throw out event

  return kStOK;
}

Int_t StPicoElecPurityMaker::FillHistograms(Int_t trig, StPicoEvent* event)
{
  if(DEBUG)std::cout << "Zach Out: In Fill Histograms trig" << trig << endl;
  //=================event selection================
  hNEvents[trig]->Fill(0);
  if(! passEventCuts(event,trig) ) return kStOK;
  if(DEBUG)std::cout << "Zach Out: After Event Cut" << endl;
  hNEvents[trig]->Fill(2);
  if(event->primaryVertex().z() < 6.)
    hNEvents[trig]->Fill(4);

  // For sample overlap
  if(trig == 0) trigCounter += 1;
  if(trig == 1) trigCounter += 2;
  if(trig == 2) trigCounter += 4;
  if(trig == 3) trigCounter += 8;

  //int triggerWord = event->triggerWord();
  //	if(triggerWord>>19 & 0x1
  Int_t triggerWORD=event->triggerWord();
  for(Int_t i=0; i<25; i++){
    if( triggerWORD>>i & 0x1 ) htriggerindex[trig]->Fill(i);

  }

  //---------------event information----------------------------
  if(DEBUG)std::cout << "Zach Out: Get Event Info" << endl;
  Double_t vzvpd=event->vzVpd();
  Double_t vztpc=event->primaryVertex().z();
  Double_t vxtpc=event->primaryVertex().x();
  Double_t vytpc=event->primaryVertex().y();
  Double_t dvz=vztpc-vzvpd;
  Double_t vr=sqrt(vxtpc*vxtpc+vytpc*vytpc);

  // Event Cuts
  if(fabs(vxtpc)<1.0e-5)  return kStOK;
  if(fabs(vytpc)<1.0e-5)  return kStOK;
  if(fabs(vztpc)<1.0e-5)  return kStOK;
  //if(fabs(fzvpd) < vZcut) return kStOK;
  //if(fabs(dvz) < dvZcut)  return kStOK;

  //do the refmult correction
  int mRunId = event->runId();
  int mZDCx = (int)event->ZDCx();
  UShort_t mGRefMult = (UShort_t)(event->grefMult());
  Int_t cent16_grefmult;
  Int_t cent9_grefmult;
  Double_t reweight;
  // NOTE: type should be double or float, not integer
  Double_t grefmultCor;


  if(runMode == AuAu){
    StRefMultCorr* grefmultCorrUtil = CentralityMaker::instance()->getgRefMultCorr() ;
    grefmultCorrUtil->init(mRunId);
    grefmultCorrUtil->initEvent(mGRefMult, vzvpd, mZDCx) ;
    cent16_grefmult = grefmultCorrUtil->getCentralityBin16();
    cent9_grefmult  = grefmultCorrUtil->getCentralityBin9();
    reweight = grefmultCorrUtil->getWeight();
    // NOTE: type should be double or float, not integer
    grefmultCor = grefmultCorrUtil->getRefMultCorr();
  }
  else
  {
    cent16_grefmult = 0;
    cent9_grefmult = 0;
    reweight = 1.;
    grefmultCor = 1.;
  }

  if(fillhistflag){
    centrality16[trig]->Fill(cent16_grefmult);
    gRefMult[trig]->Fill(event->grefMult());
    gRefMultCor[trig]->Fill(grefmultCor);
    gRefMultCorWg[trig]->Fill(grefmultCor,reweight);
    mVz_tpc[trig]->Fill(vztpc);
    mVz_vpd[trig]->Fill(vzvpd);
    mdVz[trig]->Fill(dvz);  
    mVz_vpdtpc[trig]->Fill(vztpc,vzvpd);
    mdVz_tpcVz[trig]->Fill(vztpc,dvz);
    mVxy[trig]->Fill( event->primaryVertex().x(), event->primaryVertex().y() );
    mVRvsVZ[trig]->Fill(vztpc, vr);
  }//


  //Int_t Nptrks = mPicoDst->numberOfTracks();
  Float_t Ranking = event->ranking();
  Float_t zdcx = event->ZDCx();
  Float_t bbcx = event->BBCx();
  zdcx=zdcx/1000.;
  bbcx=bbcx/1000.;
  Int_t Nmtdhits = mPicoDst->numberOfMtdHits();
  Int_t Ntofhits = mPicoDst->numberOfBTofHits();
  Int_t NRefmultPos=event->refMultPos();
  Int_t NRefmultNeg=event->refMultNeg();
  Int_t NGnremult=event->grefMult();
  Int_t NRefmult=event->refMult();
  Int_t NGtrks = event->numberOfGlobalTracks();
  Int_t Ntofmatch = event->nBTOFMatch();
  //	Int_t Nbemchits = event->
  Int_t Nbemcmatch = event->nBEMCMatch();

  //----------track information------------------------  
  Int_t numberoftracks = mPicoDst->numberOfTracks();
  if(fillhistflag){	
    mNptracks[trig]->Fill(numberoftracks);
    mNgtracks[trig]->Fill(NGtrks);
  }//

  //tVzTPC=event->primaryVertex().z()

  StThreeVectorF vertexPos;
  vertexPos = event->primaryVertex();

  Int_t Nprimarytracks=0;
  Int_t ntofmatchcount=0;
  Int_t nmtdmatchcount=0;
  Int_t nbemcmatchcount=0;
  Int_t nhftmatchcount=0;
  Int_t nmuons=0;
  Int_t ntofelecton=0;
  Int_t nbemcelectron=0;
  Float_t particleM[3]={0.938,0.140,0.494};

  // TRACK LOOP
  if(DEBUG)std::cout << "Zach Out: At Track Loop" << endl;
  for(int i=0; i<numberoftracks; i++){

    StPicoTrack* track=(StPicoTrack*) mPicoDst->track(i);

    // Check if pass track quality if fails... skip it
    hNTracks[trig]->Fill(0);
    trkHFTflag = 0;
    Bool_t isGoodTrack = passGoodTrack(event,track,trig);
    Bool_t isGoodTrack_NoEta = passGoodTrack_NoEta(event,track,trig);
    if(!isGoodTrack && !isGoodTrack_NoEta) continue;

    if(isGoodTrack_NoEta){ 
      hNTracks[trig]->Fill(2);

      Double_t meta,mpt,mphi,mcharge,mdedx;

      //change to global track
      meta=track->gMom(event->primaryVertex(),event->bField()).pseudoRapidity();
      if(track->pMom().mag()!=0) Nprimarytracks++;
      mphi=RotatePhi(track->gMom(event->primaryVertex(),event->bField()).phi());
      mpt=track->gMom(event->primaryVertex(),event->bField()).perp();
      mcharge=track->charge();
      mdedx=track->dEdx();

      if(mcharge==0||meta==0||mphi==0||mdedx==0/*||track->pMom().mag()!=0*/) continue; //remove neutral, untracked, or primary tracks
      hNTracks[trig]->Fill(4);
      if(track->isHFTTrack()){
        nhftmatchcount++;
        if(fabs(event->primaryVertex().z()) < vZcutHFT[trig])      
        {
          trkHFTflag = 1; 
          hNTracks[trig]->Fill(6);
        }
      }

      Float_t mmomentum=track->gMom(event->primaryVertex(),event->bField()).mag();
      Double_t nsigpi=track->nSigmaPion();
      Double_t nsigk=track->nSigmaKaon();
      Double_t nsigp=track->nSigmaProton();
      Double_t nsige=track->nSigmaElectron();

      mtrkpt[trig]->Fill( mpt );
      mtrketa[trig]->Fill( meta );
      mtrkphi[trig]->Fill( mphi );
      mnsigmaPI[trig]->Fill( nsigpi );
      mnsigmaP[trig] ->Fill( nsigp );
      mnsigmaK[trig] ->Fill( nsigk );
      mnsigmaE[trig] ->Fill( nsige );

      mtrketa_pt[trig]->Fill(mpt*mcharge,meta);
      mtrkphi_pt[trig]->Fill(mpt*mcharge,mphi);

      mdedx_Pt[trig]->Fill(mpt*mcharge,track->dEdx());

      double sparseFill[5] = {mpt, nsige, meta, cent16_grefmult, event->primaryVertex().z()};

      // BEMC nSig
      if(passBEMCCuts(event, track, trig))
      {
        mnSigmaE_BEMC[trig][0]->Fill(sparseFill, reweight);
        if(trkHFTflag == 1)
        {
          mnSigmaE_BEMC[trig][trkHFTflag]->Fill(sparseFill, reweight);
        } 

        // SMD and BEMC
        int checkSMD = passSMDCuts(event, track, trig);
        if(checkSMD > 0 )// if passes either: 1 = loose cuts or 2 = tight cuts
        {
          mnSigmaE_SMD[trig][0]->Fill(sparseFill, reweight);
          if(trkHFTflag == 1)
          {
            mnSigmaE_SMD[trig][trkHFTflag]->Fill(sparseFill, reweight);
          }
        }
        // Tighter SMD Cuts
        if( checkSMD == 2) // 2 = tight cuts
        {
          mnSigmaE_SMD2[trig][0]->Fill(sparseFill, reweight);
          if(trkHFTflag == 1)
          {
            mnSigmaE_SMD2[trig][trkHFTflag]->Fill(sparseFill, reweight);
          }
        }
      }
      // TOF Information
      if(passTOFCuts(event, track, trig))
      {
        Int_t tofpidid=track->bTofPidTraitsIndex();
        if(tofpidid>0){
          ntofmatchcount++;
          StPicoBTofPidTraits* btofpidtrait=(StPicoBTofPidTraits*) mPicoDst->btofPidTraits(tofpidid);

          Float_t beta=btofpidtrait->btofBeta();
          StPhysicalHelixD helix = track->helix();
          if(beta<1e-4||beta>=(USHRT_MAX-1)/20000){
            Float_t tof = btofpidtrait->btof();
            StThreeVectorF btofHitPos = btofpidtrait->btofHitPos();
            float L = tofPathLength(&vertexPos, &btofHitPos, helix.curvature()); 
            if(tof>0) beta = L/(tof*(c_light/1.0e9));
          }
          Float_t tofbeta = 1./beta;
          //Float_t tofbeta = 1./beta;
          Double_t tofm2=mmomentum*mmomentum*( 1.0/(tofbeta*tofbeta)-1.0);
          //minvsBeta_Pt[trig]->Fill(mpt,tofbeta);
          //if(tofbeta>0){
          // mtofM2_Pt[trig]->Fill(mpt,tofm2);
          //}

          // For Purity
          //mdedxvsBeta    [trig]->Fill(tofbeta, track->dEdx());
          //mnSigmaEvsBeta [trig]->Fill(tofbeta, nsige);
          //mnSigmaPIvsBeta[trig]->Fill(tofbeta, nsigpi);
          //mnSigmaKvsBeta [trig]->Fill(tofbeta, nsigk);
          //mnSigmaPvsBeta [trig]->Fill(tofbeta, nsigp);
          //mtofm2vsBeta   [trig]->Fill(tofbeta, tofm2);

          mnSigmaE_TOF[trig][0]->Fill(sparseFill, reweight);
          if(trkHFTflag == 1)
          {
            mnSigmaE_TOF[trig][trkHFTflag]->Fill(sparseFill, reweight);
          }

          Int_t tofcellid=   btofpidtrait->btofCellId();
          Int_t toftray= (int)tofcellid/192 + 1;
          Int_t tofmodule= (int)((tofcellid%192)/6.)+1;
          Float_t toflocaly = btofpidtrait->btofYLocal();
          Float_t toflocalz = btofpidtrait->btofZLocal();
          // Float_t tofhitPosx = btofpidtrait->btofHitPos().x();
          // Float_t tofhitPosy = btofpidtrait->btofHitPos().y();
          // Float_t tofhitPosz = btofpidtrait->btofHitPos().z();

        }
      }// End TOF
    }//end "isGoodTrack" with eta cut
  }//loop of all tracks
  return kStOK;
}//end of main filling fucntion

// --------------- dVz Study Loop ---------------
Int_t StPicoElecPurityMaker::dVzStudy(StPicoEvent* event){

  if(DEBUG)std::cout << "Zach Out: In dVz Study" << endl;
  Int_t trig = 3; // Only for BHT3
  if(! passEventCuts_NodVz(event,trig) ) return kStOK;
  if(DEBUG)std::cout << "Zach Out: pass noDVz Event Cut" << endl;
  Int_t numberoftracks = mPicoDst->numberOfTracks();
  StThreeVectorF vertexPos;
  vertexPos = event->primaryVertex();


  // TRACK LOOP
  for(int i=0; i<numberoftracks; i++){
    StPicoTrack* track=(StPicoTrack*) mPicoDst->track(i);
    Bool_t isGoodTrack = passGoodTrack(event,track,trig);
    Double_t meta,mpt,mphi,mcharge,mdedx;
    meta=track->gMom(event->primaryVertex(),event->bField()).pseudoRapidity();
    mphi=RotatePhi(track->gMom(event->primaryVertex(),event->bField()).phi());
    mpt=track->gMom(event->primaryVertex(),event->bField()).perp();
    mcharge=track->charge();
    mdedx=track->dEdx();
    Double_t vzvpd = event->vzVpd();
    Double_t vztpc = event->primaryVertex().z();
    Double_t dvz = vzvpd - vztpc;

    if(mcharge==0||meta==0||mphi==0||mdedx==0/*||track->pMom().mag()!=0*/) continue; //remove neutral, untracked, or primary tracks

    // BEMC nSig
    if(passBEMCCuts(event, track, trig))
    {
      int checkSMD = passSMDCuts(event, track, trig);
      if( checkSMD == 2) // 2 = tight cuts
      {
        mTPCvsVPD_Vz->Fill(vztpc,vzvpd);
        mTPCvsDVz -> Fill(vztpc,dvz);
      }
    }
  }
  return kStOK;
}

//=====================ZAOCHEN'S FUNCTION=======================================

Bool_t StPicoElecPurityMaker::checkTriggers(StPicoEvent *ev, int trigType)
{
  for(auto trg = triggers[trigType].begin(); trg < triggers[trigType].end(); ++trg)
  {
    if(ev->isTrigger(*trg))
      return true;
  }
  return false;
}


Bool_t StPicoElecPurityMaker::isBHT0(StPicoEvent *event)
{ 
  return checkTriggers(event,0);
}


Bool_t StPicoElecPurityMaker::isBHT1(StPicoEvent *event)
{ 
  return checkTriggers(event,1);
}

//-----------------------------------------                                              
Bool_t StPicoElecPurityMaker::isBHT2(StPicoEvent *event)
{
  return checkTriggers(event,2);
}

//---------------------------------------------------  
Bool_t StPicoElecPurityMaker::isBHT3(StPicoEvent *event)
{
  return checkTriggers(event,3);
}

Bool_t StPicoElecPurityMaker::isMB(StPicoEvent *event)
{ 
  return checkTriggers(event,4);
}

//------------------------------------------------------------- 

Bool_t StPicoElecPurityMaker::passGoodTrack(StPicoEvent* event, StPicoTrack* track, int trig)
{
  double fithitfrac, chargeq, fhitsdEdx, fhitsFit,feta; 
  double pt = track->gMom(event->primaryVertex(),event->bField()).perp();
  feta=track->gMom(event->primaryVertex(),event->bField()).pseudoRapidity();
  fhitsFit = track->nHitsFit();
  fithitfrac=fhitsFit/track->nHitsMax();
  fhitsdEdx = track->nHitsDedx();
  chargeq=track->charge();
  double PtCut = 0.2;

  double mdca;
  // Get DCA info
  StThreeVectorF vertexPos;
  vertexPos = event->primaryVertex();
  StDcaGeometry *dcaG = new StDcaGeometry();
  dcaG->set(track->params(),track->errMatrix());
  StPhysicalHelixD helix = dcaG->helix();
  delete dcaG;
  StThreeVectorF dcaPoint = helix.at( helix.pathLength(vertexPos.x(), vertexPos.y())  );
  double dcamag= (dcaPoint-vertexPos).mag();
  StThreeVectorF dcaP = helix.momentumAt( vertexPos.x(),vertexPos.y() );
  double dcaXY= ( (dcaPoint-vertexPos).x()*dcaP.y()-(dcaPoint-vertexPos).y()*dcaP.x() )/dcaP.perp();
  double dcaZ= dcaPoint.z() - vertexPos.z();
  mdca = dcamag;

  if(pt> PtCut && fhitsFit >= nhitsFitCut && fhitsdEdx >= nhitsdEdxCut && fithitfrac >= nhitsRatioCut && fabs(chargeq)>0 && fabs(feta) <= etaCut && mdca < dcaCut && mdca > 0.) return true;
  else return false;
}

// ----------------------------------------------------------
Bool_t StPicoElecPurityMaker::passGoodTrack_NoEta(StPicoEvent* event, StPicoTrack* track, int trig)
{
  if(DEBUG)std::cout << "Zach Out: passGoodTrack_NoEta" << endl;
  double fithitfrac, chargeq, fhitsdEdx, fhitsFit,feta; 
  double pt = track->gMom(event->primaryVertex(),event->bField()).perp();
  feta=track->gMom(event->primaryVertex(),event->bField()).pseudoRapidity();
  fhitsFit = track->nHitsFit();
  fithitfrac=fhitsFit/track->nHitsMax();
  fhitsdEdx = track->nHitsDedx();
  chargeq=track->charge();
  double PtCut = 0.2;

  double mdca;
  // Get DCA info
  StThreeVectorF vertexPos;
  vertexPos = event->primaryVertex();
  StDcaGeometry *dcaG = new StDcaGeometry();
  dcaG->set(track->params(),track->errMatrix());
  StPhysicalHelixD helix = dcaG->helix();
  delete dcaG;
  StThreeVectorF dcaPoint = helix.at( helix.pathLength(vertexPos.x(), vertexPos.y())  );
  double dcamag= (dcaPoint-vertexPos).mag();
  StThreeVectorF dcaP = helix.momentumAt( vertexPos.x(),vertexPos.y() );
  double dcaXY= ( (dcaPoint-vertexPos).x()*dcaP.y()-(dcaPoint-vertexPos).y()*dcaP.x() )/dcaP.perp();
  double dcaZ= dcaPoint.z() - vertexPos.z();
  mdca = dcamag;

  if(pt> PtCut && fhitsFit >= nhitsFitCut && fhitsdEdx >= nhitsdEdxCut && fithitfrac >= nhitsRatioCut && fabs(chargeq)>0 && mdca < dcaCut && mdca > 0.) return true;
  else return false;
}

//------------------------------------------------------------- 

Int_t StPicoElecPurityMaker::passSMDCuts(StPicoEvent* event, StPicoTrack* track, int trig)
{
  // Get SMD info

  // Get BEMC info
  Int_t emcpidtraitsid=track->emcPidTraitsIndex();
  double mpoe;
  int dsmadc = 1;
  int bemcId, btowId, nPhi,nEta;
  float zDist, phiDist,e0,adc0;
  if(emcpidtraitsid>=0){
    StPicoEmcPidTraits* emcpidtraits=(StPicoEmcPidTraits*) mPicoDst->emcPidTraits(emcpidtraitsid);

    bemcId = emcpidtraits->bemcId();
    adc0   = emcpidtraits->adc0();
    e0     = emcpidtraits->e0();
    btowId = emcpidtraits->btowId();
    zDist = emcpidtraits->zDist();
    phiDist = emcpidtraits->phiDist();
    nEta = emcpidtraits->nEta();
    nPhi = emcpidtraits->nPhi();
    mpoe = track->gMom(event->primaryVertex(),event->bField()).mag()/emcpidtraits->e();
    // get DSM Adc by finding the tower with same id as trk, then getting that ADC
    int nTrgs = mPicoDst->numberOfEmcTriggers();
    for(int j=0;j<nTrgs;j++){
      StPicoEmcTrigger *trg = (StPicoEmcTrigger*)mPicoDst->emcTrigger(j);
      if((trg->flag() & 0xf)){
        int trgId = trg->id();
        if(btowId == trgId){
          if(DEBUG)cout << "bTowId: " << btowId << " ";
          if(DEBUG)cout << "trgID: " << trgId << " ";
          dsmadc = trg->adc();
          if(DEBUG)cout << "trg->adc(): " << dsmadc << endl;
          break;
        }
      }
    }
  }
  else
  {
    mpoe = 0.0; // if no BEMC, set value = 0
    nPhi = -999;
  }
  double mpt  = track->gMom(event->primaryVertex(),event->bField()).perp();

  if(mpt > bemcPtCut && fabs(zDist) < zDistCut2 && fabs(phiDist) < phiDistCut2 && nEta >= nEtaCut2 && nPhi >= nPhiCut2)
    return 2;
  else if(mpt > bemcPtCut && fabs(zDist) < zDistCut && fabs(phiDist) < phiDistCut && nEta >= nEtaCut && nPhi >= nPhiCut)
    return 1;
  else 
    return 0;
}

Bool_t StPicoElecPurityMaker::passBEMCCuts(StPicoEvent* event, StPicoTrack* track, int trig)
{
  // Get BEMC info
  Int_t emcpidtraitsid=track->emcPidTraitsIndex();
  double mpoe;
  int dsmadc = 1;
  int bemcId, btowId, nPhi,nEta;
  float zDist, phiDist,e0,adc0;
  if(emcpidtraitsid>=0){
    StPicoEmcPidTraits* emcpidtraits=(StPicoEmcPidTraits*) mPicoDst->emcPidTraits(emcpidtraitsid);

    bemcId = emcpidtraits->bemcId();
    adc0   = emcpidtraits->adc0();
    e0     = emcpidtraits->e0();
    btowId = emcpidtraits->btowId();
    zDist = emcpidtraits->zDist();
    phiDist = emcpidtraits->phiDist();
    nEta = emcpidtraits->nEta();
    nPhi = emcpidtraits->nPhi();
    mpoe = track->gMom(event->primaryVertex(),event->bField()).mag()/emcpidtraits->e();
    // Check if hot tower. If so, return BEMC failure
    int runId = event->runId();
    if(checkHotTower(runId,btowId))
      return false;
    // get DSM Adc by finding the tower with same id as trk, then getting that ADC
    int nTrgs = mPicoDst->numberOfEmcTriggers();
    for(int j=0;j<nTrgs;j++){
      StPicoEmcTrigger *trg = (StPicoEmcTrigger*)mPicoDst->emcTrigger(j);
      if((trg->flag() & 0xf)){
        int trgId = trg->id();
        if(btowId == trgId){
          if(DEBUG)cout << "bTowId: " << btowId << " ";
          if(DEBUG)cout << "trgID: " << trgId << " ";
          dsmadc = trg->adc();
          if(DEBUG)cout << "trg->adc(): " << dsmadc << endl;
          break;
        }
      }
    }
  }
  else 
    return false;

  double mpt  = track->gMom(event->primaryVertex(),event->bField()).perp();
  //cout << "pT: " << mpt << " p/E: " << mpoe << " e0: " << e0 << " dsmadc: " << dsmadc << endl;
  if( mpt > bemcPtCut && mpoe > poeCutLow && mpoe < poeCutHigh && dsmadc > getDsmAdcCut(trig) )
    return true;
  else 
    return false;
}

Bool_t StPicoElecPurityMaker::checkHotTower(int runId, int towId)
{
  if(DEBUG)std::cout << "Zach Out: checkHotTower" << endl;
  if(runId >= 15071020 && runId <= 15086063) // range 1 of hot tower ist
  {  
    for(int i=0; i < mNHotTower1; i++)
    {
      if(towId == mHotTower1[i])
        return true;
    }
  }

  if(runId >= 15086064 && runId <= 15128024) // range 2 of hot tower ist
  {  
    for(int i=0; i < mNHotTower2; i++)
    {
      if(towId == mHotTower2[i])
        return true;
    }
  }

  if(runId >= 15128025 && runId <= 15167014) // range 3 of hot tower ist
  {  
    for(int i=0; i < mNHotTower3; i++)
    {
      if(towId == mHotTower3[i])
        return true;
    }
  }
  return false;
}

Bool_t StPicoElecPurityMaker::passTOFCuts(StPicoEvent* event, StPicoTrack* track, int trig)
{
  // Get TOF Infor
  Float_t invBeta=9999;
  Float_t toflocaly=9999;
  Float_t tofMatchFlag = -1;
  StThreeVectorF vertexPos;
  vertexPos = event->primaryVertex();
  Int_t tofpidid=track->bTofPidTraitsIndex();
  if(tofpidid>0){
    StPicoBTofPidTraits* btofpidtrait=(StPicoBTofPidTraits*) mPicoDst->btofPidTraits(tofpidid);

    //------tof information start----------
    //Float_t tofbeta=btofpidtrait->btofBeta();
    Float_t beta=btofpidtrait->btofBeta();
    StPhysicalHelixD helix = track->helix();
    if(beta<1e-4||beta>=(USHRT_MAX-1)/20000){
      Float_t tof = btofpidtrait->btof();
      StThreeVectorF btofHitPos = btofpidtrait->btofHitPos();
      float L = tofPathLength(&vertexPos, &btofHitPos, helix.curvature());
      if(tof>0) beta = L/(tof*(c_light/1.0e9));
    }
    Float_t tofbeta = 1./beta;
    invBeta = (1./tofbeta) - 1.0;
    toflocaly = btofpidtrait->btofYLocal();
    tofMatchFlag = btofpidtrait->btofMatchFlag(); 
  }
  double mpt  = track->gMom(event->primaryVertex(),event->bField()).perp();
  if(mpt < tofPtCut && fabs(invBeta) < tofInvBetaCut && tofMatchFlag > 0 && fabs(toflocaly) < toflocalyCut)
    return true;
  else return false;
}

//-------------------------------------------------------------

Bool_t StPicoElecPurityMaker::Ismuontrack(StPicoEvent* event, StPicoTrack* track)
{
  double pt;
  pt=track->gMom(event->primaryVertex(),event->bField()).perp();
  Int_t mtdpid=track->mtdPidTraitsIndex();
  if(mtdpid<=0)return false;
  StPicoMtdPidTraits* mtdpidtrait=(StPicoMtdPidTraits*) mPicoDst->mtdPidTraits(mtdpid);      
  double mtddz=mtdpidtrait->deltaZ();
  double mtddt=mtdpidtrait->deltaTimeOfFlight();
  if(track->nSigmaPion()<3.0&&track->nSigmaPion()>-1.0&&pt>1.0&&fabs(mtddz)<20.0&&mtddt>-1.0&&mtddt<1.0) return true;
  else return false;

}

//----------------------------------------------------------------
Bool_t StPicoElecPurityMaker::passEventCuts(StPicoEvent* event, int trig)
{
  Double_t vzvpd = event->vzVpd();
  Double_t vztpc = event->primaryVertex().z();
  Double_t dvz = vzvpd - vztpc;
  if(fabs(vztpc) < vZcut[trig] && fabs(dvz) < dvZcut[trig]) return true;
  else return false;
}

//----------------------------------------------------------------
Bool_t StPicoElecPurityMaker::passEventCuts_NodVz(StPicoEvent* event, int trig)
{
  Double_t vzvpd = event->vzVpd();
  Double_t vztpc = event->primaryVertex().z();
  Double_t dvz = vzvpd - vztpc;
  if(fabs(vztpc) < vZcut[trig]) return true;
  else return false;
}

//---------------------------------------------------------------
Double_t StPicoElecPurityMaker::RotatePhi(Double_t phi) const
{
  Double_t outPhi = phi;
  Double_t pi=TMath::Pi();
  while(outPhi<0) outPhi += 2*pi;
  while(outPhi>2*pi) outPhi -= 2*pi;
  return outPhi;
}

//===================================================================

void StPicoElecPurityMaker::SetDefaultCuts()
{ 
  numTrigs = 4;
  setNSigECuts(-1.5,3.0);
  setNSigPCuts(-20,20);
  setNSigKCuts(-20,20);
  setNSigPiCuts(-20,20);
  setvZCuts(0,6.0 ,3.0);  // (vZ, delVz)
  setvZCuts(1,30.0,3.0); // (vZ, delVz)
  setvZCuts(2,30.0,3.0); // (vZ, delVz)
  setvZCuts(3,30.0,30.0); // (vZ, delVz)
  setvZCutsHFT(0,6.0,4.0);  // (vZ, delVz)
  setvZCutsHFT(1,6.0,4.0); // (vZ, delVz)
  setvZCutsHFT(2,6.0,4.0); // (vZ, delVz)
  setvZCutsHFT(3,6.0,30.0); // (vZ, delVz)
  setPrimaryPtCut(3.0, 1.5); // pT < 3 (TOF), pT >1.5 (BEMC)
  setPrimaryEtaCut(1.0); // |eta| < 1.0
  setPrimaryDCACut(1.5); // eDCA < 1.5 cm
  setNhitsCuts(15.,20.,0.52); // nHitsdEdx >= 15, nHitsFit >= 20, nHitsRatio >= 0.52
  setPoECut(0.3, 1.5); // 0.3 < p/E < 1.5
  setToFBetaCut(0.03); // |1/B -1| < 0.03
  setToFLocalyCut(1.8); // |tof_localy| < 1.8
  setKaonEnhCut(0.9,1.1); // 0.9<1/B<1.1 Kaon Species Select
  setPionEnhCut(0.9,1.1); // 0.9<1/B<1.1 Pion Species Select
  setProtonEnhCut(0.9,1.1); // 0.9<1/B<1.1 Proton Species Select
  setDsmAdcCut(0,0); // dsmADC cut sets (not in MB): Use getDsmAdcCut(trig) to return value
  setDsmAdcCut(1,15); // dsmADC cut sets ()
  setDsmAdcCut(2,18); // dsmADC cut sets ()
  setDsmAdcCut(3,25); // dsmADC cut sets ()
  setSMDCuts(0,0,3.,0.8); // nEta>=, nPhi>=, zDist<, phiDist< 
  setSMDCuts2(1,1,3.,0.015); // nEta>=, nPhi>=, zDist<, phiDist< 
}

