#include "histMaker.h"

bool DEBUG = !kTRUE;
bool compareRun12 = kTRUE;
bool compareRun15p = !kTRUE; // If this is true, it will take precedence over Run12 compare

void histMaker(const char* fileName="test", int trg=1)
{
  // Set options at large
  TH1::AddDirectory(kFALSE);
  TH1::SetDefaultSumw2(kTRUE);

  trigSelect = trg; // tell analyze what sample, see getTriggerLabel()
  getTriggerLabel();

  //Get TFile
  TFile* f = new TFile(fileName, "READ");
  TFile* fC;
  TFile* fC2;
  if(compareRun15p){
    fC = new TFile(Form("/star/u/zamiller/PWGSpace/run15ppAnaTree/prod/anaTreeMaker_v3_100516/readTreeOut/out/histHadd/rootfile_temp/1117BHT%i/anaTree15_BHT%i.hists.root",trigSelect, trigSelect),"READ");
    fC2 = new TFile(Form("/star/u/zamiller/PWGSpace/run15ppAnaTree/prod/anaTreeMaker_v3_100516/readTreeOut/out/histHadd/rootfile_temp/1117BHT%i/anaTree15_BHT%i.hists_processed.root",trigSelect, trigSelect),"READ");
  }
  else
    fC= new TFile("/star/u/zamiller/PWGSpace/run15ppAnaTree/prod/anaTreeMaker_v2_080316/Run12RootFile/hist_5_2.root","READ");
  bool inFile = checkIfFileOpen(f);
  bool inFileR = checkIfFileOpen(fC);
  if(!inFile || !inFileR) exit(1);

  getHistograms(f);
  declareHistograms();
  doProjections();
  prepareCanvas();
  prepareLabels();
  drawEHadHists();
  drawEventHists();
  drawInvMassHists();
  drawQAHists();
  drawCutEfficiencyHists();
  getCorrections();
  getdNdpT(fC);
  makeInclDivUnlike(fileName);
  calculateCrossSection(fC,fC2);
  makePDF(fileName);
  writeHistsToOutFile(fileName);
  if(fC->IsOpen()) fC->Close();
  f->Close();
}

void drawEHadHists()
{
  for(int etype=0; etype<3; etype++) 
  {
    for(int ptbin=0; ptbin<numPtBins; ptbin++)
    {
      int activeCanvas = (int) ptbin/9;
      int activeBin = ptbin - activeCanvas*9; 
      drawDeltaPhi(eHadDelPhi[etype][ptbin],dPhiPt[etype][activeCanvas],trigCount[etype][ptbin],activeBin,lbl[ptbin]); 
      if(DEBUG) cout << "pt loop: " << ptbin <<endl;
    } 
    delPhiDelEta->cd(etype+1);
    gPad->SetLogz(1);
    eHadDelPhiDelEta[etype]->Draw("colz");
  }
}

bool checkIfFileOpen(TFile* f)
{
  if(!f->IsOpen())
  { std::cout << "!!! File Not Found !!!" << std::endl;
    return false;}  
  else
    std::cout<< "Opened " << f->GetName() << std::endl;
  return true;
}

void writeHistsToOutFile(const char* fileName)
{
  stringstream ss;
  TString outfile; 
  ss << fileName;
  ss >> outfile;
  outfile.ReplaceAll(".root",Form("_processed.root",trigLabel));
  TFile* oF = new TFile(outfile, "RECREATE");
  bool outFile = checkIfFileOpen(oF);
  if(!outFile){
    cout << "No Outfile made!" << endl;
    return;
  }
  for(int i=0; i<3; i++){
    ePt[i]->Write();
  }
  NPECrossSection->Write();
  npeDivPE->Write();
  oF->Close();
  return;
}  

void getCorrections(){
  // From Run 12
  TFile *infile_Trigger =new TFile("run12Corrections/TrigEfficiency_HT.root","read");
  TFile *infile_Tracking =new TFile("run12Corrections/Tracking_efficiency_HT.root","read");
  TFile *infile_PHE_re =new TFile("run12Corrections/Photonic_re_Efficiency.root","read");
  TFile *infile_purity_ht =new TFile("run12Corrections/purity_HT.root","read");
  TFile *infile_crossSec =new TFile("run12Corrections/run12_Npe_HT.root","read");

  effTrigger[0]=(TH1F *) infile_Trigger->Get("TrigEfficiency_HT0");
  effTrigger[1]=(TH1F *) infile_Trigger->Get("TrigEfficiency_HT2");
  effTracking=(TH1F *) infile_Tracking->Get("Tracking_efficiency_HT");
  effPHEReco=(TH1F *) infile_PHE_re->Get("PHE_re_efficiency");
  purityRun12=(TH1F *) infile_purity_ht->Get("purity_HT");
  xsRun12=(TH1F *) infile_crossSec->Get("run12_HT_NPE_sts");
  xsRun12sys=(TH1F *) infile_crossSec->Get("run12_HT_NPE_sys");

  //Jpsi contribution to e
  TFile* infile_jpsi = new TFile("run15Corrections/Jpsi.root");
  jpsiCorrection = rebinVariableBins((TH1F*)infile_jpsi->Get("JPsi"),numPtBinsEFF2,lowptEFF2,"jpsiCorrection");

  // From Run 15pp
  TFile* infile_purity1 = new TFile("run15Corrections/purityHists_pp15_Eta_BHT1_BEMC_processed.root","read");
  TFile* infile_purity2 = new TFile("run15Corrections/purityHists_pp15_Eta_BHT2_BEMC_processed.root","read");
  TFile *infile_Tracking15 =new TFile("run15Corrections/testSample/Tracking_efficiency_HT.root","read");
  TFile *infile_PHE_re15 =new TFile("run15Corrections/testSample/Photonic_re_Efficiency.root","read");

  effTracking15=(TH1F *) infile_Tracking15->Get("TrackingEff");
  effTrigger15[0]=(TH1F *) infile_Tracking15->Get("TrackingEffRecoEidTrg11");
  effTrigger15[1]=(TH1F *) infile_Tracking15->Get("TrackingEffRecoEidTrg18");
  effPHEReco15=(TH1F *) infile_PHE_re15->Get("PHE_re_efficiency");

  TGraphErrors* purityTG0 = (TGraphErrors*)infile_purity1->Get("drawPurityFit_0_0");
  TGraphErrors* purityTG1 = (TGraphErrors*)infile_purity2->Get("drawPurityFit_0_0");

  TH1F* purityT0 = convertTGraphErrorsToTH1(purityTG0, 100, 0, 20, "PurityT0", "Purity HT1");
  TH1F* purityT1 = convertTGraphErrorsToTH1(purityTG1, 100, 0, 20, "PurityT1", "Purity HT2");

  puritypp[0] = rebinVariableBins(purityT0,numPtBinsEFF2,lowptEFF2,"puritypp_HT1");
  puritypp[1] = rebinVariableBins(purityT1,numPtBinsEFF2,lowptEFF2,"puritypp_HT2");

  // From Run 15 pA
  TFile* infile_purity2 = new TFile("/gpfs/mnt/gpfs01/star/pwg/zamiller/run15pAuAnaTree/prod/anaTree_v2_092816/run15Corrections/purityHists_pA15_Eta_BHT1_SMD_processed.root","read");
  TFile* infile_purity3 = new TFile("/gpfs/mnt/gpfs01/star/pwg/zamiller/run15pAuAnaTree/prod/anaTree_v2_092816/run15Corrections/purityHists_pA15_Eta_BHT2_SMD_processed.root","read");

  TGraphErrors* purityTG2 = (TGraphErrors*)infile_purity2->Get("drawPurityFit_0_0");
  TGraphErrors* purityTG3 = (TGraphErrors*)infile_purity3->Get("drawPurityFit_0_0");

  TH1F* purityT2 = convertTGraphErrorsToTH1(purityTG2, 100, 0, 20, "PurityT2", "Purity HT2");
  TH1F* purityT3 = convertTGraphErrorsToTH1(purityTG3, 100, 0, 20, "PurityT3", "Purity HT3");

  purity[0] = rebinVariableBins(purityT2,numPtBinsEFF2,lowptEFF2,"purity_HT1");
  purity[1] = rebinVariableBins(purityT3,numPtBinsEFF2,lowptEFF2,"purity_HT2");

  effnSigE = rebinVariableBins(efnSigE,numPtBinsEFF2,lowptEFF2,"effnSigE");

  TH1F* temp0[6] = {effTrigger15[trigSelect-1], effTracking15, effnSigE, partECutEfficiency[1], puritypp[trigSelect-1], effPHEReco15};
  TH1F* temp[6];

  char legName[6][100] = {"Trigger Efficiency", "Tracking Efficiency", "nSigmaE Efficiency", "BEMC/SMD eID Efficiency","Electron Purity", "Photonic Rec. Efficiency"};
  TLegend* leg = new TLegend(0.45,0.68,0.88,0.89);

  efficiencies->cd();
  for(int i=0; i<6; i++){
    // Rebin so all have the same binning
    temp[i] = rebinVariableBins(temp0[i],numPtBinsEFF2,lowptEFF2,legName[i]);
    // make the total eff hist
    if(i==0)
      totalEff = (TH1F*)temp[i]->Clone();
    else if(i < 4) // combine tracking, trigger and eID efficiency 
      totalEff->Multiply(temp[i]);

    // draw individual eff
    pretty1DHist(temp[i],colors[i],20+i);
    temp[i]->GetYaxis()->SetRangeUser(0,1.5);
    temp[i]->SetTitle("Correction and Efficiencies;p_{T} (GeV/c); Efficiency");
    temp[i]->Draw((i==0)?"pe":"same pe");
    leg->AddEntry(temp[i],legName[i],"lpe");
    if(i==3)leg->AddEntry(totalEff,"Total Electron Efficiency Correction","pe");
  }
  pretty1DHist(totalEff,colors[5],20);
  totalEff->SetMarkerSize(1.4);
  totalEff->Draw("pe same");
  leg->Draw("same");
  return;
}

TH1F* convertTGraphErrorsToTH1(TGraphErrors* tg, int bins, float xlow, float xhigh, TString name, TString title)
{
  TH1F* h = new TH1F(name, title, bins, xlow, xhigh);
  int numPoints = tg->GetN();
  const int nPoints = numPoints;
  double ax[nPoints], ay[nPoints];
  for(int i=0; i<numPoints; i++)
  {
    tg->GetPoint(i,ax[i],ay[i]);
    float errorY = tg->GetErrorY(i);
    h->SetBinContent(h->GetXaxis()->FindBin(ax[i]),ay[i]);
    h->SetBinError(h->GetXaxis()->FindBin(ax[i]),errorY);
  }
  return h;
}

void calculateCrossSection(TFile* fC, TFile* fC2)
{
  TH1F* ptSpectra[3];
  TH1F* runCompPt[3];
  TString histName[3] = {"inclusivePt","USPt","LSPt"}; 
  float dEta = 1.4; // 0.7 - (-0.7)

  float eqMBXS = 42; // effective MB cross section NSD p+p: 30mb +/- 2.4  (xiaozhi Run12 slides)
  if(trigSelect == 1){
    // use NSD/0.64 for trig bias correction since VPDMB for HT1
    eqMBXS = eqMBXS*0.67; 
  }
  if(trigSelect == 2){
    // use NSD since bbcmb for HT2
    eqMBXS = eqMBXS; 
  }

  int numEvents = vertexZ->Integral();  
  float equivMBEvents = vertexZeqMB->Integral()*3.2; // uses the same event cuts as main sample
  cout << "Equiv MB Events: " << equivMBEvents << " Num Events: " << numEvents << endl;
  TH1F* vertexZ_for_mb = (TH1F*)vertexZ->Clone();
  pretty1DHist(vertexZ_for_mb,kBlack,20);
  pretty1DHist(vertexZeqMB,kRed,20);
  equivMB->cd();
  gPad->SetLogy(1);
  vertexZeqMB->GetYaxis()->SetRangeUser(1,1e10);
  vertexZeqMB->Draw("pe");
  vertexZ_for_mb->Draw("same pe");

  TH1D* runCompHists[3];  
  if(compareRun15p)
    getRun15pHists(fC,runCompHists);
  else
    getRun12Hists(fC,runCompHists);

  for(int i=0; i<3; i++)
  {
    ptSpectra[i] = rebinVariableBins(ePt[i], numPtBinsEFF2, lowptEFF2,histName[i]);
    ptSpectra[i]->SetStats(0);
    ptSpectra[i]->Scale(1.,"width");
    runCompPt[i] = (TH1F*)runCompHists[i]->Clone();
  }

  TH1D* NPECrossSection15p;
  if(compareRun15p) 
    if(fC2->IsOpen())
      NPECrossSection15p = (TH1D*)fC2->Get("NPECrossSection");
    else
    {
      cout << "Failed to Open " << fC2->GetName() << endl 
        << "Proceeding with junk histo!" << endl;
      NPECrossSection15p = new TH1D("broken", "broken", 10, 0, 10);
    }

  TH1F* ptMult = new TH1F("ptMult","ptMult",100,0,20);
  TH1F* axes = new TH1F("axes","axes",100,0,20);
  axes->SetStats(0);
  TH1F* ptMul = rebinVariableBins(ptMult, numPtBinsEFF2, lowptEFF2);
  for(int ii=0;ii<numPtBinsEFF2-1;ii++){
    float ptval = (lowptEFF2[ii+1]+lowptEFF2[ii])/2.;
    ptMul->SetBinContent(ii+1, ptval);
  }

  //////////////////////////////////////////////////////////////////
  ////  NPE = Inclusive*purity - (US-LS)/PhoRecoEff             ////
  ////  XS = 1/(2*pi*dEta*dPt)*1/N_events*1/ElecEff*NPE*eqMBXS  ////
  //////////////////////////////////////////////////////////////////
  ptSpectra[0]->Multiply(puritypp[trigSelect-1]);
  TH1F* inclpur = (TH1F*)ptSpectra[0]->Clone();
  ptSpectra[1]->Add(ptSpectra[2],-1.);
  ptSpectra[1]->Divide(effPHEReco15);
  TH1F* phey = (TH1F*)ptSpectra[1]->Clone();
  ptSpectra[0]->Add(ptSpectra[1],-1.);

  runCompPt[0]->Multiply(purityRun12);
  TH1F* inclpur12 = (TH1F*)runCompPt[0]->Clone();
  runCompPt[1]->Add(runCompPt[2],-1.);
  runCompPt[1]->Divide(effPHEReco);
  TH1F* phey12 = (TH1F*)runCompPt[1]->Clone();
  runCompPt[0]->Add(runCompPt[1],-1.);

  TLegend* legY = new TLegend(0.6,0.6,0.87,0.87);
  yieldComponents->cd(1);
  pretty1DHist(inclpur,kBlack,20);
  pretty1DHist(phey,kRed,21);
  pretty1DHist(inclpur12,kBlue,22);
  pretty1DHist(phey12,kMagenta,23);
  gPad->SetLogy(1);
  inclpur->GetYaxis()->SetRangeUser(1,1e7);
  inclpur->Draw("pe");
  phey->Draw("pe same");
  inclpur12->Draw("pe same");
  phey12->Draw("same");
  legY->AddEntry(inclpur,"Inclusive*Purity","lpe");
  legY->AddEntry(phey,"(US-LS)/#epsilon_{PHE}","lpe");
  legY->AddEntry(inclpur12,"Inclusive*Purity Run 12","lpe");
  legY->AddEntry(phey12,"(US-LS)/#epsilon_{PHE} Run 12","lpe");
  legY->Draw("same");
  yieldComponents->cd(2);
  TH1F* yieldRatio = (TH1F*)inclpur->Clone("yieldRatio");
  TH1F* yieldRatio12 = (TH1F*)inclpur12->Clone("yieldRatio12");
  yieldRatio->SetTitle("Ratio of Incl*pur/(PHE/#epsilon_{PHE});p_{T} (GeV/c);Ratio");
  yieldRatio->Divide(phey);
  yieldRatio12->Divide(phey12);
  yieldRatio->GetYaxis()->SetRangeUser(0,4);
  yieldRatio->Draw("pe");
  yieldRatio12->Draw("same pe");

  TH1F* correctNPEPHE = (TH1F*)yieldRatio->Clone();
  correctNPEPHE->Divide(yieldRatio12);

  TLegend* leg = new TLegend(0.5,0.68,0.87,0.89);
  TString histLab[2] = {"Non-Photonic Electrons","Photonic Electrons"};
  npeYield->cd(1);
  gPad->SetLogy(1);
  axes->GetYaxis()->SetRangeUser(1,1e7);
  axes->DrawClone("pe");
  for(int i=0;i<2;i++){
    axes->SetTitle("NPE Yield;p_{T} (GeV/c); dN/dpT");
    pretty1DHist(ptSpectra[i],colors[i],20+i);
    leg->AddEntry(ptSpectra[i],histLab[i],"lpe");
    ptSpectra[i]->Draw("same pe");
  }
  leg->Draw("same");
  npeYield->cd(2);
  npeDivPE = (TH1F*)ptSpectra[0]->Clone("npeDivPE");
  npeDivPE->SetTitle("NPE/PHE;p_{T} (GeV/c);Ratio NPE/PHE");
  TH1F* npeDivPE12 = (TH1F*)runCompPt[0]->Clone("npeDivPE12");
  npeDivPE->Divide(ptSpectra[1]);
  npeDivPE12->Divide(runCompPt[1]);
  axes->SetTitle("NPE/PHE;p_{T} (GeV/c);Ratio NPE/PHE");
  axes->GetYaxis()->SetRangeUser(0,3);
  axes->DrawClone("pe");
  npeDivPE->Draw("same pe");
  npeDivPE12->Draw("same pe");
  TLegend* leg2 = new TLegend(0.5,0.68,0.87,0.89);
  leg2->AddEntry(npeDivPE,"Run 15","lpe");
  if( compareRun15p )
    leg2->AddEntry(npeDivPE12,"Run 15 pp","lpe");
  else
    leg2->AddEntry(npeDivPE12,"Run 12","lpe");
  leg2->Draw("same");

  NPEYield = (TH1F*)ptSpectra[0]->Clone();
  NPECrossSection = (TH1F*)ptSpectra[0]->Clone("NPECrossSection");
  NPECrossSection->SetTitle("NPE Cross Section;p_{T} (GeV/c);E d^{3}#sigma/dp^{3} (mb GeV^{-2} c^{3})");
  axes->SetTitle("NPE Cross Section;p_{T} (GeV/c);E d^{3}#sigma/dp^{3} (mb GeV^{-2} c^{3})");

  NPECrossSection->Scale(eqMBXS/(2.*TMath::Pi()*dEta*equivMBEvents));
  NPECrossSection->Scale(1./2.); // divided by 2 for e- and e+ included
  // NPECrossSection->Scale(1.,"width"); // ALREADY DONE ABOVE, KEEP FOR REMINDER THAT IT NEEDS TO HAPPEN //divide each bin by its width
  NPECrossSection->Divide(totalEff);  // divide by efficiency for finding electron
  NPECrossSection->Divide(ptMul);
  NPECrossSection->Add(jpsiCorrection,-1); // remove j/psi contrib to e spectra
  axes->SetTitle("NPE Cross Section;p_{T} (GeV/c);E d^{3}#sigma/dp^{3} (mb GeV^{-2} c^{3})");

  /*for(int ii=1;ii<NPECrossSection->GetNbinsX();ii++){
   cout <<  NPECrossSection->GetBinError(ii)/NPECrossSection->GetBinContent(ii) << endl;
  }*/

  // get FONLL
  TFile *File_fonll=new TFile("fonll200gev.root","READ");
  TGraphErrors *gFONLLu=( TGraphErrors *) File_fonll->Get("gFONLLu");
  TGraphErrors *gFONLLc=( TGraphErrors *) File_fonll->Get("gFONLLc");
  TGraphErrors *gFONLLd=( TGraphErrors *) File_fonll->Get("gFONLLl");
  setFONLLStyle(gFONLLc,1,kBlack);
  setFONLLStyle(gFONLLu,2,kBlack);
  setFONLLStyle(gFONLLd,2,kBlack);

  TLegend* leg3 = new TLegend(0.5,0.68,0.87,0.89);
  pretty1DHistLargeMarker(NPECrossSection,colors[1],20);
  leg3->AddEntry(NPECrossSection,"Run 15 p+p","lpe");
  if(compareRun15p)
  {
    pretty1DHistLargeMarker(NPECrossSection15p,colors[2],21);
    leg3->AddEntry(NPECrossSection15p,"Run 15 p+p","lpe");
  }
  else
  {
    pretty1DHistLargeMarker(xsRun12,colors[2],21);
    leg3->AddEntry(xsRun12,"Run 12 p+p","lpe");
  }
  leg3->AddEntry(gFONLLu,"FONLL p+p 200 GeV","l");
  crossSection->cd();
  gPad->SetLogy(1);
  gPad->SetLeftMargin(0.15);
  axes->GetYaxis()->SetTitleOffset(1.5);
  axes->GetYaxis()->SetRangeUser(1e-13,1e-2);
  axes->GetXaxis()->SetRangeUser(0,14);
  axes->DrawClone("pe");
  gFONLLc->Draw("l same");
  gFONLLu->Draw("l same");
  gFONLLd->Draw("l same");
  NPECrossSection->Draw("same pe");
  if(compareRun15p)
    NPECrossSection15p->Draw("same pe");
  else
    xsRun12->Draw("same pe");
  leg3->Draw("same");

  TH1F* Rpp = (TH1F*)NPECrossSection->Clone("Rpp");
  pretty1DHist(Rpp,kBlack,26);
  Rpp->SetTitle("NPE #sigma Ratio;p_{T} (GeV/c);R15/R12 XS");
  float pAInt = Rpp->Integral(Rpp->FindBin(5.5),Rpp->FindBin(13.));
  float ppInt;
  if(compareRun15p)ppInt = NPECrossSection15p->Integral(NPECrossSection15p->FindBin(5.5),NPECrossSection15p->FindBin(13.));
  else ppInt = xsRun12->Integral(xsRun12->FindBin(5.5),xsRun12->FindBin(13.));
  //Rpp->Scale(ppInt/pAInt);
  //RpA->DrawClone("same pe");
  //leg3->AddEntry(RpA,"Run 15 p+p (Normalized High pT)","lpe");
  
  Rpp->Divide(xsRun12);
  Rpp->GetYaxis()->SetRangeUser(0,5);
  rpa->cd();
  Rpp->Draw("pe");
  //rpaLabel->Draw("same");
  return;
}

void drawQAHists()
{
  ptCompare->cd();
  gPad->SetLogz(1);
  hadPtEPt->Draw("colz");

  drawHadQA();
  drawElecQA();
  drawEEQA();
  drawPartEQA();
  drawPartECutEffic();
  drawnSigE();
}

void drawEventHists()
{
  eventHists->cd(1);
  pretty1DHist(vertexZ,kRed,20);
  vertexZ->Draw();
  eventHists->cd(2);
  gPad->SetLogy();
  pretty1DHist(refMult,kRed,20);
  refMult->GetXaxis()->SetRangeUser(0,150);
  refMult->Draw();

  zdcQA->cd(1);
  pretty2DHist(refMultvsZDCx,kRed,24);
  refMultvsZDCx->Draw();
  TF1* refZDCFit = new TF1("refZDCFit","pol2");
  refZDCFit->SetLineColor(kRed);
  refMultvsZDCx->Fit("refZDCFit");
  zdcQA->cd(2);
  pretty2DHist(refMultZDCvsRunIndex,kRed,24);
  refMultZDCvsRunIndex->Draw(); 
}

void drawPartECutEffic()
{
  TH1F* histList[3];
  char titlename[4][100] = {"EMC Matched","+ EMC eID","+ SMD Matched","+ SMD eID"};
  for(int i=0; i < 5; i++)
  {
    if(i==0) for(int l=0;l<3;l++) histList[l] = TPCTracks[l];
    if(i==1) for(int l=0;l<3;l++) histList[l] = EMCMatchedTracks[l];
    if(i==2) for(int l=0;l<3;l++) histList[l] = EMCIdTracks[l];
    if(i==3) for(int l=0;l<3;l++) histList[l] = SMDMatchedTracks[l];
    if(i==4) for(int l=0;l<3;l++) histList[l] = SMDIdTracks[l];
    pretty1DHistFill(histList[0],kRed,1001);
    pretty1DHistFill(histList[2],kAzure+1,1001);
    pretty1DHistFill(histList[1],kGreen+1,1001);
    if(i>=1){
      TLegend* leg = new TLegend(0.4,0.68,0.77,0.89);
      leg->AddEntry(histList[0],"Unlike Sign","f");
      leg->AddEntry(histList[1],"Like Sign","f");
      leg->AddEntry(histList[2],"Unlike-Like","f");
      eIDCutEffic[i-1]->cd(1);
      gPad->SetLogy(1);
      histList[0]->Draw("hist");
      histList[2]->Draw("hist same");
      histList[1]->Draw("hist same");
      leg->Draw("same");
      eIDCutEffic[i-1]->cd(2);
      gPad->SetLogy(1);
      TPCTracks[0]->Draw("hist");
      TPCTracks[2]->Draw("hist same");
      TPCTracks[1]->Draw("hist same");
      eIDCutEffic[i-1]->cd(3);

      TH1F* numerator = (TH1F*)histList[2]->Clone();
      TH1F* denominator = (TH1F*)TPCTracks[2]->Clone();
      TH1F* numRebin = rebinVariableBins(numerator,numPtBinsEFF2,lowptEFF2,"numRebin");
      TH1F* denRebin = rebinVariableBins(denominator,numPtBinsEFF2,lowptEFF2,"denRebin");
      TH1F* div = convertToHist(new TGraphAsymmErrors(numRebin,denRebin,"cl=0.683 b(1,1) mode"),"divh");
      pretty1DHist(div,kRed,20);
      div->GetYaxis()->SetRangeUser(0.,1.3);
      //div->Divide(denRebin);
      div->Draw();
      partECutEfficiency[i-1] = (TH1F*)div->Clone();
      partECutEfficiency[i-1]->SetTitle(titlename[i-1]);
    }
  }
  overlayEfficiencies();
}

void overlayEfficiencies()
{
  efficOverlay->cd();
  TH1F* axes = (TH1F*)partECutEfficiency[0]->Clone();
  axes->SetTitle("eID Cut Efficiency;p_{T} (GeV/c);Efficiency");
  axes->GetYaxis()->SetRangeUser(0.,1.4);
  axes->SetStats(kFALSE);
  axes->Draw();
  TLegend* leg = new TLegend(0.5,0.75,0.88,0.89);
  for(int i=0; i<4; i++)
  {
    pretty1DHist(partECutEfficiency[i],colors[i],20+i);
    partECutEfficiency[i]->Draw("same pe");
    leg->AddEntry(partECutEfficiency[i],partECutEfficiency[i]->GetTitle(),"lpe");
  }
  leg->Draw("same");
}

void drawPartEQA()
{
  TH2F* histList[3];
  for(int i=0; i < 6; i++)
  {
    if(i==0) for(int l=0;l<3;l++) histList[l] = nSigEPartE[l];
    if(i==1) for(int l=0;l<3;l++) histList[l] = pvePartE[l];
    if(i==2) for(int l=0;l<3;l++) histList[l] = nEtaPartE[l];
    if(i==3) for(int l=0;l<3;l++) histList[l] = nPhiPartE[l];
    if(i==4) for(int l=0;l<3;l++) histList[l] = phiDistPartE[l];
    if(i==5) for(int l=0;l<3;l++) histList[l] = zDistPartE[l];
    pElecCuts[i]->cd(1);
    histList[0]->Draw("colz");
    pElecCuts[i]->cd(2);
    histList[1]->Draw("colz");
    pElecCuts[i]->cd(3);
    histList[2]->Draw("colz");
  }
}

void drawEEQA()
{
  TString title[2] = {"Unlike Sign","Like Sign"};
  for(int i=0; i<2; i++)
  {
    pretty1DHist(eeEta[i],kRed,20);
    pretty1DHist(eePhi[i],kRed,20);
    pretty1DHist(eeDca[i],kRed,20);
    setTitleAndAxisLabels(eeEta[i],"","#eta","Counts");
    setTitleAndAxisLabels(eePhi[i],"","#phi (rad)","Counts");
    setTitleAndAxisLabels(eeDca[i],"","DCA (cm)","Counts");
    setTitleAndAxisLabels(eeEtaPhi[i],title[i],"#phi","#eta");
    eeQA[i]->cd(1);
    eeEtaPhi[i]->Draw("colz");
    eeQA[i]->cd(2);
    eePhi[i]->Draw("colz");
    eeQA[i]->cd(3);
    eeEta[i]->Draw("colz");
    eeQA[i]->cd(4);
    eeDca[i]->Draw("colz");
  }

  for(int i=0;i<6;i++)
  {
    eeOriginQA->cd(i+1);
    hPEOrigins[i]->Draw("colz");
  }

  for(int i=0;i<3;i++)
  {
    eePhivMass[i]->GetXaxis()->SetRangeUser(0,1.0);
    eePhivMass[i]->Rebin2D(8,8);
  }
  eeQAPhiv->cd(1);
  eePhivMass[0]->Draw("colz");
  eeQAPhiv->cd(2);
  eePhivMass[1]->Draw("colz");
  eeQAPhiv->cd(3);
  eePhivMass[2]->Draw("colz");
  TH2F* eePhivMassLS = (TH2F*)eePhivMass[1]->Clone();
  eePhivMassLS->SetTitle("Like Sign Combined phiv vs m_{ee}");
  eePhivMassLS->Add(eePhivMass[2]);
  eeQAPhiv->cd(4);
  eePhivMassLS->Draw("colz");
}

void drawElecQA()
{
  pretty1DHist(elecEta,kRed,20);
  pretty1DHist(elecPhi,kRed,20);
  pretty1DHist(elecDca,kRed,20);
  setTitleAndAxisLabels(elecEta,"","#eta","Counts");
  setTitleAndAxisLabels(elecPhi,"","#phi (rad)","Counts");
  setTitleAndAxisLabels(elecDca,"","DCA (cm)","Counts");
  setTitleAndAxisLabels(elecEtaPhi,"Semi-Inclusive Electrons","#phi","#eta");
  elecQA->cd(1);
  gPad->SetLogz(0);
  elecEtaPhi->Draw("colz");
  elecQA->cd(2);
  elecPhi->Draw("pe");
  elecQA->cd(3);
  elecEta->Draw("pe");
  elecQA->cd(4);
  elecDca->Draw("pe");
}

void drawHadQA()
{
  pretty1DHist(hadEta,kRed,20);
  pretty1DHist(hadPhi,kRed,20);
  pretty1DHist(hadDca,kRed,20);
  hadQA->cd(1);
  gPad->SetLogz(0);
  hadEtaPhi->Draw("colz");
  hadQA->cd(2);
  hadPhi->Draw("pe");
  hadQA->cd(3);
  hadEta->Draw("pe");
  hadQA->cd(4);
  hadDca->Draw("pe");
}

void prettyTGraph(TGraphErrors* h, int col, int style, float yrangeLow, float yrangeHigh)
{
  h->GetYaxis()->SetRangeUser(yrangeLow,yrangeHigh);
  h->SetMarkerColor(col);
  h->SetLineColor(col);
  h->SetMarkerStyle(style);
  h->SetMarkerSize(1.2);
}

void prettyTGraph(TGraph* h, int col, int style, float yrangeLow, float yrangeHigh)
{
  h->GetYaxis()->SetRangeUser(yrangeLow,yrangeHigh);
  h->SetMarkerColor(col);
  h->SetLineColor(col);
  h->SetMarkerStyle(style);
  h->SetMarkerSize(1.2);
}

void pretty1DHist(TH1* h, int col, int style)
{
  h->SetMarkerColor(col);
  h->SetLineColor(col);
  h->SetMarkerStyle(style);
  h->SetMarkerSize(1.2);
  h->SetTitleOffset(1.2,"y");
}

void pretty1DHistLargeMarker(TH1* h, int col, int style)
{
  h->SetMarkerColor(col);
  h->SetLineColor(col);
  h->SetMarkerStyle(style);
  h->SetMarkerSize(1.6);
  h->SetTitleOffset(1.2,"y");
}

void pretty2DHist(TH2* h, int col, int style)
{
  //h->SetMarkerColor(col);
  h->SetMarkerStyle(style);
  h->SetMarkerSize(1.2);
  h->SetTitleOffset(1.2,"y");
}

void pretty1DHistFill(TH1* h, int col, int style)
{
  h->SetLineColor(kBlack);
  h->SetMarkerStyle();
  h->SetFillStyle(style);
  h->SetFillColor(col);
  h->SetTitleOffset(1.2,"y");
}

void prepareLabels()
{
  getTriggerLabel();
  char textLabel[100];
  int numEvents = vertexZ->Integral();  
  float nEvents = (float)numEvents/1e6;
  sampleLabel = new TPaveText(.11,.83,.5,.89,"NB NDC");
  sprintf(textLabel,"Run 15, 200 GeV p+p Collisions");
  sampleLabel->AddText(textLabel);
  sprintf(textLabel,"%.2fM %s Events",nEvents,trigLabel);
  sampleLabel->AddText(textLabel);
  sampleLabel->SetFillColorAlpha(kWhite,0);
  for(int ptbin=0; ptbin<numPtBins; ptbin++)
  {
    lbl[ptbin] = new TPaveText(.12,.83,.45,.88,Form("NB NDC%i",ptbin));
    sprintf(textLabel,"%.2f < P_{T,e} < %.2f",lowpt[ptbin],highpt[ptbin]);
    lbl[ptbin]->AddText(textLabel);
    lbl[ptbin]->SetFillColor(kWhite);
  }
  rpaLabel = new TPaveText(.11,.76,.7,.89,"NB NDC");
  sprintf(textLabel,Form("Run 15 #sqrt{s} = 200 GeV; %s Trigger", trigLabel));
  rpaLabel->AddText(textLabel);
  sprintf(textLabel,"p+p Collisions");
  rpaLabel->AddText(textLabel);
  rpaLabel->SetFillColorAlpha(kWhite,0.3);
  
  fitLabel = new TPaveText(.62,.75,.84,.87,"NB NDC");
  fitLabel->SetFillColor(kWhite);
}

void getTriggerLabel()
{
  if(trigSelect == 0)
    sprintf(trigLabel,"BHT0");
  if(trigSelect == 1)
    sprintf(trigLabel,"BHT1");
  if(trigSelect == 2)
    sprintf(trigLabel,"BHT2");
}

void drawInvMassHists()
{
  TLegend* leg = new TLegend(0.35,0.65,0.75,0.85);
  TString legLab[3] = {"Unlike Sign","Like Sign","US - LS"};
  TString legLab2[3] = {"Unlike Sign","Like Sign Pos","Like Sign Neg"};
  TH2F* tempMassSpec = (TH2F*)eeInvMassPt[0]->Clone();
  for(int pairType=0; pairType<3; pairType++) //0 = UnlikeSign, 1 = LikeSign ++&-- Combined, 3 = US - LS
  {
    invMass->cd();
    gPad->SetLogy(0);
    pretty1DHist(eeInvMassAll[pairType],colors[pairType],20+pairType);
    eeInvMassAll[pairType]->GetXaxis()->SetRangeUser(0.,0.25);
    eeInvMassAll[pairType]->Draw((pairType==0) ? "pe" : "same pe");
    leg->AddEntry(eeInvMassAll[pairType],legLab[pairType],"lpe"); 
    if(pairType==2) leg->Draw("same");

    for(int ptbin=0; ptbin<numPtBins; ptbin++)
    {
      int activeCanvas = (int) ptbin/9;
      int activeBin = ptbin - activeCanvas*9; 
      invMassPt[activeCanvas]->cd(activeBin+1);
      gPad->SetLogy(0);
      pretty1DHist(eeInvMass[pairType][ptbin],colors[pairType],20+pairType);
      eeInvMass[pairType][ptbin]->GetXaxis()->SetRangeUser(0.,0.25);
      eeInvMass[pairType][ptbin]->Draw((pairType==0) ? "pe" : "same pe");
      lbl[ptbin]->Draw("same");
      if(pairType==2 && ptbin == 0) leg->Draw("same");
    }
    invMassVsPt->cd(pairType+1);
    gPad->SetLogz();
    eeInvMassPt[pairType]->SetTitle(legLab2[pairType]);
    eeInvMassPt[pairType]->GetYaxis()->SetRangeUser(0,0.5);
    eeInvMassPt[pairType]->Draw("colz");
  }
  invMassVsPt->cd(4);
  gPad->SetLogy(0);
  TH1D* massSpectrum = tempMassSpec->ProjectionY();
  setTitleAndAxisLabels(massSpectrum,"Unlike Sign Mass Spectrum","M_{ee} (GeV/c^{2})","Counts");
  pretty1DHist(massSpectrum,kRed,20);
  massSpectrum->Draw();
}

void drawnSigE()
{
  Double_t pT[numPtBins], pTErr[numPtBins], meanErr[numPtBins], sigmaErr[numPtBins], mean[numPtBins], sigma[numPtBins];
  TLegend* leg = new TLegend(0.55,0.75,0.98,0.92);
  TString legLab[3] = {"Unlike Sign","Like Sign","US - LS"};
  //TF1 *gaus = new TF1("Gaus","[0]*exp(-0.5*pow((x-[1])/[2],2))/sqrt(2.*TMath::Pi())/[2]",-3.5,3.5);
  TF1 *gaus = new TF1("gausn","gausn",-3.5,3.5);
  for(int pairType=0; pairType<3; pairType++) //0 = UnlikeSign, 1 = LikeSign ++&-- Combined, 2 = US - LS
  {
    leg->AddEntry(nSigE[pairType][0],legLab[pairType],"lpe"); 
    for(int ptbin=0; ptbin<numPtBins; ptbin++)
    {
      int activeCanvas = (int) ptbin/9;
      int activeBin = ptbin - activeCanvas*9; 
      nSigEPt[activeCanvas]->cd(activeBin+1);
      gPad->SetLogy(0);
      pretty1DHist(nSigE[pairType][ptbin],colors[pairType],20+pairType);
      nSigE[pairType][ptbin]->GetXaxis()->SetRangeUser(-4.,4.);
      nSigE[pairType][ptbin]->Rebin(3);
      nSigE[pairType][ptbin]->Draw((pairType==0) ? "pe" : "same pe");
      lbl[ptbin]->Draw("same");
      if(pairType==2 && ptbin == 0) leg->Draw("same");
      if(pairType==2) 
      {
        nSigE[pairType][ptbin]->Fit("gausn","IN0R","",-3.0,3.0);
        double *par= gaus->GetParameters();
        gaus->SetParameter(0,par[0]);
        gaus->SetParameter(1,par[1]);
        gaus->SetParameter(2,par[2]);
        nSigE[pairType][ptbin]->Fit("gausn","R","",-3.0,3.0);

        //Taken from Shenghui run14 analysis
        TVirtualFitter * fitter = TVirtualFitter::GetFitter();
        if(fitter == NULL)
        {
          cout << "No fitter found for TVirtualFitter" << endl;
          return;
        }
        double * cov = fitter->GetCovarianceMatrix();
        if(cov == 0)
        {
          for(int iii=0; iii<9; iii++)
            cov[iii] = 0.;
        }
        cout << "Covariance Matrix Found: " << endl;
        cout << cov[0] << " " << cov[1] << " " << cov[2] << endl;
        cout << cov[3] << " " << cov[4] << " " << cov[5] << endl;
        cout << cov[6] << " " << cov[7] << " " << cov[8] << endl << endl;
        meanSigE[ptbin]=gaus->GetParameter(1);
        meanerrSigE[ptbin]=cov[4];
        sigmaSigE[ptbin]=gaus->GetParameter(2);
        sigmaerrSigE[ptbin]=cov[8];
        corSigE[ptbin]=cov[7];
        // For TGraph
        pT[ptbin] = (highpt[ptbin]+lowpt[ptbin])/2.;
        pTErr[ptbin] = (highpt[ptbin]-lowpt[ptbin])/2.;
        mean[ptbin]=gaus->GetParameter(1);
        meanErr[ptbin]=gaus->GetParError(1);
        sigma[ptbin]=gaus->GetParameter(2);
        sigmaErr[ptbin]=gaus->GetParError(2);
      }
    }
  }
  drawnSigMeanSig(pT,pTErr,mean,meanErr,sigma,sigmaErr);
  getnSigEeff();
}

void drawnSigMeanSig(const double* pT, const double* pTErr, const double* mean, const double* meanErr, const double* sigma, const double* sigmaErr)
{
  TGraphErrors* mn = new TGraphErrors(numPtBins,pT,mean,pTErr,meanErr);
  mn->SetName("nSigmaMeanVsPt");
  TGraphErrors* sg = new TGraphErrors(numPtBins,pT,sigma,pTErr,sigmaErr);
  sg->SetName("nSigmaSigmaVsPt");
  prettyTGraph(mn,colors[0],20,-1.,3.);
  prettyTGraph(sg,colors[1],21,-1.,3.);
  TLegend* leg = new TLegend(.12,.75,.45,.87);
  leg->AddEntry(mn,"Mean","lpe");
  leg->AddEntry(sg,"Sigma","lpe");
  sg->SetTitle("n#sigma_{e} Fit Values;P_{T} (GeV/c);");
  nSigMeanSig->cd();
  sg->Draw("APE");
  mn->Draw("SAME PE");
  TF1* fSg = new TF1("fSg","[0]",2.5,12);
  fSg->SetLineColor(kRed);
  fSg->SetParameter(0,0.95);
  TF1* fMn = new TF1("fMn","[0]",2.5,12);
  fMn->SetParameter(0,-0.2);
  sg->Fit(fSg,"R");
  mn->Fit(fMn,"R");
  float sgAv = fSg->GetParameter(0); 
  float mnAv = fMn->GetParameter(0); 
  float sgSg = fSg->GetParError(0); 
  float mnSg = fMn->GetParError(0); 
  fitLabel->AddText(Form("Mean: %f #pm %f",mnAv,mnSg));
  fitLabel->AddText(Form("Sigma: %f #pm %f",sgAv,sgSg));
  fitLabel->Draw("SAME");
  leg->Draw("SAME");
}

void getnSigEeff()
{
  Double_t pT[numPtBins], pTErr[numPtBins], cutEff[numPtBins], cutEffErr[numPtBins];
  for(int ptbin=0; ptbin<numPtBins; ptbin++)
  {
    double meant=meanSigE[ptbin];
    double meanerrt=sqrt(meanerrSigE[ptbin]);
    double sigmat=sigmaSigE[ptbin];
    double sigmaerrt=sqrt(sigmaerrSigE[ptbin]);
    double cort=corSigE[ptbin]/meanerrt/sigmaerrt;
    cout << ptbin << ": " << meant << " " << meanerrt << " " << sigmat << " " << sigmaerrt << " " << cort << endl;
    twogaus[ptbin] = new TF2(Form("twogaus_%i",ptbin),"1./2./3.14/[0]/[1]/sqrt(1-[2]*[2])*exp(-1./2./(1-[2]*[2])*(pow((x[0]-[3])/[0],2)-2*[2]*(x[0]-[3])*(x[1]-[4])/[0]/[1]+pow((x[1]-[4])/[1],2)))",meant-3*meanerrt,meant+3*meanerrt,sigmat-3*sigmaerrt,sigmat+3*sigmaerrt);
    twogaus[ptbin]->SetParameter(0,meanerrt);
    twogaus[ptbin]->SetParameter(1,sigmaerrt);
    twogaus[ptbin]->SetParameter(2,cort);
    twogaus[ptbin]->SetParameter(3,meant);
    twogaus[ptbin]->SetParameter(4,sigmat);

    TF1 *Gaus=new TF1("Gaus","gausn(0)",-4,4);
    for(int j=0;j<15000;j++){
      if(j%5000==0) cout << "begin " << j << "th entry...." << endl;
      double mean1,sigma1;
      twogaus[ptbin]->GetRandom2(mean1,sigma1);
      Gaus->SetParameter(0,1.0);
      Gaus->SetParameter(1,mean1);
      Gaus->SetParameter(2,sigma1);
      nSigEeff[ptbin]->Fill((1.0*Gaus->Integral(-1.,3.))/(1.0*Gaus->Integral(-3.5,3.5)));
    }
    int activeCanvas = (int) ptbin/9;
    int activeBin = ptbin - activeCanvas*9; 
    TF1 *gaus2=new TF1("gaus2","gausn",0,1);
    gaus2->SetLineColor(kBlack);
    nSigEff[activeCanvas]->cd(activeBin+1);
    pretty1DHist(nSigEeff[ptbin],kRed,20);
    nSigEeff[ptbin]->Draw("pe");
    nSigEeff[ptbin]->Fit("gaus2");
    lbl[ptbin]->Draw("same");
    pT[ptbin] = (highpt[ptbin]+lowpt[ptbin])/2.;
    pTErr[ptbin] = (highpt[ptbin]-lowpt[ptbin])/2.;
    cutEff[ptbin] = nSigEeff[ptbin]->GetMean(1);
    cutEffErr[ptbin] = nSigEeff[ptbin]->GetRMS(1);
    efnSigE->SetBinContent(efnSigE->FindBin(pT[ptbin]),cutEff[ptbin]);
    efnSigE->SetBinError(efnSigE->FindBin(pT[ptbin]),cutEffErr[ptbin]);

    twoGaus[activeCanvas]->cd(activeBin+1);
    twogaus[ptbin]->SetTitle("2D Gaussian of #mu and #sigma for nSigE Cut");
    twogaus[ptbin]->GetXaxis()->SetTitle("#mu");
    twogaus[ptbin]->GetYaxis()->SetTitle("#sigma");
    twogaus[ptbin]->Draw("SURF2");
    lbl[ptbin]->Draw("same");
  }
  TGraphErrors* grnSigCut = new TGraphErrors(numPtBins,pT,cutEff,pTErr,cutEffErr);
  grnSigCut->SetName("nSigmaEfficVsPt");
  prettyTGraph(grnSigCut,kRed,20,0.,1.2);
  nSigCutPlot->cd();
  grnSigCut->SetTitle("n#sigma_{E} Cut Efficiency;P_{T} (GeV/c);Efficiency");
  grnSigCut->Draw("APE");
  grnSigCut->Fit("pol0","R","",2,12.);
}

void declareHistograms(){

  for(int ptbin=0; ptbin<numPtBins; ptbin++)
  {
    nSigEeff[ptbin] = new TH1F(Form("nSigEeff_%i",ptbin),"nSigE Cut Eff",1000,0,1);
  }
  efnSigE = new TH1F("efnSigE","nSigE Eff. vs p_{T};p_{T};nSigE Cut Eff",1000,0,20);

}

void doProjections()
{
  char name[3][100] = {"Inclusive","Unlike Sign Pair", "Like Sign Pair"}; 
  for(int etype=0; etype<3; etype++) // 0 = inclusive, 1 = Unlike Sign Pairs, 2 = Like Sign Pairs (for delphi), 0= +-, 1= ++, 2= -- for invMass 
  {
    for(int ptbin=0; ptbin<numPtBins; ptbin++)
    {
      eHadDelPhi[etype][ptbin] = (TH1D*)eHadDelPhiPt[etype]->ProjectionY(Form("eHadDelPhi_%i_%i",etype,ptbin),eHadDelPhiPt[etype]->GetXaxis()->FindBin(lowpt[ptbin]),eHadDelPhiPt[etype]->GetXaxis()->FindBin(highpt[ptbin])-1);
      eHadDelPhi[etype][ptbin]->SetTitle((ptbin==0) ? Form("%s Electrons",name[etype]) : "");
      trigCount[etype][ptbin] = ePt[etype]->Integral(ePt[etype]->GetXaxis()->FindBin(lowpt[ptbin]),ePt[etype]->GetXaxis()->FindBin(highpt[ptbin])-1);
      eeInvMass[etype][ptbin] = (TH1D*)eeInvMassPt[etype]->ProjectionY(Form("eeInvMass_%i_%i",etype,ptbin),eeInvMassPt[etype]->GetXaxis()->FindBin(lowpt[ptbin]),eeInvMassPt[etype]->GetXaxis()->FindBin(highpt[ptbin])-1); 

      if(etype==2)
      {
        eeInvMass[1][ptbin]->Add(eeInvMass[2][ptbin]); // combine ++ and -- to form LS
        eeInvMass[2][ptbin]->Add(eeInvMass[0][ptbin],eeInvMass[1][ptbin],1.,-1.);
      }
    }

    // No pT dependence projections
    eeInvMassAll[etype] = (TH1D*)eeInvMassPt[etype]->ProjectionY(Form("eeInvMassAll_%i",etype)); 
    eeEta[etype] = (TH1D*)eeEtaPhi[etype]->ProjectionY(Form("eeEta_%i",etype));
    eePhi[etype] = (TH1D*)eeEtaPhi[etype]->ProjectionX(Form("eePhi_%i",etype));
    eeDca[etype] = (TH1D*)eeDcaPt[etype] ->ProjectionY(Form("eeDca_%i",etype));
    if(etype==2)
    {
      eeInvMassAll[1]->Add(eeInvMassAll[2]); // combine ++ and -- to form LS
      eeInvMassAll[2]->Add(eeInvMassAll[0],eeInvMassAll[1],1.,-1.);

      eeEta[1]->Add(eeEta[2]);
      eePhi[1]->Add(eePhi[2]);
      eeDca[1]->Add(eeDca[2]);
      eeEtaPhi[1]->Add(eeEtaPhi[2]);
      eeEtaPhi[1]->Add(eeEtaPhi[2]);
    }

  }

  for(int etype=0; etype<2; etype++)
  {
    for(int ptbin=0; ptbin<numPtBins; ptbin++)
    {
      nSigE[etype][ptbin] = (TH1D*)nSigEPartE[etype]->ProjectionY(Form("nSigE_%i_%i",etype,ptbin),nSigEPartE[etype]->GetXaxis()->FindBin(lowpt[ptbin]),nSigEPartE[etype]->GetXaxis()->FindBin(highpt[ptbin])-1); 
      if(etype==1)
      {
        nSigE[2][ptbin] = (TH1D*)nSigE[0][ptbin]->Clone();
        nSigE[2][ptbin]->Add(nSigE[1][ptbin],-1.);
      }
    }
  }

  hadEta = (TH1D*)hadEtaPhi->ProjectionX();
  setTitleAndAxisLabels(hadEta,"Hadron Eta","#eta","Counts");
  hadPhi = (TH1D*)hadEtaPhi->ProjectionY();
  setTitleAndAxisLabels(hadPhi,"Hadron Phi","#phi (rad)","Counts");

  elecEta = (TH1D*)elecEtaPhi->ProjectionY();
  setTitleAndAxisLabels(elecEta,"Semi-Inclusive Electron Eta","#eta","Counts");
  elecPhi = (TH1D*)elecEtaPhi->ProjectionX();
  setTitleAndAxisLabels(elecPhi,"Semi-Inclusive Electron Phi","#phi (rad)","Counts");
  elecDca = (TH1D*)elecDcaPt->ProjectionY();
  setTitleAndAxisLabels(elecPhi,"Semi-Inclusive Electron Phi","DCA (cm)","Counts");
}

void setTitleAndAxisLabels(TH1* h, TString title, TString xLbl, TString yLbl)
{
  h->SetTitle(title);
  h->GetXaxis()->SetTitle(xLbl);
  h->GetYaxis()->SetTitle(yLbl);
}

void setTitleAndAxisLabels(TH2* h, TString title, TString xLbl, TString yLbl)
{
  h->SetTitle(title);
  h->GetXaxis()->SetTitle(xLbl);
  h->GetYaxis()->SetTitle(yLbl);
}

void prepareCanvas()
{
  for(int etype=0; etype<3; etype++) // 0 = inclusive, 1 = Unlike Sign Pairs, 2 = Like Sign Pairs
  {
    //create canvas
    for(int q=0; q<numCanvas; q++)
    {
      dPhiPt[etype][q] = new TCanvas(Form("dPhiPt_%i_%i",etype,q),"pT Dependence of DelPhi",50,50,1050,1050);
      dPhiPt[etype][q]->Divide(3,3);
    } 
  }
  delPhiDelEta = new TCanvas("delPhiDelEta","Electron-Hadron Correlations",50,50,1050,1050);
  delPhiDelEta->Divide(2,2);

  for(int q=0; q<numCanvas; q++)
  {
    invMassPt[q] = new TCanvas(Form("invMassPt_%i",q),"Invariant Mass Spectrum",50,50,1050,1050);
    invMassPt[q]->Divide(3,3);
    nSigEPt[q] = new TCanvas(Form("nSigEPt_%i",q),"nSigmaE Spectrum",50,50,1050,1050);
    nSigEPt[q]->Divide(3,3);
    nSigEff[q] = new TCanvas(Form("nSigEff_%i",q),"nSigmaE Cut Efficiency",50,50,1050,1050);
    nSigEff[q]->Divide(3,3);
    twoGaus[q] = new TCanvas(Form("twoGaus_%i",q),"Covariance Double Gaussian",50,50,1050,1050);
    twoGaus[q]->Divide(3,3);

  }
  invMass     = new TCanvas("invMass","Invariant Mass All pT",50,50,1050,1050);
  invMassVsPt = new TCanvas("invMassVsPt","Invariant Mass Vs pT",50,50,1050,1050);
  invMassVsPt->Divide(2,2);
  ptCompare     = new TCanvas("ptCompare","pT Comparison",50,50,1050,1050);
  equivMB     = new TCanvas("equivMB","Equivalent MB Events",50,50,1050,1050);
  cutEfficiency = new TCanvas("cutEfficiency","Cut Efficiency",50,50,1050,1050);
  cutEfficiency -> Divide(2,2);
  hadQA = new TCanvas("hadQA","Hadron Based QA",50,50,1050,1050);
  hadQA -> Divide(2,2);
  elecQA = new TCanvas("elecQA","Electron Based QA",50,50,1050,1050);
  elecQA -> Divide(2,2);
  for(int i=0;i<2;i++){
    eeQA[i] = new TCanvas(Form("eeQA_%i",i),"Elec Pair Based QA",50,50,1050,1050);
    eeQA[i] -> Divide(2,2);
  }
  eeQAPhiv = new TCanvas("eeQAPhiv","Elec Pair Based QA",50,50,1050,1050);
  eeQAPhiv -> Divide(2,2);
  eeOriginQA = new TCanvas("eeOriginQA","Elec Pair Based QA",50,50,1050,1050);
  eeOriginQA -> Divide(3,3);

  for(int i=0; i<6; i++)
  {
    pElecCuts[i] = new TCanvas(Form("pElecCuts_%i",i),"Partner Electron Based QA",50,50,1050,1050);
    pElecCuts[i] -> Divide(2,2);
  }
  for(int i=0; i<4; i++)
  {
    eIDCutEffic[i] = new TCanvas(Form("eIDCutEffic_%i",i),"Partner eID Based QA",50,50,1050,1050);
    eIDCutEffic[i] -> Divide(2,2);
  }
  efficOverlay = new TCanvas("efficOverlay","Partner eID Based QA",50,50,1050,1050);
  efficiencies = new TCanvas("efficiencies","Efficiencies For Cross Section",50,50,1050,1050);
  crossSection = new TCanvas("crossSection","NPE Cross Section",50,50,1050,1050);
  rpa = new TCanvas("rpa","NPE R_{pA}",50,50,1050,1050);
  npeYield = new TCanvas("npeYield","NPE Yield",50,50,1050,1050);
  npeYield->Divide(1,2);
  nSigCutPlot = new TCanvas("nSigCutPlot","nSigmaE Cut Efficiency",50,50,1050,1050);
  nSigMeanSig = new TCanvas("nSigMeanSig","nSigmaE Fit Results",50,50,1050,1050);
  dndpt = new TCanvas("dndpt","dndpt",50,50,1050,1050);
  Run12Compare = new TCanvas("Run12Compare","Run12Compare",50,50,1050,1050);
  Run12Compare->Divide(2,2);
  Run12Ratio = new TCanvas("Run12Ratio","Run12Ratio",50,50,1050,1050);
  Run12Ratio->Divide(2,2);
  inclUSRatio = new TCanvas("inclUSRatio","inclUSRatio",50,50,1050,1050);

  eventHists = new TCanvas("eventHists","Event Level Hists",50,50,1050,1050);
  eventHists->Divide(1,2);
  zdcQA = new TCanvas("zdcQA","ZDC QA Hists",50,50,1050,1050);
  zdcQA->Divide(1,2);
  yieldComponents = new TCanvas("yieldComponent","yieldComponents",50,50,1050,1050);
  yieldComponents->Divide(1,2);
  if(DEBUG) cout << "Canvas made." << endl;
}

void drawCutEfficiencyHists()
{
  int markerStart = 20;
  TString name[4] = {"Track Quality Cuts (TQC)", "TQC + eID", "TQC + HTeID", "TQC+Trigger Cuts"};
  TLegend* leg = new TLegend(0.3,0.65,0.77,0.88);
  cutEfficiency->cd(1);
  gPad->SetLogy();
  for(int l=0; l<4; l++)
  {
    ePt_eff[l]->SetLineColor(colors[l]);
    ePt_eff[l]->SetMarkerColor(colors[l]);
    ePt_eff[l]->SetMarkerStyle(markerStart+l);
    ePt_eff[l]->SetMarkerSize(0.7);
    ePt_eff[l]->Draw((l==0)?"pe":"same pe");
    leg->AddEntry(ePt_eff[l],name[l],"lpe");
  }
  leg->Draw("same");

  for(int pad=2; pad<5; pad++)
  {
    computeEfficiencyAndDraw(pad,ePt_eff[0],ePt_eff[pad-1]);
  }
}

void computeEfficiencyAndDraw(int pad, TH1F* denom, TH1F* numer)
{
  TString title[3] = {"eID Efficiency", "HT eID Efficiency", "Trigger Efficiency"};
  TH1F* num = (TH1F*)numer->Clone();
  num->Divide(denom);
  num->SetTitle(title[pad-2]);
  num->GetYaxis()->SetTitle("Efficiency");
  num->GetYaxis()->SetRangeUser(0.,1.);
  cutEfficiency->cd(pad);
  num->Draw("pe");
}

void drawDeltaPhi(TH1D* h, TCanvas* c, float norm, int activeBin, TPaveText* tpt)
{
  c->cd(activeBin+1);
  h->Rebin(4);
  float binWidth = h->GetXaxis()->GetBinWidth(10);
  h->Scale(1./norm/binWidth);
  h->GetXaxis()->SetRangeUser(-1.5,4.4);
  h->SetMarkerStyle(20);
  h->SetMarkerColor(kBlack);
  h->SetLineColor(kBlack);
  h->Draw();
  tpt->Draw("same");
}

void getHistograms(TFile* f)
{
  //Events
  refMult = (TH1F*)f->Get("hRefMultCut");
  vertexZ = (TH1F*)f->Get("hVertexZCut");
  vertexZMB = (TH1F*)f->Get("hVertexZCut_MB");
  vertexZeqMB = (TH1F*)f->Get("hVertexZCut_eqMB");

  // electrons
  eHadDelPhiPt[0] = (TH2F*)f->Get("hHadEDelPhiPt");
  eHadDelPhiPt[1] = (TH2F*)f->Get("hHadEEDelPhiPt_US");
  eHadDelPhiPt[2] = (TH2F*)f->Get("hHadEEDelPhiPt_LS");
  eHadDelPhiDelEta[0] = (TH2F*)f->Get("hHadEDelPhiDelEta_0");
  eHadDelPhiDelEta[1] = (TH2F*)f->Get("hHadEDelPhiDelEta_1");
  eHadDelPhiDelEta[2] = (TH2F*)f->Get("hHadEDelPhiDelEta_2");
  hadPtEPt = (TH2F*)f->Get("hHadPtEPt");
  ePt[0] = (TH1F*)f->Get("hEPt");
  ePt[1] = (TH1F*)f->Get("hEEPt_US");
  ePt[2] = (TH1F*)f->Get("hEEPt_LS");
  ePt_eff[0] = (TH1F*)f->Get("hEPt_eff_0"); // Track Quality Cuts only
  ePt_eff[1] = (TH1F*)f->Get("hEPt_eff_1"); // Track Quality Cuts, eID cuts
  ePt_eff[2] = (TH1F*)f->Get("hEPt_eff_2"); // Track Quality Cuts, HT eID cuts
  ePt_eff[3] = (TH1F*)f->Get("hEPt_eff_3"); // Track Quality Cuts, ADC0 cut 
  eeInvMassPt[0] = (TH2F*)f->Get("hEENumInvMassvsPtMB"); // Unlike invmass 
  eeInvMassPt[1] = (TH2F*)f->Get("hEEDenInvMassvsPtLikePosMB"); // Like ++ invmass 
  eeInvMassPt[2] = (TH2F*)f->Get("hEEDenInvMassvsPtLikeNegMB"); // Like -- invmass 

  //Electron QA
  elecEtaPhi = (TH2F*)f->Get("hEEtavsPhi");
  elecDcaPt  = (TH2F*)f->Get("hEDcavsPt");

  //Pair Hists QA
  eePhivMass[0] = (TH2F*)f->Get("hUSphivM");
  eePhivMass[1] = (TH2F*)f->Get("hLSPosphivM");
  eePhivMass[2] = (TH2F*)f->Get("hLSNegphivM");
  eeEtaPhi[0]   = (TH2F*)f->Get("hEEUSEtavsPhi");
  eeEtaPhi[1]   = (TH2F*)f->Get("hEELSPosEtavsPhi");
  eeEtaPhi[2]   = (TH2F*)f->Get("hEELSNegEtavsPhi");
  eeDcaPt[0]   = (TH2F*)f->Get("hEEUSPairDcavsPt");
  eeDcaPt[1]   = (TH2F*)f->Get("hEELSPosPairDcavsPt");
  eeDcaPt[2]   = (TH2F*)f->Get("hEELSNegPairDcavsPt");

  hPEOrigins[0] =  (TH2F*)f->Get("hPEUSOyOx");
  hPEOrigins[1] =  (TH2F*)f->Get("hPEUSOxOz");
  hPEOrigins[2] =  (TH2F*)f->Get("hPEUSOrOz");
  hPEOrigins[3] =  (TH2F*)f->Get("hPELSOyOx");
  hPEOrigins[4] =  (TH2F*)f->Get("hPELSOxOz");
  hPEOrigins[5] =  (TH2F*)f->Get("hPELSOrOz");

  // Hadron Hists
  hadEtaPhi = (TH2F*)f->Get("hHadEtaPhi");
  hadDca =   (TH1F*)f->Get("hHadDca");

  //Partner Elec Hists
  for(int i=0; i<2; i++) //0-Unlike, 1-Like
  {
    nSigEPartE[i] =  (TH2F*)f->Get(Form("hNSigEPartElec_%i",i));
    pvePartE[i] =  (TH2F*)f->Get(Form("hPvePartElec_%i",i));
    nEtaPartE[i] =  (TH2F*)f->Get(Form("hnEtaPartElec_%i",i));
    nPhiPartE[i] =  (TH2F*)f->Get(Form("hnPhiPartElec_%i",i));
    phiDistPartE[i] =  (TH2F*)f->Get(Form("hphiDistPartElec_%i",i));
    zDistPartE[i] =  (TH2F*)f->Get(Form("hzDistPartElec_%i",i));

    TPCTracks[i] =  (TH1F*)f->Get(Form("hTPCTracks_%i",i));
    EMCMatchedTracks[i] =  (TH1F*)f->Get(Form("hEMCMatchedTracks_%i",i));
    EMCIdTracks[i] =  (TH1F*)f->Get(Form("hEMCIdTracks_%i",i));
    SMDMatchedTracks[i] =  (TH1F*)f->Get(Form("hSMDMatchedTracks_%i",i));
    SMDIdTracks[i] =  (TH1F*)f->Get(Form("hSMDIdTracks_%i",i));

    if(i == 1)
    {
      makeUnlikeMinusLikePartnerElectrons();
    }
  }

  // ZDC Hists
  refMultZDCvsRunIndex =  (TH2F*)f->Get("hgRefMultZDCvsRunIndex");
  refMultvsZDCx =  (TH2F*)f->Get("hgRefMultvsZDCx");

  if(DEBUG) cout << "Get Hist." << endl;
}

void getdNdpT(TFile* fC)
{
  TString histName[3] = {"inclusivePt","USPt","LSPt"}; 
  TString legName[3] = {"Inclusive Electron", "Unlike Sign Photonic", "Like Sign Photonic"};
  TString legName12[3] = {"Run 12 Inclusive", "Run 12 Unlike", "Run 12 Like"};
  TString legName15[3] = {"Run 15 pp Inclusive", "Run 15 pp Unlike", "Run 15 pp Like"};

  TLegend* leg = new TLegend(.53,.74,.88,.88);
  dndpt->cd();
  gPad->SetLogy(1);
  TH1F* ptDist[3];
  for(int i=0; i<3; i++)
  {
    TH1F* ptSpect = (TH1F*)ePt[i]->Clone();
    ptSpect->SetStats(0);
    pretty1DHist(ptSpect,colors[i],20+i);
    leg->AddEntry(ptSpect,legName[i],"lpe");
    TH1F* ptSpectra = rebinVariableBins(ptSpect, numPtBinsEFF2, lowptEFF2);
    ptSpectra->Scale(1.,"width");
    if(i==0){
      ptSpectra->GetYaxis()->SetRangeUser(1,1e7);
      ptSpectra->SetTitle("Electron Spectra; p_{T} (GeV/c); Raw dN/dpT");
    }
    ptSpectra->Draw((i==0)?"pe":"same pe");
    ptDist[i] = (TH1F*)ptSpectra->Clone();

    // Used for determining the size needed for embedding. Used alongside relative efficiency uncertainty
    // in Run 12 embedding and the embedding size in Run 12
    /* cout << "Relative Uncertainty for " << legName[i] << endl;
    for(int ii=1; ii<ptSpectra->GetNbinsX();ii++){
      cout << ptSpectra->GetBinError(ii)/ptSpectra->GetBinContent(ii) << endl;
    }*/
  }
  sampleLabel->Draw("SAME");

  if(compareRun12 || compareRun15p){
    TH1D* runCompHists[3];
    if(compareRun15p)
      getRun15pHists(fC,runCompHists);
    else
      getRun12Hists(fC,runCompHists);

    for(int ii=0;ii<3;ii++)
    {
      runCompHists[ii]->Draw("pe same");

      if(compareRun15p)
        leg->AddEntry(runCompHists[ii],legName15[ii],"lpe");
      else
        leg->AddEntry(runCompHists[ii],legName12[ii],"lpe");
    }
    leg->Draw("same");
    
    float normLow  = 4.;
    float normHigh = 6.;
    if(trigSelect == 2)
      {
        normLow  = 6.;
        normHigh = 8.;
      }

    for(int i=0; i<3; i++)
    {
      TH1F* r15 = (TH1F*)ptDist[i]->Clone();
      TH1F* r12 = (TH1F*)runCompHists[i]->Clone();
      double r15i = r15->Integral(r15->FindBin(normLow),r15->FindBin(normHigh));
      double r12i = r12->Integral(r12->FindBin(normLow),r12->FindBin(normHigh));
      r12->Scale(r15i/r12i);
      r15->SetTitle("Run Comparison - Scaled to Match Run 15 Yield");

      TPaveText* normLabel = new TPaveText(.11,.83,.5,.89,"NB NDC");
      normLabel->AddText(Form("Normalization Factor: %.3f",r15i/r12i));
      normLabel->SetFillColorAlpha(kWhite,0);

      Run12Compare->cd(i+1);
      gPad->SetLogy(1);
      r15->Draw("pe");
      r12->Draw("same pe");
      TLegend* leg12 = new TLegend(.53,.74,.88,.88);
      leg12->AddEntry(r15, legName[i], "lpe");
      if(compareRun15p)
        leg12->AddEntry(r12,legName15[i],"lpe");
      else
        leg12->AddEntry(r12,legName12[i],"lpe");
      leg12->Draw("same");
      normLabel->Draw("same");
      TH1F* ratio = (TH1F*)r15->Clone();
      if(i==0)ratio->Multiply(puritypp[trigSelect-1]);
      TH1F* r12Pur = (TH1F*)r12->Clone();
      if(i==0)r12Pur->Multiply(purityRun12);
      ratio->Divide(r12Pur);
      if( compareRun15p )
      {
        if(i==0) ratio->SetTitle(Form("Ratio of Run 15 to Run 15 pp %s;p_{T} (GeV/c);R15pA/R15pp (Normalized Yields * Purity)",legName[i].Data()));
        else ratio->SetTitle(Form("Ratio of Run 15 to Run 15 pp %s;p_{T} (GeV/c);R15pA/R15pp Normalized Yields",legName[i].Data()));
      }
      else
      {
        if(i==0) ratio->SetTitle(Form("Ratio of Run 15 to Run 12 %s;p_{T} (GeV/c);R15/R12 (Normalized Yields * Purity)",legName[i].Data()));
        else ratio->SetTitle(Form("Ratio of Run 15 to Run 12 %s;p_{T} (GeV/c);R15/R12 Normalized Yields",legName[i].Data()));
      }
      pretty1DHist(ratio, kRed, 24);
      Run12Ratio->cd(i+1);
      ratio->GetYaxis()->SetRangeUser(0,4);
      ratio->Draw("pe");
    }

    TH1F* run12hist[3];
    for(int ii=0;ii<3;ii++) run12hist[ii] = (TH1F*)runCompHists[ii]->Clone();
    ptDist[0]->Multiply(puritypp[trigSelect-1]);
    ptDist[1]->Add(ptDist[2],-1.);
    ptDist[1]->Divide(effPHEReco15);
    ptDist[0]->Add(ptDist[1],-1.);
    if( compareRun15p )
      run12hist[0]->Multiply(puritypp[trigSelect-1]);
    else
      run12hist[0]->Multiply(purityRun12);
    run12hist[1]->Add(run12hist[2],-1.);
    run12hist[1]->Divide(effPHEReco);
    run12hist[0]->Add(run12hist[1],-1.);
    Run12Compare->cd(4);
    TLegend* leg3 = new TLegend(0.45,0.68,0.88,0.89);
    gPad->SetLogy(1);
    if( compareRun15p )
      ptDist[0]->SetTitle("Run Comparison Run 15 to Run 15 pp NPE;p_{T} (GeV/c);dN/dpT");
    else
      ptDist[0]->SetTitle("Run Comparison Run 15 to Run 12 NPE;p_{T} (GeV/c);dN/dpT");
    ptDist[0]->GetYaxis()->SetRangeUser(1,1e7);
    ptDist[0]->DrawClone("pe");
    run12hist[0]->DrawClone("pe same");
    leg3->AddEntry((TH1F*)ptDist[0]->Clone(),"Run 15 NPE","lpe");
    if( compareRun15p )
      leg3->AddEntry(run12hist[0],"Run 15 p+p NPE","lpe");
    else
      leg3->AddEntry(run12hist[0],"Run 12 NPE","lpe");
    leg3->Draw("same");
    if( compareRun15p )
      ptDist[0]->SetTitle("Ratio of Run 15 to Run 15 pp NPE;p_{T} (GeV/c);R15/R15pp (dN/dpT)");
    else
      ptDist[0]->SetTitle("Ratio of Run 15 to Run 12 NPE;p_{T} (GeV/c);R15/R12 (dN/dpT)");
    ptDist[0]->Divide(run12hist[0]);
    pretty1DHist(ptDist[0], kRed, 24);
    Run12Ratio->cd(4);
    if( compareRun15p )
      ptDist[0]->GetYaxis()->SetRangeUser(0,4);
    else
      ptDist[0]->GetYaxis()->SetRangeUser(0,50);
    ptDist[0]->Draw("pe");
  }
}  

void getRun12Hists(TFile* f12, TH1D* h[3])
{
  TH1D* h1 = (TH1D*)f12->Get(Form("mh1electronPtTrg%i",trigSelect-1));
  TH2D* h2 = (TH2D*)f12->Get(Form("mPhi_ptUnlikeTrg%i",trigSelect-1));
  TH2D* h3 = (TH2D*)f12->Get(Form("mPhi_ptlikeTrg%i",trigSelect-1));
  h[0] = (TH1D*)h1->Clone();
  h[1] = h2->ProjectionX();
  h[2] = h3->ProjectionX();
  for(int l=0;l<3;l++) 
  {
    TH1D* temp = rebinVariableBins(h[l], numPtBinsEFF2, lowptEFF2);
    temp->Scale(1.,"width");
    h[l] = temp;
  }
  pretty1DHist(h[0],kGreen+3,24);
  pretty1DHist(h[1],kMagenta,25);
  pretty1DHist(h[2],kViolet+10,26);
}

void getRun15pHists(TFile* f15p, TH1D* h[3])
{
  bool inFileR = checkIfFileOpen(f15p);
  TH1D* h1 = (TH1D*)f15p->Get("hEPt");
  TH1D* h2 = (TH1D*)f15p->Get("hEEPt_US");
  TH1D* h3 = (TH1D*)f15p->Get("hEEPt_LS");
  h[0] = (TH1D*)h1->Clone();
  h[1] = (TH1D*)h2->Clone();
  h[2] = (TH1D*)h3->Clone();
  for(int l=0;l<3;l++) 
  {
    TH1D* temp = rebinVariableBins(h[l], numPtBinsEFF2, lowptEFF2);
    temp->Scale(1.,"width");
    h[l] = temp;
  }
  pretty1DHist(h[0],kGreen+3,24);
  pretty1DHist(h[1],kMagenta,25);
  pretty1DHist(h[2],kViolet+10,26);
}

void makeUnlikeMinusLikePartnerElectrons()
{
  nSigEPartE[2] = (TH2F*)nSigEPartE[0]->Clone();
  nSigEPartE[2]->Add(nSigEPartE[1],-1.);
  nSigEPartE[2]->SetTitle("Unlike - Like");
  pvePartE[2] = (TH2F*)pvePartE[0]->Clone();
  pvePartE[2]->Add(pvePartE[1],-1.);
  pvePartE[2]->SetTitle("Unlike - Like");
  nEtaPartE[2] = (TH2F*)nEtaPartE[0]->Clone();
  nEtaPartE[2]->Add(nEtaPartE[1],-1.);
  nEtaPartE[2]->SetTitle("Unlike - Like");
  nPhiPartE[2] = (TH2F*)nPhiPartE[0]->Clone();
  nPhiPartE[2]->Add(nPhiPartE[1],-1.);
  nPhiPartE[2]->SetTitle("Unlike - Like");
  phiDistPartE[2] = (TH2F*)phiDistPartE[0]->Clone();
  phiDistPartE[2]->Add(phiDistPartE[1],-1.);
  phiDistPartE[2]->SetTitle("Unlike - Like");
  zDistPartE[2] = (TH2F*)zDistPartE[0]->Clone();
  zDistPartE[2]->Add(zDistPartE[1],-1.);
  zDistPartE[2]->SetTitle("Unlike - Like");
  TPCTracks[2] = (TH1F*)TPCTracks[0]->Clone();
  TPCTracks[2]->Add(TPCTracks[1],-1.);
  TPCTracks[2]->SetTitle("Unlike - Like");
  EMCMatchedTracks[2] = (TH1F*)EMCMatchedTracks[0]->Clone();
  EMCMatchedTracks[2]->Add(EMCMatchedTracks[1],-1.);
  EMCMatchedTracks[2]->SetTitle("Unlike - Like");
  EMCIdTracks[2] = (TH1F*)EMCIdTracks[0]->Clone();
  EMCIdTracks[2]->Add(EMCIdTracks[1],-1.);
  EMCIdTracks[2]->SetTitle("Unlike - Like");
  SMDMatchedTracks[2] = (TH1F*)SMDMatchedTracks[0]->Clone();
  SMDMatchedTracks[2]->Add(SMDMatchedTracks[1],-1.);
  SMDMatchedTracks[2]->SetTitle("Unlike - Like");
  SMDIdTracks[2] = (TH1F*)SMDIdTracks[0]->Clone();
  SMDIdTracks[2]->Add(SMDIdTracks[1],-1.);
  SMDIdTracks[2]->SetTitle("Unlike - Like");
}

void makeInclDivUnlike(const char* fileName)
{
  TFile* f15pA = new TFile(fileName, "READ"); // pA name is carried over for naming compatibility, it just references current file.
  TFile* f15pp = new TFile(Form("/star/u/zamiller/PWGSpace/run15ppAnaTree/prod/anaTreeMaker_v3_100516/readTreeOut/out/histHadd/rootfile_temp/1117BHT%i/anaTree15_BHT%i.hists.root",trigSelect, trigSelect),"READ");
  TFile* f12pp = new TFile("/star/u/zamiller/PWGSpace/run15ppAnaTree/prod/anaTreeMaker_v2_080316/Run12RootFile/hist_5_2.root","READ");
  checkIfFileOpen(f15pA);
  checkIfFileOpen(f15pp);
  checkIfFileOpen(f12pp);

  TH1D* inclDivUS[3][2];
  inclDivUS[0][0] = (TH1D*)f15pA->Get("hEPt");
  inclDivUS[0][1] = (TH1D*)f15pA->Get("hEEPt_US");
  inclDivUS[1][0] = (TH1D*)f15pp->Get("hEPt");
  inclDivUS[1][1] = (TH1D*)f15pp->Get("hEEPt_US");
  inclDivUS[2][0] = (TH1D*)f12pp->Get(Form("mh1electronPtTrg%i",trigSelect-1));
  TH2F* tempH = (TH2F*)f12pp->Get(Form("mPhi_ptUnlikeTrg%i",trigSelect-1));
  inclDivUS[2][1] = tempH->ProjectionX();

  TString legName[3] = {"Run 15 p+p","Run 15 p+p","Run 12 p+p"};
  TLegend* legR = new TLegend(.11,.7,.6,0.87);
  inclUSRatio->cd();
  inclDivUS[0][0]->SetTitle("Inclusive Elec. / Unlike Sign Electrons; p_{T} (GeV/c); Ratio");
  inclDivUS[0][0]->GetYaxis()->SetRangeUser(0,50);
  for(int col=0; col<3; col++)
  {
    if(col == 1)
      continue;
    TH1D* temp0 = (TH1D*)inclDivUS[col][0]->Clone();
    TH1D* temp1 = (TH1D*)inclDivUS[col][1]->Clone();
    inclDivUS[col][0]= rebinVariableBins(temp0, numPtBinsEFF2, lowptEFF2);
    inclDivUS[col][1]= rebinVariableBins(temp1, numPtBinsEFF2, lowptEFF2);
    inclDivUS[col][0]->Scale(1.,"width");
    inclDivUS[col][1]->Scale(1.,"width");
    inclDivUS[col][0]->Divide(inclDivUS[col][1]);
    if(col==0)inclDivUS[col][0]->GetYaxis()->SetRangeUser(0,15);
    //inclDivUS[col][0]->Rebin(4);
    pretty1DHist(inclDivUS[col][0],colors[col],20+col);
    legR->AddEntry(inclDivUS[col][0],legName[col],"lpe");
    inclDivUS[col][0]->Draw((!col) ? "pe" : "same pe");
  }
  legR->Draw("same");
}

void makePDF(const char* fileName)
{
  //Set front page
  TCanvas* fp = new TCanvas("fp","Front Page",50,50,1050,1050);
  fp->cd();
  TBox *bLabel = new TBox(0.01, 0.88, 0.99, 0.99);
  bLabel->SetFillColor(38);
  bLabel->Draw();
  TLatex tl;
  tl.SetNDC();
  tl.SetTextColor(kWhite);
  tl.SetTextSize(0.033);
  char tlName[100];
  char tlName2[100];

  TString titlename = fileName;
  int found = titlename.Last('/');
  if(found >= 0){
    titlename.Replace(0, found+1, "");
  } 
  sprintf(tlName, "RUN 15 p+A 200 GeV NPE-Hadron");
  tl.SetTextSize(0.05);
  tl.SetTextColor(kWhite);
  tl.DrawLatex(0.05, 0.92,tlName);

  TBox *bFoot = new TBox(0.01, 0.01, 0.99, 0.12);
  bFoot->SetFillColor(38);
  bFoot->Draw();
  tl.SetTextColor(kWhite);
  tl.SetTextSize(0.05);
  tl.DrawLatex(0.05, 0.05, (new TDatime())->AsString());
  tl.SetTextColor(kBlack);
  tl.SetTextSize(0.03);
  tl.DrawLatex(0.1, 0.14, titlename);

  // Place canvases in order
  TCanvas* temp = new TCanvas();
  char name[100];
  sprintf(name, "%s.pdf[", fileName);
  temp->Print(name);
  sprintf(name, "%s.pdf", fileName);
  temp = fp; // print front page
  temp->Print(name);
  temp = eventHists;
  temp->Print(name);
  temp = dndpt;
  temp->Print(name);
  if(compareRun12){
    temp = Run12Compare;
    temp->Print(name);
    temp = Run12Ratio;
    temp->Print(name);
  }
  temp = inclUSRatio;
  temp->Print(name);
  temp = efficiencies;
  temp->Print(name);
  temp = npeYield;
  temp->Print(name);
  temp = crossSection;
  temp->Print(name);
  temp = rpa;
  temp->Print(name);
  temp = equivMB;
  temp->Print(name);
  temp = yieldComponents;
  temp->Print(name);
  for(int etype=0;etype<3;etype++)
  {
    for(int q=0; q<numCanvas; q++)
    {
      temp = dPhiPt[etype][q]; // print data canvases
      temp->Print(name);
    }
  }
  temp = delPhiDelEta;
  temp->Print(name);
  for(int q=0; q<numCanvas; q++)
  {
    temp = invMassPt[q]; // print data canvases
    temp->Print(name);
  }

  temp = zdcQA;
  temp->Print(name);
  temp = invMass;
  temp->Print(name);
  temp = invMassVsPt;
  temp->Print(name);
  temp = cutEfficiency;
  temp->Print(name);
  temp = ptCompare;
  temp->Print(name);
  temp = hadQA;
  temp->Print(name);
  temp = elecQA;
  temp->Print(name);
  for(int i=0;i<2;i++)
  {
    temp = eeQA[i];
    temp->Print(name);
  }
  for(int i=0;i<4;i++)
  {
    temp = eIDCutEffic[i];
    temp->Print(name);
  }
  temp = efficOverlay;
  temp->Print(name);
  for(int i=0;i<6;i++)
  {
    temp = pElecCuts[i];
    temp->Print(name);
  }
  for(int q=0; q<numCanvas; q++)
  {
    temp = nSigEPt[q]; // print data canvases
    temp->Print(name);
  }
  temp = nSigMeanSig;
  temp->Print(name);
  for(int q=0; q<numCanvas; q++)
  {
    temp = twoGaus[q]; 
    temp->Print(name);
  }
  for(int q=0; q<numCanvas; q++)
  {
    temp = nSigEff[q];
    temp->Print(name);
  }
  temp = nSigCutPlot;
  temp->Print(name);
  temp = eeOriginQA;
  temp->Print(name);

  sprintf(name, "%s.pdf]", fileName);
  temp->Print(name);
}

TH1F* rebinVariableBins(TH1F* h, int nbins, const float* bins, TString name)
{
  if(name.EqualTo("bob"))
    name = h->GetName();
  Double_t xbins[100];
  int segs = nbins-1;
  for(int ii=0;ii<nbins;ii++) xbins[ii] = bins[ii];
  TH1F* newH = (TH1F*)h->Rebin(segs,name,xbins);
  return newH;
}

TH1D* rebinVariableBins(TH1D* h, int nbins, const float* bins, TString name)
{
  if(name.EqualTo("bob"))
    name = h->GetName();
  Double_t xbins[100];
  int segs = nbins-1;
  for(int ii=0;ii<nbins;ii++) xbins[ii] = bins[ii];
  TH1D* newH = (TH1D*)h->Rebin(segs,name,xbins);
  return newH;
}

TH1F* convertToHist(TGraphAsymmErrors* gr,TString hname)
{
  TH1D* h1 = new TH1D(hname,"",200,0,20);
  TH1F* h = (TH1F*)rebinVariableBins(h1,numPtBinsEFF2,lowptEFF2,hname);
  for(Int_t i=0;i<numPtBinsEFF2;i++)
  {
    Double_t x=0, y=0, y_err=0;
    gr->GetPoint(i,x,y);
    y_err=gr->GetErrorY(i);

    h->SetBinContent(i+1,y);
    h->SetBinError(i+1,y_err);
  }
  return h;
}

void setFONLLStyle(TGraphErrors* g,int linestyle,int color)
{
  g->SetMarkerSize(0);
  g->SetMarkerStyle(1);
  g->SetLineStyle(linestyle);
  g->SetLineWidth(1);
  g->SetLineColor(color);
}
