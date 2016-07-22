#include "histMaker.h"

bool DEBUG = !kTRUE;

void histMaker(const char* fileName="test")
{
  //Get TFile
  TFile* f = new TFile(fileName, "READ");
  if(!f->IsOpen())
  { std::cout << "!!! File Not Found !!!" << std::endl;
    exit(1); }  
  else
    std::cout<< "Opened " << fileName << std::endl;

  getHistograms(f);
  doProjections();
  prepareCanvas();
  prepareLabels();
  char textLabel[100];
  for(int etype=0; etype<3; etype++) 
  {
    for(int ptbin=0; ptbin<numPtBins; ptbin++)
    {
      int activeCanvas = (int) ptbin/9;
      int activeBin = ptbin - activeCanvas*9; 
      drawDeltaPhi(eHadDelPhi[etype][ptbin],dPhiPt[etype][activeCanvas],trigCount[etype][ptbin],activeBin,lbl[ptbin]); 
      if(DEBUG) cout << "pt loop: " << ptbin <<endl;
    } 
  }

  drawInvMassHists();
  drawQAHists();
  drawCutEfficiencyHists();
  makePDF(fileName);
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
}

void drawPartECutEffic()
{
  TH1F* histList[3];
  char titlename[4][100] = {"EMC Matched","EMC eID","SMD Matched","SMD eID"};
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
      numerator->Rebin(4);
      denominator->Rebin(4);
      numerator->SetTitle("Cut Efficiency;P_{T} (GeV/c);Efficiency");
      pretty1DHist(numerator,kRed,20);
      numerator->GetYaxis()->SetRangeUser(0.,1.3);
      numerator->Divide(denominator);
      numerator->Draw();
      partECutEfficiency[i-1] = (TH1F*)numerator->Clone();
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

void pretty1DHist(TH1* h, int col, int style)
{
  h->SetMarkerColor(col);
  h->SetLineColor(col);
  h->SetMarkerStyle(style);
  h->SetMarkerSize(0.8);
}

void pretty1DHistFill(TH1* h, int col, int style)
{
  h->SetLineColor(kBlack);
  h->SetMarkerStyle();
  h->SetFillStyle(style);
  h->SetFillColor(col);
}

void prepareLabels()
{
  char textLabel[100];

  for(int ptbin=0; ptbin<numPtBins; ptbin++)
  {
    lbl[ptbin] = new TPaveText(.53,.72,.75,.8,Form("NB NDC%i",ptbin));
    sprintf(textLabel,"%.2f < P_{T,e} < %.2f",lowpt[ptbin],highpt[ptbin]);
    lbl[ptbin]->AddText(textLabel);
    lbl[ptbin]->SetFillColor(kWhite);
  }
}

void drawInvMassHists()
{
  TLegend* leg = new TLegend(0.35,0.65,0.75,0.85);
  TString legLab[3] = {"Unlike Sign","Like Sign","US - LS"};
  for(int pairType=0; pairType<3; pairType++) //0 = UnlikeSign, 1 = LikeSign ++&-- Combined, 3 = US - LS
  {
    invMass->cd();
    gPad->SetLogy();
    eeInvMassAll[pairType]->SetLineColor(colors[pairType]);
    eeInvMassAll[pairType]->SetMarkerColor(colors[pairType]);
    eeInvMassAll[pairType]->SetMarkerStyle(20+pairType);
    eeInvMassAll[pairType]->SetMarkerSize(0.8);
    eeInvMassAll[pairType]->GetXaxis()->SetRangeUser(0.,0.5);
    eeInvMassAll[pairType]->Draw((pairType==0) ? "pe" : "same pe");
    leg->AddEntry(eeInvMassAll[pairType],legLab[pairType],"lpe"); 
    if(pairType==2) leg->Draw("same");

    for(int ptbin=0; ptbin<numPtBins; ptbin++)
    {
      int activeCanvas = (int) ptbin/9;
      int activeBin = ptbin - activeCanvas*9; 
      invMassPt[activeCanvas]->cd(activeBin+1);
      gPad->SetLogy();
      eeInvMass[pairType][ptbin]->SetLineColor(colors[pairType]);
      eeInvMass[pairType][ptbin]->SetMarkerColor(colors[pairType]);
      eeInvMass[pairType][ptbin]->SetMarkerStyle(20+pairType);
      eeInvMass[pairType][ptbin]->SetMarkerSize(0.8);
      eeInvMass[pairType][ptbin]->GetXaxis()->SetRangeUser(0.,0.5);
      eeInvMass[pairType][ptbin]->Draw((pairType==0) ? "pe" : "same pe");
      if(pairType==2 && ptbin == 0) leg->Draw("same");
    }
  }
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
  for(int ptbin = 0; ptbin < numPtBins; ptbin++){
    eHadNorm[ptbin] = 0;
    for(int etype=0; etype<3; etype++)
    {
      eHadNorm[ptbin] += trigCount[etype][ptbin];
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

  for(int q=0; q<numCanvas; q++)
  {
    invMassPt[q] = new TCanvas(Form("invMassPt_%i",q),"Invariant Mass Spectrum",50,50,1050,1050);
    invMassPt[q]->Divide(3,3);
  }
  invMass     = new TCanvas("invMass","Invariant Mass All pT",50,50,1050,1050);
  ptCompare     = new TCanvas("ptCompare","pT Comparison",50,50,1050,1050);
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
  eHadDelPhiPt[0] = (TH2F*)f->Get("hHadEDelPhiPt");
  eHadDelPhiPt[1] = (TH2F*)f->Get("hHadEEDelPhiPt_US");
  eHadDelPhiPt[2] = (TH2F*)f->Get("hHadEEDelPhiPt_LS");
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

  if(DEBUG) cout << "Get Hist." << endl;
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

void makePDF(const char* fileName)
{
  //Set front page
  TCanvas* fp = new TCanvas("fp","Front Page",100,0,1000,900);
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
  sprintf(tlName, "RUN 15 p+p 200 GeV NPE-Hadron");
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
  for(int etype=0;etype<3;etype++)
  {
    for(int q=0; q<numCanvas; q++)
    {
      temp = dPhiPt[etype][q]; // print data canvases
      temp->Print(name);
    }
  }
  for(int q=0; q<numCanvas; q++)
  {
    temp = invMassPt[q]; // print data canvases
    temp->Print(name);
  }

  temp = invMass;
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
  temp = eeOriginQA;
  temp->Print(name);

  sprintf(name, "%s.pdf]", fileName);
  temp->Print(name);
}
