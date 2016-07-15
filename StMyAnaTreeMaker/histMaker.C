#include "histMaker.h"

bool DEBUG = kTRUE;

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
    if(etype==2)
    {
      eeInvMassAll[1]->Add(eeInvMassAll[2]); // combine ++ and -- to form LS
      eeInvMassAll[2]->Add(eeInvMassAll[0],eeInvMassAll[1],1.,-1.);
    }
  }


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

  if(DEBUG) cout << "Get Hist." << endl;
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

  sprintf(name, "%s.pdf]", fileName);
  temp->Print(name);
}

