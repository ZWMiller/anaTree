
void checkPhiEta(TString filename = "test.root")
{
  TFile* f = new TFile(filename,"READ");
  TH2F* etaPhi[11];
  TH1D* eta[11];
  TH1D* phi[11];
  TCanvas* etaphiC = new TCanvas("etaphiC","etaphiC",50,50,1050,1050);
  TCanvas* etaC = new TCanvas("etaC","etaC",50,50,1050,1050);
  TCanvas* etaCz = new TCanvas("etaCz","etaCz",50,50,1050,1050);
  TCanvas* etaCzz = new TCanvas("etaCzz","etaCzz",50,50,1050,1050);
  TCanvas* phiC = new TCanvas("phiC","phiC",50,150,1050,1050);
  TCanvas* phiCz = new TCanvas("phiCz","phiCz",50,150,1050,1050);
  TCanvas* phiCzz = new TCanvas("phiCzz","phiCzz",50,150,1050,1050);

  char cutNames[11][20] = {"NoCut","pTCut","etaCut","nHitFitCut","nHistdEdxCut","nSigECut","pveCut","nEtaCut","nPhiCut","zDistCut","phiDistCut"};
  int colors[11] = {kBlack,kGreen+3,kBlue,kRed,kTeal,kMagenta,kOrange,kAzure+1,kYellow,kPink-6,kViolet+3};

  TLegend* leg = new TLegend(0.6,0.55,0.89,0.88);
  for(int i=0;i<11;i++)
  {
    etaPhi[i] = (TH2F*) f->Get(Form("electronEtaPhi_%i",i));
    eta[i] = etaPhi[i]->ProjectionX();
    phi[i] = etaPhi[i]->ProjectionY();
    if(i==0){
      eta[i]->SetTitle("Electron Eta After Various Cuts");
      phi[i]->SetTitle("Electron Phi After Various Cuts");
    }
    eta[i]->SetMarkerColor(colors[i]);
    eta[i]->SetMarkerStyle(20+i);
    eta[i]->SetMarkerSize(1.0);
    eta[i]->SetStats(0);
    phi[i]->SetMarkerColor(colors[i]);
    phi[i]->SetMarkerStyle(20+i);
    phi[i]->SetMarkerSize(1.0);
    phi[i]->SetStats(0);
    phi[i]->GetXaxis()->SetRangeUser(-2.0,5.);
    leg->AddEntry(phi[i],cutNames[i],"lpe");
    etaC->cd();
    eta[i]->GetYaxis()->SetRangeUser(0,2000);
    eta[i]->DrawClone((i==0)?"pe":"same pe");
    phiC->cd();
    phi[i]->GetYaxis()->SetRangeUser(0,4200);
    phi[i]->DrawClone((i==0)?"pe":"same pe");
    etaCz->cd();
    eta[i]->GetYaxis()->SetRangeUser(0,550);
    eta[i]->DrawClone((i==0)?"pe":"same pe");
    phiCz->cd();
    phi[i]->GetYaxis()->SetRangeUser(0,1100);
    phi[i]->DrawClone((i==0)?"pe":"same pe");
    etaCzz->cd();
    eta[i]->GetYaxis()->SetRangeUser(0,80);
    eta[i]->DrawClone((i==0)?"hist p":"same hist p");
    phiCzz->cd();
    phi[i]->GetYaxis()->SetRangeUser(0,150);
    phi[i]->DrawClone((i==0)?"hist p":"same hist p");
  }
  etaC->cd();
  leg->Draw("same");
  phiC->cd();
  leg->Draw("same");
  etaCz->cd();
  leg->Draw("same");
  phiCz->cd();
  leg->Draw("same");

  TString pdfName = filename;
  pdfName.ReplaceAll(".root",".pdf");
  TCanvas* temp = new TCanvas();
  char name[100];
  sprintf(name, "%s.pdf[", pdfName.Data());
  temp->Print(name);
  sprintf(name, "%s.pdf", pdfName.Data());
  temp->Print(name);
  temp = etaC;
  temp->Print(name);
  temp = etaCz;
  temp->Print(name);
  temp = etaCzz;
  temp->Print(name);
  temp = phiC;
  temp->Print(name);
  temp = phiCz;
  temp->Print(name);
  temp = phiCzz;
  temp->Print(name);
  sprintf(name, "%s.pdf]", pdfName.Data());
  temp->Print(name);
}
