
void checkPhiEta(TString filename = "test.root")
{
  TFile* f = new TFile(filename,"READ");
  TH2F* etaPhi[12];
  TH1D* eta[12];
  TH1D* phi[12];
  TCanvas* etaphiC = new TCanvas("etaphiC","etaphiC",50,50,1050,1050);
  TCanvas* etaC = new TCanvas("etaC","etaC",50,50,1050,1050);
  TCanvas* etaCz = new TCanvas("etaCz","etaCz",50,50,1050,1050);
  TCanvas* etaCzz = new TCanvas("etaCzz","etaCzz",50,50,1050,1050);
  TCanvas* phiC = new TCanvas("phiC","phiC",50,150,1050,1050);
  TCanvas* phiCz = new TCanvas("phiCz","phiCz",50,150,1050,1050);
  TCanvas* phiCzz = new TCanvas("phiCzz","phiCzz",50,150,1050,1050);

  char cutNames[12][20] = {"NoCut","pTCut","etaCut","nHitFitCut","nHistdEdxCut","nSigECut","pveCut","nEtaCut","nPhiCut","zDistCut","phiDistCut","dsmAdcCut"};
  int colors[12] = {kBlack,kGreen+3,kBlue,kRed,kTeal,kMagenta,kOrange,kAzure+1,kOrange-7,kPink-6,kViolet+3,kBlue-9};

  TLegend* leg = new TLegend(0.6,0.55,0.89,0.88);
  for(int i=0;i<12;i++)
  {
    etaPhi[i] = (TH2F*) f->Get(Form("electronEtaPhi_%i",i));
    eta[i] = etaPhi[i]->ProjectionX();
    phi[i] = etaPhi[i]->ProjectionY();
    eta[i]->Rebin(2);
    phi[i]->Rebin(2);
    if(i==0){
      eta[i]->SetTitle("Electron Eta After Various Cuts");
      phi[i]->SetTitle("Electron Phi After Various Cuts");
    }
    eta[i]->SetMarkerColor(colors[i]);
    eta[i]->SetLineColor(colors[i]);
    eta[i]->SetMarkerStyle(20+i);
    eta[i]->SetMarkerSize(1.2);
    eta[i]->SetStats(0);
    phi[i]->SetMarkerColor(colors[i]);
    phi[i]->SetLineColor(colors[i]);
    phi[i]->SetMarkerStyle(20+i);
    phi[i]->SetMarkerSize(1.2);
    phi[i]->SetStats(0);
    phi[i]->GetXaxis()->SetRangeUser(-2.0,5.);
    leg->AddEntry(phi[i],cutNames[i],"p");
    etaC->cd();
    if(i==0)eta[i]->GetYaxis()->SetRangeUser(0,eta[0]->GetMaximum()*2.);
    eta[i]->DrawClone((i==0)?"pe":"same pe");
    phiC->cd();
    if(i==0)phi[i]->GetYaxis()->SetRangeUser(0,phi[0]->GetMaximum()*2.);
    phi[i]->DrawClone((i==0)?"pe":"same pe");
    if(i>=1){
      etaCz->cd();
      eta[i]->GetYaxis()->SetRangeUser(0,eta[1]->GetMaximum()*2.);
      eta[i]->DrawClone((i==1)?"pe":"same pe");
      phiCz->cd();
      phi[i]->GetYaxis()->SetRangeUser(0,phi[1]->GetMaximum()*2.);
      phi[i]->DrawClone((i==1)?"pe":"same pe");
    }
    if(i>=5){
      etaCzz->cd();
      cout << eta[5]->GetMaximum() << " " << eta[5]->GetBinCenter(eta[5]->GetMaximumBin()) << endl;
      if(i==5)eta[i]->GetYaxis()->SetRangeUser(0,eta[5]->GetMaximum()/3.);
      eta[i]->DrawClone((i==5)?"hist p":"same hist p");
      phiCzz->cd();
      if(i==5)phi[i]->GetYaxis()->SetRangeUser(0,phi[5]->GetMaximum()/3.);
      phi[i]->DrawClone((i==5)?"hist p":"same hist p");
    }

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
  pdfName.ReplaceAll(".root",".EtaPhiCheck");
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
