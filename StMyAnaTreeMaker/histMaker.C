#include "histMaker.h"

bool DEBUG = !kTRUE;
bool compareRun12 = kTRUE;

void histMaker(const char* fileName="test", int trg=1)
{
  trigSelect = trg; // tell analyze what sample, see getTriggerLabel()
  getTriggerLabel();

  //Get TFile
  TFile* f = new TFile(fileName, "READ");
  TFile* f12 = new TFile("/star/u/zamiller/PWGSpace/run15ppAnaTree/prod/anaTreeMaker_v2_080316/Run12RootFile/hist_5_2.root","READ");
  bool inFile = checkIfFileOpen(f);
  stringstream ss;
  TString outfile; 
  ss << fileName;
  ss >> outfile;
  outfile.ReplaceAll(".root",Form("%s_processed.root",trigLabel));
  TFile* oF = new TFile(outfile, "RECREATE");
  bool outFile = checkIfFileOpen(oF);
  if(!inFile || !outFile) exit(1);


  getHistograms(f);
  declareHistograms();
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
    delPhiDelEta->cd(etype+1);
    gPad->SetLogz(1);
    eHadDelPhiDelEta[etype]->Draw("colz");
  }
  

  drawEventHists();
  drawInvMassHists();
  //drawQAHists();
  drawCutEfficiencyHists();
  getdNdpT(f12);
  makePDF(fileName);
  if(f12->IsOpen()) f12->Close();
  f->Close();
  oF->Close();
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

      Double_t xbins[numPtBinsEFF];
      int segs = numPtBinsEFF - 1;
      for(int ii=0;ii<numPtBinsEFF;ii++) xbins[ii] = lowptEFF[ii];
      TH1F* numerator = (TH1F*)histList[2]->Clone();
      TH1F* denominator = (TH1F*)TPCTracks[2]->Clone();
      TH1F* numRebin = numerator->Rebin(segs,"numRebin",xbins);
      TH1F* denRebin = denominator->Rebin(segs,"denRebin",xbins);
      //numerator->Rebin(20);
      //denominator->Rebin(20);
      numRebin->SetTitle("Cut Efficiency;P_{T} (GeV/c);Efficiency");
      pretty1DHist(numRebin,kRed,20);
      numRebin->GetYaxis()->SetRangeUser(0.,1.3);
      numRebin->Divide(denRebin);
      numRebin->Draw();
      partECutEfficiency[i-1] = (TH1F*)numRebin->Clone();
      partECutEfficiency[i-1]->SetTitle(titlename[i-1]);
    }
  }
  overlayEfficiencies();
  drawnSigE();
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
}

void pretty2DHist(TH2* h, int col, int style)
{
  //h->SetMarkerColor(col);
  h->SetMarkerStyle(style);
  h->SetMarkerSize(1.2);
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
        meanErr[ptbin]=cov[4];
        sigma[ptbin]=gaus->GetParameter(2);
        sigmaErr[ptbin]=cov[8];
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
  sg->Fit("pol1","R","",1.3,6.);
  mn->Draw("SAME PE");
  mn->Fit("pol1","R","",1.3,6.);
  leg->Draw("SAME");
  mn->Write();
  sg->Write();
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

    //TF1 *Gaus=new TF1("Gaus","exp(-0.5*pow((x-[0])/[1],2))/sqrt(2.*TMath::Pi())/[1]",-4,4);
    TF1 *Gaus=new TF1("Gaus","gausn(0)",-4,4);
    for(int j=0;j<10000;j++){
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
  grnSigCut->Fit("pol1","R","",1.3,6.);
  grnSigCut->Write();

}

void declareHistograms(){
  
  for(int ptbin=0; ptbin<numPtBins; ptbin++)
  {
    nSigEeff[ptbin] = new TH1F(Form("nSigEeff_%i",ptbin),"nSigE Cut Eff",1000,0,1);
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

 /* for(int ptbin = 0; ptbin < numPtBins; ptbin++){
    eHadNorm[ptbin] = 0;
    for(int etype=0; etype<3; etype++)
    {
      eHadNorm[ptbin] += trigCount[etype][ptbin];
    }
  }*/

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
  nSigCutPlot = new TCanvas("nSigCutPlot","nSigmaE Cut Efficiency",50,50,1050,1050);
  nSigMeanSig = new TCanvas("nSigMeanSig","nSigmaE Fit Results",50,50,1050,1050);
  dndpt = new TCanvas("dndpt","dndpt",50,50,1050,1050);
  Run12Compare = new TCanvas("Run12Compare","Run12Compare",50,50,1050,1050);
  Run12Compare->Divide(2,2);
  
  eventHists = new TCanvas("eventHists","Event Level Hists",50,50,1050,1050);
  eventHists->Divide(1,2);
  zdcQA = new TCanvas("zdcQA","ZDC QA Hists",50,50,1050,1050);
  zdcQA->Divide(1,2);
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

void getdNdpT(TFile* f12)
{
  TString histName[3] = {"inclusivePt","USPt","LSPt"}; 
  TString legName[3] = {"Inclusive Electron", "Unlike Sign Photonic", "Like Sign Photonic"};
  TString legName12[3] = {"Run 12 Inclusive", "Run 12 Unlike", "Run 12 Like"};

  Double_t xbins[numPtBins];
  int segs = numPtBins - 1;
  for(int ii=0;ii<numPtBins;ii++) xbins[ii] = lowpt[ii];

  TLegend* leg = new TLegend(.53,.74,.88,.88);
  dndpt->cd();
  gPad->SetLogy(1);
  for(int i=0; i<3; i++)
  {
    TH1F* ptSpectra = ePt[i];//->Rebin(segs,histName[i],xbins);
    ptSpectra->SetStats(0);
    pretty1DHist(ptSpectra,colors[i],20+i);
    leg->AddEntry(ptSpectra,legName[i],"lpe");
    ptSpectra->Rebin(4);
    ptSpectra->Scale(1./ptSpectra->GetXaxis()->GetBinWidth(5));
    if(i==0){
      ptSpectra->GetYaxis()->SetRangeUser(1,1e7);
      ptSpectra->SetTitle("Electron Spectra; p_{T} (GeV/c); Raw dN/dpT");
    }
    ptSpectra->Draw((i==0)?"pe":"same pe");
  }
  sampleLabel->Draw("SAME");

  if(compareRun12){
    TH1D* run12Hists[3];
    getRun12Hists(f12,run12Hists);
    run12Hists[0]->Draw("pe same");
    run12Hists[1]->Draw("pe same");
    run12Hists[2]->Draw("pe same");
    leg->AddEntry(run12Hists[0],legName12[0],"lpe");
    leg->AddEntry(run12Hists[1],legName12[1],"lpe");
    leg->AddEntry(run12Hists[2],legName12[2],"lpe");
  }
  leg->Draw("same");

  if(compareRun12){
    for(int i=0; i<3; i++)
    {
      TH1F* r15 = (TH1F*)ePt[i]->Clone();
      TH1F* r12 = (TH1F*)run12Hists[i]->Clone();
      double r15i = r15->Integral();
      double r12i = r12->Integral();
      r12->Scale(r15i/r12i);
      r15->SetTitle("Run Comparison - Scaled to Match Run 15 Yield");
      Run12Compare->cd(i+1);
      gPad->SetLogy(1);
      r15->Draw("pe");
      r12->Draw("same pe");
      TLegend* leg12 = new TLegend(.53,.74,.88,.88);
      leg12->AddEntry(r15, legName[i], "lpe");
      leg12->AddEntry(r12, legName12[i], "lpe");
      leg12->Draw("same");
    }
      f12->Close();
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
  for(int l=0;l<3;l++) h[l]->Rebin(4);
  Double_t bW = h[0]->GetXaxis()->GetBinWidth(10);
  Double_t bW2 = h[1]->GetXaxis()->GetBinWidth(10);
  Double_t bW3 = h[2]->GetXaxis()->GetBinWidth(10);
  h[0]->Scale(1./bW);
  h[1]->Scale(1./bW2);
  h[2]->Scale(1./bW3);
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
  temp = eventHists;
  temp->Print(name);
  temp = dndpt;
  temp->Print(name);
  if(compareRun12){
    temp = Run12Compare;
    temp->Print(name);
  }
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

