
#include <stdlib.h>
#include <iostream>
#ifndef __CINT__
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TPDF.h"
#include "TH2S.h"
#include "TH1D.h"
#include "TF1.h"
#include "TObject.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#endif 

void setpad(TVirtualPad *pad,float left, float right, float top, float bottom, int logy = 0, int logz = 0);
TLegend* myLeg(Double_t xlow, Double_t ylow, Double_t xup, Double_t yup,Int_t textFont=42,Double_t textSize=0.05);
bool isHFT(TString name);
bool isRunIndex(TString name);
TH1F* makeMeanHist(TH2F* h, float &currentHistStDev, TString outpdf);
void getErrorVal(TH1* h, TH2* h2, float& curretHistStDev, float& newMean, float lowLim, float highLim, TString outpdf);
TH1F* findBadRuns(TH1F* h, float lowLim, float highLim, int* runList, TPaveText* lbl, TString badRunList);
bool getRunList(int* runList);
void addToBadRunList(int run, TString badRunList);

void makeRunIndexQA(const char* FileName="test.picoHFMyAnaMaker.root")
{

    gROOT->Reset();
    setcolz();
    setstyle();   
    gStyle->SetOptStat(1);
    gStyle->SetOptFit(1);
    gStyle->SetOptTitle(1);
    TGaxis::SetMaxDigits(3);

    int runList[834];
    bool haveRunlist= getRunList(runList);
    TString badRunList = "badRunList.txt";
    TString nameCheck = FileName;
    if(nameCheck.Index("hists") > 0)
      badRunList.ReplaceAll(".txt","AnaTreeCheck.txt");
    if( remove( badRunList ) != 0 )
      perror( "Error deleting badRunList.txt" );
    else
      puts( "Old badRunList.txt successfully deleted" );

    char name[100];
    sprintf(name, "%s", FileName);
    TFile *f = new TFile(name);

    TCanvas *c1 = new TCanvas("c1", "c1",0,0,800,800);
    //setpad(c1,0.001,0.001,0.001,0.001, 0, 0);
    setpad(c1,0.1,0.1,0.05,0.12, 0, 0);
    TPad *pad = new TPad("pad","",0.00,0.00,1.00,1.00);
    pad->Draw();

    TIter nextkey(f->GetListOfKeys());
    TObject *t;

    TString outpdf = FileName;
    outpdf.ReplaceAll(".root","RunLevelQA");
   // TPDF *pdf = new TPDF(outpdf.Data());
    TCanvas* temp = new TCanvas("temp","temp",50,50,1050,1050);
    temp->Print(Form("%s.pdf[",outpdf.Data()));
    //Set front page
    c1->cd();
    //pad->cd();
    //setpad(gPad,0.1,0.1,0.05,0.12, 0, 0);
    TBox *bLabel = new TBox(0.01, 0.88, 0.99, 0.99);
    bLabel->SetFillColor(kBlack);
    bLabel->Draw();
    TLatex tl;
    tl.SetNDC();
    tl.SetTextColor(kWhite);
    tl.SetTextSize(0.063);
    char tlName[100];
    char tlName2[100];

    TString titlename = FileName;
    int found = titlename.Last('/');
    if(found >= 0){
        titlename.Replace(0, found+1, "");
    } 
    titlename.ReplaceAll(".root"," Run Level QA");
    sprintf(tlName, "%s", titlename.Data());
    tl.SetTextSize(0.05);
    tl.SetTextColor(kBlack);
    tl.DrawLatex(0.05, 0.7,tlName);

    TBox *bFoot = new TBox(0.01, 0.01, 0.99, 0.12);
    bFoot->SetFillColor(kBlack);
    bFoot->Draw();
    tl.SetTextColor(kWhite);
    tl.SetTextSize(0.05);
    tl.DrawLatex(0.05, 0.05, (new TDatime())->AsString());
    tl.SetTextColor(kBlack);
    tl.SetTextSize(0.04);
    //tl.DrawLatex(0.1, 0.55, titlename);
    temp = c1;
    temp->Print(Form("%s.pdf",outpdf.Data()));
  //  c1->cd();
  //  c1->Update();

    c1->cd();
    int nCounts = 0;
    int nFig = 0;
    TKey *key;
    TLine line;
    line.SetLineStyle(7);
    line.SetLineWidth(1);
    TString tName;
    tl.SetTextSize(0.035);
    float currentHistStDev = 0.;
    TLatex tLabel;
    tLabel.SetNDC();
    tLabel.SetTextColor(kRed);
    tLabel.SetTextSize(0.05);
    TString label = "Mean #pm 4*RMS";

    TPaveText* brLabel = new TPaveText(.40,.64,.65,.89,"NB NDC");
    brLabel->SetFillColorAlpha(kWhite,0);
    brLabel->SetTextSize(0.03);
    brLabel->SetTextColor(kMagenta);
    brLabel->SetTextAlign(11);
    

    //Loop through TH1 and TH2 histos and place them on PDF
    while ((key = (TKey*)nextkey())) {
        t = dynamic_cast<TObject*>(key->ReadObj());
        if(t){
          tName = t->GetName();
          if(isHFT(tName)||!isRunIndex(tName)) continue;
          if(!strcmp(t->ClassName(),"TH2D")||!strcmp(t->ClassName(),"TH2F")||!strcmp(t->ClassName(),"TH2S")){
                TH2F *h2 = (TH2F*)t;
                TH1F* hM1 = makeMeanHist(h2,currentHistStDev,outpdf);
                TF1* poly = new TF1("poly","pol0",0,832);
                TF1* polyH = new TF1("polyH","pol0",0,832);
                TF1* polyL = new TF1("polyL","pol0",0,832);
                poly->SetLineColor(kRed);
                polyH->SetLineColor(kRed);
                polyL->SetLineColor(kRed);
                polyH->SetLineStyle(9);
                polyL->SetLineStyle(9);
                hM1->Fit("poly","RWQ");
                float p0 = poly->GetParameter(0);
                float err = currentHistStDev;


                TH1F* hM1c = (TH1F*)hM1->Clone("hM1c");
                TCanvas *c2 = new TCanvas("c2", "c2",0,0,1000,1000);
                c2->Divide(1,2);
                c2->cd(1);
                hM1c->Draw("pe");
                polyH->Draw("same");
                polyL->Draw("same");
                tLabel.DrawLatex(0.16,0.85,label);
                TH1F* badRuns = findBadRuns(hM1,p0-4*err,p0+4*err, runList, brLabel, badRunList);
                double newMean = p0;
                for(int iter=0; iter<4; iter++)
                {
                  brLabel->Clear();
                  brLabel->AddText("Bad Runs:");
                  getErrorVal(hM1, h2, currentHistStDev, newMean, newMean-4*err, newMean+4*err,  outpdf); 
                  err = currentHistStDev;
                  badRuns = findBadRuns(hM1,newMean-4*err,newMean+4*err, runList, brLabel, badRunList);
                }
                c2->cd(1);
                badRuns->GetYaxis()->SetRangeUser(newMean-9*err,newMean+9*err);
                badRuns->Draw("same pe");
                brLabel->Draw("same");
                polyH->SetParameter(0,newMean+4*err);
                polyL->SetParameter(0,newMean-4*err);
                polyH->Draw("same");
                polyL->Draw("same");
                
                c2->cd(2);
                hM1->SetStats(kFALSE);
                hM1->GetYaxis()->SetRangeUser(p0-8*err,p0+8*err);
                hM1->Draw("pe");
                polyH->Draw("same");
                polyL->Draw("same");
                
                nCounts++;
                nFig++;
                temp = c2;
                temp->Print(Form("%s.pdf",outpdf.Data()));
                c2->Delete();
            }
            if(!strcmp(t->ClassName(),"TH1D")||!strcmp(t->ClassName(),"TH1F")||!strcmp(t->ClassName(),"TH1S")){
              continue;
            }
            if(!strcmp(t->ClassName(),"TList")){
                TIter next((TList*)t); 
                TObject *obj = 0;
                while ((obj = next())){ 
                    if(!strcmp(obj->ClassName(),"TH2D")||!strcmp(obj->ClassName(),"TH2F")||!strcmp(obj->ClassName(),"TH2S")){
                        TH2F *hh2 = (TH2F*)obj;
                        TH1F* hM1 = makeMeanHist(h2);
                        TF1* poly = new TF1("poly","pol1",0,810);
                        hM1->Fit("poly","RQ");
                        hM1->Draw("pe");
                        //c1->cd();
                        //c1->Update();
                    }
                    if(!strcmp(obj->ClassName(),"TH1D")||!strcmp(obj->ClassName(),"TH1F")||!strcmp(obj->ClassName(),"TH1S")){
                      continue;
                    }
                }
            }
        }
    }

    //pdf->Close();
    temp->Print(Form("%s.pdf]",outpdf.Data()));
    f->Close(); 
    return;
}

TH1F* findBadRuns(TH1F* h, float lowLim, float highLim, int* runList, TPaveText* lbl, TString badRunList)
{
  int maxBins = h->GetXaxis()->GetNbins();
  int lowX = h->GetXaxis()->GetBinLowEdge(1);
  int highX = h->GetXaxis()->GetBinUpEdge(maxBins);
  TH1F* newH = new TH1F("name","title",maxBins,lowX,highX);
  for(int i=1;i<maxBins+1;i++)
  {
    if(i>=832) continue;
    float binval = h->GetBinContent(i);
    if(binval <= lowLim || binval >= highLim && i <= 831)
    {
      newH->SetBinContent(i,binval+1e-6);
      //h->SetBinContent(i,1e7);
      if(i<=831 && (highLim !=0 || lowLim != 0)){
        addToBadRunList(runList[i-1], badRunList);
        if(runList[i-1] == 0) cout << "0 caused by id: " << i-1 << endl;
        lbl->AddText(Form("%i",runList[i-1]));
        cout << "badrun: " << runList[i-1] << endl;
      }
    }
    else 
      newH->SetBinContent(i,1e7);
  }
  newH->SetMarkerStyle(20);
  newH->SetMarkerColor(kMagenta);
  newH->SetMarkerSize(2.0);
  return newH;
}

void addToBadRunList(int run, TString badRunList)
{
  ifstream indata;
  indata.open(badRunList);
  if(indata.is_open()){
    Int_t runnum;
    while(indata>>runnum){
      if(run == runnum)
        return;
    }
  }
  std::ofstream ofs(badRunList,std::ofstream::app);
  ofs << run << endl;
  ofs.close();
}

bool getRunList(int* runList)
{
  ifstream indata;
  indata.open("runNumberList_run15pA_BACKUP092016");
  if(indata.is_open()){
    cout<<"read in total run number list and recode run number ...";
    Int_t runnum;
    Int_t newId=0;
    while(indata>>runnum){
      runList[newId] = runnum;
      newId++;
    }
    cout<<" [OK]"<<endl;  
    indata.close();
    return kTRUE;
  }else{
    cout<<"Failed to load the total run number list !!!"<<endl;
    indata.close();
    return kFALSE;
  }

}

TH1F* makeMeanHist(TH2F* h, float &currentHistStDev, TString outpdf)
{
  
  int maxBins = h->GetXaxis()->GetNbins();
  int lowX = h->GetXaxis()->GetBinLowEdge(1);
  int highX = h->GetXaxis()->GetBinUpEdge(maxBins);
  int maxBinsY = h->GetYaxis()->GetNbins();
  float lowY = h->GetYaxis()->GetBinLowEdge(1);
  float highY = h->GetYaxis()->GetBinUpEdge(maxBinsY);
  
  std::ostringstream oss,oss2,oss3,oss4;
  oss << h->GetName() << ";" << h->GetXaxis()->GetTitle() << ";Mean " << h->GetYaxis()->GetTitle();
  oss4 << h->GetName() << "_getErr;" << h->GetYaxis()->GetTitle() << ";Counts";
  TString title = oss.str();
  TString title2 = oss4.str();
  oss2 << h->GetName() << "_MEAN";
  oss3 << h->GetName() << "_getErr";
  TString name = oss2.str();
  TString name2 = oss3.str();
  TH1F* newH = new TH1F(name,title,maxBins,lowX,highX);
  newH->Sumw2();

  for(int i=0;i<maxBins;i++)
  {
     TH1D* tempH = h->ProjectionY("tempH",i+1,i+2);
     Double_t mean = tempH->GetMean();
     if(i < 832){
       newH->SetBinContent(i,mean);
     }
  }
  TF1* polT = new TF1("polT","pol0",0,832);
  newH->Fit("polT","RWQ");
  float p0 = polT->GetParameter(0);
  float newMean = p0;
  getErrorVal(newH,h,currentHistStDev, newMean, -1e6, 1e6, outpdf);
  newH->SetMarkerStyle(20);
  newH->SetMarkerColor(kBlack);
  newH->SetMarkerSize(0.7);
  newH->GetYaxis()->SetRangeUser(lowY,highY);
  return newH;
}

void getErrorVal(TH1* h, TH2* h2, float& currentHistStDev, float& newMean, float lowLim, float highLim, TString outpdf)
{
  TCanvas *c3 = new TCanvas("c3", "c3",0,0,1000,1000);
  int maxBins = h->GetXaxis()->GetNbins();
  int lowX = h->GetXaxis()->GetBinLowEdge(1);
  int highX = h->GetXaxis()->GetBinUpEdge(maxBins);
  int maxBinsY = h2->GetYaxis()->GetNbins();
  float lowY = h2->GetYaxis()->GetBinLowEdge(1);
  float highY = h2->GetYaxis()->GetBinUpEdge(maxBinsY);

  std::ostringstream oss,oss2,oss3,oss4;
  oss4 << h2->GetName() << "_getErr;" << h2->GetYaxis()->GetTitle() << ";Counts";
  TString title2 = oss4.str();
  oss3 << h2->GetName() << "_getErr";
  TString name2 = oss3.str();
  TH1F* getErr = new TH1F(name2,title2,10*maxBinsY,lowY,highY);

  for(int i=1;i<maxBins;i++)
  {
    double val = h->GetBinContent(i);
    if(val >= lowLim && val <= highLim && val != 0)
    //if(fabs(val) < 1e6 && val != 0)
      getErr->Fill(val);
  }
  currentHistStDev = getErr->GetRMS();
  newMean = getErr->GetMean();
  c3->cd();
  gStyle->SetOptStat(111100);
  getErr->SetStats(kTRUE);
  getErr->Draw();
  c3->Print(Form("%s.pdf",outpdf.Data()));
  c3->Delete();
 
}

bool isHFT(TString name)
{
  if(name.Index("wHFT")>0)
    return true;
  return false;
}

bool isRunIndex(TString name)
{
  if(name.Index("RunIndex")>0)
    return true;
  return false;
}

void setpad(TVirtualPad *pad,float left, float right, float top, float bottom, int logy, int logz){
    pad->SetFillColor(10);
    pad->SetBorderMode(0);
    pad->SetBorderSize(0);
    pad->SetFrameFillColor(10);
    pad->SetFrameBorderMode(0);
    pad->SetFrameBorderSize(0);
    pad->SetLeftMargin(left);
    pad->SetRightMargin(right);
    pad->SetTopMargin(top);
    pad->SetBottomMargin(bottom);
    pad->SetLogy(logy);
    pad->SetLogz(logz);
    return;
}
TLegend* myLeg(Double_t xlow, Double_t ylow, Double_t xup, Double_t yup, Int_t textFont, Double_t textSize)
{
    TLegend *leg = new TLegend(xlow,ylow,xup,yup);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetFillColor(10);
    leg->SetTextFont(textFont);
    leg->SetTextSize(textSize);
    return leg;
}
void setcolz()
{
    const Int_t NRGBs=5;
    const Int_t NCont=32;

    Double_t stops[NRGBs]={0.00,0.34,0.61,0.84,1.00};
    Double_t red[NRGBs]={0.00,0.00,0.87,1.00,0.51};
    Double_t green[NRGBs]={0.00,0.81,1.00,0.20,0.00};
    Double_t blue[NRGBs]={0.51,1.00,0.12,0.00,0.00};
    TColor::CreateGradientColorTable(NRGBs,stops,red,green,blue,NCont);
    gStyle->SetNumberContours(NCont);
}
setstyle()
{
    TStyle* myStyle = new TStyle("myStyle","Styles");
    myStyle->SetPalette(1,0); // avoid horrible default color scheme
    myStyle->SetOptStat(111100);
    myStyle->SetOptTitle(1);
    myStyle->SetOptDate(0);
    myStyle->SetLabelSize(0.045,"xyz"); // size of axis value font
    myStyle->SetTitleX(0.2f);
    myStyle->SetTitleY(0.96f);
    myStyle->SetTitleW(0.5f);
    myStyle->SetTickLength(0.01,"xyz");
    myStyle->SetTitleFont(62,"xyz"); // font option 
    myStyle->SetLabelFont(62,"xyz");
    myStyle->SetTitleOffset(0.8,"z");
    myStyle->SetTitleOffset(1.2,"y");
    myStyle->SetTitleFillColor(10);
    myStyle->SetLineWidth(2);
    myStyle->SetCanvasDefW(700);
    myStyle->SetCanvasDefH(600);
    myStyle->SetCanvasColor(0);// canvas...
    myStyle->SetCanvasBorderMode(0);
    myStyle->SetCanvasBorderSize(0);
    myStyle->SetPadColor(0);
    myStyle->SetPadBorderSize(1);
    myStyle->SetPadBorderMode(-1);
    myStyle->SetPadBottomMargin(0.14); //margins...
    myStyle->SetPadTopMargin(0.06);
    myStyle->SetPadLeftMargin(0.14);
    myStyle->SetPadRightMargin(0.04);
    myStyle->SetPadGridX(0); // grids, tickmarks
    myStyle->SetPadGridY(0);
    myStyle->SetPadTickX(1);
    myStyle->SetPadTickY(1);
    myStyle->SetFrameBorderSize(1);
    myStyle->SetFrameBorderMode(-1);
    myStyle->SetFrameFillColor(0);
    myStyle->SetFrameLineWidth(1.2);
    myStyle->SetPaperSize(20,24); // US letter size
    gROOT->SetStyle("myStyle");
    cout << "Styles are Set!" << endl;

}

