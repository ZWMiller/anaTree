#ifndef _histMaker
#define _histMaker

    const int numPtBins = 8;
    const int numCanvas = numPtBins/9 + 1;

    const float lowpt[numPtBins]  = 
     {1.3, 2.0, 3.0, 4.0, 
      5.0, 6.0, 8.0, 10.0};
    const float highpt[numPtBins] = 
     {2.0, 3.0, 4.0, 5.0,
      6.0, 8.0, 10., 20. };
    const int colors[5] = {kBlack,kRed,kAzure+1,kMagenta,kYellow};


    TCanvas* dPhiPt[3][numCanvas];
    TCanvas* invMassPt[numCanvas];
    TCanvas* invMass;
    TCanvas* ptCompare;
    TCanvas* cutEfficiency;
    TCanvas* hadQA;
    TCanvas* eeQA[2];
    TCanvas* eeOriginQA;
    TCanvas* elecQA;
    TCanvas* pElecCuts[6];
    TCanvas* eIDCutEffic[4];
    TCanvas* efficOverlay;
    TCanvas* eventHists;

    TPaveText* lbl[numPtBins];

    TH2F* eHadDelPhiPt[3];
    TH2F* eeInvMassPt[3];
    TH2F* hadPtEPt;
    TH1F* ePt[3];
    TH1F* ePt_eff[4];
    TH1D* eHadDelPhi[3][numPtBins];

    TH1D* eeInvMass[3][numPtBins];
    TH1D* eeInvMassAll[3];
    TH2F* eePhivMass[3];
    TH2F* eeEtaPhi[3];
    TH2F* eeDcaPt[3];
    TH1D* eeEta[3];
    TH1D* eePhi[3];
    TH1D* eeDca[3];

    TH2F *hPEOrigins[6];

    TH2F* hadEtaPhi;
    TH1D* hadEta;
    TH1D* hadPhi;
    TH1F* hadDca;

    TH2F* elecEtaPhi;
    TH2F* elecDcaPt;
    TH1D* elecDca;
    TH1D* elecEta;
    TH1D* elecPhi;

    TH2F* nSigEPartE[3];
    TH2F* pvePartE[3];
    TH2F* nEtaPartE[3];
    TH2F* nPhiPartE[3];
    TH2F* phiDistPartE[3];
    TH2F* zDistPartE[3];
    TH1F* TPCTracks[3];
    TH1F* EMCMatchedTracks[3];
    TH1F* EMCIdTracks[3];
    TH1F* SMDMatchedTracks[3];
    TH1F* SMDIdTracks[3];
    TH1F* partECutEfficiency[4];
    
    TH1F* refMult;
    TH1F* vertexZ;

    float trigCount[3][numPtBins];
    float eHadNorm;

    void doProjections();
    void prepareCanvas();
    void prepareLabels();
    void getHistograms(TFile*);
    void pretty1DHist(TH1*, int, int);
    void makeUnlikeMinusLikePartnerElectron();

    void drawQAHists();
    void drawEventHists();
    void drawHadQA();
    void drawElecQA();
    void drawEEQA();
    void drawInvMassHists();
    void drawPartECutEffic();
    void overlayEfficiencies();
    void drawPartEQA();
    void drawDeltaPhi(TH1D*, TCanvas*, float, int, TPaveText*);
    void drawCutEfficiencyHists();
    void computeEfficiencyAndDraw(int, TH1F*, TH1F*);
    void makePDF(const char*);

    void setTitleAndAxisLabels(TH1*,TString,TString,TString);
    void setTitleAndAxisLabels(TH2*,TString,TString,TString);

#endif

