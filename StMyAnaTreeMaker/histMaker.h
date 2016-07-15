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
    const int colors[5] = {kBlack,kRed,kBlue,kMagenta,kYellow};


    TCanvas* dPhiPt[3][numCanvas];
    TCanvas* invMassPt[numCanvas];
    TCanvas* invMass;
    TCanvas* ptCompare;
    TCanvas* cutEfficiency;
    TPaveText* lbl[numPtBins];

    TH2F* eHadDelPhiPt[3];
    TH2F* eeInvMassPt[3];
    TH2F* hadPtEPt;
    TH1F* ePt[3];
    TH1F* ePt_eff[4];
    TH1D* eHadDelPhi[3][numPtBins];
    TH1D* eeInvMass[3][numPtBins];
    TH1D* eeInvMassAll[3];

    float trigCount[3][numPtBins];

    void doProjections();
    void prepareCanvas();
    void prepareLabels();
    void getHistograms(TFile*);

    void drawQAHists();
    void drawInvMassHists();
    void drawDeltaPhi(TH1D*, TCanvas*, float, int, TPaveText*);
    void drawCutEfficiencyHists();
    void computeEfficiencyAndDraw(int, TH1F*, TH1F*);
    void makePDF(const char*);

#endif

