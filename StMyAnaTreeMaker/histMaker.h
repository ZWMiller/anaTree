#ifndef _histMaker
#define _histMaker

    const int numPtBins = 17;
    const int numCanvas = numPtBins/9 + 1;

    const float lowpt[numPtBins]  = 
     {1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.25, 2.5, 3.0, 3.5, 4.0, 
      5.0, 6.0, 8.0, 10.0};
    const float highpt[numPtBins] = 
     {1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.25, 2.5, 3.0, 3.5, 4.0, 5.0,
      6.0, 8.0, 10., 20.0 };
    const int colors[5] = {kBlack,kRed,kAzure+1,kMagenta,kYellow};


    const int numPtBinsEFF = 26;

    const float lowptEFF[numPtBinsEFF]  = 
     {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 
      1.9, 2.0, 2.25, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 13.0};

    TCanvas* dPhiPt[3][numCanvas];
    TCanvas* invMassPt[numCanvas];
    TCanvas* invMassVsPt;
    TCanvas* nSigEPt[numCanvas];
    TCanvas* invMass;
    TCanvas* ptCompare;
    TCanvas* cutEfficiency;
    TCanvas* hadQA;
    TCanvas* eeQA[2];
    TCanvas* eeQAPhiv;
    TCanvas* eeOriginQA;
    TCanvas* elecQA;
    TCanvas* pElecCuts[6];
    TCanvas* eIDCutEffic[4];
    TCanvas* efficOverlay;
    TCanvas* eventHists;
    TCanvas* nSigEff[numCanvas];
    TCanvas* twoGaus[numCanvas];
    TCanvas* nSigCutPlot;
    TCanvas* nSigMeanSig;
    TCanvas* dndpt;
    TCanvas* zdcQA;

    TPaveText* lbl[numPtBins];
    TPaveText* sampleLabel;

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

    TH1D* elecDcaForInt[numPtBins];
    TH1D* eeDcaForInt[3][numPtBins];

    TH2F* nSigEPartE[3];
    TH1D* nSigE[3][numPtBins];
    TH1D* nSigEForInt[3][numPtBins];
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

    TH1F* nSigEeff[numPtBins];
    TF2 *twogaus[numPtBins]; 
    
    TH1F* refMult;
    TH1F* vertexZ;
    TH2F* refMultZDCvsRunIndex;
    TH2F* refMultvsZDCx;

    int trigSelect;
    char trigLabel[100];
    float trigCount[3][numPtBins];
    float eHadNorm;
    double *nSigCovariance[numPtBins];
    double meanSigE[numPtBins];
    double meanerrSigE[numPtBins];
    double sigmaSigE[numPtBins];
    double sigmaerrSigE[numPtBins];
    double corSigE[numPtBins];

    void doProjections();
    void prepareCanvas();
    void prepareLabels();
    TString getTriggerLabel();
    void declareHistograms();
    bool checkIfFileOpen(TFile*);
    void getHistograms(TFile*);
    void pretty1DHist(TH1*, int, int);
    void pretty1DHistFill(TH1*, int, int);
    void prettyTGraph(TGraphErrors*, int, int, float, float);
    void prettyTGraph(TGraph*, int, int, float, float);
    void makeUnlikeMinusLikePartnerElectrons();

    void drawQAHists();
    void getnSigEeff();
    void drawEventHists();
    void drawHadQA();
    void drawElecQA();
    void drawEEQA();
    void drawInvMassHists();
    void drawnSigE();
    void drawPartECutEffic();
    void overlayEfficiencies();
    void drawPartEQA();
    void getdNdpT();
    void drawDeltaPhi(TH1D*, TCanvas*, float, int, TPaveText*);
    void drawCutEfficiencyHists();
    void computeEfficiencyAndDraw(int, TH1F*, TH1F*);
    void makePDF(const char*);
    void drawnSigMeanSig( const double*, const double*, const double*, const double*, const double*, const double*);

    void setTitleAndAxisLabels(TH1*,TString,TString,TString);
    void setTitleAndAxisLabels(TH2*,TString,TString,TString);


#endif

