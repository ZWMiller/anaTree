#ifndef StMyAnaTreeMaker_h
#define StMyAnaTreeMaker_h
//#include "StRoot/StPicoAnaTreeMaker/StPicoAnaTreeMaker.h"
#include "StMaker.h"
#include "TLorentzVector.h"

#define maxElectrons 50
#define maxMuons 50
#define nVzBin  10
#define nCenBin  9
#define nMagBin  2
#define nEveBin  12
#define nMaxEventsInBuffer 10

class StAnaTree;
class StEventHeader;
class StElectronTrack;
class StMuonTrack;
class StEEPair;
class StEMuPair;
class StMuMuPair;
class StEmcTrigger;
class StPicoAnaTreeMaker;
class StHadronTrack;

class TString;
class TH1F;
class TH2F;
class TH3F;
class TFile;
class TF1;
class TLorentzVector;

class StMyAnaTreeMaker : public StMaker {
	public:
		StMyAnaTreeMaker(const Char_t *name, StPicoAnaTreeMaker *treeMaker, const Char_t *outName, bool mixedEvent);
		virtual ~StMyAnaTreeMaker();

		virtual Int_t Init();
		virtual Int_t Make();
		virtual void  Clear(Option_t *opt="");
		virtual Int_t Finish();

      void addTrigger(int tr, int id) { triggers[id].push_back(tr); };
		void    declareHistograms();
		bool	passHTEIDCuts(StElectronTrack *);
		bool	passEIDCuts(StElectronTrack *);
		bool	passETrackQualityCuts(StElectronTrack *);
		bool	isHTTrigE(StElectronTrack *);
		bool	passMuIDCuts(StMuonTrack *);
		bool	passEEPairCuts(double );
		bool	passEMuPairCuts(double );
		bool	passMuMuPairCuts(double );
      bool  isElectronInValidPair(StElectronTrack*);
		int		getCentrality();
		void	printCuts();
		void	setTrigSelect(int val){ mTrigSelect = val; }
		void  	setVzCut(double min, double max){ mVzCut[0] = min; mVzCut[1] = max; }
		void  	setVzDiffCut(double min, double max){ mVzDiffCut[0] = min; mVzDiffCut[1] = max; }
		Double_t	mtddTCor(double dT, int channel);
		
		void  	setMassBinning(int nBins, double min, double max){ nMassBins = nBins; massMin = min, massMax = max; }
		void  	setPtBinning(int nBins, double min, double max){ nPtBins = nBins; ptMin = min, ptMax = max; }

	private:
		void 	makeEEMixedPairs(int magBufferPointer,int cenBufferPointer, int vzBufferPointer, int eveBufferPointer);
		void 	makeEMuMixedPairs(int magBufferPointer,int cenBufferPointer, int vzBufferPointer, int eveBufferPointer);
		void 	makeMuMuMixedPairs(int magBufferPointer,int cenBufferPointer, int vzBufferPointer, int eveBufferPointer);
		void 	copyCurrentToBuffer(int magBufferPointer,int cenBufferPointer, int vzBufferPointer, int eveBufferPointer);

      void determineTriggers();
      void clearTriggers();
      bool passHadronCuts(StHadronTrack*);
      double delPhiCorrect(double);
      Bool_t checkTriggers(StPicoEvent*, int);
      Bool_t isMB(StPicoEvent*);
      Bool_t isBHT0(StPicoEvent*);
      Bool_t isBHT1(StPicoEvent*);
      Bool_t isBHT2(StPicoEvent*);
      Bool_t isBHT3(StPicoEvent*);
      bool isHT0;
      bool isHT1;
      bool isHT2;
      bool isHT3;
      bool isMB;
      bool makeMixedEvent;
      vector<int> triggers[5]; //0-HT0, 1-HT1 ... 4-MB

      void fillElectronHists(StElectronTrack*);
      void fillHadronHists(StHadronTrack*);
      void fillMuonHists(StMuonTrack*);
      void fillEEHists(StEEPair*);
      void fillEMuHists(StEMuPair*);
      void fillMuMuHists(StMuMuPair*);

		StPicoAnaTreeMaker *mPicoAnaTreeMaker;
		StAnaTree          *mAnaTree;

		TString    mOutName;

		TFile*	   fout;
		TF1 	*fPhiVm;
		TF1   *fEDcaCut;
		TF1	*fMuDcaCut;
		TF1	*fMudTCutLow;
		TF1	*fMudTCutHigh;
		TH1F *hnEvents;
		TH1F *hnTracks;
		TH2F *hVzVpdVz;
		TH2F *hVzdVz;
		TH1F *hRefMultCut;
		TH1F *hVertexZCut;

		TH2F *hNe;
		TH2F *hNemu;
		TH2F *hNmu;
		
		TH1F *hEPt_eff[4];
		TH2F *hEEtavsPhi;
		TH2F *hEPhivsPt;
		TH2F *hEEtavsPt;
		TH2F *hEDcavsPt;
	
		TH2F *hEEtavsPhiwHft;
		TH2F *hEPhivsPtwHft;
		TH2F *hEEtavsPtwHft;
		TH2F *hEDcavsPtwHft;
		TH2F *hEDcaXYvsPtwHft;
		TH2F *hEDcaZvsPtwHft;
		
		TH2F *hPEUSOyOx;
		TH2F *hPEUSOxOz;
		TH2F *hPEUSOrOz;
		TH2F *hPELSOyOx;
		TH2F *hPELSOxOz;
		TH2F *hPELSOrOz;
		
		TH2F *hPEUSOyOxwHft;
		TH2F *hPEUSOxOzwHft;
		TH2F *hPEUSOrOzwHft;
		TH2F *hPELSOyOxwHft;
		TH2F *hPELSOxOzwHft;
		TH2F *hPELSOrOzwHft;

		TH2F *hPEUSDcavsPtwHft;
		TH2F *hPEUSDcaXYvsPtwHft;
		TH2F *hPEUSDcaZvsPtwHft;
		TH2F *hPELSDcavsPtwHft;
		TH2F *hPELSDcaXYvsPtwHft;
		TH2F *hPELSDcaZvsPtwHft;
		
		TH2F *hPEEvPvsPt;
		TH2F *hPEPvEvsPt;
		
		TH2F *hMuEtavsPhi;
		TH2F *hMuPhivsPt;
		TH2F *hMuEtavsPt;
		TH2F *hMuDcavsPt;
		
		TH2F *hMuEtavsPhiwHft;
		TH2F *hMuPhivsPtwHft;
		TH2F *hMuEtavsPtwHft;
		TH2F *hMuDcavsPtwHft;
		TH2F *hMuDcaXYvsPtwHft;
		TH2F *hMuDcaZvsPtwHft;

		TH2F *hEENumInvMassvsPtMB;
		TH2F *hEEDenInvMassvsPtLikePosMB;
		TH2F *hEEDenInvMassvsPtLikeNegMB;
		
		TH2F *hEMuNumInvMassvsPtMB;
		TH2F *hEMuDenInvMassvsPtLikePosMB;
		TH2F *hEMuDenInvMassvsPtLikeNegMB;

		TH2F *hMuMuNumInvMassvsPtMB;
		TH2F *hMuMuDenInvMassvsPtLikePosMB;
		TH2F *hMuMuDenInvMassvsPtLikeNegMB;
		
		TH2F *hEENumInvMassvsPtMBwHft;
		TH2F *hEEDenInvMassvsPtLikePosMBwHft;
		TH2F *hEEDenInvMassvsPtLikeNegMBwHft;
		
		TH2F *hEMuNumInvMassvsPtMBwHft;
		TH2F *hEMuDenInvMassvsPtLikePosMBwHft;
		TH2F *hEMuDenInvMassvsPtLikeNegMBwHft;

		TH2F *hMuMuNumInvMassvsPtMBwHft;
		TH2F *hMuMuDenInvMassvsPtLikePosMBwHft;
		TH2F *hMuMuDenInvMassvsPtLikeNegMBwHft;

		TH2F *hEENumInvMassvsPtMBnophiv;
		TH2F *hEEDenInvMassvsPtLikePosMBnophiv;
		TH2F *hEEDenInvMassvsPtLikeNegMBnophiv;
		
		TH2F *hUSphivM;
		TH2F *hLSPosphivM;
		TH2F *hLSNegphivM;

		TH3F *hEENumInvMassvsPt;
		TH3F *hEEDenInvMassvsPtLikePos;
		TH3F *hEEDenInvMassvsPtLikeNeg;
		
		TH3F *hEMuNumInvMassvsPt;
		TH3F *hEMuDenInvMassvsPtLikePos;
		TH3F *hEMuDenInvMassvsPtLikeNeg;
		
		TH3F *hMuMuNumInvMassvsPt;
		TH3F *hMuMuDenInvMassvsPtLikePos;
		TH3F *hMuMuDenInvMassvsPtLikeNeg;

		TH3F *hEENumInvMassvsPtwHft;
		TH3F *hEEDenInvMassvsPtLikePoswHft;
		TH3F *hEEDenInvMassvsPtLikeNegwHft;
		
		TH3F *hEMuNumInvMassvsPtwHft;
		TH3F *hEMuDenInvMassvsPtLikePoswHft;
		TH3F *hEMuDenInvMassvsPtLikeNegwHft;
		
		TH3F *hMuMuNumInvMassvsPtwHft;
		TH3F *hMuMuDenInvMassvsPtLikePoswHft;
		TH3F *hMuMuDenInvMassvsPtLikeNegwHft;

		
		TH3F *hEEDenInvMassvsPtMix;
		TH3F *hEEDenInvMassvsPtMixLikePos;
		TH3F *hEEDenInvMassvsPtMixLikeNeg;
		
		TH3F *hEMuDenInvMassvsPtMix;
		TH3F *hEMuDenInvMassvsPtMixLikePos;
		TH3F *hEMuDenInvMassvsPtMixLikeNeg;

		TH3F *hMuMuDenInvMassvsPtMix;
		TH3F *hMuMuDenInvMassvsPtMixLikePos;
		TH3F *hMuMuDenInvMassvsPtMixLikeNeg;

		TH3F *hEMuDenInvMassvsPtMixwHft;
		TH3F *hEMuDenInvMassvsPtMixLikePoswHft;
		TH3F *hEMuDenInvMassvsPtMixLikeNegwHft;
		
		TH3F *hEEDenInvMassvsPtMixwHft;
		TH3F *hEEDenInvMassvsPtMixLikePoswHft;
		TH3F *hEEDenInvMassvsPtMixLikeNegwHft;

		TH3F *hMuMuDenInvMassvsPtMixwHft;
		TH3F *hMuMuDenInvMassvsPtMixLikePoswHft;
		TH3F *hMuMuDenInvMassvsPtMixLikeNegwHft;

      //hadron branches
      TH1F *hHadPt;
      TH1F *hHadDca;
      TH1F *hHadDcaXY;
      TH1F *hHadDcaZ;
      TH1F *hHadPhi;
      TH1F *hEPt;
      TH1F *hEEPt_US;
      TH1F *hEEPt_LS;

      TH2F *hHadEtaPhi;
      TH2F *hHadEDelPhiPt;
      TH2F *hHadEDelPhiPt_high;
      TH2F *hHadMuDelPhiPt;
      TH2F *hHadDcaPt;
      TH2F *hHadEEDelPhiPt_LS;
      TH2F *hHadEEDelPhiPt_US;
      TH2F *hHadPtEPt;

		//TH2F *hDenInvMassvsPtMixLikePosMB;
		//TH2F *hDenInvMassvsPtMixLikeNegMB;
		//TH2F *hDenInvMassvsPtMixLikePosMBnophiv;
		//TH2F *hDenInvMassvsPtMixLikeNegMBnophiv;
		//TH2F *hMixphivM;
		//TH2F *hMixLikePosphivM;
		//TH2F *hMixLikeNegphivM;
		//TH3F *hDenInvMassvsPtMixLikePos;
		//TH3F *hDenInvMassvsPtMixLikeNeg;
		//

		Int_t mNBadRuns;
		Float_t     mHTth;
		Float_t     mHTAdc0th;
		Float_t     mEmcPtth;
		Int_t	    mTrigSelect; //-1 - all, 0 - MB, 1 - HT0, 2 - HT1, 3 - HT2, 4 - HT3, 5 - EMu, 6 - dimuon..

		Double_t	mMtdT0Corr[30][5];
		Float_t     mVzCut[2];
        Float_t     mVzDiffCut[2];

        Int_t       mnHitsFitCut[2];
        Int_t       mnHitsDedxCut[2];
        Float_t     mRatioCut[2];
        Bool_t      mHFTTrackCut;

        Float_t     mEPtCut[2];
        Float_t     mEEtaCut[2];
        Float_t     mEDcaCut[2];
        Float_t     mEInvBetaCut[2];
        Float_t     mELocalYCut[2];
        Float_t     mELocalZCut[2];
        Float_t     mEnSigECut[2];
        Float_t     mHTEnSigECut[2];
        
        Float_t     mEmcEPtCut[2];
        Float_t     mEmcEEtaCut[2];
        Float_t     mEmcEPveCut[2];
        Float_t     mEmcEDcaCut[2];
        
        Int_t      mEnEtaCut[2];
        Int_t      mEnPhiCut[2];
        Float_t     mEZDistCut[2];
        Float_t     mEPhiDistCut[2];

        Float_t     mMuPtCut[2];
        Float_t     mMuEtaCut[2];
        Float_t     mMuDcaCut[2];
        Float_t     mMunSigPiCut[2];
        Float_t     mMudTCut[2];
        Float_t     mMudZCut[2];
        Float_t     mMudYCut[2];

        Float_t     mDauEPtCut[2];
        Float_t     mDauEEtaCut[2];
        Float_t     mDauEDcaToVtxCut[2]; 
        Float_t     mDauEDcaDistCut[2]; //! DCA between two daughters
        
        Float_t     mDauMuPtCut[2];
        Float_t     mDauMuEtaCut[2];
        Float_t     mDauMuDcaToVtxCut[2];
        
        Float_t     mCosThetaStarCut[2]; //! cos(theta*)
        Float_t     mPointingAngleCut[2]; //! 

        Float_t     mPairDcaCut[2];
        Float_t     mPairDecayLCut[2]; 
        Float_t     mPairYCut[2];
        Float_t     mPairMassCut[2];

        Float_t     mEEYCut[2];
        Float_t     mEMuYCut[2];
        Float_t     mMuMuYCut[2];

        Float_t     mHadPtCut[2];
        Float_t     mHadEtaCut[2];
        Float_t     mHadDcaCut[2];
        
		Float_t     mPEMassCut[2];
		Short_t	    nMassBins;
		Short_t	    nPtBins;
		Float_t     massMin;
		Float_t     massMax;
		Float_t     ptMin;
		Float_t     ptMax;
      int         centrality;
		
		UChar_t	    current_nePlus;
		UChar_t	    current_neMinus;
		UChar_t	    current_nmuPlus;
		UChar_t	    current_nmuMinus;

		TLorentzVector  current_ePlus[maxElectrons];
		TLorentzVector  current_eMinus[maxElectrons];
		UChar_t 		current_ePlusFlag[maxElectrons]; //1 TPC+TOF only; 2 EMConly; 3 both; >4 EMC triggered
		UChar_t 		current_eMinusFlag[maxElectrons];
		UChar_t 		current_ePlusIsHFT[maxElectrons]; //0 not HFT, 1is HFT 
		UChar_t 		current_eMinusIsHFT[maxElectrons];

		TLorentzVector current_muPlus[maxMuons];
		TLorentzVector current_muMinus[maxMuons];
		UChar_t 		current_muPlusFlag[maxElectrons]; //0 MTD, 1 MTD triggered
		UChar_t 		current_muMinusFlag[maxElectrons];
		UChar_t 		current_muPlusIsHFT[maxElectrons]; //0 not HFT, 1is HFT 
		UChar_t 		current_muMinusIsHFT[maxElectrons];


		UChar_t 	nEventsInBuffer[nMagBin][nCenBin][nVzBin][nEveBin];
		Bool_t 		buffer_fullFlag[nMagBin][nCenBin][nVzBin][nEveBin];
		UChar_t 		buffer_nePlus[nMagBin][nCenBin][nVzBin][nEveBin][nMaxEventsInBuffer];
		UChar_t 		buffer_neMinus[nMagBin][nCenBin][nVzBin][nEveBin][nMaxEventsInBuffer];
		UChar_t 		buffer_ePlusFlag[nMagBin][nCenBin][nVzBin][nEveBin][nMaxEventsInBuffer][maxElectrons];
		UChar_t 		buffer_eMinusFlag[nMagBin][nCenBin][nVzBin][nEveBin][nMaxEventsInBuffer][maxElectrons];
		UChar_t 		buffer_ePlusIsHFT[nMagBin][nCenBin][nVzBin][nEveBin][nMaxEventsInBuffer][maxElectrons];
		UChar_t 		buffer_eMinusIsHFT[nMagBin][nCenBin][nVzBin][nEveBin][nMaxEventsInBuffer][maxElectrons];

		TLorentzVector buffer_ePlus[nMagBin][nCenBin][nVzBin][nEveBin][nMaxEventsInBuffer][maxElectrons];
		TLorentzVector buffer_eMinus[nMagBin][nCenBin][nVzBin][nEveBin][nMaxEventsInBuffer][maxElectrons];

		UChar_t 		buffer_nmuPlus[nMagBin][nCenBin][nVzBin][nEveBin][nMaxEventsInBuffer];
		UChar_t 		buffer_nmuMinus[nMagBin][nCenBin][nVzBin][nEveBin][nMaxEventsInBuffer];
		UChar_t 		buffer_muPlusFlag[nMagBin][nCenBin][nVzBin][nEveBin][nMaxEventsInBuffer][maxElectrons];
		UChar_t 		buffer_muMinusFlag[nMagBin][nCenBin][nVzBin][nEveBin][nMaxEventsInBuffer][maxElectrons];
		UChar_t 		buffer_muPlusIsHFT[nMagBin][nCenBin][nVzBin][nEveBin][nMaxEventsInBuffer][maxElectrons];
		UChar_t 		buffer_muMinusIsHFT[nMagBin][nCenBin][nVzBin][nEveBin][nMaxEventsInBuffer][maxElectrons];

		TLorentzVector buffer_muPlus[nMagBin][nCenBin][nVzBin][nEveBin][nMaxEventsInBuffer][maxElectrons];
		TLorentzVector buffer_muMinus[nMagBin][nCenBin][nVzBin][nEveBin][nMaxEventsInBuffer][maxElectrons];
		Int_t     	iran;
		Int_t 		current_centrality;

		ClassDef(StMyAnaTreeMaker, 1)
};

#endif
