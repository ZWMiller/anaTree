#ifndef StPicoAnaTreeMaker_h
#define StPicoAnaTreeMaker_h
#include "StRoot/StPicoDstMaker/StPicoTrack.h"
#include "StMaker.h"
#include <vector>
#include <utility>
#include <string>
#include "TRandom3.h"
class StPicoDst;
class StPicoEvent;
class StPicoDstMaker;
class StPicoMtdHit;
class StPicoTrack;
class StPicoEmcTrigger;
class StPicoMtdTrigger;
class TString;
class TH1F;
class TH2F;
class TH3F;
class TTree;
class TFile;
class TProfile2D;
class TChain;
class StAnaTree;
class StEventHeader;
class StElectronTrack;
class StPartElectronTrack;
class StMuonTrack;
class StHadronTrack;
class StEEPair;
class StPhoEEPair;
class StEMuPair;
class StMuMuPair;
class StEmcTrigger;
class StMtdTrigger;
class StEmcGeom;
class StRefMultCorr;
class CentralityMaker;
#include "StAnaTreeArrays.h"

#if !defined(ST_NO_NAMESPACES)
using namespace std;
#endif

enum triggerType{
	minbias=0,
	ht,
	mtd	
};

class StPicoAnaTreeMaker : public StMaker {
	public:
		StPicoAnaTreeMaker(Int_t mode, const char* fileName, StPicoDstMaker *picoMaker, const char *name="anaTreeMaker");
		virtual ~StPicoAnaTreeMaker();

		virtual Int_t Init();
		virtual Int_t Make();
		virtual void  Clear(Option_t *opt="");
		virtual Int_t Finish();

		void SetStatus(const char *arrType, int status);

		void setSplit(int);
		void setCompression(int);
		void setBufferSize(int);
		void setRunNumber(int);
		void setSaveHadron(bool);
		void addTrigger(int);

		void assignArrays();
		void streamerOff();
		void clearArrays();
		void zeroArrays();
		void createArrays();

		void setBranchAddresses();
		void setBranchAddresses(TChain*);

		void openRead();
		void openWrite();
		void closeWrite(); 
		void closeRead();
		Int_t MakeRead();
		Int_t MakeWrite();
		bool passEvent(StPicoDst *);
		void fillRecenterCor();
		void fillEventHeader();
		void fillTracks();
		void fillPairs();
		void fillEmcTrigger();
		void fillMtdTrigger();

		const float *getRecenterCor(int runId, int centrality);

		TClonesArray* clonesArray(TClonesArray*& p, const char* type, int size, int& counter);

		void setTriggerSelection(triggerType val){mTriggerSelection=val;}

		void setVzCut(const Float_t min, const Float_t max){mVzCut[0] = min; mVzCut[1] = max; }
		void setVzDiffCut(const Float_t min, const Float_t max){mVzDiffCut[0] = min; mVzDiffCut[1] = max; }
		void setPtCut(const Float_t min, const Float_t max){mPtCut[0] = min; mPtCut[1] = max; }
		void setEtaCut(const Float_t min, const Float_t max){mEtaCut[0] = min; mEtaCut[1] = max; }
		void setDcaCut(const Float_t min, const Float_t max){mDcaCut[0] = min; mDcaCut[1] = max; }
		void setnHitsFitCut(const Int_t min, const Int_t max){mnHitsFitCut[0] = min; mnHitsFitCut[1] = max; }
		void setnHitsDedxCut(const Int_t min, const Int_t max){mnHitsDedxCut[0] = min; mnHitsDedxCut[1] = max; }
		void setRatioCut(const Float_t min, const Float_t max){mRatioCut[0] = min; mRatioCut[1] = max; }
//		void setEPtCut(const Float_t min, const Float_t max){mEPtCut[0] = min; mEPtCut[1] = max; }
		void setEDcaCut(const Float_t min, const Float_t max){mEDcaCut[0] = min; mEDcaCut[1] = max; }
		void setEInvBetaCut(const Float_t min, const Float_t max){mEInvBetaCut[0] = min; mEInvBetaCut[1] = max; }
//		void setELocalYCut(const Float_t min, const Float_t max){mELocalYCut[0] = min; mELocalYCut[1] = max; }
//		void setELocalZCut(const Float_t min, const Float_t max){mELocalZCut[0] = min; mELocalZCut[1] = max; }
		void setEnSigECut(const Float_t min, const Float_t max){mEnSigECut[0] = min; mEnSigECut[1] = max; }
		
      void setPartEnSigECut(const Float_t min, const Float_t max){mPartEnSigECut[0] = min; mPartEnSigECut[1] = max; }
      void setPhoEPairDcaCut(const Float_t cut){mPhoEPairDcaCut = cut;}
      void setPhoEPairMassCut(const Float_t cut){mPhoEMassCut = cut;}

		void setEmcEPtCut(const Float_t min, const Float_t max){mEmcEPtCut[0] = min; mEmcEPtCut[1] = max; }
		void setEmcEEtaCut(const Float_t min, const Float_t max){mEmcEEtaCut[0] = min; mEmcEEtaCut[1] = max; }
		void setEmcEPveCut(const Float_t min, const Float_t max){mEmcEPveCut[0] = min; mEmcEPveCut[1] = max; }

		void setEnEtaCut(const Char_t min, const Char_t max){mEnEtaCut[0] = min; mEnEtaCut[1] = max; }
		void setEnPhiCut(const Char_t min, const Char_t max){mEnPhiCut[0] = min; mEnPhiCut[1] = max; }
		void setEZDistCut(const Float_t min, const Float_t max){mEZDistCut[0] = min; mEZDistCut[1] = max; }
		void setEPhiDistCut(const Float_t min, const Float_t max){mEPhiDistCut[0] = min; mEPhiDistCut[1] = max; }

		void setMuPtCut(const Float_t min, const Float_t max){mMuPtCut[0] = min; mMuPtCut[1] = max; }
		void setMuEtaCut(const Float_t min, const Float_t max){mMuEtaCut[0] = min; mMuEtaCut[1] = max; }
		void setMunSigPiCut(const Float_t min, const Float_t max){mMunSigPiCut[0] = min; mMunSigPiCut[1] = max; }
		void setMudTCut(const Float_t min, const Float_t max){mMudTCut[0] = min; mMudTCut[1] = max; }
		void setMudZCut(const Float_t min, const Float_t max){mMudZCut[0] = min; mMudZCut[1] = max; }
		void setMudYCut(const Float_t min, const Float_t max){mMudYCut[0] = min; mMudYCut[1] = max; }

		void setDauEPtCut(const Float_t min, const Float_t max){mDauEPtCut[0] = min; mDauEPtCut[1] = max; }
		void setDauEDcaToVtxCut(const Float_t min, const Float_t max){mDauEDcaToVtxCut[0] = min; mDauEDcaToVtxCut[1] = max; }
		void setDauEDcaDistCut(const Float_t min, const Float_t max){mDauEDcaDistCut[0] = min; mDauEDcaDistCut[1] = max; }

		void setDauMuPtCut(const Float_t min, const Float_t max){mDauMuPtCut[0] = min; mDauMuPtCut[1] = max; }
		void setDauMuEtaCut(const Float_t min, const Float_t max){mDauMuEtaCut[0] = min; mDauMuEtaCut[1] = max; }
		void setDauMuDcaToVtxCut(const Float_t min, const Float_t max){mDauMuDcaToVtxCut[0] = min; mDauMuDcaToVtxCut[1] = max; }

		void setCosThetaStarCut(const Float_t min, const Float_t max){mCosThetaStarCut[0] = min; mCosThetaStarCut[1] = max; }
		void setPointingAngleCut(const Float_t min, const Float_t max){mPointingAngleCut[0] = min; mPointingAngleCut[1] = max; }

		void setPairDcaCut(const Float_t min, const Float_t max){mPairDcaCut[0] = min; mPairDcaCut[1] = max; }
		void setPairDecayLCut(const Float_t min, const Float_t max){mPairDecayLCut[0] = min; mPairDecayLCut[1] = max; }
		void setPairYCut(const Float_t min, const Float_t max){mPairYCut[0] = min; mPairYCut[1] = max; }
		void setPairMassCut(const Float_t min, const Float_t max){mPairMassCut[0] = min; mPairMassCut[1] = max; }
		void setMaxRunId(int val);
		void setMaxCentrality(int val);
		void setDoCalcRecenter(Bool_t doIt){mCalcRecenter = doIt;}
		void setDoEvtPlane(Bool_t doIt){mDoEvtPlane = doIt;}
		void setInputRunList(const char *input){ mRunList = input; }
		void setInputRecenterFile(const char *input){ mRecenterFile = input; }
      void setPicoMode(int val) { mPicoMode = val;}

		Bool_t isGoodTrack(StPicoTrack* ); 
		Bool_t isTofElectron(StPicoTrack* );
		Bool_t isEmcElectron(StPicoTrack* );
		Bool_t isPartE(StPicoTrack* );
		Bool_t isMuon(StPicoTrack* );
		Bool_t isHadron(StPicoTrack* );
		Bool_t passEEPair(StElectronTrack *, StElectronTrack *, Int_t , Int_t);
		Bool_t passPhoEEPair(StPicoTrack *, StElectronTrack *, Int_t , Int_t);
		Bool_t passEMuPair(StElectronTrack *, StMuonTrack *, Int_t, Int_t);
		Bool_t passMuMuPair(StMuonTrack *, StMuonTrack *, Int_t, Int_t);

		void    declareHistos();
		void    writeHistos();

		StAnaTree* anaTree();
		TChain*  chain();
		TTree*   tree();
		void    printCuts();
		Int_t makeTriggerWord(StPicoEvent* ); 
      void printTriggerWords();

		enum ioMode {ioRead, ioWrite};
	private:

		TH1F *mhnEvents;
		TH1F *mhnTracks;
		TH2F *mhnSigEvsP;
		TH2F *mhnSigEvsPt;
		TH2F *mhnSigPivsP;
		TH2F *mhnSigPivsPt;

		TH2F *mhTofEnSigEvsPCut;
		TH2F *mhTofEnSigEvsPtCut;
		TH2F *mhEmcEnSigEvsPCut;
		TH2F *mhEmcEnSigEvsPtCut;
		TH2F *mhSmdEnSigEvsPCut;
		TH2F *mhSmdEnSigEvsPtCut;

		TH2F *mhInvBetavsP;
		TH2F *mhInvBetavsPt;

		TH2F *mhM2vsP;
		TH2F *mhM2vsPt;

		TH2F *mhTofLocalYvsTray;
		TH2F *mhTofLocalZvsTray;

		TH2F *mhnEtavsnPhi;
		TH2F *mhZDistvsPt;
		TH2F *mhPhiDistvsPt;

		TH2F *mhEvPvsPt;
		TH2F *mhPvEvsPt;

		TH2F *mhMtdMunSigPivsP;
		TH2F *mhMtdMunSigPivsPt;
		TH2F *mhMtdMunSigPivsPCut;
		TH2F *mhMtdMunSigPivsPtCut;

		TH2F *mhMtddZvsPt;
		TH2F *mhMtddYvsPt;
		TH2F *mhMtddTvsPt;

		TH2F *mhMtddZvsMod;
		TH2F *mhMtddYvsMod;
		TH2F *mhMtddTvsMod;

		TH2F *mhnSigE2PivsP;
		TH2F *mhnSigEPivsP;
		TH2F *mhnSigEKvsP;
		TH2F *mhnSigEPvsP;

		TH2F *mhnSigE2PivsPt;
		TH2F *mhnSigEPivsPt;
		TH2F *mhnSigEKvsPt;
		TH2F *mhnSigEPvsPt;

		TH2F *mhnSigEPivsPwTOF;
		TH2F *mhnSigEKvsPwTOF;
		TH2F *mhnSigEPvsPwTOF;

		TH2F *mhnSigEPivsPtwTOF;
		TH2F *mhnSigEKvsPtwTOF;
		TH2F *mhnSigEPvsPtwTOF;

		//dsmADC vs ADC0, dsmADC vs p, dsmADC vs E, ADC0 vs p, ADC0 vs E
		TH2F *mhDsmAdcvsAdc0;
		TH2F *mhDsmAdcvsP;
		TH2F *mhDsmAdcvsE;
		TH2F *mhAdc0vsP;
		TH2F *mhAdc0vsE;

		TProfile2D *cosfarwest_correction;
		TProfile2D *sinfarwest_correction;
		TProfile2D *coswest_correction;
		TProfile2D *sinwest_correction;
		TProfile2D *coseast_correction;
		TProfile2D *sineast_correction;
		TProfile2D *cosfareast_correction;
		TProfile2D *sinfareast_correction;
		TProfile2D*	etapluszplusQx; 
		TProfile2D*	etapluszplusQy;
		TProfile2D*	etapluszminusQx; 
		TProfile2D*	etapluszminusQy;
		TProfile2D*	etaminuszplusQx; 
		TProfile2D*	etaminuszplusQy;
		TProfile2D*	etaminuszminusQx; 
		TProfile2D*	etaminuszminusQy;



		StPicoDstMaker *mPicoDstMaker;
		StPicoDst      *mPicoDst;
		StAnaTree      *mAnaTree;
		StEmcGeom	   *mEmcGeom;

		Int_t      mIoMode;         //! I/O mode:  0: - read,   1: - write
		TFile*    mOutputFile;
		TFile*    mOutputHistFile;
		TChain*   mChain;
		TTree*    mTTree;

		TString   mInputFileName;        //! *.list - MuDst or picoDst
		TString   mOutputFileName;       //! FileName

		Int_t     mRunNumber;
		Int_t     mEventCounter;
		Int_t     mSplit;
		Int_t     mCompression;
		Int_t     mBufferSize;
		Bool_t    mSaveHadron;
		Int_t	    mNMaxRunId;
		Int_t	    mNMaxCentrality;
		Int_t 	 mTriggerSelection; 	//! 0 = minbias; 1 = ht; 2 = st_mtd;
      Int_t     mPicoMode;          //! 0 = Run14 pico, 1 = Run15 pico;
      Int_t     mTriggerWord;

		Bool_t	 mCalcRecenter; //! calculate Recenter 
		Bool_t	 mDoEvtPlane; //! do Eventplane 
		TString	 mRunList;
		TString	 mRecenterFile;


		Float_t     mVzCut[2];
		Float_t     mVzDiffCut[2];

		Float_t     mPtCut[2];
		Float_t     mEtaCut[2];
		Float_t     mDcaCut[2];
		Int_t       mnHitsFitCut[2];
		Int_t       mnHitsDedxCut[2];
		Float_t     mRatioCut[2];

		//Float_t     mEPtCut[2];
		Float_t     mEDcaCut[2];
		Float_t     mEInvBetaCut[2];
		//Float_t     mELocalYCut[2];
		//Float_t     mELocalZCut[2];
		Float_t     mEnSigECut[2];
		
      Float_t     mPartEnSigECut[2];
      Float_t     mPhoEPairDcaCut;
      Float_t     mPhoEMassCut;

		Float_t     mEmcEPtCut[2];
		Float_t     mEmcEEtaCut[2];
		Float_t     mEmcEPveCut[2];

		Char_t      mEnEtaCut[2];
		Char_t      mEnPhiCut[2];
		Float_t     mEZDistCut[2];
		Float_t     mEPhiDistCut[2];

		Float_t     mMuPtCut[2];
		Float_t     mMuEtaCut[2];
		Float_t     mMunSigPiCut[2];
		Float_t     mMudTCut[2];
		Float_t     mMudZCut[2];
		Float_t     mMudYCut[2];

		Float_t     mDauEPtCut[2];
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

		Float_t     mRecenterCor[8];
      vector<Int_t> partEidx;
		Float_t		mSizeAll;
		Float_t		mSizeBranch[__NANATREEARRAYS__];

		TRandom3       *mRandom;
	protected:
		friend class StAnaTree;

		TClonesArray*   mAnaTreeAllArrays[__NANATREEARRAYS__];   
      TClonesArray**  mAnaTreeArrays;   //[__NANATREEARRAYS__]
      char            mStatusArrays[__NANATREEARRAYS__];
      std::vector<int> triggers;   // vpdmb-5-p-nobsmd 

		ClassDef(StPicoAnaTreeMaker, 1)
};

#endif
inline StAnaTree* StPicoAnaTreeMaker::anaTree() { return mAnaTree; }
inline TChain* StPicoAnaTreeMaker::chain() { return mChain; }
inline TTree* StPicoAnaTreeMaker::tree() { return mTTree; }
inline void StPicoAnaTreeMaker::setSplit(int split) { mSplit = split;}
inline void StPicoAnaTreeMaker::setCompression(int comp) { mCompression = comp;}
inline void StPicoAnaTreeMaker::setBufferSize(int buf) { mBufferSize = buf; }
inline void StPicoAnaTreeMaker::setRunNumber(int run) { mRunNumber = run; }
inline void StPicoAnaTreeMaker::setSaveHadron(bool val) { mSaveHadron = val; }
inline void StPicoAnaTreeMaker::setMaxRunId(int val) { mNMaxRunId = val; }
inline void StPicoAnaTreeMaker::setMaxCentrality(int val) { mNMaxCentrality = val; }
inline void StPicoAnaTreeMaker::addTrigger(int val) { triggers.push_back((unsigned int)val); }
