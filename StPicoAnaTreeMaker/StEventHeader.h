#ifndef StEventHeader_hh
#define StEventHeader_hh

class StMuDst;
class StPicoDst;
class StPicoEvent;
class TClonesArray;
class StMuPrimaryVertex;
class StPicoAnaTree;
class StPicoTrack;
class StPicoAnaTreeMaker;
class StBTofHeader;
#include "StThreeVectorF.hh"
#include "TVector2.h"
#include "TRandom3.h"

class StEventHeader : public TObject {
	public:
		StEventHeader();
		~StEventHeader();
		StEventHeader(const StPicoDst& picoDst) ;
		StEventHeader(const StPicoDst& picoDst, const int trigWord) ;
      void prepareEventInfo(const StPicoDst& picoDst) ; // shared by different constructors
		void Clear(const Option_t*) {}

		Int_t    runId() const            { return mRunId; }
		Int_t    eventId() const          { return mEventId; }
		Float_t  bField() const           { return mBField; }
		StThreeVectorF primaryVertex() const { return StThreeVectorF(mVx,mVy,mVz); }
		Int_t    triggerWord() const      { return mTriggerWord; }
		Int_t    triggerWordMtd() const   { return mTriggerWordMtd; }
		Int_t    refMultPos() const       { return (Int_t)mRefMult-mRefMultNeg; }
		Int_t    refMultNeg() const       { return (Int_t)mRefMultNeg; }
		Int_t    refMult() const          { return (Int_t)(mRefMult); }
		Int_t    grefMult() const           { return (Int_t)(mGRefMult); }
	
	
		Float_t     grefMultCorr() const      { return (Float_t)mGRefMultCorr;}
		Float_t     ReweightGRefMult() const  { return (Float_t)mReweight_GRefMult;}
		UShort_t        Centrality9() const       { return (UShort_t)mCentrality9;}
		UShort_t        Centrality16() const      { return (UShort_t)mCentrality16;}


		// event plane
		Float_t	eventplane0() const	{return (Float_t) mEventplane0/100.;}
		Float_t	eventPlane() const	{return (Float_t) mEventplane/100.;}
		Float_t	eventPlane1() const	{return (Float_t) mEventplane1/100.;}
		Float_t	eventPlane2() const	{return (Float_t) mEventplane2/100.;}

		Float_t   Q0x() const {return (Float_t) mQ0x/100.;}
		Float_t   Q0y() const {return (Float_t) mQ0y/100.;}
	
		Float_t   Qx() const {return (Float_t) mQx/100.;}
		Float_t   Qy() const {return (Float_t) mQy/100.;}

		Float_t   Q1x() const {return (Float_t) mQ1x/100.;}
		Float_t   Q1y() const {return (Float_t) mQ1y/100.;}
		Float_t   Q2x() const {return (Float_t) mQ2x/100.;}
		Float_t   Q2y() const {return (Float_t) mQ2y/100.;}
		Float_t   Qplusx()   const {return (Float_t) mQplusx/100.;}
		Float_t   Qplusy()   const {return (Float_t) mQplusy/100.;}
		Float_t   Qminusx()  const {return (Float_t) mQminusx/100.;}
		Float_t   Qminusy()  const {return (Float_t) mQminusy/100.;}

		Float_t   Qetaplusx() const {return (Float_t) mQetaplusx/100.;}
		Float_t   Qetaplusy() const {return (Float_t) mQetaplusy/100.;}
		Float_t   Qetaminusx() const {return (Float_t) mQetaminusx/100.;}
		Float_t   Qetaminusy() const {return (Float_t) mQetaminusy/100.;}

		Int_t    numberOfPxlInnerHits() const { return (Int_t)(mNHitsHFT[0]); }
		Int_t    numberOfPxlOuterHits() const { return (Int_t)(mNHitsHFT[1]); }
		Int_t    numberOfIstHits() const      { return (Int_t)(mNHitsHFT[2]); }
		Int_t    numberOfSsdHits() const      { return (Int_t)(mNHitsHFT[3]); }

		Float_t  vzVpd() const            { return mVzVpd; }

		Float_t  ZDCx() const             { return (Float_t)mZDCx; }
		Float_t  BBCx() const             { return (Float_t)mBBCx; }

		Float_t backgroundRate() const             { return mBackgroundRate; }

		UShort_t btofTrayMultiplicity() const { return mbTofTrayMultiplicity ; }
		UShort_t numberOfGlobalTracks() const { return mNumberOfGlobalTracks ; }

		Float_t ranking() const { return mRanking ; }
		UShort_t nBEMCMatch() const { return mNBEMCMatch ; }
		UShort_t nBTOFMatch() const { return mNBTOFMatch ; }

		Int_t   ht_th(const Int_t i) { return mHT_Th[i]; }

		// other user's functions
		int      year() const;
		int      day() const;
		float    energy() const;
		bool     isMinBias() const;
		bool     isMBSlow() const;
		bool     isCentral() const;
		bool     isHT() const;    
		bool     isHT11() const; 
		bool     isHT15() const; 
		bool     isHT18() const;
		bool     isHT25() const;
		bool     isMtdTrig() const;
		bool     isDiMuon() const;
		bool     isDiMuonHFT() const;
		bool     isSingleMuon() const;
		bool     isEMuon() const;

		// set functions for trigger thresholds
		void     setHT_Th(const Int_t i, const Int_t th) { mHT_Th[i] = (UChar_t)th; }
		
	protected: //these are written out

		Int_t          mRunId;           // run number
		Int_t          mEventId;         // event number
		Float_t        mBField;          // B field in kilogauss 
		Float_t        mVx;              // primary vertex x
		Float_t        mVy;              // primary vertex y
		Float_t        mVz;              // primary vertex z
		Float_t        mVzVpd;           // vpd vz
		UInt_t         mTriggerWord;     // self-defined trigger word - see code for details
		UInt_t         mTriggerWordMtd;  // Added in Run14, st_mtd trigger
		UShort_t       mRefMultNeg;      // TPC refMult neg
		UShort_t       mRefMult;// TPC refMult neg+pos 
		UShort_t       mGRefMult;// 

		UShort_t       mGRefMultCorr;
		UChar_t        mCentrality9;
		UChar_t        mCentrality16;
		Float_t        mReweight_GRefMult; //weight of gRefMult

		UInt_t        mZDCx;           // zdcX
		UInt_t        mBBCx;
		Float_t       mBackgroundRate;
		UShort_t      mbTofTrayMultiplicity ; // BTOF tray multiplicity
		UShort_t      mNumberOfGlobalTracks ; // # of global tracks
		UShort_t      mNHitsHFT[4];

		// From StMuPrimaryVertex
		Float_t mRanking ;
		UShort_t mNBEMCMatch ;
		UShort_t mNBTOFMatch ;

		// Online HT/JP thresholds
		UChar_t mHT_Th[4];

		Short_t   mQ0x; // *100
		Short_t   mQ0y; // *100
		Short_t   mQx; // *100
		Short_t   mQy; // *100
		Short_t   mEventplane; // *100
		Short_t   mEventplane0; // *100
		Short_t   mEventplane1; // *100
		Short_t   mEventplane2; // *100
		Short_t   mQ1x; // *100
		Short_t   mQ1y; // *100
		Short_t   mQ2x; // *100
		Short_t   mQ2y; // *100
		Short_t   mQplusx; // *100
		Short_t   mQplusy; // *100
		Short_t   mQminusx; // *100
		Short_t   mQminusy; // *100
		Short_t   mQetaplusx; // *100
		Short_t   mQetaplusy; // *100
		Short_t   mQetaminusx; // *100
		Short_t   mQetaminusy; // *100

		friend class StPicoAnaTree;

		ClassDef(StEventHeader,1)
};

#endif
