
#include <TSystem>

class StMaker;
class StChain;
class StPicoDstMaker;
class StMuDstMaker;


StChain *chain;
void makeAnaTree(const Int_t runnumber=16086001,
		//    const Char_t *inputFile="/star/data54/reco/AuAu200_production_2011/FullField/P11id/2011/169/12169026/st_physics_adc_12169026_raw_4510001.MuDst.root",
		//    const Char_t *inputFile="/star/data78/reco/pp200_production_2012/ReversedFullField/P12id/2012/040/13040016/st_physics_13040016_raw_1010001.MuDst.root",
		//    const Char_t *inputFile="/star/data43/reco/pp500_production_2013/ReversedFullField/P14ia/2013/115/14115072/st_physics_14115072_raw_3690004.MuDst.root",
		//    const Char_t *inputFile="root://xrdstar.rcf.bnl.gov:1095//home/starlib/home/starreco/reco/pp500_production_2011/ReversedFullField/P11id/2011/042/12042026/st_physics_adc_12042026_raw_2500001.MuDst.root",
		//    const Char_t *inputFile="/star/u/xgn1992/work/14PicoDst/PdsfData/15087019/0/st_physics_15087019_raw_0000050.MuDst.root",
		//    const Char_t *inputFile="st_physics_15166010_raw_3000054.MuDst.root",
		//    const Char_t *inputFile="/star/data79/reco/AuAu_200_production_low_2014/ReversedFullField/P15ic/2014/145/15145024/st_physics_15145024_raw_1000048.MuDst.root",
		//    const Char_t *inputFile="root://xrdstar.rcf.bnl.gov:1095//home/starlib/home/starreco/reco/AuAu_200_production_low_2014/ReversedFullField/P15ic/2014/166/15166010/st_physics_15166010_raw_4500060.MuDst.root",
	//	const Char_t *inputFile="root://xrdstar.rcf.bnl.gov:1095//home/starlib/home/starreco/reco/AuAu_200_production_mid_2014/ReversedFullField/P15ic/2014/094/15094070/st_physics_15094070_raw_0000007.MuDst.root",
		const Char_t *inputFile="testFiles/st_physics_adc_16086001_raw_4000030.MuDst.root",
		const bool creatingPhiWgt = kFALSE, const int prodMod = 1, const int emcMode=1, const int prodType = 0
		){
	Int_t nEvents = 10000000;
	//Int_t nEvents = 10;
	//Int_t nEvents = 1000;	
	//Load all the System libraries

	gSystem->Load("libTable");
	gSystem->Load("libPhysics");
	gSystem->Load("St_base");
	gSystem->Load("StChain");
	gSystem->Load("St_Tables");
	gSystem->Load("StUtilities");        // new addition 22jul99
	gSystem->Load("StTreeMaker");
	gSystem->Load("StIOMaker");
	gSystem->Load("StarClassLibrary");
	gSystem->Load("StTriggerDataMaker"); // new starting from April 2003
	gSystem->Load("StBichsel");
	gSystem->Load("StEvent");
	gSystem->Load("StEventUtilities");
	gSystem->Load("StDbLib");
	gSystem->Load("StEmcUtil");
	gSystem->Load("StTofUtil");
	gSystem->Load("StPmdUtil");
	gSystem->Load("StPreEclMaker");
	gSystem->Load("StStrangeMuDstMaker");
	gSystem->Load("StMuDSTMaker");

	if(!creatingPhiWgt&&emcMode) {
		gSystem->Load("StTpcDb");
		gSystem->Load("StMcEvent");
		gSystem->Load("StMcEventMaker");
		gSystem->Load("StDaqLib");
		gSystem->Load("libgen_Tables");
		gSystem->Load("libsim_Tables");
		gSystem->Load("libglobal_Tables");
		gSystem->Load("StEmcTriggerMaker");
		gSystem->Load("StEmcUtil");//mine
		gSystem->Load("StEmcRawMaker");
		gSystem->Load("StEmcADCtoEMaker");
		gSystem->Load("StPreEclMaker");
		gSystem->Load("StEpcMaker");
		gSystem->Load("StEmcSimulatorMaker");
		gSystem->Load("StEmcUtil");
		gSystem->Load("StDbBroker");
		gSystem->Load("StDetectorDbMaker");
		gSystem->Load("StDbUtilities");
		gSystem->Load("StEEmcUtil");
		gSystem->Load("StEEmcDbMaker");
		gSystem->Load("St_db_Maker");
		gSystem->Load("StTriggerUtilities");
	}

	// New additions for MTD fixes in Run14
	gSystem->Load("StMagF");
	gSystem->Load("StMtdUtil");
	gSystem->Load("StMtdMatchMaker");
	//gSystem->Load("StMtdCalibMaker");

	gSystem->Load("StPicoDstMaker");
	gSystem->Load("StPicoAnaTreeMaker");
	gSystem->Load("StPicoQAMaker");
	gSystem->Load("StRefMultCorr");
   gSystem->Load("StPicoElecPurityMaker");

	chain = new StChain();

	StMuDstMaker *MuDstMaker = new StMuDstMaker(0,0,"",inputFile,"MuDst",100);
	MuDstMaker->SetStatus("*",0);
	MuDstMaker->SetStatus("MuEvent",1);
	MuDstMaker->SetStatus("PrimaryVertices",1);
	MuDstMaker->SetStatus("PrimaryTracks",1);
	MuDstMaker->SetStatus("GlobalTracks",1);
	MuDstMaker->SetStatus("CovGlobTrack",1);
	MuDstMaker->SetStatus("BTof*",1);
	MuDstMaker->SetStatus("Emc*",1);
	MuDstMaker->SetStatus("MTD*",1);


	//StMagFMaker *magfMk = new StMagFMaker; 
	//StMtdMatchMaker *mtdMatchMaker = new StMtdMatchMaker();
	//StMtdCalibMaker *mtdCalibMaker = new StMtdCalibMaker("mtdcalib"); 

	if(!creatingPhiWgt&&emcMode) {
		St_db_Maker *dbMk = new St_db_Maker("db","MySQL:StarDb","$STAR/StarDb","StarDb");

		// Endcap database
		StEEmcDbMaker* eemcDb = new StEEmcDbMaker;

		StEmcADCtoEMaker *adc2e = new StEmcADCtoEMaker();
		adc2e->setPrint(false);
		//adc2e->setFillHisto(false);
		//adc2e->setDebug(false); //more histograms
		//adc2e->setSMDRmsCut(0,0);
		adc2e->saveAllStEvent(true);
		//adc2e->setRemoveGhostEvent(false);
		//adc2e->setTallyHist(mTally);
		//adc2e->setDbName("Calibrations/emc/y3");

		StPreEclMaker *pre_ecl=new StPreEclMaker();
		pre_ecl->setPrint(kFALSE);
		StEpcMaker *epc=new StEpcMaker();
		epc->setPrint(kFALSE);

      if(prodMod==1){
		   // Trigger simulator
		   StTriggerSimuMaker* trigSimu = new StTriggerSimuMaker;
		   trigSimu->setMC(false);
		   trigSimu->useBemc();
		   trigSimu->useEemc();
		   trigSimu->useOnlineDB();
		   trigSimu->bemc->setConfig(StBemcTriggerSimu::kOffline);
      }
	}


	StPicoDstMaker *picoMaker = new StPicoDstMaker(1,inputFile,"picoDst");
	picoMaker->setRunNumber(runnumber);
	picoMaker->setProdMode(prodMod); // 0-mb, 1-ht, 4-only e or mu  
	picoMaker->setEmcMode(emcMode); // 0-No EMC, 1-EMC ON
	//picoMaker->SetDebug(1);

	/*StPicoMtdCalibMaker *calibMaker;
	if(prodMod==0 || prodMod==1){
		calibMaker = new StPicoMtdCalibMaker();
		calibMaker->setInitFromFile(kTRUE);
		calibMaker->setApplyT0(kTRUE);
		calibMaker->setCalibFileT0("StRoot/StPicoDstMaker/Run14_AuAu200_CalibDtof.offline.dat");
	}*/
	
	TString inputFileName = inputFile;
	Int_t index = inputFileName.Index("st_");
	TString mInputFileName="";
	for(int i=index;i<(int)inputFileName.Length();i++) {
		mInputFileName.Append(inputFileName(i));
	}

	TString outputFile,outQAFile;
	outQAFile=mInputFileName;
	outQAFile.ReplaceAll("MuDst.root","qa.root");
   StPicoQAMaker *qaMaker = new StPicoQAMaker("ana",picoMaker,outQAFile);
   qaMaker->setRunList("./runNumberList_run15pp_transverse");
   
   TString outPurityFile=mInputFileName;
   outPurityFile.ReplaceAll("MuDst.root","purity.root");
   StPicoElecPurityMaker *ePurMaker = new StPicoElecPurityMaker("purity",picoMaker,outPurityFile);
   ePurMaker->setRunMode(0); // 0 - pp (no centrality), 1- AuAu (centrality on)
   ePurMaker->addTrigger(480201,0); //0 -BHT0, 1-BHT1, 2-BHT2, 3-BHT3, 4-MB
   ePurMaker->addTrigger(470211,0);
   ePurMaker->addTrigger(490201,0);
   ePurMaker->addTrigger(480203,0);
   ePurMaker->addTrigger(480202,1);
   ePurMaker->addTrigger(470202,1);
   ePurMaker->addTrigger(490202,1);
   ePurMaker->addTrigger(480204,1);
   ePurMaker->addTrigger(470204,1);
   ePurMaker->addTrigger(490204,1);
   ePurMaker->addTrigger(480206,1);
   ePurMaker->addTrigger(480205,2);
   ePurMaker->addTrigger(470205,2);
   ePurMaker->addTrigger(490205,2);
	
   outputFile=mInputFileName;
	outputFile.ReplaceAll("MuDst.root","anaTree.root");
	
	StPicoAnaTreeMaker *treeMaker = new StPicoAnaTreeMaker(1,outputFile,picoMaker);
	treeMaker->setTriggerSelection(prodMod); //0-mb, 1-ht, 2-mtd
   treeMaker->setPicoMode(1); //0-Run14 PicoDstMaker, 1-Run15+ PicoDstMaker
	if(prodMod==0){
		treeMaker->setVzCut(-8,8);
		treeMaker->setVzDiffCut(-4,4);
      treeMaker->setInputRunList("./runNumberList_run15pp_transverse");
      treeMaker->addTrigger(480202);
   //   treeMaker->setInputRecenterFile("./recenter_correction.root");
      treeMaker->setPhoEPairMassCut(0.2);
	}	
	if(prodMod==1){
		treeMaker->setVzCut(-100,100);
      treeMaker->setPhoEPairMassCut(0.24);
      treeMaker->setSaveHadron(true);
      //treeMaker->setDoEvtPlane(false); //default is true
      if(prodType==0){ // prod low and mid
        treeMaker->setInputRunList("./runNumberList_run15pp_transverse");
        treeMaker->addTrigger(480201); 
        treeMaker->addTrigger(470211);
        treeMaker->addTrigger(490201);
        treeMaker->addTrigger(480203);
        treeMaker->addTrigger(480202);
        treeMaker->addTrigger(470202);
        treeMaker->addTrigger(490202);
        treeMaker->addTrigger(480204);
        treeMaker->addTrigger(470204);
        treeMaker->addTrigger(490204);
        treeMaker->addTrigger(480206);
        treeMaker->addTrigger(480205);
        treeMaker->addTrigger(470205);
        treeMaker->addTrigger(490205);
         //    treeMaker->setInputRecenterFile("./recenter_correction.root");
         treeMaker->setMaxRunId(1700);
      }
      if(prodType==1){ // prod high
      treeMaker->setInputRunList("./runNumberList_run15pp_transverse");
      //   treeMaker->setInputRecenterFile("./recenter_correction_ht_high.root");
         treeMaker->setMaxRunId(1000); //need to check
      }
	}
	if(prodMod==2){
		treeMaker->setVzCut(-100,100);
      treeMaker->addTrigger(490600);
      treeMaker->addTrigger(490602);
      //treeMaker->setDoEvtPlane(false); //default is true
      //treeMaker->setInputRunList("./runNumberList_run14AuAu200");
      //treeMaker->setInputRecenterFile("./recenter_correction.root");
	}
   //treeMaker->setDoCalcRecenter(false);  //default is false
   treeMaker->setPartEnSigECut(-3.5,3);
   treeMaker->setPhoEPairDcaCut(1);


	chain->Init();
	cout<<"chain->Init();"<<endl;
	int total = 0;
	for (Int_t i=0; i<nEvents; i++){
		if(i%100==0)
			cout << "Working on eventNumber " << i << endl;

		chain->Clear();
		int iret = chain->Make(i);

		if (iret) { cout << "Bad return code!" << iret << endl; break;}

		total++;

	}

	cout << "****************************************** " << endl;
	cout << "Work done... now its time to close up shop!"<< endl;
	cout << "****************************************** " << endl;
	chain->Finish();
	cout << "****************************************** " << endl;
	cout << "total number of events  " << total << endl;
	cout << "****************************************** " << endl;

//	delete chain;


}
