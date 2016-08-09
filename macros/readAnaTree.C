#include <TSystem>

class StMaker;
class StChain;
class StPicoDstMaker;


StChain *chain;
void readAnaTree(Int_t nEvents = 20000000, const Char_t *inputFile="processedRunsShort.list", const Char_t *outputFile="anaTree.hists.root", int trigSelect = 2, bool mixedEvent=false)
{

	//nEvents = 1000;	
	//Load all the System libraries

	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();

	gSystem->Load("StPicoDstMaker");
	gSystem->Load("StPicoAnaTreeMaker");
	gSystem->Load("StMyAnaTreeMaker");

	chain = new StChain();
   
   StPicoAnaTreeMaker *treeMaker = new StPicoAnaTreeMaker(0,inputFile,0);
	StMyAnaTreeMaker *anaMaker = new StMyAnaTreeMaker("ana",treeMaker,outputFile,mixedEvent);
	//-1 - all, 0 - MB, 1 - HT0, 2 - HT1, 3 - HT2, 4 - HT3, 5 - EMu, 6 - dimuon..
	cout<<"Trigger chosen: "<<trigSelect<<endl;
	anaMaker->setTrigSelect(trigSelect);
   anaMaker->addTrigger(480201,0); //0 -BHT0, 1-BHT1, 2-BHT2, 3-BHT3, 4-MB
   anaMaker->addTrigger(470211,0);
   anaMaker->addTrigger(490201,0);
   anaMaker->addTrigger(480203,0);
   anaMaker->addTrigger(480202,1);
   anaMaker->addTrigger(470202,1);
   anaMaker->addTrigger(490202,1);
   anaMaker->addTrigger(480204,1);
   anaMaker->addTrigger(470204,1);
   anaMaker->addTrigger(490204,1);
   //anaMaker->addTrigger(480206,1);
   anaMaker->addTrigger(480205,2);
   anaMaker->addTrigger(470205,2);
   anaMaker->addTrigger(490205,2);
	if(trigSelect==0){
		anaMaker->setVzCut(-6,6);
		anaMaker->setVzDiffCut(-3,3);
	}	
	if(trigSelect==2||trigSelect==3){
		anaMaker->setVzCut(-30,30);
		anaMaker->setVzDiffCut(-3,3);
	}	

	if(trigSelect==-1||trigSelect==4||trigSelect==5||trigSelect==6){
		anaMaker->setVzCut(-100,100);
		anaMaker->setVzDiffCut(-3,3);
	}
	chain->Init();
	cout<<"chain->Init();"<<endl;
	int total = treeMaker->chain()->GetEntries();
	cout << " Total entries = " << total << endl;
	if(nEvents>total) nEvents = total;
	for (Int_t i=0; i<nEvents; i++){
		if(i%100000==0)
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
	cout << "total number of events  " << nEvents << endl;
	cout << "****************************************** " << endl;

	delete chain;


}
