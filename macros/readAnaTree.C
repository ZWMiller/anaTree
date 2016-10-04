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
   gSystem->Load("StAnaTreeQAMaker");

   chain = new StChain();

   TString runList = "./runNumberList_completed_run15pAu";
   int number_of_runs = 0;
   std::string line;
   std::ifstream myfile(runList);
   while (std::getline(myfile, line))
     ++number_of_runs;

   StPicoAnaTreeMaker *treeMaker = new StPicoAnaTreeMaker(0,inputFile,0);
   treeMaker->setInputRunList(runList);
   StMyAnaTreeMaker *anaMaker = new StMyAnaTreeMaker("ana",treeMaker,outputFile,mixedEvent);
   //-1 - all, 0 - MB, 1 - HT0, 2 - HT1, 3 - HT2, 4 - HT3, 5 - EMu, 6 - dimuon..
   cout<<"Trigger chosen: "<<trigSelect<<endl;
   anaMaker->setTrigSelect(trigSelect);
   anaMaker->setRunList(runList);
   anaMaker->setNumberOfRuns(number_of_runs);
   anaMaker->addTrigger(500203,0); //0 -BHT0, 1-BHT1, 2-BHT2, 3-BHT3, 4-MB
   anaMaker->addTrigger(500213,0);
   anaMaker->addTrigger(500201,0);
   anaMaker->addTrigger(500204,1);
   anaMaker->addTrigger(500214,1);
   anaMaker->addTrigger(500202,1);
   //anaMaker->addTrigger(500206,1);
   anaMaker->addTrigger(500205,2);
   anaMaker->addTrigger(500215,2);

   if(trigSelect==-1||trigSelect==4||trigSelect==5||trigSelect==6){
     anaMaker->setVzCut(-100,100);
     anaMaker->setVzDiffCut(-3,3);
   }

   TString qaOutName = outputFile;
   qaOutName.ReplaceAll(".hists.root",".qa.hists.root");
   StAnaTreeQAMaker* qaMaker = new StAnaTreeQAMaker("qa",treeMaker,qaOutName.Data());
   qaMaker->setRunList(runList);
   qaMaker->setNumberOfRuns(number_of_runs);
   
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
