//==================================================================
//    Automatic script for doing analysis with ana.cc
//==================================================================
{
        TStopwatch watch;
	watch.Start();

	gSystem->CompileMacro("./cms_bristol_ana_v3.cc","k");
   	ana* myana  = new ana();


	//Switches
	myana->SetData(0); //1 for data, 0 for MC
	myana->SetGoodRuns(0); //0 for all runs, 1 for good runs
	myana->CheckTrigger(1); //0 for no trigger check, 1 for single electron trigger

	myana->Validation(0);
	myana->PlotRelisoNES(1);
	myana->SetDebug(0);
        myana->ConversionStudySwitch(0);

	//myana->SetJetAlgo("Default"); //Default, SC5, or pfjet
	//myana->SetMETAlgo("Default"); //Default, calomet_mujes, tcmet, or pfmet


        // Set cuts
        myana->SetEleETcut(30.0);
        myana->SetMuonPTcut(20.0);
        myana->SetJetETcut(30.0);
        myana->SetMETcut(20.0);


	// SD
        myana->SetRunOnSD( 1 ); //0 or 1
        myana->SetAESMETcut(  15.0 );
        myana->SetAESHTcut(  200.0 );
        myana->SetAESZveto_TwoEle( 1 ); //0 or 1


	//Input files
	myana->SetInputFile("/storage/OctEx09/SD/ttbar/*.root");
	
	myana->SetInputFile("/storage/OctEx09/SD/wenu/*.root");
	myana->SetInputFile("/storage/OctEx09/SD/zee/*.root");

	myana->SetInputFile("/storage/OctEx09/SD/enri1/*.root");
	myana->SetInputFile("/storage/OctEx09/SD/enri2/*.root");	
	myana->SetInputFile("/storage/OctEx09/SD/enri3/*.root");
	myana->SetInputFile("/storage/OctEx09/SD/bce1/*.root");
	myana->SetInputFile("/storage/OctEx09/SD/bce2/*.root");
	myana->SetInputFile("/storage/OctEx09/SD/bce3/*.root");
	

	myana->SetOutputFirstName("test");

	myana->SetLimit(-10);
	myana->EventLoop();
	
	//myana->EstimateQCD();
       	
        watch.Stop();
        watch.Print();
}
