//==================================================================
//    Automatic script for doing analysis with ana.cc
//==================================================================
{
        TStopwatch watch;
	watch.Start();

	gSystem->CompileMacro("cms_bristol_ana_v3.cc","k");
   	ana* myana  = new ana();


	//Switches
	myana->SetData(0); //1 for data, 0 for MC
	myana->SetGoodRuns(0); //0 for all runs, 1 for good runs
	myana->CheckTrigger(1); //0 for no trigger check, 1 for single electron trigger

	myana->Validation(0);
	myana->PlotRelisoNES(1);
	myana->SetDebug(0);
        myana->ConversionStudySwitch(0);

	//myana->SetJetAlgo("SC5");


        // Set cuts
        myana->SetEleETcut( 30.0 );
        myana->SetMuonPTcut( 20.0 );
        myana->SetJetETcut( 30.0 );
        myana->SetMETcut( 20.0 );


	//Input files
	myana->SetInputFile("~phtlc/store/summer09_10TeV_2/ttbar/*.root");
	myana->SetInputFile("~phtlc/store/summer09_10TeV_2/enri1/*.root");
	myana->SetInputFile("~phtlc/store/summer09_10TeV_2/enri2/*.root");
	myana->SetInputFile("~phtlc/store/summer09_10TeV_2/enri3/*.root");
	myana->SetInputFile("~phtlc/store/summer09_10TeV_2/bce1/*.root");
	myana->SetInputFile("~phtlc/store/summer09_10TeV_2/bce2/*.root");
	myana->SetInputFile("~phtlc/store/summer09_10TeV_2/bce3/*.root");
	//myana->SetInputFile("~phtlc/store/summer09_10TeV_3/wenu/*.root");
	//myana->SetInputFile("~phtlc/store/summer09_10TeV_3/zee/*.root");
	
	myana->SetOutputFirstName("test");

	myana->SetLimit(-10);
	myana->EventLoop();

	//myana->EstimateQCD();
	//myana->EstimateWjets();
       	
        watch.Stop();
        watch.Print();
}
