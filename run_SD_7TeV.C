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
        myana->SetLHCEnergyInTeV( 7 ); //LHC E_cm

	myana->Validation(0);
	myana->PlotRelisoNES(1);
	myana->SetDebug(0);
        myana->ConversionStudySwitch(0);

	//	myana->SetJetAlgo("pfjet"); //Default, SC5, pfjet
	//	myana->SetMETAlgo("calomet_mujes"); //Default, calomet_mujes, tcmet, or pfmet


        // Set cuts
        myana->SetEleETcut( 30.0 );
        myana->SetMuonPTcut( 20.0 );
        myana->SetJetETcut( 30.0 );
        myana->SetMETcut( 20.0 );


	// SD
        myana->SetRunOnSD( 1 ); //0 or 1
        myana->SetAESMETcut(  15.0 );
        myana->SetAESHTcut(  200.0 );
        myana->SetAESZveto_TwoEle( 1 ); //0 or 1


	//Input files
	myana->SetInputFile("/storage/OctEx09/SD_7TeV_JEC/ttbar7/*.root");

	myana->SetInputFile("/storage/OctEx09/SD_7TeV_JEC/wenu7/*.root");
	myana->SetInputFile("/storage/OctEx09/SD_7TeV_JEC/zee7/*.root");

	myana->SetInputFile("/storage/OctEx09/SD_7TeV_JEC/enri1_7/*.root");
	myana->SetInputFile("/storage/OctEx09/SD_7TeV_JEC/enri2_7/*.root");
	myana->SetInputFile("/storage/OctEx09/SD_7TeV_JEC/enri3_7/*.root");
	myana->SetInputFile("/storage/OctEx09/SD_7TeV_JEC/bce1_7/*.root");
	myana->SetInputFile("/storage/OctEx09/SD_7TeV_JEC/bce2_7/*.root");
	myana->SetInputFile("/storage/OctEx09/SD_7TeV_JEC/bce3_7/*.root");

	myana->SetOutputFirstName("test");

	myana->SetLimit(-10);
	myana->EventLoop();
	
	myana->EstimateQCD();
       	
        watch.Stop();
        watch.Print();
}
