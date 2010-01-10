//==================================================================
//    Automatic script for doing analysis with ana.cc
//==================================================================
{
        TStopwatch watch;
	watch.Start();

	gSystem->CompileMacro("./cms_bristol_ana_v3.cc","k");
   	ana* myana  = new ana();

	//using namespace ana;

	//Switches
	myana->SetData(  1  ); //1 for data, 0 for MC
	myana->SetGoodRuns(  0  ); //0 for all runs, 1 for good runs
	myana->CheckTrigger( 0 ); //0 for no trigger check, 1 for single electron trigger
        myana->SetLHCEnergyInTeV( 7 ); //LHC E_cm

	myana->Validation(0);
	myana->PlotRelisoNES(0);
	myana->SetDebug(0);
        myana->ConversionStudySwitch(0);
	myana->UseMissLayers(true);

        // jet options: Default (=type1), SC5, KT4, pfjet
        // met options: Default (=muon), calomet_mujes, calomet_muL2, SC5, KT4, tcmet, or pfmet
        myana->SetJetAlgo("Default");
        myana->SetMETAlgo("Default");


        // Set cuts
        myana->SetEleETcut( 30.0 );
        myana->SetMuonPTcut( 20.0 );
        myana->SetJetETcut( 10.0 );
        myana->SetMETcut( 20.0 ); //<-----

	// valid options: robustTight (Def), robustLoose, loose, tight, none
        myana->SetEleID( ana::robustTight );

	myana->SetRunPlanB( false );  //bool

	// Set cuts for QCD control sample
        myana->SetAESMETcut(  15.0 );
        myana->SetAESHTcut(  200.0 );
        myana->SetAESZveto_TwoEle( true ); //bool


	// SD
        myana->SetRunOnSD( false ); //bool


	//Input files
	// PD
	//myana->SetInputFile("/storage/top/data/MinimumBias_BeamCommissioning09-PromptReco-v2_RECO_V1_L1Bit0and4041not36to39/*.root");
	// SD
	myana->SetInputFile("/storage/top/data/MinimumBias_BeamCommissioning09-SD_AllMinBias-PromptSkimCommissioning_v2_RAW-RECO_V1_L1Bit0and4041not36to39/*.root");

	myana->SetOutputFirstName("test");

	myana->SetLimit( 100000 );
	myana->EventLoop();
	
	//myana->EstimateQCD();
	//myana->EstimateWjets();
       	
        watch.Stop();
        watch.Print();
}
