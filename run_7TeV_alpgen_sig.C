//==================================================================
//    Automatic script for doing analysis with ana.cc
//
// Example of running on  7 Tev ALPGEN signal only
//==================================================================
{
        TStopwatch watch;
	watch.Start();


	gSystem->SetIncludePath(" -I/software/cms/slc4_ia32_gcc345/lcg/roofit/5.25.02-cms/include");

	gSystem->CompileMacro("cms_bristol_ana_v3.cc","k");
   	ana *myana  = new ana();


	//Switches
	myana->SetData(0); //1 for data, 0 for MC
	myana->SetGoodRuns(0); //0 for all runs, 1 for good runs
        myana->SetLHCEnergyInTeV( 7 ); //LHC E_cm

	myana->CheckTrigger(1); //0 for no trigger check, 1 for single electron trigger
	//myana->CheckTrigger(1,"HLT_Ele15_SW_L1R"); //use non-default trigger


	myana->Validation(0);
	myana->PlotRelisoNES(1);
	myana->SetDebug(0);
        myana->ConversionStudySwitch(0);
        myana->UseMissLayers( false );


        // jet options: Default (=type1), SC5, KT4, pfjet
        // met options: Default (=muon), calomet_mujes, calomet_muL2, SC5, KT4, tcmet, or pfmet
        myana->SetJetAlgo("Default");
        myana->SetMETAlgo("calomet_mujes");


        // Set cuts
        myana->SetEleETcut( 30.0 );
        myana->SetMuonPTcut( 20.0 );
        myana->SetJetETcut( 30.0 );
        myana->SetMETcut( 30.0 ); //<-----

	// valid options: robustTight (Def), robustLoose, loose, tight, none
        myana->SetEleID( ana::robustTight );

	myana->ApplyMETcut(true); //D=true
	myana->RejectEndcapEle(true); //D=true

	// Set cuts for QCD control sample
        myana->SetAESMETcut(  15.0 );
        myana->SetAESHTcut(  200.0 );
        myana->SetAESZveto_TwoEle( true );


	// PDF unc
	//myana->StudyPDFunc(true);
	

	// SD
        //myana->SetRunOnSD( true );
	//bool use_MG_HLT_skims = true;
        //myana->SetRunOnMyHLTskim( use_MG_HLT_skims ); //MG HLTskim



        myana->SignalIsAlpgen(true);//D=false
        myana->SignalAlpgenThreshold(40); //D=40



	//Input files
	// alpgen ttbar signal
	myana->SetInputFile("/storage/top/mc/summer09_7TeV/AlpGen/tt0j_40Thres/*.root");
	myana->SetInputFile("/storage/top/mc/summer09_7TeV/AlpGen/tt1j_40Thres/*.root");
	myana->SetInputFile("/storage/top/mc/summer09_7TeV/AlpGen/tt2j_40Thres/*.root");
	myana->SetInputFile("/storage/top/mc/summer09_7TeV/AlpGen/tt3j_40Thres/*.root");
	myana->SetInputFile("/storage/top/mc/summer09_7TeV/AlpGen/tt4j_40Thres/*.root");

	myana->SetOutputFirstName("test");

	myana->SetLimit( -1 );
	myana->EventLoop();
	
	//myana->EstimateQCD();
	//myana->EstimateWjets();
       	
        watch.Stop();
        watch.Print();
}
