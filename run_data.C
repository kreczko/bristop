//==================================================================
//    Automatic script for doing analysis with ana.cc
//==================================================================
{
        TStopwatch watch;
	watch.Start();

	// include roofit
	gSystem->AddIncludePath("-I/software/cms/slc4_ia32_gcc345/lcg/roofit/5.25.02-cms/include");

	gSystem->CompileMacro("./cms_bristol_ana_v3.cc","k");
   	ana *myana  = new ana();



	//Switches
	myana->SetData( 1 ); //1 for data, 0 for MC
	myana->SetGoodRuns(0); //0 for all runs, 1 for good runs
        myana->SetLHCEnergyInTeV( 7 ); //LHC E_cm
        myana->SetIntLuminosity( 10.0 ); //integrated lumi in pb-1

	// HLT, two options
	// (a) CheckTrigger(true); to use default HLT_Ele15_SW_L1R
	// (b) CheckTrigger(true,"HLT_NEW"); to use other HLT bit
	myana->CheckTrigger(false); //0 for no trigger check, 1 for single electron trigger
	//myana->CheckTrigger(true,"HLT_Ele15_SW_L1R"); //<-- uncomment this if want a diff trigger

	myana->Validation(0);
	myana->PlotRelisoNES( 0 );
	myana->SetDebug( 0 );
        myana->ConversionStudySwitch(0);
        myana->UseMissLayers( true );

	// 35X 
	myana->SetRunOn35Xntuples(1);
	myana->CleanEvents(1); //apply scraping filter and PV filter.

	//use d0/d0error - default behaviour is false. Note, can only be used on Spring10 ntuples (RunOn35X must be true)
	myana->UseD0Significance( true );
	myana->NtupleHasD0PVBS( true );
	// choose which d0 to use: BS, PV or PVwithBS
	myana->D0ReferencePoint("PV");
	myana->Setd0Cut(3); //<-------- disable d0 cut by setting a huge upper limit


        // jet options: Default (=type1), SC5, KT4, pfjet
        // met options: Default (=muon), calomet_mujes, calomet_muL2, SC5, KT4, tcmet, or pfmet
        myana->SetJetAlgo("Default");
        myana->SetMETAlgo("calomet_mujes");


        // Set cuts
        myana->SetEleETcut( 30.0 );
        myana->SetMuonPTcut( 10.0 );
        myana->SetJetPTcut( 30.0 );
        myana->SetMETcut( 0.0 ); //<-----

	// valid options: robustTight (Def), robustLoose, loose, tight
        myana->SetEleID( ana::robustTight );

	myana->ApplyMETcut(false); //Def=true
	myana->RejectEndcapEle(true); //Def=false. Note this cut is applied after the conversion veto.

	// Set cuts for QCD control sample
        myana->SetAESMETcut(  15.0 );
        myana->SetAESHTcut(  200.0 );
        myana->SetAESZveto_TwoEle( true );

	// PDF unceratinty
	//myana->StudyPDFunc(false);
	

	// SD - default is false. These are relevant for 314 samples
        //myana->SetRunOnSD( false );

	//bool use_MG_HLT_skims = true;
	//Note - for 35X samples, this bool applies to pythia samples as well

        myana->SetRunOnMyHLTskim( true ); //bristol HLTskim


	//Input files
	//35X samples
	myana->SetInputFile("/storage/top/data/MinimumBias_Commissioning10-GR_R_35X_V7A_SD_EG-v2_RECO__TL_18May_361_10LS_json_v2/*root");
	myana->SetInputFile("/storage/top/data/MinimumBias_Commissioning10-May6thPDSkim2_SD_EG-v1_RECO__TL_18May_361_10LS_json/*root");
	myana->SetInputFile("/storage/top/data/MinimumBias_Commissioning10-SD_EG-v9_RECO__TL_18May_361_20k/*root");
	// v9, json is 132440-135175
	//myana->SetInputFile("/storage/top/data/MinimumBias_Commissioning10-SD_EG-v9_RECO__TL_19May_361_20LS_json/*root");
	

	myana->SetOutputFirstName("test");

	myana->SetLimit( -100 );
	
	myana->Begin();
	myana->EventLoop();
	
        watch.Stop();
        watch.Print();
}
