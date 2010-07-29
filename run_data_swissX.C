//==================================================================
//    Automatic script for doing analysis with ana.cc
//==================================================================
{

  cout << "\n\n EE misalignment corrected, Hcal noise filter, no L1TechTrig requirement.\n\n" << endl;

        TStopwatch watch;
	watch.Start();

	// include roofit
	gSystem->AddIncludePath("-I/software/cms/slc4_ia32_gcc345/lcg/roofit/5.25.02-cms/include");

	gSystem->CompileMacro("./cms_bristol_ana_v3.cc","k");
   	ana *myana  = new ana();



	//Switches
	myana->SetData( true ); //1 for data, 0 for MC
	myana->SetGoodRuns(0); //0 for all runs, 1 for good runs
        myana->SetLHCEnergyInTeV( 7 ); //LHC E_cm
        myana->SetIntLuminosity( 10.0 ); //integrated lumi in pb-1

	// HLT, two options
	// (a) CheckTrigger(true); to use default HLT_Ele15_SW_L1R
	// (b) CheckTrigger(true,"HLT_NEW"); to use other HLT bit
	//myana->CheckTrigger(true); //0 for no trigger check, 1 for single electron trigger
	// myana->CheckTrigger(true,"HLT_Ele15_SW_L1R"); //<-- uncomment this if want a diff trigger
	myana->CheckTrigger(true,"HLT_Photon15_Cleaned_L1R"); //<-- uncomment this if want a diff trigger

	myana->Validation(0);
	myana->AfterCutsPlots( 0 );//make "performance" plots after selection
	myana->PlotRelisoNES( 0 );
	myana->SetDebug( 0 );
        myana->ConversionStudySwitch(0);
        myana->UseMissLayers( true );
        myana->ApplyConversionVeto( 0 );


	// 35X 
	myana->SetRunOn35Xntuples(1);
	myana->CleanEvents(true); //apply scraping filter and PV filter.
	myana->RemoveScraping(true); //apply scraping filter.
	myana->CutOnSwissCross(1); //apply spike cleaning, swiss cross<0.95

	//use d0/d0error - default behaviour is false. Note, can only be used on Spring10 ntuples (RunOn35X must be true)
	myana->NtupleHasD0PVBS( true );
	myana->NtupleHasSwissCross( true );

	myana->UseD0Significance( false );
	myana->D0ReferencePoint("BS"); 	// choose which d0 to use: BS, PV or PVwithBS
	myana->Setd0Cut(0.02); //<-------- disable d0 cut by setting a huge upper limit


        // jet options: Default (=type1), SC5, KT4, pfjet
        // met options: Default (=muon), calomet_mujes, calomet_muL2, SC5, KT4, tcmet, or pfmet
        myana->SetJetAlgo("Default");
        myana->SetMETAlgo("calomet_mujes");


        // Set cuts
        myana->SetEleETcut( 30.0 );
        myana->SetMuonPTcut( 10.0 );
        myana->SetJetPTcut( 30.0 );
        myana->SetMETcut( 0.0 ); //<-----

	// valid options: robustTight, robustLoose, loose, tight, VBTF_W70 (Def)
        myana->SetEleID( ana::VBTF_W70);

	myana->ApplyMETcut(false); //Def=true
	myana->RejectEndcapEle(false); //Def=false. Note this cut is applied after the conversion veto.

	// Set cuts for QCD control sample
        myana->SetAESMETcut(  15.0 );
        myana->SetAESHTcut(  200.0 );
        myana->SetAESZveto_TwoEle( true );

	// PDF unceratinty
	//myana->StudyPDFunc(false);
	
	myana->SetRunOnV4ntuples(1);
	// SD - default is false. These are relevant for 314 samples
        //myana->SetRunOnSD( false );

	//bool use_MG_HLT_skims = true;
	//Note - for 35X samples, this bool applies to pythia samples as well

	//  myana->SetRunOnMyHLTskim( true ); //bristol HLTskim



	//data ntuples
        myana->SetInputFile("/storage/top/data/250nb/EG_Run2010A-Jun14thReReco_v1_RECO/*root");
        myana->SetInputFile("/storage/top/data/250nb/MinimumBias_Commissioning10-SD_EG-Jun14thSkim_v1_RECO/*root");
	myana->SetInputFile("/storage/top/data/250nb/EG_Run2010A-PromptReco-v4_RECO_137437_139558/*root");
	myana->SetInputFile("/storage/top/data/250nb/EG_Run2010A-Jul16-v4_RECO_139559_140159/*root");
	myana->SetInputFile("/storage/top/data/250nb/EG_Run2010A-PromptReco-v4_RECO_140160_140399/*root");



	/*
	//ntuples of e+4jet events
	myana->SetInputFile("/storage/top/data/reco_edm/nTuple_data_EG_140385_lumi_101_ev_90009543.root");
	myana->SetInputFile("/storage/top/data/reco_edm/nTuple_data_EG_140331_440601613.root");
	myana->SetInputFile("/storage/top/data/reco_edm/nTuple_run_139195_lumi_77_ev_69244083.root");
	*/



	


	myana->SetOutputFirstName("test");

	myana->SetLimit( -2000 );
	
	myana->Begin();
	myana->EventLoop();
	
        watch.Stop();
        watch.Print();
}
