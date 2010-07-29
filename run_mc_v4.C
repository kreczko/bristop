//==================================================================
//    Automatic script for doing analysis with ana.cc
//==================================================================
{
        TStopwatch watch;
	watch.Start();


	cout << "\n\n Run on  v4 e20 mc ntuples.\n" << endl;
	cout << " Use updated cross section as of 15 July" << endl;
	//cout << " all processes excpet enri1 (e20 not done yet, and it is small)" << endl;
	//cout << "  e20 skim: tt, wj, zj, QCD, t-chan, s-chan, photon+jets." << endl;

	// include roofit
	gSystem->AddIncludePath("-I/software/cms/slc4_ia32_gcc345/lcg/roofit/5.25.02-cms/include");

	gSystem->CompileMacro("./cms_bristol_ana_v3.cc","k");
   	ana *myana  = new ana();





	//Switches
	myana->SetData( 0 ); //1 for data, 0 for MC
	myana->SetGoodRuns(0); //0 for all runs, 1 for good runs
        myana->SetLHCEnergyInTeV( 7 ); //LHC E_cm
        myana->SetIntLuminosity( 10.0 ); //integrated lumi in pb-1

	// HLT, two options
	// (a) CheckTrigger(true); to use default HLT_Ele15_SW_L1R
	// (b) CheckTrigger(true,"HLT_NEW"); to use other HLT bit
	//myana->CheckTrigger(true); //0 for no trigger check, 1 for single electron trigger
	//myana->CheckTrigger(true,"HLT_Ele15_SW_L1R"); //<-- uncomment this if want a diff trigger
	myana->CheckTrigger(true,"HLT_Photon15_L1R"); //<-- uncomment this if want a diff trigger

	myana->Validation(0);
	myana->AfterCutsPlots( 0 );
	myana->PlotRelisoNES( 1 );
	myana->SetDebug( 0 );
        myana->ConversionStudySwitch(0);
        myana->UseMissLayers( true );
        myana->ApplyConversionVeto( true );

	// 35X 
	myana->SetRunOn35Xntuples(1);
	myana->CleanEvents(1); //apply scraping filter and PV filter.
	myana->RemoveScraping(0); //<-- take out scraping filter

	//use d0/d0error - default behaviour is false. Note, can only be used on Spring10 ntuples (RunOn35X must be true)
	myana->UseD0Significance( false );
	myana->D0ReferencePoint("BS");	// choose which d0 to use: BS, PV or PVwithBS
	myana->Setd0Cut(0.02); //<-------- disable d0 cut by setting a huge upper limit


        // jet options: Default (=type1), SC5, KT4, pfjet
        // met options: Default (=muon), calomet_mujes, calomet_muL2, SC5, KT4, tcmet, or pfmet
        myana->SetJetAlgo("Default");
        myana->SetMETAlgo("calomet_mujes");


        // Set cuts
        myana->SetEleETcut( 30.0 );
        myana->SetMuonPTcut( 10.0 );
        myana->SetJetPTcut( 30.0 );

	myana->ApplyMETcut(false); //Def=true
        myana->SetMETcut( 0.0 ); //<-----


	// valid options: robustTight (Def), robustLoose, loose, tight, VBTF_W70
        myana->SetEleID( ana::VBTF_W70 );
        //myana->SetEleID( ana::robustTight );

	myana->RejectEndcapEle(false); //Def=false. Note this cut is applied after the conversion veto.

	// Set cuts for QCD control sample
        myana->SetAESMETcut(  15.0 );
        myana->SetAESHTcut(  200.0 );
        myana->SetAESZveto_TwoEle( true );

	// PDF unceratinty
	//myana->StudyPDFunc(false);
	

	// SD - default is false. These are relevant for 314 samples
        //myana->SetRunOnSD( false );


	//Note - for 35X samples, this bool applies to pythia samples as well

        myana->SetRunOnMyHLTskim( false ); // HLT skim
        myana->SetRunOnMyE20skim( true ); // e20 skim


	myana->SetRunOnV4ntuples( true );

	//Input files	
	myana->SetInputFile("/storage/top/mc/spring10_7TeV_v4/MG/e20skim_ttjet/*11*.root");
	/*
	myana->SetInputFile("/storage/top/mc/spring10_7TeV_v4/MG/e20skim_ttjet/*.root");
	myana->SetInputFile("/storage/top/mc/spring10_7TeV_v4/MG/e20skim_wjet/*.root");
	myana->SetInputFile("/storage/top/mc/spring10_7TeV_v4/MG/e20skim_zjet/*.root");

	myana->SetInputFile("/storage/top/mc/spring10_7TeV_v4/pythia/e20skim_bce1/*.root");
	myana->SetInputFile("/storage/top/mc/spring10_7TeV_v4/pythia/e20skim_bce2/*.root");
	myana->SetInputFile("/storage/top/mc/spring10_7TeV_v4/pythia/e20skim_bce3/*.root");
	myana->SetInputFile("/storage/top/mc/spring10_7TeV_v4/pythia/e20skim_enri1/*.root");
	myana->SetInputFile("/storage/top/mc/spring10_7TeV_v4/pythia/e20skim_enri2/*.root");
	myana->SetInputFile("/storage/top/mc/spring10_7TeV_v4/pythia/e20skim_enri3/*.root");

	myana->SetInputFile("/storage/top/mc/spring10_7TeV_v4/MG/e20skim_tchan/*.root");
	myana->SetInputFile("/storage/top/mc/spring10_7TeV_v4/MG/e20skim_tW/*.root");

	myana->SetInputFile("/storage/top/mc/spring10_7TeV_v4/MG/e20skim_pj1/*.root");
	myana->SetInputFile("/storage/top/mc/spring10_7TeV_v4/MG/e20skim_pj2/*.root");
	myana->SetInputFile("/storage/top/mc/spring10_7TeV_v4/MG/e20skim_pj3/*.root");
	*/


	myana->SetOutputFirstName("test");

	myana->SetLimit( -10000 );
	
	myana->Begin();
	myana->EventLoop();
	
        watch.Stop();
        watch.Print();
}
