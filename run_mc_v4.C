//==================================================================
//    Automatic script for doing analysis with ana.cc
//==================================================================
{
        TStopwatch watch;
	watch.Start();


	cout << "\n\n Testing v4 mc ntuples.\n" << endl;
	//cout << "\n\n Run on MC ntuples.\n" << endl;
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
	myana->CheckTrigger(true,"HLT_Photon10_L1R"); //<-- uncomment this if want a diff trigger

	myana->Validation(0);
	myana->PlotRelisoNES( 1 );
	myana->SetDebug( 0 );
        myana->ConversionStudySwitch(0);
        myana->UseMissLayers( true );
        myana->ApplyConversionVeto( false );

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

	myana->RejectEndcapEle(false); //Def=false. Note this cut is applied after the conversion veto.

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

        myana->SetRunOnMyHLTskim( false ); //bristol HLT skim
        myana->SetRunOnMyE20skim( false ); //bristol e20 skim


	myana->SetRunOnV4ntuples( true );

	//Input files
	myana->SetInputFile("/storage/top/mc/spring10_7TeV_v4/link_ttjet/*.root");

	myana->SetOutputFirstName("test");

	myana->SetLimit( -10 );
	
	myana->Begin();
	myana->EventLoop();
	
        watch.Stop();
        watch.Print();
}
