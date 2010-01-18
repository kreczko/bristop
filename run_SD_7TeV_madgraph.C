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


	// jet options: Default (=type1), SC5, KT4, pfjet
        // met options: Default (=muon), calomet_mujes, calomet_muL2, SC5, KT4, tcmet, or pfmet
        myana->SetJetAlgo("Default");
        myana->SetMETAlgo("Default");


        // Set cuts
        myana->SetEleETcut( 30.0 );
        myana->SetMuonPTcut( 20.0 );
        myana->SetJetETcut( 30.0 );
        myana->SetMETcut( 20.0 );

        // valid options: robustTight (Def), robustLoose, loose, tight, none
        myana->SetEleID( ana::robustTight );

	// QCD estimation plan B (no met cut)
        myana->SetRunPlanB( false );

	// Define cuts for QCD control sample
        myana->SetAESMETcut(  15.0 );
        myana->SetAESHTcut(  200.0 );
        myana->SetAESZveto_TwoEle( true );

        // SD
        myana->SetRunOnSD( true );
        myana->SetRunOnMyHLTskim( true ); //MG HLTskim, ttj, wj, zj


	//Input files
	// Madgraph (HLT-skim)
        myana->SetInputFile("/storage/top/mc/summer09_7TeV/MG/HLTskim_ttjet_mg_7TeV_aod/*.root");
        myana->SetInputFile("/storage/top/mc/summer09_7TeV/MG/HLTskim_wjets_mg_7TeV/*.root");
        myana->SetInputFile("/storage/top/mc/summer09_7TeV/MG/HLTskim_zjet_mg_7TeV/*.root");
	// Madgraph (unskimmed)
        myana->SetInputFile("/storage/top/mc/summer09_7TeV/MG/sTop_tchan/*.root");
        myana->SetInputFile("/storage/top/mc/summer09_7TeV/MG/sTop_schan/*.root");
        myana->SetInputFile("/storage/top/mc/summer09_7TeV/MG/sTop_tWchan/*.root");
	// Pythia SD
	myana->SetInputFile("/storage/top/mc/OctEx09/SD_7TeV_JEC/enri1_7/*.root");
	myana->SetInputFile("/storage/top/mc/OctEx09/SD_7TeV_JEC/enri2_7/*.root");
	myana->SetInputFile("/storage/top/mc/OctEx09/SD_7TeV_JEC/enri3_7/*.root");
	myana->SetInputFile("/storage/top/mc/OctEx09/SD_7TeV_JEC/bce1_7/*.root");
	myana->SetInputFile("/storage/top/mc/OctEx09/SD_7TeV_JEC/bce2_7/*.root");
	myana->SetInputFile("/storage/top/mc/OctEx09/SD_7TeV_JEC/bce3_7/*.root");

	myana->SetOutputFirstName("test");

	myana->SetLimit(-10);
	myana->EventLoop();
	
	myana->EstimateQCD();
        myana->EstimateWjets();
       	
        watch.Stop();
        watch.Print();
}
