//==================================================================
//    Automatic script for doing analysis with ana.cc
//==================================================================
{        
	gSystem->CompileMacro("../../cms_bristol_ana_v3.cc","f");
   	ana* myana  = new ana();

	//Switches
	myana->SetData(0); //1 for data, 0 for MC
	//myana->SetGoodRuns(0); //0 for all runs, 1 for good runs
	//myana->CheckTrigger(1); //0 for no trigger check, 1 for single electron trigger

	//myana->SetIntLuminosity( 20.0 ); //pb-1

	//myana->SetIntLumiForM3( 10 ); //wanted/original int lumi

	myana->SetNtoyForM3Fit( 10000 );
	myana->EstimateWjets("FullRun_btag_new_mc_mixture.root");

}
