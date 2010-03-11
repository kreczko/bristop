/*
 * Run.cc
 *
 *  Created on: Mar 11, 2010
 *      Author: lkreczko
 */
#include "cms_bristol_ana_v3.hh"
#include "TSystem.h"
#include "TStopwatch.h"

int main(int argc, char **argv) {
	   TStopwatch watch;
		watch.Start();

		//gSystem->SetIncludePath(" -I/software/cms/slc4_ia32_gcc345/lcg/roofit/5.25.02-cms/include");
	//	gROOT->ProcessLine(".exception");
		//FIXME: this needs to be in
//		gSystem->CompileMacro("./cms_bristol_ana_v3.cc","k");
	   	ana *myana  = new ana();


		//Switches
		myana->SetData(0); //1 for data, 0 for MC
		myana->SetGoodRuns(0); //0 for all runs, 1 for good runs
	        myana->SetLHCEnergyInTeV( 7 ); //LHC E_cm

		// HLT
		myana->CheckTrigger(1); //0 for no trigger check, 1 for single electron trigger
	        //myana->CheckTrigger(1,"HLT_Ele15_SW_L1R"); //<-- uncomment this if want a diff trigger


		myana->Validation(0);
		myana->PlotRelisoNES(1);
		myana->SetDebug(0);
	        myana->ConversionStudySwitch(0);
	        myana->UseMissLayers( false );


	        // jet options: Default (=type1), SC5, KT4, pfjet
	        // met options: Default (=muon), calomet_mujes, SC5, KT4, tcmet, or pfmet
	        myana->SetJetAlgo("Default");
	        myana->SetMETAlgo("calomet_mujes");


	        // Set cuts
	        myana->SetEleETcut( 30.0 );
	        myana->SetMuonPTcut( 20.0 );
	        myana->SetJetETcut( 30.0 );
	        myana->SetMETcut( 20.0 ); //<-----

		// valid options: robustTight (Def), robustLoose, loose, tight, none
	        myana->SetEleID( ana::robustTight );

		myana->ApplyMETcut(true); //Def=true
		myana->RejectEndcapEle(true); //Def=true

		// Set cuts for QCD control sample
	        myana->SetAESMETcut(  15.0 );
	        myana->SetAESHTcut(  200.0 );
	        myana->SetAESZveto_TwoEle( true );


	        // PDF unceratinty
	        //myana->StudyPDFunc(true); //D=false



		// SD
	        myana->SetRunOnSD( true );
		bool use_MG_HLT_skims = true;
	        myana->SetRunOnMyHLTskim( use_MG_HLT_skims ); //MG HLTskim


		//Input files
		bool test_run = true;

	        if(test_run) {
		    myana->SetInputFile("/storage/top/mc/summer09_7TeV/MG/HLTskim_ttjet_mg_7TeV_aod_INVALID/*.root");
	//	    myana->SetInputFile("/storage/top/mc/Zprime_M4TeV_W400GeV-madgraph_Summer09-MC_31X_V3_7TeV-v1_GEN-SIM-RECO_nTuple/*.root");
	//	    myana->SetInputFile("/storage/top/mc/Zprime_M500GeV_W50GeV-madgraph_Summer09-MC_31X_V3_7TeV-v1_GEN-SIM-RECO_nTuple/*.root");
		    myana->SetInputFile("/storage/top/mc/Zprime_M500GeV_W5GeV-madgraph_Summer09-MC_31X_V3_7TeV-v3_GEN-SIM-RECO_nTuple/*.root");
	        }
	        else{

		    // madgraph: ttj, wj, zj
		    if( use_MG_HLT_skims ) {
			myana->SetInputFile("/storage/top/mc/summer09_7TeV/MG/HLTskim_ttjet_mg_7TeV_aod_INVALID/*.root");
			myana->SetInputFile("/storage/top/mc/summer09_7TeV/MG/HLTskim_wjets_mg_7TeV_INVALID/*.root");
			myana->SetInputFile("/storage/top/mc/summer09_7TeV/MG/HLTskim_zjet_mg_7TeV_INVALID/*.root");
		    } else {
			myana->SetInputFile("/storage/top/mc/summer09_7TeV/MG/ttjet_mg_7TeV_aod_INVALID/*.root");
			myana->SetInputFile("/storage/top/mc/summer09_7TeV/MG/wjets_mg_7TeV_INVALID/*.root");
			myana->SetInputFile("/storage/top/mc/summer09_7TeV/MG/zjet_mg_7TeV_INVALID/*.root");
		    }

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
		    myana->SetInputFile("/storage/top/mc/Zprime_M500GeV_W5GeV-madgraph_Summer09-MC_31X_V3_7TeV-v3_GEN-SIM-RECO_nTuple/*.root");
		}

		myana->SetOutputFirstName("test2");

		myana->SetLimit(-1);
		myana->EventLoop();

		if(!test_run){
		    myana->EstimateQCD();
		    myana->EstimateWjets();
	       	}
	        watch.Stop();
	        watch.Print();
	return 0;
}
