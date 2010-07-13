/*
 * Run.cc
 *
 *  Created on: Mar 11, 2010
 *      Author: lkreczko
 */
//#include "cms_bristol_ana_v3.hh"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "Analysis.h"
#include <iostream>
using namespace ROOT;
using namespace std;
int main(int argc, char **argv) {
	gROOT->ProcessLine("#include <vector>");
	//	gROOT->ProcessLine("gErrorIgnoreLevel = 2001;");
	TStopwatch watch;
	watch.Start();
	Analysis* myAnalysis = new Analysis();

	//	ana *myana = new ana();

	//Switches
	//	myana->SetData(0); //1 for data, 0 for MC
	//	myana->SetGoodRuns(0); //0 for all runs, 1 for good runs
	//	myana->SetLHCEnergyInTeV(7); //LHC E_cm

	// HLT
	//	myana->CheckTrigger(1); //0 for no trigger check, 1 for single electron trigger

	//	myana->Validation(0);
	//	myana->PlotRelisoNES(1);
	//	myana->SetDebug(0);
	//	myana->ConversionStudySwitch(0);
	//	myana->UseMissLayers(false);

	// jet options: Default (=type1), SC5, KT4, pfjet
	// met options: Default (=muon), calomet_mujes, SC5, KT4, tcmet, or pfmet
	//	myana->SetJetAlgo("Default");
	//	myana->SetMETAlgo("calomet_mujes");

	// Set cuts
	//	myana->SetEleETcut(30.0);
	//	myana->SetMuonPTcut(20.0);
	//	myana->SetJetETcut(30.0);
	//	myana->SetMETcut(20.0); //<-----

	// valid options: robustTight (Def), robustLoose, loose, tight, none
	//	myana->SetEleID(ana::robustTight);

	//	myana->ApplyMETcut(true); //Def=true
	//	myana->RejectEndcapEle(true); //Def=true

	// Set cuts for QCD control sample
	//	myana->SetAESMETcut(15.0);
	//	myana->SetAESHTcut(200.0);
	//	myana->SetAESZveto_TwoEle(true);
	//	myana->SetNumberOfJetsUsedInReco(-1);

	// PDF unceratinty
	//myana->StudyPDFunc(true); //D=false


	// SD
	//	myana->SetRunOnSD(true);
	//	bool use_MG_HLT_skims = true;
	//	myana->SetRunOnMyHLTskim(use_MG_HLT_skims); //MG HLTskim

	cout << "adding files" << endl;
	//Input files
	myAnalysis->addInputFile("/storage/top/mc/spring10_7TeV_new/MG/e20skim_ttjet/e20skim_nTuple_ttjet_99_1.root");
//	myAnalysis->addInputFile("/storage/top/mc/spring10_7TeV_new/MG/e20skim_vqq/*.root");
	//	myana->SetInputFile("/storage/top/mc/spring10_7TeV_new/MG/e20skim_zjet/*.root");

//	myAnalysis->addInputFile("/storage/top/mc/spring10_7TeV_new/pythia/e20skim_bce1/*.root");
//	myAnalysis->addInputFile("/storage/top/mc/spring10_7TeV_new/pythia/e20skim_bce2/*.root");
//	myAnalysis->addInputFile("/storage/top/mc/spring10_7TeV_new/pythia/e20skim_bce3/*.root");

	//	myana->SetInputFile("/storage/top/mc/spring10_7TeV/pythia/e20skim_enri1/*.root");//old
	//	myana->SetInputFile("/storage/top/mc/spring10_7TeV_new/pythia/e20skim_enri2/*.root");
	//	myana->SetInputFile("/storage/top/mc/spring10_7TeV_new/pythia/e20skim_enri3/*.root");

	//	myana->SetInputFile("/storage/top/mc/spring10_7TeV_new/MG/e20skim_tchan/*.root");
	//	myana->SetInputFile("/storage/top/mc/spring10_7TeV_new/MG/e20skim_tW/*.root");

	//myana->SetInputFile("/storage/top/mc/link_spring10_tchan/*.root");
	//myana->SetInputFile("/storage/top/mc/link_spring10_tW/*.root");

	//	myana->SetOutputFirstName("Release");


	//	myana->SetLimit(-1);
	//	myana->EventLoop();

	//	if (!test_run) {
	//	myana->EstimateQCD("FullRun_btag_new_mc_mixture.root");
	//		myana->EstimateWjets();
	//	}
	cout << "starting analysis" << endl;
	myAnalysis->analyze();
	watch.Stop();
	watch.Print();
	return 0;
}
