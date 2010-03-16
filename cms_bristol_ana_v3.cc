//#====================================================#
//# Last update:
//
// 15 Mar 2010: - take out redundant weight parameter in fillHistoDataAndMC().
//              - revise. review QCD plot making.
// 14 Mar 2010: - change histo array to vector. Add BookHistograms() suite.
//              - Revert all short to int to optimize code.
//              - remove doValidaion(), replaced with m_doValidation.
//              - add switch StudyZveto (D=on).
//
// 12 Mar 2010: fix small bug.
// 11 Mar 2010: replace 3D arrays (e_plus_jet & e_plus_jets_weighted) with 
//              private vectors. Simplify and imporve codes.
//
//  9 Mar 2010: - Adapt to accommodate tt0j-tt4j alpgen signal samples. 
//              - Add methods: DefineCrossSectionAlpgen7TeV, DrawTTnjetTable.
//              - Consistent use of GetNinit() & GetCrossSection().
//              - Remove old logs, keep only this year's.
//
//  3 Mar 2010: add DPhiMetJet, genMET plots.
// 25 Feb 2010: add Data_njetsVcuts TH2.
// 23 Feb 2010: fix nstage (16->13).
// 22 Feb 2010: add barrel/endcap plots for iso-met scatter plot (in NES).
// 17 Feb 2010: fix btag plot. Replace nGoodJet & Njet() with m_nGoodJet.
// 16 Feb 2010: Add SC6 and KT7. Add btag plots. Take out HT cut. Fix ("nfile =" -> "nfile +=").
//
// 15 Feb 2010: - MCTruthMAtch, change iniital value of ii from -1 to 0. 
//              - Small update in DrawSignalBGTable. Add global const ncutshown.
//
// 12 Feb 2010: - Add PDF weights. Add method to set HLT bit from script. Clean up.
//              - Remove unused mc_electron and mc_muons from branch list.
//              - Remove debug(),jetAlgo() etc (replaced with m_debug etc).
//              - Add method CheckAvailableJetMET().
//              - Update 7 TeV cross section.
//        
// 03 Feb 2010: Add function to check MC truth of reco particle MCTruthMatch(double, double)
//              and extended conversion study table to include charged pions
//
// 29 Jan 2010: Added reliso histo for trial control region.
// 26 Jan 2010: Added met plot for barrel-ele only.
//
// 22 Jan 2010: Added option ApplyMETcut() and RejectEndcapEle(), to replace SetRunPlanB().
// 20 Jan 2010: minor update.
// 19 Jan 2010: fix single top cross sections at 7 TeV.
// 15 Jan 2010: a) Added single-top xsec at 7 TeV. NB: s-chan NLO value not available.
//              b) Moved Init() into header.
//              c) Changed m3 histo to 0-960 & 960 bins; Added a copy of m3 with 0-1000 & 1000-bins.
//              d) Added skim eff for my summer09-7TeV-madgraph HLT skims.
//         
// 06 Jan 2010: a) Updated conversion algorithm routine. Separated out conversion study section in a separate function
//              b) Added in an additional electron ID option of "none" i.e. no ID applied
// 
//#====================================================#
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <time.h>

using namespace std;

#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChainElement.h"
#include "TGraphAsymmErrors.h"

// needed by m3 fit
#include "TRandom.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooHistPdf.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooConstVar.h"
using namespace RooFit;

#include "cms_bristol_ana_v3.hh" // defines ana class, including branches and leaves for tree


// Global variables/constants
//-----------------------------
const int nstage(13);
const int ntjet(5);
//const int nmctype(23); //extend to include wj, zj, QCD, VQQ, single top
const int nmctype(23+5); //extend to include tt0j-tt4j
const int nclass(16); //data+mc (incl vqq)
int ntype(1);//renamed from nhisto //FIXME

const int myprec(1); //no of decimal point for weighted nEvent
const int ncutshown(13); //11:BARREL, 12:1BTag, 13:2BTag (incl 4j)
const int nbm3 = 960;
const bool m3_use_1000_bins = false;
const bool use_old_Z_veto = false; //TEMPORARY

const string jetname[7]  = {"0j","1j","2j","3j","4j", "4mj", "allj"};
const string jetlabel[7] = {"0j","1j","2j","3j","4j", ">=4j","allj"};

const string mcname[16]  = { "data", "ttbar", "QCD", "enri1", "enri2" ,"enri3", "bce1","bce2","bce3",
			     "wj", "zj","vqq", "singleTop","tW","tchan","schan" };
const string mclabel[16] = { "data", "signal","QCD","enri1","enri2","enri3","bce1","bce2","bce3",
			     "W+jets","Z+jets","VQQ", "singleTop","tW","t-chan","s-chan" };
const string Fourjets = "$\\ge$4JETS"; //used in table





void ana::SetInputFile(const char* fname) {
  
  //cout << "start of SetInputFile" << endl;  

  // check the cut values before reading input
  static bool first_time=true;
  if(first_time) PrintCuts();

  nfile += chain->Add(fname);
  if( GetTrigger() )  chain2->Add(fname);
  //cout << "end of SetInputFile" << endl;
  first_time = false;
}

//------------------------------------------------------------------------------------

// Output file name: <first name>_<second name>.root/log
// User set in script.C the front bit of output file name ("first name"), the hist and text 
// output files will be set automatically according to type of input files:
//   data:  data
//   MC:    mc_[ttbar,wjets, bce1,2,3, enri1,2,3] (individual sample) or mc_mixture (mixture)
//
void ana::SetOutputFirstName(const string name) {

  //cout << "start of SetOutputFirstName" << endl;
  cout << "Total number in chain: "<< chain->GetEntries() << " events\n\n";

  string secondname = "noname";

  if (IsData()) {

    secondname = "data";

  } else {

    // ttnj, Alpgen (9-Mar-2010)
    mc_names.push_back("tt0j");
    mc_names.push_back("tt1j");
    mc_names.push_back("tt2j");
    mc_names.push_back("tt3j");
    mc_names.push_back("tt4j");

    mc_names.push_back("ttbar");
    mc_names.push_back("ttjet");
    mc_names.push_back("wjet");
    mc_names.push_back("wenu");
    mc_names.push_back("zjet");
    mc_names.push_back("zee");
    mc_names.push_back("enri1");
    mc_names.push_back("enri2");
    mc_names.push_back("enri3");
    mc_names.push_back("bce1");
    mc_names.push_back("bce2");
    mc_names.push_back("bce3");
    mc_names.push_back("vqq");
    mc_names.push_back("tW");
    mc_names.push_back("tchan");
    mc_names.push_back("schan");


    vector<int> nfiles(mc_names.size());
    vector<long> nevents(mc_names.size());

    cout << "List of input files:"<< endl;

    // loop over all files in chain   
    TObjArray *fileElements = chain->GetListOfFiles();
    TIter next(fileElements);
    TChainElement *chEl=0;

    map<string,bool> mc_seen;

    while (( chEl=(TChainElement*)next() )) {

      const string fname( chEl->GetTitle() );
      cout << fname << ":  "<< chEl->GetEntries() << " events";

      // try to find mc type i
      size_t i=0;
      while ( i<mc_names.size() ) {
	const string mc = mc_names.at(i);
	if ( fname.find( mc ) != string::npos ) { //match
	  //cout << "  match i="<< i<< "  mc_names=" << mc << endl;	
	  cout << "  " << mc << endl;
	  nfiles[i]++;
	  nevents[i] += chEl->GetEntries();
	  // if first time seeing this mc, then add this mc to map.
	  if ( mc_seen[mc] == 0 ) mc_seen[mc] = true;
	  break;
	}
	++i;
      }
    }
    //  cout << "mc_seen.size():" << mc_seen.size() << endl;

    //store the event count in private variable with key=mc
    for(size_t i=0; i<mc_names.size(); ++i){
      nInitialEventMC[mc_names.at(i)] = nevents.at(i);
    }

    // Set event weight map, and MC flag
    DefineCrossSection();
    SetEventWeightMap();
    SetMCFlag();

    // print out
    cout << "\n--------------------------------------------------------------------------------------------\n";
    cout << setw(10) << "mc"
	 << setw(14) << "N(files)"
	 << setw(18) << "N(events)"  
	 << setw(12) << "weight"
	 << setw(15) << "xsec"
	 << setw(14) << "Nexp (L=" << intlumi<< "/pb)";
    cout << "\n--------------------------------------------------------------------------------------------\n";
    for(size_t i=0; i<mc_names.size(); ++i){

      cout << setw(10) << mc_names.at(i) 
	   << setw(8) << nfiles.at(i) << " files"
	   << setw(11) << nevents.at(i) << " events";
      cout << setw(11);

      if (nfiles.at(i)==0) cout << "-" ;
      else cout << setprecision(4) << GetWeight(mc_names.at(i));

      // print out cross-sections
      cout << setw(14) << GetCrossSection(mc_names.at(i)) << " pb";
      cout << setw(14) << GetCrossSection(mc_names.at(i))*intlumi << endl;      
      cout << setprecision(6);
    }
    cout << "--------------------------------------------------------------------------------------------\n";

    cout << "\n MC: Listing types of MC found in input\n";
    for ( map<string,bool>::iterator i=mc_seen.begin(); i!=mc_seen.end(); ++i ){
      cout << "  " << i->first << endl;
    }

    if (mc_seen.size()==1) secondname = "mc_" + mc_seen.begin()->first; //one MC only
    else if (mc_seen.size()>1) secondname="mc_mixture";  //there are more than 1 kind of MC in input
  }// data or MC

  const string outname1 = name + "_" + secondname + ".root";
  const string outname2 = name + "_" + secondname + ".log";
  cout << "\n Analysis outputs are: " << outname1 << endl
       << "                       " << outname2 << endl << endl;
  SetOutputHistFile( outname1 );
  SetOutputTextFile( outname2 );

}//end SetOutputFirstName

void ana::SetOutputHistFile( const string name, const string mode) {

  histf = new TFile( name.c_str(), mode.c_str(), "test file" );
  outputHistFileName = name;
}

void ana::SetOutputTextFile(const string name, const string mode) {

  outfile = fopen( name.c_str(), mode.c_str() );
  outputTextFileName = name;
}

void ana::CheckTrigger(bool val) {
  checkTrig = val;
}

void ana::CheckTrigger(bool val, const string hlt) {
  checkTrig = val;
  HLTBit = hlt;
}

void ana::SetLimit(int val) {
  if(val>0) {
    numlimit = val;
  }else{
    numlimit = -1;
  }
}

//set kinematics cuts
void ana::SetEleETcut(float cut){
  ELE_ETCUT = cut;
  nCutSetInScript++;
}
void ana::SetMuonPTcut(float cut){
  MU_PTCUT = cut;
  nCutSetInScript++;
}
void ana::SetJetETcut(float cut){
  JET_PTCUT = cut;
  nCutSetInScript++;
}
void ana::SetMETcut(float cut){
  METCUT = cut;
  nCutSetInScript++;
}
void ana::SetHTcut(float cut){
  HTCUT = cut;
  nCutSetInScript++;
}

void ana::PrintCuts() const {
  
  if(IsData()) cout << " This run uses   REAL   Data\n\n";
  else         cout << " This run uses   MONTE CARLO   Data\n\n";
  cout << "***********************************************\n\n";
  cout << "  LHC c.o.m Energy :  "<< m_LHCEnergyInTeV << " TeV\n" << endl;
  cout << "***********************************************\n" << endl;
  if(checkTrig) cout << "  Trigger required :   " << HLTBit << endl;
  else          cout << "  Trigger not required" << endl;
  cout << "\n***********************************************" << endl;
  if(nCutSetInScript<4) {
    cout << "\n  WARNING! WARNING! YOU HAVE NOT SET ALL 4 CUTS!!!" << endl;
    cout << "\n  USING DEFAULT CUTS!!!" << endl;
  } else {
    cout << "\n Okay, all 4 cuts are set." << endl;
  }
  cout << "\n This run will use the following cut values:\n"<< endl;
  cout << " ELECTRON et cut =  " << ELE_ETCUT << "  GeV" << endl;
  cout << "     MUON pt cut =  " << MU_PTCUT  << "  GeV" << endl;
  cout << "      JET pt cut =  " << JET_PTCUT << "  GeV" << endl;
  cout << "      MET    cut =  " << METCUT    << "  GeV" << endl;
  //  cout << " HT cut     = " << HTCUT  << " GeV" << endl;
  cout << " ELECTRON ID requirement  =  " << printEleID() << endl;
  cout << "\n***********************************************" << endl;
  cout << "\n  HT cut in AES:  " << AES_HT_cut << " GeV";
  cout << "\n MET cut in AES:  " << AES_MET_cut << " GeV" << endl;
  if( AES_useSimpleZveto ) {
    cout << "\n  In AES, exclude events with 2 reco-electrons" << endl;
  }
  cout << "\n***********************************************" << endl;
  cout << "\n Jet collection:  " << m_jetAlgo << endl;
  cout << "\n MET collection:  " << m_metAlgo << endl;
  cout << "\n***********************************************" << endl;
  cout << "\n use_old_Z_veto : "<<  use_old_Z_veto << endl;
  cout << "\n***********************************************" << endl;
  if( !m_applyMETcut && m_rejectEndcapEle ) {
    cout << "\n  Run with Plan B:  MET is not used, add Eta(e) <1.442 as last cut" << endl;
    cout << "\n***********************************************" << endl;
  }
  if( m_applyMETcut && m_rejectEndcapEle ) {
    cout << "\n  Run with Plan C:  Both MET and Eta(e) <1.442 cuts are used" << endl;
    cout << "\n***********************************************" << endl;
  }

}//PrintCuts
//-------------------------------------------------------------------------------------------

long ana::GetNinit(const string mc) const {

  // if found this MC
  map<string,long>::const_iterator iter = nInitialEventMC.find( mc );
  if( iter != nInitialEventMC.end() )  // this MC found in input
    return iter->second;
  else
    return 0;
}
//-------------------------------------------------------------------------------------------

void ana::SetMCFlag(){
  
  //  if(nInitialEventMC["ttbar"]>0 || nInitialEventMC["ttjet"]>0) mc_sample_has_ttbar = true;
  if((GetNinit("ttbar")+GetNinit("ttjet"))>0) mc_sample_has_ttbar = true;

  // ttnj alpgen
  if( (GetNinit("tt0j")+GetNinit("tt1j")+GetNinit("tt2j")+GetNinit("tt3j")+GetNinit("tt4j")) >0 )
    mc_sample_has_ttbar = true;

  if( (GetNinit("wjet") + GetNinit("wenu")) > 0 )  mc_sample_has_Wjet  = true;
  if( (GetNinit("zjet") + GetNinit("zee"))  > 0 )  mc_sample_has_Zjet  = true;
  if((GetNinit("bce1")+ GetNinit("bce2")+ GetNinit("bce3")+
      GetNinit("enri1")+GetNinit("enri2")+GetNinit("enri3"))>0) mc_sample_has_QCD = true;
  if(GetNinit("enri1")>0)  mc_sample_has_enri1 = true;
  if(GetNinit("enri2")>0)  mc_sample_has_enri2 = true;
  if(GetNinit("enri3")>0)  mc_sample_has_enri3 = true;
  if(GetNinit("bce1")>0)   mc_sample_has_bce1  = true;
  if(GetNinit("bce2")>0)   mc_sample_has_bce2  = true;
  if(GetNinit("bce3")>0)   mc_sample_has_bce3  = true;
  if(GetNinit("vqq")>0)    mc_sample_has_VQQ   = true;
  if( (GetNinit("tW") + GetNinit("tchan") + GetNinit("schan")) > 0 )  mc_sample_has_singleTop = true;
  if(GetNinit("tW")>0)     mc_sample_has_tW    = true;
  if(GetNinit("tchan")>0)  mc_sample_has_tchan = true;
  if(GetNinit("schan")>0)  mc_sample_has_schan = true;
}
//-------------------------------------------------------------------------------------------

void ana::DefineCrossSection(){
  
  if ( m_LHCEnergyInTeV==10 ) { //10TeV = default

    cout << "\nNOTE: 10 TeV cross section not to date!!!" << endl;
    cout << "\n check https://twiki.cern.ch/twiki/bin/viewauth/CMS/CrossSections_3XSeries\n" << endl;

    cross_section["ttbar"] =          414. ;  //xs 414 (NLO+NLL at 10TeV)   241.7 pb (LO)
    cross_section["ttjet"] =          414. ;

    cross_section["wjet"]  =   (40e3*1.14) ;  //xs 40 nb  (K-factor=1.14)
    cross_section["zjet"]  =  (3.7e3*1.14) ;  //xs 3.7 nb (K-factor=1.14)
    
    cross_section["wenu"]  =  (11840*0.738*1.14) ; //xs 11840 pb (eff=0.738) (K-factor=1.14) // TO-BE-CHECKED
    cross_section["zee"]   =  (1944*1.14) ;        //xs 1944 pb (K-factor=1.14)             // TO-BE-CHECKED
    
    cross_section["enri1"] = 0.40e9 * 0.008 ;  //xs 0.4 mb (*filter efficiency=0.008)
    cross_section["enri2"] = 0.10e9 * 0.047 ;  //xs 0.1 mb
    cross_section["enri3"] =  1.9e6 * 0.15  ;  //xs 1.9e-3 mb
  
    cross_section["bce1"]  = 0.40e9 * 4.8e-4;  //xs 0.4 mb
    cross_section["bce2"]  = 0.10e9 * 2.4e-3;  //xs 0.1 mb
    cross_section["bce3"]  =  1.9e6 * 0.012 ;  //xs 1.9e-3 mb
    
    cross_section["vqq"]   =           290. ;   //xs 290 pb (W->lep)(LO)
    cross_section["tW"]    =            29. ;   //xs  29 pb (NLO) inclusive t,W decay
    cross_section["tchan"] =  130.0 * 0.324 ;   //xs 130 pb (NLO) * BR(t->blnu)
    cross_section["schan"] =    5.0 * 0.324 ;   //xs   5 pb (NLO) * BR(t->blnu)
  } 
  else if (m_LHCEnergyInTeV==7) { //7TeV

    cross_section["ttbar"] =          165. ;  //was 187 pb
    cross_section["ttjet"] =          165. ;

    cross_section["wjet"]  =  28000.; //was 24e3 * 1.14
    cross_section["zjet"]  =   2800.; //was 2.3e3 * 1.14
    
    cross_section["wenu"]  =  7899 * 0.779 * 1.14 ; //xs 7899 pb (pythia LO) (eff=0.779) (K-fac=1.14)
    cross_section["zee"]   =  1300 * 1.14 ;        //xs 1300 pb (pythia LO) (K-fac=1.14)
    
    cross_section["enri1"] = 0.2355e9 * 0.0073 ;  //xs 0.2355 mb (filter efficiency=0.0073)
    cross_section["enri2"] = 0.0593e9 * 0.059  ;  //xs 0.0593 mb
    cross_section["enri3"] =  0.906e6 * 0.148  ;  //xs 0.906e-3 mb
  
    cross_section["bce1"]  = 0.2355e9 * 0.00046;  //xs 0.2355 mb (filter efficiency=0.00046)
    cross_section["bce2"]  = 0.0593e9 * 0.00234;  //xs 0.0593 mb
    cross_section["bce3"]  =  0.906e6 * 0.0104 ;  //xs 0.906e-3 mb
    
    //cross_section["tW"]    =  11.0  ;   //xs  11 pb (NLO MCFM) inclusive t,W decay
    //cross_section["tchan"] =  20.7  ;   //xs  64 pb (NLO MCFM) * 0.324 (Br(t->blnu)) ~ 20.7
    //cross_section["schan"] =   0.99 ;   //0.99 pb is Madgraph value incl. BR(B->e,mu,tau), take from 'ProductionSummer2009at7TeV' Twiki

    cross_section["tW"]    =  10.6  ;         //xs  11 pb (NLO MCFM) inclusive t,W decay
    cross_section["tchan"] =   63. * 0.324;   //xs  63 pb (NLO MCFM) * 0.324 (Br(t->blnu)) = 20.412
    cross_section["schan"] =   4.6 * 0.324 ;  //4.6 pb x 0.324 (4.6pb is NNNLO) = 1.4904

    DefineCrossSectionAlpgen7TeV();

  } else {
    if(!IsData()) cout << "WARNING: Cross section values are not defined!" << endl;
  }

}//end DefineCrossSection
//-------------------------------------------------------------------------------------------

void ana::DefineCrossSectionAlpgen7TeV(){

  //--------------------
  // ttnj alpgen (7TeV)
  //--------------------

  if(signal_is_Alpgen) {

    cout << "\n Running on Alpgen signal"<< endl;
    cout << "NOTE: at the moment the ttbar decay table is disabled when running on Alpgen signal," << endl;
    cout << "      because mctype of 23-27 is used to identify Alpgen tt0j-tt4j, and it clashes" << endl;
    cout << "      with mctype of 1-10 corresponding to the ttbar decay modes." << endl;
    cout << "      Also, the ttbar column in S/B and MC tables is not filled." <<endl;
    cout << "\n\n  **  Alpgen Et " << signal_Alpgen_matching_threshold << "  **\n\n";
      
    vector<double> Alpgen_xsec;
  
    if (signal_Alpgen_matching_threshold==30) {     // Alpgen, Et30
      Alpgen_xsec.push_back(  79.5  ); //tt0j
      Alpgen_xsec.push_back(  58.5  ); //tt1j
      Alpgen_xsec.push_back(  27.0  ); //tt2j
      Alpgen_xsec.push_back(   9.17 ); //tt3j
      Alpgen_xsec.push_back(   2.52 ); //tt4j
    }
    else if (signal_Alpgen_matching_threshold==40) {   // Alpgen, Et40
      Alpgen_xsec.push_back(  79.6   );  //tt0j
      Alpgen_xsec.push_back(  42.6   );  //tt1j
      Alpgen_xsec.push_back(  14.8   );  //tt2j
      Alpgen_xsec.push_back(   3.87  );  //tt3j
      Alpgen_xsec.push_back(   0.828 );  //tt4j
    }
    else if (signal_Alpgen_matching_threshold==50) {    // Alpgen, Et50
      Alpgen_xsec.push_back(  79.5   );  //tt0j
      Alpgen_xsec.push_back(  32.5   );  //tt1j
      Alpgen_xsec.push_back(   8.92  );  //tt2j
      Alpgen_xsec.push_back(   1.88  );  //tt3j
      Alpgen_xsec.push_back(   0.327 );  //tt4j
    }

    // compute the total LO cross section
    double total_alpgen_signal_xsec = 0;
    for(size_t i=0; i<Alpgen_xsec.size(); i++){
      cout << "LO xsec( tt" << i << "j ):  " << Alpgen_xsec.at(i) << " pb" << endl;
      total_alpgen_signal_xsec += Alpgen_xsec.at(i);
    }
    cout << endl;
    // normalize cross section to NLO
    const double k_factor_alpgen = 165.0 / total_alpgen_signal_xsec;
      
    cross_section["tt0j"] = Alpgen_xsec.at(0) * k_factor_alpgen;
    cross_section["tt1j"] = Alpgen_xsec.at(1) * k_factor_alpgen;
    cross_section["tt2j"] = Alpgen_xsec.at(2) * k_factor_alpgen;
    cross_section["tt3j"] = Alpgen_xsec.at(3) * k_factor_alpgen;
    cross_section["tt4j"] = Alpgen_xsec.at(4) * k_factor_alpgen;
  }//ttnj alpgen

}// DefineCrossSectionAlpgen7TeV
//-------------------------------------------------------------------------------------------

double ana::GetCrossSection( const string mc ) const {
  return cross_section.find(mc)->second;
}
//-------------------------------------------------------------------------------------------

// Declare the event weights
void ana::SetEventWeightMap(){ //only if run on MC
  
   weightMap["data"]  = 1.0;

   for(size_t i=0 ; i<mc_names.size(); ++i){
     long n = GetNinit( mc_names.at(i) ) ;
     if( n>0 )
       weightMap[ mc_names.at(i) ] = GetCrossSection( mc_names.at(i) ) * intlumi / n;
     else
       weightMap[ mc_names.at(i) ] = 0;
   }
   
   //OctX for SD
   if( m_runOnSD ){
     if( m_LHCEnergyInTeV==10 ){
       cout << "\nSummer09 10 TeV, OctX SD skim efficiencies x prescale factor:" << endl;
       cout << "(note: nexp = Ninit * w / skim eff)" << endl;
       // rescale weight by filter efficiency (N_SD/N_ori), and prescale factor
       //                                  N_SD / N_ori     * pres
       const double skimEff_ttbar  =       3270 /   529750. * 100;
       const double skimEff_wenu   =      87514 /  2142960. * 20;
       const double skimEff_zee    =      15092 /  2682355. * 100;
       const double skimEff_enri1  =    6242601 / 33638282. ;
       const double skimEff_enri2  =   10792435 / 38360886. ;
       const double skimEff_enri3  =    2538711 /  5729547. ;
       const double skimEff_bce1   =     393019 /  2383833. ;
       const double skimEff_bce2   =     769808 /  2035108. ;
       const double skimEff_bce3   =     641522 /  1038080. ;

       weightMap["ttbar"] =  weightMap["ttbar"] * skimEff_ttbar;
       weightMap["wenu"]  =  weightMap["wenu"]  * skimEff_wenu ;
       weightMap["zee"]   =  weightMap["zee"]   * skimEff_zee  ;
       weightMap["enri1"] =  weightMap["enri1"] * skimEff_enri1;
       weightMap["enri2"] =  weightMap["enri2"] * skimEff_enri2;
       weightMap["enri3"] =  weightMap["enri3"] * skimEff_enri3;
       weightMap["bce1"]  =  weightMap["bce1"]  * skimEff_bce1 ;
       weightMap["bce2"]  =  weightMap["bce2"]  * skimEff_bce2 ;
       weightMap["bce3"]  =  weightMap["bce3"]  * skimEff_bce3 ;

       cout << "  skim eff   ttbar    " << skimEff_ttbar << endl;
       cout << "  skim eff   wenu     " << skimEff_wenu  << endl;
       cout << "  skim eff   zee      " << skimEff_zee   << endl;
       cout << "  skim eff   enri1    " << skimEff_enri1 << endl;
       cout << "  skim eff   enri2    " << skimEff_enri2 << endl;
       cout << "  skim eff   enri3    " << skimEff_enri3 << endl;
       cout << "  skim eff   bce1     " << skimEff_bce1  << endl;
       cout << "  skim eff   bce2     " << skimEff_bce2  << endl;
       cout << "  skim eff   bce3     " << skimEff_bce3  << endl;
       cout << endl;
     }
     else if (m_LHCEnergyInTeV==7 ){
       cout << "\nSummer09 7 TeV, OctX SD skim efficiencies x prescale factor:" << endl;
       cout << "(note: nexp = Ninit * w / skim eff)" << endl;
       // rescale weight by filter efficiency (N_SD/N_ori), and prescale factor
       //                                 N_SD / N_ori     * pres
       const double skimEff_ttbar  =      3841 /   626610. * 100 ;
       const double skimEff_wenu   =     85505 /  2078361. * 20  ;
       const double skimEff_zee    =     14683 /  2538855. * 100 ;
       const double skimEff_enri1  =   6169999 / 33505929. ;
       const double skimEff_enri2  =   9054696 / 32168675. ;
       const double skimEff_enri3  =   2463429 /  5551386. ;
       const double skimEff_bce1   =    432380 /  2752942. ;
       const double skimEff_bce2   =    840100 /  2261916. ; 
       const double skimEff_bce3   =    682720 /  1097829. ;

       weightMap["ttbar"] =  weightMap["ttbar"] * skimEff_ttbar; 
       weightMap["wenu"]  =  weightMap["wenu"]  * skimEff_wenu; 
       weightMap["zee"]   =  weightMap["zee"]   * skimEff_zee; 
       weightMap["enri1"] =  weightMap["enri1"] * skimEff_enri1; 
       weightMap["enri2"] =  weightMap["enri2"] * skimEff_enri2; 
       weightMap["enri3"] =  weightMap["enri3"] * skimEff_enri3;
       weightMap["bce1"]  =  weightMap["bce1"]  * skimEff_bce1;
       weightMap["bce2"]  =  weightMap["bce2"]  * skimEff_bce2; 
       weightMap["bce3"]  =  weightMap["bce3"]  * skimEff_bce3; 

       cout << "  skim eff   ttbar    " << skimEff_ttbar << endl;
       cout << "  skim eff   wenu     " << skimEff_wenu  << endl;
       cout << "  skim eff   zee      " << skimEff_zee   << endl;
       cout << "  skim eff   enri1    " << skimEff_enri1 << endl;
       cout << "  skim eff   enri2    " << skimEff_enri2 << endl;
       cout << "  skim eff   enri3    " << skimEff_enri3 << endl;
       cout << "  skim eff   bce1     " << skimEff_bce1  << endl;
       cout << "  skim eff   bce2     " << skimEff_bce2  << endl;
       cout << "  skim eff   bce3     " << skimEff_bce3  << endl;
       cout << endl;
     }
   }

   // for summer09 7TeV madgraph HLTskim
   if( m_runOnMyHLTskim ) {
     cout << "\nSummer09 7 TeV, Madgraph my HLT skim efficiency:" << endl;
     cout << "(note: nexp = Ninit * w / skim eff)" << endl;
     ///                         N_skim / N_ori    * pres
     const double skimEff_ttj =  610804 /  983964. ;
     const double skimEff_wj  = 2081537 / 8109289. ;
     const double skimEff_zj  =  381190 / 1068735. ;

     weightMap["ttjet"] =  weightMap["ttjet"] * skimEff_ttj;
     weightMap["wjet"]  =  weightMap["wjet"]  * skimEff_wj; 
     weightMap["zjet"]  =  weightMap["zjet"]  * skimEff_zj;

     cout << "  skim eff    ttjet      " << skimEff_ttj << endl;
     cout << "  skim eff    wjet       " << skimEff_wj  << endl;
     cout << "  skim eff    zjet       " << skimEff_zj  << endl;
     cout << endl;
   }

}//End SetEventWeightMap
//-------------------------------------------------------------------------------------------
double ana::GetWeight(string mc) const {
  return weightMap.find(mc)->second;
}
//-------------------------------------------------------------------------------------------

ana::ana(){

   // Print out code compilation time and current time
   cout << "\nCompilation Date/Time = " << __DATE__ << "  " << __TIME__ << endl;
   TDatime now;
   cout << "Current local ";
   now.Print();
   cout << endl;

   cout << "\n***********************************************************************";
   cout << "\n*                                                                     *";
   cout << "\n*        C M S     B R I S T O L     T O P     A N A L Y S I S        *";
   cout << "\n*                                                                     *";
   cout << "\n***********************************************************************";
   cout << endl << endl; 

  
   chain = new TChain("configurableAnalysis/eventB");
   chain2 = new TChain("configurableAnalysis/eventV");


   // Initialize private variables
   datafile                  = true;
   keepgood                  = true;
   checkTrig                 = true; 
   HLTBit                    = "HLT_Ele15_LW_L1R";
   numlimit                  = -1;
   this_weight               = 1.0;
   m_nGoodJet                = 0;
   m_QCDest_reliso_bin_width = 0.01;
   m_doValidation            = false;
   m_studyZveto              = true;
   m_plotRelisoNES           = true;
   m_debug                   = false;
   m_ConversionStudies       = false;
   m_jetAlgo                 = "Default";
   m_metAlgo                 = "calojet_mujes";
   m_LHCEnergyInTeV          = 7.0; //Default is 7 TeV
   m_runOnSD                 = false;
   m_runOnMyHLTskim          = false;
   m_useMisslayers           = false;
   m_ntoy                    = 2;
   m_muonCutNum              = 0;

   // Default values of kinematic cuts
   ELE_ETCUT                = 30.0;
   MU_PTCUT                 = 20.0;
   JET_PTCUT                = 30.0;
   METCUT                   = 30.0;
   HTCUT                    = 0.0;
   m_applyMETcut            = true;
   m_rejectEndcapEle        = true;
   nCutSetInScript          = 0;
   AES_HT_cut               = 200.0;
   AES_MET_cut              = 15.0;
   AES_useSimpleZveto       = true;
   m_eID                    = robustTight; //enum

   // integrated luminosity
   intlumi                 = 20.0; //pb-1
   useNewReliso            = true;
   doSystematics           = "";
   sysSample               = "";
   m_studyPDFunc           = false;
   ScientificNotation      = true;

   // initialize private variables (MC study)
   mc_sample_has_ttbar = false;
   mc_sample_has_Wjet  = false;
   mc_sample_has_Zjet  = false;
   mc_sample_has_QCD   = false;
   mc_sample_has_enri1 = false;
   mc_sample_has_enri2 = false;
   mc_sample_has_enri3 = false;
   mc_sample_has_bce1  = false;
   mc_sample_has_bce2  = false;
   mc_sample_has_bce3  = false;
   mc_sample_has_VQQ   = false;
   mc_sample_has_singleTop = false;
   mc_sample_has_tW    = false;
   mc_sample_has_tchan = false;
   mc_sample_has_schan = false;
   isTTbar = false;
   isWjets = false;
   isZjets = false;
   isQCD   = false;
   isEnri1 = false;
   isEnri2 = false;
   isEnri3 = false;
   isBce1  = false;
   isBce2  = false;
   isBce3  = false;
   isVQQ   = false;
   isSingleTop = false;
   isTW    = false;
   isTchan = false;
   isSchan = false;
   signal_is_Alpgen = false;
   signal_Alpgen_matching_threshold = 40;

   //856
   ConversionCounter = 0;
   for(int k=0;k<23; ++k){
     for(int i=0;i<2; ++i){
       for(int j=0;j<6; ++j){
         ConversionArray[k][i][j] = 0;
       }
     }
   }
   // end 856

   mycounter = 0;
}

void ana::ReadSelectedBranches() const {

   chain->SetBranchStatus("*",0); //disable all branches
   chain->SetBranchStatus("run",1);
   chain->SetBranchStatus("event",1);
   chain->SetBranchStatus("lumiBlock",1);
   chain->SetBranchStatus("beamSpot_x",1); //beam spot
   chain->SetBranchStatus("beamSpot_y",1);
   chain->SetBranchStatus("Nels",1); //electrons
   chain->SetBranchStatus("els_px",1);
   chain->SetBranchStatus("els_py",1);
   chain->SetBranchStatus("els_pz",1);
   chain->SetBranchStatus("els_energy",1);
   chain->SetBranchStatus("els_et",1);
   chain->SetBranchStatus("els_eta",1);
   chain->SetBranchStatus("els_phi",1);//z study
   chain->SetBranchStatus("els_looseId",1);
   chain->SetBranchStatus("els_tightId",1);
   chain->SetBranchStatus("els_robustLooseId",1);
   chain->SetBranchStatus("els_robustTightId",1);
   chain->SetBranchStatus("els_sigmaIEtaIEta",1); //sigma i eta i eta
   chain->SetBranchStatus("els_hadOverEm",1);
   chain->SetBranchStatus("els_dEtaIn",1);
   chain->SetBranchStatus("els_dPhiIn",1);
   chain->SetBranchStatus("els_eOverPIn",1);
   chain->SetBranchStatus("els_tIso",1);
   chain->SetBranchStatus("els_dr04EcalRecHitSumEt",1);
   chain->SetBranchStatus("els_dr04HcalTowerSumEt",1);
   chain->SetBranchStatus("els_d0dum",1);
   chain->SetBranchStatus("els_vx",1);
   chain->SetBranchStatus("els_vy",1);
   chain->SetBranchStatus("els_vpx",1);
   chain->SetBranchStatus("els_vpy",1);
   chain->SetBranchStatus("els_closestCtfTrackRef",1);
   chain->SetBranchStatus("els_tk_pt",1);
   chain->SetBranchStatus("els_tk_phi",1);
   chain->SetBranchStatus("els_tk_eta",1);
   chain->SetBranchStatus("els_tk_charge",1);
   chain->SetBranchStatus("els_tk_theta",1);
   chain->SetBranchStatus("els_shFracInnerHits",1);
   if(m_useMisslayers) chain->SetBranchStatus("els_innerLayerMissingHits",1);
   chain->SetBranchStatus("Nmus",1); //muons
   chain->SetBranchStatus("mus_cm_px",1); //global muon
   chain->SetBranchStatus("mus_cm_py",1);
   chain->SetBranchStatus("mus_cm_pz",1);
   chain->SetBranchStatus("mus_energy",1);
   chain->SetBranchStatus("mus_cm_pt",1);
   chain->SetBranchStatus("mus_cm_eta",1);
   chain->SetBranchStatus("mus_cm_chi2",1);
   chain->SetBranchStatus("mus_cm_ndof",1);
   chain->SetBranchStatus("mus_cm_d0dum",1);
   chain->SetBranchStatus("mus_tk_vx",1);
   chain->SetBranchStatus("mus_tk_vy",1);
   chain->SetBranchStatus("mus_tk_px",1);
   chain->SetBranchStatus("mus_tk_py",1);
   chain->SetBranchStatus("mus_tkHits",1);
   chain->SetBranchStatus("mus_tIso",1);
   chain->SetBranchStatus("mus_cIso",1);
   chain->SetBranchStatus("mus_id_AllGlobalMuons",1); //new
   chain->SetBranchStatus("mus_ecalvetoDep",1); //new
   chain->SetBranchStatus("mus_hcalvetoDep",1); //new

   if(m_jetAlgo=="Default"){ //using default jet-met
     chain->SetBranchStatus("Njets",1); //jets
     chain->SetBranchStatus("jets_px",1);
     chain->SetBranchStatus("jets_py",1);
     chain->SetBranchStatus("jets_pz",1);
     chain->SetBranchStatus("jets_energy",1);
     chain->SetBranchStatus("jets_eta",1);
     chain->SetBranchStatus("jets_pt",1);
     chain->SetBranchStatus("jets_emf",1);
     chain->SetBranchStatus("jets_btag_TC_highPur",1);//btag
     chain->SetBranchStatus("jets_btag_TC_highEff",1);
     chain->SetBranchStatus("jets_btag_secVertex",1);

   } else if (m_jetAlgo=="pfjet") { //PFJet
     chain->SetBranchStatus(Form("NPFJets"),1); //jets
     chain->SetBranchStatus("PFJets_px",1);
     chain->SetBranchStatus("PFJets_py",1);
     chain->SetBranchStatus("PFJets_pz",1);
     chain->SetBranchStatus("PFJets_energy",1);
     chain->SetBranchStatus("PFJets_eta",1);
     chain->SetBranchStatus("PFJets_pt",1);

   } else { //if not using default jets
     chain->SetBranchStatus(Form("Njets%s",      m_jetAlgo.c_str()),1); //jets
     chain->SetBranchStatus(Form("jets%s_px",    m_jetAlgo.c_str()),1);
     chain->SetBranchStatus(Form("jets%s_py",    m_jetAlgo.c_str()),1);
     chain->SetBranchStatus(Form("jets%s_pz",    m_jetAlgo.c_str()),1);
     chain->SetBranchStatus(Form("jets%s_energy",m_jetAlgo.c_str()),1);
     chain->SetBranchStatus(Form("jets%s_eta",   m_jetAlgo.c_str()),1);
     chain->SetBranchStatus(Form("jets%s_pt",    m_jetAlgo.c_str()),1);
     chain->SetBranchStatus(Form("jets%s_emf",   m_jetAlgo.c_str()),1);
   } 

   // 30-10-09
   chain->SetBranchStatus("Nmets",1);
   chain->SetBranchStatus("mets_et_muonCor",1);
   chain->SetBranchStatus("mets_phi_muonCor",1);
   chain->SetBranchStatus("mets_et",1);
   chain->SetBranchStatus("mets_phi",1);
   chain->SetBranchStatus("mets_gen_et",1);
   chain->SetBranchStatus("mets_gen_phi",1);
   /*
   if(m_metAlgo=="Default") { //muCor caloMET
     chain->SetBranchStatus("Nmets",1);
     chain->SetBranchStatus("mets_et_muonCor",1);
     chain->SetBranchStatus("mets_phi_muonCor",1);

   } else if(m_metAlgo=="calomet_mujes") {
     chain->SetBranchStatus("Nmets",1);
     chain->SetBranchStatus("mets_et",1);
     chain->SetBranchStatus("mets_phi",1);
   */
   if (m_metAlgo=="tcmet"){  //tcMET
     chain->SetBranchStatus("Ntcmets",1);
     chain->SetBranchStatus("tcmets_et",1);
     chain->SetBranchStatus("tcmets_phi",1);

   } else if (m_metAlgo=="pfmet"){  //PFMET
     chain->SetBranchStatus("NPFMets",1);
     chain->SetBranchStatus("PFMets_et",1);
     chain->SetBranchStatus("PFMets_phi",1);
 
   } else if (m_metAlgo=="SC5" || m_metAlgo=="SC7" ||
	      m_metAlgo=="KT4" || m_metAlgo=="KT6") { //CaloMET: SC5/7, KT4/6
     chain->SetBranchStatus(Form("Nmets%s",   m_metAlgo.c_str()),1);
     chain->SetBranchStatus(Form("mets%s_et", m_metAlgo.c_str()),1);
     chain->SetBranchStatus(Form("mets%s_phi",m_metAlgo.c_str()),1);;
   }

   chain->SetBranchStatus("Ntracks",1); //tracks
   chain->SetBranchStatus("tracks_pt",1);
   chain->SetBranchStatus("tracks_phi",1);
   chain->SetBranchStatus("tracks_eta",1);
   chain->SetBranchStatus("tracks_chg",1);
   chain->SetBranchStatus("tracks_theta",1);
   chain->SetBranchStatus("tracks_vx",1);
   chain->SetBranchStatus("tracks_vy",1);
   chain->SetBranchStatus("tracks_px",1);
   chain->SetBranchStatus("tracks_py",1);
   if(m_useMisslayers) chain->SetBranchStatus("tracks_innerLayerMissingHits",1);

   chain->SetBranchStatus("Nphotons",1); //z study
   chain->SetBranchStatus("photons_eta",1); //z study
   chain->SetBranchStatus("photons_et",1); //z study
   chain->SetBranchStatus("photons_phi",1); //z study
   chain->SetBranchStatus("photons_px",1); //z study
   chain->SetBranchStatus("photons_py",1); //z study
   chain->SetBranchStatus("photons_pz",1); //z study
   chain->SetBranchStatus("photons_energy",1); //z study
   ///-------------------------- MC Truth info ------------------------------
   if ( !IsData() ) {
     chain->SetBranchStatus("Nmc_doc",1);   // generator particles
     chain->SetBranchStatus("mc_doc_id",1);
     chain->SetBranchStatus("mc_doc_status",1);
     chain->SetBranchStatus("mc_doc_mother_id",1);
     chain->SetBranchStatus("mc_doc_mass",1);
     chain->SetBranchStatus("mc_doc_px",1);
     chain->SetBranchStatus("mc_doc_py",1);
     chain->SetBranchStatus("mc_doc_pz",1);
     chain->SetBranchStatus("mc_doc_pt",1);
     chain->SetBranchStatus("mc_doc_eta",1);
     chain->SetBranchStatus("mc_doc_phi",1); //z study
     chain->SetBranchStatus("mc_doc_energy",1); //z study

     if( m_studyPDFunc ){
       for(int i=0; i<=44; i++){
	 string pdf = Form("PDFWcteq66_%u",i);//unsigned %u
	 if(m_debug) cout <<" pdf  " << pdf << endl;
	 chain->SetBranchStatus( pdf.c_str(), 1 );
       }
     }     
   }
   ///------------------------  MC Truth info (END)  ---------------------------

   if(GetTrigger()) {
     //chain->SetBranchStatus("HLT_Ele15_LW_L1R",1); //trigger (8e29)
     //chain->SetBranchStatus("HLT_Ele10_SW_L1R",1); //trigger (1e30)
     chain->SetBranchStatus(HLTBit.c_str(),1);
   }   
}//end ReadSelectedBranches
//------------------------------------------------------------------------------------

void ana::CheckAvailableJetMET(){
  
   vector<string> availableJET;
   vector<string> availableMET;

   // check if various jet/met is in the ntuple
   TObjArray *list = chain->GetListOfBranches(); 
   //int m = list->GetEntries(); 
   //cout << "Number of branches: " << m << endl;    

   // -- Loop over all, and draw their variables into TCanvas c1
   for ( int i = 0; i < list->GetEntries(); ++i) {
     string name = ((TBranch*)(*list)[i])->GetName() ;
     // jet
     if     (name=="NjetsSC5")    availableJET.push_back("SC5");
     else if(name=="NjetsSC7")    availableJET.push_back("SC7");
     else if(name=="NjetsKT4")    availableJET.push_back("KT4");
     else if(name=="NjetsKT6")    availableJET.push_back("KT6");
     else if(name=="NPFJets")     availableJET.push_back("PF(AK5)");
     else if(name=="NjetsJPTAK5") availableJET.push_back("JPT(AK5)");
     // met
     if     (name=="NmetsSC5")    availableMET.push_back("SC5");
     else if(name=="NmetsSC7")    availableMET.push_back("SC7");
     else if(name=="NmetsKT4")    availableMET.push_back("KT4");
     else if(name=="NmetsKT6")    availableMET.push_back("KT6");
     else if(name=="NPFMets")     availableMET.push_back("PF(AK5)");
     else if(name=="Ntcmets")     availableMET.push_back("tcmet");
   }
   cout << " Available jet: ";
   for(unsigned int i=0; i<availableJET.size(); ++i){
     cout  << " " << availableJET[i];
   }
   cout << endl;
   cout << " Available met: ";
   for(unsigned int i=0; i<availableMET.size(); ++i){
     cout << " " << availableMET[i];
   }
   cout << endl;
   availableJET.clear();
   availableMET.clear();

}// end CheckAvailableJetMET()
//------------------------------------------------------------------------------------





//====================================================================================
//
//                           Histogram   Booking
//
//====================================================================================
// 13 Mar 2010
void ana::BookHistograms(){
  if(m_debug) cout << "Starting BookHistograms()"<< endl;
  if(!IsData()) ntype = 16; //MC

  if(m_doValidation) BookHistograms_valid();
  BookHistograms_basicKin();
  BookHistograms_explore();
  BookHistograms_nEle();
  //  BookHistograms_nEle_eid();
  BookHistograms_ed0();
  if(m_studyZveto) BookHistograms_zVeto();
  BookHistograms_met();
  BookHistograms_HT();
  BookHistograms_mtw();
  BookHistograms_DRemu();
  BookHistograms_DPhiEmet();
  BookHistograms_DPhimetjet();
  BookHistograms_conv(); //ask Frankie if this can be removed
  BookHistograms_QCD();
  BookHistograms_QCD_planA();
  BookHistograms_QCD_planB();
  BookHistograms_wj();
  BookHistograms_event_table();
  BookHistograms_btag();
  BookHistograms_PDFunc();//to be tested

}//end BookHistograms()
//------------------------------------------------------------------------------------


void ana::BookHistograms_valid() {

   //==================================
   //
   //  Validation Plots
   //
   //==================================
   if(m_debug) cout << "Starting BookHistograms_valid()"<< endl;

   TDirectory *dir_Valid = histf->mkdir("Valid","Validation plots");
   dir_Valid->cd();
   valid_mkHisto_cut_njet(valid_HT,        "ht",       "ht",            1000, 0, 1000);
   valid_mkHisto_cut_njet(valid_jetsEt,    "jetsEt",   "all_jetsEt",    1000, 0 ,1000);
   valid_mkHisto_cut_njet(valid_jetsEta,   "jetsEta",  "all_jetsEta",   1000, -5, 5);
   valid_mkHisto_cut_njet(valid_jetsPhi,   "jetsPhi",  "all_jetsPhi",   1000, -5, 5);
   valid_mkHisto_cut_njet(valid_jets1stEt, "jets1st",  "clean_jets1st", 1000, 0, 1000);
   valid_mkHisto_cut_njet(valid_jets2ndEt, "jets2nd",  "clean_jets2nd", 1000, 0, 1000);
   valid_mkHisto_cut_njet(valid_jets3rdEt, "jets3rd",  "clean_jets3rd", 1000, 0, 1000);
   valid_mkHisto_cut_njet(valid_jets4thEt, "jets4th",  "clean_jets4th", 1000, 0, 1000);
   valid_mkHisto_cut_njet(valid_eleEt,     "eleEt",    "eleEt",         1000, 0, 1000);
   valid_mkHisto_cut_njet(valid_eleEta,    "eleEta",   "eleEta",        1000, -5, 5);
   valid_mkHisto_cut_njet(valid_elePhi,    "elePhi",   "elePhi",        1000, -5, 5);
   valid_mkHisto_cut_njet(valid_eleCalIso, "eleCalIso","eleCalIso",     1000, 0, 1);
   valid_mkHisto_cut_njet(valid_eleTrkIso, "eleTrkIso","eleTrkIso",     1000, 0, 1);
   valid_mkHisto_cut_njet(valid_eleRelIso, "eleRelIso","eleRelIso",     1000, 0, 1);
   valid_mkHisto_cut_njet(valid_eled0,     "eled0",    "eled0",         2000, 0, 2);
   valid_mkHisto_cut_njet(valid_metEt,     "metEt",    "metEt",         1000, 0, 1000);
   valid_mkHisto_cut_njet(valid_metPhi,    "metPhi",   "metPhi",        1000, -5, 5);
   valid_mkHisto_cut_njet(valid_genTT_pt,  "genTT_pt", "genTT_pt",      1000, 0, 1000);
   valid_mkHisto_cut_njet(valid_genT_pt,   "genT_pt",  "genT_pt",       1000, 0, 1000);
   valid_mkHisto_cut_njet(valid_recoM3,    "recoM3",   "recoM3",        1000, 0, 1000);
   valid_mkHisto_cut_njet(valid_mass_ee,   "mass_ee",  "mass_ee",       1000, 0, 500);
   valid_mkHisto_cut_njet(valid_recoM3_PTMax, "recoM3_PTMax","recoM3_PTMax",1000, 0, 1000);
   valid_mkHisto_cut_njet(valid_numberTracks, "numberTracks","numberTracks", 500, 0, 500); 
   valid_mkHisto_cut_njet(valid_trackPt, "  trackPt",  "trackPt",        500, 0,500); 

}//end BookHistograms_valid()
//------------------------------------------------------------------------------------


void ana::BookHistograms_basicKin(){

   //==================================
   //
   //  Basic Kinematic Plots (Aug 09)
   //
   //==================================
   if(m_debug) cout << "Starting BookHistograms_basicKin()"<< endl;

   TDirectory *dir_basic = histf->mkdir("basic","Basic kinematics");
   dir_basic->cd();

   //-------------
   //  Electrons
   //-------------
   h_ele_ET.resize(4);
   h_ele_eta.resize(4);
   h_ele_phi.resize(4);
   h_ele_iso.resize(4);
   addHistoDataAndMC( h_nele,       "nele",     "N(e) (no cut)",                   6,0,6);
   addHistoDataAndMC( h_ele_ET[0],  "ele_ET",   "E_{T}(all e) (nocut)",            50,0,100);
   addHistoDataAndMC( h_ele_ET[1],  "ele1_ET",  "E_{T}(leading e: no cut)",        50,0,100);
   addHistoDataAndMC( h_ele_ET[2],  "ele2_ET",  "E_{T}(2^{nd} leading e: no cut)", 50,0,100);
   addHistoDataAndMC( h_ele_ET[3],  "ele3_ET",  "E_{T}(3^{rd} leading e: no cut)", 50,0,100);
   addHistoDataAndMC( h_ele_eta[0], "ele_eta",  "#eta(all e) (no cut)",           50,-2.5,2.5);
   addHistoDataAndMC( h_ele_eta[1], "ele1_eta", "#eta(leading e: no cut)",        50,-2.5,2.5);
   addHistoDataAndMC( h_ele_eta[2], "ele2_eta", "#eta(2^{nd} leading e: no cut)", 50,-2.5,2.5);
   addHistoDataAndMC( h_ele_eta[3], "ele3_eta", "#eta(3^{rd} leading e: no cut)", 50,-2.5,2.5);
   addHistoDataAndMC( h_ele_phi[0], "ele_phi",  "#phi(all e) (no cut)",           50,-3.2,3.2);
   addHistoDataAndMC( h_ele_phi[1], "ele1_phi", "#phi(leading e: no cut)",        50,-3.2,3.2);
   addHistoDataAndMC( h_ele_phi[2], "ele2_phi", "#phi(2^{nd} leading e: no cut)", 50,-3.2,3.2);
   addHistoDataAndMC( h_ele_phi[3], "ele3_phi", "#phi(3^{rd} leading e: no cut)", 50,-3.2,3.2);
   addHistoDataAndMC( h_ele_iso[0], "ele_iso",  "RelIso(all e) (no cut)",           200,0,2);
   addHistoDataAndMC( h_ele_iso[1], "ele1_iso", "RelIso(leading e: no cut)",        200,0,2);
   addHistoDataAndMC( h_ele_iso[2], "ele2_iso", "RelIso(2^{nd} leading e: no cut)", 200,0,2);
   addHistoDataAndMC( h_ele_iso[3], "ele3_iso", "RelIso(3^{rd} leading e: no cut)", 200,0,2);   

   //-----------
   //  J E T S 
   //----------- 
   // 14-8-09: jets after jet-cleaning
   h_jet_PT.resize(5);
   h_jet_eta.resize(5);
   h_jet_phi.resize(5);
   addHistoDataAndMC( h_njet,       "njet",     "jet multiplicity (ET/eta)",                   10,0,10);
   addHistoDataAndMC( h_jet_PT[0],  "jet_PT",   "P_{T}(all jets) (pass P_{T}/#eta)",           100,0,200);
   addHistoDataAndMC( h_jet_PT[1],  "jet1_PT",  "P_{T}(leading jet) (pass P_{T}/#eta)",        100,0,200);
   addHistoDataAndMC( h_jet_PT[2],  "jet2_PT",  "P_{T}(2^{nd} leading jet) (pass P_{T}/#eta)", 100,0,200);
   addHistoDataAndMC( h_jet_PT[3],  "jet3_PT",  "P_{T}(3^{rd} leading jet) (pass P_{T}/#eta)", 100,0,200);
   addHistoDataAndMC( h_jet_PT[4],  "jet4_PT",  "P_{T}(4^{th} leading jet) (pass P_{T}/#eta)", 100,0,200);
   addHistoDataAndMC( h_jet_eta[0], "jet_eta",  "#eta(all jets) (pass P_{T}/#eta)",           50,-2.5,2.5);
   addHistoDataAndMC( h_jet_eta[1], "jet1_eta", "#eta(leading jet) (pass P_{T}/#eta)",        50,-2.5,2.5);
   addHistoDataAndMC( h_jet_eta[2], "jet2_eta", "#eta(2^{nd} leading jet) (pass P_{T}/#eta)", 50,-2.5,2.5);
   addHistoDataAndMC( h_jet_eta[3], "jet3_eta", "#eta(3^{rd} leading jet) (pass P_{T}/#eta)", 50,-2.5,2.5);
   addHistoDataAndMC( h_jet_eta[4], "jet4_eta", "#eta(4^{rd} leading jet) (pass P_{T}/#eta)", 50,-2.5,2.5);
   addHistoDataAndMC( h_jet_phi[0], "jet_phi",  "#phi(all jets) (pass P_{T}/#eta)",           50,-3.2,3.2);
   addHistoDataAndMC( h_jet_phi[1], "jet1_phi", "#phi(leading jet) (pass P_{T}/#eta)",        50,-3.2,3.2);
   addHistoDataAndMC( h_jet_phi[2], "jet2_phi", "#phi(2^{nd} leading jet) (pass P_{T}/#eta)", 50,-3.2,3.2);
   addHistoDataAndMC( h_jet_phi[3], "jet3_phi", "#phi(3^{rd} leading jet) (pass P_{T}/#eta)", 50,-3.2,3.2);
   addHistoDataAndMC( h_jet_phi[4], "jet4_phi", "#phi(4^{rd} leading jet) (pass P_{T}/#eta)", 50,-3.2,3.2);

   //-----
   // MET
   //-----
   addHistoDataAndMC( h_metAlone,     "metAlone",     "Missing ET (no cut)",  200, 0, 200);
   addHistoDataAndMC( h_metAlone_phi, "metAlone_phi", "MET #phi (no cut)",   50,-3.2,3.2);   
  
}//end BookHistogram_basic()
//--------------------------------------------------------------------------------------


void ana::BookHistograms_explore(){

   //=======================
   //
   // Explore new variables
   //
   //=======================
   if(m_debug) cout << "Starting BookHistograms_explore()"<< endl;

   // 4-11-09: explore var after all but MET cuts
   TDirectory *dir_explore = histf->mkdir("explore","explore new var after all but MET cuts (>=4j)");
   dir_explore->cd();
   addHistoDataAndMC( h_exp_ele_et,  "ele_et",  "ele E_{T}",   50, 20, 120 );
   addHistoDataAndMC( h_exp_ele_eta, "ele_eta", "ele #eta",    50, -2.5, 2.5 );
   addHistoDataAndMC( h_exp_j0_pt,   "j0_pt",   "leading jet p_{T}", 100, 0, 200 );
   addHistoDataAndMC( h_exp_j1_pt,   "j1_pt",   "2nd leading jet p_{T}", 100, 0, 200 );
   addHistoDataAndMC( h_exp_DRej,    "DRej",    "#DeltaR(e,j0)",     60, 0, 6 );
   addHistoDataAndMC( h_exp_DPhiej,  "DPhiej",  "#Delta#Phi(e,j0)",  64, -3.2, 3.2 );
   addHistoDataAndMC( h_exp_DRjj,    "DRjj",    "#DeltaR(j0,j1)",    60, 0, 6 );
   addHistoDataAndMC( h_exp_DPhijj,  "DPhijj",  "#Delta#Phi(j0,j1)", 64, -3.2, 3.2 );
   addHistoDataAndMC( h_exp_met_v_eeta,  "met_v_eeta", "met v #eta(e)", 50, 0, 100, 50, -2.5, 2.5 );

}//end BookHistograms_explore
//------------------------------------------------------------------------------------


void ana::BookHistograms_nEle() {

  //--------------------
  // Electron counting
  //--------------------
  if(m_debug) cout << "Starting BookHistograms_nEle()"<< endl;

  TDirectory *dir_nEle = histf->mkdir("electron_count");
  dir_nEle->cd();
      
  // count electron per event, for each njet and mctype
  addHisto_Njet_DataAndMC( h_nEle_all,        "nEle_all",              "N(e) all", 6,0,6);
  addHisto_Njet_DataAndMC( h_nEle_s1,         "nEle_s1_EtEta",         "N(e) pass ET,eta cuts", 6,0,6);
  addHisto_Njet_DataAndMC( h_nEle_s2,         "nEle_s2_d0",            "N(e) pass ET,eta,d0 cuts", 6,0,6);
  addHisto_Njet_DataAndMC( h_nEle_s3_idLoose, "nEle_s3_idLoose",       "N(e) pass ET,eta,d0,looseId cuts", 6,0,6);
  addHisto_Njet_DataAndMC( h_nEle_s3_idTight, "nEle_s3_idTight",       "N(e) pass ET,eta,d0,tightId cuts", 6,0,6);
  addHisto_Njet_DataAndMC( h_nEle_s3_idRL,    "nEle_s3_idRobustLoose", "N(e) pass ET,eta,d0,robustLooseId cuts", 6,0,6);
  addHisto_Njet_DataAndMC( h_nEle_s3_idRT,    "nEle_s3_idRobustTight", "N(e) pass ET,eta,d0,robustTightId cuts", 6,0,6);

  // sum of electron in all events
  addHistoDataAndMC( h_eid, "eid", "Total N(ele) pass each cut", 7,0,7);//all, ET+eta, d0, loose, tight, RL, RT
  SetHistoLabelEleID( h_eid );

}//End BookHistograms_nEle
//------------------------------------------------------------------------------------


void ana::BookHistograms_eid() {
   
  //----------------------------
  //  Electron ID variables
  //----------------------------
  if(m_debug) cout << "Starting BookHistograms_eid()"<< endl;

  TDirectory *dir_eid = histf->mkdir("electron_id");
  dir_eid->cd();

}
//------------------------------------------------------------------------------------


void ana::BookHistograms_ed0() {

  //----------------------------
  //  Electron d0
  //----------------------------
   if(m_debug) cout << "Starting BookHistograms_ed0()"<< endl;

   // Electron d0 plots (bin: 0.5cm for 5000 bins, ie 1 micron per bin)
   TDirectory *dir_ed0_unCor = histf->mkdir("electron_d0_unCor", "electron d0 (uncorrected)");
   dir_ed0_unCor->cd();
   addHistoDataAndMC( h_ed0_unCor, "ed0_unCor", "ele |d0| (pass ET,eta cut)", 3000, 0, 0.3);

   TDirectory *dir_ed0 = histf->mkdir("electron_d0","electron d0 (wrt BeamSpot)");
   dir_ed0->cd();
   addHistoDataAndMC( h_ed0,      "ed0",      "ele |d0| (pass ET,eta cut)",     3000, 0, 0.3);
   addHistoDataAndMC( h_ed0_pass, "ed0_pass", "ele |d0| (pass all cuts, >=4j)", 3000, 0, 0.3);
   histf->cd();

}//end BookHistograms_ed0()
//------------------------------------------------------------------------------------


void ana::BookHistograms_zVeto() {

   //---------------------
   //  Study of Z veto
   //---------------------
   if(m_debug) cout << "Starting BookHistograms_zVeto()"<< endl;
   TDirectory *dir_zveto = histf->mkdir("z_veto","Z veto");
   dir_zveto->cd();

   // This only for Z+jet MC
   if(IsData()==false && mc_sample_has_Zjet) {
     h_nGenBasicEle_Zee_allj = new TH1F("nGenBasicEle_Zee_allj","N(GenEle) (Z->ee, allj) pass ET,eta cut",6,0,6);
     h_nGenBasicEle_Zee_allj->SetDrawOption("hist text0");
     h_nGenBasicEle_Zee_allj->SetMarkerSize(3);
     h_Zee_eta        = new TH1F("Zee_eta","eta of gen-ele from Z", 100,-5,5);
     h_Zee_pt         = new TH1F("Zee_pt","pt of gen-ele from Z", 100,0,100);
     h_Z_photon_eta   = new TH1F("Z_photon_eta",  "eta of reco-photons in Z events", 100,-4,4);
     h_Z_photon_et    = new TH1F("Z_photon_et ",  "ET of reco-photons in Z events", 100,0,100);
     h_Zee_photon_eta = new TH1F("Zee_photon_eta","eta of reco-photons in Z->ee events", 100,-4,4);
     h_Zee_photon_et  = new TH1F("Zee_photon_et","ET of reco-photons in Z->ee events", 100,0,100);
     h_Zee_photon_eteta_2D  = new TH2F("Zee_photon_eteta_2D",
				       "ET vs eta, of reco-photons in Z->ee events", 
				       100,0,100,100,-4,4);
     h_Z_Nphotons   = new TH1F("Z_Nphotons",  "N(reco-photons) in Z events (e,mu,tau)", 6,0,6);
     h_Zee_Nphotons = new TH1F("Zee_Nphotons","N(reco-photons) in Z->ee events", 6,0,6);
   }

   TDirectory *dir_z1 = dir_zveto->mkdir("NES","Normal Event Selection");
   dir_z1->cd();
   addHistoDataAndMC( h_mass_diele,     "mass_diele",     "M(e,e) (sel,rloose)",     100,0,180);   
   addHistoDataAndMC( h_mass_diele_new, "mass_diele_new", "M(e,e) (sel,rloose) NEW", 100,0,180);   

   TDirectory *dir_z2 = dir_zveto->mkdir("AES","Anti Event Selection");
   dir_z2->cd();
   addHistoDataAndMC( h_mass_diele_lowMet_1j,  "mass_diele_lowMet_1j", "M(e,e) (1 rtight, 1 rloose) lowMET 1j", 100,0,180);   
   addHistoDataAndMC( h_mass_ephoton_lowMet_1j, "mass_ephoton_lowMet_1j", "M(e,#gamma) (1 rtight, 1 #gamma(nocut)) lowMET 1j", 100,0,180);
   addHistoDataAndMC( h_Nele_lowMet_1j,   "Nele_lowMet_1j","N(reco-ele) lowMET 1j", 6,0,6);
   addHistoDataAndMC( h_Nphoton_lowMet_1j,"Nphoton_lowMet_1j","N(reco-photons) lowMET 1j", 6,0,6);
   addHistoDataAndMC( h_photon_eta_lowMet_1j, "photon_eta_lowMet_1j", "#eta(reco-photons) lowMET 1j", 100,-4,4);
   addHistoDataAndMC( h_photon_et_lowMet_1j,  "photon_et_lowMet_1j",  "E_{T}(reco-photons) lowMET 1j", 100,0,100);
   addHistoDataAndMC( h_photon1_eta_lowMet_1j,"photon1_eta_lowMet_1j","#eta(leading reco-photon) lowMET 1j", 100,-4,4);
   addHistoDataAndMC( h_photon1_et_lowMet_1j, "photon1_et_lowMet_1j", "E_{T}(leading reco-photon) lowMET 1j", 100,0,100);

}// end BookHistograms_zVeto()
//------------------------------------------------------------------------------------


void ana::BookHistograms_met() {

  //---------------------------------
  //  MET
  //---------------------------------
   if(m_debug) cout << "Starting BookHistograms_met()"<< endl;
   TDirectory *dir_MET = histf->mkdir("MET","MET after all but MET cut");
   dir_MET->cd();
   addHisto_Njet_DataAndMC( h_met,    "met",    "MET (n-1)",         200, 0, 200);
   addHisto_Njet_DataAndMC( h_met_mu, "met_mu", "#muMET (n-1)",      200, 0, 200);
   addHisto_Njet_DataAndMC( h_met_t1, "met_t1", "Type1MET (n-1)",    200, 0, 200);
   addHisto_Njet_DataAndMC( h_met_BA, "met_BA", "MET (n-1 barrel)",  200, 0, 200);
   addHisto_Njet_DataAndMC( h_met_mu_BA, "met_mu_BA", "#muMET (n-1 barrel)",   200, 0, 200);
   addHisto_Njet_DataAndMC( h_met_t1_BA, "met_t1_BA", "Type1MET (n-1 barrel)", 200, 0, 200);

   addHistoDataAndMC( h_met_ante_ISO,    "met_ante_ISO",    "Missing ET (1 GoodEle, before ISO cut)", 200, 0, 200);
   addHistoDataAndMC( h_met_ante_ISO_mu, "met_ante_ISO_mu", "muMET (1 GoodEle, before ISO cut)",      200, 0, 200);
   addHistoDataAndMC( h_met_ante_ISO_t1, "met_ante_ISO_t1", "Type1MET (1 GoodEle, before ISO cut)",   200, 0, 200);

   // NEW: 3-3-10: GenMET
   TDirectory *dir_genMET = dir_MET->mkdir("gen","gen MET");
   dir_genMET->cd();
   addHisto_Njet_DataAndMC( h_met_gen,            "met_gen",            "genMET (n-1)",                            200,    0, 200);
   addHisto_Njet_DataAndMC( h_met_gen_diff_t1,    "met_gen_diff_t1",    "MET^{t1}_{reco}-MET_{gen} (n-1)",         100, -100, 100);
   addHisto_Njet_DataAndMC( h_met_gen_diff_mu,    "met_gen_diff_mu",    "MET^{#mu}_{reco}-MET_{gen} (n-1)",        100, -100, 100);
   addHisto_Njet_DataAndMC( h_met_gen_dphi_t1,    "met_gen_dphi_t1",    "#phi(t1met)-#phi(genmet) (n-1)",  64, -3.2, 3.2);
   addHisto_Njet_DataAndMC( h_met_gen_dphi_mu,    "met_gen_dphi_mu",    "#phi(mumet)-#phi(genmet) (n-1)",  64, -3.2, 3.2);
   addHisto_Njet_DataAndMC( h_met_gen_BA,         "met_gen_BA",         "genMET (n-1 barrel)",                     200,    0, 200);
   addHisto_Njet_DataAndMC( h_met_gen_diff_t1_BA, "met_gen_diff_t1_BA", "MET^{t1}_{reco}-MET_{gen} (n-1 barrel)",  100, -100, 100);
   addHisto_Njet_DataAndMC( h_met_gen_diff_mu_BA, "met_gen_diff_mu_BA", "MET^{#mu}_{reco}-MET_{gen} (n-1 barrel)", 100, -100, 100);
   addHisto_Njet_DataAndMC( h_met_gen_dphi_t1_BA, "met_gen_dphi_t1_BA", "#phi(t1met)-#phi(genmet) (n-1 barrel)",  64, -3.2, 3.2);
   addHisto_Njet_DataAndMC( h_met_gen_dphi_mu_BA, "met_gen_dphi_mu_BA", "#phi(mumet)-#phi(genmet) (n-1 barrel)",  64, -3.2, 3.2);

}//end BookHistograms_met()
//------------------------------------------------------------------------------------


void ana::BookHistograms_HT() {

   if(m_debug) cout << "Starting BookHistograms_HT()"<< endl;

   TDirectory *dir_HT = histf->mkdir("HT","HT");
   dir_HT->cd();
   addHisto_Njet_DataAndMC( h_HT,  "HT",  "HT (after all but HT cut)", 100, 0, 1000);
}
//------------------------------------------------------------------------------------


void ana::BookHistograms_mtw() {

   //----------------------------------
   // Transverse mass of W->enu, MTW
   //----------------------------------
   if(m_debug) cout << "Starting BookHistograms_mtw()"<< endl;
   TDirectory *dir_mtw = histf->mkdir("mtw","mT(W) W transverse mass");
   dir_mtw->cd();
   addHistoDataAndMC( h_mtw_mu_incl, "mtw_mu_incl", "mT(W) #mu-cor caloMet (incl, leading iso ele)",  100, 0, 200);
   addHistoDataAndMC( h_mtw_t1_incl, "mtw_t1_incl", "mT(W) type1 caloMet (incl, leading iso ele)",   100, 0, 200);
   addHisto_Njet_DataAndMC( h_mtw_mu, "mtw_mu", "mT(W) #mu-cor caloMet (after all but MET cut)",100, 0, 200);
   addHisto_Njet_DataAndMC( h_mtw_t1, "mtw_t1", "mT(W) type1 caloMet (after all but MET cut)", 100, 0, 200);
}
//------------------------------------------------------------------------------------


void ana::BookHistograms_DRemu(){

   //----------------------
   //  Delta R (ele, muon)
   //----------------------

   if(m_debug) cout << "Starting BookHistograms_DRemu()"<< endl;
   TDirectory *dir_DRemu = histf->mkdir("DRemu","Delta R(e,mu)");
   dir_DRemu->cd();
   addHistoDataAndMC( h_DRemu_selE_GoodMu,      "DRemu_selE_GoodMu",      "#DeltaR(e,#mu) (selE,GoodMu)", 60,0,6);
   addHistoDataAndMC( h_DRemu_selE_GoodMu_pass, "DRemu_selE_GoodMu_pass", "#DeltaR(e,#mu) (selE,GoodMu) (passAllCut)", 60,0,6);
}
//------------------------------------------------------------------------------------


void ana::BookHistograms_DPhiEmet(){

   //----------------------------
   //  Delta phi (iso ele, met)
   //-----------------------------
   if(m_debug) cout << "Starting BookHistograms_DPhiEmet()"<< endl;
   TDirectory *dir_DPhiEmet = histf->mkdir("DPhi_ele_met","Delta phi(isoele,met)");
   dir_DPhiEmet->cd();
   addHistoDataAndMC( h_DPhiEmet_mu_incl, "DPhiEmet_mu_incl", "#Delta#Phi(e,met) #mu-cor caloMet (incl, leading isolated ele)",   64, -3.2, 3.2);
   addHistoDataAndMC( h_DPhiEmet_t1_incl, "DPhiEmet_t1_incl", "#Delta#Phi(e,met) type 1 caloMet (incl, leading isolated ele)",   64, -3.2, 3.2);
   addHisto_Njet_DataAndMC( h_DPhiEmet_mu, "DPhiEmet_mu", "#Delta#Phi(e,met) #mu-cor caloMet (after all but MET cut)", 64, -3.2, 3.2);
   addHisto_Njet_DataAndMC( h_DPhiEmet_t1, "DPhiEmet_t1", "#Delta#Phi(e,met) type 1 caloMet (after all but MET cut)", 64, -3.2, 3.2);
}
//------------------------------------------------------------------------------------


void ana::BookHistograms_DPhimetjet(){

   //---------------------------------------
   // 2-3-10: Delta phi (met, nearest jet) 
   //---------------------------------------
   if(m_debug) cout << "Starting BookHistograms_DPhimj()"<< endl;
   TDirectory *dir_DPhimj = histf->mkdir("DPhi_met_jet","Delta phi(met,nearest jet)");
   dir_DPhimj->cd();
   addHisto_Njet_DataAndMC( h_DPhiMetJet_mu_goodE,  "DPhiMetJet_mu_goodE",  "#Delta#phi(mumet,nearest jet) pass HLT,goodEle",   64,0,3.2);
   addHisto_Njet_DataAndMC( h_DPhiMetJet_t1_goodE,  "DPhiMetJet_t1_goodE",  "#Delta#phi(t1met,nearest jet) pass HLT,goodEle",   64,0,3.2);
   addHisto_Njet_DataAndMC( h_DPhiMetJet_gen_goodE, "DPhiMetJet_gen_goodE", "#Delta#phi(genmet,nearest jet) pass HLT,goodEle",  64,0,3.2);
   addHisto_Njet_DataAndMC( h_DPhiMetJet_mu,        "DPhiMetJet_mu",     "#Delta#phi(mumet,nearest jet) pass all but MET cut",  64,0,3.2);
   addHisto_Njet_DataAndMC( h_DPhiMetJet_t1,        "DPhiMetJet_t1",     "#Delta#phi(t1met,nearest jet) pass all but MET cut",  64,0,3.2);
   addHisto_Njet_DataAndMC( h_DPhiMetJet_gen,       "DPhiMetJet_gen",    "#Delta#phi(genmet,nearest jet) pass all but MET cut", 64,0,3.2);
}
//------------------------------------------------------------------------------------


void ana::BookHistograms_conv(){

   //--------------
   //  Conversion
   //--------------
   if(m_debug) cout << "Starting BookHistograms_conv()"<< endl;
   TDirectory *dir_conv = histf->mkdir("conversion","conversion");
   dir_conv->cd();
   //81FB
   Conv_Opti.reserve(2);
   Conv_OptiL.reserve(2);
   Conv_Optis.reserve(2);
   Conv_Opti_extragran.reserve(2);
   Conv_Opti[0]  = new TH2D("ConvOptimization_ttbar",  "ConvOptimization_ttbar",  100,  -1,  1, 100,  -1,  1);
   Conv_OptiL[0] = new TH2D("ConvOptimization_ttbar_l","ConvOptimization_ttbar_l",100, -10, 10, 100, -10, 10);
   Conv_Optis[0] = new TH2D("ConvOptimization_ttbar_s","ConvOptimization_ttbar_s",100,-0.1,0.1, 100,-0.1,0.1);
   Conv_Opti[1]  = new TH2D("ConvOptimization_enri3",  "ConvOptimization_enri3",  100,  -1,  1, 100,  -1,  1);
   Conv_OptiL[1] = new TH2D("ConvOptimization_enri3_l","ConvOptimization_enri3_l",100, -10, 10, 100, -10, 10);
   Conv_Optis[1] = new TH2D("ConvOptimization_enri3_s","ConvOptimization_enri3_s",100,-0.1,0.1, 100,-0.1,0.1);
   Conv_Opti_extragran[0] = new TH2D("ConvOptimization_ttbar_eg","ConvOptimization_ttbar_eg",2000,-0.2,0.2,2000,-0.2,0.2);
   Conv_Opti_extragran[1] = new TH2D("ConvOptimization_enri3_eg","ConvOptimization_enri3_eg",2000,-0.2,0.2,2000,-0.2,0.2);

   //Conv_CheckDelR_GSFTk_ctfTk = new TH1D("Conv_CheckDelR_GSFTk_ctfTk","Conv_CheckDelR_GSFTk_ctfTk",1000,0,10);
}
//------------------------------------------------------------------------------------



void ana::BookHistograms_QCD(){

   //====================
   //
   //  QCD estimation
   //
   //====================
   if(m_debug) cout << "Starting BookHistograms_QCD()"<< endl;

   TDirectory *dir_QCD = histf->mkdir("QCD_estimation");
   dir_QCD->cd();

   //----------------------
   // TL: 16-8-09
   // Simplify:  replace 100 lines with 4 lines
   //----------------------
   addHisto_Njet_DataAndMC( h_QCDest_CombRelIso,     "QCDest_CombRelIso",     "RelIso",     1000,0,10);

   //--------------------------------
   // AES (Anti Event Selection)
   //--------------------------------   
   TDirectory *dir_AES = dir_QCD->mkdir("AES","anti event selection");
   dir_AES->cd();
   addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES,               "QCDest_CombRelIso_AES",               "RelIso (AES)",               1000,0,10);
   addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_minusMET,      "QCDest_CombRelIso_AES_minusMET",      "RelIso (AES-met)",           1000,0,10);
   addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_minusHT,       "QCDest_CombRelIso_AES_minusHT",       "RelIso (AES-HT)",            1000,0,10);
   addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_minusTighterZ, "QCDest_CombRelIso_AES_minusTighterZ", "RelIso (AES-TighterZ)",      1000,0,10);
   addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_before,        "QCDest_CombRelIso_AES_before",        "RelIso (Before AES)",        1000,0,10);
   addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_justMET,       "QCDest_CombRelIso_AES_justMET",       "RelIso (AES: just MET)",     1000,0,10);
   addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_justHighMET,   "QCDest_CombRelIso_AES_justHighMET",   "RelIso (AES: just highMET)", 1000,0,10);
   addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_justHT,        "QCDest_CombRelIso_AES_justHT",        "RelIso (AES: just HT)",      1000,0,10);
   addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_justZ,         "QCDest_CombRelIso_AES_justZ",         "RelIso (AES: just tightZ)",  1000,0,10);
		
   //-------------------------------
   // NES (Normal Event Selection)
   //-------------------------------
   // 11 levels: 
   //-------------------
   //   1: L1    after HLT(all electrons); nGoodEle>0; muZveto; conv
   //   2: L1b   ele et,eta
   //   3: L1c   d0
   //   4: L1d1  eID H/E
   //   5: L1d2  eID Delta sigma in
   //   6: L1d3  eID Delta phi in
   //   7: L1d4  eID sigma ieta ieta
   //   8: L1d5  eID all (out-of-box)
   //   9: L2    GoodEle
   //  10: L3    !mu, !z
   //  11: L4    !conv
   //-------------------
   // notes: [11][7]=>[nLevel][njet]
   // scatter plot, iso:met


   if( m_plotRelisoNES ) {

     cout << "plotRelIsoNES" << endl;

     TDirectory *dir_NES = dir_QCD->mkdir("NES_new","normal event selection");
     dir_NES->cd();
     //   hh_QCDest_isoVmet_NES.reserve(11);//NEW 11=nLevel
     const int sc_nbIso=70; //iso in scatter plot
     const int sc_maxIso=70; //iso in scatter plot
     const int sc_nbMet=50; //met in scatter plot
     const int sc_maxMet=50; //met in scatter plot
     // weighted
     iso_addHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES,        "QCDest_isoVmet_NES",        "RelIso v MET (NES)",  sc_nbIso,0,sc_maxIso, sc_nbMet,0,sc_maxMet);
     iso_addHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES_barrel, "QCDest_isoVmet_NES_barrel", "RelIso v MET (NES barrel)",  sc_nbIso,0,sc_maxIso, sc_nbMet,0,sc_maxMet);
     iso_addHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES_endcap, "QCDest_isoVmet_NES_endcap", "RelIso v MET (NES endcap)",  sc_nbIso,0,sc_maxIso, sc_nbMet,0,sc_maxMet);
     // unweighted
     iso_addHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES_uw,        "QCDest_isoVmet_NES_uw",        "RelIso v MET (NES uw)",   sc_nbIso,0,sc_maxIso, sc_nbMet,0,sc_maxMet);
     iso_addHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES_uw_barrel, "QCDest_isoVmet_NES_uw_barrel", "RelIso v MET (NES uw barrel)",   sc_nbIso,0,sc_maxIso, sc_nbMet,0,sc_maxMet);
     iso_addHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES_uw_endcap, "QCDest_isoVmet_NES_uw_endcap", "RelIso v MET (NES uw endcap)",  sc_nbIso,0,sc_maxIso, sc_nbMet,0,sc_maxMet);

     // iso plots
     const int nbIso = 140; //nbin reliso
     const float upIso = 1.4; //max
     iso_addHisto_nlevel_nj_nmc(h_QCDest_iso_NES,        "QCDest_iso_NES",        "RelIso (NES)",        nbIso, 0,upIso );
     iso_addHisto_nlevel_nj_nmc(h_QCDest_iso_NES_barrel, "QCDest_iso_NES_barrel", "RelIso (NES barrel)", nbIso, 0,upIso );
     iso_addHisto_nlevel_nj_nmc(h_QCDest_iso_NES_endcap, "QCDest_iso_NES_endcap", "RelIso (NES endcap)", nbIso, 0,upIso );

     iso_addHisto_nlevel_nj_nmc(h_QCDest_iso_NES_hiMET,        "QCDest_iso_NES_hiMET",        "RelIso (NES hiMET)",        nbIso, 0,upIso );
     iso_addHisto_nlevel_nj_nmc(h_QCDest_iso_NES_hiMET_barrel, "QCDest_iso_NES_hiMET_barrel", "RelIso (NES hiMET barrel)", nbIso, 0,upIso );
     iso_addHisto_nlevel_nj_nmc(h_QCDest_iso_NES_hiMET_endcap, "QCDest_iso_NES_hiMET_endcap", "RelIso (NES hiMET endcap)", nbIso, 0,upIso );

     iso_addHisto_nlevel_nj_nmc(h_QCDest_iso_NES_loMET,        "QCDest_iso_NES_loMET",        "RelIso (NES loMET)",        nbIso, 0,upIso );
     iso_addHisto_nlevel_nj_nmc(h_QCDest_iso_NES_loMET_barrel, "QCDest_iso_NES_loMET_barrel", "RelIso (NES loMET barrel)", nbIso, 0,upIso );
     iso_addHisto_nlevel_nj_nmc(h_QCDest_iso_NES_loMET_endcap, "QCDest_iso_NES_loMET_endcap", "RelIso (NES loMET endcap)", nbIso, 0,upIso );


     /*
     dir_NES->cd();
     const string levelno[nLevel] = {"L1","L1b","L1c","L1d1","L1d2","L1d3","L1d4","L1d5","L2","L3","L4"};
     const string levelinfo[nLevel] = {"L1 HLT","L1b E_{T},#eta","L1c d0",
				       "L1d1 eID H/E",
				       "L1d2 eID |#Delta#sigma_{in}|",
				       "L1d3 eID |#Delta#phi_{in}|",
				       "L1d4 eID #sigma_{i#eta i#eta}",
				       "L1d5 eID",
				       "L2 GoodEle","L3 MuZVeto","L4 Conv"};

     for (int iLevel=0; iLevel<nLevel; ++iLevel) {

       const string aname_sc   = "QCDest_isoVmet_NES_"                 + levelno[iLevel];
       const string aname_scBA = "QCDest_isoVmet_NES_barrel_"          + levelno[iLevel];
       const string aname_scEN = "QCDest_isoVmet_NES_endcap_"          + levelno[iLevel];
       const string aname_sc2  = "QCDest_isoVmet_NES_uw_"              + levelno[iLevel];//unweighted
       const string aname_sc2BA = "QCDest_isoVmet_NES_uw_barrel_"      + levelno[iLevel];//unweighted
       const string aname_sc2EN = "QCDest_isoVmet_NES_uw_endcap_"      + levelno[iLevel];//unweighted
       const string aname      = "QCDest_CombRelIso_NES_"              + levelno[iLevel];
       const string anameBA    = "QCDest_CombRelIso_NES_barrel_"       + levelno[iLevel];
       const string anameEN    = "QCDest_CombRelIso_NES_endcap_"       + levelno[iLevel];
       const string aname_lo   = "QCDest_CombRelIso_NES_loMET_"        + levelno[iLevel];
       const string aname_loBA = "QCDest_CombRelIso_NES_loMET_barrel_" + levelno[iLevel];
       const string aname_loEN = "QCDest_CombRelIso_NES_loMET_endcap_" + levelno[iLevel];
       const string aname_hi   = "QCDest_CombRelIso_NES_hiMET_"        + levelno[iLevel];
       const string aname_hiBA = "QCDest_CombRelIso_NES_hiMET_barrel_" + levelno[iLevel];
       const string aname_hiEN = "QCDest_CombRelIso_NES_hiMET_endcap_" + levelno[iLevel];

       const string info_sc   = Form("RelIso v MET (NES_%s)",        levelinfo[iLevel].c_str() );
       const string info_scBA = Form("RelIso v MET (NES_%s barrel)", levelinfo[iLevel].c_str() );
       const string info_scEN = Form("RelIso v MET (NES_%s endcap)", levelinfo[iLevel].c_str() );
       const string info_sc2  = Form("RelIso v MET (NES_%s uw)",     levelinfo[iLevel].c_str() );//unweighted
       const string info_sc2BA = Form("RelIso v MET (NES_%s uw barrel)",  levelinfo[iLevel].c_str() );//unweighted
       const string info_sc2EN = Form("RelIso v MET (NES_%s uw endcap)",  levelinfo[iLevel].c_str() );//unweighted
       const string info      = Form("RelIso (NES_%s)",              levelinfo[iLevel].c_str() );
       const string infoBA    = Form("RelIso (NES_%s barrel)",       levelinfo[iLevel].c_str() );
       const string infoEN    = Form("RelIso (NES_%s endcap)",       levelinfo[iLevel].c_str() );
       const string info_lo   = Form("RelIso (NES_%s loMET)",        levelinfo[iLevel].c_str() );
       const string info_loBA = Form("RelIso (NES_%s loMET barrel)", levelinfo[iLevel].c_str() );
       const string info_loEN = Form("RelIso (NES_%s loMET endcap)", levelinfo[iLevel].c_str() );
       const string info_hi   = Form("RelIso (NES_%s hiMET)",        levelinfo[iLevel].c_str() );
       const string info_hiBA = Form("RelIso (NES_%s hiMET barrel)", levelinfo[iLevel].c_str() );
       const string info_hiEN = Form("RelIso (NES_%s hiMET endcap)", levelinfo[iLevel].c_str() );
   
       const int nb=100; //scatter plot
       const int nBin = 200; //reliso
       const float xUp = 2.0; //max
       cout << "iLevel : "<< iLevel  << endl;
     
       cout << "\n   aa"<< endl;
       addHistoNjet_TH2(h_QCDest_isoVmet_NES[iLevel],                 aname_sc,   "", info_sc,   nb,0,2,nb,0,100);
       cout << "\n   bb"<< endl;

       addHistoNjet_TH2(h_QCDest_isoVmet_NES_barrel[iLevel],          aname_scBA, "", info_scBA, nb,0,2,nb,0,100);
       cout << "\n   cc"<< endl;
       addHistoNjet_TH2(h_QCDest_isoVmet_NES_endcap[iLevel],          aname_scEN, "", info_scEN, nb,0,2,nb,0,100);
       cout << "\n   dd"<< endl;
       addHistoNjet(h_QCDest_CombRelIso_NES[iLevel],              aname,      "", info,      nBin,0,xUp);
       cout << "'n   ee"<< endl;
       addHistoNjet(h_QCDest_CombRelIso_NES_barrel[iLevel],       anameBA,    "", infoBA,    nBin,0,xUp);
       cout << "\n   ff"<< endl;
       addHistoNjet(h_QCDest_CombRelIso_NES_endcap[iLevel],       anameEN,    "", infoEN,    nBin,0,xUp);
       addHistoNjet(h_QCDest_CombRelIso_NES_loMET[iLevel],        aname_lo,   "", info_lo,   nBin,0,xUp);
       addHistoNjet(h_QCDest_CombRelIso_NES_loMET_barrel[iLevel], aname_loBA, "", info_loBA, nBin,0,xUp);
       addHistoNjet(h_QCDest_CombRelIso_NES_loMET_endcap[iLevel], aname_loEN, "", info_loEN, nBin,0,xUp);
       addHistoNjet(h_QCDest_CombRelIso_NES_hiMET[iLevel],        aname_hi,   "", info_hi,   nBin,0,xUp);
       addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_barrel[iLevel], aname_hiBA, "", info_hiBA, nBin,0,xUp);
       addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_endcap[iLevel], aname_hiEN, "", info_hiEN, nBin,0,xUp);

       if(!IsData()){
	 // scatter plot, iso:met (weighted)
	 cout << "\n   QCD"<< endl;
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_QCD[iLevel],   aname_sc, "__QCD",   info_sc+" (QCD)", nb,0,2,nb,0,100);
	 cout << "\n   bce1"<< endl;
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_bce1[iLevel],  aname_sc, "__bce1",  info_sc+" (bce1)",nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_bce2[iLevel],  aname_sc, "__bce2",  info_sc+" (bce2)",nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_bce3[iLevel],  aname_sc, "__bce3",  info_sc+" (bce3)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_enri1[iLevel], aname_sc, "__enri1", info_sc+" (enri1)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_enri2[iLevel], aname_sc, "__enri2", info_sc+" (enri2)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_enri3[iLevel], aname_sc, "__enri3", info_sc+" (enri3)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_ttbar[iLevel], aname_sc, "__ttbar", info_sc+" (ttbar)",  nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_Wjet[iLevel],  aname_sc, "__Wjet",  info_sc+" (W+jets)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_Zjet[iLevel],  aname_sc, "__Zjet",  info_sc+" (Z+jets)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_singleTop[iLevel],aname_sc, "__singleTop", info_sc+" (single top)", nb,0,2,nb,0,100);

	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_barrel_QCD[iLevel],   aname_scBA, "__QCD",   info_scBA+" (QCD)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_barrel_bce1[iLevel],  aname_scBA, "__bce1",  info_scBA+" (bce1)",nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_barrel_bce2[iLevel],  aname_scBA, "__bce2",  info_scBA+" (bce2)",nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_barrel_bce3[iLevel],  aname_scBA, "__bce3",  info_scBA+" (bce3)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_barrel_enri1[iLevel], aname_scBA, "__enri1", info_scBA+" (enri1)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_barrel_enri2[iLevel], aname_scBA, "__enri2", info_scBA+" (enri2)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_barrel_enri3[iLevel], aname_scBA, "__enri3", info_scBA+" (enri3)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_barrel_ttbar[iLevel], aname_scBA, "__ttbar", info_scBA+" (ttbar)",  nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_barrel_Wjet[iLevel],  aname_scBA, "__Wjet",  info_scBA+" (W+jets)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_barrel_Zjet[iLevel],  aname_scBA, "__Zjet",  info_scBA+" (Z+jets)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_barrel_singleTop[iLevel],aname_scBA, "__singleTop", info_scBA+" (single top)", nb,0,2,nb,0,100);

	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_endcap_QCD[iLevel],      aname_scEN, "__QCD",   info_scEN+" (QCD)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_endcap_bce1[iLevel],     aname_scEN, "__bce1",  info_scEN+" (bce1)",nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_endcap_bce2[iLevel],     aname_scEN, "__bce2",  info_scEN+" (bce2)",nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_endcap_bce3[iLevel],     aname_scEN, "__bce3",  info_scEN+" (bce3)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_endcap_enri1[iLevel],    aname_scEN, "__enri1", info_scEN+" (enri1)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_endcap_enri2[iLevel],    aname_scEN, "__enri2", info_scEN+" (enri2)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_endcap_enri3[iLevel],    aname_scEN, "__enri3", info_scEN+" (enri3)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_endcap_ttbar[iLevel],    aname_scEN, "__ttbar", info_scEN+" (ttbar)",  nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_endcap_Wjet[iLevel],     aname_scEN, "__Wjet",  info_scEN+" (W+jets)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_endcap_Zjet[iLevel],     aname_scEN, "__Zjet",  info_scEN+" (Z+jets)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_endcap_singleTop[iLevel],aname_scEN, "__singleTop", info_scEN+" (single top)", nb,0,2,nb,0,100);

	 // scatter plot, iso:met (unweighted, for each MC)
	 //addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_QCD[iLevel],      aname_sc2, "__QCD",   info_sc2+" (QCD)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_bce1[iLevel],     aname_sc2, "__bce1",  info_sc2+" (bce1)",nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_bce2[iLevel],     aname_sc2, "__bce2",  info_sc2+" (bce2)",nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_bce3[iLevel],     aname_sc2, "__bce3",  info_sc2+" (bce3)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_enri1[iLevel],    aname_sc2, "__enri1", info_sc2+" (enri1)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_enri2[iLevel],    aname_sc2, "__enri2", info_sc2+" (enri2)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_enri3[iLevel],    aname_sc2, "__enri3", info_sc2+" (enri3)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_ttbar[iLevel],    aname_sc2, "__ttbar", info_sc2+" (ttbar)",  nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_Wjet[iLevel],     aname_sc2, "__Wjet",  info_sc2+" (W+jets)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_Zjet[iLevel],     aname_sc2, "__Zjet",  info_sc2+" (Z+jets)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_singleTop[iLevel],aname_sc2, "__singleTop", info_sc2+" (single top)", nb,0,2,nb,0,100);

	 //	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_barrel_QCD[iLevel],   aname_sc2BA, "__QCD",   info_sc2BA+" (QCD)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_barrel_bce1[iLevel],  aname_sc2BA, "__bce1",  info_sc2BA+" (bce1)",nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_barrel_bce2[iLevel],  aname_sc2BA, "__bce2",  info_sc2BA+" (bce2)",nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_barrel_bce3[iLevel],  aname_sc2BA, "__bce3",  info_sc2BA+" (bce3)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_barrel_enri1[iLevel], aname_sc2BA, "__enri1", info_sc2BA+" (enri1)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_barrel_enri2[iLevel], aname_sc2BA, "__enri2", info_sc2BA+" (enri2)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_barrel_enri3[iLevel], aname_sc2BA, "__enri3", info_sc2BA+" (enri3)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_barrel_ttbar[iLevel], aname_sc2BA, "__ttbar", info_sc2BA+" (ttbar)",  nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_barrel_Wjet[iLevel],  aname_sc2BA, "__Wjet",  info_sc2BA+" (W+jets)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_barrel_Zjet[iLevel],  aname_sc2BA, "__Zjet",  info_sc2BA+" (Z+jets)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_barrel_singleTop[iLevel],aname_sc2BA, "__singleTop", info_sc2BA+" (single top)", nb,0,2,nb,0,100);

	 //	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_endcap_QCD[iLevel],   aname_sc2EN, "__QCD",   info_sc2EN+" (QCD)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_endcap_bce1[iLevel],  aname_sc2EN, "__bce1",  info_sc2EN+" (bce1)",nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_endcap_bce2[iLevel],  aname_sc2EN, "__bce2",  info_sc2EN+" (bce2)",nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_endcap_bce3[iLevel],  aname_sc2EN, "__bce3",  info_sc2EN+" (bce3)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_endcap_enri1[iLevel], aname_sc2EN, "__enri1", info_sc2EN+" (enri1)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_endcap_enri2[iLevel], aname_sc2EN, "__enri2", info_sc2EN+" (enri2)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_endcap_enri3[iLevel], aname_sc2EN, "__enri3", info_sc2EN+" (enri3)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_endcap_ttbar[iLevel], aname_sc2EN, "__ttbar", info_sc2EN+" (ttbar)",  nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_endcap_Wjet[iLevel],  aname_sc2EN, "__Wjet",  info_sc2EN+" (W+jets)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_endcap_Zjet[iLevel],  aname_sc2EN, "__Zjet",  info_sc2EN+" (Z+jets)", nb,0,2,nb,0,100);
	 addHistoNjet_TH2(h_QCDest_isoVmet_NES_uw_endcap_singleTop[iLevel],aname_sc2EN, "__singleTop", info_sc2EN+" (single top)", nb,0,2,nb,0,100);


	 // a) no met cut
	 addHistoNjet(h_QCDest_CombRelIso_NES_QCD[iLevel],   aname, "__QCD",   info+" (QCD)", nBin,0,xUp);       
	 addHistoNjet(h_QCDest_CombRelIso_NES_bce1[iLevel],  aname, "__bce1",  info+" (bce1)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_bce2[iLevel],  aname, "__bce2",  info+" (bce2)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_bce3[iLevel],  aname, "__bce3",  info+" (bce3)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_enri1[iLevel], aname, "__enri1", info+" (enri1)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_enri2[iLevel], aname, "__enri2", info+" (enri2)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_enri3[iLevel], aname, "__enri3", info+" (enri3)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_ttbar[iLevel], aname, "__ttbar", info+" (ttbar)",  nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_Wjet[iLevel],  aname, "__Wjet",  info+" (W+jets)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_Zjet[iLevel],  aname, "__Zjet",  info+" (Z+jets)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_singleTop[iLevel],aname, "__singleTop", info+" (single top)", nBin,0,xUp);            
	 // no met cut, barrel (QCD)
	 addHistoNjet(h_QCDest_CombRelIso_NES_barrel_QCD[iLevel],   anameBA, "__QCD",   infoBA+" (QCD)", nBin,0,xUp);       
	 addHistoNjet(h_QCDest_CombRelIso_NES_barrel_bce1[iLevel],  anameBA, "__bce1",  infoBA+" (bce1)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_barrel_bce2[iLevel],  anameBA, "__bce2",  infoBA+" (bce2)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_barrel_bce3[iLevel],  anameBA, "__bce3",  infoBA+" (bce3)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_barrel_enri1[iLevel], anameBA, "__enri1", infoBA+" (enri1)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_barrel_enri2[iLevel], anameBA, "__enri2", infoBA+" (enri2)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_barrel_enri3[iLevel], anameBA, "__enri3", infoBA+" (enri3)", nBin,0,xUp);
	 // no met cut, endcap (QCD)
	 addHistoNjet(h_QCDest_CombRelIso_NES_endcap_QCD[iLevel],   anameEN, "__QCD",   infoEN+" (QCD)", nBin,0,xUp);       
	 addHistoNjet(h_QCDest_CombRelIso_NES_endcap_bce1[iLevel],  anameEN, "__bce1",  infoEN+" (bce1)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_endcap_bce2[iLevel],  anameEN, "__bce2",  infoEN+" (bce2)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_endcap_bce3[iLevel],  anameEN, "__bce3",  infoEN+" (bce3)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_endcap_enri1[iLevel], anameEN, "__enri1", infoEN+" (enri1)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_endcap_enri2[iLevel], anameEN, "__enri2", infoEN+" (enri2)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_endcap_enri3[iLevel], anameEN, "__enri3", infoEN+" (enri3)", nBin,0,xUp);

	 // b) lowMET
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_QCD[iLevel],   aname_lo, "__QCD",   info_lo+" (QCD)", nBin,0,xUp);       
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_bce1[iLevel],  aname_lo, "__bce1",  info_lo+" (bce1)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_bce2[iLevel],  aname_lo, "__bce2",  info_lo+" (bce2)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_bce3[iLevel],  aname_lo, "__bce3",  info_lo+" (bce3)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_enri1[iLevel], aname_lo, "__enri1", info_lo+" (enri1)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_enri2[iLevel], aname_lo, "__enri2", info_lo+" (enri2)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_enri3[iLevel], aname_lo, "__enri3", info_lo+" (enri3)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_ttbar[iLevel], aname_lo, "__ttbar", info_lo+" (ttbar)",  nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_Wjet[iLevel],  aname_lo, "__Wjet",  info_lo+" (W+jets)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_Zjet[iLevel],  aname_lo, "__Zjet",  info_lo+" (Z+jets)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_singleTop[iLevel],aname_lo, "__singleTop", info_lo+" (single top)", nBin,0,xUp);     
	 // loMET, barrel (QCD)
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_barrel_QCD[iLevel],   aname_loBA, "__QCD",   info_loBA+" (QCD)", nBin,0,xUp);       
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_barrel_bce1[iLevel],  aname_loBA, "__bce1",  info_loBA+" (bce1)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_barrel_bce2[iLevel],  aname_loBA, "__bce2",  info_loBA+" (bce2)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_barrel_bce3[iLevel],  aname_loBA, "__bce3",  info_loBA+" (bce3)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_barrel_enri1[iLevel], aname_loBA, "__enri1", info_loBA+" (enri1)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_barrel_enri2[iLevel], aname_loBA, "__enri2", info_loBA+" (enri2)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_barrel_enri3[iLevel], aname_loBA, "__enri3", info_loBA+" (enri3)", nBin,0,xUp);
	 // loMET, endcap (QCD)
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_endcap_QCD[iLevel],   aname_loEN, "__QCD",   info_loEN+" (QCD)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_endcap_bce1[iLevel],  aname_loEN, "__bce1",  info_loEN+" (bce1)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_endcap_bce2[iLevel],  aname_loEN, "__bce2",  info_loEN+" (bce2)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_endcap_bce3[iLevel],  aname_loEN, "__bce3",  info_loEN+" (bce3)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_endcap_enri1[iLevel], aname_loEN, "__enri1", info_loEN+" (enri1)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_endcap_enri2[iLevel], aname_loEN, "__enri2", info_loEN+" (enri2)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_loMET_endcap_enri3[iLevel], aname_loEN, "__enri3", info_loEN+" (enri3)", nBin,0,xUp);

	 // c) highMET
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_QCD[iLevel],   aname_hi, "__QCD",   info_hi+" (QCD)", nBin,0,xUp);       
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_bce1[iLevel],  aname_hi, "__bce1",  info_hi+" (bce1)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_bce2[iLevel],  aname_hi, "__bce2",  info_hi+" (bce2)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_bce3[iLevel],  aname_hi, "__bce3",  info_hi+" (bce3)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_enri1[iLevel], aname_hi, "__enri1", info_hi+" (enri1)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_enri2[iLevel], aname_hi, "__enri2", info_hi+" (enri2)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_enri3[iLevel], aname_hi, "__enri3", info_hi+" (enri3)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_ttbar[iLevel], aname_hi, "__ttbar", info_hi+" (ttbar)",  nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_Wjet[iLevel],  aname_hi, "__Wjet",  info_hi+" (W+jets)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_Zjet[iLevel],  aname_hi, "__Zjet",  info_hi+" (Z+jets)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_singleTop[iLevel],aname_hi, "__singleTop", info_hi+" (single top)", nBin,0,xUp);            
	 // hiMET, barrel (QCD)
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_barrel_QCD[iLevel],   aname_hiBA, "__QCD",   info_hiBA+" (QCD)", nBin,0,xUp);       
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_barrel_bce1[iLevel],  aname_hiBA, "__bce1",  info_hiBA+" (bce1)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_barrel_bce2[iLevel],  aname_hiBA, "__bce2",  info_hiBA+" (bce2)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_barrel_bce3[iLevel],  aname_hiBA, "__bce3",  info_hiBA+" (bce3)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_barrel_enri1[iLevel], aname_hiBA, "__enri1", info_hiBA+" (enri1)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_barrel_enri2[iLevel], aname_hiBA, "__enri2", info_hiBA+" (enri2)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_barrel_enri3[iLevel], aname_hiBA, "__enri3", info_hiBA+" (enri3)", nBin,0,xUp);
	 // hiMET, endcap (QCD)
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_endcap_QCD[iLevel],   aname_hiEN, "__QCD",   info_hiEN+" (QCD)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_endcap_bce1[iLevel],  aname_hiEN, "__bce1",  info_hiEN+" (bce1)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_endcap_bce2[iLevel],  aname_hiEN, "__bce2",  info_hiEN+" (bce2)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_endcap_bce3[iLevel],  aname_hiEN, "__bce3",  info_hiEN+" (bce3)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_endcap_enri1[iLevel], aname_hiEN, "__enri1", info_hiEN+" (enri1)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_endcap_enri2[iLevel], aname_hiEN, "__enri2", info_hiEN+" (enri2)", nBin,0,xUp);
	 addHistoNjet(h_QCDest_CombRelIso_NES_hiMET_endcap_enri3[iLevel], aname_hiEN, "__enri3", info_hiEN+" (enri3)", nBin,0,xUp);

       }//MC
     }//4 levels of cut     
     */
   }//if m_plotRelisoNES

}//end
//------------------------------------------------------------------------------------


void ana::BookHistograms_QCD_planA(){

  //-----------------------------------------
  //  Trial definitions of AES for Plan A
  //-----------------------------------------
  if(m_debug) cout << "Starting BookHistogrmas_QCD_planA" << endl;
  TDirectory *dir_QCD_planA = histf->mkdir("QCD_planA");
  dir_QCD_planA->cd();
  addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planA1_e20, "QCDest_CombRelIso_AES_planA1_e20", "RelIso (AES A1, no conv cut, E_{T}>20)", 1000,0,10);
  addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planA1_e30, "QCDest_CombRelIso_AES_planA1_e30", "RelIso (AES A1, no conv cut, E_{T}>30)", 1000,0,10);
  addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planA2_e20, "QCDest_CombRelIso_AES_planA2_e20", "RelIso (AES A2, conv, E_{T}>20)", 1000,0,10);
  addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planA2_e30, "QCDest_CombRelIso_AES_planA2_e30", "RelIso (AES A2, conv, E_{T}>30)", 1000,0,10);
  addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planA3_e20, "QCDest_CombRelIso_AES_planA3_e20", "RelIso (AES A3 |d_{0}|>200um, E_{T}>20)", 1000,0,10);
  addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planA3_e30, "QCDest_CombRelIso_AES_planA3_e30", "RelIso (AES A3 |d_{0}|>200um, E_{T}>30)", 1000,0,10);
}//end
//------------------------------------------------------------------------------------


void ana::BookHistograms_QCD_planB(){

  //-----------------------------------------
  //  Trial definitions of AES for Plan B
  //-----------------------------------------
  if(m_debug) cout << "Starting BookHistogrmas_QCD_planB" << endl;
  TDirectory *dir_QCD_planB = histf->mkdir("QCD_planB");
  dir_QCD_planB->cd();
  addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB1_e20, "QCDest_CombRelIso_AES_planB1_e20", "RelIso (AES B1 E_{T}>20)", 1000,0,10);
  addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB1_e30, "QCDest_CombRelIso_AES_planB1_e30", "RelIso (AES B1 E_{T}>30)", 1000,0,10);
  addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB2_e20, "QCDest_CombRelIso_AES_planB2_e20", "RelIso (AES B2 E_{T}<20)", 1000,0,10);
  addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB2_e30, "QCDest_CombRelIso_AES_planB2_e30", "RelIso (AES B2 E_{T}<30)", 1000,0,10);
  addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB3_e20, "QCDest_CombRelIso_AES_planB3_e20", "RelIso (AES B3 RT=0 E_{T}>20)", 1000,0,10);
  addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB3_e30, "QCDest_CombRelIso_AES_planB3_e30", "RelIso (AES B3 RT=0 E_{T}>30)", 1000,0,10);
  addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB3b_e20, "QCDest_CombRelIso_AES_planB3b_e20", "RelIso (AES B3 RT=0 EB E_{T}>20)", 1000,0,10);//new
  addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB3b_e30, "QCDest_CombRelIso_AES_planB3b_e30", "RelIso (AES B3 RT=0 EB E_{T}>30)", 1000,0,10);//new
  addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB4_e20, "QCDest_CombRelIso_AES_planB4_e20", "RelIso (AES B4 RL=0 E_{T}>20)", 1000,0,10);
  addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB4_e30, "QCDest_CombRelIso_AES_planB4_e30", "RelIso (AES B4 RL=0 E_{T}>30)", 1000,0,10);
  addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB5_e20, "QCDest_CombRelIso_AES_planB5_e20", "RelIso (AES B5 cL=0 E_{T}>20)", 1000,0,10);
  addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB5_e30, "QCDest_CombRelIso_AES_planB5_e30", "RelIso (AES B5 cL=0 E_{T}>30)", 1000,0,10);
  addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB6_e20, "QCDest_CombRelIso_AES_planB6_e20", "RelIso (AES B6 cT=0 E_{T}>20)", 1000,0,10);
  addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB6_e30, "QCDest_CombRelIso_AES_planB6_e30", "RelIso (AES B6 cT=0 E_{T}>30)", 1000,0,10);
  addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB7_e20, "QCDest_CombRelIso_AES_planB7_e20", "RelIso (AES B7 |d_{0}|>200um RT=1 E_{T}>20)", 1000,0,10);
  addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB7_e30, "QCDest_CombRelIso_AES_planB7_e30", "RelIso (AES B7 |d_{0}|>200um RT=1 E_{T}>30)", 1000,0,10);
  addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB8_e20, "QCDest_CombRelIso_AES_planB8_e20", "RelIso (AES B8 |d_{0}|>200um RT=0 E_{T}>20)", 1000,0,10);
  addHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB8_e30, "QCDest_CombRelIso_AES_planB8_e30", "RelIso (AES B8 |d_{0}|>200um RT=0 E_{T}>30)", 1000,0,10);

}//end BookHistograms_QCD_planB()
//------------------------------------------------------------------------------------


void ana::BookHistograms_wj(){

   //--------------------------------------------------
   //  W+jets estimation: m3
   //--------------------------------------------------
   if(m_debug) cout << "Starting BookHistograms_wj()"<< endl;
   TDirectory *dir_wjets = histf->mkdir("Wjets_estimation","m3 plots");
   dir_wjets->cd();
   // m3 for signal region (isolated ele)
   h_hadTop_maxPT_mass_4j = new TH1D("hadTop_maxPT_mass_4j", "had top mass (m3) (>=4j)",  nbm3,0,nbm3);
   h_hadTop_maxPT_pt_4j   = new TH1D("hadTop_maxPT_pt_4j",   "had top highest PT (>=4j)", 100,0,1000);
   // m3 for background region (non-isolated ele)
   h_hadTop_maxPT_mass_nonIso_4j = new TH1D("hadTop_maxPT_mass_nonIso_4j", "had top mass (m3) (>=4j) nonIsoE", nbm3,0,nbm3);
   h_hadTop_maxPT_pt_nonIso_4j   = new TH1D("hadTop_maxPT_pt_nonIso_4j",   "had top highest PT (>=4j) nonIsoE", 100,0,1000);
   h_hadTop_maxPT_mass_4j->Sumw2();
   h_hadTop_maxPT_pt_4j  ->Sumw2();
   h_hadTop_maxPT_mass_nonIso_4j->Sumw2();
   h_hadTop_maxPT_pt_nonIso_4j  ->Sumw2();

   TDirectory *dir_wj2 = histf->mkdir("Wjets_estimation_nB1000","m3 plots (1000 bins, 0-1000)");
   dir_wj2->cd();
   h_hadTop_maxPT_mass_4j_1000        = new TH1D("hadTop_maxPT_mass_4j",        "had top mass (m3) (>=4j)",         1000,0,1000);
   h_hadTop_maxPT_mass_nonIso_4j_1000 = new TH1D("hadTop_maxPT_mass_nonIso_4j", "had top mass (m3) (>=4j) nonIsoE", 1000,0,1000);
   h_hadTop_maxPT_mass_4j_1000->Sumw2();
   h_hadTop_maxPT_mass_nonIso_4j_1000->Sumw2();

   // m3 for each type of MC (if running on MC)
   if (IsData()==false) {
     // (i) signal region
     //
     // TODO: To-be-reviewed: Number of bins used for M3
     //
     dir_wjets->cd();
     h_m3_tt  = new TH1D("m3_tt", "m3 (ttbar MC)", nbm3,0,nbm3);
     h_m3_wj  = new TH1D("m3_wj", "m3 (W+jets MC)",nbm3,0,nbm3);
     h_m3_zj  = new TH1D("m3_zj", "m3 (Z+jets MC)",nbm3,0,nbm3);
     h_m3_qcd = new TH1D("m3_qcd","m3 (QCD MC)",   nbm3,0,nbm3);
     h_m3_vqq = new TH1D("m3_vqq","m3 (VQQ MC)",   nbm3,0,nbm3);
     h_m3_singletop = new TH1D("m3_singletop","m3 (single top MC)",nbm3,0,nbm3);
     h_m3_bce[0]  = new TH1D("m3_bce1", "m3 (QCD bce1 MC)",   nbm3,0,nbm3);
     h_m3_bce[1]  = new TH1D("m3_bce2", "m3 (QCD bce2 MC)",   nbm3,0,nbm3);
     h_m3_bce[2]  = new TH1D("m3_bce3", "m3 (QCD bce3 MC)",   nbm3,0,nbm3);
     h_m3_enri[0] = new TH1D("m3_enri1","m3 (QCD enri1 MC)",  nbm3,0,nbm3);
     h_m3_enri[1] = new TH1D("m3_enri2","m3 (QCD enri2 MC)",  nbm3,0,nbm3);
     h_m3_enri[2] = new TH1D("m3_enri3","m3 (QCD enri3 MC)",  nbm3,0,nbm3);
     // (ii) control region
     h_m3_tt_control  = new TH1D("m3_tt_control", "m3 (ttbar MC, control reg) ", nbm3,0,nbm3);
     h_m3_wj_control  = new TH1D("m3_wj_control", "m3 (W+jets MC, control reg)", nbm3,0,nbm3);
     h_m3_zj_control  = new TH1D("m3_zj_control", "m3 (Z+jets MC, control reg)", nbm3,0,nbm3);
     h_m3_qcd_control = new TH1D("m3_qcd_control","m3 (QCD MC, control reg)",    nbm3,0,nbm3);
     h_m3_vqq_control = new TH1D("m3_vqq_control","m3 (VQQ MC, control reg)",    nbm3,0,nbm3);
     h_m3_singletop_control = new TH1D("m3_singletop_control","m3 (single top MC, control reg)", nbm3,0,nbm3);
     h_m3_bce_control[0]  = new TH1D("m3_bce1_control","m3 (QCD bce1 MC, control reg)",   nbm3,0,nbm3);
     h_m3_bce_control[1]  = new TH1D("m3_bce2_control","m3 (QCD bce2 MC, control reg)",   nbm3,0,nbm3);
     h_m3_bce_control[2]  = new TH1D("m3_bce3_control","m3 (QCD bce3 MC, control reg)",   nbm3,0,nbm3);
     h_m3_enri_control[0] = new TH1D("m3_enri1_control","m3 (QCD enri1 MC, control reg)", nbm3,0,nbm3);
     h_m3_enri_control[1] = new TH1D("m3_enri2_control","m3 (QCD enri2 MC, control reg)", nbm3,0,nbm3);
     h_m3_enri_control[2] = new TH1D("m3_enri3_control","m3 (QCD enri3 MC, control reg)", nbm3,0,nbm3);

     // 1000 bins
     dir_wj2->cd();
     // (i) signal region
     h_m3_tt_1000  = new TH1D("m3_tt", "m3 (ttbar MC)", 1000,0,1000);
     h_m3_wj_1000  = new TH1D("m3_wj", "m3 (W+jets MC)",1000,0,1000);
     h_m3_zj_1000  = new TH1D("m3_zj", "m3 (Z+jets MC)",1000,0,1000);
     h_m3_qcd_1000 = new TH1D("m3_qcd","m3 (QCD MC)",   1000,0,1000);
     h_m3_vqq_1000 = new TH1D("m3_vqq","m3 (VQQ MC)",   1000,0,1000);
     h_m3_singletop_1000 = new TH1D("m3_singletop","m3 (single top MC)",1000,0,1000);
     h_m3_bce_1000[0]  = new TH1D("m3_bce1", "m3 (QCD bce1 MC)",   1000,0,1000);
     h_m3_bce_1000[1]  = new TH1D("m3_bce2", "m3 (QCD bce2 MC)",   1000,0,1000);
     h_m3_bce_1000[2]  = new TH1D("m3_bce3", "m3 (QCD bce3 MC)",   1000,0,1000);
     h_m3_enri_1000[0] = new TH1D("m3_enri1","m3 (QCD enri1 MC)",  1000,0,1000);
     h_m3_enri_1000[1] = new TH1D("m3_enri2","m3 (QCD enri2 MC)",  1000,0,1000);
     h_m3_enri_1000[2] = new TH1D("m3_enri3","m3 (QCD enri3 MC)",  1000,0,1000);
     // (ii) control region
     h_m3_tt_control_1000  = new TH1D("m3_tt_control", "m3 (ttbar MC, control reg) ", 1000,0,1000);
     h_m3_wj_control_1000  = new TH1D("m3_wj_control", "m3 (W+jets MC, control reg)", 1000,0,1000);
     h_m3_zj_control_1000  = new TH1D("m3_zj_control", "m3 (Z+jets MC, control reg)", 1000,0,1000);
     h_m3_qcd_control_1000 = new TH1D("m3_qcd_control","m3 (QCD MC, control reg)",    1000,0,1000);
     h_m3_vqq_control_1000 = new TH1D("m3_vqq_control","m3 (VQQ MC, control reg)",    1000,0,1000);
     h_m3_singletop_control_1000 = new TH1D("m3_singletop_control","m3 (single top MC, control reg)", 1000,0,1000);
     h_m3_bce_control_1000[0]  = new TH1D("m3_bce1_control","m3 (QCD bce1 MC, control reg)",   1000,0,1000);
     h_m3_bce_control_1000[1]  = new TH1D("m3_bce2_control","m3 (QCD bce2 MC, control reg)",   1000,0,1000);
     h_m3_bce_control_1000[2]  = new TH1D("m3_bce3_control","m3 (QCD bce3 MC, control reg)",   1000,0,1000);
     h_m3_enri_control_1000[0] = new TH1D("m3_enri1_control","m3 (QCD enri1 MC, control reg)", 1000,0,1000);
     h_m3_enri_control_1000[1] = new TH1D("m3_enri2_control","m3 (QCD enri2 MC, control reg)", 1000,0,1000);
     h_m3_enri_control_1000[2] = new TH1D("m3_enri3_control","m3 (QCD enri3 MC, control reg)", 1000,0,1000);

     h_m3_tt ->Sumw2();
     h_m3_wj ->Sumw2();
     h_m3_zj ->Sumw2();
     h_m3_qcd->Sumw2();
     h_m3_vqq->Sumw2();
     h_m3_singletop->Sumw2();
     h_m3_tt_control ->Sumw2();  
     h_m3_wj_control ->Sumw2();  
     h_m3_zj_control ->Sumw2(); 
     h_m3_qcd_control->Sumw2(); 
     h_m3_vqq_control->Sumw2();
     h_m3_singletop_control->Sumw2();
     h_m3_tt_1000 ->Sumw2(); //1000 bins
     h_m3_wj_1000 ->Sumw2();
     h_m3_zj_1000 ->Sumw2();
     h_m3_qcd_1000->Sumw2();
     h_m3_vqq_1000->Sumw2();
     h_m3_singletop_1000->Sumw2();
     h_m3_tt_control_1000 ->Sumw2();  
     h_m3_wj_control_1000 ->Sumw2();  
     h_m3_zj_control_1000 ->Sumw2(); 
     h_m3_qcd_control_1000->Sumw2(); 
     h_m3_vqq_control_1000->Sumw2();
     h_m3_singletop_control_1000->Sumw2();

     for(int i=0; i<3; ++i){
       h_m3_bce[i]->Sumw2();
       h_m3_enri[i]->Sumw2();
       h_m3_bce_control[i]->Sumw2();
       h_m3_enri_control[i]->Sumw2();
       h_m3_bce_1000[i]->Sumw2();
       h_m3_enri_1000[i]->Sumw2();
       h_m3_bce_control_1000[i]->Sumw2();
       h_m3_enri_control_1000[i]->Sumw2();
     }

   }//if MC

}//end BookHistograms_wj()
//------------------------------------------------------------------------------------


void ana::BookHistograms_event_table(){
   //--------------------------------------------------
   //  Event Table (cut vs njet)
   //--------------------------------------------------
   if(m_debug) cout << "Starting BookHistograms_event_table()"<< endl;
   histf->cd();
   // Collect event count after selection
   Signal_njetsVcuts    = new TH2D("Signal_njetsVcuts",    "Events V Cuts (signal)",    ntjet+1, 0, ntjet+1, nstage, 0, nstage);
   QCD_njetsVcuts       = new TH2D("QCD_njetsVcuts",       "Events V Cuts (QCD)",       ntjet+1, 0, ntjet+1, nstage, 0, nstage);
   Wjets_njetsVcuts     = new TH2D("Wjets_njetsVcuts",     "Events V Cuts (W+jets)",    ntjet+1, 0, ntjet+1, nstage, 0, nstage);
   Zjets_njetsVcuts     = new TH2D("Zjets_njetsVcuts",     "Events V Cuts (Z+jets)",    ntjet+1, 0, ntjet+1, nstage, 0, nstage);
   VQQ_njetsVcuts       = new TH2D("VQQ_njetsVcuts",       "Events V Cuts (VQQ)",       ntjet+1, 0, ntjet+1, nstage, 0, nstage);
   SingleTop_njetsVcuts = new TH2D("SingleTop_njetsVcuts", "Events V Cuts (Single top)",ntjet+1, 0, ntjet+1, nstage, 0, nstage);
   Data_njetsVcuts      = new TH2D("Data_njetsVcuts",      "Events V Cuts (Data)",      ntjet+1, 0, ntjet+1, nstage, 0, nstage);
}
//------------------------------------------------------------------------------------


void ana::BookHistograms_btag(){

   //--------------------------------------------------
   //  BTAG
   //--------------------------------------------------
   if(m_debug) cout << "Starting BookHistograms_btag()"<< endl;
   TDirectory *dir_btag = histf->mkdir("btag","btag multiplicity (after all cuts)");
   dir_btag->cd();
   addHistoDataAndMC( h_nbtag_TCHE, "nbtag_TCHE", "nbtag (TCHEm)", 5,0,5);
   addHistoDataAndMC( h_nbtag_TCHP, "nbtag_TCHP", "nbtag (TCHPm)", 5,0,5);
   addHistoDataAndMC( h_nbtag_SSV,  "nbtag_SSV",  "nbtag (SSVm)",  5,0,5);
}
//------------------------------------------------------------------------------------


void ana::BookHistograms_PDFunc(){

   //--------------------------------------------------
   // [MC] PDF uncertainties - signal accep.eff
   //--------------------------------------------------
   if(m_debug) cout << "Starting BookHistograms_PDFunc()"<< endl;
   TDirectory *dir_pdfunc(0);
   if(!IsData() && mc_sample_has_ttbar && m_studyPDFunc ){ //for signal only
     dir_pdfunc = histf->mkdir("pdf_unc","PDF uncertainty (signal accep.eff)");
     dir_pdfunc->cd();
     h_pdf_total = new TH1F("pdf_total_sig","N(Signal event) Total",46,-1,45);
     h_pdf_pass  = new TH1F("pdf_pass_sig", "N(Signal event) Pass", 46,-1,45);
     h_pdf_eff   = new TH1F("pdf_eff_sig",  "Signal Accep.eff",     46,-1,45); 
     h_pdf_total->GetXaxis()->SetTitle("PDF set");
     h_pdf_pass->GetXaxis()->SetTitle("PDF set");
     h_pdf_eff->GetXaxis()->SetTitle("PDF set");
     h_pdf_total->SetFillColor(kGray);
     h_pdf_pass->SetFillColor(kGray);
     h_pdf_eff->SetMarkerStyle(20);
     h_pdf_total->Sumw2();
     h_pdf_pass->Sumw2();
     h_pdf_eff->Sumw2();
   }
}// end BookHistograms_PDFunc()
//------------------------------------------------------------------------------------
//====================================================================================
//
//                          End   of   Histogram   Booking
//
//====================================================================================





//====================================================================================
//
//                                MAIN   EVENT   LOOP
//
//====================================================================================
bool ana::EventLoop(){ 


   if(m_debug) cout << "Starting EventLoop()"<< endl;

   if(nfile==0) { cout << "No input file found, stop."<< endl; return false; }

   CheckAvailableJetMET();

   Init(); //initialize branch

   // Get the number of events/entries in the file chain
   Long64_t nEvents = chain->GetEntries(); 
   Long64_t nEventsAvail = nEvents;
   if(nEvents==0) { cout << "No input event found, stop." << endl; return false; }


   if(GetTrigger()) { chain->AddFriend(chain2); }


   // Read only selected branches to reduce cpu time
   ReadSelectedBranches();


   BookHistograms();








   histf->cd();




//////// Plots below to be reviewed


   //ttbar decay code from PDG codes
   TH1F *sig_all_mctype   = new TH1F("sig_all_mctype","signal decay modes (all, unweighted)",10,1,11);
   //signal_acceptance  = new TH1F("signal_acceptance","signal acceptance (unweighted) ttbar decay modes",10,1,11);

   //Electron Plots
   TH1F *h_hadOverEm_barrel = new TH1F("h_hadOverEm_barrel","H/E - Barrel",            100, 0, 0.5);
   TH1F *h_EOverPIn_barrel  = new TH1F("h_EOverPIn_barrel", "E/P (In) - Barrel",       100, 0, 10.0);
   TH1F *h_dEtaIn_barrel    = new TH1F("h_dEtaIn_barrel",   "Delta Eta (In) - Barrel", 100, -0.05, 0.05);
   TH1F *h_dPhiIn_barrel    = new TH1F("h_dPhiIn_barrel",   "Delta Phi (In) - Barrel", 100, -0.2, 0.2);
   TH1F *h_hadOverEm_endcap = new TH1F("h_hadOverEm_endcap","H/E - Endcap",            100, 0, 0.5);
   TH1F *h_EOverPIn_endcap  = new TH1F("h_EOverPIn_endcap", "E/P (In) - Endcap",       100, 0, 10.0);
   TH1F *h_dEtaIn_endcap    = new TH1F("h_dEtaIn_endcap",   "Delta Eta (In) - Endcap", 100, -0.05, 0.05);
   TH1F *h_dPhiIn_endcap    = new TH1F("h_dPhiIn_endcap",   "Delta Phi (In) - Endcap", 100, -0.2, 0.2);
   TH1F *h_cIso_barrel = new TH1F("h_cIso_barrel","Cal Iso - Barrel",20,0,10);
   TH1F *h_tIso_barrel = new TH1F("h_tIso_barrel","Trk Iso - Barrel",20,0,10);
   TH1F *h_cIso_endcap = new TH1F("h_cIso_endcap","Cal Iso - Endcap",20,0,10);
   TH1F *h_tIso_endcap = new TH1F("h_tIso_endcap","Trk Iso - Endcap",20,0,10);
   h_hadOverEm_barrel->Sumw2();
   h_EOverPIn_barrel ->Sumw2();
   h_dEtaIn_barrel   ->Sumw2();
   h_dPhiIn_barrel   ->Sumw2();
   h_hadOverEm_endcap->Sumw2();
   h_EOverPIn_endcap ->Sumw2();
   h_dEtaIn_endcap   ->Sumw2();
   h_dPhiIn_endcap   ->Sumw2();
   h_cIso_barrel->Sumw2();
   h_tIso_barrel->Sumw2();
   h_cIso_endcap->Sumw2();
   h_tIso_endcap->Sumw2();


   //Muon Plots
   TH1F *h_muon_chi2 = new TH1F("h_muon_chi2","muon chi2 (normalized) (pass pt,eta cut)",100,0,20);
   TH1F *h_muon_d0_unCor = new TH1F("h_muon_d0_unCor",  "muon d0 (uncorrected )(pass pt,eta cut) ",  1000,-0.5, 0.5);
   TH1F *h_muon_d0   = new TH1F("h_muon_d0",  "muon d0 (corrected using BeamSpot) (pass pt,eta cut)", 1000,-0.5, 0.5);
   TH1F *h_muon_hits = new TH1F("h_muon_hits","muon hits (pass pt,eta cut)",40,0,40);
   h_muon_chi2->Sumw2();
   h_muon_d0_unCor->Sumw2();
   h_muon_d0  ->Sumw2();
   h_muon_hits->Sumw2();











   


   // Histogram to store signal eff.acc
   TH1D *SignalVars = new TH1D("SignalVars","Signal acceptance and efficiency",4,0,4);







   //-----------
   // Cut names
   //-----------
   //vector<string> ve;
   ve.push_back("INITIAL   ");
   ve.push_back("GOODRUN   ");
   ve.push_back("TRIGGER   ");
   ve.push_back("$\\ge$1TEle"); //NEW
   ve.push_back("$\\ge$1TISO");
   ve.push_back("=1TISO    ");
   ve.push_back("!MUON");
   ve.push_back("MET       ");
   ve.push_back("!Z        ");
   ve.push_back("!CONV     ");
   if( !m_rejectEndcapEle ) ve.push_back("!DIFFZ    ");
   else ve.push_back("BARREL    ");
   //ve.push_back("TAGGABLE  ");
   ve.push_back("$\\ge$1+BTAG");
   ve.push_back("$\\ge$2+BTAG");
   //ve.push_back("$\\ge$1-BTAG");

   // 15 Feb 2010: find muon cut stage number
   for(int i=0; i<(int)ve.size(); i++) {
     if(ve.at(i).find("MUON")!=string::npos) m_muonCutNum=i;
   }

   // ve2 is list of cuts including njet
   ve2 = ve; //copy then insert njet cut
   ve2.insert( ve2.begin()+m_muonCutNum+1, Fourjets );


   // create a vector of nmctype elements with initial value of 0
   vector<int>               v1Dint(nmctype,0);
   vector<double>            v1Ddou(nmctype,0);
   // create a vector of ntjet elements, each elemet is a copy of v1D
   vector<vector<int> >      v2Dint(ntjet, v1Dint);
   vector<vector<double> >   v2Ddou(ntjet, v1Ddou);
   // resize to have nstage elements, each element is a copy of v2D
   e_plus_jet.resize(          nstage, v2Dint );
   e_plus_jet_weighted.resize( nstage, v2Ddou );
   v1Dint.clear();
   v1Ddou.clear();
   v2Dint.clear();
   v2Ddou.clear();



   //Vectors to store 4-vectors for selected objects
   vector<TLorentzVector> electrons;
   vector<unsigned int>   ii_electrons; //index of GoodEle
   vector<TLorentzVector> electrons_barrel;
   vector<TLorentzVector> electrons_endcap;
   vector<TLorentzVector> iso_electrons;
   vector<TLorentzVector> iso_electrons_barrel;
   vector<TLorentzVector> iso_electrons_endcap;
   vector<TLorentzVector> muons;
   vector<TLorentzVector> iso_muons;
   vector<TLorentzVector> jets;
   vector<TLorentzVector> met;
   //vector<TLorentzVector> tagged_jets; //not used
   vector<float>          electrons_isoval;
   vector<float>          electrons_isoval2;

   unsigned int index_selected_ele = 0;   //index of the selected isolated electron


   int nZ_TEST = 0; //TEST
   int nZ_ORIG = 0; //TEST



   int counter = 0;
   int counter_pass = 0;
   int printLevel = 0;


   cout << endl << " Starting! nEvents is " << nEvents << " counter is " << counter << endl;
   if(GetLimit() > 0) nEvents = GetLimit();
   cout << " Limited nEvents to " << nEvents << endl;



   // MC information (declared before event loop)
   string this_mc = "data"; //default
   int    mctype  = 0;


   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //
   //                         Loop over the events
   //
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   fprintf(outfile,"%7s %6s %10s %10s %8s %8s %6s %13s %6s %8s %8s   %s\n",
	   "no","run","event","lumiSec","tree","entry","ntj","leadingJetET","met","HT","mc", "file");
   fprintf(outfile,"---------------------------------------------------------");
   fprintf(outfile,"---------------------------------------------------------\n");
   if(m_debug) cout << "\n Will start main Event Loop"<< endl;


   // Z stuff
   int   nEvent_Z = 0;
   int   nEvent_Zee = 0;
   int   nEvent_Zmm = 0;
   int   nEvent_Ztt = 0;
   int   nEvent_2orMoreEle = 0;
   int   nEvent_EleMatch = 0;
   int   nEvent_Z_et = 0;
   int   nEvent_Z_eta = 0;
   int   nEvent_Z_eteta = 0;
   int   nEvent_Z_eteta_d0 = 0;
   int   nEvent_Z_eteta_d0_L = 0;
   int   nEvent_Z_eteta_d0_T = 0;
   int   nEvent_Z_eteta_d0_RL = 0;
   int   nEvent_Z_eteta_d0_RT = 0;
   int   nEvent_Zee_pass = 0;
   int   nEvent_Zee_highMET = 0;
   int   nEvent_Zee_lowMET = 0;
   int   nEvent_Zmm_pass = 0;
   int   nEvent_Zmm_highMET = 0;
   int   nEvent_Zmm_lowMET = 0;
   int   nEvent_Ztt_pass = 0;
   int   nEvent_Ztt_highMET = 0;
   int   nEvent_Ztt_lowMET = 0;

   // Monitor number of events with DR(e,mu) < 0.1
   int nEvent_DR_ele_muo_less_than_01 = 0;

   // record list of "interesting" files that contain selected events
   set<string> interestingFiles;


   //----------------------------------------------------------------------------
   //                   Beginning  of  Event  Loop
   //----------------------------------------------------------------------------

   for(Long64_t ev=0; ev<nEvents; ++ev) {


     // NB: LoadTree will complain (harmlessly) about unknown branch if SetBranchAddress  
     //     was called for a branch not present in the ntuple
     Long64_t lflag = chain->LoadTree(ev);
     if (lflag < 0) break;

     ++counter;
     //if(counter<900) continue;


     int nbytes = chain->GetEntry(ev);
     if(ev>100) m_debug = false; // turn off m_debug message after 100 events


     // 1- Check for good run (only applicable to data)
     bool goodrun = true;
     
     if(IsData() && KeepGoodRuns()) {
       //goodrun = mytool->good_run(run);
     }

     if(m_debug || ev<10 || ev%5000==0) {
       if(m_debug||ev<10) cout << "\nBegin processing event no. " << ev+1;
       else   	 cout << "\nBegin processing event no. " << (ev+1)/1000 << " k";
       cout << ". GoodRun=" << goodrun ;
       cout << "  << Run "<< run << ", Event "<< event << ", LumiSection " << lumiBlock << " >>  ";
       cout << printTimeNow() << endl;
       cout << " entry " << lflag 
	    << ", tree # " << chain->GetTreeNumber()  
	    << ", file " << chain->GetCurrentFile()->GetName() 
	    << ", EvtSize " << nbytes << endl;
     }
     
     // 2- Check Trigger 
     if(m_debug) cout << " Checking Trigger"<< endl;

     bool fired_single_em = true;

     if( GetTrigger() ){

       fired_single_em = false; //reset default value to "fail" if checking trigger

       if (ev==0) cout <<  "Check HLT trigger: " << HLTBit << endl; //nonIso electron trigger
       // fired_single_em = (bool)HLT_Ele15_LW_L1R;
       fired_single_em = passHLT();

       if (ev<10) { //print first 10 events
	 if (fired_single_em==true) cout << " PASS Single Electron Trigger"<< endl;
	 else  cout << " FAIL Single Electron Trigger"<< endl;
       }
     }
     

     // If running on MC, check type of MC from input filename
     // note: it is sufficient to check once per tree
     // However, if run on soup, then need to check every event
     bool RunOnMixedMC = false; //29-9-09

     if ( IsData()==false ) {

       if ( RunOnMixedMC ||    //check once per event
	    (!RunOnMixedMC && lflag==0) ) { //check once per input file

	 // Reset flags
	 isTTbar = false;
	 isWjets = false;
	 isZjets = false;
	 isQCD   = false;
	 isEnri1 = false;
	 isEnri2 = false;
	 isEnri3 = false;
	 isBce1  = false;
	 isBce2  = false;
	 isBce3  = false;
	 isVQQ   = false;
	 isSingleTop = false;
	 isTW    = false;
	 isTchan = false;
	 isSchan = false;
	 
	 if ( !RunOnMixedMC ) {
	   if(m_debug) cout << " [MC] Checking MC type from input filename"<< endl;
	   
	   //cout << "\nThis is the 1st event in the tree" << endl;
	   //cout << "Global event no. "<< ev <<". Entry number in this tree: " << lflag << endl; //TL: 14-1-09
	   const string fname( chain->GetCurrentFile()->GetName() );
	   //cout << "file: " << fname << endl;
	   
	   // figure out what type of MC from the file name
	   
	   int nmatch = 0;
	   
	   for ( size_t i=0; i < mc_names.size(); ++i ) {
	     if ( fname.find( mc_names[i] )!=string::npos ) {
	       this_mc = mc_names[i];
	       nmatch++;
	       break; //exit for loop when a match is found
	     }
	   }
	   if(nmatch==0) cout << "ERROR. File name does not match to any type of MC" << endl;
	   
	 } else {

	   //29-9-09
	   if(m_debug) cout << " [MC] Checking Event type from MC truth"<< endl;
	   //this_mc = CheckEventTypeFromMcTruth(); //TEMPORARY
	   //cout << "this_mc = " << this_mc << endl;
	   cout << "ERROR: Should not use this piece of code!" << endl;
	 }

	 
	 // set event weight of current event
	 this_weight = GetWeight(this_mc);
	 
	 
	 // set mctype code
	 if      (this_mc=="ttbar"||
		  this_mc=="ttjet"||
		  // ttnj Alpgen
		  this_mc=="tt0j"|| this_mc=="tt1j" || this_mc=="tt2j" ||
		  this_mc=="tt3j"|| this_mc=="tt4j")  {  isTTbar = true;  }  //mctype 1-10

	 // tt+j Alpgen 
	 if      (this_mc=="tt0j") mctype=23;
	 else if (this_mc=="tt1j") mctype=24;
	 else if (this_mc=="tt2j") mctype=25;
	 else if (this_mc=="tt3j") mctype=26;
	 else if (this_mc=="tt4j") mctype=27;

	 else if (this_mc=="wjet"||
		  this_mc=="wenu")  {  isWjets = true;                   mctype = 11;  }
	 else if (this_mc=="zjet"||
		  this_mc=="zee")   {  isZjets = true;                   mctype = 12;  }
	 else if (this_mc=="enri1") {  isQCD   = true;  isEnri1 = true;  mctype = 13;  }
	 else if (this_mc=="enri2") {  isQCD   = true;  isEnri2 = true;  mctype = 14;  }
	 else if (this_mc=="enri3") {  isQCD   = true;  isEnri3 = true;  mctype = 15;  }
	 else if (this_mc=="bce1")  {  isQCD   = true;  isBce1  = true;  mctype = 16;  }
	 else if (this_mc=="bce2")  {  isQCD   = true;  isBce2  = true;  mctype = 17;  }
	 else if (this_mc=="bce3")  {  isQCD   = true;  isBce3  = true;  mctype = 18;  }
	 else if (this_mc=="vqq")   {  isVQQ   = true;                   mctype = 19;  }
	 else if (this_mc=="tW")    {  isSingleTop = true; isTW    = true; mctype = 20;  }
	 else if (this_mc=="tchan") {  isSingleTop = true; isTchan = true; mctype = 21;  }
	 else if (this_mc=="schan") {  isSingleTop = true; isSchan = true; mctype = 22;  }
	 
	 if(m_debug) {
	   cout << "+++> Now starts running on  " << this_mc << "  events" << endl;
	   cout << " current entry: num " << lflag << ", tree # " << chain->GetTreeNumber()
		<< "\n filename: " <<  chain->GetCurrentFile()->GetName() << endl << endl;
	 }

       }//end if check MC
     }//end if run on MC


     
     //--------------------
     //  For TTbar MC
     //--------------------
     // FIXME (9 Mar 10): if signal is alpgen, do not check ttbar decays.
     if(!IsData() && isTTbar && !signal_is_Alpgen){ //If it isn't data figure out what kind of MC (for ttbar)
       
       if(m_debug) cout << " [MC] signal, checking type of ttbar decay." <<endl;
	 
       if(Nmc_doc==0) {
	 mctype = 1;
	 static bool first_time = true;
	 if(first_time) {
	   cout << " [MC] Nmc_doc = 0. No MC info in ntuple! Set mctype to 1 so that all ttbar events appear in ev column" << endl;
	   first_time = false;
	 }
       } else {

	 //Check ttbar decay 
	 int nlep=0;  
	 int nt=0;
	 int nW=0;
	 int nb=0;
	 int nWfromt = 0;
	 int nbfromt = 0;	 

	 if(printLevel>0) PrintGenParticles();
       
	 //Loop over mc documentation lines
	 for(unsigned int i = 0; i<Nmc_doc; ++i) { 
	   
	   if(fabs(mc_doc_id->at(i)) ==6) nt++;
	   if(fabs(mc_doc_id->at(i)) ==24) nW++;
	   if(fabs(mc_doc_id->at(i)) ==24 && fabs(mc_doc_mother_id->at(i))==6) nWfromt++;
	   if(fabs(mc_doc_id->at(i)) ==5) nb++;
	   if(fabs(mc_doc_id->at(i)) ==5 && fabs(mc_doc_mother_id->at(i))==6) nbfromt++;
	   
	   if(fabs(mc_doc_mother_id->at(i))==24){	  //Look at W daughters
	     if(fabs(mc_doc_id->at(i))==11) nlep += 1;  //add 1 for an electron
	     if(fabs(mc_doc_id->at(i))==13) nlep += 3;  //add 3 for a muon
	     if(fabs(mc_doc_id->at(i))==15) nlep += 7;  //add 7 for a tau	      
	   }
	 }
	 
	 //some reassignment for easy print-out in table form
	 if     (nlep==1)  mctype = 1; //evqq
	 else if(nlep==3)  mctype = 2; //mvqq
	 else if(nlep==7)  mctype = 3; //tvqq
	 else if(nlep==2)  mctype = 4; //evev
	 else if(nlep==6)  mctype = 5; //mvmv
	 else if(nlep==14) mctype = 6; //tvtv
	 else if(nlep==4)  mctype = 7; //evmv
	 else if(nlep==8)  mctype = 8; //evtv
	 else if(nlep==10) mctype = 9; //mvtv
	 else if(nlep==0)  mctype = 10; //qqqq
	 else cout << " ERROR:  nlep not defined " << nlep << endl;
	 	 
// 	 if(mctype >10) {
// 	   mctype=10;
// 	   cout << " ERROR: mctype is greater than 10! " << endl;
// 	 }	 

	 sig_all_mctype->Fill(mctype); //TL note: unweighted. fill with weight?

	 if(printLevel>0) cout << " ttbar code " << nlep << " " << mctype << endl;
       }//if Nmc_doc>0
     }//if ttbar MC

     //-------------------------------
     // Inspecting Z+jets M C events
     //-------------------------------
     if(!IsData() && isZjets && Nmc_doc> 0 && m_studyZveto){
       if(m_debug) cout << " Starting Inspection of Z+jets MC" << endl;

       isZee = false; //init private
       isZmm = false;
       isZtt = false;
       //vector<TLorentzVector> Zele;

       //       if(!IsData() && isZjets && Nmc_doc > 0){

	 //if(m_debug) cout << "checking Z+jets MC"<< endl;
	 //Loop over mc documentation lines, look for Z->ee
	 for(unsigned int i = 0; i<Nmc_doc; ++i) { 

	   if( mc_doc_status->at(i)!=3 ) continue; //documentation line only

	   //if( fabs(mc_doc_id->at(i))==11 && fabs(mc_doc_mother_id->at(i))==23 ) {

	   // Look for electron
	   if( fabs(mc_doc_id->at(i))==11 ) { //don't check Z parent as in some events, the Z is missing

	     isZee = true;
	     Zele.push_back( TLorentzVector( mc_doc_px->at(i), mc_doc_py->at(i), mc_doc_pz->at(i), mc_doc_energy->at(i) ) );
	     //plot eta/ET of Z electrons
	     h_Zee_eta->Fill(mc_doc_eta->at(i));
	     h_Zee_pt->Fill(mc_doc_pt->at(i));
	   
	     if(ev<10) {
	       cout << "\nZj MC:  event #" << ev << endl;
	       cout << "generator ele:  id = " << mc_doc_id->at(i) << "  pt = "<< mc_doc_pt->at(i) 
		    << "  eta = "<< mc_doc_eta->at(i) <<endl;
	     }
	   
	   }//is Z->ee
	   //if( fabs(mc_doc_id->at(i))==13 && fabs(mc_doc_mother_id->at(i))==23 ) isZmm=true;
	   //if( fabs(mc_doc_id->at(i))==15 && fabs(mc_doc_mother_id->at(i))==23 ) isZtt=true;
	   if( fabs(mc_doc_id->at(i))==13 ) isZmm=true;
	   if( fabs(mc_doc_id->at(i))==15 ) isZtt=true;       
	 }//loop over genPar
       
	 // loop over reco-photons in Z+jets events       
	 for (unsigned int i = 0; i<Nphotons; ++i) {
	   h_Z_photon_eta->Fill( photons_eta->at(i) );
	   h_Z_photon_et->Fill( photons_et->at(i) );
	   if(isZee) {
	     h_Zee_photon_eta->Fill( photons_eta->at(i) );
	     h_Zee_photon_et->Fill( photons_et->at(i) );
	     h_Zee_photon_eteta_2D->Fill( photons_et->at(i), photons_eta->at(i) );
	   }
	 }
	 h_Z_Nphotons->Fill(Nphotons);
	 if(isZee)  h_Zee_Nphotons->Fill(Nphotons); //no cut
     
	 nEvent_Z++;
	 if     (isZee) nEvent_Zee++;
	 else if(isZmm) nEvent_Zmm++;
	 else if(isZtt) nEvent_Ztt++;
	 else {
	   cout << "\nThis Z decay is not e,mu,tau! What is it?" << endl;
	   cout << "Print all MC info"<<endl;
	   PrintGenParticles();
	 }
	 if(isZee){
	   int nGenBasicEle_Zee = 0;
	   for(unsigned int i=0; i<Zele.size(); ++i){
	     if( Zele[i].Et() >30.0 && fabs(Zele[i].Eta()) < 2.5 ) nGenBasicEle_Zee++;
	   }
	   h_nGenBasicEle_Zee_allj->Fill(nGenBasicEle_Zee);
	 }
	 //}// Zjets MC
     }//if m_studyZveto


     // 2b - Get Beam Spot
     /*
       if(m_debug) cout << " Get Beam Spot" << endl;
       //     cout << "number of beam spot = " << NbeamSpot << endl;
       const float bsx = beamSpot_x->at(0);
       const float bsy = beamSpot_y->at(0);
       //     cout << "beam spot (x,y): "<< bsx << ", " << bsy << endl;
       */


     // 3- Find number of of electrons, muons, jets - after energy correctionss

     // 3A - Find number of electrons
     if(m_debug) cout << " Starting ELECTRON section"<< endl;

     int nGoodEle = 0;
     int nGoodEle_barrel = 0;
     int nGoodEle_endcap = 0;
     int nGoodIsoEle = 0;
     int nGoodIsoEle_barrel = 0;
     int nGoodIsoEle_endcap = 0;
     int nGoodNonIsoEle = 0;
     int nGoodNonIsoEle_barrel = 0;
     int nGoodNonIsoEle_endcap = 0;     
     float this_isoele_d0 = 0;

     int nHighETEle = 0; //pass et
     int nBasicEle = 0; //pass et,eta
     int nBasicD0Ele = 0; //pass et,eta,d0
     int nLooseEle = 0; //pass et,eta,d0,id
     int nTightEle = 0;
     int nRobustLooseEle = 0;
   

     //cout << "\n==> There are  " << Nels << "  reco-electrons in this event."<<endl;

     // basic
     fillHistoDataAndMC( h_nele, Nels );


     for (unsigned int i=0; i<Nels; ++i) {

       //cout << "\n--> electron no." << i+1 << endl;
       //cout << "  this electron ET/eta = " << els_et->at(i) << " / " << els_eta->at(i) << endl;

       //Compute Relative isolations
       //float RelCalIso = els_cIso->at(i)/els_et->at(i);
       //float RelTrkIso = els_tIso->at(i)/els_et->at(i);
	 
       //Compute Combined RelIso for electrons
       //       float CombRelIso = (els_dr04EcalRecHitSumEt->at(i) + els_dr04HcalTowerSumEt->at(i) + els_tIso->at(i))/els_et->at(i);// 10-9-09
       float CombRelIso = getRelIso(i); // 11-11-09
       float CombRelIso2 = els_et->at(i)/(els_et->at(i) + els_dr04EcalRecHitSumEt->at(i) + els_dr04HcalTowerSumEt->at(i) + els_tIso->at(i) ); //norm

       /*
       if(CombRelIso<0) {
	 cout << "electron no " << i << " has negative RelIso (" << CombRelIso << ")\n";
	 cout << " et  = " << els_et->at(i) << endl;
	 cout << " cIso  = " << els_dr04EcalRecHitSumEt->at(i) + els_dr04HcalTowerSumEt->at(i) << endl;
	 cout << " tIso  = " << els_tIso->at(i) << endl;
	 cout << " RelCalIso = " << RelCalIso << endl;
	 cout << " RelTrkIso = " << RelTrkIso << endl << endl;
       }
       */

       if(m_debug) cout << "-> fill electron id histogram" << endl;
       fillHistoDataAndMC( h_eid, 0.0 ); //count all electrons

       //cout << "electron ET ="<< els_et->at(i) << endl;
       //cout << "electron ET cut = "<< ELE_ETCUT << endl;
       //cout << "nHighETEle (b4) = "<< nHighETEle << endl;
       if(els_et->at(i)> ELE_ETCUT) nHighETEle++;
       //cout << "nHighETEle (after) = "<< nHighETEle << endl;


       
       // Basic kinematics (14-8-09)
       // - all electrons
       fillHistoDataAndMC( h_ele_ET[0],  els_et->at(i)  );
       fillHistoDataAndMC( h_ele_eta[0], els_eta->at(i) );
       fillHistoDataAndMC( h_ele_phi[0], els_phi->at(i) );
       fillHistoDataAndMC( h_ele_iso[0], CombRelIso     );
       if(i<3){ //first 3 electrons
         fillHistoDataAndMC( h_ele_ET[i+1],   els_et->at(i)  );
         fillHistoDataAndMC( h_ele_eta[i+1],  els_eta->at(i) );
         fillHistoDataAndMC( h_ele_phi[i+1],  els_phi->at(i) );
         fillHistoDataAndMC( h_ele_iso[i+1],  CombRelIso     );
       }
       
       /*
       // If want plots per njet, need to move the code-block to after jet cleaning
       // - all electrons
       fillHisto_Njet_DataAndMC( h_ele_ET[0],  els_et->at(i),  this_weight );
       fillHisto_Njet_DataAndMC( h_ele_eta[0], els_eta->at(i), this_weight );
       fillHisto_Njet_DataAndMC( h_ele_phi[0], els_phi->at(i), this_weight );
       fillHisto_Njet_DataAndMC( h_ele_iso[0], CombRelIso,     this_weight );
       if(i<3){ //first 3 electrons
       fillHisto_Njet_DataAndMC( h_ele_ET[i+1],   els_et->at(i),  this_weight );
       fillHisto_Njet_DataAndMC( h_ele_eta[i+1],  els_eta->at(i), this_weight );
       fillHisto_Njet_DataAndMC( h_ele_phi[i+1],  els_phi->at(i), this_weight );
       fillHisto_Njet_DataAndMC( h_ele_iso[i+1],  CombRelIso,     this_weight );
       }
       */



       //Make Et and Eta cuts 
       if (els_et->at(i) > ELE_ETCUT &&
	   fabs(els_eta->at(i)) < 2.5) {

	 if(m_debug) cout << "-> this electron has passed ET and ETA cut" << endl;

	 nBasicEle++;
	 
	 // Count how many electron passing each type of electron ID
	 if(m_debug) cout << "-> filling electron id histogram (pass Et,eta cuts)" << endl;
	 fillHistoDataAndMC( h_eid, 1.0 ); // count all ele passing et,eta cut



	 // Calculate d0 w.r.t beam spot position
	 /*
	   float vx  = els_vx->at(i);
	   float vy  = els_vy->at(i);
	   float vpx = els_vpx->at(i);
	   float vpy = els_vpy->at(i);
	   float vpt = TMath::Hypot(vpx, vpy);
	   float d0_corrected = - (-(vx - 0.0322)*vpy + vy*vpx) / vpt ; //d0 = -dxy
	   //cout << "ele d0:  uncorrected = " << els_d0->at(i) << "  corrected = "<< d0_corrected << endl;
	   */
	 float d0_corrected = compute_d0("electron",i); 



	 // Store electron d0 information
	 if(m_debug) cout << "-> filling electron d0 histograms" << endl;
	 fillHistoDataAndMC( h_ed0_unCor, fabs(els_d0dum->at(i)) );
	 fillHistoDataAndMC( h_ed0,       fabs(d0_corrected)     );


	 // (6 Mar 09) Apply d0 cut on electron, 200 micron = 0.02cm
	 if( fabs(d0_corrected) < 0.02 ) {

	   nBasicD0Ele++;
	   // Count how many electron passing each type of electron ID	  
	   if(m_debug) cout << "-> filling electron id histograms (pass d0 cut)" << endl;
	   fillHistoDataAndMC( h_eid, 2.0 );
	   if ( els_looseId->at(i) > 0 )       fillHistoDataAndMC( h_eid, 3.0 );
	   if ( els_tightId->at(i) > 0 )       fillHistoDataAndMC( h_eid, 4.0 );
	   if ( els_robustLooseId->at(i) > 0 ) fillHistoDataAndMC( h_eid, 5.0 );
	   if ( els_robustTightId->at(i) > 0 ) fillHistoDataAndMC( h_eid, 6.0 );



	   // Check isolation of electron
	   bool isIsolated = false;
	   if(useNewReliso) isIsolated = CombRelIso < 0.1; //new reliso
	   else             isIsolated = CombRelIso2 > 0.85; //old reliso
	   if (m_debug) {
	     cout << " ele  Old RelIso = " << CombRelIso2 ;
	     cout << "  New RelIso = " << CombRelIso ;
	     cout << "  isIsolated = " << isIsolated << endl;
	   }

	   
	   // Count number of good electron passing each eid (ignore those in eta gap: 1.442-1.56)
	   if ( fabs(els_eta->at(i))<1.442 || fabs(els_eta->at(i))>1.560 ) {
	     if( els_looseId->at(i) > 0 )        nLooseEle++;
	     if( els_tightId->at(i) > 0 )        nTightEle++;
	     if( els_robustLooseId->at(i) > 0 )  nRobustLooseEle++;
	   }
	   

	   //Separate into endcap and barrel
	   	   
	   if (fabs(els_eta->at(i)) < 1.442) { //Barrel

	     //Find out what isolation looks like (out of the box)
	     h_cIso_barrel->Fill(els_dr04EcalRecHitSumEt->at(i) + els_dr04HcalTowerSumEt->at(i), this_weight); //<-- add weight
	     h_tIso_barrel->Fill(els_tIso->at(i), this_weight);
	     
	     //Plot electron ID quantities	 
	     h_hadOverEm_barrel->Fill(els_hadOverEm->at(i), this_weight);
	     h_EOverPIn_barrel->Fill(els_eOverPIn->at(i), this_weight);
	     h_dEtaIn_barrel->Fill(els_dEtaIn->at(i), this_weight);
	     h_dPhiIn_barrel->Fill(els_dPhiIn->at(i), this_weight);
	   
	     //Apply "Robust Tight" Electron ID
	     
	     //	     if (els_robustTightId->at(i) == true){
	     if ( passEleID(i) ) {
	       

	       //Store 4 vector for "good" electron and increment counters
	       TLorentzVector ele(els_px->at(i),els_py->at(i),els_pz->at(i),els_energy->at(i));
	       nGoodEle_barrel++;	       
	       nGoodEle++;
	       electrons.push_back(ele);
	       ii_electrons.push_back(i); //index
	       electrons_barrel.push_back(ele);
	       electrons_isoval.push_back(CombRelIso); //needed by QCD estimation
	       electrons_isoval2.push_back(CombRelIso2); //needed by QCD estimation

	       if (isIsolated) { //Isolated
		 
		 //Store 4 vector for isolated electron and increment counters
		 nGoodIsoEle_barrel++;	       
		 nGoodIsoEle++;	       
		 iso_electrons.push_back(ele);
		 iso_electrons_barrel.push_back(ele);
		 this_isoele_d0 = fabs(d0_corrected);
		 index_selected_ele = i;
	       }
	       else { //Non-Isolated
		 
		 nGoodNonIsoEle_barrel++;	       
		 nGoodNonIsoEle++;	       
	       }

	     }	     
	   }// end barrel
	   else if (fabs(els_eta->at(i)) > 1.560 && fabs(els_eta->at(i)) < 2.5) { //Endcap

	     //Find out what isolation looks like (out of the box)
	     h_cIso_endcap->Fill(els_dr04EcalRecHitSumEt->at(i) + els_dr04HcalTowerSumEt->at(i), this_weight);
	     h_tIso_endcap->Fill(els_tIso->at(i), this_weight);
	       
	     //Plot electron ID quantities	 
	     h_hadOverEm_endcap->Fill(els_hadOverEm->at(i), this_weight);
	     h_EOverPIn_endcap->Fill(els_eOverPIn->at(i), this_weight);
	     h_dEtaIn_endcap->Fill(els_dEtaIn->at(i), this_weight);
	     h_dPhiIn_endcap->Fill(els_dPhiIn->at(i), this_weight);
	     
	     //Apply "Robust Tight" Electron ID	       
	     // if (els_robustTightId->at(i) > 0 ) {
	     if ( passEleID(i) ) {

	       TLorentzVector ele(els_px->at(i),els_py->at(i),els_pz->at(i),els_energy->at(i));
	       nGoodEle_endcap++;	       
	       nGoodEle++;
	       electrons.push_back(ele);
	       ii_electrons.push_back(i); //index
	       electrons_endcap.push_back(ele);	       
	       electrons_isoval.push_back(CombRelIso); //needed by QCD estimation
	       electrons_isoval2.push_back(CombRelIso2); //needed by QCD estimation

	       //if (CombRelIso2 > 0.85) { //Isolated (endcap)
	       if (isIsolated) { //Isolated (endcap)

		 //Store 4 vector for isolated electron and increment counters
		 nGoodIsoEle_endcap++;	       
		 nGoodIsoEle++;	       
		 iso_electrons.push_back(ele);
		 iso_electrons_endcap.push_back(ele);
		 this_isoele_d0 = fabs(d0_corrected);
		 index_selected_ele = i;
	       }
	       else { //Non-Isolated

		 nGoodNonIsoEle_endcap++;	       
		 nGoodNonIsoEle++;	       
	       }	       
	       
	     }//end if pass eidRobustTight
	     
	   }// end endcap	  	   
	   else { //Electrons in Gap
	   }
	 }// end d0 cut
       }// end et,eta cuts
     }// end electrons loop 

     if(m_debug){
       cout << "\n-- summary -- " << endl;
       cout << "Nels   = " << Nels << endl;
       cout << "nHighETEle = " << nHighETEle << endl;
       cout << "nBasicEle  = " << nBasicEle << endl;
       cout << "nBasicD0Ele = " << nBasicD0Ele << endl;
       cout << "nLooseEle  = " << nLooseEle << endl;
       cout << "nTightEle  = " << nTightEle << endl;
       cout << "nRobustLooseEle = " << nRobustLooseEle << endl;
       cout << "nGoodEle = " << nGoodEle << endl;
     }

     // (Z->ee only) for event with exactly 2 reco-ele, match them to gen Ele from Z decay
     if(isZjets && isZee) {
       int nEleMatch = 0;
       vector<int> recoElectronMatchToZee_index;
       if( Nels >= 2 ){
	 
	 nEvent_2orMoreEle++;
	 for (unsigned int i=0; i<Nels; ++i) { //loop over reco-ele
	   //cout << " reco-ele no." << i <<  endl;
	   
	   TLorentzVector this_reco_ele( els_px->at(i), els_py->at(i), els_pz->at(i), els_energy->at(i) );
	   
	   for (unsigned int j=0; j< Zele.size(); ++j) { //loop over gen-ele from Z
	     //cout << " gen-ele no." << j <<  endl;
	     
	     float delR_recoE_genE = calcDeltaR( this_reco_ele, Zele.at(j));
	     //cout << "DR(reco " << i << ",gen " << j << ") = "<< delR_recoE_genE << endl;
	     if( delR_recoE_genE < 0.1 ) {
	       nEleMatch ++;
	       recoElectronMatchToZee_index.push_back(i);
	     }
	   }
	 }// loop over ele
       }// 2 or more ele

       int nEle_Z_et = 0;
       int nEle_Z_eta = 0;
       int nEle_Z_eteta = 0;
       int nEle_Z_eteta_d0 = 0;
       int nEle_Z_eteta_d0_RL = 0;
       int nEle_Z_eteta_d0_L = 0;
       int nEle_Z_eteta_d0_T = 0;
       int nEle_Z_eteta_d0_RT = 0;

       if (nEleMatch==2) {
	 nEvent_EleMatch++;      
	 //cout << "info: this event has 2 reo-electrons matched to Z->ee"<< endl;
	 for(int i=0; i<2; ++i){
	   int k = recoElectronMatchToZee_index.at(i);
	 
	   if( els_et->at(k) > 30.0 )	  nEle_Z_et++; //et cut only
	   if( fabs(els_eta->at(k)) < 2.5 ) nEle_Z_eta++; //eta cut only

	   if ( els_et->at(k) > 30.0 && fabs(els_eta->at(k)) < 2.5 ){ //et and eta cut
	     nEle_Z_eteta++;
	     
	     // Calculate d0 w.r.t beam spot position
	     float d0_corrected = compute_d0("electron",k);
	     
	     if( fabs(d0_corrected) < 0.02 ){
	       nEle_Z_eteta_d0++;
	     
	       if(els_robustLooseId->at(k)) nEle_Z_eteta_d0_RL++;
	       if(els_robustTightId->at(k)) nEle_Z_eteta_d0_RT++;
	       if(els_tightId->at(k))       nEle_Z_eteta_d0_T++;
	       if(els_looseId->at(k))       nEle_Z_eteta_d0_L++;
	     }//d0
	   }//et,eta
	 }//ele loop
       }//if nEleMatch==2

       // count how many events with 2 or more electrons (pass each cut)
       if(nEle_Z_et>=2)          nEvent_Z_et++;
       if(nEle_Z_eta>=2)         nEvent_Z_eta++;
       if(nEle_Z_eteta>=2)       nEvent_Z_eteta++;
       if(nEle_Z_eteta_d0>=2)    nEvent_Z_eteta_d0++;
       if(nEle_Z_eteta_d0_L>=2)  nEvent_Z_eteta_d0_L++;
       if(nEle_Z_eteta_d0_T>=2)  nEvent_Z_eteta_d0_T++;
       if(nEle_Z_eteta_d0_RT>=2) nEvent_Z_eteta_d0_RT++;
       if(nEle_Z_eteta_d0_RL>=2) nEvent_Z_eteta_d0_RL++;
     }//Z->ee MC



     // 3B - Find number of muons        
     if(m_debug) cout << " Starting MUON section"<< endl;

     int nGoodMu = 0;
     int nGoodIsoMu = 0;
     int nGoodNonIsoMu = 0;
       

     for(unsigned int i = 0; i<Nmus; ++i) {
       
       //Make Pt and Eta cuts (consider global muons only)
       if ( mus_id_AllGlobalMuons->at(i) > 0 && 
	    mus_cm_pt->at(i) > MU_PTCUT &&
	    fabs(mus_cm_eta->at(i)) < 2.1 ) {
	 
	 // Correct muon d0 using BeamSpot
	 float mu_d0_corrected = compute_d0("muon",i);


	 //Plot Muon ID quantities	 
	 h_muon_chi2->Fill( mus_cm_chi2->at(i)/mus_cm_ndof->at(i), this_weight);  // "cm" means combinedMuon
	 h_muon_d0_unCor->Fill(mus_cm_d0dum->at(i), this_weight );
	 h_muon_d0->Fill(     mu_d0_corrected,   this_weight );
	 h_muon_hits->Fill(   mus_tkHits->at(i), this_weight );
	   
	 //Apply V+jets Muon ID (consider only Global Muons)
	 if ( (mus_cm_chi2->at(i)/mus_cm_ndof->at(i)) < 10 &&
	      fabs(mu_d0_corrected) < 0.02 &&
	      mus_tkHits->at(i) >= 11 &&
	      mus_hcalvetoDep->at(i) < 6.0 &&  // veto cone energy
	      mus_ecalvetoDep->at(i) < 4.0 ) {  //what else?

	   //Store 4 vector for "good" muon and increment counters
	   TLorentzVector muo(mus_cm_px->at(i),mus_cm_py->at(i),mus_cm_pz->at(i),mus_energy->at(i)); //global muon

	   nGoodMu++;	       
	   muons.push_back(muo);
	     
	     
	   //Compute Combined RelIso for muons
	   float CombRelIso = (mus_tIso->at(i) + mus_cIso->at(i) ) / mus_cm_pt->at(i);  //new reliso
	   float CombRelIso2 = 1/(1+CombRelIso); //norm (old)

	   // Apply isolation cut
	   bool isIsolated = false;
	   if(useNewReliso) isIsolated = CombRelIso < 0.05; //new reliso
	   else             isIsolated = CombRelIso2 > 0.85; //old reliso
	   if (m_debug) {
	     cout << " muon  Old RelIso = "<< CombRelIso2 ;
	     cout << "  New RelIso = " << CombRelIso ;
	     cout << "  isIsolated = " << isIsolated << endl;
	   }

	   if (isIsolated) { //Isolated

	     //Store 4 vector for isolated muon and increment counters
	     nGoodIsoMu++;	       
	     iso_muons.push_back(muo);
	   }
	   else { //Non-Isolated
	     nGoodNonIsoMu++;
	   }
	     
	 }// end muon ID cut
       }// end pt,eta cut
     }// end muons loop

     

     // 3C - Set muon, Z veto, and conversion flags
     if(m_debug) cout << " Starting MUON Veto section"<< endl;
     bool isMuon = false;
     if (nGoodIsoMu > 0)
       isMuon = true;

     if(m_debug) cout << " Starting CONVERSION Veto section"<< endl;
     bool isConversion = false;
     if (nGoodIsoEle == 1) { // If 1 electron only, check if this electron is a conversion ... if more than one ignore (will be rejected anyawy)
       isConversion = ConversionFinder(iso_electrons.at(0), mctype, index_selected_ele);
       if( DoConversionStudies() ){ConversionMCMatching(iso_electrons.at(0), mctype, isConversion);}
       //OptimiseConversionFinder(iso_electrons.at(0), mctype);
     }

     // Is the fist isolated electron in Barrel?
     bool isBarrel = false;
     if( iso_electrons.size()>0 && fabs(iso_electrons.at(0).Eta())<1.442  ) isBarrel = true;


     //-----------------------------------------------------------------------------------
     //  Z veto (normal selection)
     //-----------------------------------------------------------------------------------
     // Two components:
     //  (1) If event contain 2 or more isolated electrons, flag as Z.
     //  (2) If event has just 1 isolated electron, plus a robustLoose electron with 20 GeV,
     //      and if the di-electron invariant mass falls in window 75-106 GeV, flag as Z.
     //-----------------------------------------------------------------------------------
     if(m_debug) cout << " Starting Z Veto section"<< endl;

     // Component (1)
     bool isZ_twoE = false;
     if (nGoodIsoEle >= 2) {
       isZ_twoE = true;       
       if (nGoodIsoEle >= 3) cout << "*** Multilepton (>=3) Event!!! ***" << endl;
     }
     

     // Component (2)
     bool isZ_mee_NEW = false;

     if ( nGoodIsoEle==1 && Nels > 1  ) {

       //cout << "TEST: N(e) = " << Nels << endl;
       for (unsigned int j=0; j<Nels; ++j){ //e loop
	 
	 if(j==index_selected_ele) continue;
	   
	 //consider only 2nd electron with ET>20 GeV, eta<2.5, d0<200um, RobustLoose
	 if( els_et->at(j) < 20.0 ) continue;
	 if( fabs( els_eta->at(j) ) > 2.5 ) continue;
	 if( fabs(compute_d0("electron",j)) > 0.02 ) continue;
	 if( els_robustLooseId->at(j) < 1 ) continue;
	 
	 TLorentzVector loose(els_px->at(j),els_py->at(j),els_pz->at(j),els_energy->at(j));
	 float mass = ( iso_electrons.at(0) + loose ).M();
	 //cout << "   m(e,e) = " << mass << endl;
	 if(m_studyZveto) fillHistoDataAndMC( h_mass_diele_new, mass );
	 if ( mass > 76  &&  mass < 106 ) isZ_mee_NEW = true; //within window, flag
	 
       }//e loop

     }//if more than 2 electron
     
     if(isZ_mee_NEW==true) nZ_TEST++;
     


     //----------------
     ///1310 TEST  13-10-09
     //----------------     
     /*
     int  test_nRL = 0;
     for ( unsigned int i=0; i<Nels; ++i) {
       if ( els_et->at(i) > 20.0 && 
            fabs(els_eta->at(i)) < 2.5 && 
            els_robustLooseId->at(i) > 0 ) ++test_nRL;
     }
     bool isZ_twoRL = (test_nRL >= 2) ;
     */
     //bool isZ = isZ_twoE;
     //if(test_nRL>=2) isZ = true;
     
     //float mass_ee = -1;
     
     //1310 ----------


     
     //-----------------------------------------------------------------------------------
     //  Z veto (normal selection)   [ Original ]
     //-----------------------------------------------------------------------------------
     // 13 Mar 09: ask for 1 robustTight electron ("good"), plus a reco-electron (not cut)
     // ie ask for 2 reco-el, and then require one to be tight (also pass et,eta,d0 cuts)
     //-----------------------------------------------------------------------------------

     if(m_debug) cout << " Apply Z veto" << endl;

     bool isZ_mee_ORI = false; //mee = mass(e,e)   
     float mass_ee = -1;

     
     if ( nGoodIsoEle==1  &&  Nels > 1 ) { //better (because there are Z events with 3 ele, one is low pt (243evts/32995 

       // (a) if there are only 1 RT electron, take another electron (most energetic one if there 
       //     are more than one.
       if ( nGoodEle==1 ) {
	 // check index of the RT electron	 
	 TLorentzVector e1(electrons.at(0)); //RT
	 int index = -1;
	 for (unsigned int i=0; i<Nels; ++i){ //e loop
	   if( fabs(e1.Px() - els_px->at(i)) < 1e-6 ) index = i;//match
	 }
	 if(index<0) { 
	   cout << "WARNING: cannot find index of the RobustTight electron!"<<endl;
	 }
	 if(m_debug) { 
	   cout << "electrons.size = " << electrons.size() << endl;
	   cout << "electrons.at(0) px py pz E = " 
		<< electrons.at(0).Px() << ", "
		<< electrons.at(0).Py() << ", "
		<< electrons.at(0).Pz() << ", "
		<< electrons.at(0).E() << endl;
	   cout << "e1.px py pz E = " << e1.Px() << ", "<< e1.Py() << ", "<< e1.Pz() <<", "<< e1.E() << endl;
	   for (unsigned int i=0; i<Nels; ++i){ //e loop
	     cout <<  "   els_px(py,pz,E)->at(" << i << ") = " <<els_px->at(i) ;
	     cout << ", " << els_py->at(i) << ", " << els_pz->at(i) << ", " << els_energy->at(i) << endl;
	   }
	 }
	 vector<TLorentzVector> loose_electrons;
	 if(index>=0) { //bug-fix Z veto (Aug09)
	   //find 2nd electron passing ET>20, |eta|<2.5, d0, RobustLoose
	   for(unsigned int i=0; i<Nels; ++i){
	     if( (i-index)==0 ) continue; //skip the Good Electron
	     if( els_et->at(i) > 20  &&  fabs(els_eta->at(i)) < 2.5  &&
                 fabs(compute_d0("electron",i)) < 0.02 &&
		 els_robustLooseId->at(i)==true ) {
	       loose_electrons.push_back(TLorentzVector(els_px->at(i),els_py->at(i),els_pz->at(i),els_energy->at(i)));
	     }
	   }
	 }
	 if( loose_electrons.size()>0 ) {
	   float mass = (e1 + loose_electrons.at(0) ).M();
	   mass_ee = mass;
	 }
       }

       // (b) if there are >=2 RT electrons, use them to calculate inv mass(e,e)
       //     - if there are more than 2, use the first 2.
       else if ( nGoodEle>=2 ) {
	 float mass = ( electrons.at(0) + electrons.at(1) ).M();
	 mass_ee = mass;
       }
     }// 2 reco-e, of which 1 is RT
     if(mass_ee>0 && m_studyZveto) fillHistoDataAndMC( h_mass_diele, mass_ee );
     if( mass_ee > 76  &&  mass_ee < 106 ) isZ_mee_ORI = true;
     //if(isZ_mee) isZ = true;
     
     if(isZ_mee_ORI) nZ_ORIG++;



     // Choose which Z-veto to use
     //----------------------------
     bool isZ_mee ;
     if( use_old_Z_veto ) isZ_mee = isZ_mee_ORI; //old
     else                 isZ_mee = isZ_mee_NEW; //new     
     
     const bool isZ = isZ_twoE || isZ_mee;
    






     // 3D - Find number of good Primary Vertices
     // Ignore for now, PV currently not in ntuple?
     
     // ---->  Add PV here?  [ Aug 09 ]


     //cout << " Quality " << zv->quality << " Vertex z = " << zv->z_pos << " +- " << zv->z_err 
     //    << " Ntrk " << zv->n_trk << " Sum_pt " << zv->sum_pt << endl;	  




       

     // 3E - Find number of jets
     if(m_debug) cout << " Starting JET section"<< endl;
     // int nGoodJet =0;
     m_nGoodJet = 0; //reset
     int ntj = 0;


     if (m_jetAlgo=="Default") {
       for(unsigned int i = 0; i<Njets; ++i) {
	 if ( jets_pt->at(i) > JET_PTCUT   &&
	      fabs(jets_eta->at(i)) < 2.4  &&
	      jets_emf->at(i) > 0.01 )  { //<-- NEW: jets_emf
	   TLorentzVector jt(jets_px->at(i),jets_py->at(i),jets_pz->at(i),jets_energy->at(i));
	   m_nGoodJet++;
	   jets.push_back(jt);

	 }
       } 
     } else if (m_jetAlgo=="SC5") {
       for(unsigned int i = 0; i<NjetsSC5; ++i) {
	 if (jetsSC5_pt->at(i) > JET_PTCUT &&
	     fabs(jetsSC5_eta->at(i)) < 2.4 &&
	     jetsSC5_emf->at(i) > 0.01 )  {
	   TLorentzVector jt(jetsSC5_px->at(i),jetsSC5_py->at(i),jetsSC5_pz->at(i),jetsSC5_energy->at(i));
	   m_nGoodJet++;
	   jets.push_back(jt);
	 }
       }
     } else if (m_jetAlgo=="SC7") {
       for(unsigned int i = 0; i<NjetsSC7; ++i) {
	 if (jetsSC7_pt->at(i) > JET_PTCUT &&
	     fabs(jetsSC7_eta->at(i)) < 2.4 &&
	     jetsSC7_emf->at(i) > 0.01 )  {
	   TLorentzVector jt(jetsSC7_px->at(i),jetsSC7_py->at(i),jetsSC7_pz->at(i),jetsSC7_energy->at(i));
	   m_nGoodJet++;
	   jets.push_back(jt);
	 }
       }
     } else if (m_jetAlgo=="KT4") {
       for(unsigned int i = 0; i<NjetsKT4; ++i) {
	 if (jetsKT4_pt->at(i) > JET_PTCUT &&
	     fabs(jetsKT4_eta->at(i)) < 2.4 &&
	     jetsKT4_emf->at(i) > 0.01 )  {
	   TLorentzVector jt(jetsKT4_px->at(i),jetsKT4_py->at(i),jetsKT4_pz->at(i),jetsKT4_energy->at(i));
	   m_nGoodJet++;
	   jets.push_back(jt);
	 }
       }
     } else if (m_jetAlgo=="KT6") {
       for(unsigned int i = 0; i<NjetsKT6; ++i) {
	 if (jetsKT6_pt->at(i) > JET_PTCUT &&
	     fabs(jetsKT6_eta->at(i)) < 2.4 &&
	     jetsKT6_emf->at(i) > 0.01 )  {
	   TLorentzVector jt(jetsKT6_px->at(i),jetsKT6_py->at(i),jetsKT6_pz->at(i),jetsKT6_energy->at(i));
	   m_nGoodJet++;
	   jets.push_back(jt);
	 }
       }
     } else if (m_jetAlgo=="JPTAK5") {
       for(unsigned int i = 0; i<NjetsJPTAK5; ++i) {
	 if (jetsJPTAK5_pt->at(i) > JET_PTCUT &&
	     fabs(jetsJPTAK5_eta->at(i)) < 2.4 &&
	     jetsJPTAK5_emf->at(i) > 0.01 )  {
	   TLorentzVector jt(jetsJPTAK5_px->at(i),jetsJPTAK5_py->at(i),jetsJPTAK5_pz->at(i),jetsJPTAK5_energy->at(i));
	   m_nGoodJet++;
	   jets.push_back(jt);
	 }
       }
     } else if (m_jetAlgo=="pfjet") {
       for(unsigned int i = 0; i<NPFJets; ++i) {
	 if (PFJets_pt->at(i) > JET_PTCUT &&
	     fabs(PFJets_eta->at(i)) < 2.4 ) {// NB: EMF cannot be calculated for PFJets
	   TLorentzVector jt(PFJets_px->at(i),PFJets_py->at(i),PFJets_pz->at(i),PFJets_energy->at(i));
	   m_nGoodJet++;
	   jets.push_back(jt);
	 }
       }
     } else {
       cout << "ERROR: wrong jet collection." << endl;      
     }



     // 3F - Do "jet cleaning"
     if(m_debug) cout << " Do jet cleaning" << endl;
	
     float delR_jet_ele = 9999;
	
     //Loop over selected jets
     for(unsigned int i = 0 ; i < jets.size(); ++i) {

       //Loop over selected (isolated) electrons
       //for(unsigned int j = 0 ; j < iso_electrons.size(); ++j) {	  
       for(unsigned int j = 0 ; j < electrons.size(); ++j) {

	 //For each jet, compute dR to *all* selected electrons
	 delR_jet_ele = calcDeltaR(jets.at(i),electrons.at(j));
	 //delR_jet_ele = calcDeltaR(jets.at(i),iso_electrons.at(j));
    
	 //If delta R is within a cone of 0.3, remove this jet from selected jet list, since it's an electron
	 if (delR_jet_ele < 0.3) {
	   m_nGoodJet--;
	   jets.erase(jets.begin() + i);
	   i--;
	   break;
	 }	    
       }
     }// end jets loop
	

     // (TL) NEW: 14 Aug 09
     // basic plot: jet collection after cleaning
     fillHistoDataAndMC( h_njet, float(m_nGoodJet) );

     //m_nGoodJet = m_nGoodJet; //set private var for this event
     // if( Njet()<0 ) { cout << "ERROR!!! negative njet: " << Njet() << endl; }


     for(unsigned int i=0; i < jets.size(); ++i) {

       fillHistoDataAndMC( h_jet_PT[0],  jets.at(i).Pt() );
       fillHistoDataAndMC( h_jet_eta[0], jets.at(i).Eta() );
       fillHistoDataAndMC( h_jet_phi[0], jets.at(i).Phi() );
       if (i<4) {// fill for first 4 jets
         fillHistoDataAndMC( h_jet_PT[i+1],  jets.at(i).Pt() );
         fillHistoDataAndMC( h_jet_eta[i+1], jets.at(i).Eta() );
         fillHistoDataAndMC( h_jet_phi[i+1], jets.at(i).Phi() );
       }       
     }


     // fill electron-counting histo
     fillHisto_Njet_DataAndMC( h_nEle_all,        Nels,            1.0 ); //all e
     fillHisto_Njet_DataAndMC( h_nEle_s1,         nBasicEle,       1.0 ); //pass et,eta
     fillHisto_Njet_DataAndMC( h_nEle_s2,         nBasicD0Ele,     1.0 ); //also pass d0
     fillHisto_Njet_DataAndMC( h_nEle_s3_idLoose, nLooseEle,       1.0 ); //also pass eidL
     fillHisto_Njet_DataAndMC( h_nEle_s3_idTight, nTightEle,       1.0 ); //also pass eidT
     fillHisto_Njet_DataAndMC( h_nEle_s3_idRL,    nRobustLooseEle, 1.0 ); //also pass eidRL
     fillHisto_Njet_DataAndMC( h_nEle_s3_idRT,    nGoodEle,        1.0 ); //also pass eidRT



     // 3G - Order jets in Et
	
     // IGNORE FOR NOW - not necessary
	
     // 3H - Check z of lepton within to z of PV (within 2cm)?
	
     // IGNORE FOR NOW - no PV info in ntuples
	
     bool isDifferentInteraction = false;
	

     // 4 - Find numbers of taggable and tagged jets  
     if(m_debug) cout << " Starting BTAG section" << endl;

     //int ntaggable = 0;
     //int nbtagP    = 0;
     //int nbtagN    = 0;  
     DoBTagging(electrons); //ele needed for jet-cleaning


	
     // 5 Get MET

     if(m_debug) cout << " Starting MET section" << endl;

     if (Nmets > 1) cout << "WARNING: more than 1 met object, using index = 0 one" << endl;
	
     // NOTE, filling met4v (MET_x,MET_y,0, sqrt(MET_x^2 + MET_y^2))
     // e.g. to get |met| = met4v.Pt();
     // e.g to get met_phi, do arccos(met4v.px()/met4v.E()))
	
     //TLorentzVector met4v(mets_ex->at(0),mets_ey->at(0),0,mets_et->at(0));
     //met.push_back(met4v);
     //const float this_met = met.at(0).Pt();

     // Use muon-corrected MET
     //const float this_met     = mets_et_muonCor->at(0);
     //const float this_met_phi = mets_phi_muonCor->at(0);
     float this_met     = 0;
     float this_met_phi = 0;


     if (m_metAlgo=="Default") {
       this_met     = mets_et_muonCor->at(0);
       this_met_phi = mets_phi_muonCor->at(0);


     } else if (m_metAlgo=="calomet_mujes") {
       this_met     = mets_et->at(0);
       this_met_phi = mets_phi->at(0);

     } else if (m_metAlgo=="SC5") {
       this_met     = metsSC5_et->at(0);
       this_met_phi = metsSC5_phi->at(0);

     } else if (m_metAlgo=="SC7") {
       this_met     = metsSC7_et->at(0);
       this_met_phi = metsSC7_phi->at(0);

     } else if (m_metAlgo=="KT4") {
       this_met     = metsKT4_et->at(0);
       this_met_phi = metsKT4_phi->at(0);

     } else if (m_metAlgo=="KT6") {
       this_met     = metsKT6_et->at(0);
       this_met_phi = metsKT6_phi->at(0);

     } else if (m_metAlgo=="tcmet") {
       this_met     = tcmets_et->at(0);
       this_met_phi = tcmets_phi->at(0);

     } else if(m_metAlgo=="pfmet") {
       this_met     = PFMets_et->at(0);
       this_met_phi = PFMets_phi->at(0);
     }

     TVector2 met2v_mu;
     TVector2 met2v_t1;

     met2v_mu.SetMagPhi( mets_et_muonCor->at(0), mets_phi_muonCor->at(0) );
     met2v_t1.SetMagPhi( mets_et->at(0),         mets_phi->at(0) );
     
     // genMET
     float this_genmet = 0;
     TVector2 met2v_gen;
     if(!IsData()) { 
       this_genmet = mets_gen_et->at(0); 
       met2v_gen.SetMagPhi( mets_gen_et->at(0), mets_gen_phi->at(0) );
     }
     // MET eta? Not meaningful.
     //TLorentzVector met4p(this_met*cos(this_met_phi),this_met*sin(this_met_phi),0,this_met);     
     //float this_met_eta = met4p.Eta();
     //cout << "MET px = " << this_met*cos(this_met_phi) << endl;
     //cout << "MET py = " << this_met*sin(this_met_phi) << endl;
     //cout << "MET eta = " << this_met_eta << endl;

     fillHistoDataAndMC( h_metAlone,     this_met     );
     fillHistoDataAndMC( h_metAlone_phi, this_met_phi );




     // 6 - Calculate some kinematic variables
	
     if(m_debug) cout << " Starting HT calculation section" << endl;

     // HT =scalar sum of MET, electron PT, muon PT and jet PT
     float ht = this_met;
     for (unsigned int i=0; i<iso_electrons.size(); ++i) {
       ht += iso_electrons.at(i).Pt();
     }
     for (unsigned int i=0; i<iso_muons.size(); ++i) {
       ht += iso_muons.at(i).Pt(); 
     }
     for (unsigned int i=0; i<jets.size(); ++i) {
       ht += jets.at(i).Pt();
     }
     
     // compute mT(W)
     //     float this_mtw = -1;
     float this_mu_mtw = -1;
     float this_t1_mtw = -1;
     float this_mu_DPhiEmet = -1;
     float this_t1_DPhiEmet = -1;

     if(iso_electrons.size()>0) {

       TVector2 e2v( iso_electrons.at(0).Px(), iso_electrons.at(0).Py() );

       this_mu_mtw = compute_mtw( e2v, met2v_mu );      // muon-corr. caloMET
       this_t1_mtw = compute_mtw( e2v, met2v_t1 );      // type 1 caloMET
       fillHistoDataAndMC( h_mtw_mu_incl, this_mu_mtw );
       fillHistoDataAndMC( h_mtw_t1_incl, this_t1_mtw ); 

       /*
       cout << "\n\nele:     px = " << e2v.X()  << ",  py = " << e2v.Y() << ",  phi = " << e2v.Phi()<< endl;
       cout << "miss (mu): px = " << met2v_mu.X() << ",  py = " << met2v_mu.Y() << ",  phi = " << met2v_mu.Phi() << endl;
       float dphi = e2v.Phi() - met2v_mu.Phi() ;
       cout << "phi(ele) - phi(miss) = " << dphi << endl;
       if(fabs(dphi) > TMath::Pi() ) dphi = 2*TMath::Pi() - fabs(dphi);
       cout << " correct = " << dphi << endl;
       */
       //       TVector2 eee(iso_electrons.at(0).P() ;

       this_mu_DPhiEmet = ( e2v ).DeltaPhi( met2v_mu );
       this_t1_DPhiEmet = ( e2v ).DeltaPhi( met2v_t1 );
       /*
       cout  << "\nthis_mu_DPhiEmet = " << this_mu_DPhiEmet << endl;
       cout  << " ele phi = " << iso_electrons.at(0).Phi() << endl;
       cout  << " miss phi (mu) = "<< met2v_mu.Phi() << endl;
       cout  << " miss phi (t1) = "<< met2v_t1.Phi() << endl;
       */
       fillHistoDataAndMC( h_DPhiEmet_mu_incl, this_mu_DPhiEmet );
       fillHistoDataAndMC( h_DPhiEmet_t1_incl, this_t1_DPhiEmet );

     }


     ntj = TMath::Min(m_nGoodJet,4); //max of ntj is 4
     //if(m_nGoodJet > 4) ntj = 4;



     //------------------------------------------------------------------------------------
     //81FB  produce validation plots
     //-------------------------------
     if( m_doValidation ) {
       if(m_debug) cout << " Produce validation plots" << endl;
       const bool Boolcuts[9] = {1,fired_single_em,(nGoodEle>0),(nGoodIsoEle > 0),(nGoodIsoEle == 1),
				 !isMuon,(this_met > METCUT),!isZ,!isConversion};
       
       valid_fillHisto(valid_HT, Boolcuts, ht);
       
       for(unsigned int i = 0 ; i < jets.size(); ++i) {
	 valid_fillHisto(valid_jetsEt,  Boolcuts, jets.at(i).Et());	
	 valid_fillHisto(valid_jetsEta, Boolcuts, jets.at(i).Phi());
	 valid_fillHisto(valid_jetsPhi, Boolcuts, jets.at(i).Eta());
       }
              
       if(ntj>0) valid_fillHisto(valid_jets1stEt, Boolcuts, jets.at(0).Et());
       if(ntj>1) valid_fillHisto(valid_jets2ndEt, Boolcuts, jets.at(1).Et());
       if(ntj>2) valid_fillHisto(valid_jets3rdEt, Boolcuts, jets.at(2).Et());
       if(ntj>3) valid_fillHisto(valid_jets4thEt, Boolcuts, jets.at(3).Et());
       

       
       vector<TLorentzVector> valid_eles;
       for (unsigned int i = 0; i<Nels; ++i) {
	 bool eleBoolcuts[9] = {1,fired_single_em,(nGoodEle>0),(nGoodIsoEle > 0),(nGoodIsoEle == 1),
				!isMuon,(this_met > METCUT),!isZ,!isConversion};
	 
	 float d0_corrected = fabs(compute_d0("electron",i)); //abs
	 float RelCalIso = (els_dr04EcalRecHitSumEt->at(i) + els_dr04HcalTowerSumEt->at(i))/els_et->at(i);
	 float RelTrkIso = els_tIso->at(i)/els_et->at(i);
	 
	 if( (RelCalIso+RelTrkIso) > 0.1 ) { eleBoolcuts[3] = 0; }
	 if(els_et->at(i) < ELE_ETCUT || d0_corrected > 0.02) { eleBoolcuts[2] = 0; }
	 //if(els_robustTightId->at(i) == false || fabs(els_eta->at(i))>2.5 ) { eleBoolcuts[2] = 0; }
	 if ( passEleID(i)==false || fabs(els_eta->at(i))>2.5 ) { eleBoolcuts[2] = 0; }

	 valid_fillHisto(valid_eleEt,     eleBoolcuts,  els_et->at(i));
	 valid_fillHisto(valid_eleEta,    eleBoolcuts,  els_eta->at(i));      
	 valid_fillHisto(valid_elePhi,    eleBoolcuts,  els_phi->at(i));
	 valid_fillHisto(valid_eled0,     eleBoolcuts,  d0_corrected);
	 valid_fillHisto(valid_eleCalIso, eleBoolcuts,  RelCalIso);
	 valid_fillHisto(valid_eleTrkIso, eleBoolcuts,  RelTrkIso);
	 valid_fillHisto(valid_eleRelIso, eleBoolcuts,  (RelCalIso+RelTrkIso));

	 TLorentzVector eles(els_px->at(i),els_py->at(i),els_pz->at(i),els_energy->at(i));
	 valid_eles.push_back(eles);
       }//end of Nels loop 
       
       if(nGoodEle>0){
	 for(size_t i=0;i<valid_eles.size();++i){
	   double mass_mee = ( electrons.at(0)+valid_eles.at(i) ).M();
	   valid_fillHisto(valid_mass_ee, Boolcuts, mass_mee);
	 }
       }

       valid_fillHisto(valid_metEt,  Boolcuts,  this_met );
       valid_fillHisto(valid_metPhi, Boolcuts,  this_met_phi );

       valid_fillHisto(valid_numberTracks, Boolcuts, Ntracks);
       for(unsigned int i=0; i<Ntracks; ++i){
	 valid_fillHisto(valid_trackPt, Boolcuts, tracks_pt->at(i));
       }
       
       if(m_nGoodJet>2){
	 pair<double,double> res = compute_M3(jets);
	 valid_fillHisto(valid_recoM3,       Boolcuts, res.first );
	 valid_fillHisto(valid_recoM3_PTMax, Boolcuts, res.second );
       }
       
       if(!IsData() && isTTbar){ //check some gen. quantities for top quark
	 double genttbarpt = 0;
	 double genttbarpx = 0;
	 double genttbarpy = 0;
	 for(unsigned int g=0; g<Nmc_doc; g++){
	   if(fabs(mc_doc_id->at(g))==6){
	     valid_fillHisto(valid_genT_pt, Boolcuts, mc_doc_pt->at(g) );
	     genttbarpx+=mc_doc_px->at(g);
	     genttbarpy+=mc_doc_py->at(g);
	   }
	 }
	 genttbarpt = sqrt(genttbarpx*genttbarpx + genttbarpy*genttbarpy);
	 valid_fillHisto(valid_genTT_pt, Boolcuts, genttbarpt);
       }
     }//end m_doValidation
     //81FB end
     //------------------------------- End validation ---------------------------   




     //*******************************************************************************************
     //*******************************************************************************************
     // Now select events based on quantities filled above
     //*******************************************************************************************
     //*******************************************************************************************
	
	
     //Event counts
     //stage  0 : initial count from files
     //stage  1 : good run checked
     //stage  2 : trigger checked
     //stage  3 : >= 1 lepton passing cuts w/o isolation
     //stage  4 : >= 1 lepton passing cuts incl isolation
     //stage  5 : == 1 lepton passing cuts incl isolation
     //stage  6 : pass muon veto
     //stage 6b : at least 4 high-pt jets
     //stage  7 : pass MET cut
     //stage  8 : pass Z veto
     //stage  9 : not a conversion
     //stage 10 : letpton close to PV / barrel ele
     //stage 12 : >= 1+ tagged jet
     //stage 13 : >= 2+ tagged jet
	
     //Fill counters	
     bool e_plus_jet_pass = false;
      
     if(m_debug) { 
       cout << "this event:  weight: " << this_weight ;
       //cout << "wmap[this_mc]: "<< weightMap[this_mc] 
       cout << "  mc: " << this_mc;
       cout << "  type: " << mctype;
       cout << endl;
     }

     // INITIAL BEFORE ANY CUT
     //e_plus_jet[0][ntj][mctype]++;
     //e_plus_jet_weighted[0][ntj][mctype] += this_weight; //1 event x weight
     FillEventCounter(0, ntj, mctype); //TEST

     if(m_debug) cout << " Applying event selection" << endl;

     if(goodrun) { //Event in GOOD run list
	  
       //e_plus_jet[1][ntj][mctype]++;
       //e_plus_jet_weighted[1][ntj][mctype] += this_weight;
       FillEventCounter(1, ntj, mctype); //TEST

       if(fired_single_em) {  //Trigger
	 //e_plus_jet[2][ntj][mctype]++;	  
	 //e_plus_jet_weighted[2][ntj][mctype] += this_weight;	
	 FillEventCounter(2, ntj, mctype); //TEST

	 // Electron checks
	 if(nGoodEle > 0){
	   FillEventCounter(3, ntj, mctype); //TEST
	   
	   // Isolated electron
	   if(nGoodIsoEle > 0){
	     FillEventCounter(4, ntj, mctype); //TEST
	     
	     if(nGoodIsoEle == 1){
	       FillEventCounter(5, ntj, mctype); //TEST	     
  
	       if(!isMuon){ // Muon Veto 
		 FillEventCounter(6, ntj, mctype); //TEST

		 if( this_met > METCUT ){  // MET
		   FillEventCounter(7, ntj, mctype); //TEST

		   if(!isZ){ // Z Veto
		     FillEventCounter(8, ntj, mctype); //TEST

		     if(!isConversion){  //Conversion Veto
		       FillEventCounter(9, ntj, mctype); //TEST

		       if( ( m_rejectEndcapEle==false && !isDifferentInteraction ) ||  // PV check (DIFFZ)
			   ( m_rejectEndcapEle==true  && isBarrel ) ) {     // ele eta cut

			 FillEventCounter(10, ntj, mctype); //TEST		 
			 e_plus_jet_pass = true;

			 if(m_nbtag_SSV >= 1){ //at least one +tag
			   FillEventCounter(11, ntj, mctype); //TEST

			   if(m_nbtag_SSV >= 2){ //at least two +tag			     
			     FillEventCounter(12, ntj, mctype); //TEST
			   }
			   /*
			   if(ntaggable > 0){ // taggable
			     e_plus_jet[12][ntj][mctype]++;
			     e_plus_jet_weighted[12][ntj][mctype] += this_weight;
			     
			     if(nbtagP >= 1){ //at least one +tag
			       e_plus_jet[13][ntj][mctype]++;
			       e_plus_jet_weighted[13][ntj][mctype] += this_weight;
			       
			       if(nbtagP >= 2){ // at least two +tags
			       e_plus_jet[14][ntj][mctype]++;
			       e_plus_jet_weighted[14][ntj][mctype] += this_weight;
			       }
			     }			   
			     if(nbtagN >= 1){ //at least one -tag
			       e_plus_jet[15][ntj][mctype]++;
			       e_plus_jet_weighted[15][ntj][mctype] += this_weight;
			     }
			   }
			   */
			 }
		       }
		     }
		   }
		 }
	       }
	     }
	   }
	 }
       }
     }



     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     //+++++++++++++++++++++++++  Analysis!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     
     if(m_debug) cout << " Starting << ANALYSIS >>" << endl;

     if(goodrun && fired_single_em && nGoodEle > 0){
       fillHistoDataAndMC( h_met_ante_ISO,    this_met       ); //user-chosen MET
       fillHistoDataAndMC( h_met_ante_ISO_mu, met2v_mu.Mod() );
       fillHistoDataAndMC( h_met_ante_ISO_t1, met2v_t1.Mod() );
     }

     //E+jets Analysis
     if(e_plus_jet_pass && m_nGoodJet >= 4) {
	  
       ++counter_pass;

       //Sample hist
       // pass_ht->Fill(ht, this_weight);
       // pass_met->Fill(this_met, this_weight);
       fillHistoDataAndMC( h_ed0_pass, this_isoele_d0 );

       //Print out for each selected event
       fprintf(outfile,"%7d %6d %10d %10d %8d %8lld %5dj %11.2f %8.2f %8.2f %8s   %s\n",
	       counter_pass, run, event, lumiBlock, chain->GetTreeNumber(), 
	       lflag, m_nGoodJet, jets.at(0).Pt(), this_met, ht, this_mc.c_str(),
	       chain->GetCurrentFile()->GetName());
       interestingFiles.insert( chain->GetCurrentFile()->GetName() );
     }


     // (21 Feb 09) make some kinematics plots for events passing N-1 cuts (HT,MET)
     if(m_debug) cout << "[DEBUG] Filling N-1 histograms"<< endl;
     if( goodrun  &&  fired_single_em  &&  nGoodIsoEle==1  &&
     	 !isMuon  &&  !isZ  &&  !isConversion  &&  !isDifferentInteraction ) {

       if ( this_met > METCUT )    //after all but HT cut
	 fillHisto_Njet_DataAndMC( h_HT, ht, this_weight );
	 
       if ( ht >= HTCUT ) {        //after all but MET cut
	 fillHisto_Njet_DataAndMC( h_met,    this_met,       this_weight ); //user-chosen MET
	 fillHisto_Njet_DataAndMC( h_met_mu, met2v_mu.Mod(), this_weight );
	 fillHisto_Njet_DataAndMC( h_met_t1, met2v_t1.Mod(), this_weight );
	 if(!IsData()){//MC
	   fillHisto_Njet_DataAndMC( h_met_gen,                        this_genmet,   this_weight );//NEW 3-3-10
	   fillHisto_Njet_DataAndMC( h_met_gen_diff_mu, met2v_mu.Mod()-this_genmet,   this_weight );//NEW 3-3-10
	   fillHisto_Njet_DataAndMC( h_met_gen_diff_t1, met2v_t1.Mod()-this_genmet,   this_weight );//NEW 3-3-10
	   fillHisto_Njet_DataAndMC( h_met_gen_dphi_mu, met2v_mu.DeltaPhi(met2v_gen), this_weight );//NEW 3-3-10
	   fillHisto_Njet_DataAndMC( h_met_gen_dphi_t1, met2v_t1.DeltaPhi(met2v_gen), this_weight );//NEW 3-3-10
	 }
	 if(isBarrel){ //26-1-10
	   fillHisto_Njet_DataAndMC( h_met_BA,    this_met,       this_weight );
	   fillHisto_Njet_DataAndMC( h_met_mu_BA, met2v_mu.Mod(), this_weight );
	   fillHisto_Njet_DataAndMC( h_met_t1_BA, met2v_t1.Mod(), this_weight );
	   if(!IsData()){//MC
	     fillHisto_Njet_DataAndMC( h_met_gen_BA,                        this_genmet,   this_weight );//NEW 3-3-10
	     fillHisto_Njet_DataAndMC( h_met_gen_diff_mu_BA, met2v_mu.Mod()-this_genmet,   this_weight );//NEW 3-3-10
	     fillHisto_Njet_DataAndMC( h_met_gen_diff_t1_BA, met2v_t1.Mod()-this_genmet,   this_weight );//NEW 3-3-10
	     fillHisto_Njet_DataAndMC( h_met_gen_dphi_mu_BA, met2v_mu.DeltaPhi(met2v_gen), this_weight );//NEW 3-3-10
	     fillHisto_Njet_DataAndMC( h_met_gen_dphi_t1_BA, met2v_t1.DeltaPhi(met2v_gen), this_weight );//NEW 3-3-10
	   }
	 }
       }
       // after all but MET/HT/Njet cuts
       fillHisto_Njet_DataAndMC( h_mtw_mu, this_mu_mtw, this_weight );
       fillHisto_Njet_DataAndMC( h_mtw_t1, this_t1_mtw, this_weight );       
       fillHisto_Njet_DataAndMC( h_DPhiEmet_mu, this_mu_DPhiEmet, this_weight );
       fillHisto_Njet_DataAndMC( h_DPhiEmet_t1, this_t1_DPhiEmet, this_weight );

       // inspect distributions of the selected events
       if( m_nGoodJet >= 4 ) {
	 fillHistoDataAndMC( h_exp_ele_et,  iso_electrons.at(0).Et()  );
	 fillHistoDataAndMC( h_exp_ele_eta, iso_electrons.at(0).Eta() );
	 fillHistoDataAndMC( h_exp_j0_pt,   jets.at(0).Pt()           );
	 fillHistoDataAndMC( h_exp_j1_pt,   jets.at(1).Pt()           );
	 fillHistoDataAndMC( h_exp_DRej,    iso_electrons.at(0).DeltaR(   jets.at(0) ) );
	 fillHistoDataAndMC( h_exp_DPhiej,  iso_electrons.at(0).DeltaPhi( jets.at(0) ) );
	 fillHistoDataAndMC( h_exp_DRjj,    jets.at(0).DeltaR(   jets.at(1) ) );
	 fillHistoDataAndMC( h_exp_DPhijj,  jets.at(0).DeltaPhi( jets.at(1) ) );
	 fillHistoDataAndMC( h_exp_met_v_eeta,  this_met, iso_electrons.at(0).Eta(), this_weight );
       }
     }
     if(m_debug) cout << "[DEBUG] After N-1"<< endl;




     //----------------------------------------------------  
     //  Make isolation plots needed by QCD Estimation
     // TL: 1D method: plot Reliso distribution for events passing all cuts except electron isolation and nj
     // NOTE: at the moment still using the old definition of RelIso = et/(et+SumIso), to be changed.
     // notes: all cuts except ISO and Njet
     //cout << "200) isZ = " << isZ << endl;

     if(m_debug) cout << " Starting << Making isolation plot >>" << endl;

     float CombRelIso = -1.0;
     //float NormCombRelIso = -1.0;
     int  ii_GoodEle_mostIso = -1;
   
     //cout << "(B4) ii_GoodEle_mostIso = " << ii_GoodEle_mostIso << endl;

     if(nGoodEle>0) {
       if (nGoodEle==1) {
	 // need isolation values of the GoodEle but not available!!
	 // can resolve this by creating a companion vector to store the CombRelIso values of the GoodEles
	 //  <electrons> ++ <electrons_isoval>
	 CombRelIso = electrons_isoval.at(0);
	 //NormCombRelIso = electrons_isoval2.at(0);
	 ii_GoodEle_mostIso = ii_electrons.at(0); 
       } else {
	 // pick the most isolated one, loop over GoodEle, because if we happen to have 2 GoodEle,
	 // and one pass isolation cut and the other fails, then if we choose randomly we may choose the 
	 // non-isolated one, and contradict our event selection.
	 // bug fix (17 Feb 09)
	 float mostIso = 9999999.9; //most-isolated = smallest Reliso (not normalized)
	 float mostIso2 = -9999999.9; //most-isolated = largest Reliso (normalized)
	 for (int i=0; i<nGoodEle; ++i) {
	   float thisIso = electrons_isoval.at(i);
	   float thisIso2 = electrons_isoval2.at(i);
	   if(thisIso < mostIso) { //RelIso
	     mostIso = thisIso; 
	     ii_GoodEle_mostIso = ii_electrons.at(i); 
	   }
	   if(thisIso2 > mostIso2) { mostIso2 = thisIso2; } //NormReliso
	 }
	 CombRelIso = mostIso;
	 //NormCombRelIso = mostIso2;
       }
       //cout << "ii_GoodEle_mostIso = " << ii_GoodEle_mostIso << endl;
       
       // To identify the index of the most isolated GoodEle in the Nels collection
//        
//        for(unsigned int i=0; i<Nels; ++i){
// 	 float tempComIso = (els_tIso->at(i) + els_dr04EcalRecHitSumEt->at(i) + els_dr04HcalTowerSumEt->at(i))/els_et->at(i);
// 	 if( (tempComIso - CombRelIso)/CombRelIso  < 1e-5 ){ ii_GoodEle_mostIso = i; break;}
// 	 //changed 280909 to make relative comparison due to diff of 1.01e-6 for one event.e-5 = 0.001% 
//        }

       if(ii_GoodEle_mostIso < 0){ //bugfix: 9-9-09
	 cout <<"\nERROR: PROBLEM, COULDNT FIND MOST_ISO_ELECTRON AGAIN"<<endl;
	 cout << "event number is" << ev << endl;
	 cout << " index of most isolated GoodEle (ii_GoodEle_mostIso): " << ii_GoodEle_mostIso << endl;
	 cout << " Nels: "<< Nels << ",  nGoodEle: "<< nGoodEle << endl;
	 printf(" smallest isolation: %12.10f\n" ,CombRelIso);
	 cout << " Isolation of all electrons: " << endl;	 
	 for(unsigned int j=0; j<Nels; ++j){
	   //float tmpIso = (els_tIso->at(j) + els_dr04EcalRecHitSumEt->at(j) + els_dr04HcalTowerSumEt->at(j))/els_et->at(j);
	   float tmpIso = getRelIso(j);
	   printf("ele %d : %12.10f\n",j, tmpIso);
	 }
       }
	 
     }//nGoodEle > 0
     
     //if electron is not isolated, apply conversion algo on it so that all electrons in the distribution 
     //are treated the same wrt conversion algo
     //987
     
     if(m_debug) cout << " Apply << Conversion algo >> on all 'most-iso' good electrons" << endl;
     
     if(CombRelIso>0.1 && nGoodEle>0 ){
       //cout <<"ii_GoodEle_mostIso = "<< ii_GoodEle_mostIso  << endl;
       TLorentzVector eles_temp(els_px->at(ii_GoodEle_mostIso),els_py->at(ii_GoodEle_mostIso),els_pz->at(ii_GoodEle_mostIso),els_energy->at(ii_GoodEle_mostIso));
       isConversion = ConversionFinder(eles_temp, mctype, ii_GoodEle_mostIso);
     }
   
   
     //-------------------------------------------------------------------
     //
     //  QCD estimation: Normal & Anti Event Selection
     //
     //-------------------------------------------------------------------
     if(m_debug) cout << " Starting << QCD estimation NES AES >>" << endl;

     /// bool isZ_mee_AES = false;
     /// bool isZ_mep_AES = false;
     float mass_ep = -1; //NEW
   
     if( goodrun  &&  fired_single_em  &&  nGoodEle>0  && 
	 !isMuon  &&  !isZ  &&  !isConversion  &&  !isDifferentInteraction  &&  ht >= HTCUT ) { // <-- HT cut


       // Apply missing ET cut (Normal Selection)
       bool passALL = true;
       if ( m_applyMETcut && this_met < METCUT ) passALL = false;

       // Apply eta cut if specified (optional)
       if ( m_rejectEndcapEle && fabs( els_eta->at(ii_GoodEle_mostIso) ) > 1.442 ) passALL = false;
     
       //if( (RunPlanB()==false && this_met > METCUT)  || //plan A
       //    (RunPlanB()==true  && fabs( els_eta->at(ii_GoodEle_mostIso) ) < 1.442) ) { //plan B

       if (passALL) {
	 if(m_debug) cout << "QCDest: event pass all cuts (except isol and nj), make isolation plots" << endl;
	 //fill histo (with weight) accord to m_nGoodJet (16-8-09)
	 fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso,      CombRelIso,     this_weight );
	 //fillHisto_Njet_DataAndMC( h_QCDest_NormCombRelIso,  NormCombRelIso, this_weight );
       }
   


       
       //=====================================
       // Tighter Z veto in Anti Selection
       //================= 8-6-09 ==========
       // move the following block to outside MET cut
       //=====================================
       // - M(e,e): require 1 GoodEle (no isolation requirement), and 1 "robustloose" electron
       // - M(e,e) veto window: 65-106 GeV
       // - M(e,pho): take the 1st GoodEle, and the 1st photon (no cut)
       // - M(e,pho) veto window: 40-110 GeV
       //=====================================


       //------------------------
       // 8-6-09: Define flags for AES
       bool flag_AES_pass_metcut =  this_met < AES_MET_cut;
       bool flag_AES_pass_HTcut  = false;
       bool flag_AES_pass_tighterZveto_mee = true;
       bool flag_AES_pass_tighterZveto_mep = true;
       //------------------------
       

     
       //1510  Simple Z veto: exclude event with 2 Nels from AES
       //       if(AES_useSimpleZveto_TwoRL) flag_AES_pass_tighterZveto_mee = !isZ_twoRL;
       if(AES_useSimpleZveto) {
	 if ( Nels >=2 ) flag_AES_pass_tighterZveto_mee = false; //exclude
       }


       //======================================
       //  M(e,e): di-electron inv mass [AES]
       //======================================

       if ( AES_useSimpleZveto==false  &&  nGoodEle > 0  &&  Nels > 1 ) {

	 flag_AES_pass_tighterZveto_mee = false;

	 //cout << "TEST: N(e) = " << Nels << endl;
	 for (unsigned int i=0; i<Nels; ++i){ //1st e loop
	   
	   //cout << "\n TEST: e[" << i << "]"  ;
	   // Require 1 GoodEle
	   if( els_et->at(i) < 30.0 ) continue; 
	   if( fabs( els_eta->at(i) ) > 2.5 ) continue;
	   if( fabs( els_eta->at(i) ) > 1.442 &&
	       fabs( els_eta->at(i) ) < 1.56 ) continue; //ignore gap
	   if( fabs(compute_d0("electron",i)) > 0.02 ) continue;
	   //if( els_robustTightId->at(i) < 1 ) continue;
	   if ( passEleID(i) ==false )  continue;

	   for (unsigned int j=0; j<Nels; ++j){ //2nd e loop
	     if(j==i) continue;
	   
	     // Require "loose" electron: consider only electrons with ET>20 GeV, eta<2.5, d0<200um, RobustLoose
	     if( els_et->at(j) < 20.0 ) continue;
	     if( fabs( els_eta->at(j) ) > 2.5 ) continue;
	     if( fabs(compute_d0("electron",j)) > 0.02 ) continue; //d0 cut
	     if( els_robustLooseId->at(j) < 1 ) continue;
	     
	     TLorentzVector tight(els_px->at(i),els_py->at(i),els_pz->at(i),els_energy->at(i));
	     TLorentzVector loose(els_px->at(j),els_py->at(j),els_pz->at(j),els_energy->at(j));

	     float mass = ( tight + loose ).M();
	     //cout << "   m(e,e) = " << mass << endl;
	     // if fall in the window, flag as Z
	     if ( mass < 65. || mass > 106. )  flag_AES_pass_tighterZveto_mee = true; // not Z, keep event
	   
	   }//2nd e loop
	 }//1st e loop
	 //cout << "isZ_mee_TEST = " << isZ_mee_TEST << endl;
       }//if more than 2 electron
       

       
       //=============================================
       //  M(e,pho): ele-photon invariant mass [AES]
       //=============================================
       // revise this...
       // at the moment, it only consider the inv mass of the first goodEle and the first photon

       if( nGoodEle>=1 && Nphotons>=1 ){

	 flag_AES_pass_tighterZveto_mep = false;
	 TLorentzVector photon1( photons_px->at(0),photons_py->at(0),photons_pz->at(0),photons_energy->at(0) );
	 mass_ep = (electrons.at(0) + photon1).M(); //GoodEle + photon
	 if(mass_ep < 40. || mass_ep > 110.) flag_AES_pass_tighterZveto_mep = true; //outside window, keep event
	 
       }//require 1 GoodEle, 1 photon

       // tighter Z veto for AES
       //// if(!isZ_mee_AES) flag_AES_pass_tighterZveto_mee = true;
       //// if(!isZ_mep_AES) flag_AES_pass_tighterZveto_mep = true;
       //================= 8-6-09 ==========
     
     
       // Fill histograms for Z veto (AES), lowMET, 1j (already pass the normal Z veto)
       if(m_studyZveto) {
	 if( flag_AES_pass_metcut && m_nGoodJet==1 ) {

	   if(mass_ee>0) fillHistoDataAndMC( h_mass_diele_lowMet_1j, mass_ee ); //plot m(ee) before cut
	   if( flag_AES_pass_tighterZveto_mee ) { //event survive mee cut
	     if(mass_ep>0) fillHistoDataAndMC( h_mass_ephoton_lowMet_1j, mass_ep ); //plot m(ep) before cut
	   }
	   // record number of ele/photon in this event
	   fillHistoDataAndMC( h_Nele_lowMet_1j,    Nels     );
	   fillHistoDataAndMC( h_Nphoton_lowMet_1j, Nphotons );
	   for(unsigned int i=0; i<Nphotons; ++i){
	     fillHistoDataAndMC( h_photon_eta_lowMet_1j, photons_eta->at(i) );
	     fillHistoDataAndMC( h_photon_et_lowMet_1j,  photons_et->at(i)  );
	     if(i==0) {
	       fillHistoDataAndMC( h_photon1_eta_lowMet_1j, photons_eta->at(i) );
	       fillHistoDataAndMC( h_photon1_et_lowMet_1j,  photons_et->at(i)  );
	     }
	   }
	 }
       }//m_studyZveto
       
     
       //===============================
       // AES:  HT cut (HT<100)
       //===============================

       int ii_AES = ii_GoodEle_mostIso ; // 15-6-09, move to front

       float ht_AES = this_met;
       ht_AES += els_et->at(ii_AES);
       for (unsigned int i=0; i<jets.size(); ++i) {
	 ht_AES += jets.at(i).Pt();
       }
       
       if(ht_AES < AES_HT_cut) flag_AES_pass_HTcut = true;
       //===============================




       //987
       //----------------------------------------------
       // 8-6-09: fill histo for n-1 AES reliso plots
       //----------------------------------------------
       // Before AES cut
       fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_before, CombRelIso, this_weight );

       
       // AES: all cut
       if( flag_AES_pass_metcut && 
	   flag_AES_pass_HTcut && 
	   flag_AES_pass_tighterZveto_mee && 
	   flag_AES_pass_tighterZveto_mep ) {

	 fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES, CombRelIso, this_weight ); //All AES cut
	 //cout << "AES pass all" << endl;
       }
              

       // AES: N-1 (met)
       if( flag_AES_pass_HTcut ){

	 fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_justHT, CombRelIso, this_weight );

	 if( flag_AES_pass_tighterZveto_mee && flag_AES_pass_tighterZveto_mep ){
	   //cout << "amtb AES - met" << endl;
	   fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_minusMET, CombRelIso, this_weight );
	 }
       }
       // AES: N-1 (HT)
       if( flag_AES_pass_metcut ) {

	 fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_justMET, CombRelIso, this_weight );

	 if (flag_AES_pass_tighterZveto_mee && flag_AES_pass_tighterZveto_mep ){
	   //cout << "amtb AES - ht" << endl;
	   fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_minusHT, CombRelIso, this_weight );
	 }
	 
	 // AES: N-1 (tighter Z, both mee, mep)
	 if( flag_AES_pass_HTcut ){
	   //cout << "amtb AES - Z" << endl;
	   fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_minusTighterZ, CombRelIso, this_weight );
	 }
       }//pass AES_MET
       else { //event has high MET
	 fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_justHighMET, CombRelIso, this_weight );
       }
       // AES: just tighter Z
       if( flag_AES_pass_tighterZveto_mee && flag_AES_pass_tighterZveto_mep ) 
	 fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_justZ, CombRelIso, this_weight );

     }//pass some cuts
     // TL: end ------------------------------


     //=====================================================================
     //
     //  
     //    QCD Control Region (AES)        plan A      4-12-09
     //
     //--------------------------------------------------------
     // Trial AES definition for Plan A: MET>20.
     //
     // Common: met<15, HT<150,  exactly 1 GoodEle (RT).
     //
     // Definitions: A1 : leave out !conv cut
     //              A2 : invert !conv cut, ie plot only conversion ele.
     //              A3 : invert d0 cut = d0 > 200um.
     //--------------------------------------------------------
     // First make sure selected events do not enter control sample
     //     if(m_debug) cout << << endl;
     if( e_plus_jet_pass == false ) {

       if( goodrun  &&  fired_single_em  &&  nGoodEle==1  &&
	   !isMuon  &&  !isZ  &&  !isDifferentInteraction  ) {
 
	 // the Good Electron
	 unsigned int theGE = ii_electrons.at(0);  //first good ele

	 float ht_AES = this_met + els_et->at( theGE );
	 for (unsigned int i=0; i<jets.size(); ++i) {
	   ht_AES += jets.at(i).Pt();
	 }

	 if ( this_met < AES_MET_cut  &&  ht_AES < AES_HT_cut ) {
	   // satisfy criteria for plan A control sample
                                                      
	   const double this_et  = els_et->at( theGE );
	   const double this_iso = getRelIso( theGE );


	   // A1: no conv cuts
	   //cout << "[DEBUG] A1" << endl;
	   if ( this_et > 20. ) {  fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planA1_e20, this_iso, this_weight );
	     if ( this_et > 30. )  fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planA1_e30, this_iso, this_weight );	
	   }

	   // A2: conv only
	   //if(m_debug) cout << "   calling ConversionFinder()" << endl;
	   //cout << "[DEBUG] A2" << endl;
	   //	   TLorentzVector eles_temp(els_px->at(0),els_py->at(0),els_pz->at(0),els_energy->at(0));
	   bool this_is_conv = ConversionFinder( electrons.at(0), mctype, 0); //1st good ele

	   if ( this_is_conv ) {
	     if ( this_et > 20. ) { fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planA2_e20, this_iso, this_weight );
	       if ( this_et > 30. ) fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planA2_e30, this_iso, this_weight );
	     }
	   }

	   // A3: invert d0 (>200 um)
	   //cout << "[DEBUG] A3" << endl;
	   if ( fabs(compute_d0("electron",theGE)) > 0.02 ) {
	     if ( this_et > 20. ) { fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planA3_e20, this_iso, this_weight );
	       if ( this_et > 30. ) fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planA3_e30, this_iso, this_weight );
	     }
	   }
	 } 
       }
     }
     // END: plan A control sample
     //=====================================================================
     if(m_debug) cout << "End of planA control sample" << endl;



     //=====================================================================
     //
     //    QCD Control Region (AES)        plan B      10-11-09
     //
     //----------------------------------------------------
     // Trial AES definition for Plan B: noMET, eta(e)<1.5
     //
     // Definitions: B1: - 1 electron with ET>20, |eta|<1.442, d0<200um, RT, not conversion.
     //                  - 1 jet
     //                  - simple Z veto: no 2 reco ele     
     //                  - std cuts: HLT, !mu
     //
     // Definitions: B2: - 1 electron with ET<30, |eta|<1.442, d0<200um, RT, not conversion.
     //                  - 1 jet
     //                  - simple Z veto: no 2 reco ele     
     //                  - std cuts: HLT, !mu
     // Definitions: B3: - 1 electron withh ET>20/30, |eta|<1.442, d0<200um, fail RT, not conversion.
     //                  - 1 jet
     //                  - simple Z veto: no 2 reco ele
     //                  - std cuts: HTL, !mu
     // Definitions: B4: same as B3 but failing RL ID.
     // Definitions: B5: same as B3 but failing category-based Loose ID.
     // Definitions: B6: same as B3 but failing category-based Tight ID.
     // Definitions: B7: same as B1 but inverting d0 cut, req d0>200um..
     // Definitions: B8: same as B7 but fail RT.
     //-----------------------------------------------------
     // First make sure selected events do not enter control sample
     if( e_plus_jet_pass == false ) {
       if( goodrun  &&  fired_single_em  &&    Nels==1   &&  !isMuon  &&  !isDifferentInteraction ) {//AAAAA

	 if(m_debug) cout << "  passing selections" << endl;
	 // common criteria for B1 & 2
	 //if(  fabs(els_eta->at(0))           < 1.442  &&
	 if(  fabs(els_eta->at(0))           < 2.5  &&
	      fabs(compute_d0("electron",0)) < 0.02     ) {
	   

	   if(m_debug) cout << "  passing ele eta and d0 cuts" << endl;
	   //TLorentzVector eles_temp(els_px->at(0),els_py->at(0),els_pz->at(0),els_energy->at(0));
	   //isConversion = ConversionFinder(eles_temp, mctype, 0);
	   
	   //if ( !isConversion ) {
	   if ( true ) {
	     
	     const double this_iso = getRelIso(0);
	     
	     if ( els_robustTightId->at(0)  > 0 ) { //pass RT ID
	       
	       // Definition B1:  lowered EleET cut
	       //------------------------------------
	       if(m_debug) cout << "    B1" << endl;
	       if ( els_et->at(0) > 20. )  { fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB1_e20, this_iso, this_weight );
		 if ( els_et->at(0) > 30. )  fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB1_e30, this_iso, this_weight );
	       }
	       // Definition B2:  inverted EleET cut
	       //------------------------------------
	       if(m_debug) cout << "    B2" << endl;
	       if( els_et->at(0) < 30. ) {  fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB2_e30, this_iso, this_weight );
		 if( els_et->at(0) < 20. )  fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB2_e20, this_iso, this_weight );
	       }
	       
	     } else {
	       // Definition B3: fail RT ID
	       //---------------------------
	       if(m_debug) cout << "    B3" << endl;
	       if ( els_et->at(0) > 20. ) { fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB3_e20, this_iso, this_weight );
		 if ( els_et->at(0) > 30. )  fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB3_e30, this_iso, this_weight );
	       } 
	       if( fabs(els_eta->at(0))  < 1.442 ){ //BARREL
		 if ( els_et->at(0) > 20. ) { fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB3b_e20, this_iso, this_weight );
		   if ( els_et->at(0) > 30. )  fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB3b_e30, this_iso, this_weight );
		 }
	       }
	     }
	     // 2-12-09
	     if(m_debug) cout << "    B4,5,6" << endl;
	     if( els_et->at(0) > 20. ) {
	       if( els_robustLooseId->at(0) == 0 ) fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB4_e20, this_iso, this_weight );
	       if( els_looseId->at(0) == 0 )       fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB5_e20, this_iso, this_weight );
	       if( els_tightId->at(0) == 0 )       fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB6_e20, this_iso, this_weight );
	       
	       if( els_et->at(0) > 30. ) {
		 if( els_robustLooseId->at(0) == 0 ) fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB4_e30, this_iso, this_weight );
		 if( els_looseId->at(0) == 0 )       fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB5_e30, this_iso, this_weight );
		 if( els_tightId->at(0) == 0 )       fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB6_e30, this_iso, this_weight );
	       }
	     }
	     
	   }
	 }
	 // B7: invert d0 cut
	 if(m_debug) cout << "    B7" << endl;
	 if(  els_et->at(0) > 20.  && 
	      fabs(els_eta->at(0)) < 2.5  &&
	      fabs(compute_d0("electron",0)) > 0.02 ) { //<---
	   
	   const double this_iso = getRelIso(0);
	   
	   if( els_robustTightId->at(0)  > 0 ) { //pass RT ID
	     fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB7_e20, this_iso, this_weight );
	     if( els_et->at(0) > 30. )   fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB7_e30, this_iso, this_weight );
	   }
	   else { //fail RT ID
	     if(m_debug) cout << "    B8" << endl;
	     fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB8_e20, this_iso, this_weight );
	     if( els_et->at(0) > 30. )   fillHisto_Njet_DataAndMC( h_QCDest_CombRelIso_AES_planB8_e30, this_iso, this_weight );
	   }//RT ID
	 }
       }//pass selections
     }// plan B control sample
     //=====================================================================
     if(m_debug) cout << "End of planB control sample" << endl;



     //---------------------------------------------------------------
     // 9 Jun 09: inspect isolation/met at each key stage of cuts
     //           to study correlation of iso:met for QCD events
     // level 1: HLT
     //       2: nGoodEle>0
     //       3: mu, Z veto
     //       4: conversion cut
     //---------------------------------------------------------------       
     //cout << "* CombRelIso = " << CombRelIso << endl;

     if( m_plotRelisoNES ) {

       if( goodrun  &&  fired_single_em ) {


	 string metside = "loMET"; //low met
	 if( this_met > METCUT ) metside = "hiMET"; //high met
	 //cout << "this_met ="<< this_met << endl;
	 //cout << "metside = " << metside << endl;


	 // Pass L1
	 if(m_debug) cout << "-> Filling Reliso NES histograms, L1 HLT" << endl; 
	 // fill for all electrons
	 for(unsigned int ie=0; ie<Nels; ie++){

	   //double tmpRelIso = (els_tIso->at(ie) + els_dr04EcalRecHitSumEt->at(ie) + els_dr04HcalTowerSumEt->at(ie)) / els_et->at(ie);
	   float tmpRelIso = getRelIso(ie);

	   // barrel or endcap?
	   string etaside = "barrel";
	   bool inBarrel = true;
	   if( fabs(els_eta->at(ie)) > 1.5 ) { etaside = "endcap"; inBarrel=false;}

	   iso_fillHisto_NES( 0, inBarrel, tmpRelIso, this_met );
// 	   iso_fillHisto_nlevel_nj_nmc( h_QCDest_isoVmet_NES[0],    tmpRelIso, this_met, this_weight ); //NEW 0: L1
// 	   iso_fillHisto_nlevel_nj_nmc( h_QCDest_isoVmet_NES_uw[0], tmpRelIso, this_met, 1           ); //NB: unweighted
// 	   iso_fillHisto_nlevel_nj_nmc( h_QCDest_CombRelIso_NES[0], tmpRelIso, this_met              );
// 	   if(inBarrel) {
// 	     iso_fillHisto_nlevel_nj_nmc( h_QCDest_isoVmet_NES_barrel[0],    tmpRelIso, this_met );
// 	     iso_fillHisto_nlevel_nj_nmc( h_QCDest_CombRelIso_NES_barrel[0], tmpRelIso, this_met );
// 	   } else {
// 	     iso_fillHisto_nlevel_nj_nmc( h_QCDest_isoVmet_NES_endcap[0],    tmpRelIso, this_met );
// 	     iso_fillHisto_nlevel_nj_nmc( h_QCDest_CombRelIso_NES_endcap[0], tmpRelIso, this_met );
// 	   }

// 	   fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_L1",                tmpRelIso, this_met, this_weight );
// 	   fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_"+etaside+"_L1",    tmpRelIso, this_met, this_weight );
// 	   fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_uw_L1",             tmpRelIso, this_met, 1 );//unweighted
// 	   fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_uw_"+etaside+"_L1", tmpRelIso, this_met, 1 );//unweighted
// 	   fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_L1",             tmpRelIso, this_weight ); //no met cut
// 	   fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+metside+"_L1", tmpRelIso, this_weight );
// 	   fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+etaside+"_L1", tmpRelIso, this_weight );
// 	   fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+metside+"_"+etaside+"_L1", tmpRelIso, this_weight );



	   // fill only for electrons with at least 30 GeV ET, and pass eta cut
	   if( els_et->at(ie) > ELE_ETCUT && fabs(els_eta->at(ie)) < 2.5 &&
	       ( fabs(els_eta->at(ie)) < 1.442 || fabs(els_eta->at(ie)) > 1.56 ) ) {

	     iso_fillHisto_NES( 1, inBarrel, tmpRelIso, this_met );
	     /*
	     iso_fillHisto_nlevel_nj_nmc( h_QCDest_isoVmet_NES[1],    tmpRelIso, this_met, this_weight ); //NEW 1: L1b
	     iso_fillHisto_nlevel_nj_nmc( h_QCDest_isoVmet_NES_uw[1], tmpRelIso, this_met, 1           ); //NB: unweighted
	     iso_fillHisto_nlevel_nj_nmc( h_QCDest_CombRelIso_NES[1], tmpRelIso, this_met              );
	     if(isInBarrel) {
	       iso_fillHisto_nlevel_nj_nmc( h_QCDest_isoVmet_NES_barrel[1],    tmpRelIso, this_met );
	       iso_fillHisto_nlevel_nj_nmc( h_QCDest_CombRelIso_NES_barrel[1], tmpRelIso, this_met );
	     } else {
	       iso_fillHisto_nlevel_nj_nmc( h_QCDest_isoVmet_NES_endcap[1],    tmpRelIso, this_met );
	       iso_fillHisto_nlevel_nj_nmc( h_QCDest_CombRelIso_NES_endcap[1], tmpRelIso, this_met );
	     }
	     */

// 	     fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_L1b",                tmpRelIso, this_met, this_weight );
// 	     fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_"+etaside+"_L1b",    tmpRelIso, this_met, this_weight );
// 	     fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_uw_L1b",             tmpRelIso, this_met, 1 );//unweighted
// 	     fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_uw_"+etaside+"_L1b", tmpRelIso, this_met, 1 );//unweighted
// 	     fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_L1b",             tmpRelIso, this_weight ); //no met cut
// 	     fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+metside+"_L1b", tmpRelIso, this_weight );
// 	     fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+etaside+"_L1b", tmpRelIso, this_weight );
// 	     fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+metside+"_"+etaside+"_L1b", tmpRelIso, this_weight );
	   	   

	     // d0 cut
	     // Calculate d0 w.r.t beam spot
	     float d0_corrected  = compute_d0("electron",ie);

	     if( fabs(d0_corrected) < 0.02 ){
	     
	       iso_fillHisto_NES( 2, inBarrel, tmpRelIso, this_met );

// 	       fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_L1c",                tmpRelIso, this_met, this_weight );
// 	       fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_"+etaside+"_L1c",    tmpRelIso, this_met, this_weight );
// 	       fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_uw_L1c",             tmpRelIso, this_met, 1 );//unweighted
// 	       fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_uw_"+etaside+"_L1c", tmpRelIso, this_met, 1 );//unweighted
// 	       fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_L1c",             tmpRelIso, this_weight );
// 	       fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+metside+"_L1c", tmpRelIso, this_weight ); 
// 	       fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+etaside+"_L1c", tmpRelIso, this_weight );
// 	       fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+metside+"_"+etaside+"_L1c", tmpRelIso, this_weight ); 
	     
	       // eID (barrel)
	       //bool pass_eid_c0 =  els_robustTightId->at(ie) > 0; //out-of-box eID variable
	       bool pass_eid_c0 =  passEleID(ie); //out-of-box eID variable
	       bool pass_eid_c1 = false;
	       bool pass_eid_c2 = false;
	       bool pass_eid_c3 = false;
	       bool pass_eid_c4 = false;
	       bool pass_eid = false;

	       // electron ID variables
	       // 
	       // 31X: eID cuts have changed!!! values need to be updated!  TODO
	       //

	       if( fabs(els_eta->at(ie)) < 1.442 ) {  // barrel

                 // Updated eID cut values: tight Fixed Threshold "RobustTight"
		 // ref: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideCutBasedElectronID
		 if(     els_hadOverEm->at(ie) < 0.01   ) pass_eid_c1 = true;
		 if(  fabs(els_dEtaIn->at(ie)) < 0.0040 ) pass_eid_c2 = true;
		 if(  fabs(els_dPhiIn->at(ie)) < 0.025  ) pass_eid_c3 = true;	       
		 if( els_sigmaIEtaIEta->at(ie) < 0.0099 ) pass_eid_c4 = true;
	       
	       }
	       else if ( fabs(els_eta->at(ie)) > 1.56 ) { // endcap

		 if(     els_hadOverEm->at(ie) < 0.01   ) pass_eid_c1 = true;
		 if(  fabs(els_dEtaIn->at(ie)) < 0.0066 ) pass_eid_c2 = true;
		 if(  fabs(els_dPhiIn->at(ie)) < 0.020  ) pass_eid_c3 = true;
		 if( els_sigmaIEtaIEta->at(ie) < 0.028  ) pass_eid_c4 = true;

	       }
	     
	       pass_eid = pass_eid_c0 && pass_eid_c1 &&  pass_eid_c2 && pass_eid_c3 && pass_eid_c4; //all 5
	     
	       if(pass_eid_c1) {


		 iso_fillHisto_NES( 3, inBarrel, tmpRelIso, this_met );

// 		 iso_fillHisto_nlevel_nj_nmc( hh_QCDest_isoVmet_NES[3], tmpRelIso, this_met ) ; //NEW 3:L1d1
// 		 if(isInBarrel) iso_fillHisto_nlevel_nj_nmc( hh_QCDest_isoVmet_NES_barrel[3], tmpRelIso, this_met );
// 		 else           iso_fillHisto_nlevel_nj_nmc( hh_QCDest_isoVmet_NES_endcap[3], tmpRelIso, this_met );

// 		 fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_L1d1",                tmpRelIso, this_met, this_weight );
// 		 fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_"+etaside+"_L1d1",    tmpRelIso, this_met, this_weight );
// 		 fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_uw_L1d1",             tmpRelIso, this_met, 1 );
// 		 fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_uw_"+etaside+"_L1d1", tmpRelIso, this_met, 1 );
// 		 fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_L1d1",             tmpRelIso, this_weight );
// 		 fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+metside+"_L1d1", tmpRelIso, this_weight );
// 		 fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+etaside+"_L1d1", tmpRelIso, this_weight );
// 		 fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+metside+"_"+etaside+"_L1d1", tmpRelIso, this_weight );

		 if(pass_eid_c2) {
		   iso_fillHisto_NES( 4, inBarrel, tmpRelIso, this_met );

// 		   iso_fillHisto_nlevel_nj_nmc( hh_QCDest_isoVmet_NES[4], tmpRelIso, this_met ) ; //NEW 4: L1d2
// 		   if(isInBarrel) iso_fillHisto_nlevel_nj_nmc( hh_QCDest_isoVmet_NES_barrel[4], tmpRelIso, this_met );
// 		   else           iso_fillHisto_nlevel_nj_nmc( hh_QCDest_isoVmet_NES_endcap[4], tmpRelIso, this_met );

// 		   fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_L1d2",                tmpRelIso, this_met, this_weight );
// 		   fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_"+etaside+"_L1d2",   tmpRelIso, this_met, this_weight );
// 		   fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_uw_L1d2",             tmpRelIso, this_met, 1 );
// 		   fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_uw_"+etaside+"_L1d2", tmpRelIso, this_met, 1 );
// 		   fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_L1d2",             tmpRelIso, this_weight );
// 		   fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+metside+"_L1d2", tmpRelIso, this_weight );
// 		   fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+etaside+"_L1d2", tmpRelIso, this_weight );
// 		   fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+metside+"_"+etaside+"_L1d2", tmpRelIso, this_weight );

		   if(pass_eid_c3) {

		     iso_fillHisto_NES( 5, inBarrel, tmpRelIso, this_met );


// 		     iso_fillHisto_nlevel_nj_nmc( hh_QCDest_isoVmet_NES[5], tmpRelIso, this_met ) ; //NEW 5:L1d3
// 		     if(isInBarrel) iso_fillHisto_nlevel_nj_nmc( hh_QCDest_isoVmet_NES_barrel[5], tmpRelIso, this_met );
// 		     else           iso_fillHisto_nlevel_nj_nmc( hh_QCDest_isoVmet_NES_endcap[5], tmpRelIso, this_met );

// 		     fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_L1d3",                tmpRelIso, this_met, this_weight );
// 		     fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_"+etaside+"_L1d3",    tmpRelIso, this_met, this_weight );
// 		     fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_uw_L1d3",             tmpRelIso, this_met, 1 );
// 		     fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_uw_"+etaside+"_L1d3", tmpRelIso, this_met, 1 );
// 		     fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_L1d3",             tmpRelIso, this_weight );
// 		     fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+metside+"_L1d3", tmpRelIso, this_weight );
// 		     fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+etaside+"_L1d3", tmpRelIso, this_weight );
// 		     fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+metside+"_"+etaside+"_L1d3", tmpRelIso, this_weight ); 

		     if(pass_eid_c4) {

		       iso_fillHisto_NES( 3, inBarrel, tmpRelIso, this_met );


// 		       iso_fillHisto_nlevel_nj_nmc( hh_QCDest_isoVmet_NES[6], tmpRelIso, this_met ) ; //NEW 6:L1d4
// 		       if(isInBarrel) iso_fillHisto_nlevel_nj_nmc( hh_QCDest_isoVmet_NES_barrel[6], tmpRelIso, this_met );
// 		       else           iso_fillHisto_nlevel_nj_nmc( hh_QCDest_isoVmet_NES_endcap[6], tmpRelIso, this_met );

// 		       fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_L1d4",                tmpRelIso, this_met, this_weight );
// 		       fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_"+etaside+"_L1d4",    tmpRelIso, this_met, this_weight );
// 		       fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_uw_L1d4",             tmpRelIso, this_met, 1 );
// 		       fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_uw_"+etaside+"_L1d4", tmpRelIso, this_met, 1 );
// 		       fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_L1d4",             tmpRelIso, this_weight );
// 		       fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+metside+"_L1d4", tmpRelIso, this_weight );
// 		       fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+etaside+"_L1d4", tmpRelIso, this_weight );
// 		       fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+metside+"_"+etaside+"_L1d4", tmpRelIso, this_weight );

		       if(pass_eid_c0) {

			 iso_fillHisto_NES( 7, inBarrel, tmpRelIso, this_met );


// 			 iso_fillHisto_nlevel_nj_nmc( hh_QCDest_isoVmet_NES[7], tmpRelIso, this_met ) ; NEW 7:L1d5
// 			 if(isInBarrel) iso_fillHisto_nlevel_nj_nmc( hh_QCDest_isoVmet_NES_barrel[7], tmpRelIso, this_met );
// 			 else           iso_fillHisto_nlevel_nj_nmc( hh_QCDest_isoVmet_NES_endcap[7], tmpRelIso, this_met );


// 			 fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_L1d5",                tmpRelIso, this_met, this_weight );
// 			 fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_"+etaside+"_L1d5",    tmpRelIso, this_met, this_weight );
// 			 fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_uw_L1d5",             tmpRelIso, this_met, 1 );
// 			 fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_uw_"+etaside+"_L1d5", tmpRelIso, this_met, 1 );
// 			 fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_L1d5",             tmpRelIso, this_weight );
// 			 fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+metside+"_L1d5", tmpRelIso, this_weight );
// 			 fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+etaside+"_L1d5", tmpRelIso, this_weight );
// 			 fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+metside+"_"+etaside+"_L1d5", tmpRelIso, this_weight );

		       }//pass c0
		     }//pass c4
		   }//pass c3
		 }//pass c2
	       }//pass c1
	     }//pass d0
	   }//pass et,eta
	 }// loop over Nels
	 

	 if(  nGoodEle>0  ) {


	   // barrel or endcap?
	   if(ii_GoodEle_mostIso < 0) cout << "error/warning: could not find most isolated GoodEle." << endl;
	   string etaside = "barrel";
	   bool inBarrel = true;
	   if( fabs(els_eta->at(ii_GoodEle_mostIso))>1.56 ) {
	     etaside="endcap"; inBarrel=false;
	   }
	   /*
	     cout << "\n" << endl;
	     cout << "ele eta: " << els_eta->at(ii_GoodEle_mostIso) << "  (" << etaside<< ")" << endl;
	     cout << "\n" << endl;
	     cout << " CombRelIso = "<< CombRelIso << endl;
	     cout << " Check      = "<< (els_tIso->at(ii_GoodEle_mostIso)+els_cIso->at(ii_GoodEle_mostIso))/els_et->at(ii_GoodEle_mostIso) << endl;
	     cout << endl;
	   */

	   // cout << "** CombRelIso = " << CombRelIso << endl;

	   if(CombRelIso<0) {
	     static int mm = 0;
	     if(mm==0) cout << "(Printing the first 20 occurances.)"<< endl;
	     if(mm<10) {
	       cout << " info: negative CombRelIso: " << CombRelIso << endl;
	       ++mm;
	     }
	   }
	 
	   if(m_debug) cout << "-> Filling Reliso NES histograms, L2" << endl;


	   iso_fillHisto_NES( 8, inBarrel, CombRelIso, this_met );


// 	   iso_fillHisto_nlevel_nj_nmc( hh_QCDest_isoVmet_NES[8], CombRelIso, this_met ) ; //NEW 8: L2
// 	   if(inBarrel) iso_fillHisto_nlevel_nj_nmc( hh_QCDest_isoVmet_NES_barrel[8], tmpRelIso, this_met );
// 	   else           iso_fillHisto_nlevel_nj_nmc( hh_QCDest_isoVmet_NES_endcap[8], tmpRelIso, this_met );

// 	   fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_L2",                CombRelIso, this_met, this_weight );
// 	   fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_"+etaside+"_L2",    CombRelIso, this_met, this_weight );
// 	   fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_uw_L2",             CombRelIso, this_met, 1 );//unweighted
// 	   fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_uw_"+etaside+"_L2", CombRelIso, this_met, 1 );//unweighted
// 	   fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_L2",             CombRelIso, this_weight );
// 	   fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+metside+"_L2", CombRelIso, this_weight );
// 	   fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+etaside+"_L2", CombRelIso, this_weight );
// 	   fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+metside+"_"+etaside+"_L2", CombRelIso, this_weight );


	   if( !isMuon  &&  !isZ ) {

	     if(m_debug) cout << "-> Filling Reliso NES histograms, L3" << endl;

	     //	     iso_fillHisto_nlevel_nj_nmc( hh_QCDest_isoVmet_NES[9], CombRelIso, this_met ) ; //NEW 9: L3
	     iso_fillHisto_NES(9, inBarrel, CombRelIso, this_met);

// 	     fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_L3",                CombRelIso, this_met, this_weight );
// 	     fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_"+etaside+"_L3",    CombRelIso, this_met, this_weight );
// 	     fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_uw_L3",             CombRelIso, this_met, 1 );//unweighted
// 	     fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_uw_"+etaside+"_L3", CombRelIso, this_met, 1 );//unweighted
// 	     fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_L3",             CombRelIso, this_weight );
// 	     fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+metside+"_L3", CombRelIso, this_weight );
// 	     fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+etaside+"_L3", CombRelIso, this_weight );
// 	     fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+metside+"_"+etaside+"_L3", CombRelIso, this_weight );

	     if( !isConversion  &&  !isDifferentInteraction ) { //NB: no HT cut

	       if(m_debug) cout << "-> Filling Reliso NES histograms, L4" << endl;

	       //iso_fillHisto_nlevel_nj_nmc( hh_QCDest_isoVmet_NES[10], CombRelIso, this_met ) ; //NEW 10: L4
	       iso_fillHisto_NES(10, inBarrel, CombRelIso, this_met);

// 	       fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_L4",                CombRelIso, this_met, this_weight );
// 	       fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_"+etaside+"_L4",    CombRelIso, this_met, this_weight );
// 	       fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_uw_L4",             CombRelIso, this_met, 1 );//unweighted
// 	       fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_isoVmet_NES_uw_"+etaside+"_L4", CombRelIso, this_met, 1 );//unweighted
// 	       fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_L4",             CombRelIso, this_weight );
// 	       fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+metside+"_L4", CombRelIso, this_weight );
// 	       fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+etaside+"_L4", CombRelIso, this_weight );
// 	       fillHistoNjet_DataAndMC( "QCD_estimation/NES/QCDest_CombRelIso_NES_"+metside+"_"+etaside+"_L4", CombRelIso, this_weight );
	     }
	   }
	 }//nGoodEle > 0
       }//goodrun, trigger
     }//m_plotRelisoNES switch







     // study Z veto ---------
     if( m_studyZveto &&
	 goodrun  &&  fired_single_em  &&  nGoodEle>0  && 
     	 !isMuon  &&  !isZ  &&  !isConversion  &&  !isDifferentInteraction  &&  ht >= HTCUT ) {

       if(m_debug) cout << "study z veto" << endl;
       if(isZee) {
	 nEvent_Zee_pass++;
	 if( this_met > METCUT ) nEvent_Zee_highMET++;
	 else                    nEvent_Zee_lowMET++;
       }	 
       else if(isZmm) {
	 nEvent_Zmm_pass++;
	 if( this_met > METCUT ) nEvent_Zmm_highMET++;
	 else                    nEvent_Zmm_lowMET++;
       }
       else if(isZtt) {
	 nEvent_Ztt_pass++;
	 if( this_met > METCUT ) nEvent_Ztt_highMET++;
	 else                    nEvent_Ztt_lowMET++;
       }      
     }// study z veto
   

     // Add m3 plots for W+jets estimation
     //---------------------------------------
     // Plot for event with 4 or more jets only
     // Require at least 3 jest to reconstruct t->jjj

     if(m_debug) cout << "compute m3" << endl;

     // notes: all cuts except ISO
     if ( goodrun  &&  fired_single_em  &&  nGoodEle>0  &&  this_met > METCUT &&
	  !isMuon  &&  !isZ  &&  !isConversion  &&  !isDifferentInteraction  &&  ht >= HTCUT &&
	  m_nGoodJet >=4 ) { 
       reco_hadronicTop_highestTopPT( jets, nGoodIsoEle );
     }



     //-------------------
     // Add Delta R(e,mu)
     //-------------------
     if(m_debug) cout << "DR(e,mu)" << endl;
     // Get isolated electrons     
     for ( size_t e=0; e < iso_electrons.size(); e++ ) {       
       
       //TLorentzVector ele( els_px->at(e), els_py->at(e), els_pz->at(e), els_energy->at(e) );       
       //cout << "ele pt : "<< iso_electrons.at(e).Pt()  << endl;
       if(iso_electrons.at(e).Pt()<1)  cout << "iso ele pt < 1, ev:"<<ev<< endl;
       // Get all GoodMuon
       for ( size_t m=0; m < muons.size(); m++ ) {
	 //TLorentzVector muo( mus_cm_px->at(m), mus_cm_py->at(m), mus_cm_pz->at(m), mus_energy->at(m) );
	 //cout << "mu pt : "<< muons.at(m).Pt()  << endl;
	 if(muons.at(m).Pt()<1)  cout << "mu pt < 1, ev:"<<ev<< endl;
	 float DRemu = (iso_electrons.at(e)).DeltaR( muons.at(m) );
	 if(ev<100) cout << "ev: " << ev << ", DRemu: " << DRemu << endl;
	 fillHistoDataAndMC( h_DRemu_selE_GoodMu, DRemu );
	 if ( e_plus_jet_pass && m_nGoodJet>=4 ) { // pass ALL cuts
	   fillHistoDataAndMC( h_DRemu_selE_GoodMu_pass, DRemu );
	   if( DRemu < 0.1 ) nEvent_DR_ele_muo_less_than_01 ++;
	 }
       }//mu
     }//e




     //------------------------------------
     // Delta phi (met,nearest jet) 2-3-10
     //------------------------------------
     if(m_debug) cout << "Dphi(met,jet)" << endl;
     if( m_nGoodJet > 0 ) {
       if ( goodrun && fired_single_em && nGoodEle>0 ) {
	 Double_t min_dphi_mu = 10;
	 Double_t min_dphi_t1 = 10;
	 Double_t min_dphi_gen = 10;
	 for(size_t j=0; j<jets.size(); j++){
	   // compute dphi(m,j)
	   TVector2 j2v( jets.at(j).Px(), jets.at(j).Py() );
	   min_dphi_mu = TMath::Min( min_dphi_mu, TMath::Abs(met2v_mu.DeltaPhi(j2v)) );
	   min_dphi_t1 = TMath::Min( min_dphi_t1, TMath::Abs(met2v_t1.DeltaPhi(j2v)) );
	   if(!IsData()) { min_dphi_gen = TMath::Min( min_dphi_gen, TMath::Abs(met2v_gen.DeltaPhi(j2v)) ); }
	   if(min_dphi_mu>3.2) cout << "INFO: min_dphi_mu > 3.2!" << endl;
	   if(min_dphi_t1>3.2) cout << "INFO: min_dphi_t1 > 3.2!" << endl;
	 }
	 fillHisto_Njet_DataAndMC( h_DPhiMetJet_mu_goodE, min_dphi_mu, this_weight );
	 fillHisto_Njet_DataAndMC( h_DPhiMetJet_t1_goodE, min_dphi_t1, this_weight );
	 if(!IsData()) { fillHisto_Njet_DataAndMC( h_DPhiMetJet_gen_goodE, min_dphi_gen, this_weight ); }

	 // apply all but MET cut (N-1)
	 if(  nGoodIsoEle==1  &&  !isMuon  &&  !isZ  &&  !isConversion 
	      && ( ( m_rejectEndcapEle==false && !isDifferentInteraction ) ||  // PV check (DIFFZ)                        
		   ( m_rejectEndcapEle==true  && isBarrel) )                   // or ele eta cut  
	      ) {
	   fillHisto_Njet_DataAndMC( h_DPhiMetJet_mu,  min_dphi_mu,  this_weight );
	   fillHisto_Njet_DataAndMC( h_DPhiMetJet_t1,  min_dphi_t1,  this_weight );
	   if(!IsData()) { fillHisto_Njet_DataAndMC( h_DPhiMetJet_gen, min_dphi_gen, this_weight ); }
	 }//met N-1

       }//one good ele
     }//one good jet
     //------------------------------------
     


     if(m_debug) {
       cout << "\n---------------------------";
       cout << "\n End: Clearing vectors";
       cout << "\n---------------------------\n";
     }
     electrons.clear();
     ii_electrons.clear();
     electrons_barrel.clear();
     electrons_endcap.clear();
     iso_electrons.clear();
     iso_electrons_barrel.clear();
     iso_electrons_endcap.clear();
     muons.clear();
     iso_muons.clear();
     jets.clear();
     met.clear();
     //tagged_jets.clear();
     electrons_isoval.clear();
     electrons_isoval2.clear();
     Zele.clear();

     //---------------------------
     // PDF uncertainties (begin)
     //---------------------------
     if( !IsData() && isTTbar && m_studyPDFunc ){
       if(ev==1) {
	 cout << "PDF uncertainties" << endl;
	 cout << "PDFWcteq66_0   " << PDFWcteq66_0 << endl;
	 cout << "PDFWcteq66_1   " << PDFWcteq66_1 << endl;
	 cout << "PDFWcteq66_2   " << PDFWcteq66_2 << endl;
	 cout << "PDFWcteq66_3   " << PDFWcteq66_3 << endl;
	 cout << "PDFWcteq66_4   " << PDFWcteq66_4 << endl;
	 cout << "PDFWcteq66_10  " << PDFWcteq66_10 << endl;
	 cout << "PDFWcteq66_20  " << PDFWcteq66_20 << endl;
	 cout << "PDFWcteq66_30  " << PDFWcteq66_30 << endl;
	 cout << "PDFWcteq66_40  " << PDFWcteq66_40 << endl;
	 cout << "PDFWcteq66_44  " << PDFWcteq66_44 << endl;
       }
       fillHisto_PDF_weights(h_pdf_total);
       if(e_plus_jet_pass) fillHisto_PDF_weights(h_pdf_pass);
     }
     //-------------------------
     // PDF uncertainties (end)
     //-------------------------     


     //--------------
     //  B-tagging
     //--------------
     if(e_plus_jet_pass && m_nGoodJet>=4) { //pass all cuts
       fillHistoDataAndMC( h_nbtag_TCHE, m_nbtag_TCHE );
       fillHistoDataAndMC( h_nbtag_TCHP, m_nbtag_TCHP );
       fillHistoDataAndMC( h_nbtag_SSV,  m_nbtag_SSV  );
     }
     //--------------



   }// loop over all the events

   //----------------------------------------------------------------------------
   //                         End  of  Event  Loop
   //----------------------------------------------------------------------------


   


   if(!IsData()){
     cout << "\n (num of event flagged as Z, di-e inv mass)" << endl;
     cout << " nZ_ORIG (old) = " << nZ_ORIG << endl;
     cout << " nZ_TEST (new) = " << nZ_TEST << endl << endl;
   }


   const bool study_z_veto = false;
   if(study_z_veto){
     cout << endl;
     cout <<  "nEvent looped over (Z+jets)            = "<<  nEvent_Z << endl;
     cout <<  "|->nEvent which is really Z->ee        = "<<  nEvent_Zee << endl;
     cout <<  "|->nEvent which is really Z->mumu      = "<<  nEvent_Zmm << endl;
     cout <<  "|->nEvent which is really Z->tautau    = "<<  nEvent_Ztt << endl;
     cout << endl;
     cout <<  "|->nEvent which is really Z->ee        = "<<  nEvent_Zee << endl;
     cout <<  "  |-> 2 or more reco-ele (no cut)      = "<<  nEvent_2orMoreEle << endl;
     cout <<  "    |-> subset: match to genEle Z->ee  = "<<  nEvent_EleMatch << " (" 
	  << (float)nEvent_EleMatch/nEvent_2orMoreEle*100.0 << " %)"<< endl;
     cout <<  "       |-> pass cut (2 or more):"<<  endl;
     cout <<  "       |-> subset: et only             = "<<  nEvent_Z_et << endl; 
     cout <<  "       |-> subset: eta only            = "<<  nEvent_Z_eta << endl; 
     cout <<  "       |-> subset: et,eta              = "<<  nEvent_Z_eteta << endl; 
     cout <<  "       |-> subset: et,eta,d0           = "<<  nEvent_Z_eteta_d0 << endl; 
     cout <<  "       |-> subset: et,eta,d0,idL       = "<<  nEvent_Z_eteta_d0_L<< endl;  
     cout <<  "       |-> subset: et,eta,d0,idT       = "<<  nEvent_Z_eteta_d0_T<< endl;  
     cout <<  "       |-> subset: et,eta,d0,idRL      = "<<  nEvent_Z_eteta_d0_RL<< endl; 
     cout <<  "       |-> subset: et,eta,d0,idRT      = "<<  nEvent_Z_eteta_d0_RT<< endl;
     cout <<  "|-> true Z, reco-ele, pass et,eta:"<< endl;
     cout <<  "    nEvent     |-> 0 reco-e            = "<<h_nEle_s1[6][10]->GetBinContent(1)<<endl; //allj
     cout <<  "               |-> 1 reco-e            = "<<h_nEle_s1[6][10]->GetBinContent(2)<<endl;
     cout <<  "               |-> 2 reco-e            = "<<h_nEle_s1[6][10]->GetBinContent(3)<<endl;
     cout <<  "               |-> 3 reco-e            = "<<h_nEle_s1[6][10]->GetBinContent(4)<<endl;
     cout <<  "               |-> 4 reco-e            = "<<h_nEle_s1[6][10]->GetBinContent(5)<<endl;
     cout <<  "|-> true Z->ee, GenEle pass et,eta:"<< endl;
     cout <<  "    nEvent     |-> 0 gen-e             = "<<h_nGenBasicEle_Zee_allj->GetBinContent(1)<<endl;
     cout <<  "               |-> 1 gen-e             = "<<h_nGenBasicEle_Zee_allj->GetBinContent(2)<<endl;
     cout <<  "               |-> 2 gen-e             = "<<h_nGenBasicEle_Zee_allj->GetBinContent(3)<<endl;
     cout <<  "               |-> 3 gen-e             = "<<h_nGenBasicEle_Zee_allj->GetBinContent(4)<<endl;
     cout <<  "               |-> 4 gen-e             = "<<h_nGenBasicEle_Zee_allj->GetBinContent(5)<<endl;
     cout <<  "Zee --> pass QCDest scope (reliso)     = "<<  nEvent_Zee_pass << endl;
     cout <<  "    |-> high met                       = "<<  nEvent_Zee_highMET << endl;
     cout <<  "    |-> low met                        = "<<  nEvent_Zee_lowMET << endl;
     cout <<  "Zmm --> pass QCDest scope (reliso)     = "<<  nEvent_Zmm_pass << endl;
     cout <<  "    |-> high met                       = "<<  nEvent_Zmm_highMET << endl;
     cout <<  "    |-> low met                        = "<<  nEvent_Zmm_lowMET << endl;
     cout <<  "Ztt --> pass QCDest scope (reliso)     = "<<  nEvent_Ztt_pass << endl;
     cout <<  "    |-> high met                       = "<<  nEvent_Ztt_highMET << endl;
     cout <<  "    |-> low met                        = "<<  nEvent_Ztt_lowMET << endl;
     cout << endl;
   }
   //---------
   fprintf(outfile, "\nList of files that contain selected events:\n\n");
   for(set<string>::const_iterator it = interestingFiles.begin(); it != interestingFiles.end(); ++it) {
     fprintf(outfile, "%s\n", it->c_str());
   }
   fclose(outfile);

   // MC - make bin zero the total bin for data-like print-out of numbers of events
   if(!IsData()){
     for(int i=0; i<nstage; ++i) {
       for(int j=0; j<ntjet; ++j){
	 for(int k=1; k<nmctype; ++k){
	   e_plus_jet[i][j][0]          += e_plus_jet[i][j][k];
	   e_plus_jet_weighted[i][j][0] += e_plus_jet_weighted[i][j][k];
	 }
       }
     }
   }

   // Fill histograms (event count tables)
   fillHisto_event_tables();


   // Fill Signal accpetance and efficiency histograms
   if(mc_sample_has_ttbar){
     float total = 0;
     float totalb = 0;
     for(int k=1;k<11;++k){
       totalb += e_plus_jet[9][4][k];
       for(int j=0;j<ntjet;++j){
	 total += e_plus_jet[0][j][k];
       }
     }
     
     const double SignalAcc = totalb/total;
     const double SignalAccUnc = sqrt( SignalAcc*(1-SignalAcc)/totalb );     
 
     SignalVars->SetBinContent(1,SignalAcc);
     SignalVars->SetBinContent(2,SignalAccUnc);

     SignalVars->SetBinContent(3,totalb); //signal passed event selection
     SignalVars->SetBinContent(4,e_plus_jet[9][4][1]);//evqq channel passed event selection
   }


   // Histogram for PDF uncertainties
   if(m_studyPDFunc){
     h_pdf_eff->Divide(h_pdf_pass,h_pdf_total,1.,1.,"B");
   }


   //Summary Info
   cout << "\n\\begin{verbatim}";
   cout << "\n%------------------------------------------------------------------------------------\n";
   cout << "                             Summary Information";
   cout << "\n%------------------------------------------------------------------------------------\n";

   cout << "Integrated luminosity assumed = " << intlumi << "/pb" << endl;
   if( GetTrigger() ) cout << "Applied HLT trigger: " << HLTBit << endl;
   cout << " Electron ET cut =  " << ELE_ETCUT << "  GeV" << endl;
   cout << "     Muon PT cut =  " << MU_PTCUT  << "  GeV" << endl;
   cout << "      Jet PT cut =  " << JET_PTCUT << "  GeV" << endl;
   cout << "      MET    cut =  " << METCUT    << "  GeV" << endl;
   //cout << "       HT    cut =  " << HTCUT     << "  GeV" << endl;
   cout << " Electron ID = " << printEleID() << endl;
   cout << " Electron RelIso formula = " ;
   if(useNewReliso) cout << "new"; else cout << "old";
   if(!m_applyMETcut && m_rejectEndcapEle) 
     cout << "\n Using plan B, does not use MET, use barrel ele only (eta<1.442)" << endl;
   if(m_applyMETcut && m_rejectEndcapEle) 
     cout << "\n Using plan C, both MET & ele eta cuts are used (eta<1.442)" << endl;
   cout << "\n\\end{verbatim}\n"<< endl;


   cout << "\n Monitor:" << endl;
   cout << " Number of event with DR(e,mu) less than 0.1  =  " << nEvent_DR_ele_muo_less_than_01 << endl;
   cout << " e: isolated electron" << endl;
   cout << " mu: all GoodMuon\n" << endl;     



   DrawEventPerNjetTable();  //unweighted number
   if(IsData()==false) {
     DrawEventPerNjetTable();  //weighted number
   }


   // FB+TL (21-1-09)
   // When running on MC, produce tables of event count according to type of MC
   if(IsData()==false) {

     cout << endl;
     DrawSignalBGTable(); //weighted table

     DrawMCTypeTable( "Events Table (per MC type; unweighted)" ); //unweighted table
     DrawMCTypeTable( "Events Table (per MC type; weighted)  " ); //weighted table


     // Break down table for tt+j Alpgen
     if( mc_sample_has_ttbar && signal_is_Alpgen ) {
       DrawTTnjetTable( "TT+j Table (unweighted)" );
       DrawTTnjetTable( "TT+j Table (weighted)  " );
     }

     // Break down table for QCD
     if( mc_sample_has_QCD ) {
       //       DrawQCDTable( "QCD Table (weighted)  ", ve );
       //       DrawQCDTable( "QCD Table (unweighted)", ve );
       DrawQCDTable( "QCD Table (unweighted)" );
       DrawQCDTable( "QCD Table (weighted)  " );

     }

     // Break down table for single top
     if( mc_sample_has_singleTop ) {
       DrawSingleTopTable( "Single Top Table (unweighted)" );
       DrawSingleTopTable( "Single Top Table (weighted)  " );
     }


     // Acceptance note printout (ttbar)
     // >= 3 tight jets, >= 4 tight jets

     // TL (21-1-09):  only run this if we've ran over some ttbar MC in input
     if( mc_sample_has_ttbar && !signal_is_Alpgen ) { //FIXME

       DrawSignalAcceptanceTable(ve);

       // Set labels for all_mctype histogram (ttbar decay modes)
       sig_all_mctype->GetXaxis()->SetBinLabel(1, "evqq"); //semilep
       sig_all_mctype->GetXaxis()->SetBinLabel(2, "mvqq");
       sig_all_mctype->GetXaxis()->SetBinLabel(3, "tvqq");
       sig_all_mctype->GetXaxis()->SetBinLabel(4, "evev"); //dilepton
       sig_all_mctype->GetXaxis()->SetBinLabel(5, "mvmv");
       sig_all_mctype->GetXaxis()->SetBinLabel(6, "tvtv");
       sig_all_mctype->GetXaxis()->SetBinLabel(7, "evmu");
       sig_all_mctype->GetXaxis()->SetBinLabel(8, "evtv");
       sig_all_mctype->GetXaxis()->SetBinLabel(9, "mvtv");
       sig_all_mctype->GetXaxis()->SetBinLabel(10,"qqqq"); //fully hadronic

     }//end mc_sample_has_ttbar


     // (19 Feb 09) print how many QCD events we have after all cuts except reliso & njet
     if (mc_sample_has_QCD){
       const int QCD_bc = 2;
       cout << "\n QCD" << endl;
       cout << " New RelIso  mc" << setw(10) << intlumi <<"/pb" << endl;
       cout << "   1j" 
	    << setw(10) << h_QCDest_CombRelIso[1][QCD_bc]->GetEntries() 
	    << setw(10) << h_QCDest_CombRelIso[1][QCD_bc]->Integral() << endl;       
       cout << "   2j" 
	    << setw(10) << h_QCDest_CombRelIso[2][QCD_bc]->GetEntries() 
	    << setw(10) << h_QCDest_CombRelIso[2][QCD_bc]->Integral() << endl;
       cout << "   3j" 
	    << setw(10) << h_QCDest_CombRelIso[3][QCD_bc]->GetEntries() 
	    << setw(10) << h_QCDest_CombRelIso[3][QCD_bc]->Integral() << endl;
       cout << " >=4j" 
	    << setw(10) << h_QCDest_CombRelIso[5][QCD_bc]->GetEntries() 
	    << setw(10) << h_QCDest_CombRelIso[5][QCD_bc]->Integral() << endl;
     }

   }//end MC

      
   cout.precision(myprec); //reset precision

   //Print event tables with errors 
   PrintErrorTables( ve );
 
   //Print event tables with errors, but without applying scientific notation
   ScientificNotation = false;
   PrintErrorTables( ve );


   // Close the histogram file
   histf->Write();
   histf->Close();

   //check histo vector (is it worth clearing the histo vectors???)
   cout << "h_nEle_all.size(): "<< h_nEle_all.size() << endl;
   cout << "h_exp_ele_et.size(): "<< h_exp_ele_et.size() << endl;
   cout << "h_ele_ET.size(): "<< h_ele_ET.size() << endl;
   

   //856
   if( DoConversionStudies() ) {PrintConversionTable();}
   //end 856

   cout << "\n***********************************************************************";
   cout << "\n*                                                                     *";
   cout << "\n*           A N A L Y S I S       C O M P L E T E D                   *";
   cout << "\n*                                                                     *";
   cout << "\n***********************************************************************";
   cout << "\n*   Available events:  " << left << setw(47) << nEventsAvail << "*";
   cout << "\n*   Limit:             " << left << setw(47) << nEvents      << "*";
   cout << "\n*   Processed events:  " << left << setw(47) << counter      << "*";
   cout << "\n*   Passed events:     " << left << setw(47) << counter_pass << "*";
   cout << "\n***********************************************************************\n";

   // Print current time
   TDatime now;
   cout << endl << "Current local ";
   now.Print();
   cout << endl;

   return true;

}// end EventLoop()
//=============================================================================================
//
//                                END   of   Event   Loop
//
//=============================================================================================





float ana::calcDeltaPhi(const TLorentzVector& p1,const TLorentzVector& p2) const {
  float delPhi(999.0);
  delPhi = fabs( p1.Phi() - p2.Phi() );
  if (delPhi >  TMath::Pi()) delPhi = 2.0*TMath::Pi() - delPhi;
  return delPhi;
}

float ana::calcDeltaEta(const TLorentzVector& p1,const TLorentzVector& p2) const {
  float delEta(999.0);
  delEta = p1.PseudoRapidity() - p2.PseudoRapidity();
  return delEta;
}

// Calculate Delta R
float ana::calcDeltaR(const TLorentzVector& p1, const TLorentzVector& p2) const {
  
  float phi1 = p1.Phi();
  float eta1 = p1.PseudoRapidity();
  
  float phi2 = p2.Phi();
  float eta2 = p2.PseudoRapidity();
  
  float delR=  calcDeltaR(phi1,eta1,phi2,eta2);
  
  //cout << " phi1= " << phi1 << " phi2= " << phi2
  //<< " eta1= " << eta1 << " eta2= " << eta2
  //<< " delR= " << delR << endl;
  
  return delR;
}   
 
//Calculate Delta R
float ana::calcDeltaR(const float phi1, const float eta1, const float phi2, const float eta2) const {
  
  // Calculate delta phi and delta eta
  float delPhi = fabs(phi1 - phi2);
  if (delPhi >  TMath::Pi()) delPhi = 2.0*TMath::Pi() -delPhi;
  float delEta = eta1-eta2;
  
  // Calculate the delta R
  float delR = sqrt(delPhi*delPhi+delEta*delEta);             
  
  return delR;      
}


//========================================================================================
//
//                                QCD Estimation (start)
//
//========================================================================================

void ana::EstimateQCD() {
  EstimateQCD( outputHistFileName );
}

bool ana::EstimateQCD( const string inputFile ) {

  cout << "--------------------------------------------------------" << endl;
  cout << "                    QCD Estimation" << endl;
  cout << "--------------------------------------------------------" << endl;

  string isoVariable = "QCDest_CombRelIso"; //new reliso
  if(!useNewReliso) isoVariable = "QCDest_NormCombRelIso"; //old reliso

  // Open the histogram file to add histogram
  TFile f( inputFile.c_str(), "UPDATE" );
  outfile = fopen( outputTextFileName.c_str(), "a" ); //append

  TCanvas c1("c1","Electron isolation",980,700);
  c1.Divide(2,2);


  //f.Print();
  const int nrange = 1;// number of trial fit ranges
  //  const double fit_range_low[nrange] = {    0,  0.5 };
  //  const double fit_range_up[nrange]  = { 0.85, 0.85 };

  const int nj = 4;

  double sig_from = 0; //new reliso
  double sig_upto = 0.1;
  if(!useNewReliso) { //old reliso
    sig_from = 0.85;
    sig_upto = 1.0;
  }
  if(useNewReliso)  cout << "Using New Reliso formulation" << endl;
  else              cout << "Using Old Reliso formulation" << endl;
  cout << "Signal region: " << sig_from <<  " to " << sig_upto << endl;


  // use single function with fixed range
  const double fit_range_low_oldReliso[nj] = { 0.4, 0.4, 0.4, 0.4 };
  const double fit_range_up_oldReliso[nj]  = { 0.8, 0.8, 0.8, 0.8 };

  const double fit_range_low_newReliso[nj] = { 0.2, 0.2, 0.2, 0.2 }; //<-----
  const double fit_range_up_newReliso[nj]  = { 1.0, 1.0, 1.0, 1.0 }; //<-----

  double fit_range_low[nj] ;
  double fit_range_up[nj]  ;

  int    rebin_oldReliso[nj] = { 1, 2,  5, 10 };
  //  int    rebin_newReliso[nj] = { 2, 2,  5, 10 };  //<---- 
  int    rebin_newReliso[nj] = { 10, 10,  10, 10 };  //<---- 
  int    rebin[nj] ;

  for (int i=0; i<nj; ++i){ //set the fit range according to reliso choice
    if (useNewReliso) {
      fit_range_low[i] = fit_range_low_newReliso[i];
      fit_range_up[i]  = fit_range_up_newReliso[i];
      rebin[i]         = rebin_newReliso[i];
    }else{
      fit_range_low[i] = fit_range_low_oldReliso[i];
      fit_range_up[i]  = fit_range_up_oldReliso[i];
      rebin[i]         = rebin_oldReliso[i];
    }
  }

  string func = "gaus";
  //string func = "landau";
  //string func = "pol3";
  if(!useNewReliso) func = "gaus"; //old reliso
  cout << "Fit function: " << func << endl;

  //double nall_actual_sig[nrange][nj]; //all=s+b, for 0j,1j,2j,3j,4j,>=4j
  double nall_actual_ctr[nrange][nj];
  double nqcd_actual_sig[nrange][nj]; //qcd only
  double nqcd_actual_ctr[nrange][nj];
  double nqcd_actual_all[nj];
  double nqcd_actual_all_unweighted[nj];
  double n_extrap[nrange][nj];
  double n_extrap_err_plus[nrange][nj];
  double n_extrap_err_minus[nrange][nj];
  double fit_chi2[nrange][nj];
  int    fit_ndf[nrange][nj];

  // initialize arrays to zero
  for (int i=0; i<nrange; ++i){
    for (int j=0; j<nj; ++j){
      //nall_actual_sig[i][j] = 0;
      nall_actual_ctr[i][j] = 0;
      nqcd_actual_sig[i][j] = 0;
      nqcd_actual_ctr[i][j] = 0;
      nqcd_actual_all[j]    = 0;
      nqcd_actual_all_unweighted[j] = 0;
      n_extrap[i][j]           = 0;
      n_extrap_err_plus[i][j]  = 0;
      n_extrap_err_minus[i][j] = 0;
      fit_chi2[i][j] = 0;
      fit_ndf[i][j]  = 0;
    }
  }

  // take CombRelIso distributions, fit in range 0-0.85, extrapolate to 0.85-1.0
  const string myjetbin[4] = { "1 jet", "2 jets", "3 jets", "#geq4 jets"}; //root label
  const string jetlabel[4] = { "1j", "2j", "3j", "$\\ge$4j"}; //latex label
  const string jettext[4]  = { "1j", "2j", "3j", ">=4j"}; //text label

  for (int i=0; i<nrange; ++i) {

    /*
      cout << "-------------------- " << endl;
      cout << " Range: " << low << " to "<< up << endl;
      cout << "-------------------- " << endl;    
    */

    // 27 Feb 09
    vector<double> mpv;

    for (int j=0; j<nj; ++j) { //nj: 0,1,2,3,4,>=4j -> 1,2,3,>=4j


      const double fit_from = fit_range_low[j];
      const double fit_upto = fit_range_up[j];
      cout << "\n\n-----------------------" << endl;
      cout << " " << jettext[j] << "  fit range: " << fit_from << " to " << fit_upto << endl;


      //char *hname = Form("QCDest_NormCombRelIso_%dj",j+1);
      char *hname = Form("%s_%dj", isoVariable.c_str(), j+1);
      // if (j+1==4) hname = Form("QCDest_NormCombRelIso_4mj";
      if (j+1==4) hname = Form( "%s_4mj", isoVariable.c_str() );
      TH1F *allEvent = (TH1F*)f.Get(Form("QCD_estimation/%s__data",hname)); //all events
      TH1F *qcdEvent = (TH1F*)f.Get(Form("QCD_estimation/%s__QCD", hname)); //QCD only
      if (allEvent==0 || qcdEvent==0) {
	cout << "Warning: zero pointer to reliso histograms, exit EstimateQCD()" <<endl;
	return false;
      }



      const double binw = 0.01; //original binw = 0.01 = 1.1/110 (old) = 10/1000 (new)
      // signal bin (86, 100) for old reliso
      // signal bin ( 1, 100) for new reliso
      const int ctr_bin_low = (int)(fit_from/binw + 1);
      const int ctr_bin_up  = (int)(fit_upto/binw);
      const int sig_bin_low = (int)(sig_from/binw + 1);
      const int sig_bin_up  = (int)(sig_upto/binw);
      cout << "     signal bin:  " << sig_bin_low  << " to " << sig_bin_up ;
      cout << "  control bin:  " << ctr_bin_low  << " to " << ctr_bin_up << endl;

      //int sig_bin_low = 86;
      //int sig_bin_up =  100;
      //      n_actual_ctr[i][j] = allEvent->Integral( (int)(low*100), (int)(up*100) ); //actual number of _all_ event in control region
      nall_actual_ctr[i][j] = allEvent->Integral( ctr_bin_low, ctr_bin_up ); //actual number of _all_ event in control region
      //nall_actual_sig[i][j] = allEvent->Integral( sig_bin_low, sig_bin_up ); //actual number of _all_ event in signal region

      nqcd_actual_ctr[i][j] = qcdEvent->Integral( ctr_bin_low, ctr_bin_up ); //actual number of _QCD_ event in control region
      nqcd_actual_sig[i][j] = qcdEvent->Integral( sig_bin_low, sig_bin_up ); //actual number of _QCD_ event in signal region
      /*
	cout <<  "qcdEvent->Integral( ctr_bin_low, ctr_bin_up ) " <<qcdEvent->Integral( ctr_bin_low, ctr_bin_up ) << endl; 
	cout <<  "qcdEvent->Integral( 41, 80 )                  " << qcdEvent->Integral( 41, 85 ) << endl; 
	cout <<  "qcdEvent->Integral( sig_bin_low, sig_bin_up ) " <<qcdEvent->Integral( sig_bin_low, sig_bin_up ) << endl; 
	cout <<  "qcdEvent->Integral( 86, 100                 ) " << qcdEvent->Integral( 86, 100 ) << endl; 

	cout <<  "allEvent->Integral( ctr_bin_low, ctr_bin_up ) " <<allEvent->Integral( ctr_bin_low, ctr_bin_up ) << endl; 
	cout <<  "allEvent->Integral( 41,80 )                   " << allEvent->Integral( 41, 80 ) << endl; 
	cout <<  "allEvent->Integral( sig_bin_low, sig_bin_up ) " <<allEvent->Integral( sig_bin_low, sig_bin_up ) << endl; 
	cout <<  "allEvent->Integral( 86, 100 )                 " << allEvent->Integral( 86, 100 ) << endl; 
      */

      nqcd_actual_all_unweighted[j] = qcdEvent->GetEntries();
      nqcd_actual_all[j]            = qcdEvent->Integral();
      nqcd_actual_ctr[i][j]         = qcdEvent->Integral( ctr_bin_low, ctr_bin_up );


      const int this_rebin = rebin[j];


      // set value of global variable
      //QCDest_reliso_bin_width = binw*this_rebin; //binw = original bin width = 0.01
      Set_Reliso_bin_width( binw*this_rebin ); //binw = original bin width = 0.01
      cout << "     rebin = " << this_rebin << ", binw = " << Get_Reliso_bin_width() << endl;

      TH1F* all = (TH1F*)allEvent->Clone();
      TH1F* qcd = (TH1F*)qcdEvent->Clone();

      all->Rebin(this_rebin);  //all events
      qcd->Rebin(this_rebin);  //qcd only


      // Draw
      c1.cd(j+1);
      //all->SetFillColor(kGray);
      qcd->SetFillColor(kAzure-9); //blue
      qcd->SetLineColor(kAzure+2); //blue (darker)
      //h1->SetFillStyle(3001);
      all->Draw("ahist");
      if(useNewReliso) all->GetXaxis()->SetRangeUser(0,1.6);
      all->SetMinimum(0);
      all->SetTitle( Form("RelIso distribution (%s)", myjetbin[j].c_str()) );
      all->GetXaxis()->SetTitle("RelIso");
      //if(!useNewReliso) all->SetTitle( Form("Old RelIso distribution (%s)", myjetbin[j].c_str()) );



      // Do the fit if there are non-0 event in control region

      if ( nall_actual_ctr[i][j] > 0 ) {

	//cout << "do fit" << endl;
	//	h1->Fit("gaus","Q0","goff",0.,0.85); //"Q": Quiet, "0": do not draw
	//	h1->Fit("gaus","Q","h", low, up); //"Q": Quiet, "0": do not draw

	// For 3j, >=4j, set limit on MPV (for Landau only)
	//---------------------------------
	if( useNewReliso && j>1 && func=="landau") { //j=2 is 3j
	  
	  // MPV-Constrained Landau Fit	  
	  TF1 *fitf = new TF1("landau","landau",0,2);
	  
	  //fitf->SetParameters( 1, mpv.at(0), 0.1 ); //set initial values of parameters
	  //fitf->FixParameter( 1, mpv.at(0) );   // fix 2nd parameter is MPV
	  
	  // Set limit on MPV parameter for 3,>=4 jet bins
	  if(mpv.size()>=2){

	    cout << "mpv fitted: 1j= " << mpv.at(0) << ", 2j=" << mpv.at(1) << endl;
	    if( mpv.at(0) <= mpv.at(1) ) {
	      fitf->SetParLimits(1, mpv.at(0), mpv.at(1) );
	      //fitf->SetParLimits(1, 0.1, 1000 ); 
	    } else {
	      fitf->SetParLimits(1, mpv.at(1), mpv.at(0) );
	      //fitf->SetParLimits(1, 0.1, 1000 );
	    }
	  } else {
	    cout << "\n could not find 2 mpv values from 1j,2j fits.\n" << endl;
	  }	
	  
	  all->Fit(fitf,"V","ah", fit_from, fit_upto); //"Q": Quiet, "0": do not draw, "ah": (no axis)
	  //all->Fit(fitf,"BV","ah", fit_from, fit_upto); //"BV": if fixing parameter value
	  delete fitf;//test
	}
	else{
	  // Free Fit
	  all->Fit(func.c_str(),"V","ah", fit_from, fit_upto); //"Q": Quiet, "0": do not draw, "ah": (no axis)
	}

	qcd->Draw("ahist same"); //qcd part, no axis
        qcd->Draw("ae same");
        all->Draw("ae same");
        all->Draw("axis same");


	TF1 *myf = all->GetFunction(func.c_str());
	myf->SetLineColor(kRed);
	myf->Draw("same");

	TF1 *myf2 = (TF1*)myf->Clone();
	myf2->SetLineColor(kBlue);
	myf2->SetRange( sig_from, sig_upto ); //new
	//if(useNewReliso) myf2->SetRange( 0,    0.1 ); //new
	//else             myf2->SetRange( 0.85, 1.0 ); //old
	myf2->Draw("same");

	// resize the stat box
	gPad->Update();
	TPaveStats *st = (TPaveStats*)all->FindObject("stats");
	st->SetX1NDC(0.6);
	st->SetY1NDC(0.5);
	st->SetX2NDC(0.97);
	st->SetY2NDC(0.97);
	st->Draw();

	// Get fit results
	//------------------
	fit_chi2[i][j] = myf->GetChisquare();
	fit_ndf[i][j]  = myf->GetNDF();
	// Extrapolate fitted function

	//double bin_width_now = 1.1/110*this_rebin;  //xxxxxxxxxxx
	//cout << "QCDest_reliso_bin_width (2) "<< QCDest_reliso_bin_width << endl;
	
	//n_extrap[i][j] = myf->Integral( sig_from, sig_upto ) / QCDest_reliso_bin_width;
	n_extrap[i][j] = myf->Integral( sig_from, sig_upto ) / Get_Reliso_bin_width();
	//if(useNewReliso)  n_extrap[i][j] = myf->Integral( 0,    0.1 ) / QCDest_reliso_bin_width;
	//else              n_extrap[i][j] = myf->Integral( 0.85, 1.0 ) / QCDest_reliso_bin_width;


	//cout << "n (signal region) extrapolated (est): " << n_extrap[i][j] << endl;
	//cout << "n (signal region) actual: " << n_actual_sig[i][j] << endl;	
	//cout << "n (control/fit region) actual: " << n_actual_ctr[i][j] << endl;

	// calculate uncertainty on n_extrap by varying the fit function
	if (func=="gaus"||func=="landau") {
	  pair<double,double> estimate_error = estimateQCD_computeFitResultUncertainty( n_extrap[i][j], myf );
	  n_extrap_err_plus[i][j]  = estimate_error.first ;
	  n_extrap_err_minus[i][j] = estimate_error.second ;
	}

	// Get fit parameters
	if(func=="landau") {
	  mpv.push_back(myf->GetParameter(1));
	}


	// Add legend
	//-------------
	TLegend *leg = new TLegend(0.22,0.55,0.6,0.88);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->AddEntry(all, "All events (s+b)");
	leg->AddEntry(qcd, "QCD events","f");
	leg->AddEntry(myf, "Fit to all events");
	leg->AddEntry(myf2,"Extrapolation");
	leg->Draw();
      
      }// end fit

      /*
      // jet bin label
      TLatex *la = new TLatex(0.25, 0.73, myjetbin[j].c_str());
      la->SetNDC(); //use NDC
      la->SetTextSize(0.1);
      la->Draw();
      */

    }// end nj loop: 0 to 4 jets

    cout << "-----------------------------------------------"<< endl;
    cout << "  " << intlumi << "/pb     rb      True    Estimate     Diff" << endl;
    printf("   1 jet:   %2d %10.1f  %10.1f  %6.1f %%\n", rebin[0], nqcd_actual_sig[i][0], n_extrap[i][0], (n_extrap[i][0]/nqcd_actual_sig[i][0]-1)*100 );
    printf("   2 jet:   %2d %10.1f  %10.1f  %6.1f %%\n", rebin[1], nqcd_actual_sig[i][1], n_extrap[i][1], (n_extrap[i][1]/nqcd_actual_sig[i][1]-1)*100 ); 
    printf("   3 jet:   %2d %10.1f  %10.1f  %6.1f %%\n", rebin[2], nqcd_actual_sig[i][2], n_extrap[i][2], (n_extrap[i][2]/nqcd_actual_sig[i][2]-1)*100 );
    printf(" >=4 jet:   %2d %10.1f  %10.1f  %6.1f %%\n", rebin[3], nqcd_actual_sig[i][3], n_extrap[i][3], (n_extrap[i][3]/nqcd_actual_sig[i][3]-1)*100 );
    cout << "-----------------------------------------------"<< endl;

    for(int j=0; j<4; ++j){
      cout << "Unc of QCD estimate (" << jettext[j] << "):  +" 
	   << n_extrap_err_plus[i][j] << " / -" <<  n_extrap_err_minus[i][j] << " events" << endl;
    }

    cout << "\nEstimated number (>=4 jets): " << n_extrap[0][3] << " events";
    cout << " (+" << n_extrap_err_plus[i][3] 
	 << "/-" <<  n_extrap_err_minus[i][3] << " events)" << endl;

    // save canvas (for each fit range)
    c1.SaveAs(Form("c1_reliso_fit_%d.ps", i+1)); // hard-coded filename for now, can easily be changed if we want

  }// end nrange loop
  c1.Close();


  // Print result table
  //---------------------
  fprintf(outfile, "\n%%- - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  fprintf(outfile, "%%    QCD Estimation Results\n");
  fprintf(outfile, "%%- - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  fprintf(outfile, "%%\\scriptsize\n"); //small font size
  fprintf(outfile, "\\addtolength{\\tabcolsep}{-0.9mm}\n");
  fprintf(outfile, "\\begin{tabular}{|r|rr|rr|rr|rcccc|}\\hline\n");
  fprintf(outfile, " & \\multicolumn{2}{c|}{Total QCD}\n");
  fprintf(outfile, " & \\multicolumn{2}{c|}{Control region}\n");
  fprintf(outfile, " & \\multicolumn{2}{c|}{Signal region (QCD)}\n");
  fprintf(outfile, " & \\multicolumn{5}{c|}{Fit results} \\\\\n");
  fprintf(outfile, "         &   True &     True &  True (All) &  True (QCD) & True");
  fprintf(outfile, " & Estimate & Diff & $\\chi^2$/dof & Function & Range & Bin\\\\\n");
  fprintf(outfile, "         & (MC) & (%2.0f/pb) & (%2.0f/pb) & (%2.0f/pb) & (%2.0f/pb) & (%2.0f/pb) &&&&&width \\\\ \\hline\n",
	  intlumi,intlumi,intlumi,intlumi,intlumi);

  for (int i=0; i<nrange; ++i) { //range

    for (int j=0; j<nj; ++j) {
      double diff = -1;
      if ( nqcd_actual_sig[i][j] > 0 )
	diff = ( n_extrap[i][j]/nqcd_actual_sig[i][j] -1 ) *100;
      
      //      fprintf(outfile, Form( "%8s & ", jetlabel[j].c_str() ) ); // %10.4g = 4 sig figures
      fprintf(outfile, "%8s & ", jetlabel[j].c_str() ); // %10.4g = 4 sig figures
      fprintf(outfile, "%6.0f & %9.1f & ",  nqcd_actual_all_unweighted[j], nqcd_actual_all[j] ); //total QCD
      fprintf(outfile, "%8.1f & %8.1f & ",  nall_actual_ctr[i][j], nqcd_actual_ctr[i][j] );  //control region: all, QCD
      fprintf(outfile, "%8.1f & ",          nqcd_actual_sig[i][j] );  //signal region: actual QCD

      if ( nall_actual_ctr[i][j] > 0 ) {
	fprintf(outfile, "%8.1f & ",  n_extrap[i][j] );  //QCD estimate
	fprintf(outfile, "%5.1f \\%% & %6.2f & ", diff, fit_chi2[i][j]/fit_ndf[i][j] );
	fprintf(outfile, "%6s & %2.1f--%2.2f & %5.2f\\\\\n", func.c_str(), fit_range_low[j], fit_range_up[j], rebin[j]/(110/1.1) );
	
      } else {
	//if zero actual event, print "-" in the extrapolation result column
	fprintf(outfile, "%8s & %8s & %6s & %6s & %9s & -\\\\\n", 
		"-", "-", "-","-","-");
      }
      
    }// end nj loop
    fprintf(outfile, "\\hline\n");
  }// end fit-ranges loop

  fprintf(outfile, "\\end{tabular}\n");
  fprintf(outfile, "%%- - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  fclose(outfile);

  cout << "\n  Actual number (>=4 jets): " << nqcd_actual_sig[0][3] << " +/- " << sqrt(nqcd_actual_sig[0][3]) 
       << " events.  (unc=sqrt(N))" << endl;
  

  //add QCD estimate results to root file
  TH1D *QCDEstimate = new TH1D("QCDEstimate","QCD estimate results (>=4j)",8,0,8);
  // (i) range 1: 0 to 0.85
  QCDEstimate->SetBinContent( 1, n_extrap[0][3]           );
  QCDEstimate->SetBinContent( 2, n_extrap_err_plus[0][3]  );
  QCDEstimate->SetBinContent( 3, n_extrap_err_minus[0][3] );
  //average error = (plus + minus)/2
  QCDEstimate->SetBinContent( 4, (n_extrap_err_plus[0][3]+n_extrap_err_minus[0][3])/2 ); //average
  // TL 12-2-09: no longer using gaussian for >=4j bin, error on fit, take 30% (ad-hoc value)
  //  QCDEstimate->SetBinContent( 4, n_extrap[0][3]*0.3 ); // assume 30% error

  /*
  // (ii) range 2: 0.5 to 0.85
  QCDEstimate->SetBinContent(5,n_extrap[1][3]);
  QCDEstimate->SetBinContent(6,n_extrap_err_plus[1][3]);
  QCDEstimate->SetBinContent(7,n_extrap_err_minus[1][3]);
  QCDEstimate->SetBinContent(8,(n_extrap_err_plus[1][3]+n_extrap_err_minus[1][3])/2); //average
  */

  QCDEstimate->Write();
  f.Close();

  cout << "--------------------------------------------------------" << endl;
  cout << "                    QCD Estimation (end)" << endl;
  cout << "--------------------------------------------------------" << endl;
  return true;

}// end of EstimateQCD()

//-------------------------------------------------------------------------------------------
// Method to calculate QCD estimate uncertainty, by varying each of the 3 parameters
// in the Gaussian fit (ie constant, mean, sigma), by 1sigma plus and minus.
// The asymmetric uncertainty is calculated by adding in quadrature individual deviations.
// -> adapt for the case when we use landau to fit (also 3 parameters: constant, MPV, sigma).
//-------------------------------------------------------------------------------------------
pair<double,double> ana::estimateQCD_computeFitResultUncertainty( const double est, TF1 *myf ) const {

  // get fit results
  double constant = myf->GetParameter(0); //1st parameter: constant
  double mean     = myf->GetParameter(1); //2nd parameter: mean (gaus) or MPV (landau)
  double sigma    = myf->GetParameter(2); //3rd parameter: sigma
  double constant_err = myf->GetParError(0);
  double mean_err     = myf->GetParError(1);
  double sigma_err    = myf->GetParError(2);
  double est_const_pos;
  double est_const_neg;
  double est_mean_pos;
  double est_mean_neg;
  double est_sigma_pos;
  double est_sigma_neg;
  double new_est[2];

  // 1) vary constant, get new estimate
  new_est[0] = estimateQCD_newEstimateByVaryingFitFunction( constant + constant_err, mean, sigma );
  new_est[1] = estimateQCD_newEstimateByVaryingFitFunction( constant - constant_err, mean, sigma );   
  pair<double,double> est_mean_pair = estimateQCD_assign_pos_neg_estimate( est, new_est );
  est_const_pos = est_mean_pair.first;
  est_const_neg = est_mean_pair.second;

  // 2) vary mean(or MPV), get new estimate
  new_est[0] = estimateQCD_newEstimateByVaryingFitFunction( constant, mean + mean_err, sigma );
  new_est[1] = estimateQCD_newEstimateByVaryingFitFunction( constant, mean - mean_err, sigma );
  est_mean_pair = estimateQCD_assign_pos_neg_estimate( est, new_est );
  est_mean_pos = est_mean_pair.first;
  est_mean_neg = est_mean_pair.second;

  // 3) vary sigma, get new estimate
  new_est[0] = estimateQCD_newEstimateByVaryingFitFunction( constant, mean, sigma + sigma_err );
  new_est[1] = estimateQCD_newEstimateByVaryingFitFunction( constant, mean, sigma - sigma_err );
  est_mean_pair = estimateQCD_assign_pos_neg_estimate( est, new_est );
  est_sigma_pos = est_mean_pair.first;
  est_sigma_neg = est_mean_pair.second;

  // 4) Compute total uncertainty (sum in quadrature)
  double unc_pos =  sqrt( TMath::Power( (est_const_pos - est), 2) + 
			  TMath::Power( ( est_mean_pos - est), 2) + 
			  TMath::Power( (est_sigma_pos - est), 2 ) ) ;
  double unc_neg =  sqrt( TMath::Power( (est_const_neg - est), 2) + 
			  TMath::Power( ( est_mean_neg - est), 2) + 
			  TMath::Power( (est_sigma_neg - est), 2 ) );
  //cout << "total unc +ve =  " << unc_pos << endl;
  //cout << "total unc -ve =  " << unc_neg << endl;
  pair<double,double> unc(unc_pos,unc_neg);
  //  unc.first  = unc_pos;
  //  unc.second = unc_neg;
  return unc;
}

double ana::estimateQCD_newEstimateByVaryingFitFunction( const double this_constant, const double this_mean, const double this_sigma ) const {
        
  // 1) create new function using each error of fit parameters
  TF1 *f1a;
  if(useNewReliso) f1a = new TF1("f1a","landau",0,1.1);
  else             f1a = new TF1("f1a","gaus", 0,1);

  f1a->SetParameter("Constant", this_constant );
  if(useNewReliso) f1a->SetParameter("MPV",    this_mean  ); //landau
  else             f1a->SetParameter("Mean",   this_mean  ); //gaus
  f1a->SetParameter("Sigma",    this_sigma    );

  double new_est = 0;
//   if(useNewReliso) new_est =  f1a->Integral(0,   0.1) / QCDest_reliso_bin_width;
//   else             new_est =  f1a->Integral(0.85,1.0) / QCDest_reliso_bin_width;
  if(useNewReliso) new_est =  f1a->Integral(0,   0.1) / Get_Reliso_bin_width();
  else             new_est =  f1a->Integral(0.85,1.0) / Get_Reliso_bin_width();


  //cout << "QCD new est = " << new_est << endl;//"  diff = " << new_est-est << endl;
  delete f1a;
  return new_est;
}
// Method to check whether the new estimate is larger or smaller than the nominal estimate
pair<double,double> ana::estimateQCD_assign_pos_neg_estimate( const double est, const double new_est[2] ) const {

  double pos = est;
  double neg = est;
  
  if ( new_est[0] > est && new_est[1] < est ) {      
    pos = new_est[0]; 
    neg = new_est[1]; 
  }
  else if ( new_est[0] < est && new_est[1] > est ) { 
    pos = new_est[1];  
    neg = new_est[0]; 
  }
  else {
    //cout << "--> Info: both new estimates are on the same side (nom=" << est << "): "
    //     << new_est[0]<< ", " <<new_est[1] << endl;     
    if ( new_est[0] > est && new_est[1] > est ) { //+ve side
      if ( new_est[0] >= new_est[1] ) { pos = new_est[0]; }
      else {                            pos = new_est[1]; }
    }
    if ( new_est[0] < est && new_est[1] < est ) { //-ve side
      if ( new_est[0] <= new_est[1] ) { neg = new_est[0]; }
      else {                            neg = new_est[1]; }
    }
  }
  pair<double,double> retval(pos,neg);
  return retval;
}
//========================================================================================
//
//                                  QCD Estimation (end)
//
//========================================================================================







//========================================================================================
//
//                                 Historam Manipulation
//
//========================================================================================



/*
// arrays

// ------------ my method to add a set of 1D histograms acc to njet ----------------------
void ana::addHistoNjet( TH1F* h[], const string name, const string ext, const string title,
			const int nbin, const float xlow, const float xup ){
  const string jetname[7]  = {"0j","1j","2j","3j","4j", "4mj", "allj"};
  const string jetlabel[7] = {"0j","1j","2j","3j","4j", ">=4j","allj"};
  for (unsigned int j=0; j<7; ++j) {
    char hname[70];
    char htitle[100];
    sprintf( hname,  "%s_%s%s",  name.c_str(), jetname[j].c_str(), ext.c_str() );
    sprintf( htitle, "%s (%s)",  title.c_str(), jetlabel[j].c_str() );
    h[j] = new TH1F(hname, htitle, nbin, xlow, xup);
    // Call Sumw2 to store sum of square of weights per bin
    // Important to get histogram error to reflect statistical error of MC events, rather than sqrt(nexp)
    h[j]->Sumw2();
  }
}

// ------------ my method to add a set of 2D histograms acc to njet ----------------------
void ana::addHistoNjet( TH2F* h[], const string name, const string ext, const string title,
			const int nbinx, const float xlow, const float xup,
			const int nbiny, const float ylow, const float yup ){
  const string jetname[7]  = {"0j","1j","2j","3j","4j", "4mj", "allj"};
  const string jetlabel[7] = {"0j","1j","2j","3j","4j", ">=4j","allj"};
  for (unsigned int j=0; j<7; ++j) {
    char hname[70];
    char htitle[100];
    sprintf( hname,  "%s_%s%s",  name.c_str(), jetname[j].c_str(), ext.c_str() );
    sprintf( htitle, "%s (%s)",  title.c_str(), jetlabel[j].c_str() );
    h[j] = new TH2F(hname, htitle, nbinx, xlow, xup, nbiny, ylow, yup);
    h[j]->Sumw2();
  }
}

//-------------- my method to fill 1D histograms acc to njet -----------------------------
void ana::fillHistoNjet(TH1F* h[], const float value, const double w) {
  
  h[6]->Fill(value, w); //allj
  if( m_nGoodJet < 5 ) h[ m_nGoodJet ]->Fill(value, w); //=0,1,2,3,4j
  if( m_nGoodJet > 3 ) h[5]->Fill(value, w); //4 or more jets
}


//-------------- my method to fill 2D histograms acc to njet -----------------------------
void ana::fillHistoNjet(TH2F* h[], const float v1, const float v2, const double w) {
  
  h[6]->Fill(v1, v2, w); //allj
  if( m_nGoodJet < 5 ) h[ m_nGoodJet ]->Fill(v1, v2, w); //=0,1,2,3,4j
  if( m_nGoodJet > 3 ) h[5]->Fill(v1, v2, w); //4 or more jets
}

//-------------- my method to fill 1D histograms acc to njet, GIVEN eventClass ------------------
void ana::fillHistoNjet2D(TH1F* h[][16], const int ec, const float value, const double w ) {

  h[6][ec]->Fill(value, w); //allj  
  if( m_nGoodJet < 5 )  h[ m_nGoodJet ][ec]->Fill(value, w); //0-4j
  if( m_nGoodJet > 3 )  h[5][ec]->Fill(value, w); //>=4j
}
//----------------------------------------------------------------------------------------







//----- my method to add a set of 1D histograms (for each type of MC when running on MC) ------
void ana::addHistoDataAndMC( TH1F* h[], const string name, const string title,
			     const int nbin, const float xlow, const float xup ) const {

  //  int nhisto = 1; //1 for real data
  //  if(!IsData()) nhisto = 16; //MC

  for (int i=0; i<ntype; ++i) {
    // only do for mc present in input
    if( i>0 && is_mc_present(i)==false ) continue;
    char hname[70];
    char htitle[100];
    sprintf( hname,  "%s__%s",  name.c_str(),  mcname[i].c_str() );
    sprintf( htitle, "%s (%s)", title.c_str(), mclabel[i].c_str() );
    h[i] = new TH1F(hname, htitle, nbin, xlow, xup);
    h[i]->Sumw2();
  }
}

//----- my method to add a set of 2D histograms (for each type of MC when running on MC) ------
void ana::addHistoDataAndMC( TH2F* h[], const string name, const string title,
			     const int nbin, const float xlow, const float xup,
			     const int nbiny, const float ylow, const float yup ) const {

  //  int nhisto = 1; //1 for real data
  //  if(!IsData()) nhisto = 16; //MC

  for (int i=0; i<ntype; ++i) {
    // only do for mc present in input
    if( i>0 && is_mc_present(i)==false ) continue;
    char hname[70];
    char htitle[100];
    sprintf( hname,  "%s__%s",  name.c_str(),  mcname[i].c_str() );
    sprintf( htitle, "%s (%s)", title.c_str(), mclabel[i].c_str() );
    h[i] = new TH2F(hname, htitle, nbin, xlow, xup, nbiny, ylow, yup);
    h[i]->Sumw2();
  }
}






//-------- my method to fill 1D histograms (acc to MC type when running on MC) ---------
void ana::fillHistoDataAndMC(TH1F* h[], const float value, const double w ) const {

  if(m_debug) cout << "\nStarting << fillHistoDataAndMC >>: " << h[0]->GetName() << endl;
  h[0]->Fill(value, w); //all events (data)
  if(!IsData()) { //run on MC
    if      (isTTbar)       h[1]->Fill(value,w);
    else if (isQCD) {       h[2]->Fill(value,w); //all QCD
      if      (isEnri1)     h[3]->Fill(value,w);
      else if (isEnri2)     h[4]->Fill(value,w);
      else if (isEnri3)     h[5]->Fill(value,w);
      else if (isBce1)      h[6]->Fill(value,w);
      else if (isBce2)      h[7]->Fill(value,w);
      else if (isBce3)      h[8]->Fill(value,w);
    } else if (isWjets)     h[9]->Fill(value,w);
    else if (isZjets)       h[10]->Fill(value,w);
    else if (isVQQ)         h[11]->Fill(value,w);
    else if (isSingleTop) { h[12]->Fill(value,w); //all single top
      if      (isTW)        h[13]->Fill(value,w);
      else if (isTchan)     h[14]->Fill(value,w);
      else if (isSchan)     h[15]->Fill(value,w);
    }
  }
  if(m_debug) cout << "End of << fillHistoDataAndMC >>" << endl;
}
//--------------------------------------------------------------------------------------
//-------- my method to fill 2D histograms (acc to MC type when running on MC) ---------
void ana::fillHistoDataAndMC(TH2F* h[], const float v1, const float v2, const double w ) const {

  if(m_debug) cout << "\nStarting << fillHistoDataAndMC 2D >>: " << h[0]->GetName() << endl;
  h[0]->Fill(v1,v2, w); //all events (data)
  if(!IsData()) { //run on MC
    if      (isTTbar)       h[1]->Fill(v1,v2,w);
    else if (isQCD) {       h[2]->Fill(v1,v2,w); //all QCD
      if      (isEnri1)     h[3]->Fill(v1,v2,w);
      else if (isEnri2)     h[4]->Fill(v1,v2,w);
      else if (isEnri3)     h[5]->Fill(v1,v2,w);
      else if (isBce1)      h[6]->Fill(v1,v2,w);
      else if (isBce2)      h[7]->Fill(v1,v2,w);
      else if (isBce3)      h[8]->Fill(v1,v2,w);
    } else if (isWjets)     h[9]->Fill(v1,v2,w);
    else if (isZjets)       h[10]->Fill(v1,v2,w);
    else if (isVQQ)         h[11]->Fill(v1,v2,w);
    else if (isSingleTop) { h[12]->Fill(v1,v2,w); //all single top
      if      (isTW)        h[13]->Fill(v1,v2,w);
      else if (isTchan)     h[14]->Fill(v1,v2,w);
      else if (isSchan)     h[15]->Fill(v1,v2,w);
    }
  }
  if(m_debug) cout << "End of << fillHistoDataAndMC 2D >>" << endl;
}
//--------------------------------------------------------------------------------------
*/




// //-------------- my method to fill 1D histograms acc to njet & mctype -----------------------------
// void ana::fillHistoNjet_DataAndMC(const string name, const float value, const double w ) {
  
//   if(m_debug) cout << "\nStarting << fillHistoNjet_DataAndMC >>" << endl;

//   //  const string jetname[7]  = {"0j","1j","2j","3j","4j", "4mj", "allj"};

//   // first fill for "allj", then fill for "njet (0-4)", then fill for ">=4j"
//   for (int i=0; i<3; ++i) {

//     string jetbin = "allj";
//     TH1F *h;
//     char hname[70];

//     if(i==1) {
//       if( m_nGoodJet>=0 && m_nGoodJet<5 )  jetbin = jetname[ m_nGoodJet ]; //0j to 4j
//       else continue;
//     }else if(i==2){
//       if( m_nGoodJet>=4) jetbin = "4mj"; //4 or more
//       else continue;
//     }

//     sprintf( hname,  "%s_%s",  name.c_str(), jetbin.c_str() );
//     if(m_debug) cout <<  "Filling histo: " << hname << endl;

//     h = (TH1F*)histf->Get( Form("%s", hname) );
//     if( m_debug && h==0 ) cout << "histo not found" << endl;

//     if(h>0) h->Fill(value, w);  

//     // 2a) all-jet, MC type
//     if(!IsData()) { //MC
//       if(isTTbar) {//ttbar
// 	h = (TH1F*)histf->Get( Form("%s__ttbar", hname) );
// 	if(h>0) h->Fill(value, w);

//       } else if(isWjets) { //wj
// 	h = (TH1F*)histf->Get( Form("%s__Wjet", hname) );
// 	if(h>0) h->Fill(value, w);

//       } else if(isZjets) { //zj
// 	h = (TH1F*)histf->Get( Form("%s__Zjet", hname ) );
// 	if(h>0) h->Fill(value, w);
    
//       } else if(isQCD) {
// 	h = (TH1F*)histf->Get( Form("%s__QCD", hname) );
// 	if(h>0) h->Fill(value, w);
// 	if(isEnri1){
// 	  h = (TH1F*)histf->Get( Form("%s__enri1", hname) );
// 	  if(h>0) h->Fill(value, w);
// 	} else if(isEnri2){
// 	  h = (TH1F*)histf->Get( Form("%s__enri2", hname) );
// 	  if(h>0) h->Fill(value, w);
// 	} else if(isEnri3){
// 	  h = (TH1F*)histf->Get( Form("%s__enri3", hname) );
// 	  if(h>0) h->Fill(value, w);
// 	} else if(isBce1){
// 	  h = (TH1F*)histf->Get( Form("%s__bce1", hname) );
// 	  if(h>0) h->Fill(value, w);
// 	} else if(isBce2){
// 	  h = (TH1F*)histf->Get( Form("%s__bce2", hname) );
// 	  if(h>0) h->Fill(value, w);
// 	} else if(isBce3){
// 	  h = (TH1F*)histf->Get( Form("%s__bce3", hname) );
// 	  if(h>0) h->Fill(value, w);
// 	}
//       }else if(isSingleTop) {
// 	h = (TH1F*)histf->Get( Form("%s__singleTop", hname) );
// 	if(h>0) h->Fill(value, w);
// 	if(isTW) {
// 	  h = (TH1F*)histf->Get( Form("%s__tW", hname) );
// 	  if(h>0) h->Fill(value, w);
// 	}else if(isTchan){
// 	  h = (TH1F*)histf->Get( Form("%s__tchan", hname) );
// 	  if(h>0) h->Fill(value, w);
// 	}else if(isSchan){
// 	  h = (TH1F*)histf->Get( Form("%s__schan", hname) );
// 	  if(h>0) h->Fill(value, w);
// 	}
//       }
//     }//if run on MC
//   }// loop: "allj", "0-4j",">=4j"
  
//   if(m_debug) cout << "End of << fillHistoNjet_DataAndMC >>" << endl;
// }

// //-------------- my method to fill 2D histograms acc to njet & mctype -----------------------------
// void ana::fillHistoNjet_DataAndMC(const string name, const float v1, const float v2, const double w ) {
  
//   if(m_debug) cout << "\nStarting << fillHistoNjet_DataAndMC (h2D) >>" << endl;

//   //  const string jetname[7]  = {"0j","1j","2j","3j","4j", "4mj", "allj"};

//   // first fill for "allj", then fill for "njet (0-4)", then fill for ">=4j"
//   for (int i=0; i<3; ++i) {

//     string jetbin = "allj";
//     TH2F *h;
//     char hname[70];

//     if(i==1) {
//       if(m_nGoodJet>=0 && m_nGoodJet<5)  jetbin = jetname[ m_nGoodJet ]; //0j to 4j
//       else continue;
//     }else if(i==2){
//       if( m_nGoodJet>=4) jetbin = "4mj"; //4 or more
//       else continue;
//     }

//     sprintf( hname,  "%s_%s",  name.c_str(), jetbin.c_str() );
//     if(m_debug) cout <<  "Filling histo: " << hname << endl;

//     h = (TH2F*)histf->Get( Form("%s", hname) );
//     if( m_debug && h==0 ) cout << "histo not found" << endl;
    
//     if(h>0) h->Fill(v1, v2, w);  

//     // 2a) all-jet, MC type
//     if(!IsData()) { //MC
//       if(isTTbar) {//ttbar
// 	h = (TH2F*)histf->Get( Form("%s__ttbar", hname) );
// 	if(h>0) h->Fill(v1, v2, w);

//       } else if(isWjets) { //wj
// 	h = (TH2F*)histf->Get( Form("%s__Wjet", hname) );
// 	if(h>0) h->Fill(v1, v2, w);

//       } else if(isZjets) { //zj
// 	h = (TH2F*)histf->Get( Form("%s__Zjet", hname) );
// 	if(h>0) h->Fill(v1, v2, w);
    
//       } else if(isQCD) {
// 	h = (TH2F*)histf->Get( Form("%s__QCD", hname) );
// 	if(h>0) h->Fill(v1, v2, w);
// 	if(isEnri1){
// 	  h = (TH2F*)histf->Get( Form("%s__enri1", hname) );
// 	  if(h>0) h->Fill(v1, v2, w);
// 	} else if(isEnri2){
// 	  h = (TH2F*)histf->Get( Form("%s__enri2", hname) );
// 	  if(h>0) h->Fill(v1, v2, w);
// 	} else if(isEnri3){
// 	  h = (TH2F*)histf->Get( Form("%s__enri3", hname) );
// 	  if(h>0) h->Fill(v1, v2, w);
// 	} else if(isBce1){
// 	  h = (TH2F*)histf->Get( Form("%s__bce1", hname) );
// 	  if(h>0) h->Fill(v1, v2, w);
// 	} else if(isBce2){
// 	  h = (TH2F*)histf->Get( Form("%s__bce2", hname) );
// 	  if(h>0) h->Fill(v1, v2, w);
// 	} else if(isBce3){
// 	  h = (TH2F*)histf->Get( Form("%s__bce3", hname) );
// 	  if(h>0) h->Fill(v1, v2, w);
// 	}
//       }else if(isSingleTop) {
// 	h = (TH2F*)histf->Get( Form("%s__singleTop", hname) );
// 	if(h>0) h->Fill(v1, v2, w);
// 	if(isTW) {
// 	  h = (TH2F*)histf->Get( Form("%s__tW", hname) );
// 	  if(h>0) h->Fill(v1, v2, w);
// 	}else if(isTchan){
// 	  h = (TH2F*)histf->Get( Form("%s__tchan", hname) );
// 	  if(h>0) h->Fill(v1, v2, w);
// 	}else if(isSchan){
// 	  h = (TH2F*)histf->Get( Form("%s__schan", hname) );
// 	  if(h>0) h->Fill(v1, v2, w);
// 	}
//       }
//     }//if run on MC
//   }// loop: "allj", "0-4j",">=4j"
  
//   if(m_debug) cout << "End of << fillHistoNjet_DataAndMC (h2D) >>" << endl;
// }

/*
//----- my method to add a set of 1D histograms (for each type of MC when running on MC) ------
void ana::addHisto_Njet_DataAndMC( TH1F* h[7][16], const string name, const string title,
				   const int nbin, const float xlow, const float xup ) {
  // example: xxx_[1-4j]__[data,ttbar,wj..]

  const string jetname[7]  = {"0j","1j","2j","3j","4j", "4mj", "allj"};
  const string jetlabel[7] = {"0j","1j","2j","3j","4j", ">=4j","allj"};

  //  int ntype = 1; //1 for real data
  //  if(!IsData()) ntype = 16; //MC

  for (int j=0; j<7; ++j) {
    for (int i=0; i<ntype; ++i) {
      // only do for mc present in input
      char hname[70];
      char htitle[100];
      if( i>0 && is_mc_present(i)==false ) continue;
      sprintf( hname,  "%s_%s__%s",    name.c_str(),  jetname[j].c_str(), mcname[i].c_str()  );
      sprintf( htitle, "%s (%s, %s)",  title.c_str(), jetlabel[j].c_str(), mclabel[i].c_str() );       
      h[j][i] = new TH1F(hname, htitle, nbin, xlow, xup);
      h[j][i]->Sumw2();      
    }
  }
}

//-------------- my method to fill  1D  histograms acc to njet & mctype ------------------
void ana::fillHisto_Njet_DataAndMC( TH1F* h[7][16], const float value, const double w ) {

  if(m_debug) cout << "\nStarting << fillHisto_Njet_DataAndMC >>: " << h[0][0]->GetName() << endl;

  if(h[0][0]==0) cout << "[ERROR] histo " << h[0][0]->GetName() << " not found!"<< endl;

  // ALL data (2nd dimention = eventClass = 0)  
  fillHistoNjet2D( h, 0, value, w );

  if(!IsData()) { //MC
    if(isTTbar)           fillHistoNjet2D( h, 1,  value, w );
    else if(isQCD) {      fillHistoNjet2D( h, 2,  value, w );
      if     (isEnri1)    fillHistoNjet2D( h, 3,  value, w );
      else if(isEnri2)    fillHistoNjet2D( h, 4,  value, w );
      else if(isEnri3)    fillHistoNjet2D( h, 5,  value, w );
      else if(isBce1)     fillHistoNjet2D( h, 6,  value, w );
      else if(isBce2)     fillHistoNjet2D( h, 7,  value, w );
      else if(isBce3)     fillHistoNjet2D( h, 8,  value, w );
    } else if(isWjets)    fillHistoNjet2D( h, 9,  value, w );
    else if(isZjets)      fillHistoNjet2D( h, 10, value, w );
    else if(isWjets)      fillHistoNjet2D( h, 11, value, w );
    else if(isSingleTop){ fillHistoNjet2D( h, 12, value, w );
      if     (isTW)       fillHistoNjet2D( h, 13, value, w );
      else if(isTchan)	  fillHistoNjet2D( h, 14, value, w );
      else if(isSchan)	  fillHistoNjet2D( h, 15, value, w );
    }
  }
  if(m_debug) cout << "End of << fillHisto_Njet_DataAndMC >>" << endl;
}
*/







//========================================================================================
//
//                                   histogram vectors
//
//========================================================================================

// // ------------ my method to add a set of 1D histograms acc to njet ----------------------
// void ana::addHistoNjet( vector<TH1*>& h, const string name, const string ext, const string title,
// 			const int nbin, const float xlow, const float xup ){
//   //  const string jetname[7]  = {"0j","1j","2j","3j","4j", "4mj", "allj"};
//   //  const string jetlabel[7] = {"0j","1j","2j","3j","4j", ">=4j","allj"};
//   h.resize(7);
//   for (unsigned int j=0; j<7; ++j) {
//     char hname[70];
//     char htitle[100];
//     sprintf( hname,  "%s_%s%s",  name.c_str(), jetname[j].c_str(), ext.c_str() );
//     sprintf( htitle, "%s (%s)",  title.c_str(), jetlabel[j].c_str() );
//     h[j] = new TH1F(hname, htitle, nbin, xlow, xup);
//     // Call Sumw2 to store sum of square of weights per bin
//     // Important to get histogram error to reflect statistical error of MC events, rather than sqrt(nexp)
//     h[j]->Sumw2();
//   }
// }

// // ------------ my method to add a set of 2D histograms acc to njet ----------------------
// void ana::addHistoNjet_TH2( vector<TH2*>& h, const string name, const string ext, const string title,
// 			    const int nbinx, const float xlow, const float xup,
// 			    const int nbiny, const float ylow, const float yup ){
//   if(m_debug) cout << "Starting << addHistoNjet_TH2 >> " << name << endl;
//   //  const string jetname[7]  = {"0j","1j","2j","3j","4j", "4mj", "allj"};
//   //  const string jetlabel[7] = {"0j","1j","2j","3j","4j", ">=4j","allj"};
//   /*
//   cout << "\n\n   addHistoNjet_TH2:  " << name << endl;
//   cout << "name=" << name << ",  title="<< title << endl;
//   cout << "ext=" << ext << endl;
//   cout << "amtb 0: h.resize(7)" << endl;
//   */
//   //  h.resize(7);//err
//   h.reserve(7);
//   //  cout << "amtb 1" << endl;

//   for (unsigned int j=0; j<7; ++j) {
//     //    cout << " amtb 2   j=" << j << endl;  
//     char hname[70];
//     char htitle[100];
//     sprintf( hname,  "%s_%s%s",  name.c_str(), jetname[j].c_str(), ext.c_str() );
//     sprintf( htitle, "%s (%s)",  title.c_str(), jetlabel[j].c_str() );
//     //cout << " amtb 3  " << hname << "  "<< htitle <<endl;
//     h[j] = new TH2F(hname, htitle, nbinx, xlow, xup, nbiny, ylow, yup);
//     h[j]->Sumw2();
//     //cout << " amtb after new, sumw2" << endl;
//   }
//   cout << "  end" << endl;
// }

// //-------------- my method to fill 1D histograms acc to njet -----------------------------
// void ana::fillHistoNjet(vector<TH1*>& h, const float& value, const double& w) {
  
//   h[6]->Fill(value, w); //allj
//   if( m_nGoodJet < 5 ) h[ m_nGoodJet ]->Fill(value, w); //=0,1,2,3,4j
//   if( m_nGoodJet > 3 ) h[5]->Fill(value, w); //4 or more jets
// }


// //-------------- my method to fill 2D histograms acc to njet -----------------------------
// void ana::fillHistoNjet(vector<TH2*>& h, const float& v1, const float& v2, const double& w) {
  
//   h[6]->Fill(v1, v2, w); //allj
//   if( m_nGoodJet < 5 ) h[ m_nGoodJet ]->Fill(v1, v2, w); //=0,1,2,3,4j
//   if( m_nGoodJet > 3 ) h[5]->Fill(v1, v2, w); //4 or more jets
// }


//--------- method to fill TH1D histograms acc to njet, GIVEN eventClass ------------------
void ana::fillHistoNjet2D(v2D_TH1& h, const int& ec, const float& value, const double& w ) {

  if(m_debug) cout << "fillHistoNjet2D"<< endl;
  h[6][ec]->Fill(value, w); //allj  
  if( m_nGoodJet < 5 )  h[ m_nGoodJet ][ec]->Fill(value, w); //0-4j
  if( m_nGoodJet > 3 )  h[5][ec]->Fill(value, w); //>=4j
}
//----------------------------------------------------------------------------------------


//----- method to add a set of TH1 histograms (for each type of MC when running on MC) ------
void ana::addHistoDataAndMC( vector<TH1*>& h, 
			     const string& name, const string& title,
			     const int& nbin, const float& xlow, const float& xup ) {
  h.reserve(ntype);

  for (int i=0; i<ntype; ++i) {
    // only do for mc present in input
    if( i>0 && is_mc_present(i)==false ) continue;
    char hname[70];
    char htitle[100];
    sprintf( hname,  "%s__%s",  name.c_str(),  mcname[i].c_str() );
    sprintf( htitle, "%s (%s)", title.c_str(), mclabel[i].c_str() );
    h[i] = new TH1F(hname, htitle, nbin, xlow, xup);
    h[i]->Sumw2();
  }
}
//----- method to add a set of TH2 histograms (for each type of MC when running on MC) ------
void ana::addHistoDataAndMC( vector<TH2*>& h, 
			     const string& name, const string& title,
			     const int& nbin, const float& xlow, const float& xup,
			     const int& nbiny, const float& ylow, const float& yup ) {
  h.reserve(ntype);

  for (int i=0; i<ntype; ++i) {
    // only do for mc present in input
    if( i>0 && is_mc_present(i)==false ) continue;
    char hname[70];
    char htitle[100];
    sprintf( hname,  "%s__%s",  name.c_str(),  mcname[i].c_str() );
    sprintf( htitle, "%s (%s)", title.c_str(), mclabel[i].c_str() );
    h[i] = new TH2F(hname, htitle, nbin, xlow, xup, nbiny, ylow, yup);
    h[i]->Sumw2();
  }
}


//-------- method to fill TH1 histograms (acc to MC type when running on MC) ---------
//void ana::fillHistoDataAndMC(vector<TH1*>& h, const float& value, const double& w ) {
void ana::fillHistoDataAndMC(vector<TH1*>& h, const float& value ) {

  if(m_debug) cout << "\nStarting << fillHistoDataAndMC >>: " << h[0]->GetName() << endl;
  h[0]->Fill(value, this_weight); //all events (data)
  if(!IsData()) { //run on MC
    if      (isTTbar)       h[1]->Fill( value, this_weight);
    else if (isQCD) {       h[2]->Fill( value, this_weight); //all QCD
      if      (isEnri1)     h[3]->Fill( value, this_weight);
      else if (isEnri2)     h[4]->Fill( value, this_weight);
      else if (isEnri3)     h[5]->Fill( value, this_weight);
      else if (isBce1)      h[6]->Fill( value, this_weight);
      else if (isBce2)      h[7]->Fill( value, this_weight);
      else if (isBce3)      h[8]->Fill( value, this_weight);
    } else if (isWjets)     h[9]->Fill( value, this_weight);
    else if (isZjets)       h[10]->Fill(value, this_weight);
    else if (isVQQ)         h[11]->Fill(value, this_weight);
    else if (isSingleTop) { h[12]->Fill(value, this_weight); //all single top
      if      (isTW)        h[13]->Fill(value, this_weight);
      else if (isTchan)     h[14]->Fill(value, this_weight);
      else if (isSchan)     h[15]->Fill(value, this_weight);
    }
  }
  if(m_debug) cout << "  End" << endl;
}
//-------- my method to fill TH2 histograms (acc to MC type when running on MC) ---------
void ana::fillHistoDataAndMC(vector<TH2*>& h, const float& v1, const float& v2, const double& w ) {
//void ana::fillHistoDataAndMC(vector<TH2*>& h, const float& v1, const float& v2, ) {

  if(m_debug) cout << "\nStarting << fillHistoDataAndMC 2D >>: " << h[0]->GetName() << endl;
  h[0]->Fill(v1,v2, w); //all events (data)
  if(!IsData()) { //run on MC
    if      (isTTbar)       h[1]->Fill(v1,v2,w);
    else if (isQCD) {       h[2]->Fill(v1,v2,w); //all QCD
      if      (isEnri1)     h[3]->Fill(v1,v2,w);
      else if (isEnri2)     h[4]->Fill(v1,v2,w);
      else if (isEnri3)     h[5]->Fill(v1,v2,w);
      else if (isBce1)      h[6]->Fill(v1,v2,w);
      else if (isBce2)      h[7]->Fill(v1,v2,w);
      else if (isBce3)      h[8]->Fill(v1,v2,w);
    } else if (isWjets)     h[9]->Fill(v1,v2,w);
    else if (isZjets)       h[10]->Fill(v1,v2,w);
    else if (isVQQ)         h[11]->Fill(v1,v2,w);
    else if (isSingleTop) { h[12]->Fill(v1,v2,w); //all single top
      if      (isTW)        h[13]->Fill(v1,v2,w);
      else if (isTchan)     h[14]->Fill(v1,v2,w);
      else if (isSchan)     h[15]->Fill(v1,v2,w);
    }
  }
  if(m_debug) cout << "  End" << endl;
}
//--------------------------------------------------------------------------------------


//----- my method to add a set of 1D histograms (for each type of MC when running on MC) ------
void ana::addHisto_Njet_DataAndMC( v2D_TH1& h, const string& name, const string& title,
				   const int& nbin, const float& xlow, const float& xup ) {
  // example: xxx_[1-4j]__[data,ttbar,wj..]

  if(m_debug) cout << "\nStaring << addHisto_Njet_DataAndMC >> vec:" << name << endl;

  h.resize(7);
  
  for (int j=0; j<7; ++j) {
    h[j].resize(ntype);
    for (int i=0; i<ntype; ++i) {
      // only do for mc present in input
      if( i>0 && is_mc_present(i)==false ) continue;
      char hname[70];
      char htitle[100];
      sprintf( hname,  "%s_%s__%s",    name.c_str(),  jetname[j].c_str(),  mcname[i].c_str() );
      sprintf( htitle, "%s (%s, %s)",  title.c_str(), jetlabel[j].c_str(), mclabel[i].c_str() );       
      h[j][i] = new TH1F(hname, htitle, nbin, xlow, xup);
      h[j][i]->Sumw2();      
    }
  }
}

//-------------- my method to fill  1D  histograms acc to njet & mctype ------------------
void ana::fillHisto_Njet_DataAndMC( v2D_TH1& h, const float& value, const double& w ) {

  if(m_debug) cout << "\nStarting << fillHisto_Njet_DataAndMC >> vec: " << h[0][0]->GetName() << endl;
  if(h[0][0]==0) cout << "[ERROR] histo " << h[0][0]->GetName() << " not found!"<< endl;

  // ALL data (2nd dimention = eventClass = 0)  
  fillHistoNjet2D( h, 0, value, w );

  if(!IsData()) { //MC
    if(isTTbar)           fillHistoNjet2D( h, 1,  value, w );
    else if(isQCD) {      fillHistoNjet2D( h, 2,  value, w );
      if     (isEnri1)    fillHistoNjet2D( h, 3,  value, w );
      else if(isEnri2)    fillHistoNjet2D( h, 4,  value, w );
      else if(isEnri3)    fillHistoNjet2D( h, 5,  value, w );
      else if(isBce1)     fillHistoNjet2D( h, 6,  value, w );
      else if(isBce2)     fillHistoNjet2D( h, 7,  value, w );
      else if(isBce3)     fillHistoNjet2D( h, 8,  value, w );
    } else if(isWjets)    fillHistoNjet2D( h, 9,  value, w );
    else if(isZjets)      fillHistoNjet2D( h, 10, value, w );
    else if(isWjets)      fillHistoNjet2D( h, 11, value, w );
    else if(isSingleTop){ fillHistoNjet2D( h, 12, value, w );
      if     (isTW)       fillHistoNjet2D( h, 13, value, w );
      else if(isTchan)	  fillHistoNjet2D( h, 14, value, w );
      else if(isSchan)	  fillHistoNjet2D( h, 15, value, w );
    }
  }
  if(m_debug) cout << "  End" << endl;
}


//
// book 3D histogram h[cut][nj][proc] for reliso
//                     11    7   16
// 1D
void ana::iso_addHisto_nlevel_nj_nmc(v3D_TH1& h,
				     const string& name,  const string& title,
				     const int& xnb,  const float& xlow,  const float& xup ){

  if(m_debug) cout << "Starting << iso_addisto_nlevel_nj_nmc >>: " << name << endl;
  //  cout << "amtb 1"<< endl;
  // cut
  const int nLevel=11;
  const string level_name[nLevel] = {"L1","L1b","L1c","L1d1","L1d2","L1d3","L1d4","L1d5","L2","L3","L4"};
  const string level_label[nLevel] = {"L1 HLT","L1b E_{T},#eta","L1c d0",
				      "L1d1 eID H/E",
				      "L1d2 eID |#Delta#sigma_{in}|",
				      "L1d3 eID |#Delta#phi_{in}|",
				      "L1d4 eID #sigma_{i#eta i#eta}",
				      "L1d5 eID",
				      "L2 GoodEle","L3 MuZVeto","L4 Conv"};
  h.resize(nLevel);

  for (int k=0; k<nLevel; ++k) {//nlevel

    h[k].resize(7);

    for (int j=0; j<7; ++j) {//nj

      h[k][j].resize(ntype);

      for (int i=0; i<ntype; ++i) {//nmc

	// only do for mc present in input
	if( i>0 && is_mc_present(i)==false ) continue;

	char hname[70];
	char htitle[100];
	sprintf( hname,  "%s_%s_%s__%s",  name.c_str(),  level_name[k].c_str(), jetname[j].c_str(),  mcname[i].c_str() );
	sprintf( htitle, "%s (%s,%s,%s)", title.c_str(), level_label[k].c_str(), jetlabel[j].c_str(), mclabel[i].c_str() );
	h[k][j][i] = new TH1F(hname, htitle, xnb, xlow, xup );
	h[k][j][i]->Sumw2();  
      }//i
    }//j
  }//k
  if(m_debug) cout << "  End" << endl;
}//end iso_addHisto_nlevel_nj_nmc

void ana::iso_addHisto_nlevel_nj_nmc(v3D_TH2& h,
				     const string& name,  const string& title,
				     const int& xnb,  const float& xlow,  const float& xup,
				     const int& ynb,  const float& ylow,  const float& yup){

  if(m_debug) cout << "Starting << iso_addisto_nlevel_nj_nmc >>: " << name << endl;
  //  cout << "amtb 1"<< endl;
  // cut
  const int nLevel=11;
  const string level_name[nLevel] = {"L1","L1b","L1c","L1d1","L1d2","L1d3","L1d4","L1d5","L2","L3","L4"};
  const string level_label[nLevel] = {"L1 HLT","L1b E_{T},#eta","L1c d0",
				      "L1d1 eID H/E",
				      "L1d2 eID |#Delta#sigma_{in}|",
				      "L1d3 eID |#Delta#phi_{in}|",
				      "L1d4 eID #sigma_{i#eta i#eta}",
				      "L1d5 eID",
				      "L2 GoodEle","L3 MuZVeto","L4 Conv"};
  // njet
  //  const string jetname[7]  = {"0j","1j","2j","3j","4j", "4mj", "allj"};
  //  const string jetlabel[7] = {"0j","1j","2j","3j","4j", ">=4j","allj"};
  // proc = mcname
  //  cout << "amtb 2"<< endl;
  h.resize(nLevel);
  //  cout << "amtb 3"<< endl;
  for (int k=0; k<nLevel; ++k) {//nlevel
    //cout << "amtb 4"<< endl;
    h[k].resize(7);
    //cout << "amtb 5"<< endl;
    for (int j=0; j<7; ++j) {//nj
      //cout << "amtb 6"<< endl;
      h[k][j].resize(ntype);
      ///      cout << "amtb 7"<< endl;
      for (int i=0; i<ntype; ++i) {//nmc
	//cout << "amtb 8"<< endl;
	// only do for mc present in input
	if( i>0 && is_mc_present(i)==false ) continue;
	//cout << "amtb 9"<< endl;
	char hname[70];
	char htitle[100];
	sprintf( hname,  "%s_%s_%s__%s",  name.c_str(),  level_name[k].c_str(), jetname[j].c_str(),  mcname[i].c_str() );
	sprintf( htitle, "%s (%s,%s,%s)", title.c_str(), level_label[k].c_str(), jetlabel[j].c_str(), mclabel[i].c_str() );
	h[k][j][i] = new TH2F(hname, htitle, xnb, xlow, xup, ynb, ylow, yup);
	//cout << "amtb 96"<< endl;
	//string hname = Form( "%s_%s_%s__%s",   name.c_str(),  level_name[k].c_str(), jetname[j].c_str(),  mcname[i].c_str() );
	//cout << "amtb 97  "<< hname << endl;
	//string htitle = Form( "%s (%s,%s,%s)",  title.c_str(), level_label[k].c_str(), jetlabel[j].c_str(), mclabel[i].c_str() );       
	//cout << "amtb 98  " << htitle << endl;
	//h[k][j][i] = new TH2F(Form("%s",hname.c_str()), Form("%s",htitle.c_str()), xnb, xlow, xup, ynb, ylow, yup);
	h[k][j][i]->Sumw2();
	//cout << "amtb 99"<< endl;
  
      }//i
    }//j
  }//k
  if(m_debug) cout << "  End" << endl;
}//end iso_addHisto_nlevel_nj_nmc


void ana::iso_fillHisto_NES( const int& ilevel, const bool& inBarrel,
			     const float& iso, const float& met ){
  if(m_debug) cout << "\n\n++ Starting << iso_fillHisto_NES >>" << endl;

  iso_fillHisto_nlevel_nj_nmc( h_QCDest_isoVmet_NES[ilevel],    iso, met, this_weight ); //NEW 1: L1b
  iso_fillHisto_nlevel_nj_nmc( h_QCDest_isoVmet_NES_uw[ilevel], iso, met, 1           ); //NB: unweighted
  iso_fillHisto_nlevel_nj_nmc( h_QCDest_iso_NES[ilevel],        iso );
  if(inBarrel) {
    iso_fillHisto_nlevel_nj_nmc( h_QCDest_isoVmet_NES_barrel[ilevel], iso, met, this_weight );
    iso_fillHisto_nlevel_nj_nmc( h_QCDest_iso_NES_barrel[ilevel],     iso );
  } else {
    iso_fillHisto_nlevel_nj_nmc( h_QCDest_isoVmet_NES_endcap[ilevel], iso, met, this_weight );
    iso_fillHisto_nlevel_nj_nmc( h_QCDest_iso_NES_endcap[ilevel],     iso );
  }
  if(m_debug) cout << "\n++ End of << iso_fillHist_NES >>" << endl;
}
//------------------------------------------------------------------------------------
// fill 2D vector of TH1 (weighted)
void ana::iso_fillHisto_nlevel_nj_nmc(v2D_TH1& h, const float& v1 ) {
  if(m_debug) cout << "\n   **** Starting << iso_fillHisto_nlevel_nj_nmc >>" << endl;

  vector<int> njbin;
  njbin.push_back(6); //fill [6]: allj
  if(m_nGoodJet<5)  njbin.push_back(m_nGoodJet); //[nj]:0j-4j
  if(m_nGoodJet>=4) njbin.push_back(5); //[5]:>=4j

  //  cout << "njbin.size: "<< njbin.size() << endl;
  //  cout << "goodjet: "<< m_nGoodJet << endl;

  for(unsigned int j=0; j<njbin.size(); ++j){//nj
    int thejetbin = njbin.at(j);
    fillHistoDataAndMC( h[thejetbin], v1 ); //h[j] is a vector of TH1
  }

  if(m_debug) cout << "\n   **** End of << iso_fillHisto_nlevel_nj_nmc >>" << endl;
}//end iso_addHisto_nlevel_nj_nmc
//------------------------------------------------------------------------------------


// fill 2D vector of TH2 (can be weighted or unweighted)
void ana::iso_fillHisto_nlevel_nj_nmc(v2D_TH2& h, const float& v1, const float& v2, const double& w) {
  if(m_debug) cout << "\n   +++ Starting << iso_fillHisto_nlevel_nj_nmc >>" << endl;

  vector<int> njbin;
  njbin.push_back(6); //fill [6]: allj
  if(m_nGoodJet<5) njbin.push_back(m_nGoodJet); //[nj]:0j-4j
  if(m_nGoodJet>=4) njbin.push_back(5); //[5]:>=4j

  //  cout << "njbin.size: "<< njbin.size() << endl;
  //  cout << "goodjet: "<< m_nGoodJet << endl;
  if(m_debug) cout << "will call fillHistoDataAandMC"<< endl;

  for(unsigned int j=0; j<njbin.size(); ++j){//nj
    int thejetbin = njbin.at(j);
    fillHistoDataAndMC( h[thejetbin], v1, v2, w ); //h[j] is a vector of TH2
  }

  if(m_debug) cout << "\n   +++ End of << iso_fillHisto_nlevel_nj_nmc >>" << endl;
}//end iso_addHisto_nlevel_nj_nmc
//------------------------------------------------------------------------------------






//----------------------------------------------------------------------------------------
//81FB 
//void ana::valid_mkHisto_cut_njet(TH1F* h[][7], const string name, const string title,
//                                 const int nbin, const float xlow, const float xup ){
void ana::valid_mkHisto_cut_njet(v2D_TH1& h, const string& name, const string& title,
                                 const int& nbin, const float& xlow, const float& xup ){

  const string cutname[9]  = {"Initial", "Trigger", ">=1T ele",">=1T Iso","==1T Iso", "noMuon", "MET", "Z veto", "noConv"};
  const string cutlabel[9] = {"Initial", "Trigger", ">=1Tele",">=1TIso","==1TIso", "noMuon", "MET", "Zveto", "noConv"};

  h.resize(9);
  for(int i=0; i<9; ++i){ //cut
    h[i].resize(7);
    for(int j=0; j<7; ++j){ //nj
      char hname[70];
      char htitle[100];
      sprintf( hname,  "%s_%s_%s",    name.c_str(),  cutlabel[i].c_str(), jetname[j].c_str() );
      sprintf( htitle, "%s %s (%s)",  title.c_str(), cutname[i].c_str(),  jetlabel[j].c_str() );
      h[i][j] = new TH1F(hname, htitle, nbin, xlow, xup);      
    }
  }
}
//----------------------------------------------------------------------------------------


//void ana::valid_fillHisto(TH1F* h[][7], const bool cuts[8], int nj, double value) const {
void ana::valid_fillHisto(v2D_TH1& h, const bool cuts[8], const double& value) const {

  for(int i=0; i<9; ++i){ //cut
    if(cuts[i]) {
      h[i][6]->Fill(value);//all jets
      if(m_nGoodJet < 5){ h[i][m_nGoodJet]->Fill(value); }//fill 0-4j
      if(m_nGoodJet > 3){ h[i][5]->Fill(value); }//fill >=4j
    }
    else{break;}
  }
}
//81FB end
//----------------------------------------------------------------------------------------

//=========================================================================================
//
//                                   end of histo vectors
//
//=========================================================================================




//----------------------------------------------------------------------------------------

void ana::fillHisto_event_tables() {

   for(int i=0; i<nstage; ++i) {
     for(int j=0; j<ntjet; ++j) {

       // *** All events ***
       Data_njetsVcuts->Fill(j,nstage-i-1, e_plus_jet_weighted[i][j][0]); //NEW, 3rd dimension is [0] for data.
       Data_njetsVcuts->Fill(5,nstage-i-1, e_plus_jet_weighted[i][j][0]); //sum for all jets

       // *** MC ***
       if(!IsData()) {
	 if ( mc_sample_has_ttbar ) {//MC contains Signal
	   for(int k=1; k<11; ++k){
	     Signal_njetsVcuts->Fill(j,nstage-i-1, e_plus_jet_weighted[i][j][k]);
	     Signal_njetsVcuts->Fill(5,nstage-i-1, e_plus_jet_weighted[i][j][k]);//sum for all jets
	   }	   
	   if(signal_is_Alpgen) {
	     for(int k=23; k<=27; ++k) {//mctype 23-27
               Signal_njetsVcuts->Fill(j,nstage-i-1, e_plus_jet_weighted[i][j][k]);
               Signal_njetsVcuts->Fill(5,nstage-i-1, e_plus_jet_weighted[i][j][k]);//sum for all jets
	     }
	   }//ttnj alpgen
	   
	 }
	 if ( mc_sample_has_QCD ) {//MC contains QCD
	   for(int k=13; k<19; ++k){ // mctype is 13 to 18 for QCD
	     QCD_njetsVcuts->Fill(j,nstage-i-1, e_plus_jet_weighted[i][j][k]);
	     QCD_njetsVcuts->Fill(5,nstage-i-1, e_plus_jet_weighted[i][j][k]);//sum for all jets
	   }
	 }
	 if ( mc_sample_has_Wjet ) {//MC contains Wjets (mctype=11)
	   Wjets_njetsVcuts->Fill(j,nstage-i-1, e_plus_jet_weighted[i][j][11]);
	   Wjets_njetsVcuts->Fill(5,nstage-i-1, e_plus_jet_weighted[i][j][11]);//sum for all jets
	 }
	 if ( mc_sample_has_Zjet ) {//MC contains Zjets (mctype=12)
	   Zjets_njetsVcuts->Fill(j,nstage-i-1, e_plus_jet_weighted[i][j][12]);
	   Zjets_njetsVcuts->Fill(5,nstage-i-1, e_plus_jet_weighted[i][j][12]);//sum for all jets
	 }
	 if ( mc_sample_has_VQQ ) {//MC contains VQQ (mctype=19)
	   VQQ_njetsVcuts->Fill(j,nstage-i-1, e_plus_jet_weighted[i][j][19]);
	   VQQ_njetsVcuts->Fill(5,nstage-i-1, e_plus_jet_weighted[i][j][19]);//sum for all jets
	 }
	 if ( mc_sample_has_singleTop ) {//MC contains single top (mctypes are 20-22)
	   for(int k=20; k<23; ++k){
	     SingleTop_njetsVcuts->Fill(j,nstage-i-1, e_plus_jet_weighted[i][j][k]);
	     SingleTop_njetsVcuts->Fill(5,nstage-i-1, e_plus_jet_weighted[i][j][k]);//sum for all jets
	   }
	 }
       }//end MC
     }
   }

   // Set labels of the histograms (cuts and Njet)
   SetHistoLabelCutNjet( Signal_njetsVcuts,     ve );
   SetHistoLabelCutNjet( QCD_njetsVcuts,        ve );
   SetHistoLabelCutNjet( Wjets_njetsVcuts,      ve );
   SetHistoLabelCutNjet( Zjets_njetsVcuts,      ve );
   SetHistoLabelCutNjet( VQQ_njetsVcuts,        ve );
   SetHistoLabelCutNjet( SingleTop_njetsVcuts,  ve );
   SetHistoLabelCutNjet( Data_njetsVcuts,       ve );

}//end fillHisto_event_tables
//----------------------------------------------------------------------------------------


void ana::fillHisto_PDF_weights( TH1F* h ){
  // 45 bins
  h->Fill(-1,  1);

  h->Fill(0.,  PDFWcteq66_0);
  h->Fill(1.,  PDFWcteq66_1);
  h->Fill(2.,  PDFWcteq66_2);
  h->Fill(3.,  PDFWcteq66_3);
  h->Fill(4.,  PDFWcteq66_4);
  h->Fill(5.,  PDFWcteq66_5);
  h->Fill(6.,  PDFWcteq66_6);
  h->Fill(7.,  PDFWcteq66_8);
  h->Fill(8.,  PDFWcteq66_8);
  h->Fill(9.,  PDFWcteq66_9);
  h->Fill(10., PDFWcteq66_10);

  h->Fill(11., PDFWcteq66_11);
  h->Fill(12., PDFWcteq66_12);
  h->Fill(13., PDFWcteq66_13);
  h->Fill(14., PDFWcteq66_14);
  h->Fill(15., PDFWcteq66_15);
  h->Fill(16., PDFWcteq66_16);
  h->Fill(17., PDFWcteq66_18);
  h->Fill(18., PDFWcteq66_18);
  h->Fill(19., PDFWcteq66_19);
  h->Fill(20., PDFWcteq66_20);

  h->Fill(21., PDFWcteq66_21);
  h->Fill(22., PDFWcteq66_22);
  h->Fill(23., PDFWcteq66_23);
  h->Fill(24., PDFWcteq66_24);
  h->Fill(25., PDFWcteq66_25);
  h->Fill(26., PDFWcteq66_26);
  h->Fill(27., PDFWcteq66_28);
  h->Fill(28., PDFWcteq66_28);
  h->Fill(29., PDFWcteq66_29);
  h->Fill(30., PDFWcteq66_30);

  h->Fill(31., PDFWcteq66_31);
  h->Fill(32., PDFWcteq66_32);
  h->Fill(33., PDFWcteq66_33);
  h->Fill(34., PDFWcteq66_34);
  h->Fill(35., PDFWcteq66_35);
  h->Fill(36., PDFWcteq66_36);
  h->Fill(37., PDFWcteq66_38);
  h->Fill(38., PDFWcteq66_38);
  h->Fill(39., PDFWcteq66_39);
  h->Fill(40., PDFWcteq66_40);

  h->Fill(41., PDFWcteq66_41);
  h->Fill(42., PDFWcteq66_42);
  h->Fill(43., PDFWcteq66_43);
  h->Fill(44., PDFWcteq66_44);

}//end fillHisto_PDF_weights
//----------------------------------------------------------------------------------------

//========================================================================================
//
//                            End of Historam Manipulation
//
//========================================================================================






//========================================================================================
//
//                                  Conversion Finder
//
//========================================================================================
// This was originally ConversionFinder but has been replaced (07-01-10). It is left here until the new routine has been thoroughly tested. 
bool ana::ConversionFinder2(const TLorentzVector& e1, int mctype, int index_selected_ele)  {

  if(m_debug) cout << "Starting << ConversionFinder2 >>" << endl;
  bool isthisConversion = false;

  const float bfield=3.8;
  
  float phi_ie = e1.Phi();
  float eta_ie = e1.PseudoRapidity();

  int mytrackref  = static_cast<int>( els_closestCtfTrackRef->at(index_selected_ele) );
  /*
  float tphi1;
  float teta1;
  if(mytrackref>-1){
    tphi1 = tracks_phi->at(mytrackref);
    teta1 = tracks_eta->at(mytrackref);
    float tEleTrackDelR = calcDeltaR( tphi1, teta1, phi_ie, eta_ie);
    ///    cout<<"tEleTrackDelR: "<<tEleTrackDelR <<endl;
    Conv_CheckDelR_GSFTk_ctfTk->Fill(tEleTrackDelR);
  }
  else{Conv_CheckDelR_GSFTk_ctfTk->Fill(10);}
  */


  //declare track variables
  float phi1,eta1,phi2,eta2;
  float tk1Curvature;
  float tk2Curvature;
  float tk1r;
  float tk2r;
  float d0;
  float d02;
  float tk1x;
  float tk1y;
  float tk2x;
  float tk2y;
  float distmag;
  float dist;
  float dcot;
  
 
  //first loop over tracks, trying to match to the electron
  for (unsigned int i=0; i<Ntracks; ++i){
    if(mytrackref > -1){
      i = mytrackref;
    }

    phi1 = tracks_phi->at(i);
    eta1 = tracks_eta->at(i);
    
    float EleTrackDelR = calcDeltaR( phi1, eta1, phi_ie, eta_ie);
    
    //consider only tracks that fall within a small cone around electron
    if(EleTrackDelR > 0.3 && mytrackref < 0) continue;
    
    //calculate curvature and radius of track
    tk1Curvature = -0.3*bfield*tracks_chg->at(i)/( 100*( tracks_pt->at(i) ) );
    tk1r = fabs(1/tk1Curvature);

    /*
    //  d0 = tracks_d0->at(i);
    // Calculate d0 w.r.t beam spot position (x,y)=(0.0325839, 0.0001793) (3 Mar 09)
    float vx = tracks_vx->at(i);
    float vy = tracks_vy->at(i);
    float px = tracks_px->at(i);
    float py = tracks_py->at(i);
    float pt = tracks_pt->at(i);
    d0 = - (-(vx-0.0325839)*py + (vy-0.0001793)*px) / pt ; //d0 = -dxy
    */
    // NEW (Aug 09)
    d0 = compute_d0("track",i);



    //calculate coordinates of centre of track "circle"
    tk1x = ((1/tk1Curvature) - d0)*cos(phi1);
    tk1y = ((1/tk1Curvature) - d0)*sin(phi1);
    //loop over the tracks again, try to find other track in a conversion
    for (unsigned int j=0; j<Ntracks; ++j){
      //avoid using the same track:
      if(i==j) continue;
      phi2 = tracks_phi->at(j);
      eta2 = tracks_eta->at(j);

      //For a restriction on the position of the second track with respect to the
      //electron, uncomment the line following. The main advantage is a slight
      //reduction in signal loss with very little change in the conversion
      //rejection rate
      //if( calcDeltaR( phi2, eta2, phi_ie, eta_ie) > 0.5 ) continue;

      //check the tracks have opposite charges
      if(tracks_chg->at(i)*tracks_chg->at(j) > 0.0) continue;


      tk2Curvature = -0.3*bfield*tracks_chg->at(j)/( 100*( tracks_pt->at(j) ) );
      tk2r = fabs(1/tk2Curvature);


      //d02 = tracks_d0->at(j);
      // Calculate d0 w.r.t beam spot position (x,y)=(0.0325839, 0.0001793) (3 Mar 09)
//       float vx2 = tracks_vx->at(j);
//       float vy2 = tracks_vy->at(j);
//       float px2 = tracks_px->at(j);
//       float py2 = tracks_py->at(j);
//       float pt2 = tracks_pt->at(j);
//       d02 = - (-(vx2-0.0322)*py2 + vy2*px2) / pt2 ; //d0 = -dxy
      d02 = compute_d0("track",j);


      tk2x = ((1/tk2Curvature) - d02)*cos(phi2);
      tk2y = ((1/tk2Curvature) - d02)*sin(phi2);

      //calculate distance between curves
      distmag = (tk2x - tk1x)*(tk2x -tk1x) + (tk2y -tk1y)*(tk2y-tk1y);
      distmag = sqrt(distmag);
      dist = distmag-(tk1r+tk2r);
      dcot = 1/tan(tracks_theta->at(i)) - 1/tan(tracks_theta->at(j));

      // FB/TL (13 Feb 09): update cut values for better performace
      // studies shown that using the new criteria, acceptance of conversion
      // events increase from 33% to ~50%, with littler effect on the signal 
      // rejection (1.4% to 2%)
      // if( dist > -0.02 && dist < 0.01 && fabs(dcot) < 0.02)
      if( dist > -0.04 && dist < 0.04 && fabs(dcot) < 0.03)
	{
	  //cout <<"Found a Conversion!"<<endl;
	  isthisConversion = true;
	  break;
	}

    }//end of second track loop

    if(isthisConversion) break;
    if(mytrackref > -1) break;
  }//end of first track loop

  //856
  if( DoConversionStudies() ) {

    if ( Nmc_doc==0 ){
    // print only first 20 occurances
    static int mm = 0;
    if(mm==0) cout << "(Printing first 20 occurances.)" << endl;
    if(mm<20) { cout << "PROBLEM: Nmc_doc = 0. Cannot do conversion studies (involve MC matching)." << endl; ++mm; }
    }
    else if ( Nmc_doc>0 ){

      int didConv = 0;
      if(isthisConversion){
	didConv = 1;
      }

      //  int mctype = 1;
      ConversionArray[mctype][didConv][0]++;
      
      float mc_phi;
      float mc_eta;
      int ii=0;
      float tempDelR = 100;
      
      
      for(unsigned int i=0;i<Nmc_doc;++i){
	mc_phi = mc_doc_phi->at(i);
	mc_eta = mc_doc_eta->at(i);
	float empDelR = calcDeltaR( mc_phi, mc_eta, phi_ie, eta_ie);
	if(empDelR > 0.3) continue;
	
	if(empDelR < tempDelR){tempDelR = empDelR;ii=i;}
      }//end mc particle loop  
      
      //     cout << "amtb 4c" << endl;
      //     cout << "ii="<<ii << endl;
      //     cout << "Nmc_doc="<< Nmc_doc << endl;
      //     cout << "mc_doc_id(ii)="<< mc_doc_id->at(ii) << endl;
      
      if( fabs( mc_doc_id->at(ii) ) == 11){
	if(m_debug)  cout << "It matches closest to a real electron" << endl;  
	ConversionArray[mctype][didConv][1]++;
      }
      
      else if( fabs( mc_doc_id->at(ii) ) == 22){
	if(m_debug) cout << "It matches closest to a photon" << endl;
	ConversionArray[mctype][didConv][2]++;
      }
      
      else if( fabs( mc_doc_id->at(ii) ) == 111 || fabs( mc_doc_id->at(ii) ) == 221  ){
	if(m_debug) cout << "It matches closest to a pi zero or eta" << endl;
	ConversionArray[mctype][didConv][3]++;
      }
      else{ 
	if(m_debug) cout << "It matches closest to something else" << endl;
	ConversionArray[mctype][didConv][4]++;
      }  
      
      //end 856
    }//has genPar
  }//if DoConversionStudies

  return isthisConversion;

}//end ConversionFinder2


//return generated particle match to reco object based on delR
float ana::MCTruthMatch(float eta, float phi){
  float mc_phi,mc_eta,tempDelR,finalDelR=100;
  int ii = 0;
  for(unsigned int i=0;i<Nmc_doc;++i){
    mc_phi = mc_doc_phi->at(i);
    mc_eta = mc_doc_eta->at(i);
    tempDelR = calcDeltaR( mc_phi, mc_eta, phi, eta);
    if(tempDelR > 0.3) continue;

    if(tempDelR < finalDelR){finalDelR = tempDelR;ii=i;}
  }//end mc particle loop                               
  if(ii==-1) cout << "ERROR: index for mc_doc_id is -1"<< endl;
  return mc_doc_id->at(ii);
}





bool ana::ConversionFinder(const TLorentzVector& e1, int mctype, int index_selected_ele) {

  //Should do no more than return whether the electron passes the conversion algo or not
  //This currently (07-01-10) doesn't actually use mctype or e1. For now leave in case they might be useful for 
  //later alterations 


  if(m_debug) cout << "Starting << ConversionFinder >>" << endl;
  bool isthisConversion = false;

  //if the electron has 3 or more missing layers, declare it a conversion
  if(m_useMisslayers) { if( els_innerLayerMissingHits->at(index_selected_ele) > 2 ) return true; }
  //else, continue the routine. Note, we may want to apply this separately to the routine. For now, leave it here (07-01-10)


  const float bfield=3.8;

  int mytrackref  = static_cast<int>( els_closestCtfTrackRef->at(index_selected_ele) );

  float phi1,eta1,phi2,eta2;
  float tk1Curvature,tk2Curvature;
  float tk1r,tk2r;
  float d0,d02;
  float tk1x,tk1y,tk2x,tk2y;
  float distmag,dist,dcot;
  float EleTrackDelR,theta1,charge1;
  float NumMissLayers = -1;

  //list of cut values subject to possible optimisation:
  // 1) EleTrackDelR > 0.3
  // 2) Dist
  // 3) Dcot
  // 4) Bool opposite track charges - required or not.
  // 5) EleTrackDelR between both tracks < 0.3
 

  //  mytrackref = 0;
  //if track ref is valid, AND the ctf track and gsf track have a minimum required number of matching hits
  // then use this CTF track (remember to skip in loop over track collection). Otherwise we use the GSF track 
  if(mytrackref>=0 && els_shFracInnerHits->at(index_selected_ele) > 0.45){//use CTF params
    phi1 = tracks_phi->at(mytrackref);
    eta1 = tracks_eta->at(mytrackref);
    charge1 = tracks_chg->at(mytrackref);
    d0 = compute_d0("track",mytrackref);
    tk1Curvature = -0.3*bfield*charge1/( 100*( tracks_pt->at(mytrackref) ) );
    theta1 = tracks_theta->at(mytrackref);
    if(m_useMisslayers) NumMissLayers = tracks_innerLayerMissingHits->at(mytrackref); 
  }
  else{//use GSF params
    phi1 = els_tk_phi->at(index_selected_ele);
    eta1 = els_tk_eta->at(index_selected_ele);
    charge1 = els_tk_charge->at(index_selected_ele);
    d0 = compute_d0("electron",index_selected_ele);
    tk1Curvature = -0.3*bfield*charge1/( 100*( els_tk_pt->at(index_selected_ele) ) );
    theta1 = els_tk_theta->at(index_selected_ele);
    if(m_useMisslayers) NumMissLayers = els_innerLayerMissingHits->at(index_selected_ele);
    }
   
  tk1r = fabs(1/tk1Curvature);
  tk1x = ((1/tk1Curvature) - d0)*cos(phi1);
  tk1y = ((1/tk1Curvature) - d0)*sin(phi1);
  

  for (unsigned int j=0; j<Ntracks; ++j){//loop over the track collection 
    
    if( (abs(mytrackref)-j)==0 ) continue;//works if trackref is valid or not (invalid it will be less than zero, which j is not)
    phi2 = tracks_phi->at(j);
    eta2 = tracks_eta->at(j);
 
    //require difference in gsf and ctf charges to be different   
    if(charge1*tracks_chg->at(j) > 0.0) continue;

    EleTrackDelR = calcDeltaR(phi1, eta1, phi2, eta2);
    if(EleTrackDelR > 0.3) continue;
    
    tk2Curvature = -0.3*bfield*tracks_chg->at(j)/( 100*( tracks_pt->at(j) ) );
    tk2r = fabs(1/tk2Curvature);
    
    d02 = compute_d0("track",j);
    
    
    tk2x = ((1/tk2Curvature) - d02)*cos(phi2);
    tk2y = ((1/tk2Curvature) - d02)*sin(phi2);
    
    distmag = (tk2x - tk1x)*(tk2x -tk1x) + (tk2y -tk1y)*(tk2y-tk1y);
    distmag = sqrt(distmag);
    dist = distmag-(tk1r+tk2r);
    dcot = 1/tan(theta1) - 1/tan(tracks_theta->at(j));
    
    
    if( dist > -0.04 && dist < 0.04 && fabs(dcot) < 0.03)
      {
	bool mlPass =true;
	//require similar number of missing layers for each track
	if(m_useMisslayers) {if( fabs(NumMissLayers - tracks_innerLayerMissingHits->at(j) ) >2  ){mlPass=false;}}
	if(mlPass){
	  isthisConversion = true;
	  break;
	}
      }
    
    if(isthisConversion) break;
  }//end of second track loop
  
  if(m_debug) cout << "Ending << ConversionFinder >>" << endl;
  return isthisConversion;
}//end ConversionFinder

void ana::ConversionMCMatching(const TLorentzVector& e1, int mctype, bool isthisConversion){

  //float phi_ie,eta_ie;
  float phi_ie = e1.Phi();
  float eta_ie = e1.PseudoRapidity();

  int didConv = 0;
  if(isthisConversion){
    didConv = 1;
  }

  //totals the number of electrons flagged or not flagged as conversions
  ConversionArray[mctype][didConv][0]++;

  float mc_phi,mc_eta,empDelR;
  int ii=0;
  float tempDelR = 100;

  if( Nmc_doc>0 ){//only run if there are MC particles in the file 

    for(unsigned int i=0;i<Nmc_doc;++i){
      mc_phi = mc_doc_phi->at(i);
      mc_eta = mc_doc_eta->at(i);
      empDelR = calcDeltaR( mc_phi, mc_eta, phi_ie, eta_ie);
      if(empDelR > 0.3) continue;
      
      if(empDelR < tempDelR){tempDelR = empDelR;ii=i;}
    }//end mc particle loop   
    
    if( fabs( mc_doc_id->at(ii) ) == 11){
      if(m_debug)  cout << "It matches closest to a real electron" << endl;
      ConversionArray[mctype][didConv][1]++;
    }
    
    else if( fabs( mc_doc_id->at(ii) ) == 22){
      if(m_debug) cout << "It matches closest to a photon" << endl;
      ConversionArray[mctype][didConv][2]++;
    }
    
    else if( fabs( mc_doc_id->at(ii) ) == 111 || fabs( mc_doc_id->at(ii) ) == 221  ){
      if(m_debug) cout << "It matches closest to a pi zero or eta" << endl;
      ConversionArray[mctype][didConv][3]++;
    }
    else if( fabs( mc_doc_id->at(ii) ) == 211 ){
      if(m_debug) cout << "It matches closest to charged pion" << endl;
      ConversionArray[mctype][didConv][4]++;
    }
    else{
      if(m_debug) cout << "It matches closest to something else" << endl;
      ConversionArray[mctype][didConv][5]++;
    }
    
  }

  return;

}







// 06-01-10: Should remove this or at least change it - the algo is obsolete. 
//856
void ana::OptimiseConversionFinder(const TLorentzVector& e1, int mctype){

  //  bool isthisConversion = false;

  const float bfield=3.8;

  float phi_ie = e1.Phi();
  float eta_ie = e1.PseudoRapidity();

  float mc_phi;
  float mc_eta;
  int ii=0;
  float tempDelR = 100;

  for(unsigned int i=0;i<Nmc_doc;++i){
    mc_phi = mc_doc_phi->at(i);
    mc_eta = mc_doc_eta->at(i);
    float empDelR = calcDeltaR( mc_phi, mc_eta, phi_ie, eta_ie);
    if(empDelR > 0.3) continue;

    if(empDelR < tempDelR){tempDelR = empDelR;ii=i;}
  }//end mc particle loop                                                          
  
  //if particle doesn't match a photon and isn't ttbar, return. For ttbar, want all 
  if( fabs( mc_doc_id->at(ii) ) != 22 && fabs( mc_doc_id->at(ii) ) != 111 && mctype >10){return;}


  //declare track variables     
  float phi1;
  float eta1;
  float phi2;
  float eta2;
  float tk1Curvature;
  float tk2Curvature;
  float tk1r;
  float tk2r;
  float d0;
  float d02;
  float tk1x;
  float tk1y;
  float tk2x;
  float tk2y;
  float distmag;
  float dist;
  float dcot;

  float Mindist = 100;
  float Mindcot = 100;
  float min_didc = 1000;
 
  //first loop over tracks, trying to match to the electron                      
  for (unsigned int i=0; i<Ntracks; ++i){

    phi1 = tracks_phi->at(i);
    eta1 = tracks_eta->at(i);

    float EleTrackDelR = calcDeltaR( phi1, eta1, phi_ie, eta_ie);

    //consider only tracks that fall within a small cone around electron         
    if(EleTrackDelR > 0.3) continue;

    //calculate curvature and radius of track                                    
    tk1Curvature = -0.3*bfield*tracks_chg->at(i)/( 100*( tracks_pt->at(i) ) );
    tk1r = fabs(1/tk1Curvature);

    d0 = compute_d0("track",i);

    //calculate coordinates of centre of track "circle"                            
    tk1x = ((1/tk1Curvature) - d0)*cos(phi1);
    tk1y = ((1/tk1Curvature) - d0)*sin(phi1);

    //loop over the tracks again, try to find other track in a conversion          
    for (unsigned int j=0; j<Ntracks; ++j){
      //avoid using the same track:                                                
      if(i==j) continue;

      phi2 = tracks_phi->at(j);
      eta2 = tracks_eta->at(j);

      //check the tracks have opposite charges                                     
      if(tracks_chg->at(i)*tracks_chg->at(j) > 0.0) continue;


      tk2Curvature = -0.3*bfield*tracks_chg->at(j)/( 100*( tracks_pt->at(j) ) );
      tk2r = fabs(1/tk2Curvature);

      d02 = compute_d0("track",j);


      tk2x = ((1/tk2Curvature) - d02)*cos(phi2);
      tk2y = ((1/tk2Curvature) - d02)*sin(phi2);

      //calculate distance between curves                                      
      distmag = (tk2x - tk1x)*(tk2x -tk1x) + (tk2y -tk1y)*(tk2y-tk1y);
      distmag = sqrt(distmag);
      dist = distmag-(tk1r+tk2r);
      dcot = 1/tan(tracks_theta->at(i)) - 1/tan(tracks_theta->at(j));

      //float temp_mindidc = sqrt(dist*dist/(3.2*3.2) + dcot*dcot/(2.7*2.7));
      float temp_mindidc = sqrt(dist*dist + dcot*dcot);

      if(temp_mindidc < min_didc){min_didc = temp_mindidc;Mindist = dist; Mindcot = dcot;}

      //if(mctype<11){ Conv_Opti[0]->Fill(dist,dcot);Conv_Optis[0]->Fill(dist,dcot);}
      //if(mctype==15){Conv_Opti[1]->Fill(dist,dcot);Conv_Optis[1]->Fill(dist,dcot);}

      // if( dist > -0.04 && dist < 0.04 && fabs(dcot) < 0.03)
      //{
          //cout <<"Found a Conversion!"<< endl;                      
      //  isthisConversion = true;
	  //    break;
      //}

    }//end of second track loop                                                
    //   if(isthisConversion) break;
  }//end of first track loop                                                   

  //  return isthisConversion;
  if(mctype<11){ Conv_Opti[0]->Fill(Mindist,Mindcot);Conv_Optis[0]->Fill(Mindist,Mindcot);Conv_Opti_extragran[0]->Fill(Mindist,Mindcot);}
  if(mctype==15){Conv_Opti[1]->Fill(Mindist,Mindcot);Conv_Optis[1]->Fill(Mindist,Mindcot);Conv_Opti_extragran[1]->Fill(Mindist,Mindcot);}

}






void ana::PrintConversionTable(){

  TString MySamples[14] = {"ttbar","W+jet","Z+Jet","Enri1","Enri2","Enri3","bce1","bce2","bce3","vqq","tW","tchan","schan","data"};
  TString ConvNames[13] = {"&  Iso Ele ", "&  Convers ","& Electron ","&   Photon ","&   PiZero ","&  ChrPion ","&    Other ","& Non Conv ","& Electron ","&   Photon ","&   PiZero ","&  ChrPion ","&    Other "};

  cout << endl << endl << "-------------------------------------------------------------------";
  cout << endl << "-------------------------------------------------------------------";
  cout << endl << "----------------Conversion Breakdown-------------------------------";
  cout << endl << "-------------------------------------------------------------------";
  cout << endl;
  cout << "Shows Conversions with matching particle, and non conversions with matching particle, for all isolated reco electrons"<<
    endl<<endl;
  
  cout <<"mctype ";
  for(int i=0;i<12;++i) cout<<ConvNames[i];
  cout<<ConvNames[12]<<" \\\\"  <<endl;

  //for(int k=0;k<20;k++){                                                                                                              
  int ff;
  
  for(int k=10;k!=1;++k){
    if(k==11||k==12) continue;
    if(k==10){
      for(int kkc=0;kkc<2;++kkc){
	for(int j=0;j<6;++j){
	  for(int kk=1;kk<10;++kk){ConversionArray[10][kkc][j]+=ConversionArray[kk][kkc][j];
	  }
	}
      }
    }
    if(k==0){ff=13;}
    else{ff = k-10;}
    cout<< setw(7) <<MySamples[ff] <<"&  "<<setw(7)<<(ConversionArray[k][0][0]+ConversionArray[k][1][0]);
    
    for(int i=1;i>-1;--i){
      //     cout<<ConversionArray[0][i][0]<<setw(8);                                                                         
      for(int j=0;j<6;++j){
	cout<<" & "<<setw(8)<<ConversionArray[k][i][j];
      }
    }
    cout << "  \\\\" << endl;
    if(k==22) k=-1;
  }
  cout<<endl;
  
}//end PrintConversionTable()
// end 856
//====================================================================================
//
//                           End  of  Conversion  Finder
//
//====================================================================================





// ----------- my method to compute hadronic top mass (m3) --------------
// m3 = invariant mass of 3 jets which have highest vector sum PT
//
void ana::reco_hadronicTop_highestTopPT(const vector<TLorentzVector>& jetColl, const int nGoodIsoEle ){

  pair<double,double> res = compute_M3(jetColl);
  const double max_top_PT = res.second;
  //double had_top_mass = res.first;
  const double m3 = res.first; //renamed from "had_top_mass"
  //cout << "selected Top (max pt):  pt=" << max_top_PT << "  m=" << had_top_mass << endl;
  
  // fill histo
  if ( nGoodIsoEle > 0 ) { // this is signal region
    
    if (jetColl.size()>=4) {
      h_hadTop_maxPT_pt_4j->Fill(max_top_PT, this_weight);
      h_hadTop_maxPT_mass_4j->Fill(m3, this_weight);
      h_hadTop_maxPT_mass_4j_1000->Fill(m3, this_weight);//15Jan10

      if (IsData()==false) {  // fill according to type of MC
	if (isTTbar)       {  h_m3_tt->Fill(     m3, this_weight);  h_m3_tt_1000->Fill(     m3, this_weight);  }
	else if (isWjets)  {  h_m3_wj->Fill(     m3, this_weight);  h_m3_wj_1000->Fill(     m3, this_weight);  }
	else if (isZjets)  {  h_m3_zj->Fill(     m3, this_weight);  h_m3_zj_1000->Fill(     m3, this_weight);  }
	else if (isQCD)  {    h_m3_qcd->Fill(    m3, this_weight);  h_m3_qcd_1000->Fill(    m3, this_weight); 	
	  if     (isEnri3) {  h_m3_enri[2]->Fill(m3, this_weight);  h_m3_enri_1000[2]->Fill(m3, this_weight);  }
	  else if(isEnri2) {  h_m3_enri[1]->Fill(m3, this_weight);  h_m3_enri_1000[1]->Fill(m3, this_weight);  }
	  else if(isEnri1) {  h_m3_enri[0]->Fill(m3, this_weight);  h_m3_enri_1000[0]->Fill(m3, this_weight);  }
	  else if(isBce3)  {  h_m3_bce[2]->Fill( m3, this_weight);  h_m3_bce_1000[2]->Fill( m3, this_weight);  }
	  else if(isBce2)  {  h_m3_bce[1]->Fill( m3, this_weight);  h_m3_bce_1000[1]->Fill( m3, this_weight);  }
	  else if(isBce1)  {  h_m3_bce[0]->Fill( m3, this_weight);  h_m3_bce_1000[0]->Fill( m3, this_weight);  }
	}
	else if (isSingleTop) { h_m3_singletop->Fill(m3, this_weight); h_m3_singletop_1000->Fill(m3, this_weight); }
	else if (isVQQ)       { h_m3_vqq->Fill(m3, this_weight);       h_m3_vqq_1000->Fill(m3, this_weight); }
      }//MC

    }
    //printf("selected combination = %d %d %d\n", which_3jet[0], which_3jet[1], which_3jet[2] ) ;
    //cout << "selected combination (idx) = " << which_3jet[0] << " "<< which_3jet[1] << " " << which_3jet[2] << endl;
  }// end pass_iso_ele (signal region)
  else {   // fail pass_iso_ele, this is bg region
    if ( jetColl.size() >= 4 ) {
      h_hadTop_maxPT_pt_nonIso_4j->Fill( max_top_PT, this_weight );
      h_hadTop_maxPT_mass_nonIso_4j->Fill( m3, this_weight );
      h_hadTop_maxPT_mass_nonIso_4j_1000->Fill( m3, this_weight );
      if (IsData()==false) {  // fill according to type of MC
	if (isTTbar)       {  h_m3_tt_control->Fill(     m3, this_weight);  h_m3_tt_control_1000->Fill(     m3, this_weight); }
	else if (isWjets)  {  h_m3_wj_control->Fill(     m3, this_weight);  h_m3_wj_control_1000->Fill(     m3, this_weight); }
	else if (isZjets)  {  h_m3_zj_control->Fill(     m3, this_weight);  h_m3_zj_control_1000->Fill(     m3, this_weight); }
	else if (isQCD) {     h_m3_qcd_control->Fill(    m3, this_weight);  h_m3_qcd_control_1000->Fill(    m3, this_weight); 
	  if(isEnri3)       { h_m3_enri_control[2]->Fill(m3, this_weight);  h_m3_enri_control_1000[2]->Fill(m3, this_weight); }
	  else if(isEnri2)  { h_m3_enri_control[1]->Fill(m3, this_weight);  h_m3_enri_control_1000[1]->Fill(m3, this_weight); }
	  else if(isEnri1)  { h_m3_enri_control[0]->Fill(m3, this_weight);  h_m3_enri_control_1000[0]->Fill(m3, this_weight); }
	  else if(isBce3)   { h_m3_bce_control[2]->Fill( m3, this_weight);  h_m3_bce_control_1000[2]->Fill( m3, this_weight); }
	  else if(isBce2)   { h_m3_bce_control[1]->Fill( m3, this_weight);  h_m3_bce_control_1000[1]->Fill( m3, this_weight); }
	  else if(isBce1)   { h_m3_bce_control[0]->Fill( m3, this_weight);  h_m3_bce_control_1000[0]->Fill( m3, this_weight); }
	}
	else if (isSingleTop) { h_m3_singletop_control->Fill(m3, this_weight);  h_m3_singletop_control->Fill(m3, this_weight); }
	else if (isVQQ)       { h_m3_vqq_control->Fill(m3, this_weight);        h_m3_vqq_control->Fill(m3, this_weight); }
      }
    }
  }

}// end reco_hadronicTop_highestTopPT
//-------------------------------------------------------------------------------------
// NEW
// return <M3,3j_vecsum_pt>
pair<double, double> ana::compute_M3( const vector<TLorentzVector>& jetColl ) const {

  double max_top_PT   = 0;
  double had_top_mass = 0; //ie m3
  vector<unsigned int> which_3jet(3); //3 places for int

  // Permute over all 3-jet combination, compute vector sum PT, find 3 jets which give the highest value
  for ( unsigned int a=0; a < (jetColl.size()-2); ++a ) {
        
    // const pat::Jet& j1 = jetColl[a];
    TLorentzVector j1( jetColl[a].Px(), jetColl[a].Py(), jetColl[a].Pz(), jetColl[a].Energy() );
    //cout << "j1 Pt=" << j1.Pt() << " (check: " << jetColl[a].Pt() << ")" << endl;
    
    for ( unsigned int b=a+1; b < (jetColl.size()-1); ++b ) {
      
      TLorentzVector j2( jetColl[b].Px(), jetColl[b].Py(), jetColl[b].Pz(), jetColl[b].Energy() );
      //cout << "j2 Pt=" << j2.Pt() << endl;
      
      for ( unsigned int c=b+1; c < jetColl.size(); ++c ) {
	
	TLorentzVector j3( jetColl[c].Px(), jetColl[c].Py(), jetColl[c].Pz(), jetColl[c].Energy() );
	//cout << "j3 Pt=" << j3.Pt() << endl;
	
	TLorentzVector top( j1 + j2 + j3 ); //construct a hadronic top candidate
	//cout << "top Pt=" << top.Pt() << "  mass=" << top.M() << endl;
	double this_top_PT = top.Pt();
	
	// find top quark candidate with the highest PT (ie 3-jet combination with highest PT)
	if ( this_top_PT > max_top_PT ) {
	  max_top_PT    = this_top_PT;
	  had_top_mass  = top.M();  //invariant mass of top
	  which_3jet[0] = a; 
	  which_3jet[1] = b; 
	  which_3jet[2] = c;
	}
      }//3rd jet
    }//2nd jet
  }//1st jet
  //cout << "BBB selected Top (max pt):  pt=" << max_top_PT << "  m=" << had_top_mass << endl;
  
  return pair<double,double>(had_top_mass, max_top_PT); //M3,pt
}



//----------------------------------------------------------------------------------------
//                                 W+jets Estimation
//----------------------------------------------------------------------------------------
// Descriptions:
//   o  The code are designed to run on data _or_ MC.
//   o  When running on data, it takes the m3 template shapes for ttbar and W+jets from MC.
//      The shapes for QCD is taken from data in the control region (nonIsoE), that is 
//      we assume QCD is dominating.
//   o  When running on MC only, it takes the m3 shapes from MC for ttbar, W+jets and QCD.
//      The shapes in the signal region are used to generate the pseudo-data. 
//      The number of event generated is the sum of Poisson-fluctuated n(tt), n(wj) and n(QCD).
//      To do the fit, we use template shape from MC (signal region) for ttbar and W+jets;
//      for QCD, we can choose to use the m3 distribution from control region for: 
//      (a) all events (tt+wj+QCD), or (b) only QCD events. At the moment the code assume
//      approach (b). (To use (a), uncomment lines that says "Use all events for QCD template")
//   o  In both cases, the size of QCD int the fit is constrained to be a Gaussian, with 
//      mean equals to the QCD estimate, and width equals to the uncertainty of the estimate.
//   o  The code assumes that the various m3 template histograms from MC are read in from 
//      a root file which will be made using a separate macro.
//   o  The fit produces the number of ttbar, W+jets and QCD. It also produces a root file 
//      called m3_out.root to store histograms.
//
// IMPORTANT NOTICE: the fit code requires ROOT v5.22.00.
//

void ana::EstimateWjets() {

  EstimateWjets(outputHistFileName,""); //run on MC
}

bool ana::EstimateWjets(const string inputFile_data, const string inputFile_mc) {
  
  cout << "--------------------------------------------------------" << endl;
  cout << "                    W+jets Estimation" << endl;
  cout << "--------------------------------------------------------" << endl;

  vector<double> nFitResults;


  //--------------------------
  // import m3 distributions
  //--------------------------
  // aim: to import m3 from one root file

  TFile *fdata(0); //data: real or simulated
  TFile *fmc(0);   //template from MC
  fdata = new TFile( inputFile_data.c_str(), "UPDATE" );
  if (IsData()) {
    fmc = new TFile( inputFile_mc.c_str() ); //separate file
  }else{
    fmc = fdata; //point to same file
  }

  cout << "inputFile_data = " << inputFile_data << endl;
  cout << "inputFile_mc = " << inputFile_mc << endl;



  // used as template to fit, also to generate pseudo-data
  // m3 in signal region (iso)  
  string dirW = "Wjets_estimation"; //960 bins
  if( m3_use_1000_bins ) dirW = "Wjets_estimation_1000"; //1000 bins
  cout << "--------------------------------------------"<< endl;
  cout << "Taking m3 histograms from dir: " << dirW << endl;
  cout << "--------------------------------------------"<< endl;
  TH1D *h_m3_tt   = (TH1D*)fmc->Get(Form("%s/m3_tt",dirW.c_str()));
  TH1D *h_m3_wj   = (TH1D*)fmc->Get(Form("%s/m3_wj",dirW.c_str()));

  // 21-5-09: Take wj template from Fastsim normal MC
  const bool take_wj_template_from_fastsim = false;
  if (take_wj_template_from_fastsim==true) {
    cout << "NOTE: Taking W+jets shape from wj Fastsim as the default template." << endl;
    TFile *fwj_fastsim = new TFile("test_mc_wjet_fastsim_nom.root");
    h_m3_wj = (TH1D*)fwj_fastsim->Get(Form("%s/m3_wj",dirW.c_str()));
  }



  TH1D *h_m3_zj(0);
  TH1D *h_m3_stop(0);
  //TH1D *h_m3_qcd(0);         //NB: for QCD template, we use data in control region (was called "control_all")
  TH1D *h_m3_qcd_control(0); //only used when doing "QCD template unc"

  if (IsData()==false) { //MC
    h_m3_qcd    = (TH1D*)fmc->Get(Form("%s/m3_qcd",dirW.c_str()));
    h_m3_zj     = (TH1D*)fmc->Get(Form("%s/m3_zj",dirW.c_str()));
    h_m3_stop   = (TH1D*)fmc->Get(Form("%s/m3_singletop",dirW.c_str()));
    h_m3_qcd_control = (TH1D*)fmc->Get(Form("%s/m3_qcd_control",dirW.c_str()));
  }



  TH1D *h_m3_MCdata_control = (TH1D*)fdata->Get(Form("%s/hadTop_maxPT_mass_nonIso_4j",dirW.c_str()));
  TH1D *control_all = (TH1D*)h_m3_MCdata_control->Clone("control_all");



  //-----------------------------------------------------
  // May 22: adapt code to use Thorsten's M3 histograms
  //-----------------------------------------------------
  const bool use_Thorsten_histo = false;
  const bool use_Thorsten_eventNo = false;

  if(use_Thorsten_histo) {
    cout << "\n *** Using Thorsten's M3 histograms to fit ***" << endl;
  
    TFile *ftt = new TFile("Thorsten/ttbar.root");
    TFile *fwj = new TFile("Thorsten/wjets_Fastsim.root"); //from FASTsim
    TFile *fqcd = new TFile("Thorsten/QCD_template.root");
    TFile *fstop = new TFile("Thorsten/stop.root");
    
    h_m3_tt = (TH1D*)ftt->Get("h_m3_100");
    h_m3_wj = (TH1D*)fwj->Get("h_m3_100");
    control_all = (TH1D*)fqcd->Get("h_m3_100");
    h_m3_stop = (TH1D*)fstop->Get("h_m3_100");
  }

  /*
  // 27 May 09: use KA's ttjet shape  
  TFile *ftt = new TFile("Thorsten/ttbar.root");
  h_m3_tt = (TH1D*)ftt->Get("h_m3_100");
  cout << "\n*** NOTE: For signal ttjet, use KA shape. ***" << endl;
  */
  /*
  // 27 May 09: use ttbar pythia as default shape  
  TFile *ftt = new TFile("test_mc_ttbar.root");
  h_m3_tt = (TH1D*)ftt->Get("Wjets_estimation/m3_tt");
  cout << "\n*** NOTE: For signal, use ttbar pythia shape. ***" << endl;
  */
  
  // 1 June 09: to do ttbar ISR/FSR systematics, use fastsim ttbar nominal as reference  
  if( sysSample.find("ttbar_fastsim_ISRFSR") != string::npos ) { //match "ttbar_fastsim_ISRFSR"
    TFile *ftt = new TFile("test_mc_ttbar_fastsim_ISRFSR_nom.root");
    h_m3_tt = (TH1D*)ftt->Get(Form("%s/m3_tt",dirW.c_str()));
    cout << "\n*** NOTE: For signal, use Fastsim ttbar pythia (nominal ISR/FSR) shape. ***" << endl;
  }
  
  // 18 June 09: to do ttjet scale/threshold systematics, use fastsim ttjet nominal as reference
  if( sysSample.find("ttjet_fastsim") != string::npos ) { //match "ttjet_fastsim"
    TFile *ftt = new TFile("test_mc_ttjet_fastsim_nom.root");
    h_m3_tt = (TH1D*)ftt->Get(Form("%s/m3_tt",dirW.c_str()));
    cout << "\n*** NOTE: To do ttjet scale/threshold systematics, use normal Fastsim ttjet shape as template. ***\n" << endl;
  }
  


  cout << "m3 distributions:" << endl;
  cout << "pointer tt: "<< h_m3_tt << endl;
  cout << "pointer wj: "<< h_m3_wj << endl;
  cout << "pointer zj: "<< h_m3_zj << endl;
  //cout << "pointer qcd: "<< h_m3_qcd << endl;
  if(h_m3_tt>0)  h_m3_tt->ls();
  if(h_m3_wj>0)  h_m3_wj->ls();
  //if(h_m3_qcd>0) h_m3_qcd->ls();
  if ( h_m3_tt==0 || h_m3_wj==0 || control_all==0 ) {
    cout << "Warning: zero pointer to either h_m3_tt or h_m3_wj or control_all. Exit EstimateWjets()." << endl;
    return false;
  }
  // If we don't find the m3 histograms, stop execution of code.
  if ( IsData()==false ) {
    if (h_m3_zj==0 || h_m3_stop==0)  {
      cout << "Warning: could not find some m3 histogram(s). Exit EstimateWjets()" << endl;
      return false;
    }
  }


  // Define (default) fit templates
  //--------------------------------
  TH1D *temp_tt = (TH1D*)h_m3_tt->Clone("temp_tt");
  TH1D *temp_wj = (TH1D*)h_m3_wj->Clone("temp_wj");
  TH1D *temp_stop = (TH1D*)h_m3_stop->Clone("temp_stop");
  TH1D *temp_qcd = (TH1D*)control_all->Clone("temp_qcd");




  // Define (default) pseudodata models
  //------------------------------------
  TH1D *model_tt = (TH1D*)h_m3_tt->Clone("model_tt");
  TH1D *model_wj = (TH1D*)h_m3_wj->Clone("model_wj");
  TH1D *model_zj = (TH1D*)h_m3_zj->Clone("model_zj");
  TH1D *model_stop = (TH1D*)h_m3_stop->Clone("model_stop");
  TH1D *model_qcd = (TH1D*)control_all->Clone("model_qcd");



  // pick one for QCD
  //// NOTE: switch to use control region to generate if there are very few QCD events in signal region.
  //// TO BE SWITCH to signal region if we have enough QCD events	
  //TH1D *h_qcd_model = (TH1D*)h_m3_control_qcd->Clone("model_qcd");  //<--- control region
  //TH1D *h_qcd_model = (TH1D*)h_m3_qcd->Clone("model_qcd");   //<--- signal region


  //---------------------------------------------
  // 25 Apr 09: Systematic Shape Uncertainties
  //---------------------------------------------
  // Swith: Use a different shape to evaluate systematic unc
  
  if(doSystematics=="shape"||doSystematics=="accepShape"){

    cout << "\nDoing *shape* uncertainty\n" << endl;
    
    TFile *fmc2(0);

    // (1) if vary  ttbar
    if ( sysSample.find( "tt" ) != string::npos ) { //match "tt"

      if(!use_Thorsten_histo) { //use our histo

	if(sysSample=="ttbar")                 fmc2 = new TFile("test_mc_ttbar.root");
	else if(sysSample=="ttjet_largerISR")  fmc2 = new TFile("test_mc_ttjet_largerISR.root"); //old
	else if(sysSample=="ttjet_smallerISR") fmc2 = new TFile("test_mc_ttjet_smallerISR.root"); //old
	// new pythia ttbar ISR/FSR sample (1 June)
	//else if(sysSample=="ttbar_fastsim_ISRFSR_nom")   fmc2 = new TFile("test_mc_ttbar_fastsim_ISRFSR_nom.root");
	else if(sysSample=="ttbar_fastsim_ISRFSR_large") fmc2 = new TFile("test_mc_ttbar_fastsim_ISRFSR_large.root");
	else if(sysSample=="ttbar_fastsim_ISRFSR_small") fmc2 = new TFile("test_mc_ttbar_fastsim_ISRFSR_small.root");
	// new sample 15-Jun-09
	else if(sysSample=="ttjet_fastsim_nom")  fmc2 = new TFile("test_mc_ttjet_fastsim_nom.root");
	else if(sysSample=="ttjet_fastsim_T40")  fmc2 = new TFile("test_mc_ttjet_fastsim_T40.root");
	else if(sysSample=="ttjet_fastsim_T10")  fmc2 = new TFile("test_mc_ttjet_fastsim_T10.root");
	else if(sysSample=="ttjet_fastsim_scaleUp")   fmc2 = new TFile("test_mc_ttjet_fastsim_scaleUp.root");
	else if(sysSample=="ttjet_fastsim_scaleDown") fmc2 = new TFile("test_mc_ttjet_fastsim_scaleDown.root");

	if(fmc2>0) model_tt = (TH1D*)fmc2->Get(Form("%s/m3_tt",dirW.c_str()))->Clone("model_tt");

      }else{ //use Thorsten's histo
	cout << "NOTE: Use Thorsten's histo for ttbar shape" << endl;
	if(sysSample=="ttbar")  fmc2 = new TFile("Thorsten/ttbar_Pythia.root");
	model_tt = (TH1D*)fmc2->Get("h_m3_100")->Clone("model_tt");
      }
    }
    
    // (2) if vary  wj    
    //  TFile *fmc2 = new TFile("test_mc_SETTHIS.root");    
    else if ( sysSample.find( "wj" ) != string::npos ) { //match "wj"

      if(!use_Thorsten_histo) { //use our histo

	if(sysSample=="wjet_fastsim_nom")    fmc2 = new TFile("test_mc_wjet_fastsim_nom.root");
	else if(sysSample=="wjet_thres20")   fmc2 = new TFile("test_mc_wjet_thres20.root");
	else if(sysSample=="wjet_thres5")    fmc2 = new TFile("test_mc_wjet_thres5.root");
	else if(sysSample=="wjet_scaleUp")   fmc2 = new TFile("test_mc_wjet_scaleUp.root");
	else if(sysSample=="wjet_scaleDown") fmc2 = new TFile("test_mc_wjet_scaleDown.root");
	model_wj = (TH1D*)fmc2->Get(Form("%s/m3_wj",dirW.c_str()))->Clone("model_wj");
      
      }else{ //use Thorsten's histo
	cout << "NOTE: Use Thorsten's histo for wj systematics shape" << endl;
	if(sysSample=="wjet_fastsim_nom")    fmc2 = new TFile("Thorsten/wjets_FastSim.root");
	else if(sysSample=="wjet_thres20")   fmc2 = new TFile("Thorsten/wjets_FastSim_T20.root");
	else if(sysSample=="wjet_thres5")    fmc2 = new TFile("Thorsten/wjets_FastSim_T5.root");
	else if(sysSample=="wjet_scaleUp")   fmc2 = new TFile("Thorsten/wjets_FastSim_scaleup.root");
	else if(sysSample=="wjet_scaleDown") fmc2 = new TFile("Thorsten/wjets_FastSim_scaledown.root");
	model_wj = (TH1D*)fmc2->Get("h_m3_100")->Clone("model_wj");       
      }
    }
    cout << ",  file: " << fmc2->GetName() << endl;
  }//end if vary shape

  /*
  // 27 May 09
  cout << "NOTE: Use Thorsten's pythia histo for ttbar shape (as model)" << endl;
  TFile *fmc2 = new TFile("Thorsten/ttbar_Pythia.root");
  model_tt = (TH1D*)fmc2->Get("h_m3_100")->Clone("model_tt");
  */
  

  /*
  cout << "\nNOTE: Use wj fastsim as generation model (pseudodata)\n" << endl;
  TFile *fmc2 = new TFile("test_mc_wjet_fastsim_nom.root");
  model_wj = (TH1D*)fmc2->Get("Wjets_estimation/m3_wj")->Clone("model_wj");

  cout << "\nNOTE: Use ttbar pythia as generation model (pseudodata)\n" << endl;
  TFile *fmc3 = new TFile("test_mc_ttbar.root");
  model_tt = (TH1D*)fmc3->Get("Wjets_estimation/m3_tt")->Clone("model_wj");
  */
  


  //Uncomment the next five lines if carrying out JES systematic, and change pdf inputs to *2 (i.e. h_m3_tt -> h_m3_tt2). fmc2 should point to a file with the standard templates. 
  //    TFile *fmc2 = new TFile("test_mc_mixture.root");
  //TH1D *h_m3_tt2 = (TH1D*)fmc2->Get("Wjets_estimation/m3_tt")->Clone("model_tt");
  //TH1D *h_m3_wj2 = (TH1D*)fmc2->Get("Wjets_estimation/m3_wj")->Clone("model_wj");
  //h_m3_tt2->Rebin(8);
  //h_m3_wj2->Rebin(8);


  //---------------------------
  // SYS: Unc due to QCD template 
  //---------------------------
  if( doSystematics=="QCDtemplate" ) {
    cout << "For 8.1.2, Uncertainty due to QCD template" << endl;
    cout << "Use only QCD events in control region as data model" << endl;
    model_qcd = (TH1D*)h_m3_qcd_control->Clone("model_qcd");
  }





  // 1-2-09: use more appropriate binning (1000/12.5 bins = 80GeV/bin)
  //const int m3_bin_used_in_AN = 80; //GeV
  int rB = 1;
  //  const int rB = nbm3 / m3_bin_used_in_AN; //800 bins / 10 = 80GeV/bin
  if      ( temp_tt->GetNbinsX() == 100 )  rB = 8; //12.5 bins
  else if ( temp_tt->GetNbinsX() == 960 )  rB = 80; //12 bins
  else if ( temp_tt->GetNbinsX() == 1000 )  rB = 80; //12.5 bins
      
  temp_tt->Rebin(rB); 
  temp_wj->Rebin(rB);
  temp_qcd->Rebin(rB);
  temp_stop->Rebin(rB); 
 
  model_tt->Rebin(rB); 
  model_wj->Rebin(rB);
  model_zj->Rebin(rB);
  model_qcd->Rebin(rB);
  model_stop->Rebin(rB); 

  if(temp_tt->GetNbinsX()==75)  temp_tt->Rebin(6);
  if(model_tt->GetNbinsX()==75) model_tt->Rebin(6);

  //----------------------------------------







  //----------------------------------------
  //  If running on data, nfit is set automatically to 1.
  //  If running on MC, can specify how many toy MC to generate and fit
  //
  int nfit = m_ntoy;
  if (IsData()) nfit = 1;
  //----------------------------------------

  //-----------------------------------------------------
  // true value: expected signal and bgnd events from MC

  TH2D *SigEvents   = (TH2D*)fmc->Get("Signal_njetsVcuts");
  TH2D *WjetsEvents = (TH2D*)fmc->Get("Wjets_njetsVcuts");
  TH2D *ZjetsEvents = (TH2D*)fmc->Get("Zjets_njetsVcuts");
  TH2D *QCDEvents   = (TH2D*)fmc->Get("QCD_njetsVcuts");
  TH2D *SingleTopEvents = (TH2D*)fmc->Get("SingleTop_njetsVcuts");

  /*  
  // use only when redoing m3 fit for 200/pb using 20/pb hist file
  if(1){ //multiply by 10
    SigEvents->Scale(10);
    WjetsEvents->Scale(10);
    ZjetsEvents->Scale(10);
    QCDEvents->Scale(10);
    SingleTopEvents->Scale(10);
  }
  */



  // new event yields for wj fastsim sys samples (to do accept unc)
  //const double  ntrue_wjet_fastsim_nom = 79.9; //91.77;
  double  ntrue_wjet_T20  = 108.9; //125.06;
  double  ntrue_wjet_T5   = 82.7; //95.01;
  double  ntrue_wjet_scaleUp   = 54.6; //62.73;
  double  ntrue_wjet_scaleDown = 160.3; //184.09;
  double  ntrue_ttbar          = 180.7;
  //  const double  ntrue_ttjet_largerISR   = 208.93; 
  //  const double  ntrue_ttjet_smallerISR  = 215.38; 
  //  const double  ntrue_ttbar_fastsim_ISRFSR_nom   = 160.9;
  const double  ntrue_ttbar_fastsim_ISRFSR_large = 181.95; // unscaled = 159.7;
  const double  ntrue_ttbar_fastsim_ISRFSR_small = 169.78; // unscaled = 149.0;
  // ttjet fastsim (scaled to Fullsim)
  const double  ntrue_ttjet_fastsim_T40 = 188.4;
  const double  ntrue_ttjet_fastsim_T10 = 190.1;
  const double  ntrue_ttjet_fastsim_scaleUp = 187.3;
  const double  ntrue_ttjet_fastsim_scaleDown = 191.8;


  if(use_Thorsten_eventNo){
    ntrue_wjet_T20 = 74;
    ntrue_wjet_T5  = 55;
    ntrue_wjet_scaleUp = 36;
    ntrue_wjet_scaleDown = 99;
    ntrue_ttbar = 172;
  }


  double ntt_true  = SigEvents->GetBinContent(5,5); //180.70
  double nwj_true  = WjetsEvents->GetBinContent(5,5);
  double nzj_true  = ZjetsEvents->GetBinContent(5,5);
  double nqcd_true = QCDEvents->GetBinContent(5,5);
  double nstop_true = SingleTopEvents->GetBinContent(5,5);
  /*
  ntt_true = 180.7; // pythia
  cout << "\n NOTE: normalize n(signal) to " << ntt_true << endl;

  */


  if(use_Thorsten_eventNo){
    // Use Thorsten's event yields
    cout << "NOTE: Using Thorsten's event numbers." << endl;
    ntt_true = 177;
    nwj_true = 60;
    nzj_true = 12;
    nqcd_true = 30;
    nstop_true = 8;    
  }

  const double nstop_exp = nstop_true; //used in Gaus-constrain


  //----------------------------------
  // QCD normalization uncertainty
  //----------------------------------

  // Default
  double nqcd_est     = nqcd_true;    // = mean of Gaussian constrain
  double nqcd_est_err = nqcd_true*0.5; // = width of Gaussian constrain

  // Do as Thorsten, vary the "true" mean in number of QCD evnt in pseudodata 

  // Our numbers
  if(doSystematics=="QCDnorm_genest" )  {
    nqcd_true    = 17.0;  // QCd est
    nqcd_est     = 17.0;
    nqcd_est_err = 8.5;
  } else if(doSystematics=="QCDnorm_genest15" ) {
    nqcd_true    = 25.5;  // QCD est +50%
    nqcd_est     = 17.0;
    nqcd_est_err = 8.5;
  } else if(doSystematics=="QCDnorm_genest05" ) {
    nqcd_true    = 8.5;  // QCD est -50%
    nqcd_est     = 17.0;
    nqcd_est_err = 8.5;
  }

  
  if(use_Thorsten_eventNo) {
    // Thorsten's numbers
    if(doSystematics=="QCDnorm_genest15" ) {
      nqcd_true    = 45;  // QCD est +50%
      nqcd_est     = 30;
      nqcd_est_err = 15;
    } else if(doSystematics=="QCDnorm_genest05" ) {
      nqcd_true    = 15;  // QCD est -50%
      nqcd_est     = 30;
      nqcd_est_err = 15;
    }
  }
  cout << "QCD mean of pseudodata (true) = " << nqcd_true << endl;
  cout << "QCD Gaussian-constraint, mean = " << nqcd_est;
  cout << "  width = " << nqcd_est_err << endl;


  //--------------------------------------
  // Single Top normalization uncertainty
  //--------------------------------------
  const double nstop_true_plus = nstop_true *1.46;
  const double nstop_true_minus = nstop_true *(1-0.46);
  if(doSystematics=="stop_plus" )  { //+46%
    cout << "NOTE: more single top in pseudodata: + 46%, " << nstop_true_plus << endl;
    nstop_true  = nstop_true_plus  ;
  } else if(doSystematics=="stop_minus" ) { //-46%
    cout << "NOTE: less single top in pseudodata: - 46%, " << nstop_true_minus << endl;
    nstop_true = nstop_true_minus  ;
  }
  





  //------------------------------------
  // Acceptance (rate) uncertainty
  //------------------------------------
  // 19-5-09
  if(doSystematics=="accep"||doSystematics=="accepShape") {
    cout << "\n Doing *acceptance* uncertainty\n" << endl;
    // wj
    if(sysSample=="wjet_fastsim_nom")    {
      cout << "use the same nwj_true from fullsim" << endl;
      //nwj_true = ntrue_wjet_fastsim_nom;
    }
    else if(sysSample=="wjet_thres20")   nwj_true = ntrue_wjet_T20;
    else if(sysSample=="wjet_thres5")    nwj_true = ntrue_wjet_T5;
    else if(sysSample=="wjet_scaleUp")   nwj_true = ntrue_wjet_scaleUp;
    else if(sysSample=="wjet_scaleDown") nwj_true = ntrue_wjet_scaleDown;
    // tt
    else if(sysSample=="ttbar")    ntt_true = ntrue_ttbar;
    else if(sysSample=="ttbar_fastsim_ISRFSR_nom")  {
      cout << "use the same ntt_true from fullsim" << endl;
      //ntt_true = ntrue_ttbar_fastsim_ISRFSR_nom;
    }
    else if(sysSample=="ttbar_fastsim_ISRFSR_large") ntt_true = ntrue_ttbar_fastsim_ISRFSR_large;
    else if(sysSample=="ttbar_fastsim_ISRFSR_small") ntt_true = ntrue_ttbar_fastsim_ISRFSR_small;
    else if(sysSample=="ttjet_fastsim_T40")       ntt_true = ntrue_ttjet_fastsim_T40;
    else if(sysSample=="ttjet_fastsim_T10")       ntt_true = ntrue_ttjet_fastsim_T10;
    else if(sysSample=="ttjet_fastsim_scaleUp")   ntt_true = ntrue_ttjet_fastsim_scaleUp;
    else if(sysSample=="ttjet_fastsim_scaleDown") ntt_true = ntrue_ttjet_fastsim_scaleDown;    
  }


  cout << "MC true expected number (signal region, isol)" << endl;
  cout << "True n(tt):  " <<  ntt_true << endl;
  cout << "True n(wj):  " <<  nwj_true << endl;
  cout << "True n(zj):  " <<  nzj_true << endl;
  cout << "True n(qcd): " << nqcd_true << endl;
  cout << "True n(singletop): " << nstop_true << endl;

  



  //-----------------------------------------------------
//   // constrain on QCD: QCD estimate
//   double nqcd_est     = nqcd_true;     // ideal, take true MC expectation <--- we'll put our QCD estimate here
//   double nqcd_est_err = nqcd_true*0.5;  // 50%
//   // To do QCD normalization uncertainty
//   if ( doSystematics.find( "QCDnorm" ) != string::npos ) { //match "QCDnorm"
//     cout << "Doing QCD normalizatoin uncertainy" << endl;
//     nqcd_est = 17.0;
//     nqcd_est_err = 8.5;
//     cout << "In QCD Gaussian-constraint, mean = " << nqcd_est;
//     cout << "  width = " << nqcd_est_err << endl;
//   }

  /*
    nqcd_est = 17.0; //nqcd_true;   //ideal case
    //nqcd_est = 17.0; //to do QCD norm sys //est=17.0, 25.5, 8.5
    //    nqcd_est_err = sqrt(nqcd_true);
    nqcd_est_err = 8.5; // nqcd_est*0.5; //do the same as Thorsten, 50% of the estimate.
  */
  /*

  TH1D *qcdEstimate = (TH1D*)fdata->Get("QCDEstimate");
  //qcdEstimate->Scale(10);  //use only when redoing m3 fit for 200/pb using 20/pb hist file

  //  if(false){//some statement to check we have done qcd estimate (xxxxxxxx)
  if(qcdEstimate>0){ //some statement to check we have done qcd estimate
    nqcd_est     = qcdEstimate->GetBinContent(1);
    nqcd_est_err = qcdEstimate->GetBinContent(4); //average error
    // 13 Feb 09: if error is larger than estimate, assume 100% error
    if (nqcd_est_err > nqcd_est) nqcd_est_err = nqcd_est;
  }
  */
  cout << "QCD Estimate used: n(qcd):  " <<  nqcd_est << " +/- "<< nqcd_est_err << endl;
  //-----------------------------------------------------



  //---------------
  // add  singlet top as Gaussian constraint like qcd
  const double nstop_exp_err = sqrt( nstop_exp + 0.3*0.3*nstop_exp*nstop_exp );
  //---------------

  //------------------------------------------------
  //  Create a new output root file
  //------------------------------
  //    TFile fout("out.root","recreate");
  //  TFile *fout = new TFile("m3_out.root","recreate");
  //------------------------------------------------
  
  // Book histograms to store fit results
  TH1D *h_ntt_fit     = new TH1D("h_ntt_fit",  "n(tt) fit",  100, 0, 500*intlumi/20);
  TH1D *h_nwj_fit     = new TH1D("h_nwj_fit",  "n(wj) fit",  100, -100, 450*intlumi/20);
  TH1D *h_nqcd_fit    = new TH1D("h_nqcd_fit", "n(qcd) fit", 100, 0*intlumi/20, 80*intlumi/20); //rescale (20)
  TH1D *h_nstop_fit   = new TH1D("h_nstop_fit", "n(singletop) fit", 1000, 0, 20*intlumi/20);
  // good fit only
  TH1D *h_ntt_goodfit   = new TH1D("h_ntt_goodfit",  "n(tt) goodfit",  100, 0, 500*intlumi/20);
  TH1D *h_nwj_goodfit   = new TH1D("h_nwj_goodfit",  "n(wj) goodfit",  100, -100, 450*intlumi/20);
  TH1D *h_nqcd_goodfit  = new TH1D("h_nqcd_goodfit", "n(qcd) goodfit", 100, 0*intlumi/20, 80*intlumi/20);
  TH1D *h_nstop_goodfit = new TH1D("h_nstop_goodfit", "n(singletop) goodfit", 1000, 0, 20*intlumi/20);

  TH1D *h_ntt_fiterr  = new TH1D("h_ntt_fiterr",  "error of n(ntt) fit",  1000, 0, 200);
  TH1D *h_nwj_fiterr  = new TH1D("h_nwj_fiterr",  "error of n(wj) fit",   1000, 0, 200);
  TH1D *h_nqcd_fiterr = new TH1D("h_nqcd_fiterr", "error of n(qcd) fit",  1000, 0, 200);//mod
  TH1D *h_nstop_fiterr = new TH1D("h_nstop_fiterr", "error of n(singletop) fit",  1000, 0, 200);

  // Declare histograms to store number of generated (Poisson-fluctuated) events 
  // (when running on MC).
  TH1D *h_ntt_fluc(0);
  TH1D *h_nwj_fluc(0);
  TH1D *h_nzj_fluc(0);
  TH1D *h_nqcd_fluc(0);
  TH1D *h_nstop_fluc(0);
  TH1D *h_ntt_pull(0); //pull of fit results = (fit-true)/error
  TH1D *h_nwj_pull(0);
  TH1D *h_nqcd_pull(0);
  TH1D *h_nstop_pull(0);
  TH1D *h_ntt_goodpull(0); //pull of fit results = (fit-true)/error
  TH1D *h_nwj_goodpull(0);
  TH1D *h_nqcd_goodpull(0);
  TH1D *h_nstop_goodpull(0);


  if ( IsData()==false ) {
    h_ntt_fluc  = new TH1D("h_ntt_fluc", "generated number of tt events (Poisson)", 100,0,1000);
    h_nwj_fluc  = new TH1D("h_nwj_fluc", "generated number of wj events (Poisson)", 100,0,500);
    h_nzj_fluc  = new TH1D("h_nzj_fluc", "generated number of zj events (Poisson)", 100,0,500);
    h_nqcd_fluc = new TH1D("h_nqcd_fluc","generated number of qcd events (Poisson)", 100,0,200); //rescale (20)
    h_nstop_fluc = new TH1D("h_nstop_fluc","generated number of singletop events (Poisson)", 100,0,100);

    // pull
    h_ntt_pull  = new TH1D("h_ntt_pull",  "pull of n(tt)",  100, -5, 5);
    h_nwj_pull  = new TH1D("h_nwj_pull",  "pull of n(wj)",  100, -5, 5);
    h_nqcd_pull = new TH1D("h_nqcd_pull", "pull of n(qcd)", 1200, -7, 5);
    h_nstop_pull = new TH1D("h_nstop_pull", "pull of n(singletop)", 1000, -5, 5);

    h_ntt_goodpull  = new TH1D("h_ntt_goodpull",  "pull of n(tt) good",  100, -5, 5);
    h_nwj_goodpull  = new TH1D("h_nwj_goodpull",  "pull of n(wj) good",  100, -5, 5);
    h_nqcd_goodpull = new TH1D("h_nqcd_goodpull", "pull of n(qcd) good", 1200, -7, 5);
    h_nstop_goodpull = new TH1D("h_nstop_goodpull", "pull of n(singletop) good", 1000, -5, 5);

  }
  // keep a few set of pseudo-data
  TH1D *h_m3_sum_toy_keep[5];
  //TH1D *h_m3_sum_ctr_toy_keep[5];
  
  for(int i=0; i<5; ++i){
    h_m3_sum_toy_keep[i]     = new TH1D(Form("h_m3_sum_toy_%d",i+1),     Form("m3 sum toy mc %d",i+1), 100,0,1000);
    //h_m3_sum_ctr_toy_keep[i] = new TH1D(Form("h_m3_sum_ctr_toy_%d",i+1), Form("m3 sum ctr toy mc %d",i+1), 100,0,1000);
  }
  

  // keep data (real/sim) to write to root file
  TH1D *h_m3_data_keep(0);
  TH1D *h_m3_data_control_keep(0);

  // keep a few qcd pdf (use data in control region)
  RooPlot* m3_qcd_pdf[5];

  bool printText = true;

  RooPlot* m3frame1(0);
  RooPlot* m3frame2(0);
  RooPlot* m3frame3(0);
  RooPlot* m3frame4(0);
  RooPlot* m3frame5(0);
  RooPlot* m3frame6(0);
  RooPlot* m3frame7(0);
  RooPlot* m3frame8(0);
  RooPlot* m3frame9(0);
  RooPlot* m3frame10(0);
  RooPlot* m3frame11(0);
  //RooPlot* m3frame21(0);
  //RooPlot* m3frame22(0);


  //---------------
  // perform fit
  //---------------


  for ( int n=0; n<nfit; ++n ) {


    if (n>20) printText = false;
    cout << "Fits carried out: " << n+1 << " of " << nfit << endl;
    if( (n+1)%100 == 1 || (n+1)==nfit ){
      printText = true;
    }
    
    //--------------------------------------
    // 1) Get the m3 distribution 
    //    (can be real data or pseudo-data)
    //--------------------------------------
    TH1D *h_m3_data(0);          // in signal region (Iso)
    TH1D *h_m3_data_control(0);  // in control region (nonIso) - take this to be QCD shape



    int ntt_fluc = 0;
    int nwj_fluc = 0;
    int nzj_fluc = 0;
    int nqcd_fluc = 0;
    int nstop_fluc = 0;


    if ( IsData()==true ) { // Real LHC data

      TH1D *h_m3_LHCdata         = (TH1D*)fdata->Get(Form("%s/hadTop_maxPT_mass_4j",dirW.c_str()));
      TH1D *h_m3_LHCdata_control = (TH1D*)fdata->Get(Form("%s/hadTop_maxPT_mass_nonIso_4j",dirW.c_str()));
      h_m3_data         = (TH1D*)h_m3_LHCdata->Clone("h_m3_data"); //rename
      h_m3_data_control = (TH1D*)h_m3_LHCdata_control->Clone("h_m3_data_control");
      
    } else {  // Monte Carlo

      // a) Poission fluctuated number of events
      ntt_fluc   = gRandom->Poisson(ntt_true);
      nwj_fluc   = gRandom->Poisson(nwj_true);
      nzj_fluc   = gRandom->Poisson(nzj_true);
      nqcd_fluc  = gRandom->Poisson(nqcd_true);
      nstop_fluc = gRandom->Poisson(nstop_true);
      h_ntt_fluc  ->Fill( ntt_fluc );
      h_nwj_fluc  ->Fill( nwj_fluc );
      h_nzj_fluc  ->Fill( nzj_fluc );
      h_nqcd_fluc ->Fill( nqcd_fluc );
      h_nstop_fluc->Fill( nstop_fluc );
      
      if(printText) {
	cout << "\n>>> toy experiment #" << n+1 << " of " << nfit << endl;
	cout << "-------------------------------------" << endl;
	cout << "Poisson-fluc n(tt):  " << ntt_fluc << endl;
	cout << "Poisson-fluc n(wj):  " << nwj_fluc << endl;
	cout << "Poisson-fluc n(zj):  " << nzj_fluc << endl;
	cout << "Poisson-fluc n(qcd): " << nqcd_fluc << endl;
	cout << "Poisson-fluc n(singletop): " << nstop_fluc << endl;
      }


      // b) Book histos to store pseudo-data (can use smaller binning if needed)
      //const int nbin = 12;
      //      const int max = 960;
      // TEST 15-Jan-10
      const int nbin = temp_tt->GetNbinsX();
      const float max = temp_tt->GetXaxis()->GetXmax();
      if(printText) cout << "\n m3 Pseudo-data: " << nbin << " bins, range 0-" << max << endl << endl;
      // signal region
      TH1D *h_m3_tt_toy   = new TH1D("h_m3_tt_toy",   "h_m3_tt_toy",   nbin,0,max);
      TH1D *h_m3_wj_toy   = new TH1D("h_m3_wj_toy",   "h_m3_wj_toy",   nbin,0,max);
      TH1D *h_m3_zj_toy   = new TH1D("h_m3_zj_toy",   "h_m3_zj_toy",   nbin,0,max);
      TH1D *h_m3_qcd_toy  = new TH1D("h_m3_qcd_toy",  "h_m3_qcd_toy",  nbin,0,max);
      TH1D *h_m3_stop_toy = new TH1D("h_m3_stop_toy", "h_m3_stop_toy", nbin,0,max);
      TH1D *h_m3_sum_toy  = new TH1D("h_m3_sum_toy",  "h_m3_sum_toy",  nbin,0,max);


      //-------------------------------------------------
      // c) Generate pseudo-data for this toy experiment
      //-------------------------------------------------
      for ( int i=0; i < ntt_fluc; ++i ) {	h_m3_tt_toy->Fill(   model_tt->GetRandom() );    }
      for ( int i=0; i < nwj_fluc; ++i ) {	h_m3_wj_toy->Fill(   model_wj->GetRandom() );    }
      for ( int i=0; i < nzj_fluc; ++i ) {	h_m3_zj_toy->Fill(   model_zj->GetRandom() );    }     
      for ( int i=0; i < nstop_fluc; ++i ) {    h_m3_stop_toy->Fill( model_stop->GetRandom() );  }
      for ( int i=0; i < nqcd_fluc; ++i ) {	h_m3_qcd_toy->Fill(  model_qcd->GetRandom() );   }
      



      // d) prepares the pseudo-data in signal region
      h_m3_sum_toy->Add( h_m3_tt_toy );
      h_m3_sum_toy->Add( h_m3_wj_toy );
      if (doSystematics=="addZj_fitWithWjShape"){
	h_m3_sum_toy->Add( h_m3_zj_toy );
      }
      h_m3_sum_toy->Add( h_m3_stop_toy );
      h_m3_sum_toy->Add( h_m3_qcd_toy );


      if (n<5) {
	h_m3_sum_toy_keep[n]     = (TH1D*)h_m3_sum_toy->Clone(Form("h_m3_sum_toy_%d",n+1));
	h_m3_sum_toy_keep[n]->SetTitle( Form("m3 sum toy mc %d",n+1) );
      }


      const double n_gen_event = h_m3_sum_toy->GetEntries();
      if(printText) cout << "Total gen event (Poisson-fluc): " << n_gen_event << endl;

      // prepare fake data distributions
      h_m3_data         = (TH1D*)h_m3_sum_toy->Clone("h_m3_data");
      h_m3_data_control = (TH1D*)control_all->Clone("h_m3_data_control");

      delete h_m3_tt_toy;
      delete h_m3_wj_toy;
      delete h_m3_zj_toy;
      delete h_m3_qcd_toy;
      delete h_m3_stop_toy;
      delete h_m3_sum_toy;
            
    }// end MC

    if(n==0) {
      h_m3_data_keep         = (TH1D*)h_m3_data->Clone("h_m3_data_keep");
      h_m3_data_control_keep = (TH1D*)h_m3_data_control->Clone("h_m3_data_control_keep");
    }

    // number of event observed in data (either real or fake)
    const double n_event_obs = h_m3_data->GetEntries();
    

    // export generated m3 histogram (pseudo-data) to RooDataHist object
    RooRealVar m3("m3","m3",0,960); //changed from 1000 to 960
    RooDataHist data("data", "dataset with m3", m3, h_m3_data ); //real or MC data


    //------------------------
    //  2)      f i t
    //------------------------
    // signal and background histograms to make PDFs
    RooDataHist rh_tt("rh_tt",   "tt",  m3, temp_tt);
    RooDataHist rh_wj("rh_wj",   "wj",  m3, temp_wj);
    RooDataHist rh_stop("rh_stop",  "singletop",  m3, temp_stop);
    RooDataHist rh_qcd("rh_qcd", "qcd",  m3, temp_qcd);

    

    // use 0-th order interpolation
    RooHistPdf pdf_tt ("pdf_tt",  "Signal pdf",    m3, rh_tt,  0);
    RooHistPdf pdf_wj ("pdf_wj",  "W+jets pdf",   m3, rh_wj,  0);
    RooHistPdf pdf_qcd( "pdf_qcd", "QCD pdf ",    m3, rh_qcd, 0);
    RooHistPdf pdf_singletop( "pdf_singletop", "single top pdf", m3, rh_stop, 0);


    // create fit model and parameters                       min    max   unit
    // limit on each parameter is 0 -- total number of observed events
    RooRealVar ntt ("ntt",  "number of tt signal events",   -150, 2*n_event_obs, "event");
    RooRealVar nwj ("nwj",  "number of W+jets bgnd events", -150, 2*n_event_obs, "event");

    //RooRealVar ntt ("ntt",  "number of tt signal events",   -150, n_event_obs, "event");
    //RooRealVar nwj ("nwj",  "number of W+jets bgnd events", -150, n_event_obs, "event");

    RooRealVar nqcd("nqcd", "number of QCD bgnd events",  nqcd_est, 0, n_event_obs, "event");
    RooRealVar nstop("nstop", "number of single top bgnd events", nstop_exp, 0, n_event_obs, "event");
    
   
    // add constrain to QCD PDF
    RooConstVar qcd_est("qcd_est", "qcd est", (nqcd_est) ); //mean of gaus 
    RooConstVar qcd_est_unc("qcd_est_unc", "uncertainty of qcd est", nqcd_est_err); //sigma of gaus
    RooGaussian nqcd_constraint("nqcd_constraint","nqcd_constraint", nqcd, qcd_est, qcd_est_unc );

    
    // multiply constraint with QCD PDF
    RooProdPdf pdf_qcd_constraint("pdf_qcd_constraint","constrained QCD pdf",
				  RooArgSet(pdf_qcd, nqcd_constraint) );


    // Mar 09: add constraint to single top PDF using expected value from MC
    RooConstVar singletop_exp("singletop_exp", "singletop exp", ( nstop_exp ) ); //mean of gauss
    RooConstVar singletop_exp_unc("singletop_exp_unc", "uncertainty of singletop exp", nstop_exp_err); //sigma of gaus
    RooGaussian nstop_constraint("nstop_constraint","nstop_constraint", 
				 nstop, singletop_exp, singletop_exp_unc );
    // Multiply constraint with single top PDF                                      
    RooProdPdf pdf_singletop_constraint("pdf_singletop_constraint","constrained single top pdf",
					RooArgSet(pdf_singletop, nstop_constraint) );



    // Create a model with constrained QCD and single top
    //----------------
    if(printText) cout << " ---> Create data model" << endl;
    RooAddPdf model2("model2", "sig+wj+qcd+stop",      
		     RooArgList( pdf_tt, pdf_wj,  pdf_qcd_constraint, pdf_singletop_constraint ),    
		     RooArgList(    ntt,    nwj,                nqcd,            nstop ) ) ;  
    
    


    // Fit model to data  (with internal constraint on nqcd and nstop)
    //--------------------
    if(printText) cout << " ---> Perform RooFit" << endl;
    model2.fitTo(data, Extended(kTRUE), Constrain( RooArgList(nqcd,nstop) ), PrintLevel(-1) );



    if(n==0){
      m3frame1 = m3.frame();
      m3frame2 = m3.frame();
      m3frame3 = m3.frame();
      m3frame4 = m3.frame();
      m3frame5 = m3.frame();
      m3frame6 = m3.frame();
      m3frame7 = m3.frame();
      m3frame8 = m3.frame();
      m3frame9 = m3.frame();
      m3frame10 = m3.frame();
      m3frame11 = m3.frame();      
      //m3frame21 = m3.frame();
      //m3frame22 = m3.frame();
      
      
      //-----------------
      // 15-12-09
      //-----------------
      // error type (Def=Poission for integer histogram bin content)
      // ref: http://root.cern.ch/root/html/RooAbsData.html
      //  RooAbsData::Poission
      //  RooAbsData::SumW2
      //  RooAbsData::None
      //  RooAbsData::Auto      
      RooAbsData::ErrorType etype = RooAbsData::Poisson;  //Def
      if ( !IsData() ) etype =  RooAbsData::SumW2;
      //-----------------


      rh_tt.plotOn( m3frame1, MarkerSize(1), DataError(etype) );
      rh_wj.plotOn( m3frame2, MarkerSize(1), DataError(etype) );
      rh_qcd.plotOn( m3frame3, MarkerSize(1), DataError(etype) );
      
      pdf_tt.plotOn( m3frame4 );
      pdf_wj.plotOn( m3frame5 );
      pdf_qcd.plotOn( m3frame6 );
      pdf_qcd_constraint.plotOn( m3frame7 );
            
      data.plotOn(m3frame8, MarkerSize(1)); //plot pseudo-data
      model2.plotOn(m3frame8); //plot composite pdf (s+b model)
      model2.plotOn(m3frame8,Components(pdf_tt), LineStyle(kDashed), LineColor(kBlue+1));
      model2.plotOn(m3frame8,Components(pdf_wj), LineStyle(kDashed), LineColor(kRed+1));
      model2.plotOn(m3frame8,Components(pdf_qcd_constraint), LineStyle(kDashed), LineColor(kOrange-6));
      
      // single top plots
      rh_stop.plotOn( m3frame9, MarkerSize(1), DataError(etype) );
      pdf_singletop.plotOn( m3frame10 );
      pdf_singletop_constraint.plotOn( m3frame11 );
    }
      
    // QCD pdf
    if(n<5){
      m3_qcd_pdf[n] = m3.frame();
      m3_qcd_pdf[n]->SetNameTitle(Form("m3_qcd_pdf_%d",n+1),Form("m3 pdf: qcd (control region) toy %d",n+1));
      pdf_qcd.plotOn( m3_qcd_pdf[n] );
    }
    
    if(printText) {
      cout << "\n" << n+1 << ") Fit results (with constraint):\n";
      ntt.Print();
      nwj.Print();
      nqcd.Print();
      nstop.Print();
    }
    // record fit results
    double ntt_fit   = ntt.getVal();
    double nwj_fit   = nwj.getVal();
    double nqcd_fit  = nqcd.getVal();
    double nstop_fit  = nstop.getVal();

    double ntt_fiterr  = ntt.getError();
    double nwj_fiterr  = nwj.getError();
    double nqcd_fiterr = nqcd.getError();
    double nstop_fiterr = nstop.getError();

    h_ntt_fiterr ->Fill( ntt_fiterr );
    h_nwj_fiterr ->Fill( nwj_fiterr );
    h_nqcd_fiterr->Fill( nqcd_fiterr );
    h_nstop_fiterr->Fill( nstop_fiterr );

    h_ntt_fit ->Fill( ntt_fit );
    h_nwj_fit ->Fill( nwj_fit );
    h_nqcd_fit->Fill( nqcd_fit );
    h_nstop_fit->Fill( nstop_fit );
  


    // calculate the pull when running on MC
    if(IsData()==false) {
      double ntt_pull  = (ntt_fit - ntt_true) / ntt_fiterr;
      //double nwj_pull  = (nwj_fit - nwj_true) / nwj_fiterr;
      double nwj_pull  = (nwj_fit - nwj_true  ) / nwj_fiterr;
      if(doSystematics=="addZj_fitWithWjShape"){
	nwj_pull = (nwj_fit - nwj_true - nzj_true ) / nwj_fiterr; //include also z+jets
      }
      double nqcd_pull = (nqcd_fit - nqcd_true) / nqcd_fiterr;
      double nstop_pull = (nstop_fit - nstop_true) / nstop_fiterr;

      h_ntt_pull ->Fill( ntt_pull );
      h_nwj_pull ->Fill( nwj_pull );
      h_nqcd_pull->Fill( nqcd_pull );
      h_nstop_pull->Fill( nstop_pull );

      // to exclude problematic fit
      if( ntt_fiterr<100 && nwj_fiterr>20 ) {
	//if(fabs(nwj_pull)<5. && ) {
	h_ntt_goodfit ->Fill( ntt_fit );
	h_nwj_goodfit ->Fill( nwj_fit );
	h_nqcd_goodfit->Fill( nqcd_fit );
	h_nstop_goodfit->Fill( nstop_fit );
	h_ntt_goodpull ->Fill( ntt_pull );
	h_nwj_goodpull ->Fill( nwj_pull );
	h_nqcd_goodpull->Fill( nqcd_pull );
	h_nstop_goodpull->Fill( nstop_pull );
      }
      if(printText) {
	cout << "fit ntt:  " << ntt_fit  << " pm " << ntt_fiterr ;
	cout << " (pull = " << ntt_pull << ")" << endl;
	cout << "fit nwj:  " << nwj_fit  << " pm " << nwj_fiterr ;
	cout << " (pull = " << nwj_pull << ")" << endl;
	cout << "fit nqcd: " << nqcd_fit << " pm " << nqcd_fiterr;
	cout << " (pull = " << nqcd_pull << ")" << endl;
	cout << "fit nstop: " << nstop_fit << " pm " << nstop_fiterr;
	cout << " (pull = " << nstop_pull << ")" << endl;
	cout << endl;

	//record results for the 1st toy MC only
	if(n==0){
	  nFitResults.push_back(ntt_fit);
	  nFitResults.push_back(ntt_fiterr);
	  nFitResults.push_back(nwj_fit);
	  nFitResults.push_back(nwj_fiterr);
	  nFitResults.push_back(nqcd_fit);
	  nFitResults.push_back(nqcd_fiterr);
	  nFitResults.push_back(nstop_fit);
	  nFitResults.push_back(nstop_fiterr);
	}
	// 1 Feb 09: Print result table
	printf("\n               %8s  %6s %18s  %12s\n", "true", "fluc","Gaus constraint","fit");
	printf(" ttbar           %8.2f %6d  %18s  %8.2f +/- %4.2f\n",
	       ntt_true, ntt_fluc, "-", ntt_fit, ntt_fiterr);
	// W only
	printf(" W+jets          %8.2f %6d  %18s  %8.2f +/- %4.2f\n",
	       nwj_true, nwj_fluc, "-", nwj_fit, nwj_fiterr);
	// W+Z
	if(doSystematics=="addZj_fitWithWjShape"){
	  printf(" W/Z+jets      %8.2f %6d  %18s  %8.2f +/- %4.2f\n",
		 nwj_true + nzj_true, nwj_fluc + nzj_fluc, "-", nwj_fit, nwj_fiterr);
	}
	printf(" QCD             %8.2f %6d  %8.2f +/- %4.2f  %8.2f +/- %4.2f\n", 
	       nqcd_true, nqcd_fluc,  nqcd_est, nqcd_est_err, nqcd_fit, nqcd_fiterr);	
	printf(" single top      %8.2f %6d  %8.2f +/- %4.2f  %8.2f +/- %4.2f\n", 
	       nstop_true, nstop_fluc,  nstop_exp, nstop_exp_err, nstop_fit, nstop_fiterr);
	
	printf("\n");
	printf(" W+jets          %8.2f %6d\n", nwj_true, nwj_fluc);
	printf(" Z+jets          %8.2f %6d\n", nzj_true, nzj_fluc);
	printf(" single top      %8.2f %6d\n", nstop_true, nstop_fluc);
      }
    }

    delete h_m3_data;
    delete h_m3_data_control;

  }// end nfit loop

  //add M3 method data to root file

  TH1D *Signal_fitted = new TH1D("Signal_fitted","Signal fitted",2,0,2);
  TH1D *Wjet_fitted   = new TH1D("Wjet_fitted",  "Wjet fitted",  2,0,2);
  TH1D *QCD_fitted    = new TH1D("QCD_fitted",   "QCD fitted",   2,0,2);
  TH1D *SingleTop_fitted = new TH1D("SingleTop_fitted", "Single Top fitted", 2,0,2);

  Signal_fitted->SetBinContent( 1, nFitResults.at(0) );
  Signal_fitted->SetBinContent( 2, nFitResults.at(1) );

  Wjet_fitted->SetBinContent (  1, nFitResults.at(2) );
  Wjet_fitted->SetBinContent(   2, nFitResults.at(3) );

  QCD_fitted->SetBinContent(    1, nFitResults.at(4) );
  QCD_fitted->SetBinContent(    2, nFitResults.at(5) );

  SingleTop_fitted->SetBinContent( 1, nFitResults.at(6) );
  SingleTop_fitted->SetBinContent( 2, nFitResults.at(7) );


  // remove pointer to histograms in input root file
  delete SigEvents;
  delete QCDEvents;
  delete WjetsEvents;
  delete ZjetsEvents;
  delete SingleTopEvents;

  // delete qcdEstimate;

  delete h_m3_tt; //template
  delete h_m3_wj; //template
  delete h_m3_MCdata_control;
  delete control_all;
  if (IsData()==false){
    delete h_m3_qcd; 
    delete h_m3_zj;
    delete h_m3_stop;
  }
 

  if( fdata->GetDirectory("m3_fit")>0 ) {
    //cout << "warning: Directory m3_fit exists. Exit EstimateWjets()."<< endl;
    //return false;
    fdata->rmdir("m3_fit");
  }

  TDirectory *m3_fit = fdata->mkdir("m3_fit");
  m3_fit->cd();


  // write histo to root file
  h_m3_data_keep->Write();
  h_m3_data_control_keep->SetTitle("m3 (data in control reg, used as QCD template)");
  h_m3_data_control_keep->Write();



  m3frame1->SetName("m3frame1");  m3frame1->SetTitle("m3 hist: tt (mc)");  m3frame1->Write(); //rh_tt
  m3frame2->SetName("m3frame2");  m3frame2->SetTitle("m3 hist: wj (mc)");  m3frame2->Write(); //rh_wj
  m3frame3->SetName("m3frame3");  m3frame3->SetTitle("m3 hist: qcd (data control)"); m3frame3->Write(); //rh_data_control
  m3frame4->SetName("m3frame4");  m3frame4->SetTitle("m3 pdf: tt (mc)");  m3frame4->Write(); //pdf_tt
  m3frame5->SetName("m3frame5");  m3frame5->SetTitle("m3 pdf: wj (mc)");  m3frame5->Write(); //pdf_wj
  m3frame6->SetName("m3frame6");  m3frame6->SetTitle("m3 pdf: qcd (control)");   m3frame6->Write();//pdf_qcd_control
  m3frame7->SetName("m3frame7");  m3frame7->SetTitle("m3 pdf: qcd constraint");  m3frame7->Write(); //pdf_qcd_constraint
  m3frame8->SetName("m3frame8");  m3frame8->SetTitle("m3 pseudo-data and fit");  m3frame8->Write();
  m3frame9->SetName("m3frame9");  m3frame9->SetTitle("m3 hist: single top (mc)");  m3frame9->Write();
  m3frame10->SetName("m3frame10");  m3frame10->SetTitle("m3 pdf: single top (mc)");  m3frame10->Write();
  m3frame11->SetName("m3frame11");  m3frame11->SetTitle("m3 pdf: single top constraint (mc)");  m3frame11->Write();
  //  m3frame21->SetName("m3frame21");  m3frame21->SetTitle("m3 hist: data control");  m3frame21->Write();
  //  m3frame22->SetName("m3frame22");  m3frame22->SetTitle("m3 pdf: zj (mc)");  m3frame22->Write(); //pdf_22

  h_ntt_fit->Write();
  h_nwj_fit->Write();
  h_nqcd_fit->Write();
  h_nstop_fit->Write();

  h_ntt_goodfit->Write();
  h_nwj_goodfit->Write();
  h_nqcd_goodfit->Write();
  h_nstop_goodfit->Write();

  h_ntt_fiterr->Write();
  h_nwj_fiterr->Write();
  h_nqcd_fiterr->Write();
  h_nstop_fiterr->Write();

  h_ntt_pull->Write();
  h_nwj_pull->Write();
  h_nqcd_pull->Write();
  h_nstop_pull->Write();

  h_ntt_goodpull->Write();
  h_nwj_goodpull->Write();
  h_nqcd_goodpull->Write();
  h_nstop_goodpull->Write();

  h_ntt_fluc->Write();
  h_nwj_fluc->Write();
  h_nzj_fluc->Write();
  h_nqcd_fluc->Write();
  h_nstop_fluc->Write();


  const int nn = min(nfit,5);
  for(int i=0; i<nn; ++i)  h_m3_sum_toy_keep[i]->Write();
  //  for(int i=0; i<nn; ++i)  h_m3_sum_ctr_toy_keep[i]->Write();
  for(int i=0; i<nn; ++i)  m3_qcd_pdf[i]->Write();

  Signal_fitted->Write();
  Wjet_fitted->Write();
  QCD_fitted->Write();
  SingleTop_fitted->Write();

  model_tt->Write();
  model_wj->Write();
  model_zj->Write();
  model_qcd->Write();
  model_stop->Write();

  temp_tt->Write();
  temp_wj->Write();
  temp_qcd->Write();
  temp_stop->Write();

  fdata->Close();
  fmc->Close();

  cout << "qcd est: " << nqcd_est << ", unc: " << nqcd_est_err << endl;
  cout << "--------------------------------------------------------" << endl;
  cout << "                    W+jets Estimation (end)" << endl;
  cout << "--------------------------------------------------------\n" << endl;
 
  return true;
}
//----------------------------------------------------------------------------------------
//                                 W+jets Estimation (end)
//----------------------------------------------------------------------------------------



void ana::StudySystematics(const string& name,const string& name2){ //TEST
//void ana::StudySystematics(const string name,const string name2){
  doSystematics = name;
  sysSample = name2;
  cout << "Study Systematics"<< endl;
  cout << "doSystematics = " << doSystematics << endl;
  cout << "sysSample = " << sysSample << endl;
}

//----------------------------------------------------------------------------------------
void ana::FillEventCounter(const int istage, const int& ntj, const int& mctype){
  //cout << "istage: " << istage << "   ";
  //cout << "ntj:    " << ntj << "   ";
  //cout << "mctype: " << mctype << endl;
  e_plus_jet[istage][ntj][mctype]++;
  e_plus_jet_weighted[istage][ntj][mctype] += this_weight;
}

//----------------------------------------------------------------------------------------
// Table: Event count as a function of Njet
// 1st call: unweighted
// 2nd call: weighted
void ana::DrawEventPerNjetTable() const {

  static bool first_time(true);
  if(first_time) cout.precision(0);      //unweighted table
  else           cout.precision(myprec); //weighted table

  if(first_time){    
    cout << endl << "%------------------------------------------------------------------";
    cout << endl << "%------------------------------------------------------------------";
    cout << endl << "           Event counts per jet multiplicity";
    cout << endl << "%------------------------------------------------------------------";
    cout << endl << "%------------------------------------------------------------------";
    cout << endl << "\\\\\n";
    cout << endl << "\\begin{tabular}{|l|rrrrr|r|}\n\\hline"; 
    cout << endl << "Unweighted numbers" << endl;
  }
  else {
    //  cout << "\n TEST: uses e_plus_jet (vector)\n" << endl;
    cout << "Normalized for " << intlumi << " pb$^{-1}$" << endl;
  }
  cout << setw(23) 
       << " &" << setw(13) << "0-jet"
       << " &" << setw(13) << "1-jet"
       << " &" << setw(13) << "2-jets"
       << " &" << setw(13) << "3-jets"
       << " &" << setw(13) << "$\\ge$4-jets" 
       << " &" << setw(23) << "Total \\\\\n\\hline";
  cout << endl;
  
  for(int i=0; i<ncutshown; ++i) { //nstage
    printCutStage(i,ve.at(i));
    double total = 0;

    for(int j=0; j<ntjet; ++j){
      double nevent;
           
      if(first_time)  nevent = (double)e_plus_jet[i][j][0];
      else 	      nevent = e_plus_jet_weighted[i][j][0];

      cout << " & " << setw(12) << fixed << nevent;
      total += nevent;
    }
    cout << " & " << setw(13) << fixed << total << " \\\\" << endl;
  }
  cout << "\\hline" << endl;
  if(IsData() || !first_time) cout << "\\end{tabular}\n" << endl;
  first_time = false;
}
//----------------------------------------------------------------------------------------

void ana::DrawMCTypeTable( const string title ) const {

  static bool first_time(true);
  if(first_time) cout.precision(0);      //unweighted table
  else           cout.precision(myprec); //weighted table

  if(first_time) cout << "\\newpage\n"<< endl;
  cout << "\n%------------------------------------------------------------------------\n";
  cout << "%       " << title ;
  cout << "\n%------------------------------------------------------------------------\n";
  
  if(first_time) {
    cout << "\\begin{tabular}{|l|rrrrrr|r|}"<< endl;
    cout << "\\hline"<< endl;
  }
  cout << "\\multicolumn{7}{|l|}";
  if(first_time) cout << "{Actual number of MC events passing selection}";
  else           cout << "{Expected number of events for " << intlumi << "/pb}";
  cout << "\\\\\\hline" << endl;

  cout << "           Cut       "
       << " &" << setw(13) << "\\ttbar{} "
       << " &" << setw(13) << "W+jets "
       << " &" << setw(13) << "Z+jets "
       << " &" << setw(13) << "QCD "
    // << " &" << setw(13) << "VQQ "
       << " &" << setw(13) << "Single Top "
       << " &" << setw(25) << "Total \\\\\n\\hline" << endl;

  int njbegin = 0;
  int p=0;

  for(int i=0; i<ncutshown; ++i){ //loop over cuts

    printCutStage(p, ve2.at(p));
    if(ve2.at(p)==Fourjets)  { njbegin = 4; i--; }
    p++;

    double ntt   = 0;
    double nwj   = 0;
    double nzj   = 0;
    double nqcd  = 0;
    //double nvqq  = 0;
    double nstop = 0;

    for(int j=njbegin; j<ntjet; ++j) {   //sum up jet bins
      for(int k=1;  k<nmctype; ++k) {    //nmc
	double nevent = 0;
	if(first_time) nevent = e_plus_jet[i][j][k];
	else           nevent = e_plus_jet_weighted[i][j][k];

	if( (k>= 1 && k<=10) ||
	    (k>=23 && k<=27) ) ntt += nevent; //sig: 1-10; 23-27(alpgen)
	else if(k==11)         nwj += nevent; //11
	else if(k==12)         nzj += nevent; //12
	else if(k<=18)        nqcd += nevent; //13-18
	//else if(k==19)     nvqq  += nevent; //19
	else if(k<=22)       nstop += nevent; //20-22
      }
    }
    double total_obs = ntt + nwj + nzj + nqcd + nstop; //vqq

    cout << " & " << setw(12) << fixed << ntt;
    cout << " & " << setw(12) << fixed << nwj;
    cout << " & " << setw(12) << fixed << nzj;
    cout << " & " << setw(12) << fixed << nqcd;
    // cout << " & " << setw(12) << fixed << nvqq;
    cout << " & " << setw(12) << fixed << nstop;
    
    //print total column for this cut
    cout << " & " << setw(13) << fixed << total_obs << " \\\\" <<endl;
 
  }// end loop over cut (nstage)
  cout << "\\hline" << endl;
  if(!first_time) cout << "\\end{tabular}\\\\[5mm]" << endl;
  first_time = false;
}
//end DrawMCTypeTable
//---------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
// Print weighted event numbers
//
void ana::DrawSignalBGTable() const {

  cout << "\n%---------------------------------------------------------------------\n";
  cout << "       Expected Signal and Background for " << intlumi << "/pb";
  cout << "\n%---------------------------------------------------------------------\n\n";

  cout << "\\begin{tabular}{|l|r|rr|r|}\\hline" << endl;
  cout << "           Cut       "
       << " &  " << setw(11) << "Total Events"
       << " &  " << setw(14) << "Total Signal (S)"
       << " &  " << setw(16) << "Total Background (B)"
       << " &  " << setw(10) << "S/B \\\\\n\\hline" << endl;

  int njbegin = 0;
  int p=0;

  for(int i=0; i<ncutshown; ++i){ //loop over cuts
    
    double total_sig = 0;
    double total_bkg = 0;
   
    if(ve2.at(p)==Fourjets){ njbegin = 4; i--; }
    /*
    //Total Signal
    for(int k=1; k<11; ++k){ //loop over ttbar mc types (code 1 to 10)
      for(int j=njbegin;j<ntjet;++j){ total_sig += nevent[i][j][k]; }  //sum up jet bins
    }

    //Total BG
    for(int k=11; k<23; ++k){ //loop over all bg mc types (code 11 to 22)
      for(int j=njbegin;j<ntjet;++j) {  total_bkg += nevent[i][j][k]; }  //sum up jet bins
    }
    */
    for(int j=njbegin; j<ntjet; ++j) { // sum up jet bins
      
      for(int k=1; k<11; ++k){ //loop over ttbar mc types (code 1 to 10)
	total_sig += e_plus_jet_weighted[i][j][k];
      }
      for(int k=11; k<23; ++k){ //loop over all bg mc types (code 11 to 22)
	total_bkg += e_plus_jet_weighted[i][j][k];
      }
    }

    //Print cut, S+B, S, B
    printCutStage(p, ve2.at(p));
    p++;
    cout << " & " << setw(13) << fixed << total_sig + total_bkg;
    cout << " & " << setw(13) << fixed << total_sig;
    cout << " & " << setw(13) << fixed << total_bkg;

    //Print Signal-to-background ratio (S/B)
    if( total_bkg > 0 ) 
      cout << " & " << setw(11) << fixed << total_sig/total_bkg;
    else
      cout << " & " << setw(11) << "-" ;
    cout << " \\\\" << endl;
  
  }// end loop over cut (nstage)

  cout << "\\hline\n\\end{tabular}\n" << endl;
}
//end DrawSignalBGTable
//---------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------
// Draw event count table for break down of QCD
//
void ana::DrawQCDTable(const string QCDtitle) const {

  static bool first_time(true);
  if(first_time) cout.precision(0);      //unweighted table
  else           cout.precision(myprec); //weighted table

  if(first_time) cout << "\\newpage\n" << endl;
  cout << "%------------------------------------------------------------------------" << endl;
  cout << "%                   " << QCDtitle << endl;
  cout << "%------------------------------------------------------------------------" << endl;

  if(first_time){
    cout << "\\begin{tabular}{|l|rrrrrr|r|}" << endl;
    cout << "\\hline" << endl;
  }
  cout << "\\multicolumn{8}{|l|}";
  if(first_time) cout << "{Break down of actual number of QCD events passing selection}";
  else           cout << "{Break down of expected QCD events for " << intlumi << "/pb}";
  cout << "\\\\\n\\hline" << endl;

  cout << "           Cut       ";
  cout << " &" << setw(13) << "enri1";
  cout << " &" << setw(13) << "enri2";
  cout << " &" << setw(13) << "enri3";
  cout << " &" << setw(13) << "bce1";
  cout << " &" << setw(13) << "bce2";
  cout << " &" << setw(13) << "bce3";
  cout << " &" << setw(23) << "AllQCD \\\\\\hline" << endl;

  int njbegin = 0;
  int p=0;

  for(int i=0; i<ncutshown; ++i){ //cut stage

    double totalAllQCD = 0;

    //    printCutStage(i, ve.at(i));
    printCutStage(p, ve2.at(p));
    if(ve2.at(p)==Fourjets)  { njbegin = 4; i--; }
    p++;

    for(int k=13; k<19; ++k) { //mctype (QCD): 13-18
      double totalT = 0;
      for(int j=njbegin; j<ntjet; ++j) { //njet
	// totalT += nevent[i][j][k]; 
	if(first_time) totalT += (double)e_plus_jet[i][j][k]; //unweighted
	else           totalT += e_plus_jet_weighted[i][j][k]; //weighted
      }
      totalAllQCD += totalT;
      cout << " & " << setw(12) << totalT;
    }
    cout << " & " << setw(13) << totalAllQCD << " \\\\"<< endl;

  }
  cout << "\\hline" << endl;
  if(!first_time) cout << "\\end{tabular}\\\\[5mm]\n" << endl;
  first_time = false;

}//end DrawQCDTable
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// Draw event count table for break down of Single Top
//
void ana::DrawSingleTopTable( const string title ) const {
    
  static bool first_time(true);
  if(first_time) cout.precision(0);      //unweighted table
  else           cout.precision(myprec); //weighted table
 
  if(first_time) cout << "\\newpage\n" << endl;
  cout << "%------------------------------------------------------------------------" << endl;
  cout << "%                      " << title << endl;
  cout << "%------------------------------------------------------------------------" << endl;

  if(first_time){
    cout << "\\begin{tabular}{|l|rrr|r|}" << endl;
    cout << "\\hline" << endl;
  }
  cout << "\\multicolumn{5}{|l|}";
  if(first_time) cout << "{Break down of actual number of single top events passing selection}";
  else           cout << "{Break down of expected single top events for " << intlumi << "/pb}";
  cout << "\\\\\n\\hline" << endl;

  cout << "          Cut        ";
  cout << " &" << setw(13) << "tW-chan";
  cout << " &" << setw(13) << "t-chan";
  cout << " &" << setw(13) << "s-chan";
  cout << " &" << setw(23) << "AllSingleTop \\\\\\hline" << endl;

  int njbegin = 0;
  int p=0;

  for(int i=0; i<ncutshown; ++i){ //cut stage
    double allSingleTop = 0;
 
   //    printCutStage(i, ve.at(i));
    printCutStage(p, ve2.at(p));
    if(ve2.at(p)==Fourjets)  { njbegin = 4; i--; }
    p++;

    for(int k=20; k<23; ++k) { //mctype (single top): 20-22
      double totalT = 0;
      for(int j=njbegin; j<ntjet; ++j) { //njet
	// totalT += nevent[i][j][k]; 
	if(first_time) totalT += (double)e_plus_jet[i][j][k];
	else           totalT += e_plus_jet_weighted[i][j][k];

      }
      allSingleTop += totalT;
      cout << " & " << setw(12) << totalT;
    }
    cout << " & " << setw(13) << allSingleTop << " \\\\"<< endl;
    
  }
  cout << "\\hline" << endl;
  if(!first_time) cout << "\\end{tabular}\\\\[5mm]" << endl;
  first_time = false;

}//end DrawSingleTopTable
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// Draw event count table for break down of TT + n jet
//
void ana::DrawTTnjetTable( const string title ) const {
  
  static bool first_time(true);
  if(first_time) cout.precision(0);      //unweighted table
  else           cout.precision(myprec); //weighted table

  if(first_time) cout << "\\newpage\n" << endl;
  cout << "%------------------------------------------------------------------------" << endl;
  cout << "%                         " << title << "                                " << endl;
  cout << "%------------------------------------------------------------------------" << endl;

  if(first_time){
    cout << "\\begin{tabular}{|l|rrr|r|}" << endl;
    cout << "\\hline" << endl;
  }
  cout << "\\multicolumn{7}{|l|}";
  if(first_time) cout << "{Break down of actual number of tt+j events passing selection}";
  else           cout << "{Break down of expected tt+j events for " << intlumi << "/pb}";
  cout << "\\\\\n\\hline" << endl;
  cout << "          Cut        ";
  cout << " &" << setw(13) << "tt0j";
  cout << " &" << setw(13) << "tt1j";
  cout << " &" << setw(13) << "tt2j";
  cout << " &" << setw(13) << "tt3j";
  cout << " &" << setw(13) << "tt4j";
  cout << " &" << setw(20) << "Total \\\\\\hline" << endl;

  int njbegin = 0;

  int p=0;

  for(int i=0; i<ncutshown; i++){ //loop over cuts

    printCutStage(p, ve2.at(p));
    if(ve2.at(p)==Fourjets)  { njbegin = 4; i--; }
    p++;

    double sum = 0;

    for(int k=23; k<=27; k++) { //mctype (tt+j): 23-27

      double totalT = 0;
      for(int j=njbegin; j<ntjet; ++j) { //njet
	if(first_time) totalT += (double)e_plus_jet[i][j][k];
	else           totalT += e_plus_jet_weighted[i][j][k];
      }
      sum += totalT;
      cout << " & " << setw(12) << totalT;
    }
    cout << " & " << setw(13) << sum << " \\\\"<< endl;   
  }
  cout << "\\hline" << endl;
  if(!first_time) cout << "\\end{tabular}\\\\[5mm]" << endl;
  first_time = false;

}//end DrawTTnjetTable
//---------------------------------------------------------------------------------------------


// NOTE: print unweighted event numbers
//void ana::DrawSignalAcceptanceTable(const double nevent[][5][23+5], vector<string> ve) const {
void ana::DrawSignalAcceptanceTable( vector<string> ve ) const {

  cout.precision(0); //reset precision

  vector<string> tt;
  tt.push_back("&   Total ");
  tt.push_back("&    evqq ");
  tt.push_back("&    mvqq ");
  tt.push_back("&    tvqq ");
  tt.push_back("&    evev ");
  tt.push_back("&    mvmv ");
  tt.push_back("&    tvtv ");
  tt.push_back("&    evmv ");
  tt.push_back("&    evtv ");
  tt.push_back("&    mvtv ");
  tt.push_back("&    qqqq ");    

     
  for(int h=3; h<5; h++){//jet cut: 3 or 4
    
    TString njets_text = "$\\ge$"; //"$>=$";
    njets_text += h;
    njets_text += "JETS";
      
    if(h==3) cout << endl << "\\newpage";
    cout << endl << "%------------------------------------------------------------------";
    cout << endl << "%------------------------------------------------------------------";
    cout << endl << "      Acceptance for ttbar MC with "<< njets_text<< "  (unweighted)";
    cout << endl << "%------------------------------------------------------------------";
    cout << endl << "%------------------------------------------------------------------";
    cout << endl << "\\begin{tabular}{|l|r|rrrrrrrrrr|}\\hline";

    cout << endl << "           Cut        ";
    for(unsigned int k=0; k<tt.size(); ++k) cout << tt.at(k);
    cout << "\\\\\n\\hline" << endl;
       
    int startnj=0;
    int countstage=0;

    for(int i=0; i<ncutshown; ++i) { //nstage

      printCutStage( countstage, ve.at(i) );
      countstage++;

      int total_tt=0;
      vector<int> total_mode(10,0);

      for(int k=1; k<11; ++k) {
	for(int j=startnj; j<ntjet; ++j) total_mode[k-1] += e_plus_jet[i][j][k];
	total_tt += total_mode[k-1];
      }
      // print first the sum of all ttbar
      cout << " &  " << setw(6) << total_tt ;
      // then each decay mode
      for(int k=0; k<10; ++k) {
 	cout << " &  " << setw(6) << total_mode.at(k) ;
      }
      cout << " \\\\ ";

      if(ve.at(i)=="!MUON") { // extra step for >= n jet cut
	startnj=h;
	ve.at(i) = njets_text;
	i--;
      }
      
      cout << endl;
    }
    cout << "\\hline\n\\end{tabular}"<< endl;
    
    if(h==3) {
      cout << "\\\\[5mm]" << endl;
      ve.at(6) = "!MUON"; //reset
    }
    
  }//cut on 3 or 4 jet

}//end DrawSignalAcceptanceTable
//------------------------------------------------------------------------------------





//====================================================================================
//
//                             START  of  Error  Tables
//
//====================================================================================

void ana::PrintErrorTables(vector<string> ve) const {


  const int mynstage = 11; //up to !DIFFZ

  double e_plus_jet_errors[mynstage][5][24];
  double e_plus_jet_effic[mynstage][5][24];
  double e_plus_jet_effic_unc[mynstage][5][24];

  for(int i=0; i<mynstage; ++i){
    for(int j=0; j<5; ++j){
      for(int k=0; k<24; ++k){
	e_plus_jet_errors[i][j][k]    = 0;
	e_plus_jet_effic[i][j][k]     = 0;
	e_plus_jet_effic_unc[i][j][k] = 0;
      }
    }
  }


  string ttsample = "ttbar";
  if(GetNinit("ttjet")>0) ttsample = "ttjet";
  //if(GetNinit("TTJet")>0) ttsample = "TTJet";


  string wjetSample = "wjet"; //changed from char*
  //  if(nInitialEventMC["WJET"]>0) wjetSample = "WJET";//or whatever the other wjet sample is


  string kIndexmcNames[24] = {"",ttsample,ttsample,ttsample,ttsample,ttsample,ttsample,ttsample,ttsample,ttsample,ttsample,
			      wjetSample,"zjet","enri1","enri2","enri3","bce1","bce2","bce3","vqq","tW","tchan","schan",ttsample};




  for(int i=0; i<mynstage; ++i){
    for(int j=0; j<5; ++j){
      for(int k=1; k<23; ++k){
	
	const long ni = GetNinit( kIndexmcNames[k] );

	//e_plus_jet_effic[i][j][k] = e_plus_jet[i][j][k]/GetNinit(kIndexmcNames[k]);
		
	//e_plus_jet_effic_unc[i][j][k] = sqrt(  e_plus_jet_effic[i][j][k]*(1 - e_plus_jet_effic[i][j][k])/GetNinit(kIndexmcNames[k]) );

	//if(e_plus_jet_effic[i][j][k]==0){e_plus_jet_effic_unc[i][j][k] = GetBayesUncertainty(GetNinit(kIndexmcNames[k]));}

	//----------------------
	// Ask FB to check....
	if ( ni==0 ) continue;

	if( e_plus_jet[i][j][k] > 0 ) {
	  e_plus_jet_effic[i][j][k] = e_plus_jet[i][j][k]/ni  ;
	  e_plus_jet_effic_unc[i][j][k] = sqrt(  e_plus_jet_effic[i][j][k]*(1 - e_plus_jet_effic[i][j][k]) / ni );
	} else {
	  // 0 event pass
	  //	  e_plus_jet_effic[i][j][k] = 0;
	  e_plus_jet_effic_unc[i][j][k] = GetBayesUncertainty( ni );
	}

	//e_plus_jet_errors[i][j][k] = e_plus_jet_effic_unc[i][j][k]*weightMap[kIndexmcNames[k]]*GetNinit(kIndexmcNames[k]);
	e_plus_jet_errors[i][j][k] = e_plus_jet_effic_unc[i][j][k]*GetCrossSection(kIndexmcNames[k])*intlumi; //TL 21-8-09 (ask FB to check)
	//----------------------

      }//end of k loop

      e_plus_jet_effic[i][j][23] = 0;
      for(int k=1;k<11;++k){
	e_plus_jet_effic[i][j][23] += e_plus_jet_effic[i][j][k];
      }
      e_plus_jet_effic_unc[i][j][23] = sqrt(e_plus_jet_effic[i][j][23]*(1-e_plus_jet_effic[i][j][23])/GetNinit(ttsample));
      //e_plus_jet_errors[i][j][23] =  e_plus_jet_effic_unc[i][j][23]*weightMap[ttsample]*GetNinit(ttsample);
      e_plus_jet_errors[i][j][23] =  e_plus_jet_effic_unc[i][j][23]*GetCrossSection(ttsample)*intlumi; //TL 21-8-09 (ask FB to check)

    }
  }


  // QCD sum
  double Sum_Effic2[mynstage][ntjet];
  double Sum_Effic_unc_pos[mynstage][ntjet];
  double Sum_Effic_unc_neg[mynstage][ntjet];
  double Sum_pass_weighted2[mynstage][ntjet];
  // single top sum
  double STopSum_Effic2[mynstage][ntjet];
  double STopSum_Effic_unc_pos[mynstage][ntjet];
  double STopSum_Effic_unc_neg[mynstage][ntjet];
  double STopSum_pass_weighted2[mynstage][ntjet];

  for(int i=0; i<mynstage; ++i){
    for(int j=0; j<5; ++j){
      Sum_Effic2[i][j] = 0.0;
      Sum_Effic_unc_pos[i][j] = 0.0;
      Sum_Effic_unc_neg[i][j] = 0.0;
      Sum_pass_weighted2[i][j] = 0.0;

      STopSum_Effic2[i][j] = 0.0;
      STopSum_Effic_unc_pos[i][j] = 0.0;
      STopSum_Effic_unc_neg[i][j] = 0.0;
      STopSum_pass_weighted2[i][j] = 0.0;
    }
  }


  // QCD
  for(int i=0; i<mynstage; ++i){
    for(int j=0; j<5; ++j){
      for(int k=13; k<19; ++k){
	Sum_pass_weighted2[i][j] +=e_plus_jet_weighted[i][j][k];
	double cer = e_plus_jet_errors[i][j][k];
	Sum_Effic_unc_pos[i][j] += (cer*cer);
	if(e_plus_jet_effic[i][j][k] !=0){Sum_Effic_unc_neg[i][j] += (cer*cer);}
      }
      Sum_Effic_unc_pos[i][j] = sqrt(Sum_Effic_unc_pos[i][j]);
      Sum_Effic_unc_neg[i][j] = sqrt(Sum_Effic_unc_neg[i][j]);
    }
  }

  // Single top
  for(int i=0; i<mynstage; ++i){
    for(int j=0; j<5; ++j){
      for(int k=20; k<23; ++k){
        STopSum_pass_weighted2[i][j] +=e_plus_jet_weighted[i][j][k];
        double cer = e_plus_jet_errors[i][j][k];
        STopSum_Effic_unc_pos[i][j] += (cer*cer);
        if(e_plus_jet_effic[i][j][k] !=0){STopSum_Effic_unc_neg[i][j] += (cer*cer);}
      }
      STopSum_Effic_unc_pos[i][j] = sqrt(STopSum_Effic_unc_pos[i][j]);
      STopSum_Effic_unc_neg[i][j] = sqrt(STopSum_Effic_unc_neg[i][j]);
    }
  }



  ofstream myfile;
  myfile.open("Analyzeroutput.txt",ios::app);

  myfile.setf(ios::fixed,ios::floatfield);


  myfile <<endl;
  myfile <<endl;
  myfile <<endl;
  myfile <<endl;
  myfile <<endl;

  double Allevents[12]; //mynstage
  double AlleventsUncPos[12];
  double AlleventsUncNeg[12];
  double JustSignal[12];
  double JustSignalUnc[12];
  double JustBG[12];
  double JustBGUncPos[12];
  double JustBGUncNeg[12];

  for(int i=0; i<12; ++i){
    Allevents[i] = 0;
    AlleventsUncPos[i] = 0;
    AlleventsUncNeg[i] = 0;
    JustBG[i] = 0;
    JustBGUncPos[i] = 0;
    JustBGUncNeg[i] = 0;
    JustSignal[i] = 0;
    JustSignalUnc[i] = 0;
  }

  myfile.precision(1);

  myfile <<"\\begin{tabular}{|l|ccccc|c|}"<<endl<<"\\hline"<<endl;

  myfile << setw(23) << " Cuts  &"
	 << setw(23) << " \\ttbar{}  &"
	 << setw(23) << " W+jets   &"
	 << setw(23) << " Z+jets   &"
	 << setw(23) << " QCD   &"
    //	 << setw(23) << " VQQ   &"
	 << setw(23) << " Single Top   &"
	 << setw(20) << " Total" << " \\\\ \n\\hline" << endl;
  //    &   W+jets  &  Z+jets  &   QCD   &   VQQ    &  Single Top & Total \\\\ \\hline"<<endl;
  

  // First 7 cuts (up to mu-veto)
  for(int i=0; i<7; ++i){
    double TotalEvents = 0;
    double TotalErrorPos = 0;
    double TotalErrorNeg = 0;

    //myfile << setw(11) << left << ve.at(i) << " "<< right;
    printCutStage(myfile,i,ve.at(i));

    double total = 0;
    double totalerr = 0;

    // signal
    for(int j=0; j<5; ++j){
      for(int k=1; k<11; ++k){total += e_plus_jet_weighted[i][j][k];}
      totalerr += e_plus_jet_errors[i][j][23]*e_plus_jet_errors[i][j][23];
    }
    //myfile << "&" << setw(6) << ScrNum(total) << "$\\pm$" <<setw(6) << left << ScrNum(sqrt(totalerr));
    //myfile << "&" << setw(6) << ScrNum(total) << "$\\pm$" <<setw(6) << left << ScrNum(sqrt(totalerr)) << right;
    printLine(myfile, total, sqrt(totalerr));
    TotalEvents      += total;
    TotalErrorPos    += totalerr;
    TotalErrorNeg    += totalerr;
    JustSignal[i]    += total;
    JustSignalUnc[i] += totalerr;


    // Backgrounds
    // W+jets
    total = 0;
    totalerr = 0;
    for(int j=0; j<5; ++j){
      total    += e_plus_jet_weighted[i][j][11];
      totalerr += e_plus_jet_errors[i][j][11]*e_plus_jet_errors[i][j][11];
    }
    //myfile << "&" << setw(6) << ScrNum(total) << "$\\pm$" << setw(6) << ScrNum( sqrt(totalerr));
    printLine(myfile, total, sqrt(totalerr));
    TotalEvents     += total;
    TotalErrorPos   += totalerr;
    TotalErrorNeg   += totalerr;
    JustBG[i]       += total;
    JustBGUncPos[i] += totalerr;
    JustBGUncNeg[i] += totalerr;

    // Z+jets
    total = 0;
    totalerr = 0;
    for(int j=0; j<5; ++j){
      total += e_plus_jet_weighted[i][j][12];
      totalerr += e_plus_jet_errors[i][j][12]*e_plus_jet_errors[i][j][12];
    }
    //myfile<<"&"<<setw(6) <<ScrNum(total) << "$\\pm$" <<setw(6) << ScrNum( sqrt(totalerr));
    printLine(myfile, total, sqrt(totalerr));
    TotalEvents     += total;
    TotalErrorPos   += totalerr;
    TotalErrorNeg   += totalerr;
    JustBG[i]       += total;
    JustBGUncPos[i] += totalerr;
    JustBGUncNeg[i] += totalerr;

    
    total = 0;
    totalerr = 0;
    double totalNerr = 0;
    for(int j=0; j<5; ++j){
      total += Sum_pass_weighted2[i][j];
      totalerr += Sum_Effic_unc_pos[i][j]*Sum_Effic_unc_pos[i][j];
      totalNerr += Sum_Effic_unc_neg[i][j]*Sum_Effic_unc_neg[i][j];
    }
    //myfile<<"&"<<setw(6) <<ScrNum(total)  << "$\\pm$"<< ScrNum( sqrt(totalerr));
    printLine(myfile, total, sqrt(totalerr));
    TotalEvents     += total;
    TotalErrorPos   += totalerr;
    TotalErrorNeg   += totalNerr;
    JustBG[i]       += total;
    JustBGUncPos[i] += totalerr;
    JustBGUncNeg[i] += totalNerr;


    // VQQ
    /*    
    total = 0;
    totalerr = 0;
    for(int j=0; j<5; ++j){
      total    += e_plus_jet_weighted[i][j][19];
      totalerr += e_plus_jet_errors[i][j][19]*e_plus_jet_errors[i][j][19];
    }
    //myfile<<"&"<<setw(6) <<ScrNum(total) << "$\\pm$" <<setw(6) << ScrNum( sqrt(totalerr));
    printLine(myfile, total, sqrt(totalerr));
    */
    /*
      TotalEvents     += total;
      TotalErrorPos   += totalerr;
      TotalErrorNeg   += totalerr;
      JustBG[i]       += total;
      JustBGUncPos[i] += totalerr;
      JustBGUncNeg[i] += totalerr;
    */

    // Single top
    total = 0;
    totalerr = 0;
    for(int j=0; j<5; ++j){
      total += STopSum_pass_weighted2[i][j];
      totalerr += STopSum_Effic_unc_pos[i][j]*STopSum_Effic_unc_pos[i][j];
    }
    //myfile<<"&"<<setw(6) <<ScrNum(total) << "$\\pm$" <<setw(6) << ScrNum( sqrt(totalerr));
    printLine(myfile, total, sqrt(totalerr));
    TotalEvents     += total;
    TotalErrorPos   += totalerr;
    TotalErrorNeg   += totalerr;
    JustBG[i]       += total;
    JustBGUncPos[i] += totalerr;
    JustBGUncNeg[i] += totalerr;

    // ALL
    //    myfile << "& $"<<ScrNum(TotalEvents)<<"^{+" <<ScrNum( sqrt(TotalErrorPos))<<"}_{-"<<ScrNum(sqrt(TotalErrorNeg)) <<"}$ \\\\"<<endl;
    //myfile << "& "<<ScrNum(TotalEvents)<<"$\\pm$" <<ScrNum( sqrt(TotalErrorPos))<<" \\\\"<<endl;
    printLine(myfile, TotalEvents, sqrt(TotalErrorPos));
    myfile << " \\\\" << endl;
    Allevents[i]       = TotalEvents;
    AlleventsUncPos[i] = TotalErrorPos;
    AlleventsUncNeg[i] = TotalErrorNeg;
  }


  

  // For cut 4mj to DIFFZ
  for(int i=6; i<mynstage; ++i){ //cut (up to !DIFFZ)

    double TotalEvents = 0;
    double TotalErrorPos = 0;
    double TotalErrorNeg = 0;

    if(i==6) { printCutStage(myfile, i, Fourjets); }
    else {     printCutStage(myfile, i, ve.at(i)); }

    double total = 0;
    for(int k=1; k<11; ++k){total += e_plus_jet_weighted[i][4][k];}

    //myfile << "&" << setw(6) << ScrNum(total) << "$\\pm$" << setw(6) << ScrNum( e_plus_jet_errors[i][4][23]);
    printLine(myfile, total, e_plus_jet_errors[i][4][23]);
    TotalErrorPos      += e_plus_jet_errors[i][4][23]*e_plus_jet_errors[i][4][23];
    TotalErrorNeg      += e_plus_jet_errors[i][4][23]*e_plus_jet_errors[i][4][23];
    TotalEvents        += total;
    JustSignal[i+1]    += total;
    JustSignalUnc[i+1] += e_plus_jet_errors[i][4][23]*e_plus_jet_errors[i][4][23];

    
    // W+jets
    //myfile << "&" << setw(6) <<ScrNum(e_plus_jet_weighted[i][4][11]) << "$\\pm$" << setw(6) << ScrNum( e_plus_jet_errors[i][4][11]);
    printLine(myfile, e_plus_jet_weighted[i][4][11], e_plus_jet_errors[i][4][11]);
    TotalErrorPos     += e_plus_jet_errors[i][4][11]*e_plus_jet_errors[i][4][11];
    TotalErrorNeg     += e_plus_jet_errors[i][4][11]*e_plus_jet_errors[i][4][11];
    TotalEvents       += e_plus_jet_weighted[i][4][11];
    JustBG[i+1]       += e_plus_jet_weighted[i][4][11];
    JustBGUncPos[i+1] += e_plus_jet_errors[i][4][11]*e_plus_jet_errors[i][4][11];
    JustBGUncNeg[i+1] += e_plus_jet_errors[i][4][11]*e_plus_jet_errors[i][4][11];
    
    // Z+jets
    //myfile<<"&" << setw(6) << ScrNum(e_plus_jet_weighted[i][4][12]) << "$\\pm$" << setw(6) << ScrNum( e_plus_jet_errors[i][4][12]);
    printLine(myfile, e_plus_jet_weighted[i][4][12], e_plus_jet_errors[i][4][12]);
    TotalErrorPos     += e_plus_jet_errors[i][4][12]*e_plus_jet_errors[i][4][12];
    TotalErrorNeg     += e_plus_jet_errors[i][4][12]*e_plus_jet_errors[i][4][12];
    TotalEvents       += e_plus_jet_weighted[i][4][12];
    JustBG[i+1]       += e_plus_jet_weighted[i][4][12];
    JustBGUncPos[i+1] += e_plus_jet_errors[i][4][12]*e_plus_jet_errors[i][4][12];
    JustBGUncNeg[i+1] += e_plus_jet_errors[i][4][12]*e_plus_jet_errors[i][4][12];


    // QCD
    //myfile << " &  " << setw(6) <<ScrNum(Sum_pass_weighted2[i][4]) <<"$\\pm$" << setw(6) <<ScrNum( Sum_Effic_unc_pos[i][4]);
    printLine(myfile, Sum_pass_weighted2[i][4], Sum_Effic_unc_pos[i][4]);
    TotalErrorPos     += Sum_Effic_unc_pos[i][4]*Sum_Effic_unc_pos[i][4];
    TotalErrorNeg     += Sum_Effic_unc_neg[i][4]*Sum_Effic_unc_neg[i][4];
    TotalEvents       += Sum_pass_weighted2[i][4];
    JustBG[i+1]       += Sum_pass_weighted2[i][4];
    JustBGUncPos[i+1] += Sum_Effic_unc_pos[i][4]*Sum_Effic_unc_pos[i][4];
    JustBGUncNeg[i+1] += Sum_Effic_unc_neg[i][4]*Sum_Effic_unc_neg[i][4];

    // VQQ
    //myfile<<"&"<< setw(6) <<ScrNum(e_plus_jet_weighted[i][4][19]) << "$\\pm$" << setw(6) << ScrNum( e_plus_jet_errors[i][4][19]);
    /*
      printLine(myfile, e_plus_jet_weighted[i][4][19], e_plus_jet_errors[i][4][19]);
    
      TotalErrorPos     += e_plus_jet_errors[i][4][19]*e_plus_jet_errors[i][4][19];
      TotalErrorNeg     += e_plus_jet_errors[i][4][19]*e_plus_jet_errors[i][4][19];
      TotalEvents       += e_plus_jet_weighted[i][4][19];
      JustBG[i+1]       += e_plus_jet_weighted[i][4][19];
      JustBGUncPos[i+1] += e_plus_jet_errors[i][4][19]*e_plus_jet_errors[i][4][19];
      JustBGUncNeg[i+1] += e_plus_jet_errors[i][4][19]*e_plus_jet_errors[i][4][19];
    */

    // Single top
    //myfile << " &  " <<ScrNum(STopSum_pass_weighted2[i][4]) << "$\\pm$" << ScrNum(STopSum_Effic_unc_pos[i][4]);
    printLine(myfile, STopSum_pass_weighted2[i][4], STopSum_Effic_unc_pos[i][4]);
    TotalErrorPos     += STopSum_Effic_unc_pos[i][4];
    TotalErrorNeg     += STopSum_Effic_unc_pos[i][4];
    TotalEvents       += STopSum_pass_weighted2[i][4];
    JustBG[i+1]       += STopSum_pass_weighted2[i][4];
    JustBGUncPos[i+1] += STopSum_Effic_unc_pos[i][4];
    JustBGUncNeg[i+1] += STopSum_Effic_unc_neg[i][4];

    // All
    Allevents[i+1]       += TotalEvents;
    AlleventsUncPos[i+1] += TotalErrorPos;
    AlleventsUncNeg[i+1] += TotalErrorNeg;
    //    myfile << "& "<<ScrNum(TotalEvents)<<"$\\pm$" <<ScrNum( sqrt(TotalErrorPos))<<" \\\\"<<endl;
    printLine(myfile, TotalEvents, sqrt(TotalErrorPos) );
    myfile << " \\\\" << endl;
  }
  myfile <<"  \\hline"<<endl<<"\\end{tabular}"<<endl;
  myfile<<endl<<endl<<endl;



  // Break Down of Single Top (with errors)
  //----------------------------------------
  myfile << "\\begin{tabular}{|l|rrr|r|}" << endl;
  myfile << "\\hline" << endl;
  myfile << "\\multicolumn{5}{|l|}";
  myfile << "{Break down of expected single top events for " << (int)intlumi << "/pb}";
  myfile << "\\\\\n\\hline" << endl;

  myfile << "            Cut      ";
  myfile << " &" << setw(21) << "tW-chan  ";
  myfile << " &" << setw(21) << "t-chan  ";
  myfile << " &" << setw(21) << "s-chan  ";
  myfile << " &" << setw(29) << "AllSingleTop \\\\\\hline" << endl;


  int njbegin = 0;

  for(int i=0; i<11; ++i){ //cut stage (up to HT)
    double allSingleTop = 0;
    double allSingleTopEr = 0;
    printCutStage(myfile, i, ve.at(i));
    //   myfile << " Stage " << setw(2) << i << " " << setw(11) << left << ve.at(i) << right;
    for(int k=20; k<23; ++k) { //mctype (single top): 20-22
      double totalT = 0;
      double totalEr = 0;
      for(int j=njbegin; j<ntjet; ++j) { //njet
        totalT += e_plus_jet_weighted[i][j][k];
	totalEr += e_plus_jet_errors[i][j][k]*e_plus_jet_errors[i][j][k];
      }
      totalEr = sqrt(totalEr);
      allSingleTop += totalT;
     
      //myfile << " & " << setw(12) << ScrNum(totalT) <<"$\\pm$"<< ScrNum(totalEr);
      printLine(myfile, totalT, totalEr);
    }
    for(int j=njbegin; j<ntjet; ++j) {allSingleTopEr +=STopSum_Effic_unc_pos[i][j]*STopSum_Effic_unc_pos[i][j];}
    allSingleTopEr = sqrt(allSingleTopEr);

    myfile << " & " << setw(10) << ScrNum(allSingleTop) << "$\\pm$" << setw(4) << left << ScrNum(allSingleTopEr) 
	   << right << " \\\\"<< endl;

    // insert >=4j cut
    if(ve.at(i)=="!MUON"){ 
      njbegin = 4; ve.at(i) = Fourjets; i--;
    }
  }

  myfile << "\\hline" << endl;
  myfile << "\\end{tabular}\\\\[5mm]" << endl;
  myfile << endl << endl;
  ve.at(6) = "!MUON"; //restore


  //BEGIN QCD BREAKDOWN TABLE
  //-------------------------
  myfile << "\\begin{tabular}{|l|cccccc|c|}" << endl;
  myfile << "\\hline" << endl;
  myfile << "\\multicolumn{8}{|l|}";
  myfile << "{Break down of expected QCD events for " << intlumi << "/pb}";
  myfile << "\\\\\n\\hline" << endl;
  
  myfile << "            Cut      ";
  myfile << " &" << setw(24) << "enri1 ";
  myfile << " &" << setw(24) << "enri2 ";
  myfile << " &" << setw(24) << "enri3 ";
  myfile << " &" << setw(24) << "bce1 ";
  myfile << " &" << setw(24) << "bce2 ";
  myfile << " &" << setw(24) << "bce3 ";
  myfile << " &" << setw(27) << "AllQCD \\\\\\hline" << endl;

  //ntjet = 5;
  njbegin = 0;

  for(int i=0; i<11; ++i){ //cut stage (up to HT)
    double totalAllQCD = 0;
    double totalAllQCDErPos = 0;
    double totalAllQCDErNeg = 0;
    //myfile << "" << setw(2) << i << " " << setw(11) << left << ve.at(i) << right;
    printCutStage(myfile,i,ve.at(i));

    for(int k=13; k<19; ++k) { //mctype (QCD): 13-18
      double totalT = 0;
      double totalEr = 0;
      for(int j=njbegin; j<ntjet; ++j) { //njet
        totalT  += e_plus_jet_weighted[i][j][k];
        totalEr += e_plus_jet_errors[i][j][k]*e_plus_jet_errors[i][j][k];
      }
      totalEr = sqrt(totalEr);
      totalAllQCD += totalT;

      myfile << " &" << setw(12) << ScrNum(totalT);
      if(totalT !=0) myfile <<"$\\pm$"<< setw(7) << left << ScrNum(totalEr) << right;
      else{myfile <<"$<$"<< setw(9) << left << ScrNum(totalEr) << right;}
    }
    for(int j=njbegin; j<ntjet; ++j) {
      totalAllQCDErPos +=Sum_Effic_unc_pos[i][j]*Sum_Effic_unc_pos[i][j];
      totalAllQCDErNeg +=Sum_Effic_unc_neg[i][j]*Sum_Effic_unc_neg[i][j];
    }
    
    totalAllQCDErPos = sqrt(totalAllQCDErPos);
    totalAllQCDErNeg = sqrt(totalAllQCDErNeg);
    //    myfile << " & " << setw(13) << "$"<<ScrNum(totalAllQCD) <<"^{+"<<ScrNum(totalAllQCDErPos)<<"}_{-"<<ScrNum(totalAllQCDErNeg) <<"}$ \\\\"<< endl;
    myfile << " & " << setw(13) <<ScrNum(totalAllQCD) <<"$\\pm$"<< setw(8) << left<< ScrNum(totalAllQCDErPos) << right << " \\\\"<< endl;
    
    if(ve.at(i)=="!MUON"){
      njbegin = 4;  ve.at(i) = Fourjets;  i--;
    }
  }
  myfile << "\\hline" << endl;
  myfile << "\\end{tabular}\\\\[5mm]\n" << endl<<endl;

  ve.at(6)="!MUON";//restore



  myfile << "\n%---------------------------------------------------------------------\n";
  myfile << "       Expected Signal and Background for " << intlumi << "/pb";
  myfile << "\n%---------------------------------------------------------------------\n\n";
  myfile << "\\begin{tabular}{|l|r|rr|r|}\\hline" << endl;
  myfile << "          Cut        "
	 << " &  " << setw(30) << "Total Events "
	 << " &  " << setw(24) << "Total Signal (S) "
	 << " &  " << setw(26) << "Total Background (B) "
	 << " &  " << setw(20) << "S/B   \\\\\n\\hline" << endl;

  
  // Insert >=4j cut after muon-veto
  const int muVeto_pos = 7;
  ve.insert( ve.begin()+muVeto_pos, Fourjets);

  for(int i=0; i <= 11; ++i){ //loop over cuts (up to DIFFZ)

    printCutStage(myfile, i, ve.at(i));
    
    //myfile << " & " << setw(12) << ScrNum(Allevents[i]) << "$^{+"<<ScrNum(sqrt(AlleventsUncPos[i])) << "}_{-"<< ScrNum(sqrt(AlleventsUncNeg[i]))<
    //myfile << " & " << setw(12) << ScrNum(JustSignal[i])<<"$\\pm$"<<ScrNum(sqrt(JustSignalUnc[i]));
    //myfile << " & " << setw(12) << ScrNum(JustBG[i])<< "$^{+"<<ScrNum(sqrt(JustBGUncPos[i]))  << "}_{-"<<ScrNum( sqrt(JustBGUncNeg[i]))<<"}$";

    myfile << " & " << setw(12) << right << ScrNum(Allevents[i]) 
	   << setw(24) << left << Form( "$^{+%s}_{-%s}$",ScrNum(sqrt(AlleventsUncPos[i])).c_str(), ScrNum(sqrt(AlleventsUncNeg[i])).c_str() ) << left;
    myfile << " & " << setw(10) << right << ScrNum(JustSignal[i])
	   << setw(10) << left << "$\\pm$"+ScrNum(sqrt(JustSignalUnc[i])) ;
    myfile << " & " << setw(12) << right << ScrNum(JustBG[i])
	   << setw(24) << left << Form( "$^{+%s}_{-%s}$", ScrNum(sqrt(JustBGUncPos[i])).c_str(),ScrNum( sqrt(JustBGUncNeg[i])).c_str() ) << right;
						
    //Print Signal-to-background ratio (S/B)
    if( JustBG[i] > 0 ) {
      double errSBGPos = sqrt(JustSignalUnc[i]/(JustSignal[i]*JustSignal[i]) + JustBGUncPos[i]/(JustBG[i]*JustBG[i]) );
      myfile << " & " << setw(11) << ScrNum(JustSignal[i]/JustBG[i])<< "$\\pm$" <<ScrNum(errSBGPos*JustSignal[i]/JustBG[i]);
    } else {
      myfile << " & " << setw(12) << "-" ;
    }
    myfile << " \\\\" << endl;
  
  }// end loop over cut (nstage)
  myfile << "\\hline\n\\end{tabular}\n" << endl<<endl<<endl;

  // restore (take out 4mj)
  //  ve.erase( ve.begin()+muVeto_pos);


  myfile << "%---------------------------------------"<< endl;
  myfile << "            NjetVcut table "<< endl;
  myfile << "%---------------------------------------\n"<< endl;
  PrintError_NjetVcut(myfile, e_plus_jet_errors );

  myfile.close();

}//end PrintErrorTables
//----------------------------------------------------------------------------------


double ana::GetBayesUncertainty(int Ninitial) const{

  if(Ninitial==0) return 0; //NEW, added by TL (Ask FB to check)
  
  TH1D *h_total = new TH1D("total","events that are incident",1,0,1);
  h_total->SetBinContent( 1, Ninitial );
  TH1D *h_pass = new TH1D("pass","events that pass",1,0,1);
  h_pass->SetBinContent( 1, 0 );

  TGraphAsymmErrors *tge = new TGraphAsymmErrors();
  tge->BayesDivide(h_pass,h_total);
  delete h_total; delete h_pass;
  double error1 = tge->GetErrorYhigh(0);
  delete tge;
  return error1;

}// end GetBayesUncertainty
//---------------------------------------------------------------------

void ana::printLine(ofstream &myfile, double num, double err) const {
  myfile << " &" << setw(9) << ScrNum(num) << "$\\pm$" << setw(7) << left << ScrNum(err) << right;
}

string ana::ScrNum(double num) const{

  string Number = "";
  if(ScientificNotation){
    if(num/10000.0 > 1){
      Number = Form("%.1E",num);
    }
    else{Number = Form("%.1f",num);}
  }
  else{
    Number = Form("%.1f",num);
  }
  return Number;
}
//----------------------------------------------------------------------------------


void ana::PrintError_NjetVcut(ofstream& myfile, const double e_plus_jet_errors[][5][24] ) const {

  myfile<<"\\begin{tabular}{|c|ccccc|c|}"<<endl;
  myfile<<"\\hline"<<endl;
    
  myfile << setw(21) << "Cut  "
	 << " &" << setw(26) << "0-jet "
	 << " &" << setw(26) << "1-jet "
	 << " &" << setw(26) << "2-jets "
	 << " &" << setw(26) << "3-jets "
	 << " &" << setw(26) << "$\\ge$4-jets "
	 << " &" << setw(36) << "Total   \\\\\n\\hline\n";

  for(int i=0; i<11; ++i) { //nstage

    printCutStage(myfile,i,ve.at(i));

    double total = 0;
    double PosError = 0;
    double NegError = 0;
    for(int j=0; j<ntjet; ++j){
      double myTotal = 0;
      double myPosErr = 0;
      double myNegError = 0;
      for(int k=1;k<23;++k){
	if(k==19) continue;//skip vqq
	myTotal += e_plus_jet_weighted[i][j][k];
	myPosErr += e_plus_jet_errors[i][j][k]*e_plus_jet_errors[i][j][k];
	if(e_plus_jet_weighted[i][j][k] != 0){ myNegError += e_plus_jet_errors[i][j][k]*e_plus_jet_errors[i][j][k];}
      }
      myfile << " & " << setw(12) << ScrNum(myTotal) << "$\\pm$"<< setw(8) << left << ScrNum(sqrt(myPosErr)) << right;
      total    += myTotal;
      PosError += myPosErr;
      NegError += myNegError;
    }//loop jets
    myfile << " & " << setw(13) << fixed << ScrNum(total) << "$\\pm$"<< setw(8) << left << ScrNum(sqrt(PosError)) <<" \\\\ \n" << right;
  }//loop cuts
  myfile << "\n\\hline\n\\end{tabular}" << endl;         

}//end PrintError_NjetVcut
//----------------------------------------------------------------------------------
//==================================================================================
//
//                          END  of  Error  tables
//
//==================================================================================





void ana::SetHistoLabelCutNjet( TH2D *this_njetsVcuts, vector<string>& ve ) const {

  // ve.size() = nstage
  for (size_t i=0; i < ve.size(); ++i) {
    this_njetsVcuts->GetYaxis()->SetBinLabel(i+1, ve.at( ve.size()-i-1 ).c_str() ); //i=0, bin=1, stage=13, cut=-1btag
  }
  this_njetsVcuts->GetXaxis()->SetBinLabel(1,"0-jet");
  this_njetsVcuts->GetXaxis()->SetBinLabel(2,"1-jet");
  this_njetsVcuts->GetXaxis()->SetBinLabel(3,"2-jets");
  this_njetsVcuts->GetXaxis()->SetBinLabel(4,"3-jets");
  this_njetsVcuts->GetXaxis()->SetBinLabel(5,"#geq4-jets");
  this_njetsVcuts->GetXaxis()->SetBinLabel(6,"Total");
  this_njetsVcuts->SetLabelSize(0.07,"x");
  this_njetsVcuts->GetYaxis()->SetLabelOffset(0.0); 
  this_njetsVcuts->SetOption("text");
  this_njetsVcuts->SetMarkerSize(2.0);
  this_njetsVcuts->SetStats(0);
  
}//end SetHistoLabelCutNjet
//---------------------------------------------------------------------------------------------

//void ana::SetHistoLabelEleID( TH1F *eid[] ) const {
void ana::SetHistoLabelEleID( vector<TH1*>& eid ) const {
  
  //  int nhisto = 1; //data
  //  if(!IsData()) nhisto = 16; //MC
  // ntype = nhisto
  for (int i=0; i<ntype; ++i){
    if( i>0 && is_mc_present(i)==false ) continue;
    eid[i]->GetXaxis()->SetBinLabel(1,"allEle");
    eid[i]->GetXaxis()->SetBinLabel(2,"ET/Eta");
    eid[i]->GetXaxis()->SetBinLabel(3,"d0");
    eid[i]->GetXaxis()->SetBinLabel(4,"eidLoose");
    eid[i]->GetXaxis()->SetBinLabel(5,"eidTight");
    eid[i]->GetXaxis()->SetBinLabel(6,"eidRobustLoose");
    eid[i]->GetXaxis()->SetBinLabel(7,"eidRobustTight");    
    eid[i]->SetLabelSize(0.04,"x");
    eid[i]->SetOption("hist text0");
    eid[i]->SetMarkerSize(1.5);
    eid[i]->SetStats(0);
  }
}//end SetHistoLabelEleID
//---------------------------------------------------------------------------------------------

bool ana::is_mc_present( const int& code ) const {

  switch ( code ) {
  case  1 : return mc_sample_has_ttbar; break;
  case  2 : return mc_sample_has_QCD;   break;
  case  3 : return mc_sample_has_enri1; break;
  case  4 : return mc_sample_has_enri2; break;
  case  5 : return mc_sample_has_enri3; break;
  case  6 : return mc_sample_has_bce1;  break;
  case  7 : return mc_sample_has_bce2;  break;
  case  8 : return mc_sample_has_bce3;  break;
  case  9 : return mc_sample_has_Wjet;  break;
  case 10 : return mc_sample_has_Zjet;  break;
  case 11 : return mc_sample_has_VQQ;   break;
  case 12 : return mc_sample_has_singleTop; break;
  case 13 : return mc_sample_has_tW;    break;
  case 14 : return mc_sample_has_tchan; break;
  case 15 : return mc_sample_has_schan; break;
  default : return false; break;
  }
}//end is_mc_present
//---------------------------------------------------------------------------------------------

//--------- compute (signed) d0 corrected w.r.t Beam Spot -------------------------------------
float ana::compute_d0 ( const string& lepton, const int& i ) const {
  
  float vx = 0;
  float vy = 0;
  float px = 0;
  float py = 0;

  if (lepton=="electron") {
    vx = els_vx->at(i);
    vy = els_vy->at(i);
    px = els_vpx->at(i);
    py = els_vpy->at(i);
  }
  else if (lepton=="muon") {
    vx = mus_tk_vx->at(i);
    vy = mus_tk_vy->at(i);
    px = mus_tk_px->at(i);
    py = mus_tk_py->at(i);
  }
  else { //track  
    vx = tracks_vx->at(i);
    vy = tracks_vy->at(i);
    px = tracks_px->at(i);
    py = tracks_py->at(i);
  }

  float refx = beamSpot_x->at(0);
  float refy = beamSpot_y->at(0);
  
  float d0_corrected = - (-(vx - refx)*py + (vy - refy)*px) / TMath::Hypot(px,py) ; //d0 = -dxy 
  //cout << "d0 corrected (px,py):  " << d0_corrected << endl;

  return d0_corrected;
}//end compute_d0
//---------------------------------------------------------------------------------------------

//--------- compute mT(W) given electron 2-vector and MET -------------------------------------
float ana::compute_mtw ( const TVector2& e, const TVector2& miss ) const {
  
  float leppx = e.Px();
  float leppy = e.Py();
  float lepet = e.Mod();
  //cout << "lep: px=" << leppx << "  py=" << leppy << " et=" << lepet << endl;
  
  float mspx = miss.Px();
  float mspy = miss.Py();
  float mset = miss.Mod();

  //  float mset = TMath::Hypot(mspx,mspy);
  //cout << "met: px=" << mspx << "  py=" << mspy << "  mset=" << mset << endl;
  
  // need ET, px, py
  // leppt = lepet
  //double leppx = (*electrons)[which_e].px();
  //double leppy = (*electrons)[which_e].py();
  //double mspx = (*mets)[0].px();
  //double mspy = (*mets)[0].py();
  /*      
	  cout << "lepet "  << electron_ET << endl;
	  cout << "mspx "  << mspx << endl;
	  cout << "mspy "  << mspy << endl;
	  cout << "mset "  << missing_ET << endl;
  */
  float z1 = pow( ( lepet + mset ), 2);
  float z2 = pow( (leppx + mspx), 2 ) + pow( (leppy + mspy), 2 );
  //  double z2 = 2*(leppx*mspx + leppy*mspy) ;                                               
  
  float z3 = z1 - z2;
  float mtlm = -1.0;
  //cout << "z1/z2/z3: "  << z1 << " / " << z2 << " / " << z3 << endl;                        
  
  if (z3 > 0) {
    mtlm = sqrt( z3 );
  }
  //cout << "mtlm: " << mtlm << endl;                    
  return mtlm;
  
}//end compute_mtw
//---------------------------------------------------------------------------------------------

void ana::PrintGenParticles() const {
  cout << setfill('-') << setw(105) << "" << setfill(' ') << endl;
  cout << setw(3) << ""
       << setw(8) << "PDG id"
       << setw(8) << "Status"
       << setw(10) << "Parent id"
       << setw(15) << "Mass"
       << setw(12) << "Px"
       << setw(12) << "Py"
       << setw(12) << "Pz"
       << setw(12) << "Pt"
       << setw(12) << "Energy"
       << endl;
  cout << setfill('-') << setw(105) << "" << setfill(' ') << endl;
  for(unsigned int i=0; i<Nmc_doc; ++i) { 
    cout << setw(3) << i 
	 << setw(8) << mc_doc_id->at(i) 
	 << setw(8) << mc_doc_status->at(i)
	 << setw(10) << mc_doc_mother_id->at(i) 
	 << setw(15) << mc_doc_mass->at(i) 
	 << setw(12) << mc_doc_px->at(i) 
	 << setw(12) << mc_doc_py->at(i) 
	 << setw(12) << mc_doc_pz->at(i) 
	 << setw(12) << mc_doc_pt->at(i)  
	 << setw(12) << mc_doc_energy->at(i)
	 << endl;
  }
  cout << setfill('-') << setw(105) << "" << setfill(' ') << endl;
}//end PrintGenParticles
//---------------------------------------------------------------------------------------------

void ana::printCutStage( const int& i, const string& cut ) const {
  cout << " Stage " << setw(2) << i << " " << setw(11) << left << cut << right;  
}
//---------------------------------------------------------------------------------------------

void ana::printCutStage( ofstream& os, const int& i, const string& cut ) const {
  os << " Stage " << setw(2) << i << " " << setw(11) << left << cut << right;  
}
//---------------------------------------------------------------------------------------------

// *** NOT READY FOR USED ***
/*
string ana::CheckEventTypeFromMcTruth() const {
  bool found_top = false;
  bool found_W   = false;
  bool found_Z   = false;
  // Look only at "documentation lines" (in Pythia), ie status=3 particles
  for(unsigned int i=0; i<Nmc_doc; ++i) {
    //cout << "mc " << i << ",  st "<< mc_doc_status->at(i) << endl;
    if ( mc_doc_status->at(i) != 3 ) break;
    if      ( fabs(mc_doc_id->at(i)) == 6  ) found_top = true;
    else if ( fabs(mc_doc_id->at(i)) == 24 ) found_W = true; 
    else if ( fabs(mc_doc_id->at(i)) == 22 ) found_Z = true; 
  }
  if(found_top) return "ttbar";
  if(found_W) return "wenu";
  if(found_Z) return "zee";
  return "none";
}
*/
//---------------------------------------------------------------------------------------------

float ana::getRelIso(const unsigned int& i) const {
  return (els_dr04EcalRecHitSumEt->at(i) + els_dr04HcalTowerSumEt->at(i) + els_tIso->at(i))/els_et->at(i);
}
//---------------------------------------------------------------------------------------------

bool ana::passEleID(const unsigned int& i) const {
  if ( EleID()==robustTight ) return (els_robustTightId->at(i) > 0);
  if ( EleID()==robustLoose ) return (els_robustLooseId->at(i) > 0);
  if ( EleID()==tight )       return (els_tightId->at(i) > 0);
  if ( EleID()==loose )       return (els_looseId->at(i) > 0);
  if ( EleID()==none )        return true;
  return false;
}
//---------------------------------------------------------------------------------------------

string ana::printEleID() const {
  switch ( EleID() ) {
  case robustTight: return "robustTight";
  case robustLoose: return "robustLoose";
  case tight:       return "tight";
  case loose:       return "loose"; 
  case none:        return "none"; 
  default:          return "robustTight (Default)";
  }
}
//---------------------------------------------------------------------------------------------

bool ana::passHLT() const {
  if ( HLTBit=="HLT_Ele15_LW_L1R" ) return (bool)HLT_Ele15_LW_L1R;
  if ( HLTBit=="HLT_Ele15_SW_L1R" ) return (bool)HLT_Ele15_SW_L1R;
  return false;
}
//---------------------------------------------------------------------------------------------

string ana::printTimeNow() const {
  TDatime now;
  return now.AsSQLString();
}
//---------------------------------------------------------------------------------------------

void ana::DoBTagging(const vector<TLorentzVector>& electrons){
  // .. TO DO ...
  // for each btag algo, find nbjets, fill histo
  // pick one algo for report in cut table.
  m_nbtag_TCHE = 0;
  m_nbtag_TCHP = 0;
  m_nbtag_SSV  = 0;
  //int nbtag_softMu = 0;
     
  // 16 Feb 2010: btag	   
  //--------------------
  // jets_btag_TC_highPur
  // jets_btag_TC_highEff
  // jets_btag_TC_secVertex = simpleSecondaryVertex
  // jets_btag_TC_softMuon
  // https://twiki.cern.ch/twiki/bin/view/CMS/BtagOctober09ExerciseDetails#Working_Point_definition 
  // Operating points, see Francisco's talk: http://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=69434
  // TCHE M = 3.3 (7TeV), 3.99 (10TeV)
  // TCHP M = 1.93 (7TeV), 2.17(10TeV)
  // SSV M  = 1.74 (710TeV) 2.02 (10TeV)
  //-------------------------
  float discut_TCHE = 3.3; //Medium (mistagEff=1%) 7 TeV
  float discut_TCHP = 1.93; //Medium
  float discut_SSV = 1.74;// Medium
  if(m_LHCEnergyInTeV==10){ 
    discut_TCHE = 3.99;
    discut_TCHP = 2.17;
    discut_SSV  = 2.02;
  }
 
  // Loop over all jets (Default AK5 Calo Jets)
  for(unsigned int i=0; i< Njets; ++i){
 
   if ( jets_pt->at(i) > JET_PTCUT   &&
	 fabs(jets_eta->at(i)) < 2.4  &&
	 jets_emf->at(i) > 0.01 )  {

      TLorentzVector thisjet(jets_px->at(i),jets_py->at(i),jets_pz->at(i),jets_energy->at(i));

      if( jetNotNearElectron(thisjet, electrons) ) {	 // Require DR(e,j)>0.3
	if(jets_btag_TC_highEff->at(i) > discut_TCHE) m_nbtag_TCHE++;
	if(jets_btag_TC_highPur->at(i) > discut_TCHP) m_nbtag_TCHP++;
	if(jets_btag_secVertex->at(i)  > discut_SSV)  m_nbtag_SSV++;
      }//DR
    }//quality cuts
  }
  /*
      fillHistoDataAndMC( h_nbtag_TCHE, nbtag_TCHE, this_weight );
      fillHistoDataAndMC( h_nbtag_TCHP, nbtag_TCHP, this_weight );
      fillHistoDataAndMC( h_nbtag_SSV,  nbtag_SSV,  this_weight );
     
      h_nbtag_TCHE->Fill(nbtag_TCHE); //track counting
      h_nbtag_TCHP->Fill(nbtag_TCHP);
      h_nbtag_SSV->Fill(nbtag_SSW);
  */
}
//---------------------------------------------------------------------------------------------

bool ana::jetNotNearElectron(const TLorentzVector& j, const vector<TLorentzVector>& e) const {
  for(unsigned int i=0; i<e.size(); ++i){
    if( j.DeltaR(e[i]) < 0.3 ) return false; //disregard this jet
  }
  return true;
}
//-- eof ------------------------------------------------------------------------------------------
