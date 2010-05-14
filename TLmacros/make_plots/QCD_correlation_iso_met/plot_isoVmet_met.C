//
// 18 Mar 2010
// cut: at each level of NES. Barrel/Endcap/Both.
// plot: 
// 1) 2D plot of iso:met for 1-4j, allj
// 2) project met from 2D: noISO, loISO, hiISO.
// permutation: nj: 1-4j
//              proc: enri1,2,3, bce1,2,3, QCD
//
// NES stage   2D   <------- MET ------>
//       1     2D   noISO iso<0.1 iso>0.1
//       2      x     x      x      x
//       3
//---------------------------------------------
/*
#include <iostream>
#include <string>
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TPaveLabel.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TROOT.h"

using namespace std;
*/

string etaRegion="both"; //barrel,endcap,both
string metType = "t1";

string metcut_label[2] = {"t1met>30"," #mumet>30"};
int metLabel = 0; //D=t1

string nj="1j";
bool plot_2D_unweighted = false;

const int nLevel = 11;
const int nh = 3*nLevel;

bool debug = 0;


//--------------------------------------------------------------------------------
// note: 
// plot_isoVmet_met() -> plot() -> plot_njet() -> do_plot_isoVmet()
//--------------------------------------------------------------------------------
void plot_isoVmet_met(string etaRegion_user="both",string metType_user="t1"){
  etaRegion = etaRegion_user;
  metType = metType_user;
  cout << "etaRegion :  "<< etaRegion << endl;
  cout << "met: " << metType << endl;
  plot();
}


void plot() {
  if(metType=="t1")  metLabel = 0; //t1MET
  else               metLabel = 1; //muMET

  plot_nj("1j");  
  plot_nj("2j");
  plot_nj("3j");
  plot_nj("4mj");
  plot_nj("allj");
}


void plot_nj(string nj_user="1j") {
  nj = nj_user;

  do_plot_isoVmet("enri1");
  do_plot_isoVmet("enri2");
  do_plot_isoVmet("enri3");
  do_plot_isoVmet("bce1");
  do_plot_isoVmet("bce2");
  do_plot_isoVmet("bce3");   
  do_plot_isoVmet("QCD"); //total QCD

  do_plot_isoVmet("ttbar");
}


void do_plot_isoVmet(string sel_user = "enri1") {

  string sel = sel_user;//"enri1";
  cout << "etaRegion, met, nj, sel:  " << etaRegion << ", " << metType << ", " << nj << ", "<< sel << endl;


  TFile f("../../test_mc_mixture.root");


  int    rb      = 10; // rebin of reliso (1D)
  double met_max = 50;
  double iso_max = 1.39; //1D

  //  gStyle->SetStatX(0.91);
  //  gROOT->SetStyle("Default");
  gStyle->SetPadTopMargin(0.13);
  gStyle->SetPadRightMargin(0.125);
  gStyle->SetOptStat(1110); //n,mean,rms
  gStyle->SetNumberContours(99);

  

  string opname = Form("isoVmet_met_%s_%s", nj.c_str(), sel.c_str());
  if(etaRegion!="both")
    opname = Form("isoVmet_met_%s_%s_%s", etaRegion, nj.c_str(), sel.c_str());
  //  if(sel=="QCD") opname = Form("iso_met_%s_qcd_total",nj.c_str()); 
  //  int nLevel = 7+4;
  //  int nh = 3*nLevel;

  if(debug) cout << "defining histo" << endl;

  TH1D *h[nh]; //projected 1D
  TH2F *h2[nLevel];

  const string var = "QCD_estimation/NES/QCDest_isoVmet_NES";

  //    "QCD_estimation/NES/QCDest_isoVmet_NES_L1",
  //    "QCD_estimation/NES/QCDest_isoVmet_NES_barrel_L1",
  //    "QCD_estimation/NES/QCDest_isoVmet_NES_endcap_L1",

  //    "QCD_estimation/NES/QCDest_isoVmet_NES_uw_L1",
  //    "QCD_estimation/NES/QCDest_isoVmet_NES_uw_barrel_L1",
  //    "QCD_estimation/NES/QCDest_isoVmet_NES_uw_endcap_L1",


  const string level[nLevel] ={ 
    "_L1",
    "_L1b",
    "_L1c",
    "_L1d1",
    "_L1d2",
    "_L1d3",
    "_L1d4",
    "_L1d5",
    "_L2",
    "_L3",
    "_L4"
  };



  int npadx = 4;
  int npady = nLevel;
  int npad  = 4*nLevel;
  double length = 200;
  double height = 160;

  // 1st page (6 rows)
  int nrow_on_page1 = 6;
  int npad_p1 = nrow_on_page1*npadx;


  if(debug) cout << "creating cc1" << endl;
  // Note: need to create on heap, otherwise code crashes.
  TCanvas *cc1 = new TCanvas("cc1","cc1",(npadx*length), (nrow_on_page1*height)*1.05);

  
  TPad *c0 = new TPad("c0","title pad", 0, 0.9624, 1, 1     ); 
  TPad *c1 = new TPad("c1","main pad",  0, 0,      1, 0.9524);
  c0->Draw();
  c1->Draw();

  if(debug) cout << "drawing c0,c1" << endl;

  // banner
  c0->cd();
  TPaveLabel myinfo(0,0,1,1,"");
  myinfo.SetLabel( Form("Reliso:met correlation, %s, %s, %s, %s (7 TeV)", 
			nj.c_str(), sel.c_str(), metcut_label[metLabel].c_str(), etaRegion) );
  myinfo.SetFillColor(0);
  myinfo.Draw();

  //TCanvas c1("c1","c1",10,10,npadx*length,nrow_on_page1*height);
  c1->cd();
  c1->Divide(npadx,nrow_on_page1,0.004,0.004); //smaller separation between pads

  // 2nd page (5 rows)
  if(debug) cout << "creating c2" << endl;
  TCanvas *c2= new TCanvas("c2","c2",10,10,npadx*length,(npady-6)*height);
  c2->Divide(npadx,npady-6,0.004,0.004);

  if(debug) cout << "getting 2D histo" << endl;

  if(sel=="QCD") plot_2D_unweighted = false;

  //  cout << "amtb 0" << endl;

  unsigned int npro = 0; //n projected
  for(int i=0; i<nLevel; i++){
    if(debug)  cout << "i(nLevel)=" << i << endl;

    // eg: QCDest_isoVmet_NES_[uw]_[barrel]_1j__data

    string hname = var;
    if(plot_2D_unweighted) hname += "_uw";    
    if(etaRegion=="barrel"||etaRegion=="endcap") hname = hname +"_"+etaRegion;

    hname = hname + level[i] + "_" + nj + "__" + sel;

    if(debug) cout << hname << endl;

    h2[i] = (TH2F*)f.Get( Form("%s",hname) );

    //h[i]->SetTitle( Form("%s", title[i].c_str()) );
    h2[i]->GetXaxis()->SetTitle("RelIso");
    h2[i]->GetYaxis()->SetTitle("MET (GeV)");

    // project met
    h2[i]->ProjectionY("met");
    h[npro] = (TH1D*)met->Clone();
    //    cout << "amtb 1" << endl;

    h2[i]->ProjectionY("met",1,5);//low iso
    h[npro+1] = (TH1D*)met->Clone();

    const int iso_nbin = h2[i]->GetNbinsX();

    h2[i]->ProjectionY("met",6,iso_nbin);//high iso
    h[npro+2] = (TH1D*)met->Clone();

    npro += 3;
  }
  //  cout << "amtb 2" << endl;


  // Check projected histo
  //-----------------------
  if(debug) {
    cout << "check projected histo" << endl;
    for(int i=0; i<nh; i++){
      cout << "i=" << i << endl;
      h[i]->ls();
    }
  }
  

  // Adjust binning depending on number of entries
  //-----------------------------------------------
  // 2D
  if(debug) cout << "rebin 2D histo" << endl;
  for(int i=0; i<nLevel; i++){
    int rb2d = 1;
    if(h2[i]->GetEntries() < 20000) rb2d = 2;
    if(h2[i]->GetEntries() <  6000) rb2d = 3;
    if(h2[i]->GetEntries() <  3000) rb2d = 4;
    if(h2[i]->GetEntries() <  1000) rb2d = 5; //iso=0.1
    if(h2[i]->GetEntries() <   500) rb2d = 10;
    h2[i]->Rebin2D(rb2d,rb2d);
  }

  // 1D
  if(debug) cout << "rebin 1D histo" << endl;
  for(int i=0; i<nh; i++){
    int rb = 1;
    if(h[i]->GetEntries() < 1000) rb = 2;
    h[i]->Rebin(rb);
  }
  

  //---------
  //  Draw
  //---------
  if(debug) cout << "drawing h" << endl;

  // 1D plot
  int i = 0;
  for(int n=1; n<=npad; n++){

    if(debug) cout << "n = " << n << "  " << h[i] << endl;

    c1->cd(n);
    if(n>npad_p1) c2->cd(n-npad_p1); //end page

    // scattter plot
    // skip pad 1,5,9,13,17 (for scatter plot)
    if( ( n-1)%4==0 ) // if (n-1) is multiple integer of 4
      continue;

    if(h[i]>0) {
      h[i]->GetXaxis()->SetRangeUser(0,met_max);
      h[i]->SetXTitle("MET (GeV)");
      //if(n%4==2) h[i]->SetTitle("RelIso<0.1"); //2nd column
      if(n%4==3) h[i]->SetTitle("RelIso<0.1"); //3rd column
      if(n%4==0) h[i]->SetTitle("RelIso>0.1"); //4th column
      h[i]->Draw();
    }
    i++;    
  }
  

  if(debug) cout << "drawing h2" << endl;

  // scatter plot (on pad 1,5,9,14,17)
  gStyle->SetStatX(0.91);
  i=0;
  for(int n=1; n<=npad; n=n+4){
    //cout << "n = " << n << endl;
    c1->cd(n);
    if(n>npad_p1) c2->cd(n-npad_p1); //end page

    h2[i]->GetYaxis()->SetRangeUser(0,met_max); //met
    h2[i]->GetXaxis()->SetRangeUser(0,iso_max); //iso
    
    if(h2[i]==0) cout << h2[i]->GetName() << endl;
    h2[i]->Draw("colz");
    i++;
  }

  if(debug) cout << "save as: " << opname << endl;

  // single-page plot (needed for elog)
  //c1->SaveAs(Form("%s_1.png",opname.c_str()));
  //c2->SaveAs(Form("%s_2.png",opname.c_str()));

  // multiple-page pdf
  if(debug) cout << "saving" << endl;
  cc1->Print(Form("%s.ps(",opname.c_str()));
  c2->Print(Form("%s.ps)",opname.c_str()));
    
  if(debug) cout << "closing file" << endl;
  f.Close();

  if(debug) cout << "closing canvases" << endl;
  cc1->Close();
  c2->Close();
 
  if(debug) cout << "ps2pdf" << endl;
  gROOT->ProcessLine( Form(".!ps2pdf %s.ps",opname.c_str()) );
  gROOT->ProcessLine( Form(".!rm -f %s.ps",opname.c_str()) );
  // move pdf into subdir
  string dir="both_barrel_and_endcap";
  if(etaRegion=="barrel") dir="barrel_only";
  if(etaRegion=="endcap") dir="endcap_only";
  gROOT->ProcessLine( Form(".!mv -f %s.pdf %s",opname,dir) );

  // open the pdf
  // gROOT->ProcessLine(Form(".!open %s.pdf",opname));

  if(debug) cout << "end of main code" << endl;

}
