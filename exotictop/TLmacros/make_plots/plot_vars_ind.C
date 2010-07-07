//
// 5 Mar 2010
// plot met: multiple page pdf (ps2pdf)
//
// vars: iso, met, mtw, DPhiEmet, DPhiMetJet.
//
//         1j  2j  3j  >=4j   |    (2nd page)  |   (3rd page)
//  enri1  x   x   x    x     |  qcd  x x x x  | QCD v wj
//  enri2                     |   tt           | QCD v zj
//  enri3                     |   wj           | QCD v sig
//  bce1                      |   zj           | sig v wj
//  bce2                      | data           | sig v zj
//  bce3                      | data(+mc)      |  wj v zj
//-----------------------------------------------------------------

// Globals
const int nproc = 12;
string mydir;
string mylabel;
double mylabelsize;
string myxt;
double bw[4];
bool setxrange = false;
double xmax = 0;

const bool debug = false;


//---------------------------------------------------------------
bool plot_vars_ind() {

  plot_iso_ind();
  plot_met_ind();
  plot_mtw_ind();
  plot_DPhiEmet_ind();
  plot_DPhiMetJet_ind();

  // study met
  plot_met_gen();
  plot_met_gen_diff();
  plot_met_gen_dphi();
}


//---------------------------------------------------------------
void plot_iso_ind() {

  mydir = "QCD_estimation";
  myxt  = "RelIso";
  mylabelsize = 0.06;
  //  bw[0]=0.1; bw[1]=0.1; bw[2]=0.1; bw[3]=0.1;  //bin width
  double mbw[4] = { 0.1, 0.1, 0.1, 0.1 };
  memcpy(bw, mbw, sizeof(bw));
  setxrange = true;
  xmax = 1.4;

  // NES
  mylabel = "RelIso in Signal Region (7 TeV)";
  do_plot_vars_ind("QCDest_CombRelIso");
  
  // AES (old)
  mydir = "QCD_estimation/AES";
  mylabel = "RelIso in Old Control Region (7 TeV)";
  do_plot_vars_ind("QCDest_CombRelIso_AES");

  // AES (new: planB3_e20)
  mydir = "QCD_planB";
  mylabel = "RelIso in New Control Region: B3_e20 (7 TeV)";
  do_plot_vars_ind("QCDest_CombRelIso_AES_planB3_e20");

}
//---------------------------------------------------------------

void plot_met_ind() {
  mydir = "MET";
  myxt  = "MET (GeV)";
  mylabelsize = 0.05;
  //  bw[0]=2; bw[1]=2; bw[2]=4; bw[3]=4;  //GeV (bin width)
  double mbw[4] = { 2, 2, 4, 4 };
  memcpy(bw, mbw, sizeof(bw));
  setxrange = true;
  xmax = 150;

  mylabel = "MET after all cuts except MET: #mu-MET (7 TeV)";
  do_plot_vars_ind("met_mu");

  mylabel = "MET after all cuts except MET: Type1-MET (7 TeV)";
  do_plot_vars_ind("met_t1");
}
//---------------------------------------------------------------

void plot_mtw_ind() {
  mydir = "mtw";
  myxt  = "M_{T}(W) (GeV)";
  mylabelsize = 0.05;
  double mbw[4] = { 2, 4, 4, 4 }; //orig: 2GeV/bin, 100 bins 0-200
  memcpy(bw, mbw, sizeof(bw));
  setxrange = false;

  mylabel = "M_{T}(W) after all cuts except MET: #mu-MET (7 TeV)";
  do_plot_vars_ind("mtw_mu");

  mylabel = "M_{T}(W) after all cuts except MET: Type1-MET (7 TeV)";
  do_plot_vars_ind("mtw_t1");
}
//---------------------------------------------------------------

void plot_DPhiEmet_ind() {
  mydir = "DPhi_ele_met";
  myxt  = "#Delta#Phi(e,met)";
  mylabelsize = 0.05;
  double mbw[4] = { 0.1, 0.1, 0.2, 0.2 };
  memcpy(bw, mbw, sizeof(bw));
  setxrange = false;

  mylabel = "#Delta#Phi(e_{iso},#mu-met) after all cuts except #mu-MET (7 TeV)";
  do_plot_vars_ind("DPhiEmet_mu");

  mylabel = "#Delta#Phi(e_{iso},t1-met) after all cuts except t1-MET (7 TeV)";
  do_plot_vars_ind("DPhiEmet_t1");
}
//---------------------------------------------------------------

void plot_DPhiMetJet_ind() {
  mydir = "DPhi_met_jet";
  myxt  = "#Delta#Phi(met,jet)";
  mylabelsize = 0.05;
  double mbw[4] = { 0.2, 0.2, 0.2, 0.2 };
  memcpy(bw, mbw, sizeof(bw));
  setxrange = false;

  mylabel = "#Delta#Phi(#mu-met,nearest jet) after all cuts except #muMET (7 TeV)";
  do_plot_vars_ind("DPhiMetJet_mu");

  mylabel = "#Delta#Phi(Type1-met,nearest jet) after all cuts except t1MET (7 TeV)";
  do_plot_vars_ind("DPhiMetJet_t1");
}
//---------------------------------------------------------------

void plot_met_gen() {
  mydir = "MET/gen";
  myxt  = "GenMET (GeV)";
  mylabelsize = 0.05;
  double mbw[4] = { 5, 5, 5, 5 };
  memcpy(bw, mbw, sizeof(bw));
  setxrange = true;
  xmax = 100;

  mylabel = "GenMET after all cuts except MET (7 TeV) barrel";
  do_plot_vars_ind("met_gen_BA");

}
//---------------------------------------------------------------

void plot_met_gen_diff() {
  mydir = "MET/gen";
  myxt  = "#Delta MET (GeV)";
  mylabelsize = 0.05;
  double mbw[4] = { 4, 4, 4, 4 }; //ori=2 GeV-bin
  memcpy(bw, mbw, sizeof(bw));
  setxrange = false;
  //xmax = 100;

  mylabel = "#Delta MET=Type1-Gen, after all cuts except MET (7 TeV) barrel";
  do_plot_vars_ind("met_gen_diff_t1_BA");

  mylabel = "#Delta MET=#muMET-Gen, after all cuts except MET (7 TeV) barrel";
  do_plot_vars_ind("met_gen_diff_mu_BA");

}
//---------------------------------------------------------------

void plot_met_gen_dphi() {
  mydir = "MET/gen";
  myxt  = "#Delta#Phi";
  mylabelsize = 0.05;
  double mbw[4] = { .2, .2, .2, .2 }; //ori=0.2 GeV-bin
  memcpy(bw, mbw, sizeof(bw));
  setxrange = false;
  //xmax = 100;

  mylabel = "#Delta#Phi=Type1-Gen, after all cuts except MET (7 TeV) barrel";
  do_plot_vars_ind("met_gen_dphi_t1_BA");

  mylabel = "#Delta#Phi=#muMET-Gen, after all cuts except MET (7 TeV) barrel";
  do_plot_vars_ind("met_gen_dphi_mu_BA");

}
//---------------------------------------------------------------


//--------------------------
//    m a i n   c o d e
//--------------------------

bool do_plot_vars_ind( string var = "met_mu" ) {

  //  double bw[4] = {2,2,4,4};  // GeV (bin width)
  //  int rb[4];  
  //  for(short a=0; a<4; ++a)  rb[a]= bw[a]/1; //originally bw=1 GeV (200 bins for 0-200)


  //string out = Form("met_mu_ind_w%.2f.ps",bw[0]);
  //if(plot_t1met) out = Form("met_t1_ind_w%.2f.ps",bw[0]);
  //  string mettype = "mu";
  //  if(plot_t1met) mettype = "t1";
  string out = var + "_ind_" + Form("w%.2f",bw[0]) + ".ps";
  cout << "var:     " << var << endl;
  cout << "output:  " << out << endl;


  TFile *f = new TFile("../test_mc_mixture.root");
  

  gROOT->SetStyle("Default"); //Back to default style
  gStyle->SetOptStat(1110);//off title
  gStyle->SetStatH(0.2); //bigger
  gStyle->SetStatW(0.4);
  gStyle->SetStatY(0.9);
  gStyle->SetPadLeftMargin(0.15); //bigger
  gStyle->SetPadBottomMargin(0.15); //bigger 
  gStyle->SetPadTopMargin(0.16); //bigger 
  gStyle->SetTitleH(0.08);
  gStyle->SetTitleAlign(13); //23 = left-aligned
  gROOT->UseCurrentStyle();
  TGaxis::SetMaxDigits(4);

  //--------------
  // Color (fill)
  //--------------
  Color_t Col_sig  = kRed-7;
  Color_t Col_qcd  = kAzure-9;
  Color_t Col_wj   = kGreen-7;
  Color_t Col_zj   = kYellow-9;
  Color_t Col_vqq  = kGray+1;
  Color_t Col_stop = kViolet-4;
  // new
  Color_t Col_data = kGreen-2;
  Color_t Col_enri = kBlue-9;
  Color_t Col_bce  = kRed-9;



  
  int nx = 4;
  int ny = 6;
  int npad_on_canvas = nx*ny;

  TCanvas *cc1 = new TCanvas("cc1","iso",(nx*200)*1.1,(ny*200)*1.05);
  TPad *c0 = new TPad("c0","title pad", 0.1, 0.9624, 1, 1 ); 
  TPad *c1 = new TPad("c1","main pad", 0.1, 0, 1, 0.9524 );
  c0->Draw();
  c1->Draw();

  c0->cd();

  TPaveLabel myinfo(0,0,1,1,"");
  //  if(!plot_t1met) myinfo.SetLabel("MET after all cuts except MET: #mu-MET (7 TeV)");    
  //  else            myinfo.SetLabel("MET after all cuts except MET: Type1-MET (7 TeV)");
  myinfo.SetLabel( mylabel.c_str() );
  myinfo.Draw();

  c1->cd();
  c1->Divide(nx,ny,0.003,0.002); //smaller margin

  

  string proc[nproc] ={"enri1","enri2","enri3","bce1","bce2","bce3",
		       "QCD","ttbar","wj","zj","data","data"};
  string nj[5] ={"0j","1j","2j","3j","4mj"};
  char *jlabel[4] = {"(1j)","(2j)","(3j)","(#geq4j)"};
  int nvar = nproc*4;

  int npad=1;
  int iv = 0; 

  TH1D *h[50];
  //  TH1D *hh[50];
  //  int nj_start = 1; //for NES, start from 1j
  //  if(plot_AES) nj_start = 0; //for AES, start from 0j



  for(int i=0; i<nproc; ++i){


    // when drawing second page, clear canvas
    if (npad==npad_on_canvas+1) {
      c1->Clear();
      c1->Divide(nx,ny,0.003,0.003);
      npad=1;
    }

    cout << "i="<< i << "  proc = "<< proc[i]<< endl;


    for(int j=1; j <= 4; ++j){ //nj
      //    for(int j=nj_start; j < (4+nj_start); ++j){ //nj
      c1->cd(npad);
      //  cout << " nj = "<< nj[j]<< endl;
      /*
      if(!plot_t1met) 
	h[iv] = (TH1D*)f->Get( Form( "%s/met_mu_%s__%s", nj[j].c_str(), proc[i].c_str() ) );
      else
	h[iv] = (TH1D*)f->Get( Form( "%s/met_t1_%s__%s", nj[j].c_str(), proc[i].c_str() ) );
      */
      string myhisto = mydir + "/" +  var + "_" + nj[j] + "__" + proc[i];
      h[iv] = (TH1D*)f->Get( myhisto.c_str() );
      
      //cout << "histo pointer: "<< h[iv] << endl;
      if(h[iv]==0) { 
	cout << "histo not found, exit." << endl; 
	return 0;
      }

      if( h[iv] > 0 ) {
	
	if( proc[i].find("enri") !=string::npos ) h[iv]->SetFillColor(Col_enri);
	if( proc[i].find("bce")  !=string::npos ) h[iv]->SetFillColor(Col_bce);
	if( proc[i].find("ttbar")!=string::npos ) h[iv]->SetFillColor(Col_sig);
	if( proc[i].find("wj")   !=string::npos ) h[iv]->SetFillColor(Col_wj);
	if( proc[i].find("zj")   !=string::npos ) h[iv]->SetFillColor(Col_zj);
	if( proc[i].find("data") !=string::npos ) h[iv]->SetFillColor(Col_data);
	if( proc[i].find("QCD")  !=string::npos ) h[iv]->SetFillColor(Col_qcd);

	h[iv]->SetTitle(Form("%s %s",proc[i].c_str(),jlabel[j-1]));

	//	h[iv]->SetLabelSize(0.05,"xy");//0.06
	h[iv]->SetLabelSize(mylabelsize,"xy");//0.06
	h[iv]->SetXTitle(myxt.c_str());
	h[iv]->GetXaxis()->SetTitleSize(0.06);

	//cout << "h "<< iv <<":"<< h[iv]->GetNbinsX() << endl;
	double xrange = h[iv]->GetXaxis()->GetXmax() - h[iv]->GetXaxis()->GetXmin();
	double bw_orig = xrange/h[iv]->GetNbinsX();
	int rb = int(bw[j-1]/bw_orig + 0.5); //cast to int
	if(debug){
	  cout << "xrange           " << xrange << endl;
	  cout << "bw_orig          " << bw_orig << endl;
	  cout << "bw               " << bw[j-1] << endl;
	  cout << "rb (bw/bw_orig)  " << rb << endl;
	}
	if(rb==0) return;
	if(h[iv]->GetNbinsX() > 60) h[iv]->Rebin(rb);
	//cout << " after: "<< h[iv]->GetNbinsX() << endl;

	h[iv]->SetMinimum(0); //min on y-axis
	h[iv]->Draw("hist");
	if(setxrange) h[iv]->GetXaxis()->SetRangeUser(0,xmax);

	//h[iv]->Draw("");
	//	hh[iv] = (TH1D*)h[iv]->Clone();
      }
    
      iv++;
      npad++;
    }//nj


    if(iv==npad_on_canvas)  cc1->Print(Form("%s(",out.c_str())); //open ps
    //    if(iv==npad_on_canvas*2)  c1->Print("h1.ps)");
    //    if(iv==nvar)            c1->Print("h1.ps)");  
    //        if(iv==nvar) c1->Print("h1.ps)");
    
  }//proc
  
  
  //--------------------------
  // Super-impose on 1 plot
  //--------------------------
  TH1D *qcd[4]; //1-4j
  TH1D *sig[4];
  TH1D *wj[4];
  TH1D *zj[4];

  int pad_no = 5*4+1;

  int q = 24;

  for(int k=0; k<4; ++k){
    if(debug) cout << "nj " << k+1 << endl;
    qcd[k] = (TH1D*)h[q+k]->Clone();
    sig[k] = (TH1D*)h[q+k+4]->Clone();
    wj[k]  = (TH1D*)h[q+k+8]->Clone();
    zj[k]  = (TH1D*)h[q+k+12]->Clone();

    c1->cd(pad_no);
    qcd[k]->Draw("ahist same");
    if(k+1 < 3){
      wj[k]->Draw("ahist same");
      zj[k]->Draw("ahist same");
      sig[k]->Draw("ahist same");
      sig[k]->Draw("axis same");
    }else{ //3,4j
      sig[k]->Draw("ahist same");
      wj[k]->Draw("ahist same");
      zj[k]->Draw("ahist same");
      sig[k]->Draw("axis same");     
    }

    pad_no++;
  }
  cc1->Print(Form("%s",out.c_str()));

  // page 3: compare diff processes
  c1->Clear();
  c1->Divide(nx,ny,0.003,0.003);

  pad_no = 1;

  double max_sig[4];
  double max_wj[4];
  double max_zj[4];
  double max_qcd[4];


  for(int k=0; k<4; ++k){
    wj[k]->SetFillStyle(3425);
    zj[k]->SetFillStyle(3452);
    qcd[k]->SetFillStyle(3475);
    sig[k]->SetFillStyle(3457);

    max_sig[k] = sig[k]->GetMaximum();
    max_wj[k]  = wj[k]->GetMaximum();
    max_zj[k]  = zj[k]->GetMaximum();
    max_qcd[k] = qcd[k]->GetMaximum();
  }

  TH1D *h1;
  TH1D *h2;

  // 3.1) compare QCD with wj
  for(int k=0; k<4; ++k){
    c1->cd(pad_no);
    
    if( max_qcd[k] > max_wj[k] ) {  h1 = qcd[k]; h2 = wj[k]; }
    else {                          h2 = qcd[k]; h1 = wj[k]; }

    h1->SetTitle(Form("QCD vs W+j %s",jlabel[k]));
    h1->DrawClone("ahist");
    h2->Draw("ahist same");
    h2->Draw("axis same");
    ++pad_no;
  }

  // 3.2) compare QCD with zj
  for(int k=0; k<4; ++k){
    c1->cd(pad_no);
    if( max_qcd[k] > max_zj[k] ) {  h1 = qcd[k]; h2 = zj[k]; }
    else {                          h2 = qcd[k]; h1 = zj[k]; }
    h1->SetTitle(Form("QCD vs Z+j %s",jlabel[k]));
    h1->DrawClone("ahist");
    h2->Draw("ahist same");
    h2->Draw("axis same");
    ++pad_no;
  }

  // 3.3) compare QCD with sig
  for(int k=0; k<4; ++k){
    c1->cd(pad_no);
    if( max_qcd[k] > max_sig[k] ) {  h1 = qcd[k]; h2 = sig[k]; }
    else {                           h2 = qcd[k]; h1 = sig[k]; }
    h1->SetTitle(Form("QCD vs Signal %s",jlabel[k]));
    h1->DrawClone("ahist");
    h2->Draw("ahist same");
    h2->Draw("axis same");
    ++pad_no;
  }

  // 3.4) compare sig with wj
  for(int k=0; k<4; ++k){
    c1->cd(pad_no);
    if( max_sig[k] > max_wj[k] ) {  h1 = sig[k]; h2 = wj[k]; }
    else {                          h2 = sig[k]; h1 = wj[k]; }
    h1->SetTitle(Form("Signal vs W+j %s",jlabel[k]));
    h1->DrawClone("ahist");
    h2->Draw("ahist same");
    h2->Draw("axis same");
    ++pad_no;
  }

  // 3.5) compare sig with zj
  for(int k=0; k<4; ++k){
    c1->cd(pad_no);
    if( max_sig[k] > max_zj[k] ) {  h1 = sig[k]; h2 = zj[k]; }
    else {                          h2 = sig[k]; h1 = zj[k]; }
    h1->SetTitle(Form("Signal vs Z+j %s",jlabel[k]));
    h1->DrawClone("ahist");
    h2->Draw("ahist same");
    h2->Draw("axis same");
    ++pad_no;
  }

  // 3.6) compare wj with zj
  for(int k=0; k<4; ++k){
    c1->cd(pad_no);
    if( max_wj[k] > max_zj[k] ) {  h1 = wj[k]; h2 = zj[k]; }
    else {                         h2 = wj[k]; h1 = zj[k]; }
    h1->SetTitle(Form("W+j vs Z+j %s",jlabel[k]));
    //    h1->SetTitleSize(0.2);
    h1->DrawClone("ahist");
    h2->Draw("ahist same");
    h2->Draw("axis same");
    ++pad_no;
  }

  // close ps
  cc1->Print(Form("%s)",out.c_str()));
  cc1->Close();
  gROOT->ProcessLine(Form(".!ps2pdf %s",out.c_str()));
  gROOT->ProcessLine(Form(".!rm %s",out.c_str()));
  
  return 0;
}
