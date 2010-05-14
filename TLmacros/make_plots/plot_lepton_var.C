//
// 25 Mar 2010
// plot electron variables for electrons passing ET and eta cuts: multiple page pdf (ps2pdf)
//
//  enri1,2,3 bce1,2,3
//  QCD  tt  w  z  data data
//  compare x x x x x x
//
//---------------------------------------
const int nproc = 11;
const bool debug = 0;


void plot_lepton_var( ) {
  plot("ele", "barrel");
  plot("ele", "endcap");
  plot("muon");
}

void plot( string lepton="ele", string etaRegion="barrel" ) {

  //----------------------------
  //  Specify variable to plot
  //----------------------------
  //  int var_index = 0;
  string out ;
  if(lepton=="ele") out = Form("ele_var_%s.ps",etaRegion);
  else              out = Form("muon_var.ps");

  
  //-------------
  //  Variables
  //-------------
  vector<string> var;
  if(lepton=="ele"){
    var.push_back("electron/ele_hadOverEm"); //plus '_barrel' or '_endcap'
    var.push_back("electron/ele_dEtaIn");
    var.push_back("electron/ele_dPhiIn");
    var.push_back("electron/ele_sigmaIEtaIEta");
    var.push_back("electron/ele_EoverPIn");
    //  var.push_back("electron/ele_fBrem");
    var.push_back("electron/ele_tIso");
    var.push_back("electron/ele_cIso");
    var.push_back("electron/ele_tIso_dr03"); //tk
    var.push_back("electron/ele_cIso_dr03"); //calo
    var.push_back("electron/ele_eIso_dr03"); //ecal
    var.push_back("electron/ele_hIso_dr03"); //hcal
    var.push_back("electron/ele_tIso_dr04"); //tk
    var.push_back("electron/ele_cIso_dr04"); //calo
    var.push_back("electron/ele_eIso_dr04"); //ecal
    var.push_back("electron/ele_hIso_dr04"); //hcal
  }else{
    var.push_back("muon/muon_normchi2");
    var.push_back("muon/muon_d0");
    var.push_back("muon/muon_tkHits");
  }

  vector<string> type; //linear or log in y
  if(lepton=="ele"){
    type.push_back("log"); //H/E
    type.push_back("lin"); //d eta
    type.push_back("lin"); //d phi
    type.push_back("lin"); //sigma ieta ieta
    type.push_back("lin"); //e/p
    //type.push_back("lin"); //fBrem
    type.push_back("log"); //tIso   Def
    type.push_back("lin"); //cIso
    type.push_back("log"); //tk    dr03
    type.push_back("lin"); //calo
    type.push_back("lin"); //ecal
    type.push_back("log"); //hcal
    type.push_back("log"); //tk   dr04
    type.push_back("lin"); //calo
    type.push_back("lin"); //ecal
    type.push_back("log"); //hcal
  }else{
    // muon
    type.push_back("lin"); //norm chi2
    type.push_back("log"); //d0
    type.push_back("lin"); //tkHits
  }

  //-------------
  // Rebin
  //-------------
  vector<int> rbs;
  /*
  rbs.push_back(1); //
  rbs.push_back(1); //
  rbs.push_back(1); //
  rbs.push_back(1); //
  rbs.push_back(1); //
  rbs.push_back(1); //
  rbs.push_back(1); //
  rbs.push_back(1); //
  rbs.push_back(1); //
  */

  //  for(short p=0; p<rbs.size(); ++p ) cout << "rb " << p <<" = "<< rbs[p] << endl;





  TFile *f = new TFile("../test_mc_mixture.root");
  

  gStyle->SetOptStat(1110);//off title
  gStyle->SetStatH(0.2); //bigger
  gStyle->SetStatW(0.4);
  gStyle->SetStatY(0.88);
  gStyle->SetPadLeftMargin(0.15); //bigger
  gStyle->SetPadBottomMargin(0.15); //bigger 
  gStyle->SetPadTopMargin(0.14); //bigger 
  gStyle->SetTitleH(0.08);
  gStyle->SetTitleAlign(13); //23 = left-aligned
  gStyle->SetTitleFillColor(0); //23 = left-aligned
  gROOT->UseCurrentStyle();

  //--------------
  // Color (fill)
  //--------------
  Color_t Col_sig = kRed-7;
  Color_t Col_qcd = kAzure-9;
  Color_t Col_wj  = kGreen-7;
  Color_t Col_zj = kYellow-9;
  Color_t Col_vqq = kGray+1;
  Color_t Col_stop = kViolet-4;
  // new
  Color_t Col_data = kGreen-2;
  Color_t Col_enri = kBlue-9;
  Color_t Col_bce  = kRed-9;


  //-------------
  // Canvas
  //-------------
  int nx = 6;
  int ny = 3;
  int npad_on_canvas = nx*ny;

  TCanvas *cc1 = new TCanvas("cc1","iso",(nx*200)*1,(ny*200)*1.05);
  TPad *c0 = new TPad("c0","title pad", 0, 0.96, 1, 1 ); 
  TPad *c1 = new TPad("c1","main pad", 0, 0, 1, 0.95 );
  c0->Draw();
  c1->Draw();

  c0->cd();

  TPaveLabel *myinfo = new TPaveLabel(0,0,1,1,"TEST");
  myinfo->Draw();

  c1->cd();
  c1->Divide(nx,ny,0.002,0.003); //smaller margin

  gPad->GetCanvas()->FeedbackMode(kTRUE); //TEST - magic line to overcome gPad error as follow
  /// Error: non class,struct,union object $gPad used with . or -> plot_lepton_var.C:221:
  

  string proc[nproc] ={"enri1","enri2","enri3","bce1","bce2","bce3",
		       "QCD","ttbar","wj","zj","data"};
  

  int nvar = var.size();
  
  int npad=1;
  int iv = 0; 

  TH1D *h[50];
  TH1D *hcheck;

  for(int n=0; n<nvar; ++n){
    
    cout << "var: " << var[n] << endl;

    // First check if this variable is in the root file
    if(lepton=="ele") var[n] = var[n] + "_" + etaRegion;

    hcheck = (TH1D*)f->Get(Form("%s__%s", var[n].c_str(), proc[0].c_str()));
    if(hcheck==0) { cout << var[n] << " histo not found!" << endl; continue; }
    

    string banner_info = Form("Distributions of %s for selected leptons (et/pt + #eta cut)",var[n].c_str());
    myinfo->SetLabel(banner_info.c_str());
    c0->cd();
    myinfo->Draw();
  
    //    int rb = rbs[n];
    int rb = 1;
    if(debug) cout << "rb=" << rb << endl;
    
    
    for(int i=0; i<nproc; ++i){
      
      if(debug) cout << "proc " << i << " "<< proc[i] << endl;

      h[i] = (TH1D*)f->Get(Form("%s__%s", var[n].c_str(), proc[i].c_str()));
    
      if(rb>1) {
	cout << "nbin b4: " << h[i]->GetNbinsX() << endl;
	h[i]->Rebin(rb);
	cout << "nbin after: "<< h[i]->GetNbinsX() << endl;
      }
      // Set colour
      if( h[i] > 0 ) {
	if( proc[i].find("enri") !=string::npos ) h[i]->SetFillColor(Col_enri);
	if( proc[i].find("bce")  !=string::npos ) h[i]->SetFillColor(Col_bce);
	if( proc[i].find("ttbar")!=string::npos ) h[i]->SetFillColor(Col_sig);
	if( proc[i].find("wj")   !=string::npos ) h[i]->SetFillColor(Col_wj);
	if( proc[i].find("zj")   !=string::npos ) h[i]->SetFillColor(Col_zj);
	if( proc[i].find("data") !=string::npos ) h[i]->SetFillColor(Col_data);
	if( proc[i].find("QCD")  !=string::npos ) h[i]->SetFillColor(Col_qcd);
	c1->cd(i+1);

	if(debug) cout << type[n] << ", nEntries: " << h[i]->GetEntries() << endl;
	
	h[i]->Draw("hist");
	
	if( h[i]->GetEntries() == 0 ) continue;

	//cout << "amituofo" << endl;
	
	if(type[n]=="log") {
	  //cout << "set log" << endl;
	  gPad->SetLogy(1); 
	}
	else {
	  //cout << "dont set log" << endl;
	  gPad->SetLogy(0);	
	}

      }//if histo pointer valid

    }//loop of proc


    // pad no 12: data and break down
    if(debug) cout << "go to pad 12" << endl;
    c1->cd(12);
    if(type[n]=="log") gPad->SetLogy(1);  else  gPad->SetLogy(0);
    TH1D *qcd1;
    TH1D *sig1;
    TH1D *wj1;
    TH1D *zj1;
    int q = 6; //no. 7 is QCD

    qcd1 = (TH1D*)h[6]->Clone();
    sig1 = (TH1D*)h[7]->Clone();
    wj1  = (TH1D*)h[8]->Clone();
    zj1  = (TH1D*)h[9]->Clone();
    all  = (TH1D*)h[10]->Clone();

    all->Draw("ahist");
    sig1->Draw("ahist same");
    qcd1->Draw("ahist same");
    wj1->Draw("ahist same");
    zj1->Draw("ahist same");
    qcd1->Draw("axis same");


    //cout << "amtb 2" << endl;

    // pad 13-16: compare 
    double max_sig;
    double max_qcd;
    double max_wj;
    double max_zj;

    TH1D *qcd;
    TH1D *sig;
    TH1D *wj;
    TH1D *zj;
    qcd = (TH1D*)qcd1->Clone();
    sig = (TH1D*)sig1->Clone();
    wj  = (TH1D*)wj1->Clone();
    zj  = (TH1D*)zj1->Clone();

    wj->SetFillStyle(3425);
    zj->SetFillStyle(3452);
    qcd->SetFillStyle(3475);
    sig->SetFillStyle(3457);

    max_sig = sig->GetMaximum();
    max_qcd = qcd->GetMaximum();
    max_wj  = wj->GetMaximum();
    max_zj  = zj->GetMaximum();

    TH1D *h1;
    TH1D *h2;
    //cout << "amtb 3" << endl;

    int pad_no = 13;
  
    // 3.1) compare QCD with wj
    if(debug) cout << "go to pad 13" << endl;

    c1->cd(pad_no);
    if(type[n]=="log") gPad->SetLogy(1);  else  gPad->SetLogy(0);
    if( max_qcd > max_wj ) {  h1 = qcd; h2 = wj; }
    else {                    h2 = qcd; h1 = wj; }  
    h1->SetTitle(Form("QCD vs W"));
    h1->DrawClone("ahist");
    h2->Draw("ahist same");
    h2->Draw("axis same");
    ++pad_no;
  
    // 3.2) compare QCD with zj
    c1->cd(pad_no);   
    if(type[n]=="log") gPad->SetLogy(1);  else  gPad->SetLogy(0);
    if( max_qcd > max_zj ) {  h1 = qcd; h2 = zj; }
    else {                    h2 = qcd; h1 = zj; } 
    h1->SetTitle(Form("QCD vs Z"));
    h1->DrawClone("ahist");
    h2->Draw("ahist same");
    h2->Draw("axis same");
    ++pad_no;
  
    // 3.3) compare QCD with sig
    c1->cd(pad_no);   
    if(type[n]=="log") gPad->SetLogy(1);  else  gPad->SetLogy(0);
    if( max_qcd > max_sig ) {  h1 = qcd; h2 = sig; }
    else {                     h2 = qcd; h1 = sig; }
    h1->SetTitle(Form("QCD vs Signal"));
    h1->DrawClone("ahist");
    h2->Draw("ahist same");
    h2->Draw("axis same");
    ++pad_no;

    // 3.4) compare sig with wj
    c1->cd(pad_no);   
    if(type[n]=="log") gPad->SetLogy(1);  else  gPad->SetLogy(0);
    if( max_sig > max_wj ) {  h1 = sig; h2 = wj; }
    else {                    h2 = sig; h1 = wj; }
    h1->SetTitle(Form("Signal vs W"));
    h1->DrawClone("ahist");
    h2->Draw("ahist same");
    h2->Draw("axis same");
    ++pad_no;

    // 3.5) compare sig with zj
    c1->cd(pad_no);
    if(type[n]=="log") gPad->SetLogy(1);  else  gPad->SetLogy(0);
    if( max_sig > max_wj ) {  h1 = sig; h2 = wj; }
    else {                    h2 = sig; h1 = wj; }  
    h1->SetTitle(Form("Signal vs W"));
    h1->DrawClone("ahist");
    h2->Draw("ahist same");
    h2->Draw("axis same");
    ++pad_no;

    // 3.6) compare wj with zj
    c1->cd(pad_no);
    if(type[n]=="log") gPad->SetLogy(1);  else  gPad->SetLogy(0);
    if( max_wj > max_zj ) {  h1 = wj; h2 = zj; }
    else {                   h2 = wj; h1 = zj; }  
    h1->SetTitle(Form("W vs Z"));
    h1->DrawClone("ahist");
    h2->Draw("ahist same");
    h2->Draw("axis same");


    // Print picture to file
    if(n==0)    
      cc1->Print(Form("%s(",out.c_str())); //open
    else if(n<(nvar-1)) 
      cc1->Print(Form("%s", out.c_str()));

  }//var

  // close ps
  cc1->Print(Form("%s)",out.c_str()));
  cc1->Close();
  gROOT->ProcessLine(Form(".!ps2pdf %s",out.c_str()));
  gROOT->ProcessLine(Form(".!rm %s",out.c_str()));


}
