//
// 4 Mar 2010
// plot reliso, 3 pads:
// old AES and new AES
//  1) control sample (AES), 1j, composition
//  2) control sample vs QCD in signal sample (AES vs NES) in 1j
//  3) QCD in signal sample: compare 1j,2j,3j,4mj
//--------------------------------------------------

string myAES;
string myop;
bool debug = false;

void plot_iso_aes_compo_nes_3in1(){
  plot_old_AES();
  plot_new_AES();
}

void plot_old_AES(){
  myAES = "QCD_estimation/AES/QCDest_CombRelIso_AES";
  myop  = "old_AES_e30";
  do_plot_iso_aes_compo_nes_3in1(); //old
}

void plot_new_AES(){
  myAES = "QCD_planB/QCDest_CombRelIso_AES_planB3_e20";
  myop  = "new_AES_planB3_e20";
  do_plot_iso_aes_compo_nes_3in1(); //old
}


void do_plot_iso_aes_compo_nes_3in1(){

  cout << "\n QCD Control Sample is " << myAES << endl;

  double bw = 0.1;
  int rb = bw/0.01;

  double up = 1.49;


  //---------------
  //  input file
  //---------------
  char *file_data = "../test_mc_mixture.root";
  char *file_mc   = "../test_mc_mixture.root";

  TFile *fdata = new TFile(file_data);
  TFile *fmc   = new TFile(file_mc);
  

  //---------------
  //  output file
  //---------------
  string out = Form("iso_aes_compo_nes_3in1_%s_w%.2f",myop.c_str(), bw);
  //  string out = Form("iso_aes_compo_nes_3in1_e%d_w%.2f",take_AES_et, bw);
  

  //--------------
  // Color (fill)
  //--------------
  Color_t Col_sig = kRed-7;
  Color_t Col_qcd = kAzure-9;
  Color_t Col_wj  = kGreen-7;
  Color_t Col_zj  = kYellow-9;
  Color_t Col_vqq = kGray+1;
  Color_t Col_stop = kViolet-4;



  gStyle->SetTitleH(0.07);
  gStyle->SetTitleW(0.9);
  gStyle->SetTitleAlign(13); //Def in style =11 left-bottom
  gStyle->SetPadTopMargin(0.12);
  gStyle->SetOptStat(1);


  //--------------------------
  // PLOT 1: AES composition
  //--------------------------
  cout << "Amituofo: Pad 1"<< endl;

  TCanvas *c1 = new TCanvas("c1","iso (AES, 1j)", 3*400,400);
  c1->Divide(3,1);
  c1->cd(1);
  gPad->SetGrid();
  //  gPad->SetTicks();
  
  TH1D *h_sig = (TH1D*)fmc->Get( Form("%s_1j__ttbar", myAES.c_str()) );
  TH1D *h_wj  = (TH1D*)fmc->Get( Form("%s_1j__wj",    myAES.c_str()) );
  TH1D *h_zj  = (TH1D*)fmc->Get( Form("%s_1j__zj",    myAES.c_str()) );
  TH1D *h_qcd = (TH1D*)fmc->Get( Form("%s_1j__QCD",   myAES.c_str()) );

  h_sig->SetFillColor(Col_sig);
  h_wj ->SetFillColor(Col_wj);
  h_zj ->SetFillColor(Col_zj);
  h_qcd->SetFillColor(Col_qcd);

  h_sig->Rebin(rb);
  h_wj ->Rebin(rb);
  h_zj ->Rebin(rb);
  h_qcd->Rebin(rb);
    
  h_sig->GetXaxis()->SetRangeUser(0,up);
  h_wj ->GetXaxis()->SetRangeUser(0,up);
  h_zj ->GetXaxis()->SetRangeUser(0,up);
  h_qcd->GetXaxis()->SetRangeUser(0,up);
  

  THStack *h = new THStack("stack","RelIso (AES,1j)");
  h->Add(h_qcd);
  h->Add(h_sig);
  h->Add(h_wj);
  h->Add(h_zj);


  // Do this because of usual strange ROOT behaviours!
  // To set the x range.
  h->SetTitle("(a) Composition of Control Sample (1j)");

  h->Draw("hist"); //need to draw first to get xaxis 

  h->GetXaxis()->SetRangeUser(0,up);
  h->GetXaxis()->SetTitle("RelIso");
  h->Draw("hist");
  gPad->Update();
  
  //cout << "yaxis_max: "<< c1->GetUymax() << endl;
  //cout << "xaxis_max: "<< c1->GetUxmax() << endl;

  // y axis (no label)
  float yaxis_max = gPad->GetUymax(); 
  TGaxis *yaxis = new TGaxis(gPad->GetUxmin(),gPad->GetUymin(),
                             gPad->GetUxmin(),gPad->GetUymax(),
                             0,yaxis_max,510,"-RU");
  yaxis->SetLabelFont(42);
  yaxis->Draw();
  gPad->Update();
  // x axis (no label)
  float xaxis_max = gPad->GetUxmax(); 
  TGaxis *xaxis = new TGaxis(gPad->GetUxmin(),gPad->GetUymin(),
                             gPad->GetUxmax(),gPad->GetUymin(),
                             0,xaxis_max,510,"+RU");
  xaxis->SetLabelFont(42);
  xaxis->Draw();

  TLegend leg(0.6,0.6,0.85,0.85);
  leg.AddEntry(h_qcd, "QCD",    "f");
  leg.AddEntry(h_wj,  "W+jets", "f");
  leg.AddEntry(h_zj,  "Z+jets", "f");
  leg.AddEntry(h_sig, "ttbar",  "f");
  leg.SetTextFont(42);
  leg.SetFillColor(0);
  leg.Draw();





  //----------------------
  // PLOT 2: AES vs NES
  //----------------------
  cout << "Amituofo: Pad 2"<< endl;

  bool draw_error_bar = true;

  c1->cd(2);
  
  TH1D *h0 = (TH1D*)fmc->Get("QCD_estimation/QCDest_CombRelIso_1j__QCD");
  TH1D *h1 = (TH1D*)fdata->Get( Form("%s_1j__data", myAES.c_str()) );

  h0->Rebin(rb);
  h1->Rebin(rb);

  h0->GetXaxis()->SetRangeUser(0,up);
  h1->GetXaxis()->SetRangeUser(0,up);

  h0->SetTitle("(b) Control Sample vs QCD in Signal Sample (1j)");
  h1->SetTitle("(b) Control Sample vs QCD in Signal Sample (1j)");

  h0->SetName("QCD in Signal Sample (1j)");
  h1->SetName("Control Sample");

  h0->SetXTitle("RelIso");
  h1->SetXTitle("RelIso");

  //control, light brown
  h1->SetFillColor(41);

  // QCD in signal sample
  h0->SetLineColor(kRed);
  h0->SetMarkerColor(kRed);
  h0->SetMarkerStyle(20);


  //h0->SetLineWidth(2);
  //h1->SetLineWidth(2);
  

  // Normalized to 1 (exclude under/overflows)
  h0->Scale(1.0/h0->GetSumOfWeights());
  h1->Scale(1.0/h1->GetSumOfWeights());

  // Get stat box
  h0->Draw();
  gPad->Update();
  
  if(debug) {
    string mao;
    cout << "stop ";
    cin >> mao;
    cout << "------"<< endl;
    h0->GetListOfFunctions()->ls();
    cout << "------"<< endl;
  }
  TPaveStats *s1 = (TPaveStats*)h0->FindObject("stats");
  if(debug){
    cout << "h0" << endl;
    h0->ls();
    s1->ls();
  }

  h1->Draw();
  gPad->Update();
  TPaveStats *s2 = (TPaveStats*)h1->FindObject("stats");
  if(debug){
    cout << "h1" << endl;
    h1->ls();
    s2->ls();
  }

  double xp = 0.98; // x of top righ point
  double yp = 0.9;  // y of top right point
  double bx = 0.2;  // box width
  double by = 0.2;  // box height
  double s1pos[4] = { xp-bx,   yp-by,    xp,   yp    }; 
  double s2pos[4] = { xp-bx,   yp-2*by,  xp,   yp-by }; 

    
  s1->SetX1NDC(s1pos[0]);
  s1->SetY1NDC(s1pos[1]);
  s1->SetX2NDC(s1pos[2]);
  s1->SetY2NDC(s1pos[3]);
  
  s2->SetX1NDC(s2pos[0]);
  s2->SetY1NDC(s2pos[1]);
  s2->SetX2NDC(s2pos[2]);
  s2->SetY2NDC(s2pos[3]);
  

  double max0 = h0->GetMaximum();
  double max1 = h1->GetMaximum();
  cout << "max of h0 (True QCD): " << max0 << endl;
  cout << "max of h1 (Control):  " << max1 << endl;

  if (max1 > max0) { //draw CS first
    if(draw_error_bar){

      //h1n->SetMaximum(0.3);
      h1->Draw("ah");
      h1->Draw("ae same");
      h0->Draw("ahe same");

      //h0n->Draw("ae same");
      //h1n->Draw("ah same");
      h1->Draw("axis same");
    } else {
      h1->Draw("ahist");
      h0->Draw("ahist same");
      h1->Draw("ahist same");
      h1->Draw("axis same");
      s1->Draw();
    }
  } else {
    if(draw_error_bar){
      h0->Draw("ah");
      h1->Draw("ah same");
      h0->Draw("ae same");
      h1->Draw("ah same");
      h1->Draw("axis same");
    } else {
      h0->Draw("ahist");
      h1->Draw("ahist same");
      h1->Draw("axis same");
      s1->Draw();
    }   
  }

  TLegend leg2(0.32,0.64,0.8,0.88);
  leg2.AddEntry(h1, "Control Sample",           "f"  );
  leg2.AddEntry(h0, "QCD in Signal Sample (1j)","lep");
  leg2.SetFillStyle(0);
  leg2.SetBorderSize(0);
  leg2.Draw();





  //-------------------------------------
  // PLOT 3: QCD in Signal Sample (NES)
  //-------------------------------------
  cout << "Amituofo: Pad 3" << endl;
  c1->cd(3);

  //  gStyle->SetTitleW(0.6);

  TH1D *QCD_1j  = (TH1D*)fmc->Get("QCD_estimation/QCDest_CombRelIso_1j__QCD")->Clone();
  TH1D *QCD_2j  = (TH1D*)fmc->Get("QCD_estimation/QCDest_CombRelIso_2j__QCD")->Clone();
  TH1D *QCD_3j  = (TH1D*)fmc->Get("QCD_estimation/QCDest_CombRelIso_3j__QCD")->Clone();
  TH1D *QCD_4mj = (TH1D*)fmc->Get("QCD_estimation/QCDest_CombRelIso_4mj__QCD")->Clone();
  TH1D *control = (TH1D*)h1->Clone();
 
  QCD_1j->Rebin(rb);
  QCD_2j->Rebin(rb);
  QCD_3j->Rebin(rb);
  QCD_4mj->Rebin(rb);

  QCD_1j->GetXaxis()->SetRangeUser(0,up);
  QCD_2j->GetXaxis()->SetRangeUser(0,up);
  QCD_3j->GetXaxis()->SetRangeUser(0,up);
  QCD_4mj->GetXaxis()->SetRangeUser(0,up);

  QCD_1j->SetLineColor(kRed);
  QCD_2j->SetLineColor(kBlue);
  QCD_3j->SetLineColor(kGreen+1);
  QCD_4mj->SetLineColor(kBlack);

  QCD_1j->SetMarkerColor(kRed);
  QCD_2j->SetMarkerColor(kBlue);
  QCD_3j->SetMarkerColor(kGreen+1);
  QCD_4mj->SetMarkerColor(kBlack);

  QCD_1j->SetMarkerStyle(20);
  QCD_2j->SetMarkerStyle(21);
  QCD_3j->SetMarkerStyle(22);
  QCD_4mj->SetMarkerStyle(23);

  QCD_4mj->SetTitle("QCD in Signal Sample");

  QCD_1j->SetXTitle("RelIso");
  QCD_2j->SetXTitle("RelIso");
  QCD_3j->SetXTitle("RelIso");
  QCD_4mj->SetXTitle("RelIso");
  control->SetXTitle("RelIso");

  QCD_1j->SetTitle("(c) QCD in Signal Sample");
  QCD_2j->SetTitle("(c) QCD in Signal Sample");
  QCD_3j->SetTitle("(c) QCD in Signal Sample");
  QCD_4mj->SetTitle("(c) QCD in Signal Sample");
  control->SetTitle("(c) QCD in Signal Sample");

  gStyle->SetOptStat(0);

  // control sample
  control->SetFillColor(41);

  // Normalized to 1 (exclude under/overflows)
  QCD_1j->Scale(1.0/QCD_1j->GetSumOfWeights());
  QCD_2j->Scale(1.0/QCD_2j->GetSumOfWeights());
  QCD_3j->Scale(1.0/QCD_3j->GetSumOfWeights());
  QCD_4mj->Scale(1.0/QCD_4mj->GetSumOfWeights());
  control->Scale(1.0/control->GetSumOfWeights());

  control->SetMaximum(0.4); //y max (0.3)

  control->Draw("ahist");
  control->Draw("ae same");
  QCD_4mj->Draw("ae1 same");
  QCD_3j->Draw("ae1 same");
  QCD_2j->Draw("ae1 same");
  QCD_1j->Draw("ae1 same");
  QCD_1j->Draw("axis same");


  TLegend leg3(0.4,0.5,0.98,0.9);
  leg3.SetFillColor(0);
  leg3.AddEntry(control,"Control Sample",                   "f" );
  leg3.AddEntry(QCD_1j, "QCD in Signal Sample (1 jet)",     "lp");
  leg3.AddEntry(QCD_2j, "QCD in Signal Sample (2 jets)",    "lp");
  leg3.AddEntry(QCD_3j, "QCD in Signal Sample (3 jets)",    "lp");
  leg3.AddEntry(QCD_4mj,"QCD in Signal Sample (#geq4 jets)","lp");
  leg3.Draw();


  //---------
  //  Save
  //---------

  c1->SaveAs(Form("%s.gif",out.c_str()));
  c1->Close();

  // make pdf
  /*
  c1->SaveAs(Form("%s.eps",out.c_str()));
  gROOT->ProcessLine(Form(".!ps2pdf -dEPSCrop %s.eps",out));
  gROOT->ProcessLine(Form(".!rm %s.eps",out));
  */


}
