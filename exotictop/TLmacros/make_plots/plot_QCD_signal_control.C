//
// 29 Mar 2010
//
// Comparison of True QCD in signal region to Control Region (CR):
// A) take QCD CR from EME+BCE (has isolation bias).
// B) take QCD control sample from 1j dijet_pt15 sample (no isolation bias).
//
// Plots RelIso:
// (1a) True QCD in signal region in 1,2,3,4mj
// (1b) Control shape in B3_e20 for 1,2,3,4mj
// (2a) sig 1j: Control_1j vs true
// (2b) sig 2j: Control_1j vs true
// (2c) sig 3j: Control_1j vs true
// (2d) sig 4mj: Control_1j vs true
//
//--------------------------------------------------
bool use_dijet_cs = true; //cs=control sample
bool use_1T = true;
string template_jetbin = "1j";

bool superimpose_control_234j = true;
string out;//output name

void plot_QCD_signal_control(){
  plot_using_enrich_cs();
  plot_using_dijet_cs();
}

void plot_using_enrich_cs(){
  // 4 templates
  use_dijet_cs = false;
  use_1T = false;
  out = "QCD_signal_control_enri_4T.eps";
  plot();

  // 1 template (2j)
  use_1T = true;
  template_jetbin = "2j";
  out = "QCD_signal_control_enri_1T2j.eps";
  plot();

  // 1 template (allj)
  use_1T = true;
  template_jetbin = "allj";
  out = "QCD_signal_control_enri_1Tallj.eps";
  plot();

}

void plot_using_dijet_cs(){
  use_dijet_cs = true;

  // use 4 templates (1-4j)
  use_1T = false;
  out = "QCD_signal_control_dijet_4T.eps";
  plot();

  // use 1 template (1j)
  use_1T = true;
  template_jetbin = "1j";
  out = "QCD_signal_control_dijet_1T1j.eps";
  plot();
  out = "QCD_signal_control_dijet_1T1j_b.eps";
  superimpose_control_234j = false;
  plot();

  // use 1 template (allj)
  template_jetbin = "allj";
  superimpose_control_234j = true;
  out = "QCD_signal_control_dijet_1Tallj.eps";
  plot();

}


bool plot(){

  double bw = 0.1;
  int rb = bw/0.01;

  double up = 1.5;


  TFile *f = new TFile("../test_mc_mixture.root");//NES
  TFile *f2(0);
  if(use_dijet_cs)
    f2 = new TFile("~/cms/analysis/results/2010_Feb/25Feb_7TeV_plainQCD_pt15/test_data.root");//Control
  else
    f2 = new TFile("../test_mc_mixture.root");//NES
  
  if(f->IsZombie()) return 0;
  if(f2->IsZombie()) return 0;
  

  gStyle->SetTitleH(0.07);
  gStyle->SetTitleW(0.9);
  //gStyle->SetTitleAlign(13); //Def in style =11 left-bottom
  gStyle->SetOptStat(0);
  //  gStyle->SetMarkerSize(1);


  //  TCanvas *c1 = new TCanvas("c1","RelIso for QCD", 4*300,2*300);
  TCanvas *c1 = new TCanvas("c1","RelIso for QCD", 4*300,2*300);
  c1->Divide(4,2);
  c1->cd(1);

  //-------------------------------------
  // PLOT 1: QCD in Signal Sample (NES)
  //-------------------------------------



  TH1D *QCD_1j  = (TH1D*)f->Get("QCD_estimation/QCDest_CombRelIso_1j__QCD")->Clone();
  TH1D *QCD_2j  = (TH1D*)f->Get("QCD_estimation/QCDest_CombRelIso_2j__QCD")->Clone();
  TH1D *QCD_3j  = (TH1D*)f->Get("QCD_estimation/QCDest_CombRelIso_3j__QCD")->Clone();
  TH1D *QCD_4mj = (TH1D*)f->Get("QCD_estimation/QCDest_CombRelIso_4mj__QCD")->Clone();
  
  QCD_1j->Rebin(rb);
  QCD_2j->Rebin(rb);
  QCD_3j->Rebin(rb);
  QCD_4mj->Rebin(rb);

  QCD_1j->GetXaxis()->SetRangeUser(0,up);
  QCD_2j->GetXaxis()->SetRangeUser(0,up);
  QCD_3j->GetXaxis()->SetRangeUser(0,up);
  QCD_4mj->GetXaxis()->SetRangeUser(0,up);

  Color_t c1j = kRed;
  Color_t c2j = kBlue;
  Color_t c3j = kGreen+1;
  Color_t c4j = kViolet-1;
  Color_t cAllj = kOrange+1;

  QCD_1j->SetLineColor(c1j);
  QCD_2j->SetLineColor(c2j);
  QCD_3j->SetLineColor(c3j);
  QCD_4mj->SetLineColor(c4j);//black

  QCD_1j->SetMarkerColor(c1j);
  QCD_2j->SetMarkerColor(c2j);
  QCD_3j->SetMarkerColor(c3j);
  QCD_4mj->SetMarkerColor(c4j);//black

  QCD_1j->SetMarkerStyle(20);
  QCD_2j->SetMarkerStyle(21);
  QCD_3j->SetMarkerStyle(22);
  QCD_4mj->SetMarkerStyle(23);

  QCD_1j->SetMarkerSize(0.7);
  QCD_2j->SetMarkerSize(0.7);
  QCD_3j->SetMarkerSize(0.7);
  QCD_4mj->SetMarkerSize(0.7);


  QCD_1j->SetXTitle("RelIso");
  QCD_2j->SetXTitle("RelIso");
  QCD_3j->SetXTitle("RelIso");
  QCD_4mj->SetXTitle("RelIso");

  QCD_1j->SetTitle("True QCD in Signal Sample");
  QCD_2j->SetTitle("True QCD in Signal Sample");
  QCD_3j->SetTitle("True QCD in Signal Sample");
  QCD_4mj->SetTitle("True QCD in Signal Sample");

  // Normalized to 1 (exclude under/overflows)
  QCD_1j->Scale(1.0/QCD_1j->GetSumOfWeights());
  QCD_2j->Scale(1.0/QCD_2j->GetSumOfWeights());
  QCD_3j->Scale(1.0/QCD_3j->GetSumOfWeights());
  QCD_4mj->Scale(1.0/QCD_4mj->GetSumOfWeights());
  /*
  double tallest = TMath::Max(QCD_1j->GetMaximum(),QCD_2j->GetMaximum());
  tallest = TMath::Max(tallest,QCD_3j->GetMaximum());
  tallest = TMath::Max(tallest,QCD_4mj->GetMaximum());
  QCD_1j->SetMaximum(tallest*1.2);
  */
  QCD_1j->SetMaximum(0.18);
  QCD_1j->DrawCopy("e1");
  QCD_2j->DrawCopy("e1 same");
  QCD_3j->DrawCopy("e1 same");
  QCD_4mj->DrawCopy("e1 same");


  //------------------------------
  //  PLOT 2: AES: control sample
  //------------------------------
  c1->cd(2);
  string AES = "B3_e20";

  TH1D *QCD_AES_1j = (TH1D*)f2->Get( Form("QCD_planB/QCDest_CombRelIso_AES_plan%s_1j__data",AES.c_str()));
  TH1D *QCD_AES_2j = (TH1D*)f2->Get( Form("QCD_planB/QCDest_CombRelIso_AES_plan%s_2j__data",AES.c_str()));
  TH1D *QCD_AES_3j = (TH1D*)f2->Get( Form("QCD_planB/QCDest_CombRelIso_AES_plan%s_3j__data",AES.c_str()));
  TH1D *QCD_AES_4mj = (TH1D*)f2->Get( Form("QCD_planB/QCDest_CombRelIso_AES_plan%s_4mj__data",AES.c_str()));
  TH1D *QCD_AES_allj = (TH1D*)f2->Get( Form("QCD_planB/QCDest_CombRelIso_AES_plan%s_allj__data",AES.c_str()));



  QCD_AES_1j->Rebin(rb);
  QCD_AES_2j->Rebin(rb);
  QCD_AES_3j->Rebin(rb);
  QCD_AES_4mj->Rebin(rb);
  QCD_AES_allj->Rebin(rb);

  QCD_AES_1j->GetXaxis()->SetRangeUser(0,up);
  QCD_AES_2j->GetXaxis()->SetRangeUser(0,up);
  QCD_AES_3j->GetXaxis()->SetRangeUser(0,up);
  QCD_AES_4mj->GetXaxis()->SetRangeUser(0,up);
  QCD_AES_allj->GetXaxis()->SetRangeUser(0,up);

  QCD_AES_1j->SetLineColor(c1j);
  QCD_AES_2j->SetLineColor(c2j);
  QCD_AES_3j->SetLineColor(c3j);
  QCD_AES_4mj->SetLineColor(c4j);

  QCD_AES_1j->SetMarkerColor(c1j);
  QCD_AES_2j->SetMarkerColor(c2j);
  QCD_AES_3j->SetMarkerColor(c3j);
  QCD_AES_4mj->SetMarkerColor(c4j);
  QCD_AES_allj->SetMarkerColor(cAllj);

  QCD_AES_1j->SetMarkerStyle(20);
  QCD_AES_2j->SetMarkerStyle(21);
  QCD_AES_3j->SetMarkerStyle(22);
  QCD_AES_4mj->SetMarkerStyle(23);
  QCD_AES_allj->SetMarkerStyle(20);

  QCD_AES_1j->SetMarkerSize(0.7);
  QCD_AES_2j->SetMarkerSize(0.7);
  QCD_AES_3j->SetMarkerSize(0.7);
  QCD_AES_4mj->SetMarkerSize(0.7);
  QCD_AES_allj->SetMarkerSize(0.7);


  QCD_AES_1j->SetXTitle("RelIso");
  QCD_AES_2j->SetXTitle("RelIso");
  QCD_AES_3j->SetXTitle("RelIso");
  QCD_AES_4mj->SetXTitle("RelIso");
  QCD_AES_allj->SetXTitle("RelIso");

  char *title_AES = "QCD Control Sample (all events)";
  if(use_dijet_cs) title_AES="QCD Control Sample (Dijet pt15)";

  QCD_AES_1j->SetTitle(title_AES);
  QCD_AES_2j->SetTitle(title_AES);
  QCD_AES_3j->SetTitle(title_AES);
  QCD_AES_4mj->SetTitle(title_AES);


  // Normalized to 1 (exclude under/overflows)
  QCD_AES_1j->Scale(1.0/QCD_AES_1j->GetSumOfWeights());
  QCD_AES_2j->Scale(1.0/QCD_AES_2j->GetSumOfWeights());
  QCD_AES_3j->Scale(1.0/QCD_AES_3j->GetSumOfWeights());
  QCD_AES_4mj->Scale(1.0/QCD_AES_4mj->GetSumOfWeights());
  QCD_AES_allj->Scale(1.0/QCD_AES_allj->GetSumOfWeights());

  QCD_AES_1j->SetMaximum(0.18);
  QCD_AES_2j->SetMaximum(0.18);
  QCD_AES_3j->SetMaximum(0.18);
  QCD_AES_4mj->SetMaximum(0.18);
  QCD_AES_allj->SetMaximum(0.18);
 

  QCD_AES_1j->DrawCopy("e1");
  if(superimpose_control_234j){
    QCD_AES_2j->DrawCopy("e1 same");
    QCD_AES_3j->DrawCopy("e1 same");
    QCD_AES_4mj->DrawCopy("e1 same");
  }

  

  TLegend leg(0.6,0.6,0.9,0.9);
  leg.SetFillColor(0);
  leg.AddEntry(QCD_1j, "1 jet","lp");
  leg.AddEntry(QCD_2j, "2 jets","lp");
  leg.AddEntry(QCD_3j, "3 jets","lp");
  leg.AddEntry(QCD_4mj,"#geq4 jets","lp");
  leg.Draw();
  c1->cd(1); leg.Draw();
  //c1->cd(2); leg.Draw();




  //---------------------------------
  // Plot 3: control (allj)
  //---------------------------------
  c1->cd(3);

  if(use_dijet_cs) QCD_AES_allj->SetTitle("QCD CS (Dijet pt15, allj)");
  else             QCD_AES_allj->SetTitle("QCD CS (All events, allj)");
  QCD_AES_allj->DrawCopy();
  



  //-----------------------------------------------
  // 2nd Row Plots: control vs signal (1j,2j,3j,4j)
  //------------------------------------------------

  // Formatting for control
  //------------------------
  QCD_AES_1j->SetMinimum(0);
  QCD_AES_2j->SetMinimum(0);
  QCD_AES_3j->SetMinimum(0);
  QCD_AES_4mj->SetMinimum(0);

  QCD_AES_1j->SetFillColor(41);
  QCD_AES_1j->SetLineColor(kBlack);
  QCD_AES_1j->SetMarkerColor(kBlack);
  QCD_AES_1j->SetTitle("Shape comparison: 1j");

  QCD_AES_2j->SetFillColor(41);
  QCD_AES_2j->SetLineColor(kBlack);
  QCD_AES_2j->SetMarkerColor(kBlack);
  QCD_AES_2j->SetTitle("Shape comparison: 2j");

  QCD_AES_3j->SetFillColor(41);
  QCD_AES_3j->SetLineColor(kBlack);
  QCD_AES_3j->SetMarkerColor(kBlack);
  QCD_AES_3j->SetTitle("Shape comparison: 3j");

  QCD_AES_4mj->SetFillColor(41);
  QCD_AES_4mj->SetLineColor(kBlack);
  QCD_AES_4mj->SetMarkerColor(kBlack);
  QCD_AES_4mj->SetTitle("Shape comparison: #geq4j");

  QCD_AES_allj->SetFillColor(41);
  QCD_AES_allj->SetLineColor(kBlack);
  QCD_AES_allj->SetMarkerColor(kBlack);
  QCD_AES_allj->SetTitle("Shape comparison: #geq4j");

  // If using one template
  TH1D *QCD_AES_temp;
  if(use_1T){
    if     (template_jetbin=="1j")   QCD_AES_temp = (TH1D*)QCD_AES_1j->Clone("QCD_AES_temp");
    else if(template_jetbin=="2j")   QCD_AES_temp = (TH1D*)QCD_AES_2j->Clone("QCD_AES_temp");
    else if(template_jetbin=="allj") QCD_AES_temp = (TH1D*)QCD_AES_allj->Clone("QCD_AES_temp");
  }

  // 2nd Row: Plot 1: "True QCD in 1 jet"  vs  "Control"
  c1->cd(5);
  if(use_1T){
    QCD_AES_temp->SetTitle("Shape comparison: 1j");
    QCD_AES_temp->DrawCopy("hist");
    QCD_AES_temp->DrawCopy("e same");    
  } else {
    QCD_AES_1j->DrawCopy("hist");
    QCD_AES_1j->DrawCopy("e same");
  }
  QCD_1j->Draw("he same");


  // 2nd Row: Plot 2: "True QCD in 2 jet"  vs  "Control"
  c1->cd(6);
  if(use_1T){
    QCD_AES_temp->SetTitle("Shape comparison: 1j");
    QCD_AES_temp->DrawCopy("hist");
    QCD_AES_temp->DrawCopy("e same");
  }else{
    QCD_AES_2j->Draw("hist");
    QCD_AES_2j->Draw("e same");
  }
  QCD_2j->Draw("he same");


  // 2nd Row: Plot 3: "True QCD in 3 jet"  vs  "Control"
  c1->cd(7);
  if(use_1T){
    QCD_AES_temp->SetTitle("Shape comparison: 3j");
    QCD_AES_temp->DrawCopy("hist");
    QCD_AES_temp->DrawCopy("e same");    
  }else{
    QCD_AES_3j->Draw("hist");
    QCD_AES_3j->Draw("e same");
  }
  QCD_3j->Draw("he same");

  // 2nd Row: Plot 4: "True QCD in >=4 jet"  vs  "Control"
  c1->cd(8);
  if(use_1T){
    QCD_AES_temp->SetTitle("Shape comparison: #geq4j");
    QCD_AES_temp->DrawCopy("hist");
    QCD_AES_temp->DrawCopy("e same");    
  }else{
    QCD_AES_4mj->Draw("hist");
    QCD_AES_4mj->Draw("e same");
  }
  QCD_4mj->Draw("he same");

  

  TLegend leg3(0.5,0.7,0.99,0.9);
  /*
  TLegend leg4T_1j(0.5,0.7,0.99,0.9);
  TLegend leg4T_2j(0.5,0.7,0.99,0.9);
  TLegend leg4T_3j(0.5,0.7,0.99,0.9);
  TLegend leg4T_4j(0.5,0.7,0.99,0.9);
  TLegend leg1T(0.5,0.7,0.99,0.9);

  leg4T_1j.SetFillColor(0);
  leg4T_2j.SetFillColor(0);
  leg4T_3j.SetFillColor(0);
  leg4T_4j.SetFillColor(0);
  leg1T.SetFillColor(0);
  
  // leg: 4 templates (1j)
  leg4T_1j.AddEntry(QCD_AES_1j, "Control Sample (1j)","f" );
  leg4T_1j.AddEntry(QCD_1j,     "True QCD (1j)",      "lp");

  // leg: 4 templates (2j)
  leg4T_2j.AddEntry(QCD_AES_2j, "Control Sample (2j)","f" );
  leg4T_2j.AddEntry(QCD_2j,     "True QCD (2j)",      "lp");

  // leg: 4 templates (3j)
  leg4T_3j.AddEntry(QCD_AES_3j, "Control Sample (3j)","f" );
  leg4T_3j.AddEntry(QCD_3j,     "True QCD (3j)",      "lp");

  // leg: 4 templates (4mj)
  leg4T_4j.AddEntry(QCD_AES_4mj, "Control Sample (#geq4j)","f" );
  leg4T_4j.AddEntry(QCD_4mj,     "True QCD (#geq4j)",      "lp");

  // leg: 1 template
  leg1T_1j.AddEntry(QCD_AES_1j, "Control Sample (1j)","f" );
  leg1T_1j.AddEntry(QCD_1j,     "True QCD (1j)",      "lp");

  // leg: 1 template
  leg1T_2j.AddEntry(QCD_AES_2j, "Control Sample (2j)","f" );
  leg1T_2j.AddEntry(QCD_2j,     "True QCD (2j)",      "lp");

  // leg: 1 template
  leg1T_3j.AddEntry(QCD_AES_allj, "Control Sample (2j)","f" );
  leg1T_3j.AddEntry(QCD_3j,     "True QCD (2j)",      "lp");
*/

  leg3.SetFillColor(0);

  char *oneTemplateLabel = Form("Control Sample (%s)",template_jetbin.c_str());


  c1->cd(5);
  if(use_1T) 
    leg3.AddEntry(QCD_AES_1j, oneTemplateLabel,     "f");
  else
    leg3.AddEntry(QCD_AES_1j, "Control Sample (1j)","f");
  leg3.AddEntry(QCD_1j,     "True QCD (1j)","lp");
  leg3.DrawClone();

  c1->cd(6);
  leg3.Clear();
  if(use_1T) 
    leg3.AddEntry(QCD_AES_1j, oneTemplateLabel,     "f");
  else 
    leg3.AddEntry(QCD_AES_2j, "Control Sample (2j)","f");
  leg3.AddEntry(QCD_2j,     "True QCD (2j)","lp");
  leg3.DrawClone();

  c1->cd(7);
  leg3.Clear();
  if(use_1T)
    leg3.AddEntry(QCD_AES_1j, oneTemplateLabel,     "f");
  else
    leg3.AddEntry(QCD_AES_3j, "Control Sample (3j)","f");
  leg3.AddEntry(QCD_3j,     "True QCD (3j)","lp");
  leg3.DrawClone();

  c1->cd(8);
  leg3.Clear();
  if(use_1T)
    leg3.AddEntry(QCD_AES_1j, oneTemplateLabel,     "f");
  else      
    leg3.AddEntry(QCD_AES_4mj, "Control Sample (#geq4j)","f");
  leg3.AddEntry(QCD_4mj,     "True QCD (#geq4j)","lp");
  leg3.DrawClone();


  c1->SaveAs(Form("%s",out));
  c1->Close();

  
  gROOT->ProcessLine(Form(".!ps2pdf -dEPSCrop %s",out));
  gROOT->ProcessLine(Form(".!rm -f %s",out));


}

