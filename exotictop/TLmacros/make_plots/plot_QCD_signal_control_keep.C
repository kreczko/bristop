//
// 4-2-2010
// Plots RelIso:
// (1a) True QCD in signal region in 1,2,3,4mj
// (1b) Control shape in B3_e20 for 1,2,3,4mj
// (2a) 1j: Control vs true
// (2b) 2j: Control vs true
// (2c) 3j: Control vs true
// (2d) 4mj: Control vs true
//
//--------------------------------------------------
void plot_QCD_signal_control(){

  double bw = 0.1;
  int rb = bw/0.01;

  double up = 1.5;


  TFile *f = new TFile("../test_mc_mixture.root");
  


  gStyle->SetTitleH(0.07);
  gStyle->SetTitleW(0.9);
  //gStyle->SetTitleAlign(13); //Def in style =11 left-bottom
  gStyle->SetOptStat(0);


  TCanvas *c1 = new TCanvas("c1","RelIso for QCD", 4*300,2*300);
  c1->Divide(4,2);
  c1->cd(1);

  //-------------------------------------
  // PLOT 1: QCD in Signal Sample (NES)
  //-------------------------------------



  TH1D *QCD_1j = (TH1D*)f->Get("QCD_estimation/QCDest_CombRelIso_1j__QCD")->Clone();
  TH1D *QCD_2j = (TH1D*)f->Get("QCD_estimation/QCDest_CombRelIso_2j__QCD")->Clone();
  TH1D *QCD_3j = (TH1D*)f->Get("QCD_estimation/QCDest_CombRelIso_3j__QCD")->Clone();
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


  QCD_1j->SetXTitle("RelIso");
  QCD_2j->SetXTitle("RelIso");
  QCD_3j->SetXTitle("RelIso");
  QCD_4mj->SetXTitle("RelIso");

  QCD_1j->SetTitle("QCD in Signal Sample");
  QCD_2j->SetTitle("QCD in Signal Sample");
  QCD_3j->SetTitle("QCD in Signal Sample");
  QCD_4mj->SetTitle("QCD in Signal Sample");

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

  TH1D *QCD_AES_1j = (TH1D*)f->Get( Form("QCD_planB/QCDest_CombRelIso_AES_plan%s_1j__data",AES.c_str()));
  TH1D *QCD_AES_2j = (TH1D*)f->Get( Form("QCD_planB/QCDest_CombRelIso_AES_plan%s_2j__data",AES.c_str()));
  TH1D *QCD_AES_3j = (TH1D*)f->Get( Form("QCD_planB/QCDest_CombRelIso_AES_plan%s_3j__data",AES.c_str()));
  TH1D *QCD_AES_4mj = (TH1D*)f->Get( Form("QCD_planB/QCDest_CombRelIso_AES_plan%s_4mj__data",AES.c_str()));



  QCD_AES_1j->Rebin(rb);
  QCD_AES_2j->Rebin(rb);
  QCD_AES_3j->Rebin(rb);
  QCD_AES_4mj->Rebin(rb);

  QCD_AES_1j->GetXaxis()->SetRangeUser(0,up);
  QCD_AES_2j->GetXaxis()->SetRangeUser(0,up);
  QCD_AES_3j->GetXaxis()->SetRangeUser(0,up);
  QCD_AES_4mj->GetXaxis()->SetRangeUser(0,up);

  QCD_AES_1j->SetLineColor(c1j);
  QCD_AES_2j->SetLineColor(c2j);
  QCD_AES_3j->SetLineColor(c3j);
  QCD_AES_4mj->SetLineColor(c4j);

  QCD_AES_1j->SetMarkerColor(c1j);
  QCD_AES_2j->SetMarkerColor(c2j);
  QCD_AES_3j->SetMarkerColor(c3j);
  QCD_AES_4mj->SetMarkerColor(c4j);

  QCD_AES_1j->SetMarkerStyle(20);
  QCD_AES_2j->SetMarkerStyle(21);
  QCD_AES_3j->SetMarkerStyle(22);
  QCD_AES_4mj->SetMarkerStyle(23);

  QCD_AES_1j->SetXTitle("RelIso");
  QCD_AES_2j->SetXTitle("RelIso");
  QCD_AES_3j->SetXTitle("RelIso");
  QCD_AES_4mj->SetXTitle("RelIso");

  QCD_AES_1j->SetTitle("QCD Control Sample (all events)");
  QCD_AES_2j->SetTitle("QCD Control Sample (all events)");
  QCD_AES_3j->SetTitle("QCD Control Sample (all events)");
  QCD_AES_4mj->SetTitle("QCD Control Sample (all events)");

  // Normalized to 1 (exclude under/overflows)
  QCD_AES_1j->Scale(1.0/QCD_AES_1j->GetSumOfWeights());
  QCD_AES_2j->Scale(1.0/QCD_AES_2j->GetSumOfWeights());
  QCD_AES_3j->Scale(1.0/QCD_AES_3j->GetSumOfWeights());
  QCD_AES_4mj->Scale(1.0/QCD_AES_4mj->GetSumOfWeights());

  QCD_AES_1j->SetMaximum(0.18);
  QCD_AES_2j->SetMaximum(0.18);
  QCD_AES_3j->SetMaximum(0.18);
  QCD_AES_4mj->SetMaximum(0.18);

  QCD_AES_1j->DrawCopy("e1");
  QCD_AES_2j->DrawCopy("e1 same");
  QCD_AES_3j->DrawCopy("e1 same");
  QCD_AES_4mj->DrawCopy("e1 same");


  

  TLegend leg(0.6,0.6,0.9,0.9);
  leg.SetFillColor(0);
  leg.AddEntry(QCD_1j, "1 jet","lp");
  leg.AddEntry(QCD_2j, "2 jets","lp");
  leg.AddEntry(QCD_3j, "3 jets","lp");
  leg.AddEntry(QCD_4mj,"#geq4 jets","lp");
  leg.Draw();
  c1->cd(1); leg.Draw();
  c1->cd(2); leg.Draw();

  //---------------------------------
  // Plot 3: signal vs control (1j)
  //---------------------------------

  c1->cd(5);
  //  QCD_1j->SetLineColor(kRed);
  QCD_AES_1j->SetFillColor(41);
  QCD_AES_1j->SetLineColor(kBlack);
  QCD_AES_1j->SetMarkerColor(kBlack);
  QCD_AES_1j->SetTitle("Shape comparison: 1j");
  QCD_AES_1j->Draw("hist");

  QCD_AES_1j->Draw("e same");
  QCD_1j->Draw("he same");

  c1->cd(6);
  //QCD_2j->SetLineColor(kBlue);
  QCD_AES_2j->SetFillColor(41);
  QCD_AES_2j->SetLineColor(kBlack);
  QCD_AES_2j->SetMarkerColor(kBlack);
  QCD_AES_2j->SetTitle("Shape comparison: 2j");
  QCD_AES_2j->Draw("hist");
  QCD_AES_2j->Draw("e same");
  QCD_2j->Draw("he same");

  c1->cd(7);
  //QCD_3j->SetLineColor(kGreen+1);
  QCD_AES_3j->SetFillColor(41);
  QCD_AES_3j->SetLineColor(kBlack);
  QCD_AES_3j->SetMarkerColor(kBlack);
  QCD_AES_3j->SetTitle("Shape comparison: 3j");
  QCD_AES_3j->Draw("hist");
  QCD_AES_3j->Draw("e same");
  QCD_3j->Draw("he same");

  c1->cd(8);
  //QCD_4mj->SetLineColor(kViolet+1);
  QCD_AES_4mj->SetFillColor(41);
  QCD_AES_4mj->SetLineColor(kBlack);
  QCD_AES_4mj->SetMarkerColor(kBlack);
  QCD_AES_4mj->SetTitle("Shape comparison: #geq4j");
  QCD_AES_4mj->Draw("hist");
  QCD_AES_4mj->Draw("e same");
  QCD_4mj->Draw("he same");

  
  TLegend leg3(0.55,0.7,0.9,0.9);
  leg3.SetFillColor(0);
  leg3.AddEntry(QCD_AES_1j, "Control Sample","f");
  leg3.AddEntry(QCD_1j,     "True QCD","lp");
  leg3.Draw();
  c1->cd(5); leg3.DrawClone();

  leg3.Clear();
  leg3.AddEntry(QCD_AES_2j, "Control Sample","f");
  leg3.AddEntry(QCD_2j,     "True QCD","lp");
  c1->cd(6); leg3.DrawClone();

  leg3.Clear();
  leg3.AddEntry(QCD_AES_3j, "Control Sample","f");
  leg3.AddEntry(QCD_3j,     "True QCD","lp");
  c1->cd(7); leg3.DrawClone();

  leg3.Clear();
  leg3.AddEntry(QCD_AES_4mj, "Control Sample","f");
  leg3.AddEntry(QCD_4mj,     "True QCD","lp");
  c1->cd(8); leg3.DrawClone();


  c1->SaveAs("QCD_signal_control.gif");
  c1->Close();
}

