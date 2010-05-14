void plot_iso_QCD_sig_reg_superimpose_1234j() {

  TFile f("../test_mc_mixture.root");
  
  TH1D *QCD_1j = (TH1D*)f.Get("QCD_estimation/QCDest_CombRelIso_1j__QCD");
  TH1D *QCD_2j = (TH1D*)f.Get("QCD_estimation/QCDest_CombRelIso_2j__QCD");
  TH1D *QCD_3j = (TH1D*)f.Get("QCD_estimation/QCDest_CombRelIso_3j__QCD");
  TH1D *QCD_4mj = (TH1D*)f.Get("QCD_estimation/QCDest_CombRelIso_4mj__QCD");

  int rb = 10;
  QCD_1j->Rebin(rb);
  QCD_2j->Rebin(rb);
  QCD_3j->Rebin(rb);
  QCD_4mj->Rebin(rb);

  double up = 1.5;
  QCD_1j->GetXaxis()->SetRangeUser(0,up);
  QCD_2j->GetXaxis()->SetRangeUser(0,up);
  QCD_3j->GetXaxis()->SetRangeUser(0,up);
  QCD_4mj->GetXaxis()->SetRangeUser(0,up);

  QCD_1j->SetLineColor(kBlack);
  QCD_2j->SetLineColor(kBlue);
  QCD_3j->SetLineColor(kGreen);
  QCD_4mj->SetLineColor(kRed);

  QCD_1j->SetFillColor(kBlack);
  QCD_2j->SetFillColor(kBlue);
  QCD_3j->SetFillColor(kGreen);
  QCD_4mj->SetFillColor(kRed);

  QCD_1j->SetMarkerColor(kBlack);
  QCD_2j->SetMarkerColor(kBlue);
  QCD_3j->SetMarkerColor(kGreen);
  QCD_4mj->SetMarkerColor(kRed);

  QCD_1j->SetMarkerStyle(20);
  QCD_2j->SetMarkerStyle(21);
  QCD_3j->SetMarkerStyle(22);
  QCD_4mj->SetMarkerStyle(23);

  TCanvas c1("c1","RelIso of QCD in signal region (superimpose 1,2,3,4j)",500,500);

  gStyle->SetOptStat(0);
  
  QCD_4mj->SetTitle("QCD in signal region");
  QCD_4mj->SetXTitle("RelIso");

  QCD_4mj->DrawNormalized("e1");
  QCD_3j->DrawNormalized("e1 same");
  QCD_2j->DrawNormalized("e1 same");
  QCD_1j->DrawNormalized("e1 same");

  TLegend leg(0.5,0.6,0.9,0.9);
  leg.SetFillColor(0);
  leg.AddEntry(QCD_1j,"1 jet","lp");
  leg.AddEntry(QCD_2j,"2 jets","lp");
  leg.AddEntry(QCD_3j,"3 jets","lp");
  leg.AddEntry(QCD_4mj,"#geq4 jets","lp");
  leg.Draw();

  c1.SaveAs("iso_QCD_sig_reg_superimpose_1234j.gif");
  c1.Close();
}
