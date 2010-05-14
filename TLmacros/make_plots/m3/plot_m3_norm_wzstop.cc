{
  // 9 Mar 09: make m3 plots for the note
  // plot m3 (normalized) for W+jets, Z+jets, QCD, VQQ and single top



  gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.06);
  //gStyle->SetLineScalePS(2);  
  gStyle->SetHistLineWidth(3);
  gStyle->SetMarkerSize(1.5);
  gROOT->ForceStyle();


  //  TFile f("~/cms/analysis/results/18May09_FB_new_ttjet/run_mixture_18May.root");
  TFile f("../../test_mc_mixture.root");



  Wjets_estimation.cd();


  // TH1D *hs = (TH1D*)gDirectory->Get("m3_ttbar");
  TH1D *hwj = (TH1D*)gDirectory->Get("m3__wj");
  TH1D *hqcd = (TH1D*)gDirectory->Get("m3__QCD");
  TH1D *hzj = (TH1D*)gDirectory->Get("m3__zj");
  //TH1D *hvqq = (TH1D*)gDirectory->Get("m3__vqq");
  TH1D *hstop = (TH1D*)gDirectory->Get("m3__singleTop");


  const int rb = 80;
  hwj->Rebin(rb);
  hqcd->Rebin(rb);
  hzj->Rebin(rb);
  //hvqq->Rebin(rb);
  hstop->Rebin(rb);

  // set titles
  hwj->SetTitle("M3 (stacked plot)");
  //  hwj->SetTitle("m3 (summer08)");
  //  hqcd->SetTitle("m3 (summer08)");
  
  hwj->GetXaxis()->SetTitle("M3 (GeV)");
  hzj->GetXaxis()->SetTitle("M3 (GeV)");
  /*
  hs->GetXaxis()->SetRangeUser(0,600);
  hwj->GetXaxis()->SetRangeUser(0,600);
  hqcd->GetXaxis()->SetRangeUser(0,600);
  */
  // Color scheme
  // Signal : blue
  // QCD    : red
  // W+jets : green  
  //----------------------------------------------------
  // light color (for filled histo)
  Color_t Col_sig = kRed-7;//kBlue-7; +1
  Color_t Col_qcd = kAzure-9 ;//kAzure+9;//kGreen-6
  Color_t Col_wj  = kGreen-9; //kGreen-6;//kAzure+9;
  Color_t Col_zj  = kYellow-7;
  Color_t Col_vqq = kGray+1;
  Color_t Col_stop = kViolet-4;
  // darker (for line)
  Color_t Col_sig2 = kRed-4;
  Color_t Col_qcd2 = kAzure+2;
  Color_t Col_wj2  = kGreen+2;
  Color_t Col_zj2  = kOrange-1;
  Color_t Col_vqq2 = kGray+1;
  Color_t Col_stop2 = kViolet;
  //----------------------------------------------------


  //  TCanvas c1("c1","m3 (summer08)",750,500);
  TCanvas c1("c1","M3",600,400);


 hwj->SetFillColor(0);
 hqcd->SetFillColor(0);
 hzj->SetFillColor(0);
 hstop->SetFillColor(0);

 // set titles
 hwj->SetTitle("M3 (normalized)");
 hwj->GetXaxis()->SetTitle("M3 (GeV)");
 hqcd->GetXaxis()->SetTitle("M3 (GeV)");



 hwj->SetLineColor(Col_wj2);
 hqcd->SetLineColor(Col_qcd2);
 hzj->SetLineColor(Col_zj2);
 //hvqq->SetLineColor(Col_vqq2);
 hstop->SetLineColor(Col_stop2);

 hwj->SetMarkerColor(Col_wj2);
 hqcd->SetMarkerColor(Col_qcd2);
 hzj->SetMarkerColor(Col_zj2);
 //hvqq->SetMarkerColor(Col_vqq2);
 hstop->SetMarkerColor(Col_stop2);

 int lw = 3; //line width
 hwj->SetLineWidth(lw);
 hqcd->SetLineWidth(lw);
 hzj->SetLineWidth(lw);
 //hvqq->SetLineWidth(lw);
 hstop->SetLineWidth(lw);

 //hwj->SetLineStyle(lw);
 hzj->SetLineStyle(kDashed);
 hstop->SetLineStyle(kDotted);
 //hvqq->SetLineWidth(lw);
 //hstop->SetLineWidth(lw);

 hwj->SetMarkerStyle(20);
 hzj->SetMarkerStyle(21);
 hstop->SetMarkerStyle(22);
 



 hzj->SetTitle("M3 (normalized to unit area)");

 hzj->DrawNormalized("h");
 hwj->DrawNormalized("h same");
 hstop->DrawNormalized("h same");
 //hvqq->DrawNormalized("same");
 // hstop->DrawNormalized("same");

 // add legend
 // ----------
 TLegend *leg = new TLegend(0.5,0.6,0.9,0.84);  // x1,y1,x2,y2
 //leg->SetTextSize(0.04);
 leg->SetFillStyle(0);
 leg->SetBorderSize(0);

 leg->AddEntry(hwj,  "W+jets",    "lep");
 leg->AddEntry(hzj,  "Z+jets",    "lep");
 leg->AddEntry(hstop,"Single top","lep");
 //leg->AddEntry(hqcd,"QCD","l");
 //leg->AddEntry(hvqq,"VQQ","l");

 leg->SetTextFont(42);

 leg->Draw();

 

 string out = "m3_norm_wzstop.eps";
 c1.SaveAs(Form("%s",out.c_str()));

 gROOT->ProcessLine(Form(".!ps2pdf -dEPSCrop %s",out.c_str()));
 gROOT->ProcessLine(Form(".!rm -f %s",out.c_str()));

}
