{
  // 23 Mar 2010
  // make m3 plots for the note: 2 in 1
  //  stacked (all) ,  normalized (sig,wj)
  // save plot in ps, then convert to pdf

  
  bool fillSolidColour = true;

  gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.06);
  gStyle->SetLineScalePS(2);  

  //  TFile f("~/cms/analysis/results/18May09_FB_new_ttjet/run_mixture_18May.root");
  TFile f("../../test_mc_mixture.root");

  Wjets_estimation.cd();


  TH1D *hs = (TH1D*)gDirectory->Get("m3__ttbar");
  TH1D *hwj = (TH1D*)gDirectory->Get("m3__wj");
  TH1D *hqcd = (TH1D*)gDirectory->Get("m3__QCD");
  TH1D *hzj = (TH1D*)gDirectory->Get("m3__zj");
  //TH1D *hvqq = (TH1D*)gDirectory->Get("m3__vqq");
  TH1D *hstop = (TH1D*)gDirectory->Get("m3__singleTop");

  const int rb = 80;
  hs->Rebin(rb);
  hwj->Rebin(rb);
  hqcd->Rebin(rb);
  hzj->Rebin(rb);
  //hvqq->Rebin(rb);
  hstop->Rebin(rb);

  // set titles
  hs->SetTitle("M3 (stacked plot)");
  hs->GetXaxis()->SetTitle("M3 (GeV)");
  /*
  hs->GetXaxis()->SetRangeUser(0,600);
  hwj->GetXaxis()->SetRangeUser(0,600);
  hqcd->GetXaxis()->SetRangeUser(0,600);
  */
  //----------------------------------------------------
  // Color scheme
  //---------------
  // Signal : blue
  // QCD    : red
  // W+jets : green  
  //------------------
  // light color (for filled histo)
  Color_t Col_sig = kRed-7;//kBlue-7; +1
  Color_t Col_qcd = kAzure-9 ;//kAzure+9;//kGreen-6
  Color_t Col_wj  = kGreen-9; //kGreen-6;//kAzure+9;
  Color_t Col_zj  = kYellow-6;
  Color_t Col_vqq = kGray+1;
  Color_t Col_singletop = kViolet-4;
  // darker (for line)
  Color_t Col_sig2 = kRed-4;
  Color_t Col_qcd2 = kAzure+2;
  Color_t Col_wj2  = kGreen+2;
  Color_t Col_zj2  = kYellow-6;
  Color_t Col_vqq2 = kGray+1;
  Color_t Col_singletop2 = kViolet-4;
  //----------------------------------------------------

  hs->SetFillColor(Col_sig);
  hwj->SetFillColor(Col_wj);
  hqcd->SetFillColor(Col_qcd);
  hzj->SetFillColor(Col_zj);
  hstop->SetFillColor(Col_singletop);



  if(fillSolidColour==false){



    Style_t Sty_sig = 3001;
    Style_t Sty_stop = 1001; //solid
    Style_t Sty_wj = 3354;
    Style_t Sty_zj = 3345;
    Style_t Sty_qcd = 3344;

    hs->SetFillStyle(Sty_sig);
    hwj->SetFillStyle(Sty_wj);
    hqcd->SetFillStyle(Sty_qcd);
    hzj->SetFillStyle(Sty_zj);
    hstop->SetFillStyle(Sty_stop);
  }





  TCanvas c1("c1","m3 (summer08)",1200,400); 
  c1.Divide(2,1);

  c1.cd(1);

  /*
  // draw overlayed plot: not stacking up
  hs->Draw();
  hwj->Draw("same");
  hqcd->Draw("same");
  */
  /*
  gPad->SetFrameLineWidth(1);
  gStyle->SetLineWidth(1);
  gStyle->SetHistLineWidth(1);
  */

  THStack *stack = new THStack("stack","M3 (stacked)");

  //stack->Add(m3_vqq); //first

  /*
  stack->Add(m3_qcd);
  stack->Add(m3_zj);
  stack->Add(m3_wj);
  stack->Add(m3_singletop);
  stack->Add(m3_tt);
  */
  stack->Add(hqcd);
  stack->Add(hzj);
  stack->Add(hwj);
  stack->Add(hstop);
  stack->Add(hs);

  stack->Draw("hist");  
  stack->GetXaxis()->SetTitle("M3 (GeV)");
  TLine line(0,0,700,0);
  line.Draw();
  TLine line2(0,0,0,20);
  line2.Draw();

 TH1D *hs2  = hs->Clone("hs2");
 TH1D *hwj2 = hwj->Clone("hwj2");
 TH1D *hqcd2 = hqcd->Clone("hqcd2");
 TH1D *hzj2  = hzj->Clone("hzj2");
 //TH1D *hvqq2 = hvqq->Clone("hvqq2");
 TH1D *hstop2 = hstop->Clone("hstop2");

 hs2->SetFillColor(0);
 hwj2->SetFillColor(0);
 hqcd2->SetFillColor(0);
 hzj2->SetFillColor(0);
 hstop2->SetFillColor(0);

 // set titles
 hs2->SetTitle("M3 (normalized to unit area)");
 // hwj2->SetTitle("Normalized m3 (Summer08)");
 // hqcd2->SetTitle("Normalized m3 (Summer08)");

 hs2->GetXaxis()->SetTitle("M3 (GeV)");
 hwj2->GetXaxis()->SetTitle("M3 (GeV)");
 hqcd2->GetXaxis()->SetTitle("M3 (GeV)");

 // hs2->SetLineStyle();
 // hwj2->SetLineStyle(kDashed);
 // hqcd2->SetLineStyle(kDotted);

 hs2->SetLineColor(Col_sig2);
 hwj2->SetLineColor(Col_wj2);
 hqcd2->SetLineColor(Col_qcd2);
 hzj2->SetLineColor(Col_zj2);
 //hvqq2->SetLineColor(Col_vqq2);
 hstop2->SetLineColor(Col_singletop2);

 hs2->SetMarkerColor(Col_sig2);
 hwj2->SetMarkerColor(Col_wj2);
 hqcd2->SetMarkerColor(Col_qcd2);
 hzj2->SetMarkerColor(Col_zj2);
 //hvqq2->SetMarkerColor(Col_vqq2);
 hstop2->SetMarkerColor(Col_singletop2);

 int lw = 2;
 hs2->SetLineWidth(lw);
 hwj2->SetLineWidth(lw);
 hqcd2->SetLineWidth(lw);
 hzj2->SetLineWidth(lw);
 //hvqq2->SetLineWidth(lw);
 hstop2->SetLineWidth(lw);

 hwj2->SetLineStyle(kDashed);
 hwj2->SetMarkerStyle(20);
 hs2->SetMarkerStyle(22);

 c1.cd(2);

 hs2->DrawNormalized("h");
 hwj2->DrawNormalized("h same");
 //hqcd2->DrawNormalized("same");
 //hzj2->DrawNormalized("same");
 //hvqq2->DrawNormalized("same");
 // hstop2->DrawNormalized("same");

 // add legend
 // ----------
 TLegend *leg = new TLegend(0.54,0.52,0.88,0.85);  // x1,y1,x2,y2
 //leg->SetTextSize(0.04);
 leg->SetFillStyle(0);
 leg->SetBorderSize(0);
 leg->AddEntry(hs,"t#bar{t}","f");
 leg->AddEntry(hstop,"Single top","f");
 leg->AddEntry(hwj,"W+jets","f");
 leg->AddEntry(hzj,"Z+jets","f");
 leg->AddEntry(hqcd,"QCD","f");
 //leg->AddEntry(hvqq,"VQQ","f");
 leg->SetTextFont(42);


 TLegend *leg2 = new TLegend(0.56,0.63,0.94,0.8);  // x1,y1,x2,y2
 //leg->SetTextSize(0.04);
 leg2->SetFillStyle(0);
 leg2->SetBorderSize(0);
 leg2->AddEntry(hs2,"t#bar{t}","lep");
 //leg2->AddEntry(hstop2,"Single top","p");
 leg2->AddEntry(hwj2,"W+jets","lep");
 //leg2->AddEntry(hzj2,"Z+jets","p");
 //leg2->AddEntry(hqcd2,"QCD","l");
 //leg2->AddEntry(hvqq2,"VQQ","l");
 leg2->SetTextFont(42);

 c1.cd(1);
 leg->Draw();
 c1.cd(2);
 leg2->Draw();


 
 //c1.SaveAs("m3_stack_norm.png");
 //c1.SaveAs("m3_stack_norm.pdf");

 // make filled histo
 //c1.SaveAs("m3_stack_norm.ps"); //ok
 // gROOT->ProcessLine(".! ps2pdf m3_stack_norm.ps m3_stack_norm.pdf");

 // make shaded histo
 string out = "m3_stack_norm_2in1.eps";
 c1.SaveAs(Form("%s",out)); //ok
 gROOT->ProcessLine(Form(".!ps2pdf -dEPSCrop %s",out));
 gROOT->ProcessLine(Form(".!rm -f %s",out));
 
}
