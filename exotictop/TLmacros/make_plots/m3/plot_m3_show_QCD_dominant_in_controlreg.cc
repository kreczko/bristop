{

  gStyle->SetOptStat(0);
  gStyle->SetTitleH(0.07);
  gStyle->SetTitleW(0.8);


  //  TFile f("~/cms/analysis/results/18May09_FB_new_ttjet/run_mixture_18May.root");
  TFile f("../../test_mc_mixture.root");


  TCanvas c1("c1","m3 (summer08)",1200,400); 
  c1.Divide(2,1);


  TH1D *h_m3_control   = (TH1D*)f.Get("Wjets_estimation/hadTop_maxPT_mass_nonIso_4j");
  TH1D *h_m3_qcd_control = (TH1D*)f.Get("Wjets_estimation/m3_control__QCD");
  

  Color_t Col_qcd = kAzure-9;
  Color_t Col_qcd2 = kAzure+2;

  const int rb = 80;
  h_m3_control->Rebin(rb);
  h_m3_qcd_control->Rebin(rb);


  // first plot (to scale)
  c1.cd(1);

  //  h_m3_control->SetLineWidth(3);
  h_m3_qcd_control->SetFillColor(Col_qcd);

  h_m3_control->SetTitle("M3 distribution in control region (non-isol.)");
  h_m3_control->GetXaxis()->SetTitle("M3 (GeV)");
  h_m3_control->GetYaxis()->SetTitle("N(events)");

  h_m3_control->Draw("ahist");
  h_m3_qcd_control->Draw("ahist same");
  h_m3_control->Draw("axis same");


  // normalized-to-1 plot
  c1.cd(2);
  TH1D *h_m3_control2   = h_m3_control->Clone();
  TH1D *h_m3_qcd_control2 = h_m3_qcd_control->Clone();

  h_m3_control2->SetTitle("Normalized M3 distribution in control region (non-isol.)");

  h_m3_control2->SetMarkerStyle(4);  //open circle
  //  h_m3_control2->SetMarkerSize(1.6);  //bullet
  //  h_m3_control2->DrawNormalized("");

  h_m3_qcd_control2->SetFillColor(0);
  h_m3_qcd_control2->SetLineColor(kBlue);
  h_m3_qcd_control2->SetMarkerColor(kBlue);
  h_m3_qcd_control2->SetLineStyle(kDashed);
  h_m3_qcd_control2->SetMarkerStyle(22);
  //  h_m3_qcd_control2->SetMarkerSize(1.6);

  //  h_m3_control2->SetLineWidth(2);
  //  h_m3_qcd_control2->SetLineWidth(2);

  // norm to 1
  h_m3_control2->Scale(1./h_m3_control2->Integral());
  h_m3_qcd_control2->Scale(1./h_m3_qcd_control2->Integral());

  h_m3_control2->SetMaximum(0.38);
  h_m3_qcd_control2->SetMaximum(0.38);

  h_m3_control2->Draw("");
  h_m3_qcd_control2->Draw("same");


  // Add legend
  //-------------
  c1.cd(1);
  TLegend *leg = new TLegend(0.5,0.6,0.88,0.85);  // x1,y1,x2,y2
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  //  leg->AddEntry(h_m3_control,    "t#bar{t}+all backgrounds","l");
  leg->AddEntry(h_m3_control,    "data control region","l");
  leg->AddEntry(h_m3_qcd_control,"QCD","f");
  leg->SetTextFont(42);
  leg->Draw();

  c1.cd(2);
  TLegend *leg2 = new TLegend(0.5,0.6,0.88,0.85);  // x1,y1,x2,y2
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  //  leg2->AddEntry(h_m3_control2,    "t#bar{t}+all backgrounds","lp");
  leg2->AddEntry(h_m3_control2,    "data control region","lp");
  leg2->AddEntry(h_m3_qcd_control2,"QCD","lp");
  leg2->SetTextFont(42);
  leg2->Draw();

  
  char *outname ="m3_show_QCD_dominant_in_controlreg";
  //  c1.SaveAs( Form("%s.png", outname ));
  c1.SaveAs( Form("%s.eps", outname ));

  gROOT->ProcessLine(Form(".!ps2pdf -dEPSCrop %s.eps",outname));
  gROOT->ProcessLine(Form(".!rm %s.eps",outname));

}
