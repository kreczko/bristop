//
// 13 Feb 2010
//
// scale reliso plots in signal region for "data" (all mc events)
// for 1/pb and 10/pb.
//
// to plot original histo with error reflecting MC stat error, 
// set poisson to false.
//
// to run:
//   root plot_reliso_poisson_error();
//---------------------------------------------------------------

void plot_reliso_poisson_error(){
  make_plot_reliso_poisson_error(false); //MC stat error, (n=20/pb)
  make_plot_reliso_poisson_error(true,1); //expected Poisson error (n=1/pb)
  make_plot_reliso_poisson_error(true,10);
  make_plot_reliso_poisson_error(true,20);
}

void make_plot_reliso_poisson_error(bool m_poisson=true, 
			       double m_intlumi=20){

  bool poisson = m_poisson;//false; //false = mc error

  double intlumi = m_intlumi; //pb-1

  //------------------------------
  gStyle->SetStatW(0.3);
  gStyle->SetTitleH(0.08);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleYSize(0.05);
  gStyle->SetTitleXOffset(1.0);
  gStyle->SetTitleYOffset(1.3);
  gStyle->SetPadGridY(1);
  gROOT->ForceStyle();
  //------------------------------
  Color_t Col_sig = kRed-7;
  Color_t Col_qcd = kAzure-9;
  //------------------------------


  //  TFile f("../12Feb_7TeV_planC_t1met30/test_mc_mixture.root");
  TFile f("../../test_mc_mixture.root");
  
  //Get reliso "data"
  TH1D *h1 = (TH1D*)f.Get("QCD_estimation/QCDest_CombRelIso_1j__data");
  TH1D *h2 = (TH1D*)f.Get("QCD_estimation/QCDest_CombRelIso_2j__data");
  TH1D *h3 = (TH1D*)f.Get("QCD_estimation/QCDest_CombRelIso_3j__data");
  TH1D *h4 = (TH1D*)f.Get("QCD_estimation/QCDest_CombRelIso_4mj__data");

  TH1D *s1 = (TH1D*)f.Get("QCD_estimation/QCDest_CombRelIso_1j__ttbar");
  TH1D *s2 = (TH1D*)f.Get("QCD_estimation/QCDest_CombRelIso_2j__ttbar");
  TH1D *s3 = (TH1D*)f.Get("QCD_estimation/QCDest_CombRelIso_3j__ttbar");
  TH1D *s4 = (TH1D*)f.Get("QCD_estimation/QCDest_CombRelIso_4mj__ttbar");

  TH1D *q1 = (TH1D*)f.Get("QCD_estimation/QCDest_CombRelIso_1j__QCD");
  TH1D *q2 = (TH1D*)f.Get("QCD_estimation/QCDest_CombRelIso_2j__QCD");
  TH1D *q3 = (TH1D*)f.Get("QCD_estimation/QCDest_CombRelIso_3j__QCD");
  TH1D *q4 = (TH1D*)f.Get("QCD_estimation/QCDest_CombRelIso_4mj__QCD");

  h1->Rebin(5);
  h2->Rebin(10);
  h3->Rebin(10);
  h4->Rebin(10);
  s1->Rebin(5);
  s2->Rebin(10);
  s3->Rebin(10);
  s4->Rebin(10);
  q1->Rebin(5);
  q2->Rebin(10);
  q3->Rebin(10);
  q4->Rebin(10);

  h1->GetXaxis()->SetRangeUser(0,1.59);
  h2->GetXaxis()->SetRangeUser(0,1.59);
  h3->GetXaxis()->SetRangeUser(0,1.59);
  h4->GetXaxis()->SetRangeUser(0,1.59);

  h1->SetMarkerStyle(20);
  h2->SetMarkerStyle(20);
  h3->SetMarkerStyle(20);
  h4->SetMarkerStyle(20);

  h1->SetXTitle("RelIso");
  h2->SetXTitle("RelIso");
  h3->SetXTitle("RelIso");
  h4->SetXTitle("RelIso");

  h1->SetYTitle("Entry/0.05");
  h2->SetYTitle("Entry/0.1");
  h3->SetYTitle("Entry/0.1");
  h4->SetYTitle("Entry/0.1");

  s1->SetFillColor(Col_sig);
  s2->SetFillColor(Col_sig);
  s3->SetFillColor(Col_sig);
  s4->SetFillColor(Col_sig);
  q1->SetFillColor(Col_qcd);
  q2->SetFillColor(Col_qcd);
  q3->SetFillColor(Col_qcd);
  q4->SetFillColor(Col_qcd);

  s1->SetLineColor(Col_sig);
  s2->SetLineColor(Col_sig);
  s3->SetLineColor(Col_sig);
  s4->SetLineColor(Col_sig);
  q1->SetLineColor(Col_qcd);
  q2->SetLineColor(Col_qcd);
  q3->SetLineColor(Col_qcd);
  q4->SetLineColor(Col_qcd);

  int style_sig = 3002;//3457;
  int style_qcd = 3001;//3475;
  s1->SetFillStyle(style_sig);
  s2->SetFillStyle(style_sig);
  s3->SetFillStyle(style_sig);
  s4->SetFillStyle(style_sig);
  q1->SetFillStyle(style_qcd);
  q2->SetFillStyle(style_qcd);
  q3->SetFillStyle(style_qcd);
  q4->SetFillStyle(style_qcd);


  TH1D *g1 = h1->Clone();
  TH1D *g2 = h2->Clone();
  TH1D *g3 = h3->Clone();
  TH1D *g4 = h4->Clone();

  // scale to intlumi
  //------------------
  double s = intlumi/20.0;
  g1->Scale(s);
  g2->Scale(s);
  g3->Scale(s);
  g4->Scale(s);
  s1->Scale(s);
  s2->Scale(s);
  s3->Scale(s);
  s4->Scale(s);
  q1->Scale(s);
  q2->Scale(s);
  q3->Scale(s);
  q4->Scale(s);
  

  for(int i=0; i<g1->GetNbinsX(); ++i){

    g1->SetBinError(i,sqrt(g1->GetBinContent(i))); //poisson error
    g2->SetBinError(i,sqrt(g2->GetBinContent(i))); //poisson error
    g3->SetBinError(i,sqrt(g3->GetBinContent(i))); //poisson error
    g4->SetBinError(i,sqrt(g4->GetBinContent(i))); //poisson error

  }
  g1->SetTitle(Form("1j (%.0f pb^{-1})",intlumi));
  g2->SetTitle(Form("2j (%.0f pb^{-1})",intlumi));
  g3->SetTitle(Form("3j (%.0f pb^{-1})",intlumi));
  g4->SetTitle(Form("#geq4j (%.0f pb^{-1})",intlumi));

  g1->SetMinimum(0);
  g2->SetMinimum(0);
  g3->SetMinimum(0);
  g4->SetMinimum(0);

  //980,700
  TCanvas c1("c1","reliso (poisson error)",4*300,700);
  c1.Divide(4,1,0.005,0.005);
  
  if(poisson){
    c1.cd(1);
    g1->Draw();
    s1->Draw("hist same");
    q1->Draw("hist same");
    g1->Draw("same");
    
    c1.cd(2);
    g2->Draw();
    s2->Draw("hist same");
    q2->Draw("hist same");
    g2->Draw("same");    

    c1.cd(3);
    g3->Draw();
    s3->Draw("hist same");
    q3->Draw("hist same");
    g3->Draw("same");    

    c1.cd(4);
    g4->Draw();
    s4->Draw("hist same");
    q4->Draw("hist same");
    g4->Draw("same");    

    c1.SaveAs(Form("reliso_poisson_error_%.0fpb.gif",intlumi));

  }else{

    h1->SetTitle("1j (MC stat error)");
    h2->SetTitle("2j (MC stat error)");
    h3->SetTitle("3j (MC stat error)");
    h4->SetTitle("#geq4j (MC stat error)");

    c1.cd(1);  
    h1->Draw(); 
    s1->Draw("hist same");
    q1->Draw("hist same");
    h1->Draw("same");

    c1.cd(2);  
    h2->Draw(); 
    s2->Draw("hist same"); 
    q2->Draw("hist same");
    h2->Draw("same");

    c1.cd(3);
    h3->Draw();
    s3->Draw("hist same");
    q3->Draw("hist same");
    h3->Draw("same");

    c1.cd(4);
    h4->Draw();
    s4->Draw("hist same");
    q4->Draw("hist same");
    h4->Draw("same");

    c1.SaveAs(Form("reliso_mc_stat_error_%.0fpb.gif",intlumi));
  }

}
