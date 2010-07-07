//
// 24 Mar 2010
// Macro to make plots for pull and nfit distributions
//
// outputs are named after "sample"
//-------------------------------------------------------

void plot_nfit_pull(char *sample="ref") {

  gStyle->SetTitleH(0.08);
  gStyle->SetPadTopMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetLabelSize(0.05,"xy");
  gStyle->SetTitleSize(0.06,"xy");
  gStyle->SetOptStat(1111110); //show full details (no name)
  //gStyle->SetOptStat(1110);
  //gStyle->SetStatH(0.25);
  gStyle->SetStatW(0.19);

  gROOT->ForceStyle();


  //char *sample = "wj_fastsim_nom"; //zj_nom
  //char *sample = "wj_thres20";
  //char *sample = "wj_thres5";
  //char *sample = "wj_scaleUp";
  //char *sample = "wj_scaleDown";
  //char *sample = "ttjet";
  //char *sample = "singletop";
  //char *sample = "qcd";

  //cout << "which sample: ";
  //cin >> sample;

  TFile f("FullRun_btag_new_mc_mixture.root");
  //  TFile f(Form("Norm__%s.root",sample));

  TCanvas c1("c1","m3 fit",1000,600);
  c1.Divide(2,2);

  f.cd("m3_fit");



  Color_t cp = kBlue; //pull
  Color_t cf = kRed;  //nfit


  //----------------------
  //
  // 1st plot: nfit
  //
  //-----------------------

  // increase the height of stat box for nfit
  gStyle->SetStatH(0.26);


  c1.cd(1);
  ntt_fit->Draw("e");
  ntt_fit->SetTitle("Fitted N(t#bar{t}) distribution");
  ntt_fit->GetYaxis()->SetTitle("N(toy experiment)");
  ntt_fit->GetXaxis()->SetTitle("Fitted N(event)");
  ntt_fit->Fit("gaus","Q0"); //Q=suppress message
  TF1 *f1 = ntt_fit->GetFunction("gaus");
  f1->SetLineColor(cf);
  f1->SetLineWidth(2);
  f1->Draw("same");


  c1.cd(2);
  nwj_fit->Draw("e");
  if(sample!="addZj_fitWithWjShape") 
    nwj_fit->SetTitle("Fitted N(W+jets) distribution");
  else
    nwj_fit->SetTitle("Fitted N(W/Z+jets) distribution");
  nwj_fit->SetTitle("Fitted N(W+jets) distribution");
  nwj_fit->GetYaxis()->SetTitle("N(toy experiment)");
  nwj_fit->GetXaxis()->SetTitle("Fitted N(event)");
  nwj_fit->Fit("gaus","Q0");
  TF1 *f2 = nwj_fit->GetFunction("gaus");
  f2->SetLineColor(cf);
  f2->SetLineWidth(2);
  f2->Draw("same");


  c1.cd(3);
  //nqcd_fit->Rebin(4);
  nqcd_fit->Draw("e");
  nqcd_fit->SetTitle("Fitted N(QCD) distribution");
  nqcd_fit->GetYaxis()->SetTitle("N(toy experiment)");
  nqcd_fit->GetXaxis()->SetTitle("Fitted N(event)");

  nqcd_fit->GetXaxis()->SetRangeUser(1.6,2.8);

  //nqcd_fit->GetXaxis()->SetRangeUser(22,40);
  /*
  if(sample=="QCD_norm_est")  nqcd_fit->GetXaxis()->SetRangeUser(10,30);
  if(sample=="QCD_norm_est0.5")  nqcd_fit->GetXaxis()->SetRangeUser(0,20);
  if(sample=="QCD_norm_genest")  nqcd_fit->GetXaxis()->SetRangeUser(0,35);
  if(sample=="QCD_norm_genest1.5")  nqcd_fit->GetXaxis()->SetRangeUser(0,35);
  if(sample=="QCD_norm_genest0.5")  nqcd_fit->GetXaxis()->SetRangeUser(0,35);
  */
  nqcd_fit->Fit("gaus","Q0");
  TF1 *f3 = nqcd_fit->GetFunction("gaus");
  f3->SetLineColor(cf);
  f3->SetLineWidth(2);
  f3->Draw("same");


  c1.cd(4);
  nstop_fit->Rebin(2);
  nstop_fit->Draw("e");
  nstop_fit->SetTitle("Fitted N(single top) distribution");
  nstop_fit->GetYaxis()->SetTitle("N(toy experiment)");
  nstop_fit->GetXaxis()->SetTitle("Fitted N(event)");
  //  nstop_fit->GetXaxis()->SetRangeUser(8.5,10.5);
  //  if(sample=="wj_fastsim_nom"||sample=="wj_thres20")

  nstop_fit->GetXaxis()->SetRangeUser(1.4,2.2);//10/pb

  //if(sample!="stop_plus") nstop_fit->GetXaxis()->SetRangeUser(0,22);
  //  if(sample!="stop_minus") nstop_fit->GetXaxis()->SetRangeUser(0,22);
  nstop_fit->Fit("gaus","Q0");
  TF1 *f4 = nstop_fit->GetFunction("gaus");
  f4->SetLineColor(cf);
  f4->SetLineWidth(2);
  f4->Draw("same");


  string out = Form("nfit_%s.eps",sample);

  c1.SaveAs(Form("%s",out));

  gROOT->ProcessLine(Form(".!ps2pdf -dEPSCrop %s",out));
  gROOT->ProcessLine(Form(".!rm -f %s",out));




  //--------------------------
  //
  // 2nd plot: pull
  //
  //--------------------------
  gStyle->SetStatH(0.2);
  gROOT->ForceStyle();  

  c1.cd(1);
  ntt_pull->Draw("e");
  ntt_pull->SetTitle("Pull distribution of N(t#bar{t})");
  ntt_pull->GetYaxis()->SetTitle("N(toy experiment)");
  ntt_pull->GetXaxis()->SetTitle("Pull");
  ntt_pull->Fit("gaus","Q0");
  TF1 *f1a = ntt_pull->GetFunction("gaus");
  f1a->SetLineColor(cp);
  f1a->SetLineWidth(2);
  f1a->Draw("same");


  c1.cd(2);
  nwj_pull->Draw("e");
  if(sample!="addZj_fitWithWjShape") 
    nwj_pull->SetTitle("Pull distribution of N(W+jets)");
  else
    nwj_pull->SetTitle("Pull distribution of N(W/Z+jets)");
  nwj_pull->GetYaxis()->SetTitle("N(toy experiment)");
  nwj_pull->GetXaxis()->SetTitle("Pull");
  nwj_pull->Fit("gaus","Q0");
  TF1 *f2a = nwj_pull->GetFunction("gaus");
  f2a->SetLineColor(cp);
  f2a->SetLineWidth(2);
  f2a->Draw("same");


  c1.cd(3);
  nqcd_pull->Draw("e");
  //  nqcd_pull->Rebin(5);
  nqcd_pull->SetTitle("Pull distribution of N(QCD)");
  nqcd_pull->GetYaxis()->SetTitle("N(toy experiment)");
  nqcd_pull->GetXaxis()->SetTitle("Pull");
  nqcd_pull->GetXaxis()->SetRangeUser(-0.5,0.5);
  if(sample=="QCD_norm_est0.5")     nqcd_pull->GetXaxis()->SetRangeUser(-5,5);
  if(sample=="QCD_norm_genest1.5")  nqcd_pull->GetXaxis()->SetRangeUser(-5,5);
  if(sample=="QCD_norm_genest0.5")  nqcd_pull->GetXaxis()->SetRangeUser(-5,5);
  nqcd_pull->Fit("gaus","Q0");
  TF1 *f3a = nqcd_pull->GetFunction("gaus");
  f3a->SetLineColor(cp);
  f3a->SetLineWidth(2);
  f3a->Draw("same");


  c1.cd(4);
  nstop_pull->Draw("e");
  nstop_pull->SetTitle("Pull distribution of N(single top)");
  nstop_pull->GetXaxis()->SetTitle("Pull");
  nstop_pull->GetYaxis()->SetTitle("N(toy experiment)");
  if(sample!="stop_plus" && sample!="stop_minus") 
    nstop_pull->GetXaxis()->SetRangeUser(-0.5,0.5);
  nstop_pull->Fit("gaus","Q0");
  TF1 *f4a = nstop_pull->GetFunction("gaus");
  f4a->SetLineColor(cp);
  f4a->SetLineWidth(2);
  f4a->Draw("same");


  out = Form("pull_%s.eps",sample);

  c1.SaveAs(Form("%s",out));
  c1.Close();

  gROOT->ProcessLine(Form(".!ps2pdf -dEPSCrop %s",out));
  gROOT->ProcessLine(Form(".!rm -f %s",out));


}
