// B1-B8
void plot_planB( string nj = "1j" ) {

  //string nj = "1j"; //0j, 1j, 2j, 3j, 4mj
  cout << "plot njet (0j,1j,2j,3j,4mj,allj) = " << nj << endl;
  //  cin >> nj;

  gStyle->SetOptStat(1110);
  gStyle->SetStatH(0.08);

  TFile f("../test_mc_mixture.root");


  TCanvas c1("c1","plan B, define control region",2*350,3*350);
  c1.Divide(2,3);

  vector<string> var;
  var.push_back("QCD_planB/QCDest_CombRelIso_AES_planB1_e20");
  var.push_back("QCD_planB/QCDest_CombRelIso_AES_planB1_e30");

  var.push_back("QCD_planB/QCDest_CombRelIso_AES_planB2_e20");
  var.push_back("QCD_planB/QCDest_CombRelIso_AES_planB2_e30");

  var.push_back("QCD_planB/QCDest_CombRelIso_AES_planB3_e20");
  var.push_back("QCD_planB/QCDest_CombRelIso_AES_planB3_e30");

  var.push_back("QCD_planB/QCDest_CombRelIso_AES_planB4_e20");
  var.push_back("QCD_planB/QCDest_CombRelIso_AES_planB4_e30");

  var.push_back("QCD_planB/QCDest_CombRelIso_AES_planB5_e20");
  var.push_back("QCD_planB/QCDest_CombRelIso_AES_planB5_e30");

  var.push_back("QCD_planB/QCDest_CombRelIso_AES_planB6_e20");
  var.push_back("QCD_planB/QCDest_CombRelIso_AES_planB6_e30");

  var.push_back("QCD_planB/QCDest_CombRelIso_AES_planB7_e20");
  var.push_back("QCD_planB/QCDest_CombRelIso_AES_planB7_e30");

  var.push_back("QCD_planB/QCDest_CombRelIso_AES_planB8_e20");
  var.push_back("QCD_planB/QCDest_CombRelIso_AES_planB8_e30");

  // plan A
  var.push_back("QCD_planA/QCDest_CombRelIso_AES_planA1_e20");
  var.push_back("QCD_planA/QCDest_CombRelIso_AES_planA1_e30");

  var.push_back("QCD_planA/QCDest_CombRelIso_AES_planA2_e20");
  var.push_back("QCD_planA/QCDest_CombRelIso_AES_planA2_e30");


  // compare control_nj to signal_nj
  TH1D *QCDinSig_nj = (TH1D*)f.Get(Form("QCD_estimation/QCDest_CombRelIso_%s__QCD",nj.c_str()));
  int rb_smallpad = 10;
  QCDinSig_nj->Rebin(rb_smallpad);
  QCDinSig_nj->GetXaxis()->SetRangeUser(0,1.5);
  QCDinSig_nj->SetStats(0);
  QCDinSig_nj->Scale(1.0/QCDinSig_nj->GetSumOfWeights());

  int npad = 1;

  for(int i=0; i<var.size(); ++i){

    //cout << "Plotting " << var[i] << endl;
    TH1D *all = (TH1D*)f.Get( Form("%s_%s__data", var[i].c_str(), nj.c_str() ) );
    TH1D *qcd = (TH1D*)f.Get( Form("%s_%s__QCD",  var[i].c_str(), nj.c_str() ) );
    TH1D *bce1 = (TH1D*)f.Get( Form("%s_%s__bce1",  var[i].c_str(), nj.c_str() ) );
    TH1D *bce2 = (TH1D*)f.Get( Form("%s_%s__bce2",  var[i].c_str(), nj.c_str() ) );
    TH1D *bce3 = (TH1D*)f.Get( Form("%s_%s__bce3",  var[i].c_str(), nj.c_str() ) );

    if(all==0) {
      cout << "\n error: histo 'all' not found, go to next one.\n"<< endl;
      continue;
    }

    TH1D *bce = bce1->Clone();
    bce->Add(bce2);
    bce->Add(bce3);


    all->GetXaxis()->SetRangeUser(0,1.5);
    qcd->GetXaxis()->SetRangeUser(0,1.5);
    bce->GetXaxis()->SetRangeUser(0,1.5);
    
    all->SetFillColor(42); //light brown
    qcd->SetFillColor(kAzure-9);
    bce->SetFillColor(kGreen+1);//kSpring-4);
    bce->SetLineColor(kGreen+1);//kSpring-4);
    bce->SetFillStyle(3452); // slopy lines
    /*
    if(var[i]=="QCD_planB/QCDest_CombRelIso_AES_planB3_e20") {
      all->SetTitle(Form("RelIso (AES B3 RT=0 ET>20) (%s,data)",nj.c_str()));
    }
    if(var[i]=="QCD_planB/QCDest_CombRelIso_AES_planB3_e30") {
      all->SetTitle(Form("RelIso (AES B3 RT=0 ET>30) (%s,data)",nj.c_str()));
    }
    */    
    //cout << "find A1_e20: " << bool(var[i].find("A1_e20")!=string::npos) << endl;
    if(npad>6 || var[i].find("A1_e20")!=string::npos) {
      //cout << "reset canvas"<< endl;
      npad=1; //reset and clear pads
      c1.Clear();
      c1.Divide(2,3);
    }
    c1.cd(npad); 
    npad++;
    all->Draw("ahist");
    qcd->Draw("ahist same");
    bce->Draw("ahist same");
    qcd->Draw("axis same");

    

    // Compare "control sample" with "QCD in signal sample (1j)"
    TPad *p1 = new TPad("p1","p1",0.5,0.35,1.0,0.89);
    p1->Draw();
    p1->cd();
    // allcopy = all = control sample
    TH1D *allcopy = all->Clone();
    allcopy->Rebin(rb_smallpad);
    allcopy->GetXaxis()->SetRangeUser(0,1.5);
    allcopy->SetStats(0);
    //    allcopy->SetLineColor(kRed);
    allcopy->SetLineWidth(2);
    allcopy->SetFillColor(41);
    //
    //  allcopy->SetMarkerStyle(20);
    double integral = allcopy->GetSumOfWeights();
    if(integral>0) 
      allcopy->Scale(1.0/integral);
    
    //    QCD_1j->Draw();
    //    allcopy->Draw("same");
    QCDinSig_nj->SetLineWidth(2);
    QCDinSig_nj->SetLineColor(kRed);    
    
    if(QCDinSig_nj->GetMaximum() > allcopy->GetMaximum()) {
      QCDinSig_nj->Draw("ah");
      allcopy->Draw("ahist same");
      allcopy->Draw("ae same");
      QCDinSig_nj->Draw("ah same");
      QCDinSig_nj->Draw("axis same");
    }else{
      allcopy->Draw("ahist");
      allcopy->Draw("ae same");
      QCDinSig_nj->Draw("ah same");
      QCDinSig_nj->Draw("axis same");
    }
    TLegend *leg = new TLegend(0.5,0.6,0.9,0.9);
    leg->AddEntry(QCDinSig_nj,"QCD in Sig.","le"); 
    leg->AddEntry(allcopy,"Control Sample","lef");
    leg->Draw();
    //    QCDinSig_nj->DrawNormalized(); //excluding under/overflows
    //    allcopy->DrawNormalized("same");

    if(i+1==6)  c1.SaveAs(Form("test_planB_%s_page1.gif",nj.c_str()));
    if(i+1==12) c1.SaveAs(Form("test_planB_%s_page2.gif",nj.c_str()));
    if(i+1==16) c1.SaveAs(Form("test_planB_%s_page3.gif",nj.c_str()));
    if(i+1==var.size()) 
      c1.SaveAs(Form("test_planA_%s.gif",nj.c_str()));

  }

  //  c1.SaveAs(Form("test3_planB_%s.gif",nj.c_str())); 
  //  c1.SaveAs(Form("test_planB_%s.gif",nj.c_str()));

}
