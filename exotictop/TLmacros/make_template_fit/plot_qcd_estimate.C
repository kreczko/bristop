//
// 17-Feb 2010
// Plots (signed) Delta = est-true, for QCD estimate
//-----------------------------------------------------

const int n = 4;
const int nrange = 9;

bool plot_qcd_estimate(string func="template",
		       string control_sample ="B3_e20",
		       int ntemp = 1,
		       string tempbin ="2j"
		       ){

  if(ntemp==1) 
    cout << ">> using 1 template from "<< tempbin << endl;
  else 
    cout << ">> using 4 templates" << endl;


  double x[4] = {1,2,3,4};
  //  double xx[2] = {3,4};
  /*
  string func = "gaus";
  cout << "enter function: ";
  cin >> func;
  //   string sfunc = "pol3";
  
  */
  // Gaus (Free-fits)
  ifstream result;
  if(ntemp==4){
    result.open(Form("est_%s_%s_4T.txt",
		     func.c_str(),
		     control_sample.c_str()) );
  }else if (ntemp==1) {
    result.open(Form("est_%s_%s_1T_%s.txt",
		     func.c_str(),
		     control_sample.c_str(),
		     tempbin.c_str()) );
  }
  else{
    cout << "stop, result file not found." << endl;
    return;
  }

  int m = nrange*(4);
  double est[500];
  double max_dev = 0;
  double max_dev_1j = 0;
  double max_dev_2j = 0;
  double max_dev_3j = 0;
  double max_dev_4j = 0;

  for(int i=0; i<m; ++i){

    result >> est[i];

    if( fabs(est[i]) > max_dev ) max_dev = fabs(est[i]);


    if( (i+1)%4 == 1 ) { //1j
      max_dev_1j = TMath::Max( max_dev_1j, fabs(est[i]) );
    }
    else if( (i+1)%4 == 2 ) { //2j
      max_dev_2j = TMath::Max( max_dev_2j, fabs(est[i]) );
    }
    else if( (i+1)%4 == 3 ) { //3j
      max_dev_3j = TMath::Max( max_dev_3j, fabs(est[i]) );
    }
    else if( (i+1)%4 == 4 ) { //4j
      max_dev_4j = TMath::Max( max_dev_4j, fabs(est[i]) );
    }
    
  }//read result

  result.close();
  max_dev = TMath::Max( max_dev_1j, max_dev_2j);
  max_dev = TMath::Max( max_dev, max_dev_3j);
  max_dev = TMath::Max( max_dev, max_dev_4j);

  cout << "\n For 1 j, |max| dev is " << max_dev_1j ;
  cout << "\n For 2 j, |max| dev is " << max_dev_2j ;
  cout << "\n For 3 j, |max| dev is " << max_dev_3j ;
  cout << "\n For 4 j, |max| dev is " << max_dev_4j ;
  cout << "\n For allj, all ranges, |max| dev is " << max_dev << endl;


  TCanvas c1("c1","QCD estimates",600,600);

  double y[nrange][4];

  int index=0;

  for(int i=0; i<nrange; ++i){
    for(int j=0; j<4; ++j){ //4 values
      y[i][j] = est[index]; //read in
      index++;
      //cout << "index="<<index<< endl;
      //cout << y[i][j]<<endl;;
    }
  }

  double yy[nrange][4];

  for(int i=0; i<nrange; ++i){
    for(int j=0; j<2; ++j){ //2 values
      yy[i][j] = est[index]; //read in
      index++;
      //cout << "index="<<index<< endl;
    }
  }



  gStyle->SetMarkerSize(1.7);
  gStyle->SetMarkerStyle(20);
  c1.SetTopMargin(0.1);
  c1.SetLeftMargin(0.12);
  c1.SetRightMargin(0.35);

  TGraph *gr1 = new TGraph(n,x,y[1-1]);
  TGraph *gr2 = new TGraph(n,x,y[2-1]);
  TGraph *gr3 = new TGraph(n,x,y[3-1]);
  TGraph *gr4 = new TGraph(n,x,y[4-1]);
  TGraph *gr5 = new TGraph(n,x,y[5-1]);
  TGraph *gr6 = new TGraph(n,x,y[6-1]);
  TGraph *gr7 = new TGraph(n,x,y[7-1]);
  TGraph *gr8 = new TGraph(n,x,y[8-1]);
  TGraph *gr9 = new TGraph(n,x,y[9-1]);

  TGraph *gr_1 = gr1;
  TGraph *gr_2 = gr2;
  TGraph *gr_3 = gr3;
  TGraph *gr_4 = gr4;
  TGraph *gr_5 = gr5;
  TGraph *gr_6 = gr6;
  TGraph *gr_7 = gr7;
  TGraph *gr_8 = gr8;
  TGraph *gr_9 = gr9;

  gr1->SetMarkerColor(kGreen+1);
  gr2->SetMarkerColor(kGreen+2);
  gr3->SetMarkerColor(kGreen+3);
  gr4->SetMarkerColor(kAzure+7);
  gr5->SetMarkerColor(kAzure-3);
  gr6->SetMarkerColor(kBlue);
  gr7->SetMarkerColor(kOrange);
  gr8->SetMarkerColor(kOrange-1);
  gr9->SetMarkerColor(kOrange-6);

  gr_1->SetMarkerColor(kRed);
  gr_2->SetMarkerColor(kRed+1);
  gr_3->SetMarkerColor(kRed+2);
  gr_4->SetMarkerColor(kViolet-4);
  gr_5->SetMarkerColor(kViolet);
  gr_6->SetMarkerColor(kViolet-1);
  gr_7->SetMarkerColor(kGray+1);
  gr_8->SetMarkerColor(kGray+2);
  gr_9->SetMarkerColor(kGray+3);

  gr_1->SetMarkerStyle(22);
  gr_2->SetMarkerStyle(22);
  gr_3->SetMarkerStyle(22);
  gr_4->SetMarkerStyle(22);
  gr_5->SetMarkerStyle(22);
  gr_6->SetMarkerStyle(22);
  gr_7->SetMarkerStyle(22);
  gr_8->SetMarkerStyle(22);
  gr_9->SetMarkerStyle(22);




  // To get desired x range, draw blank histo
  gStyle->SetTitleW(0.9);
  if(func=="gaus")
    h = new TH1F("h","Variation of QCD estimates with fit range (Gaussian)",4,0.5,4.5);
  if(func=="pol3")
    h = new TH1F("h","Variation of QCD estimates with fit range (Pol3)",4,0.5,4.5);
  if(func=="pol4")
    h = new TH1F("h","Variation of QCD estimates with fit range (Pol4)",4,0.5,4.5);
  if(func=="template")
    h = new TH1F("h","Variation of QCD estimates with range (template fit)",4,0.5,4.5);



  h->SetStats(kFALSE); // no statistics
  h->Draw();
  h->SetYTitle("Deviation = (Est-True)/True");

  double show_range = int(max_dev)+1;

  h->GetYaxis()->SetRangeUser( 0-show_range, show_range );
  h->GetXaxis()->SetRangeUser(0.5,5.5);
  h->GetXaxis()->SetBinLabel(1.,"1j");
  h->GetXaxis()->SetBinLabel(2.,"2j");
  h->GetXaxis()->SetBinLabel(3.,"3j");
  h->GetXaxis()->SetBinLabel(4.,"#geq4j");
  h->GetXaxis()->SetLabelSize(0.07);
  h->GetYaxis()->SetTitleOffset(1.3);


  // Constrained fits
  gr_1->Draw("P");
  gr_2->Draw("P");
  gr_3->Draw("P");
  gr_4->Draw("P");
  gr_5->Draw("P");
  gr_6->Draw("P");
  gr_7->Draw("P");
  gr_8->Draw("P");
  gr_9->Draw("P");


  c1.SetGrid(1,1);



  TLegend leg(0.65,0.1,0.98,0.9);
  leg.SetHeader("Fit Range");
  leg.SetFillColor(0);

  leg.AddEntry(gr_1,"0.2-1.0","p");
  leg.AddEntry(gr_2,"0.2-1.2","p");
  leg.AddEntry(gr_3,"0.2-1.4","p");
  leg.AddEntry(gr_4,"0.3-1.1","p");
  leg.AddEntry(gr_5,"0.3-1.3","p");
  leg.AddEntry(gr_6,"0.3-1.5","p");
  leg.AddEntry(gr_7,"0.4-1.2","p");
  leg.AddEntry(gr_8,"0.4-1.4","p");
  leg.AddEntry(gr_9,"0.4-1.6","p");

  leg.Draw();

  
  //-----------
  // Save plot
  //-----------
  string out = Form("%s_%s_%dT", func.c_str(),  control_sample.c_str(),   ntemp);
  if(ntemp==1) out = out + "_" + tempbin; // + "_1j"
    
  cout << "save as: " << out << endl;

  c1.SaveAs( Form("qcd_estimate_%s.gif", out.c_str()) );

  // zoom in/out
  h->GetYaxis()->SetRangeUser( -1, 1 );
  c1.SaveAs( Form("qcd_estimate_%s_zoom_11.gif", out.c_str()) );

  h->GetYaxis()->SetRangeUser( -2, 2 );
  c1.SaveAs( Form("qcd_estimate_%s_zoom_22.gif", out.c_str()) );

  h->GetYaxis()->SetRangeUser( -4, 4 );
  c1.SaveAs( Form("qcd_estimate_%s_zoom_44.gif", out.c_str()) );

  h->GetYaxis()->SetRangeUser( -6, 6 );
  c1.SaveAs( Form("qcd_estimate_%s_zoom_66.gif", out.c_str()) );

  h->GetYaxis()->SetRangeUser( -8, 8 );
  c1.SaveAs( Form("qcd_estimate_%s_zoom_88.gif", out.c_str()) );


  gDirectory->DeleteAll();
  return;
}
