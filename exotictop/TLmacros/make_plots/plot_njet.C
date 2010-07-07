//
// 25 Mar 2010 
// macro to make njet plots
// adapted from Frankie's macro xSecCalc.cc
//  2 in 1: linear, log
// 
// adapt for "data" mode.
// when we have data, specify the name of data file in FileNameData. 
//-------------------------------------------
bool isData = true;
bool debug = false; 
int intlumi = 100; //in pb-1

void plot_njet() {

  // set to true if want solid fill colour
  bool wantSolidColour = true;

  double mytopmargin = 0.12;

  gStyle->SetLineScalePS(2); //D=3
  gStyle->SetLabelSize(0.06,"y");
  gStyle->SetPadTopMargin(mytopmargin);
  gStyle->SetMarkerStyle(20); //big data point
  //  gStyle->SetTitleW(0.); //new


  // all samples
  //  TString FileName = "~/cms/analysis/results/2009_Oct/8oct_octx_SD/test_mc_mixture.root";
  //  TString FileName = "../test_mc_mixture.root";

  TString FileNameMC   = "../test_mc_mixture.root"; //mc
  TString FileNameData = "../test_mc_mixture.root"; //data


  TFile *RootFileMC;
  TFile *RootFileData;
  //TString SampleName[5] = {"mixture","ttbar","wjet","zjet","qcd"};
  RootFileMC   = new TFile(Form("%s",(char*)FileNameMC));
  RootFileData = new TFile(Form("%s",(char*)FileNameData));


  // Read intlumi from file (MC)
  TH1D *info = (TH1D*)RootFileMC->Get("info");
  if(info>0) {
    intlumi = info->GetBinContent(3);
    cout << "reading in int lumi of " << intlumi << "/pb"<< endl;
  }

  // title
  string myTitle = Form("e+jets (7 TeV, %d pb^{-1})",intlumi);
  cout << myTitle << endl;
  //string myTitle ="CMS Preliminary (20 pb^{-1})";



  //-------------------------------------
  Color_t Col_sig  = kRed-7;   //kBlue-7; +1
  Color_t Col_qcd  = kAzure-9; //kGreen-6
  Color_t Col_wj   = kGreen-9; //kAzure+9;
  Color_t Col_zj   = kYellow-7;
  Color_t Col_vqq  = kGray+1;
  Color_t Col_stop = kViolet-4;
  //-------------------------------------





  //create a stack for events from each sample as a function of njets
  njetEventPass = new THStack("njetEventPass","Events per sample passing event selection");
  

  //Either access a 1-D histogram that is already an njet distribution,
  //or if we make a multiple-dimensional histogram based on the e_plus_jet[i][j][k]
  //we can extract the events passing


  const float Jetend = 4;//5 for 1->=4 jets, 6 to include total for all jets
  int Jetenda = Jetend;
  //declare histograms for stacked plot
  TH1D *DataEventsPass       = new TH1D("njetsData",      "",Jetenda,0,Jetenda);
  TH1D *SignalEventsPass     = new TH1D("njetsSignal",    "",Jetenda,0,Jetenda);
  TH1D *QCDEventsPass        = new TH1D("njetsQCD",       "",Jetenda,0,Jetenda);
  TH1D *WjetsEventsPass      = new TH1D("njetsWjets",     "",Jetenda,0,Jetenda);
  TH1D *ZjetsEventsPass      = new TH1D("njetsZjets",     "",Jetenda,0,Jetenda);
  //TH1D *VQQEventsPass      = new TH1D("njetsVQQ",       "",Jetenda,0,Jetenda);
  TH1D *SingleTopEventsPass  = new TH1D("njetsSingleTop", "",Jetenda,0,Jetenda);

  // Get 2D histogram from root file
  //----------------------------------
  SignalEvents = (TH2D*)RootFileMC->Get(Form("Signal_njetsVcuts"));
  QCDEvents    = (TH2D*)RootFileMC->Get(Form("QCD_njetsVcuts"));
  WjetsEvents  = (TH2D*)RootFileMC->Get(Form("Wjets_njetsVcuts"));
  ZjetsEvents  = (TH2D*)RootFileMC->Get(Form("Zjets_njetsVcuts"));
  //VQQEvents    = (TH2D*)RootFileMC->Get(Form("VQQ_njetsVcuts"));
  SingleTopEvents = (TH2D*)RootFileMC->Get(Form("SingleTop_njetsVcuts"));

  TH2D* DataEvents(0);
  if(isData)  DataEvents = (TH2D*)RootFileData->Get("Data_njetsVcuts"); 
  if(isData && debug) DataEvents->ls();

  //  int endCut = 5; //5th from below, ie Stage 10: HT -> endCut=5
  int endCut = 3; //3th from below, ie Stage 11 (not counting 4j): BARREL -> endCut=3
  int startbin = 0;
  float Neventsin4j = 0;

  for(int j=1;j<=Jetend;j++){

    //second index is for cut/stage.
    float binvalue_sig       =     SignalEvents->GetBinContent(j+1,endCut);
    float binvalue_qcd       =        QCDEvents->GetBinContent(j+1,endCut);
    float binvalue_wj        =      WjetsEvents->GetBinContent(j+1,endCut);
    float binvalue_zj        =      ZjetsEvents->GetBinContent(j+1,endCut);
    //float binvalue_vqq     =       VQQEvents->GetBinContent(j+1,endCut);
    float binvalue_singletop = SingleTopEvents->GetBinContent(j+1,endCut);
    float binvalue_data;
    float binerror_data;
    if(isData) binvalue_data =      DataEvents->GetBinContent(j+1,endCut);
    if(isData) binerror_data =      DataEvents->GetBinError(  j+1,endCut); //=sqrt(binvalue)

    SignalEventsPass->SetBinContent(j+startbin,binvalue_sig);
    QCDEventsPass   ->SetBinContent(j+startbin,binvalue_qcd);
    WjetsEventsPass ->SetBinContent(j+startbin,binvalue_wj);
    ZjetsEventsPass ->SetBinContent(j+startbin,binvalue_zj);
    //VQQEventsPass   ->SetBinContent(j+startbin,binvalue_vqq);
    SingleTopEventsPass->SetBinContent(j+startbin,binvalue_singletop);

    if(isData) {
      DataEventsPass->SetBinContent(j+startbin, binvalue_data);
      DataEventsPass->SetBinError(  j+startbin, binerror_data );
    }

  }//end nj loop: 1,2,3,>=4j
  


  SignalEventsPass->SetFillColor(Col_sig);
  QCDEventsPass   ->SetFillColor(Col_qcd);
  WjetsEventsPass ->SetFillColor(Col_wj);
  ZjetsEventsPass ->SetFillColor(Col_zj);
  SingleTopEventsPass->SetFillColor(Col_stop);



  if (!wantSolidColour) {
    gStyle->SetHatchesSpacing(1);
    gStyle->SetHatchesLineWidth(1);

    SignalEventsPass->SetFillStyle(3001);
    //SingleTopEventsPass->SetFillStyle(3002);

    QCDEventsPass->SetFillStyle(3344);
    WjetsEventsPass->SetFillStyle(3354);
    ZjetsEventsPass->SetFillStyle(3345);
  }
  

  // legend
  //---------
  TLegend *leg1 = new TLegend(0.6,0.6,1,0.88); 

  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->AddEntry(DataEventsPass,      " Data",           "LPE");
  leg1->AddEntry(SignalEventsPass,    " Signal t#bar{t}","f");
  leg1->AddEntry(SingleTopEventsPass, " Single top",     "f");
  leg1->AddEntry(WjetsEventsPass,     " W+jets",         "f");
  leg1->AddEntry(ZjetsEventsPass,     " Z+jets",         "f");
  leg1->AddEntry(QCDEventsPass,       " QCD",            "f"); 
  //  leg1->SetTextSize(0.06);


  // add to stacked histogram
  //-------------------------
  //njetEventPass->Add(VQQEventsPass); //bottom
  njetEventPass->Add(QCDEventsPass);
  njetEventPass->Add(ZjetsEventsPass);
  njetEventPass->Add(WjetsEventsPass);
  njetEventPass->Add(SingleTopEventsPass);
  njetEventPass->Add(SignalEventsPass); //top


  //---------
  // x axis
  //---------
  string jetlabels[5] = {"0-jet","1-jet","2-jets","3-jets","#geq4-jets"};

  TAxis *xaxis  = SignalEventsPass ->GetXaxis();
  
  for(int i=1;i<5;i++){
    xaxis->SetBinLabel(i,Form("%s",jetlabels[i].c_str())); //starts at 1j
  }


  //-----------------
  //  Finally, draw
  //-----------------

  //  TCanvas *c1 = new TCanvas("c1","c1",10,10,800,410);
  TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,520);
  c1->Divide(2,1);



  //----------------
  // (1) Plot 1-4j
  //----------------
  c1->cd(1);
  gPad->SetLeftMargin(0.17);
  gPad->SetRightMargin(0.03);
  njetEventPass->SetTitle( myTitle.c_str() );
  njetEventPass->Draw();
  njetEventPass->GetXaxis()->SetLabelSize(0.1);
  if(isData) {
    cout << "Drawing data points on pad 1."<< endl;
    //    DataEventsPass->Draw("same");
    DataEventsPass->SetMarkerSize(1);
    DataEventsPass->DrawCopy("e0 same");
  }


  for(int i=1;i<5;i++){
    njetEventPass->GetXaxis()->SetBinLabel(i,Form("%s",jetlabels[i].c_str())); //starts at 1j
  }
  //    njetEventPass->Draw("axis same");
  
  //njetEventPass->GetYaxis()->SetTicks("+-");
  
  // vertical axis
  // "U": no label
  float yaxis_max = gPad->GetUymax(); 
  TGaxis *yaxis = new TGaxis(gPad->GetUxmin(),gPad->GetUymin(),
			     gPad->GetUxmin(),gPad->GetUymax(),
			     0,yaxis_max,510,"-RU");
  yaxis->SetLabelFont(42);
  yaxis->SetLabelSize(0.06);
  yaxis->Draw();
  
  // x axis (no label)
  //    float xaxis_max = gPad->GetUxmax(); 
  TGaxis *x2 = new TGaxis(0,gPad->GetUymin(),4,gPad->GetUymin(),
			  0,4,400,"+RU");
  x2->Draw();

  // horizontal axis
  //  TLine l(0,0,4,0);
  //  l.Draw(); 
  
  //TLine lb(4,0,4,2000);
  //lb.Draw();
  
  leg1->Draw();



  //--------------------------------
  // (2) Plot 1-4j (log scale in y)
  //--------------------------------
  cout << "Drawing log scale plot on pad 2."<< endl;

  c1->cd(2);

  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.1);
  gPad->SetLogy(1);
    
  THStack *njetlog = (THStack*)njetEventPass->Clone("njetlog");
   
  njetlog->SetTitle( myTitle.c_str() );
  njetlog->GetYaxis()->SetTicks("-");//+-
  njetlog->Draw("hist");
    
  njetlog->SetMinimum(20); //30, 50
  njetlog->SetMaximum(2e4); //30, 50
  if(intlumi==20) { //10/pb
    njetlog->SetMinimum(4.0);
    njetlog->SetMaximum(1e4);
  }
  else if(intlumi==10) { //10/pb
    njetlog->SetMinimum(2.0);
    njetlog->SetMaximum(1e4);
  }

  if(isData) {      
    DataEventsPass->DrawCopy("e0 same");
  }

  //    njetlog->GetYaxis()->SetRangeUser(2,2e4);
  //njetlog->SetMaximum(1e5);
  leg1->Draw();
        
  // horizontal line
  //l.Draw();
  // vertical line
  //  float yup = 12;
  //  float ylow = 1.6e4;
  float yup = 1.6e4;
  float ylow = 12;
  if(intlumi==10) yup = 5e3;
  TLine *l3  = new TLine(0,ylow,0,yup);  //left vertical line   
  TLine *l3b = new TLine(4,ylow,4,yup); //right vertical line
  //  l3->Draw();
  //  l3b->Draw();

  TLine *line = new TLine(0,0,1,1);
  line->DrawLineNDC(0.9, 0.1, 0.9, 1-mytopmargin); //right vertical
  //line->DrawLineNDC(0.1, 0.1, 0.1, 1-mytopmargin); //left vertical



  /*
  // Draw horizontal axis for log plot
  TF1 *f1=new TF1("f1","+x",0,4);
  float yymin = gPad->GetUymin();
  TGaxis *A1 = new TGaxis(0,yymin,4,yymin,"f1",400,"+");
  A1->SetTitle("");
  A1->SetLabelSize(0);
  A1->Draw();
  */

  // x axis (no label)
  //    float xaxis_max = gPad->GetUxmax(); 
  //TGaxis *x2 = new TGaxis(0,gPad->GetUymin(),4,gPad->GetUymin(),
  //			    0,4,400,"+RU");
  x2->Draw();
  


  //-----------------------------------
  // check if the numbers are correct
  //-----------------------------------

  //printf("%3s (%10s) %10s %10s %10s %10s %10s %10s %12s\n","bin","label",
  //	 "n(tt)","n(wj)","n(zj)","n(qcd)","n(vqq)","n(singletop)","Total");

  printf("%3s (%10s) %10s %10s %10s %10s %10s %12s\n","bin","label",
	 "n(tt)","n(wj)","n(zj)","n(qcd)","n(singletop)","Total");

  float nTotal = 0;

  for(int i=1; i<=Jetenda; i++){

    nTotal = 0;
    nTotal += SignalEventsPass->GetBinContent(i);
    nTotal += QCDEventsPass->GetBinContent(i);
    nTotal += WjetsEventsPass->GetBinContent(i);
    nTotal += ZjetsEventsPass->GetBinContent(i);
    //nTotal += VQQEventsPass->GetBinContent(i);
    nTotal += SingleTopEventsPass->GetBinContent(i);

    printf("%3d (%10s) %10.2f %10.2f %10.2f %10.2f %10.2f   %12.2f\n",i,
	   SignalEventsPass->GetXaxis()->GetBinLabel(i),
	   SignalEventsPass->GetBinContent(i),
	   WjetsEventsPass->GetBinContent(i),
	   ZjetsEventsPass->GetBinContent(i),
	   QCDEventsPass->GetBinContent(i),
	   //VQQEventsPass->GetBinContent(i),
	   SingleTopEventsPass->GetBinContent(i),
	   nTotal);
  }


  //---------------
  //  Save plots
  //---------------  
  c1->SaveAs("njet_2in1.eps");
  c1->cd(2);
  gPad->SaveAs("njet_log.eps");
  c1->Close();

  gROOT->ProcessLine(".!ps2pdf -dEPSCrop njet_2in1.eps");
  gROOT->ProcessLine(".!ps2pdf -dEPSCrop njet_log.eps");
  gROOT->ProcessLine(".!rm njet_2in1.eps");
  gROOT->ProcessLine(".!rm njet_log.eps");

}
