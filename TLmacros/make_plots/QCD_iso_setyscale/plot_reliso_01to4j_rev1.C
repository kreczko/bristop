//
// 3 Mar 2010: reliso plot, 4 in 1 (separate pdf)
//
// op: NES, old AES (e30), new AES (planB3_e30)
//
// sel = QCD, data, dijet
//
// (1) plot_reliso_QCD_01to4j()
//     1j 2j 3j 4mj (NES)  reliso_QCD_1to4j.pdf
//     1j 2j 3j 4mj (AES)  reliso_QCD_1to4j_AES.pdf
//     1j 2j 3j 4mj (AES)  reliso_QCD_1to4j_AES_planB3.pdf
//
// (2) plot_reliso_data_01to4j()
//     1j 2j 3j 4mj (NES)  reliso_data_1to4j.pdf
//     1j 2j 3j 4mj (AES)  reliso_data_1to4j_AES.pdf
//     1j 2j 3j 4mj (AES)  reliso_data_1to4j_AES_planB3.pdf
//
// (3) plot_reliso_dijet_01to4j(15)
//     1j 2j 3j 4mj (NES)  reliso_dijet_pt( )_1to4j.pdf
//     1j 2j 3j 4mj (AES)  reliso_dijet_pt( )_1to4j_AES.pdf
//     1j 2j 3j 4mj (AES)  reliso_dijet_pt( )_1to4j_AES_planB3.pdf
//
//----------------------------------------------------------------
const bool plot_also_0j = false;
const bool setyscale = true;

string old_AES = "AES_e30";//old AES ele et cut
string planB   = "planB3_e20";
string plotTitle;
int myTitleW = 0.5;


//-----------------------
//  main
//-----------------------
void plot_reliso_01to4j_rev1(string sel="QCD"){
  if     (sel=="QCD")   plot_reliso_QCD_01to4j();
  else if(sel=="data")  plot_reliso_data_01to4j();
  else if(sel=="dijet") plot_reliso_dijet_01to4j();
}

//-----------------------
//  Total QCD
//-----------------------
void plot_reliso_QCD_01to4j() {
  plot_reliso_01to4j("QCD");
}

//-----------------------
//  data
//-----------------------
void plot_reliso_data_01to4j() {
  plotTitle = "Data";
  plot_reliso_01to4j("data");
}

//----------------------
//  pythia dijet
//----------------------
void plot_reliso_dijet_01to4j(){
  plot_reliso_dijet_01to4j(15);
  plot_reliso_dijet_01to4j(30);
  plot_reliso_dijet_01to4j(80);
  plot_reliso_dijet_01to4j(170);
  plot_reliso_dijet_01to4j(300);
  plot_reliso_dijet_01to4j(470);
  plot_reliso_dijet_01to4j(800);
  plot_reliso_dijet_01to4j(1400);
}

int dijet_pt;
void plot_reliso_dijet_01to4j(int dijet_pt_user=15) {
  dijet_pt = dijet_pt_user;
  cout << "dijet pt: " << dijet_pt << endl;
  plotTitle = "QCD Dijets";
  myTitleW = 0.8;
  plot_reliso_01to4j("dijet");
}

//----------------
//  actual code
//----------------
bool plot_reliso_01to4j(string sel="QCD") {

  if(sel!="QCD" && sel!="data" && sel!="dijet") {
    cout << "wrong input, sel must be QCD, data, or dijet" << endl;
    return 0;
  }

  int  this_rebin = 10;

  //  int rebin[8] = {1,2,2,2,
  //		  1,2,2,2};


  int nplot = 4; //1j-4j
  if(plot_also_0j) nplot = 5;



  //  int run_eleET_cut = input;
  //  cout <<  "enter ele et cut = ";
  //  cin >> run_eleET_cut;



  gStyle->SetOptStat(0);
  gStyle->SetTitleH(0.08);
  gStyle->SetTitleW(myTitleW);//@@@
  gStyle->SetLabelSize(0.05,"xy");
  gStyle->SetTitleSize(0.05); //axis title
  gStyle->SetLineScalePS(2);
  gStyle->SetPadTopMargin(0.15) ;///new
  gStyle->SetPadLeftMargin(0.12) ;///new
  gROOT->ForceStyle();
  TGaxis::SetMaxDigits(4);//D=5



  TCanvas c1("c1","reliso in QCD", nplot*300,300);
  c1.Divide(nplot,1,0.002,0.002);

  //  TFile f(Form("~/cms/analysis/results/2009_Oct/9oct_all/e%d/test_mc_mixture.root",run_eleET_cut));
  string file = "../test_mc_mixture.root";
  if(setyscale) file = "../../test_mc_mixture.root";

  TFile f;
  //  TFile f("../25Feb_7TeV_plainQCD_pt15/test_data.root");
  //  TFile f("../25Feb_7TeV_plainQCD_pt30/test_data.root");
  //  TFile f("../26Feb_7TeV_plainQCD_pt80/test_data.root");
  if(sel=="dijet"){
    if     (dijet_pt==15)   file = "../25Feb_7TeV_plainQCD_pt15/test_data.root";
    else if(dijet_pt==30)   file = "../25Feb_7TeV_plainQCD_pt30/test_data.root";
    else if(dijet_pt==80)   file = "../26Feb_7TeV_plainQCD_pt80/test_data.root";
    else if(dijet_pt==170)  file = "../26Feb_7TeV_plainQCD_pt170/test_data.root";
    else if(dijet_pt==300)  file = "../26Feb_7TeV_plainQCD_pt300/test_data.root";
    else if(dijet_pt==470)  file = "../26Feb_7TeV_plainQCD_pt470/test_data.root";
    else if(dijet_pt==800)  file = "../26Feb_7TeV_plainQCD_pt800/test_data.root";
    else if(dijet_pt==1400) file = "../26Feb_7TeV_plainQCD_pt1400/test_data.root";
  }
  f.Open(file.c_str());
  //  cout << "file: "<< f.GetName() << endl;
  cout << "file: "<< file << endl;



  
  //-------------------------  
  Color_t Col_sig = kRed-7;//kBlue-7; +1
  Color_t Col_qcd = kAzure-9;//kGreen-6
  Color_t Col_wj  = kGreen-9;//kAzure+9;
  Color_t Col_zj  = kYellow;
  // darker
  Color_t Col_qcd2 = kAzure+2;
  Color_t Col_wj2  = kGreen+2;
  //-------------------------


  const char *jetlabel[15] = {"0-jet","1-jet", "2-jets", "3-jets","#geq4-jets",
			      "0-jet","1-jet", "2-jets", "3-jets","#geq4-jets",
			      "0-jet","1-jet", "2-jets", "3-jets","#geq4-jets"};

  int npad=1;

  for (int i=1; i<=15; ++i ) {

    if(nplot==4 && (i==1||i==6||i==11) ) continue;//skip 0j


    char *var = "QCD_estimation/QCDest_CombRelIso";
    if(i>5) var = "QCD_estimation/AES/QCDest_CombRelIso_AES";
    if(i>10) var = Form("QCD_planB/QCDest_CombRelIso_AES_%s", planB.c_str());
 

    string nj = "0j";
    if      (i==2||i== 7||i==12)  nj = "1j";
    else if (i==3||i== 8||i==13)  nj = "2j";
    else if (i==4||i== 9||i==14)  nj = "3j";
    else if (i==5||i==10||i==15)  nj = "4mj";
    


    // All events/data
    TH1D *data = (TH1D*)gDirectory->Get( Form("%s_%s__data",var,nj.c_str()) );

    // All QCD
    TH1D *qcd = (TH1D*)gDirectory->Get( Form("%s_%s__QCD",var,nj.c_str()) );
    //qcd->ls();

    // All bce 
    TH1D *bce1 = (TH1D*)gDirectory->Get( Form("%s_%s__bce1",var,nj.c_str()) );
    TH1D *bce2 = (TH1D*)gDirectory->Get( Form("%s_%s__bce2",var,nj.c_str()) );
    TH1D *bce3 = (TH1D*)gDirectory->Get( Form("%s_%s__bce3",var,nj.c_str()) );

    // All enri
    TH1D *enri1 = (TH1D*)gDirectory->Get( Form("%s_%s__enri1",var,nj.c_str()) );
    TH1D *enri2 = (TH1D*)gDirectory->Get( Form("%s_%s__enri2",var,nj.c_str()) );
    TH1D *enri3 = (TH1D*)gDirectory->Get( Form("%s_%s__enri3",var,nj.c_str()) );
    
    TH1D *bce;
    TH1D *enri;

    if(sel=="QCD"){
      if(bce1==0||bce2==0||bce3==0||
	 enri1==0||enri2==0||enri3==0) {
	cout << "couldn't plot as QCD reliso histograms not found."<< endl;
	return 0;
      }
      bce = (TH1D*)bce1->Clone("bce");
      bce->Add(bce2);
      bce->Add(bce3);
      enri = (TH1D*)enri1->Clone("enri");
      enri->Add(enri2);
      enri->Add(enri3);

      qcd->Rebin(this_rebin);
      bce->Rebin(this_rebin);
      enri->Rebin(this_rebin);
    }
    else{
      data->Rebin(this_rebin);
    }
    //    int  this_rebin = rebin[i-1];
    



    //---------------------
    // Rescale (for dijet) 
    //---------------------
    if(sel=="dijet") {
      if     (dijet_pt==15)   data->Scale( 3571 );
      else if(dijet_pt==30)   data->Scale(  288 );
      else if(dijet_pt==80)   data->Scale( 14.6 );
      else if(dijet_pt==170)  data->Scale( 0.29 );
      else if(dijet_pt==300)  data->Scale( 1.6E-2 );
      else if(dijet_pt==470)  data->Scale( 9.7e-4 );
      else if(dijet_pt==800)  data->Scale( 2.4e-5 );
      else if(dijet_pt==1400) data->Scale( 2.0e-7);      
    }

    //-----------------
    // Formatting
    //-----------------

    //if(i<5) qcd->SetTitle(Form("QCD, #slash{E}_{T}<20 GeV (%s)",jetlabel[i-1]));
    //else    qcd->SetTitle(Form("QCD, #slash{E}_{T}>20 GeV (%s)",jetlabel[i-1]));
    if(sel=="QCD"){
      qcd->SetTitle(Form("QCD (%s)",jetlabel[i-1]));
      if(i>5) qcd->SetTitle(Form("QCD (%s) (AES)",jetlabel[i-1]));
      
      qcd->GetXaxis()->SetTitle("RelIso");
      qcd->GetXaxis()->SetRangeUser(0,1);    
      qcd->SetMinimum(0);
      if(i<=5 && setyscale) {//NES
	if     (nj=="0j")  qcd->SetMaximum( 600);
	else if(nj=="1j")  qcd->SetMaximum(1100);
	else if(nj=="2j")  qcd->SetMaximum( 350);
	else if(nj=="3j")  qcd->SetMaximum(  65);
	else if(nj=="4mj") qcd->SetMaximum(  15);
      }
    }
    else {
      
      data->SetTitle(Form("%s (%s)", plotTitle.c_str(),jetlabel[i-1]));
      if(i>5) data->SetTitle(Form("%s (%s) (AES)", plotTitle.c_str(),jetlabel[i-1]));

      data->GetXaxis()->SetTitle("RelIso");
      data->GetXaxis()->SetRangeUser(0,1);    
      data->SetMinimum(0);
    }



    if(npad==nplot+1) npad=1; //reset
    c1.cd(npad);



    //-------------
    // Draw plots
    //-------------
    if(sel=="QCD"){
      qcd->Draw();

      // components
      bce->SetLineColor(kBlue);
      bce->SetFillColor(kBlue);
      //bce->SetLineStyle(kDashed);
      enri->SetLineColor(kRed);
      enri->SetFillColor(kRed);
      //enri->SetLineStyle(kDotted);
      bce->SetFillStyle(3005);
      enri->SetFillStyle(3004);
      
      enri->Draw("hist same");
      bce->Draw("hist same");
    }
    else{
      data->Draw();
    }



    //------
    // Save
    //------
    string out;
    if(sel=="QCD"){
      if(nplot==4) out = "reliso_QCD_1to4j";
      if(nplot==5) out = "reliso_QCD_0to4j";
    }
    else if(sel=="data"){
      if(nplot==4) out = "reliso_data_1to4j";
      if(nplot==5) out = "reliso_data_0to4j";
    }
    else{
      if(nplot==4) out = Form("reliso_dijet_pt%d_1to4j",dijet_pt) ;
      if(nplot==5) out = Form("reliso_dijet_pt%d_0to4j",dijet_pt) ;
    }

    if(i>5&&i<=10) {
      out = out + "_" + old_AES;
      //out = out + Form("_e%d",old_run_eleET_cut);
    }else{
      if(i>10) out = out + "_AES_" + planB;
    }


    
    //    char *outname = out.c_str();

    if(npad==nplot) {
      
      c1.SaveAs(Form("%s.gif",out.c_str()));
      // save as pdf
      /*
      c1.SaveAs(Form("%s.eps",outname));
      //convert eps to pdf
      gROOT->ProcessLine(Form(".!ps2pdf -dEPSCrop %s.eps",outname));
      gROOT->ProcessLine(Form(".!rm %s.eps",outname));
      */
    }
    npad++;    
  }
  c1.Close();

  return 0;
}
