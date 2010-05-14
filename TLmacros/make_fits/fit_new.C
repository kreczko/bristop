//
// 31 Mar 2010
//
// * Perform fit on RelIso of all (s+b) to estimate QCD
// * for a given range.
//
// (a) free functional fits in all jet bins
// (b) free fit in 1,2j, but constrained fits in 3,4mj, using 
//     fitted results from 1,2j.
// (c) plot RelIso without fitting.
//
//----------------------------------------------------
#include <iomanip>
#include <vector>

bool limit_gaus_mean_12j = false;
double gaus_mean_min = 0.3;
const double gaus_mean_max = 0.6;


const int nj = 4;
const double intlumi = 20; //pb-1

const bool debug =0;// true;

const double sig_from = 0; //new reliso
const double sig_upto = 0.1;
double fit_from;
double fit_upto;

double this_bw;

double nqcd_actual_sig[nj];
double n_extrap[nj];
double fit_chi2[nj];  
int    fit_ndf[nj];

string do_fit;

int nFreePara ;
double fit_par0[nj];
double fit_par1[nj];
double fit_par2[nj];
double fit_par3[nj];

int fit_new(){
  fit();
}
int fit(string do_fit_user="fix",
	string func_user="gaus", 
	double fit_from_user=0.2, 
	double fit_upto_user=1.4,
	bool   zoom_in = true) {
  

  double bw = 0.1;
  int rb = bw/0.01;
  
  
  string func = func_user;
  //  cout << "enter func (gaus/pol3): ";
  //  cin >> func; 
  //string func = "pol3";
  //string func = "pol4";
  //string func = "landau";

  //double fit_from = fit_from_user; //0.3;  // <------
  //double fit_upto = fit_upto_user; //1.0;  // <------

  do_fit = do_fit_user;

  fit_from = fit_from_user; //0.3;  // <------
  fit_upto = fit_upto_user; //1.0;  // <------


  const int rebin[4] = { 5, 10,  10, 10 };  //<---- 

  gStyle->SetTitleH(0.1);
  gStyle->SetStatH(0.22); //0.24);
  gStyle->SetStatW(0.22); //0.26);
  gStyle->SetOptStat(1); //on stat
  //  gStyle->SetTitleOffset(-0.05,"X"); //on stat
  //  gStyle->SetTitleSize(0.06); //axis title
  
  if(do_fit=="none")  gStyle->SetOptStat(0);//off stat

  //  gStyle->SetOptStat(1110);//off title


  TFile f("../test_mc_mixture.root");



  // fixPara: Fit parameters
  if(do_fit=="fix"){
    TF1 ff(func.c_str(),func.c_str(),0,1);
    nFreePara = ff.GetNumberFreeParameters();
  }


  TCanvas *c1 = new TCanvas("c1","c1",980,700);
  c1->Divide(2,2);


  char *njet[4] = {"1j","2j","3j","4mj"} ;
  char *njlabel[4] = {"1 jet","2 jets","3 jets","#geq4 jets"} ;


  if(debug) cout << "About to enter nj loop"<< endl;



  for (int j=0; j<4; j++){ 

    if(debug) {
      cout << "\n Njet: " << j+1 << endl;
      cout << "----------" << endl;
    }
    TH1D *all   = (TH1D*)f.Get(Form("QCD_estimation/QCDest_CombRelIso_%s__data", njet[j]));
    TH1D *qcd   = (TH1D*)f.Get(Form("QCD_estimation/QCDest_CombRelIso_%s__QCD",  njet[j]));
    TH1D *enri1 = (TH1D*)f.Get(Form("QCD_estimation/QCDest_CombRelIso_%s__enri1",njet[j]));
    TH1D *enri2 = (TH1D*)f.Get(Form("QCD_estimation/QCDest_CombRelIso_%s__enri2",njet[j]));
    TH1D *enri3 = (TH1D*)f.Get(Form("QCD_estimation/QCDest_CombRelIso_%s__enri3",njet[j]));
    TH1D *bce1  = (TH1D*)f.Get(Form("QCD_estimation/QCDest_CombRelIso_%s__bce1", njet[j]));
    TH1D *bce2  = (TH1D*)f.Get(Form("QCD_estimation/QCDest_CombRelIso_%s__bce2", njet[j]));
    TH1D *bce3  = (TH1D*)f.Get(Form("QCD_estimation/QCDest_CombRelIso_%s__bce3", njet[j]));


    const double binw = 0.01; //original binw = 0.01 = 1.1/110 (old) = 10/1000 (new)
      
    const int sig_bin_low = (int)(sig_from/binw + 1);
    const int sig_bin_up  = (int)(sig_upto/binw);

    nqcd_actual_sig[j] = qcd->Integral( sig_bin_low, sig_bin_up ); //actual number of _QCD_ event in signal region
    
    if(1||debug) {
      printf("sig_bin_low  =  %.2f\n", sig_bin_low );
      printf("sig_bin_up  =  %.2f\n", sig_bin_up );
      printf("qcd->Integral( sig_bin_low, sig_bin_up )  =  %.2f\n", nqcd_actual_sig[j] );
    }


    if(debug) cout << "amtb 1"<< endl;

    all->Rebin(rebin[j]);
    qcd->Rebin(rebin[j]);
    enri1->Rebin(rebin[j]);
    enri2->Rebin(rebin[j]);
    enri3->Rebin(rebin[j]);
    bce1->Rebin(rebin[j]);
    bce2->Rebin(rebin[j]);
    bce3->Rebin(rebin[j]);

    // QCD component
    THStack *qcdcomp = new THStack("qcdcomp","QCD components");
    enri1->SetFillColor(kAzure+2); //normal
    enri2->SetFillColor(kAzure+1);
    enri3->SetFillColor(kAzure-9);
    bce1->SetFillColor(kSpring+3); //kGreen+2 darkest
    bce2->SetFillColor(kSpring+2); //kGreen 8
    bce3->SetFillColor(kSpring+1); //kSpring-9
    qcdcomp->Add(bce1);
    qcdcomp->Add(bce2);
    qcdcomp->Add(bce3);
    qcdcomp->Add(enri1);
    qcdcomp->Add(enri2);
    qcdcomp->Add(enri3);

    //    double this_bw = binw*rebin[j];
    this_bw = binw*rebin[j];



    if(debug) cout << "amtb 2"<< endl;

    all->GetXaxis()->SetRangeUser(0,1.6-0.01);
    qcd->GetXaxis()->SetRangeUser(0,1.6-0.01);

    all->SetMarkerStyle(20);
    //    all->SetMarkerSize(0.8); //if run in batch mode

    all->SetLineColor(kRed);
    qcd->SetLineColor(kAzure+2);

    qcd->SetFillColor(kAzure-9);

    all->SetLineWidth(2);
    qcd->SetLineWidth(2);

    all->SetTitle(Form("%s",njlabel[j]));
    all->SetXTitle("RelIso");
    //    all->GetXaxis()->SetTitleSize(0.2);

    //all->SetName(func.c_str());
    all->SetName(Form("%s (%.1f-%.1f)", func.c_str(), fit_from, fit_upto));

    all->SetLineColor(kBlack);

    c1->cd(j+1);

    gStyle->SetOptFit(112);

    // To zoom in on y-axis
    if(zoom_in) {

      // find tallest bin from 3rd to 9th (0.2-1.0)
      double highest = 0;
      unsigned int scanfrom = all->FindBin(0.2);
      unsigned int scanto   = all->FindBin(1.0) - 1;//do not count bin 1.0
      //cout << "find highest bin from relios 0.2 to 1.0, bin: "<< scanfrom << "-"<< scanto << endl;
      for (short i=scanfrom; i<scanto; i++) {
	double bin_height  = all->GetBinContent(i);
	//cout << "bin " <<  i << "  " << bin_height << endl;
	if( bin_height > highest ) highest = bin_height;
      }
      //cout << "highest: " << highest <<endl;
      all->GetYaxis()->SetRangeUser( 0, highest*2.3  );
    
      /*
      if(j+1==1) all->GetYaxis()->SetRangeUser(0,3000); //1j
      if(j+1==2) all->GetYaxis()->SetRangeUser(0,800); //2j
      if(j+1==3) all->GetYaxis()->SetRangeUser(0,200); //3j
      if(j+1==4) all->GetYaxis()->SetRangeUser(0,60); //4mj
      */
    }

    if(debug) cout << "amtb: draw histo"<< endl;

    //    all->Draw("ahist");
    gPad->SetBottomMargin(0.14); //TEST
    all->GetXaxis()->SetTitleOffset(0.8); //on stat

    all->GetXaxis()->SetTitleSize(0.07);

    all->Draw("ahist"); 
    qcdcomp->Draw("ahist same"); //add
    qcd->Draw("ae same");
    all->Draw("ae same");
    all->Draw("axis same");
    


    if(debug) cout << "amtb 4"<< endl;




    if(do_fit=="none"){
      continue;
    }

    if(debug) cout << "perform the fit" << endl;
    //--------------------
    //  Perform the Fit
    //--------------------
    //    all->Fit(func.c_str(),"Q0","ah", fit_from, fit_upto);

    //    TF1 *myf;

    if (do_fit=="free") {
      
      //----------
      // Free Fit
      //----------
      if(debug) cout << "free fit"<< endl;

      //    TF1 *fitf = new TF1("landau","landau",0,2);
      //1st opt: "V": verbose, "Q": Quiet, "0": do not draw
       all->Fit(func.c_str(),"Q0","ah", fit_from, fit_upto);

       //cout << "after fit:  all->ls()"<< endl;
       //all->ls();


       getFitRes(j, gPad, func,all);

    }
    else if(do_fit=="fix") {

      //-------------------
      //  Constrained Fit
      //-------------------
      if(debug) cout << "constrained fit" << endl;
    
      TF1 *fitf;
      if(func=="gaus") {
	fitf = new TF1("gaus","gaus",0,2);
      } 
      else if (func=="pol3") {
	// pol3 = A + Bx + Cx^2+ Dx^3 = a(1 + b x + c x2 + d x3)
	//
	fitf = new TF1("pol3","[0] * ( 1 + [1]*x + [2]*x^2 + [3]*x^3 )",0,2);
      }
      else if (func=="landau"){
	fitf = new TF1("landau","landau",0,2);
      }

      // (a) free-fit in 1,2 jet bins
      //
      // For gaussian, require mean to be whithin 0.2-0.6, the flat or peak region.
      //
      if( j<2 ) { //1,2j - Free Fit OR gaus-mean-constrained

	if(func=="gaus" && limit_gaus_mean_12j){
	  
	  cout << "constraining gaus mean in12j to "<< gaus_mean_min << "-" << gaus_mean_max << endl;
	  fitf->SetParLimits( 1, gaus_mean_min, gaus_mean_max ); //mean of gaus (1j)

	}

	//      all->Fit(func.c_str(),"V","ah", fit_from, fit_upto); //Free Fit
	// 1st opt: "V"=verbose; "Q"=quite
	all->Fit(fitf,"Q0","ah", fit_from, fit_upto); //Free Fit (Q)

      } else { //3,4mj - Constrained Fit
	
	string fopt ;
	if ( func=="gaus" ) {
	  //fitf = new TF1(func.c_str(),func.c_str(),0,2);
	  fitf->FixParameter( 1, (fit_par1[0]+fit_par1[1])/2 ); //mean or p1
	  fitf->FixParameter( 2, (fit_par2[0]+fit_par2[1])/2 ); //sigma or p2      
	  fopt="BQ0";
	} else if ( func=="pol3" ) {
	  // fix b,c,d
	  double avg_b = (fit_par1[0] + fit_par1[1]) / 2;
	  double avg_c = (fit_par2[0] + fit_par2[1]) / 2;
	  double avg_d = (fit_par3[0] + fit_par3[1]) / 2;
	  
	  fitf->FixParameter( 1, avg_b );
	  fitf->FixParameter( 2, avg_c );
	  fitf->FixParameter( 3, avg_d );
	  fopt="BQ0";
	}    
	else if(func=="landau"){
	  double mpv_min = min(fit_par1[0], fit_par1[1]);//mpv or p1
	  double mpv_max = max(fit_par1[0], fit_par1[1]);
	  cout << "mpv_min " << mpv_min << endl;
	  cout << "mpv_max " << mpv_max << endl;
	  //cout << "mpv_avg " << (mpv_min+mpv_max)/2 << endl;
	  //cout << "sigma_avg " << (fit_par2[0]+fit_par2[0])/2 << endl;
	  //fitf->SetParameters( 10, (mpv_min+mpv_max)/2, (fit_par2[0]+fit_par2[0])/2);
	  fitf->SetParLimits( 1, mpv_min, mpv_max );
	  //fitf->SetParLimits( 1, 0.37, 0.44 );
	  fopt="Q0"; //"B" gives wrong answer
	}
  
	// "B"=use fixed parameter value
	all->Fit(fitf,fopt.c_str(),"ah", fit_from, fit_upto);
      }
      delete fitf;

      if(debug) cout << "getFitRes 1" << endl;
      getFitRes(j, gPad, func,all);

    }
  

    if(debug) cout << "amtb go to next nj"<< endl;
    
    
  }// njet loop

  if(debug) cout << "amtb: after nj loop"<< endl;

  // Legend
  //--------
  if(debug) cout << "amtb: creating legend"<< endl;

  double coord[4] = {0.33,0.5,0.584,0.9};
  TLegend leg( coord[0], coord[1], coord[2], coord[3] );

  
  TF1 *blue=new TF1("blue","pol0",0,1);
  TF1 *red=new TF1("red","pol0",0,1);
  blue->SetLineColor(kBlue);
  red->SetLineColor(kRed);
  
  leg.SetFillColor(0);
  leg.AddEntry(all,  "All events (S+B)","LP"); //Line+Marker, E(error)
  
  leg.AddEntry(enri1,"EME pt 20-30","f");
  leg.AddEntry(enri2,"EME pt 30-80","f");
  leg.AddEntry(enri3,"EME pt 80-170","f");
  leg.AddEntry(bce1, "BCE pt 20-30","f");
  leg.AddEntry(bce2, "BCE pt 30-80","f");
  leg.AddEntry(bce3, "BCE pt 80-170","f");
  
  if(do_fit != "none"){
    leg.AddEntry(red,  "Fit","l");
    leg.AddEntry(blue, "Extrapolation","l");
  }else{
    leg.SetX1(coord[0]+0.3);//move to right
    leg.SetX2(coord[2]+0.3);
  }
  if(debug) cout << "amtb: draw legend"<< endl;
  c1->cd(1); leg.Draw();
  c1->cd(2); leg.Draw();
  c1->cd(3); leg.Draw();
  c1->cd(4); leg.Draw();  
  

  if(do_fit=="none"){

    if(debug) cout << "amtb: save plot"<< endl;
    if(!zoom_in)  {
      c1->SaveAs("QCD_reliso.gif");
      c1->SaveAs("QCD_reliso.pdf");
    }
    else {
      c1->SaveAs("QCD_reliso_zoom.gif");
      c1->SaveAs("QCD_reliso_zoom.pdf");
    }
    c1->Close();
    return 0;
  }



  if(debug) cout << "amtb: print fit results"<< endl;

  //TDatime now;
  //now.Print();
  cout << "\n-----------------------------------------------"<< endl;
  cout << "Fit " << func << " (range: " << fit_from << "-"<< fit_upto << ")" <<endl;
  cout << "-----------------------------------------------"<< endl;
  cout << "  " << intlumi << "/pb     rb      True    Estimate     Diff" << endl;
  printf("   1 jet:   %2d %10.1f  %10.1f  %6.1f %%\n", rebin[0], nqcd_actual_sig[0], n_extrap[0], (n_extrap[0]/nqcd_actual_sig[0]-1)*100 );
  printf("   2 jet:   %2d %10.1f  %10.1f  %6.1f %%\n", rebin[1], nqcd_actual_sig[1], n_extrap[1], (n_extrap[1]/nqcd_actual_sig[1]-1)*100 ); 
  printf("   3 jet:   %2d %10.1f  %10.1f  %6.1f %%\n", rebin[2], nqcd_actual_sig[2], n_extrap[2], (n_extrap[2]/nqcd_actual_sig[2]-1)*100 );
  printf(" >=4 jet:   %2d %10.1f  %10.1f  %6.1f %%\n", rebin[3], nqcd_actual_sig[3], n_extrap[3], (n_extrap[3]/nqcd_actual_sig[3]-1)*100 );
  cout << "-----------------------------------------------\n"<< endl;
  
  if(debug) cout << "amtb: print constrained-fit results"<< endl;

  if (do_fit=="fix"){
     
    cout << " Fitted parameter values" << endl;
    cout << "-----------------------------------------------"<< endl;
    if     (func=="gaus")   cout << " Gaus         const        mean       sigma" << endl;
    else if(func=="landau") cout << " Landau       const        mpv        sigma" << endl;
    else                    cout << "           p0           p1          p2         p3" << endl;
    cout << "-----------------------------------------------"<< endl;
    char *jjlabel[4] = {"1","2","3",">=4"};
    for(short i=0;i<4;++i){
      if(nFreePara==3)    
	printf(" %3s jet: %10.1f %10.2f %10.2f\n",jjlabel[i],fit_par0[i],fit_par1[i],fit_par2[i]);
      else if(nFreePara==4)    
	printf(" %3s jet: %10.1f %10.2f %10.2f %10.2f\n",jjlabel[i],fit_par0[i],fit_par1[i],fit_par2[i],fit_par3[i]);
    }
    cout << "-----------------------------------------------"<< endl;
    
  }


  if(debug) cout << "amtb: save picture"<< endl;

  //c1->ls();
  //c1->SaveAs("test.gif");
  string out = Form("%s_r%.1fto%.1f",func.c_str(), fit_from,fit_upto);
  if(do_fit=="fix") out = Form("fixPara_%s_r%.1fto%.1f",func.c_str(), fit_from,fit_upto);

  if(!zoom_in) {
    c1->SaveAs(Form("fit_%s.gif",     out));
    c1->SaveAs(Form("fit_%s.pdf",     out));
  } else {
    c1->SaveAs(Form("fit_zoom_%s.gif",out));
    c1->SaveAs(Form("fit_zoom_%s.pdf",out));
  }
  c1->Close();

  //  cout << "amtb 123 b"<< endl;

  if(debug) cout << "amtb: write result to text file"<< endl;

  ofstream myfile;

  string outText;
  outText = "est_";
  if(zoom_in) outText += "zoom_";
  outText = outText + do_fit + "_" + func;
  
  if(debug) cout << "amtb: " << outText << endl;

  myfile.open(Form("%s.txt",outText.c_str()), ios::app);
  
  //if(!zoom_in) myfile.open(Form("est_%s.txt",func.c_str()),ios::app);
  //else         myfile.open(Form("est_zoom_%s.txt",func.c_str()),ios::app);
  myfile.setf(ios::fixed,ios::floatfield);
  myfile << endl;
  // for free fit, write out 1-4j  
    for(int k=0; k<4; ++k){    myfile << n_extrap[k]/nqcd_actual_sig[k]-1 << endl;  }
    /*
  if(do_fit=="free") {
    for(int k=0; k<4; ++k){    myfile << n_extrap[k]/nqcd_actual_sig[k]-1 << endl;  }
  }
  else {
    // for constrained fit, write out 3,4j only
    for(int k=2; k<4; ++k){    myfile << n_extrap[k]/nqcd_actual_sig[k]-1 << endl;  }
  }
    */
  if(debug) cout << "amtb: close text file" << endl; 
  myfile.close();

  if(debug) cout << "amtb: complete"<< endl;
  
  return 0;
}

//-------------------------------------------------------------------------------
// just plot without fitting
void fit_none(bool zoom_in=true){
  fit( "none", "", 0, 1, zoom_in );
}

//-------------------------------------------------------------------------------
// free fit
void fit(string func_user="gaus", 
	 double fit_from_user=0.2, 
	 double fit_upto_user=1.6,
	 bool   zoom_in = true) {
  fit("free",func_user,fit_from_user,fit_upto_user,zoom_in);
}

//-------------------------------------------------------------------------------
// constrained fit in 3,4j (Free in 1,2j)
void fit_fixPara(string func_user="gaus", 
		 double fit_from_user=0.2, 
		 double fit_upto_user=1.6,
		 bool   zoom_in = true) {
  fit("fix",func_user,fit_from_user,fit_upto_user,zoom_in);
}
//-------------------------------------------------------------------------------
// constrained fit in 3,4j (gaus-mean-constrained in 1,2j > 0.3)
void fit_fixPara_mean12j_geq03(string func_user="gaus", 
			       double fit_from_user=0.2, 
			       double fit_upto_user=1.6,
			       bool   zoom_in = true) {
  limit_gaus_mean_12j = true;
  gaus_mean_min = 0.3;
  fit("fix",func_user,fit_from_user,fit_upto_user,zoom_in);
}
//-------------------------------------------------------------------------------
// constrained fit in 3,4j (gaus-mean-constrained in 1,2j > 0.4)
void fit_fixPara_mean12j_geq04(string func_user="gaus", 
			       double fit_from_user=0.2, 
			       double fit_upto_user=1.6,
			       bool   zoom_in = true) {
  limit_gaus_mean_12j = true;
  gaus_mean_min = 0.4;
  fit("fix",func_user,fit_from_user,fit_upto_user,zoom_in);
}

//-------------------------------------------------------------------------------
void getFitRes(int j, TPad *gpad, string func, TH1D *all){

  if(debug) cout << "\n start of getFitRes\n"<< endl;

  // Get fitted function
  //---------------------
  if(debug) cout << "get fitted function"<< endl;

  TF1 *myf = all->GetFunction(func.c_str());
  myf->SetLineColor(kRed);

  if(debug) cout << "clone myf"<< endl;
  TF1 *myf2 = (TF1*)myf->Clone(); //range 0-0.1
  myf2->SetLineColor(kBlue);
  myf2->SetRange( 0, 0.1 );

  TF1 *myf3 = (TF1*)myf->Clone(); //range 0.1 to 0.2
  myf3->SetLineColor(kBlue);
  myf3->SetLineStyle(kDashed);
  myf3->SetRange( 0.1, fit_from );

  if(debug) cout << "amtb: go to pad"<< endl;
  gpad->cd();

  if(debug) cout << "amtb: draw functions"<< endl;
  myf->Draw("same");
  myf2->Draw("same");
  myf3->Draw("same");


  // Get fit results
  //----------------
  if(debug) cout << "amtb: get fit results"<< endl;
  fit_chi2[j] = myf->GetChisquare();
  fit_ndf[j]  = myf->GetNDF();
  n_extrap[j] = myf->Integral( sig_from, sig_upto ) / this_bw;


  if(do_fit=="fix") getParametersForConstrainedFit(j,myf);

  if(debug) cout << "\n end of getFitRes\n"<< endl;
    
}

//-------------------------------------------------------------------------------
void getParametersForConstrainedFit(int j, TF1 *myf){
 
  // fixPara: Get fit results
  //--------------------------
  fit_par0[j] = myf->GetParameter(0);//1st para p0
  fit_par1[j] = myf->GetParameter(1);//2nd para p1
  if(nFreePara>=3) {
    fit_par2[j] = myf->GetParameter(2);//3rd para p2
    if(nFreePara>=4) {
      fit_par3[j] = myf->GetParameter(3);//4th para p3
    }
  }
}
//-------------------------------------------------------------------------------
