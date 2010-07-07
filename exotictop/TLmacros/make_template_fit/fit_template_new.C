//
// 19 Mar 2010
// Plan B, QCD template fit.
//  
// - Use control sample as template.
// - Then fit to signal region
// - template from 1j-control to fit 1j-signal, etc
//---------------------------------------------------

#include <iomanip>

using namespace RooFit;

const int nj = 4;
//string jetbin_to_use_if_using_one_template = "4j";

//------------------------------------------------------------------------
// [#0] WARNING:InputArguments -- RooAbsPdf::fitTo(model) WARNING: a likelihood fit is request of what appears to be weighte
// d data. 
// While the estimated values of the parameters will always be calculated taking the weights into account, 
//        there are multiple ways to estimate the errors on these parameter values. You are advised to make an 
//   explicit choice on the error calculation: 
//   - Either provide SumW2Error(kTRUE), to calculate a sum-of-weights corrected HESSE error matrix 
//   (error will be proportional to the number of events)
//   - Or provide SumW2Error(kFALSE), to return errors from original HESSE error matrix
//   (which will be proportional to the sum of the weights)
//   If you want the errors to reflect the information contained in the provided dataset, choose kTRUE. 
//        If you want the errors to reflect the precision you would be able to obtain with an unweighted dataset 
//   with 'sum-of-weights' events, choose kFALSE.
//------------------------------------------------------------------------
bool use_dijet_template = false;

// Fit using 1 template from control region (use dijet pt15)
void fit_template_1T_dijet( string control_sample="B3_e20", string jetbin="1j", 
			    double fit_from, double fit_upto ){
  use_dijet_template = true;
  fit_template_new(control_sample, true, jetbin, fit_from, fit_upto);
}

// Fit using 1 template from control region (use usual s+b (enri+bce))
void fit_template_1T( string control_sample="B3_e20", string jetbin="2j", 
		      double fit_from, double fit_upto ){
  fit_template_new(control_sample, true, jetbin, fit_from, fit_upto);
}

// Fit using 4 templates from control region
void fit_template_4T( string control_sample="B3_e20", double fit_from, double fit_upto ){
  fit_template_new(control_sample, false, "", fit_from, fit_upto);
}


// MAIN code
void fit_template_new( string control_sample ="B3_e20",
		       bool use_one_template = true,
		       string jetbin_to_use_if_using_one_template = "2j",
		       double fit_from_user = 0.2, 
		       double fit_upto_user = 1.2 ){

  
  cout << "---------------------------"<< endl;
  cout << "control sample:  " << control_sample << endl;
  cout << "---------------------------"<< endl;
  if(use_one_template)  
    cout << "Using one template from " << jetbin_to_use_if_using_one_template<< endl;
  else                  
    cout << "Using four templates from 1-4j." << endl;  
  cout << "---------------------------"<< endl;

  //-----------------
  // RooFit setting
  //-----------------
  // error type (Def=Poission for integer histogram bin content)
  // ref: http://root.cern.ch/root/html/RooAbsData.html
  //  RooAbsData::Poission
  //  RooAbsData::SumW2
  //  RooAbsData::None
  //  RooAbsData::Auto
  RooAbsData::ErrorType etype = RooAbsData::SumW2; 
  //-----------------

  gStyle->SetOptStat(1111111);
  gStyle->SetStatW(0.35);
  
  double fit_from = fit_from_user;
  double fit_upto = fit_upto_user;

  int    rebin[4] = { 10, 10,  10, 10 };  //<---- 


  cout << "\n\n\n====================\n";
  cout << "\n Range:  " << fit_from << " to " << fit_upto;
  //  cout << "\n\n Fit option:  " << fopt ;
  cout << "\n\n====================\n\n"<< endl;




  TFile f("../test_mc_mixture.root");
  //TFile f2("../test_mc_mixture.root"); 
  // use QCD template from dijet_pt15
  TFile *f2;
  if(use_dijet_template){
    f2 = new TFile("~/cms/analysis/results/2010_Feb/25Feb_7TeV_plainQCD_pt15/test_data.root");
  }else{
    f2 = new TFile("../test_mc_mixture.root");
  }


  TCanvas c1("c1","QCD template fit",980,700);
  c1.Divide(2,2);
 
  const double tiny = 1e-6;

  string njetA[4] = {"1j","2j","3j","4mj"} ;
  char *njlabel[4] = {"1 jet","2 jets","3 jets","#geq4 jets"} ;


  double nQCDInSignal_true[nj]; //double nqcd_actual_sig[nj];
  double nQCDInSignal_est[nj]; //n_extrap[nj];
  double nQCDInSignal_estErr[nj]; //fit error on est
  double Deviation[nj];

  for(int a=0; a<4; a++){
    nQCDInSignal_true[a]   = 0;
    nQCDInSignal_est[a]    = 0;
    nQCDInSignal_estErr[a] = 0;
    Deviation[a]           = 0;
  }


  for(int j=0; j<4; j++) {


    cout << "\n-------------" << endl;
    cout << " fitting " <<  njlabel[j];
    cout << "\n-------------" << endl;



    char *njet = njetA[j].c_str();




    ///----------------------
    ///  Fit Control sample
    ///----------------------

    cout << "\n(a) Create template from control sample:" << endl;

    string templateBin = jetbin_to_use_if_using_one_template.c_str(); // use this template throughout
   
    if(!use_one_template) { //use 4 templates
      templateBin = njet;
    }

    TH1F *control = (TH1F*)f2->Get(Form("QCD_planB/QCDest_CombRelIso_AES_plan%s_%s__data",
				      control_sample.c_str(), templateBin.c_str()));

    cout << " integral of control = " << control->Integral() << endl;
    control->GetXaxis()->SetRangeUser(0,fit_upto-tiny);
    cout << " integral of control (after setrange) = " << control->Integral() << endl;
    //  control->DrawCopy();



    if( control->GetNbinsX() > 100) {
      control->Rebin(rebin[0]); 
    }

    cout << "integral of control (after rebin) = " << control->Integral() << endl;
    // after rebinning
  
    control->GetXaxis()->SetRangeUser(0,fit_upto-tiny);
    cout << "integral of control (after rebin+setrange) = " << control->Integral() << endl;
    //control->DrawCopy("h");



    // create pdf
    //c1.cd(3);
    //control->GetXaxis()->SetRangeUser(fit_from,fit_upto);
    //    RooRealVar reliso("reliso","reliso", 0, fit_upto);
    RooRealVar reliso("reliso","reliso", 0, 1.6);
    reliso.setBins(12);
    //control->GetXaxis()->SetRangeUser(0,fit_upto-tiny);
    //cout << "integral ori histo qcd = " << control->Integral() << endl;
  
  
    RooDataHist rh_qcd("rh_qcd","qcd", reliso, control);
    cout << "RooDataHist control rh_qcd" << endl;
    cout << "sumEntries, reliso 0-0.1: " << rh_qcd.sumEntries("reliso>0&&reliso<0.1") << endl;
    cout << "sumEntries, reliso 0.1-0.2: " << rh_qcd.sumEntries("reliso>0.1&&reliso<0.2") << endl;
  

    RooHistPdf pdf_qcd( "pdf_qcd", "QCD pdf ", reliso, rh_qcd, 0);




    ///----------------------
    ///  Fit Signal sample
    ///----------------------


    cout << "\n(b) Fitting Signal sample" << endl;


    /*
    // create pdf
    //c1.cd(3);
    //control->GetXaxis()->SetRangeUser(fit_from,fit_upto);
    //    RooRealVar reliso("reliso","reliso", 0, fit_upto);
    RooRealVar reliso("reliso","reliso", 0, 1.6);
    reliso.setBins(12);
    //control->GetXaxis()->SetRangeUser(0,fit_upto-tiny);
    //cout << "integral ori histo qcd = " << control->Integral() << endl;
    
  
    RooDataHist rh_qcd("rh_qcd","qcd", reliso, control);
    cout << "RooDataHist control rh_qcd" << endl;
    cout << "sumEntries, reliso 0-0.1: " << rh_qcd.sumEntries("reliso>0&&reliso<0.1") << endl;
    cout << "sumEntries, reliso 0.1-0.2: " << rh_qcd.sumEntries("reliso>0.1&&reliso<0.2") << endl;
    */



    const double binw = 0.01;// original bin width;
    double this_bw = binw*rebin[0];
    
    const double sig_from = 0;
    const double sig_upto = 0.1;
    const int sig_bin_low = (int)(sig_from/this_bw+0.1);
    const int sig_bin_up  = (int)(sig_upto/this_bw+0.1);
    const int fit_bin_low = (int)(fit_from/this_bw+0.1);
    const int fit_bin_up  = (int)(fit_upto/this_bw+0.1);
    
    cout << "\n sig_bin_low    " << sig_bin_low << endl;
    cout << " sig_bin_up     " << sig_bin_up << endl;
    cout << " fit_bin_low    " << fit_bin_low << endl;
    cout << " fit_bin_up     " << fit_bin_up << endl << endl;
    
    /*
      double fracLOW  = control_norm_subrange->Integral(sig_bin_low+1, sig_bin_up);
      double fracGAP  = control_norm_subrange->Integral(sig_bin_up +1, fit_bin_low);
      double fracHIGH = control_norm_subrange->Integral(fit_bin_low+1, fit_bin_up);
      double fracALL  = control_norm_subrange->Integral(sig_bin_low+1, fit_bin_up);
    */
    
    
    ///---------------------
    ///  Fit signal sample
    ///---------------------
    

    // Prepare signal sample
    bool do_closure_test = 0;//true;
    TH1F *ho_all;
    TH1F *ho_qcd;
    if(do_closure_test) {
      ho_all = (TH1F*)control->Clone();
      ho_qcd = (TH1F*)control->Clone();
    } else {
      ho_all = (TH1F*)f.Get(Form("QCD_estimation/QCDest_CombRelIso_%s__data",njet)); 
      ho_qcd = (TH1F*)f.Get(Form("QCD_estimation/QCDest_CombRelIso_%s__QCD",njet)); 
    }
    
    //  all->GetXaxis()->SetRangeUser(0,fit_upto-tiny);
    //  hqcd->GetXaxis()->SetRangeUser(0,fit_upto-tiny);
    
    ho_qcd->SetFillColor(kAzure-9);
    ho_qcd->SetFillStyle(3001);


    // Rebin 
    if(!do_closure_test) {
      ho_all->Rebin(rebin[0]);
      ho_qcd->Rebin(rebin[0]);
    }
    
    
    TH1F *all = (TH1F*)ho_all->Clone();
    TH1F *hqcd = (TH1F*)ho_qcd->Clone();
    /*
      TH1F *all = new TH1F("all",ho_all->GetTitle(),   fit_bin_up,0,fit_upto);
      TH1F *hqcd = new TH1F("hqcd",ho_qcd->GetTitle(), fit_bin_up,0,fit_upto);
      for(int i=1; i<=fit_bin_up; ++i){
      all->SetBinContent(i,ho_all->GetBinContent(i));
      hqcd->SetBinContent(i,ho_qcd->GetBinContent(i));
      all->SetBinError(i,ho_all->GetBinError(i));
      hqcd->SetBinError(i,ho_qcd->GetBinError(i));
      }
      hqcd->SetFillColor(kAzure-9);
      hqcd->SetFillStyle(3001);
    */
    
    // after rebin
    
    // build model template p.d.f.
    cout << "\n Build model p.d.f\n"<< endl;
    //  RooHistPdf model("model","model", reliso, rh_qcd, 0);
    
    const double nevent_observed_in_fit_range = all->Integral(fit_bin_low+1,fit_bin_up);
    cout << "\n fit_bin_low = " << fit_bin_low << endl;
    cout << "\n fit_bin_up  = " << fit_bin_up << endl;
    cout << "\n nevent_observed_in_fit_range = " << nevent_observed_in_fit_range << endl ;
    cout << "\n setting limit on nqcd to (0,2*nevent_observed_in_fit_range)" << endl;
    // set upper limit to twice nevent
    RooRealVar nqcd("nqcd", "number of QCD events", 0, nevent_observed_in_fit_range*2, "event");
    
    RooAddPdf model("model","model", pdf_qcd, nqcd);
    
    
    // Define reference ranges in reliso
    reliso.setRange("whole",       0,  fit_upto);
    reliso.setRange("signal",      0,    0.1   );
    reliso.setRange("gap",       0.1,  fit_from);
    reliso.setRange("extrapolate", 0,  fit_from);
    reliso.setRange("fit",   fit_from, fit_upto);
    
    cout << "\n----------------\n"<< endl;
    cout << "\n Find model integrals: model.createIntegral\n"<< endl;
    cout << "\n----------------\n"<< endl;
    
    double m0 = (model.createIntegral(reliso, NormSet(reliso), Range("whole")))->getVal();
    double m1 = (model.createIntegral(reliso, NormSet(reliso), Range("signal")))->getVal();
    double m2 = (model.createIntegral(reliso, NormSet(reliso), Range("gap")))->getVal();
    double m3 = (model.createIntegral(reliso, NormSet(reliso), Range("fit")))->getVal();
    double m4 = (model.createIntegral(reliso, NormSet(reliso), Range("extrapolate")))->getVal();
    double mm4 = (model.createIntegral(reliso, Range("fit")))->getVal();
    
    cout << " integral 0 to 1.2 = "<< m0 <<endl; 
    cout << " integral 0 to 0.1 = "<< m1 <<endl;
    cout << " integral 0.1 to 0.2 = "<< m2 <<endl;
    cout << " integral 0 to 0.2 = "<< m3 <<endl;
    cout << " integral 0.2 to 1.2 = "<< m4 <<endl;
    cout << " integral 0.2 to 1.2 = "<< mm4 <<endl;
    
    
    cout << "sum 0 to 1.2 = "  << m1+m2+m4 << endl;
    cout << "\n----------------\n\n"<< endl;
    
    
    /*
    double nAll_all = all->Integral() ;
    double nAll_low = all->Integral(1,1) ; //(0,0.1)->(1,1) 
    double nAll_gap = all->Integral(2,2) ; //(0.1,0.2)->(2,2) = int(0.1/binw+0.1)+1, same
    double nAll_high = all->Integral(3,12) ; //(0.3,1.2)->(3,12)
    double nAll_sub = all->Integral(1,12) ; //(0.1,1.2)->(1,2)
    
    cout << endl;
    cout << "nEvent in Signal sample, true, [  0,0.1] = " << nAll_low  << endl;
    cout << "nEvent in Signal sample, true, [0.1,0.2] = " << nAll_gap << endl;
    cout << "nEvent in Signal sample, true, [0.2,1.2] = " << nAll_high << endl;
    cout << "nEvent in Signal sample, true, [  0,1.2] = " << nAll_sub  << endl;
    */

    double nAll_low  = all->Integral(sig_bin_low+1,sig_bin_up) ; //(0,0.1)->(1,1) 
    double nAll_gap  = all->Integral(sig_bin_up+1, fit_bin_low) ; //(0.1,0.2)->(2,2) = int(0.1/binw+0.1)+1, same
    double nAll_high = all->Integral(fit_bin_low+1,fit_bin_up) ; //(0.3,1.2)->(3,12)
    double nAll_sub  = all->Integral(sig_bin_low+1,fit_bin_up) ; //(0.1,1.2)->(1,2)

    
    printf("nTotal event in Signal sample (sig) [ %.1f, %.1f ] = %6g\n", sig_from, sig_upto, nAll_low );
    printf("nTotal event in Signal sample (gap) [ %.1f, %.1f ] = %6g\n", sig_upto, fit_from, nAll_gap );
    printf("nTotal event in Signal sample (fit) [ %.1f, %.1f ] = %6g\n", fit_from, fit_upto, nAll_high);
    printf("nTotal event in Signal sample (sub) [ %.1f, %.1f ] = %6g\n", sig_from, fit_upto, nAll_sub );


    double nQCDInSignal_true_all = hqcd->Integral() ;    
    /*
    double nQCDInSignal_true_all = hqcd->Integral() ;
    double nQCDInSignal_true_low = hqcd->Integral(1,1) ;
    double nQCDInSignal_true_gap = hqcd->Integral(2,2) ;
    double nQCDInSignal_true_high = hqcd->Integral(3,12) ;
    double nQCDInSignal_true_sub = hqcd->Integral(1,12) ;    

    cout << endl;
    cout << "QCD in Signal sample, true, [  0,0.1] = " << nQCDInSignal_true_low  << endl;
    cout << "QCD in Signal sample, true, [0.1,0.2] = " << nQCDInSignal_true_gap << endl;
    cout << "QCD in Signal sample, true, [0.2,1.2] = " << nQCDInSignal_true_high << endl;
    cout << "QCD in Signal sample, true, [  0,1.2] = " << nQCDInSignal_true_sub  << endl;
    */

    double nQCDInSignal_true_low  = hqcd->Integral(sig_bin_low+1,sig_bin_up) ; //(0,0.1)->(1,1) 
    double nQCDInSignal_true_gap  = hqcd->Integral(sig_bin_up+1, fit_bin_low) ; //(0.1,0.2)->(2,2) = int(0.1/binw+0.1)+1, same
    double nQCDInSignal_true_high = hqcd->Integral(fit_bin_low+1,fit_bin_up) ; //(0.3,1.2)->(3,12)
    double nQCDInSignal_true_sub  = hqcd->Integral(sig_bin_low+1,fit_bin_up) ; //(0.1,1.2)->(1,2)
    
    printf("True nQCD in Signal sample (sig) [ %.1f, %.1f ] = %6g\n", sig_from, sig_upto, nQCDInSignal_true_low );
    printf("True nQCD in Signal sample (gap) [ %.1f, %.1f ] = %6g\n", sig_upto, fit_from, nQCDInSignal_true_gap );
    printf("True nQCD in Signal sample (fit) [ %.1f, %.1f ] = %6g\n", fit_from, fit_upto, nQCDInSignal_true_high);
    printf("True nQCD in Signal sample (sub) [ %.1f, %.1f ] = %6g\n", sig_from, fit_upto, nQCDInSignal_true_sub );

    
    nQCDInSignal_true[j] = hqcd->Integral( sig_bin_low+1, sig_bin_up );
    cout << "\n\n AMITUOFO \n\n"<< endl;
    cout <<  " nQCDInSignal_true[j]    =   " <<   nQCDInSignal_true[j] << endl;
    cout << "\n\n AMITUOFO \n\n"<< endl;

    /*
    c1.cd(10);
    all->GetXaxis()->SetRangeUser(0,fit_upto-tiny);
    hqcd->GetXaxis()->SetRangeUser(0,fit_upto-tiny);
    all->DrawCopy("h");
    hqcd->DrawCopy("hist same");
    */
    
    // Prepare data in *fit range*
    all->GetXaxis()->SetRangeUser(fit_from,fit_upto);
    hqcd->GetXaxis()->SetRangeUser(fit_from,fit_upto);
    cout << "\n Prepare RooDataHist data and qcd\n"<< endl;
    //    RooRealVar reliso_0_to_16();
    RooDataHist data("data","data in signal sample", reliso, all);
    RooDataHist qcd("qcd","QCD in signal sample", reliso, hqcd);
    
    
    cout << "\n Print entry in RooDataHist data"<< endl;
    double sume = data.sumEntries("reliso<0.1");
    cout << "sumEntry(reliso<0.1) = " << sume << endl;
    //  cout << " 'data': nEvent in signal sample [  0 - 0.1]:  " << data.createIntegral(reliso, NormSet(reliso), Range(0,0.1)))->getVal() << endl;
    //  cout << " 'data': nEvent in signal sample [  0 - 0.1]:  " << data->createIntegral(0,0.1)->getVal() << endl;
    
    
    // Binned ML fit
    //fit to 1-jet range 0.2 to 1.2
    //  pdf_qcd.fitTo( data, Range( 0.2, 1.2 ) );y

    // model.fitTo( data, Range("fit")); //for ROOT 5.22

    model.fitTo( data, Range("fit"), SumW2Error(kFALSE) ); //for ROOT 5.26

    
    /*
    // chi2 fit
    cout << "\n Performe chi2 fit" << endl;
    RooChi2Var chi2("chi2","chi2",model,data,Extended(),Range("fit"));
    RooMinuit m(chi2);
    m.fit("mhv");
    */
    /*
    c1.cd(7);
    cout << "pad 7" << endl;
    cout << "\n model.createHistogram\n"<< endl;
    TH1F* thefit = (TH1F*)model.createHistogram("thefit",reliso,Binning(12));
    thefit->ls();
    thefit->SetMarkerSize(2);
    thefit->Draw("hist text");
    */

    cout << "amtb 3" << endl;
    nqcd.Print();
    



    c1.cd(j+1);

    // NB: when drawing QCD after data, the fit is wrong! The fit is to QCD, not data!
    
    double highest = 0;
    for(short i=3; i<10; i++){
      double this_bin = all->GetBinContent(i);
      if(this_bin > highest) highest = this_bin;
    }
    double yup = highest * 2.3;
    /*
    double yup = 3000;
    if(njetA[j]=="2j")   yup = 700;
    if(njetA[j]=="3j")   yup = 150;
    if(njetA[j]=="4mj")  yup = 30;
    */


    all->SetStats(0);
    all->GetYaxis()->SetRangeUser(0,yup);
    all->GetXaxis()->SetRangeUser(0,1.5);
    all->SetFillColor(kGray);
    all->SetLineColor(kGray);
    all->SetTitle(Form("Binned ML Fit: %s",njlabel[j]));
    all->SetXTitle("RelIso");
    all->Draw("hist");

    RooPlot *isoF4 = reliso.frame();
    cout << "amtb 4" << endl;
    qcd.plotOn(isoF4, 
	       DataError(etype), 
	       MarkerColor(kBlue) );
    cout << "amtb 5" << endl;
    data.plotOn(isoF4, DataError(etype));
    
    cout << "amtb 6" << endl;
    //  data.plotOn(isoF4, DataError(etype));
    //cout << "amtb 6b" << endl;
    const double nfit = nqcd.getVal();
    const double nfiterr = nqcd.getError();
    cout << "\n nfit = " << nfit << "\n"<< endl;
    cout << "\n nfiterr = " << nfiterr << "\n"<< endl;

    model.plotOn(isoF4, 
		 Range("extrapolate"), 
		 LineStyle(kDashed),
		 LineWidth(3),
		 Normalization(nfit, RooAbsReal::NumEvent));
    
    cout << "amtb 7" << endl;
    model.plotOn(isoF4, 
		 Range("fit"), 
		 LineColor(kRed),
		 LineWidth(3),
		 Normalization(nfit, RooAbsReal::NumEvent));	       
    /*
      model.paramOn(isoF4, 
      Label("fit result"),
      Layout(0.5,0.98,0.93));
    */
    
    if(!do_closure_test)  isoF4->SetAxisRange(0,yup,"y");
    isoF4->SetTitleOffset(1.25,"Y");
    isoF4->Draw("same");
    
  
    //-----------------------------
    // Results
    //-----------------------------
    //c1.cd(12);
    
    
    double frac_QCD_allIso = m0;
    double frac_QCD_lowIso = m1;
    double frac_QCD_gapIso = m2;
    double frac_QCD_highIso = m3;
    
    /*
    TH1F *hfit = (TH1F*)model.createHistogram("hfit",reliso);
    
    //  hfit->Scale( nfit/hfit->Integral(fit_bin_low+1, fit_bin_up) );
    double scale_factor = nfit/ frac_QCD_highIso ; //use model.createIntegral(0
    cout << "scale_factor = "<< scale_factor<< endl;
    hfit->Scale( scale_factor );
    hfit->GetYaxis()->SetRangeUser(0,yup);
    
    //  hfit->SetLineColor(kBlue);
    //  hfit->SetLineStyle(kDashed);
    hfit->SetMarkerSize(2);
    hfit->DrawCopy("hist text");
    cout << "\nhfit: \n";
    cout << " bin #1:  "<< hfit->GetBinContent(1) << endl;
    cout << " bin #2:  "<< hfit->GetBinContent(2) << endl;
    cout << " int(bin #3,13):  "<< hfit->Integral(3,12) << endl;
    */






  /* Requie ROOT 5.25
  // cumulative distribution function
  RooAbsReal *cdf = model.createCdf(reliso);
  RooPlot *isoFF = reliso.frame();
  cdf->plotOn(isoFF,LineWidth(1));
  c1.cd(8);
  isoFF->Draw();
  */

    
    cout << endl;
    cout << "QCD in Signal sample, true, [all]     = " << nQCDInSignal_true_all  << endl;
    cout << "QCD in Signal sample, true, [  0,1.2] = " << nQCDInSignal_true_sub  << endl;
    cout << "QCD in Signal sample, true, [  0,0.1] = " << nQCDInSignal_true_low  << endl;
    cout << "QCD in Signal sample, true, [0.1,0.2] = " << nQCDInSignal_true_gap  << endl;
    cout << "QCD in Signal sample, true, [0.2,1.2] = " << nQCDInSignal_true_high << endl;
    
    
    cout << endl << "fraction in ranges in control sample (roofit):" << endl;
    cout << "frac QCD allIso  = " << frac_QCD_allIso << endl;
    cout << "frac QCD lowIso  = " << frac_QCD_lowIso << endl;
    cout << "frac QCD gapIso  = " << frac_QCD_gapIso << endl;
    cout << "frac QCD highIso = " << frac_QCD_highIso << endl;

    /*    
    cout << endl << "fraction in ranges in control sample (root):" << endl;
    cout << "frac QCD allIso  = " << fracALL << endl;
    cout << "frac QCD lowIso  = " << fracLOW << endl;
    cout << "frac QCD gapIso  = " << fracGAP << endl;
    cout << "frac QCD highIso = " << fracHIGH << endl;
    */

    cout << "\n predicted QCD numbers (roofit):"<< endl;
    
    double estimate_0     = nfit / frac_QCD_highIso * frac_QCD_lowIso ;
    double estimate_error = nfiterr / nfit * estimate_0;

    nQCDInSignal_est[j] = estimate_0;
    nQCDInSignal_estErr[j] = estimate_error;



    cout << "\n\n AMITUOFO \n\n"<< endl;
    cout <<  " nQCDInSignal_est[j]    =   " <<   nQCDInSignal_est[j] << endl;
    cout <<  " nQCDInSignal_estErr[j]    =   " <<   nQCDInSignal_estErr[j] << endl;
    cout << "\n\n AMITUOFO \n\n"<< endl;



    cout << "\n QCD estimate [0 to 0.1] = " ;
    cout << nfit << " / " << frac_QCD_highIso << " * " << frac_QCD_lowIso ;
    cout << " = "<< estimate_0 << endl;
    cout << "   = nfit / frac_QCD_highIso * frac_QCD_lowIso " << endl;
    
    double estGAP_0 = nfit / frac_QCD_highIso * frac_QCD_gapIso ;
    cout << "\n QCD estimate [0.1 to 0.2] = ";
    cout << nfit << " / " << frac_QCD_highIso << " * " << frac_QCD_gapIso ;
    cout << " = "<< estGAP_0 << endl;
    cout << "   = nfit / frac_QCD_highIso * frac_QCD_gapIso " << endl;
    
    double estHIGH_0 = nfit;
    cout << "\n QCD estimate [0.2 to 1.2] = " ;
    cout << " = " << estHIGH_0 << endl;
    cout << "   = nfit / frac_QCD_highIso * frac_QCD_highIso " << endl;
    
    double estALL_0 = nfit / frac_QCD_highIso * frac_QCD_allIso ;
    cout << "\n QCD estimate [all] = " ;
    cout << nfit << " / " << frac_QCD_highIso << " * " << frac_QCD_allIso ;
    cout << " = "<< estALL_0 << endl;
    cout << "   = nfit / frac_QCD_highIso * frac_QCD_allIso " << endl;
    
    /*    
    cout << "\n predicted QCD numbers (root):"<< endl;
    
    double estimate = nfit / fracHIGH * fracLOW;
    cout << "\n QCD estimate [0 to 0.1] = " ;
    cout << nfit << " / " << fracHIGH << " * " << fracLOW ;
    cout << " = "<< estimate << endl;
    cout << "   = nfit / frac_QCD_highIso * frac_QCD_lowIso " << endl;
    
    double estGAP = nfit / fracHIGH * fracGAP ;
    cout << "\n QCD estimate [0.1 to 0.2] = ";
    cout << nfit << " / " << fracHIGH << " * " << fracGAP;
    cout << " = "<< estGAP << endl;
    cout << "   = nfit / frac_QCD_highIso * frac_QCD_gapIso " << endl;
    
    double estHIGH = nfit;
    cout << "\n QCD estimate [0.2 to 1.2] = " ;
    cout << nfit << " = "<< estHIGH << endl;
    cout << "   = nfit / frac_QCD_highIso * frac_QCD_highIso " << endl;
    
    double estALL = nfit / fracHIGH * fracALL;
    cout << "\n QCD estimate [all] = " ;
    cout << nfit << " / " << fracHIGH << " * " << fracALL ;
    cout << " = "<< estALL << endl;
    cout << "   = nfit / frac_QCD_highIso * frac_QCD_allIso " << endl;
    



    double estHIGH = nfit;

    double estALL = estHIGH / fracHIGH * fracALL ;
  
    cout << "\n--------------------------------------------------"<< endl;
    cout << "  Signal sample: R "<< fit_from << "-" << fit_upto << "  (" << njet << ")" << endl;
    cout << "--------------------------------------------------"<< endl;
    printf("Iso region      Fraction    Estimated QCD    True QCD        Dev\n");
    printf(" Low (sig)   %8.3f  %12.2f  %12.2f  %12.2f\n", fracLOW,  estimate,  nQCDInSignal_true_low, 
	   estimate/nQCDInSignal_true_low-1);
    printf(" Gap         %8.3f  %12.2f  %12.2f  %12.2f\n", fracGAP,  estGAP  ,  nQCDInSignal_true_gap, 
	 estGAP/nQCDInSignal_true_gap-1);
    printf(" High (fit)  %8.3f  %12.2f  %12.2f  %12.2f\n", fracHIGH, estHIGH ,  nQCDInSignal_true_high,
	   estHIGH/nQCDInSignal_true_high-1);
    printf(" Total       %8.3f  %12.2f  %12.2f  %12.2f\n", fracALL,  estALL  ,  nQCDInSignal_true_sub, 
	   estALL/nQCDInSignal_true_sub-1);
    cout << "--------------------------------------------------"<< endl;
    */
    
    cout << "\n Using fraction from integral of RooFit's model.createIntegral() method" <<endl;
    cout << "--------------------------------------------------"<< endl;
    cout << "  Signal sample: R "<< fit_from << "-" << fit_upto << "  (" << njet << ")" << endl;
    cout << "--------------------------------------------------"<< endl;
    printf("Iso region      Fraction    Estimated QCD    True QCD        Dev\n");
    printf(" Low (sig)   %8.3f  %12.2f  %12.2f  %12.2f\n", frac_QCD_lowIso,  estimate_0,  nQCDInSignal_true_low, 
	 estimate_0/nQCDInSignal_true_low-1);
    printf(" Gap         %8.3f  %12.2f  %12.2f  %12.2f\n", frac_QCD_gapIso,  estGAP_0  ,  nQCDInSignal_true_gap, 
	   estGAP_0/nQCDInSignal_true_gap-1);
    printf(" High (fit)  %8.3f  %12.2f  %12.2f  %12.2f\n", frac_QCD_highIso, estHIGH_0 ,  nQCDInSignal_true_high,
	 estHIGH_0/nQCDInSignal_true_high-1);
    printf(" Total       %8.3f  %12.2f  %12.2f  %12.2f\n", frac_QCD_allIso,  estALL_0  ,  nQCDInSignal_true_sub, 
	   estALL_0/nQCDInSignal_true_sub-1);
    cout << "--------------------------------------------------"<< endl;
    

    
    Deviation[j] = nQCDInSignal_est[j]/nQCDInSignal_true[j] -1;

    // Show results in plot
    TPaveText *res= new TPaveText(0.55,0.6,0.9,0.9,"NDC");
    res->SetBorderSize(1);      
    res->SetFillColor(0);
    res->AddText(Form("Range (%.1f,%.1f)", fit_from,fit_upto));
    res->AddText(Form("True = %.1f",   nQCDInSignal_true[j] ));
    res->AddText(Form("Est = %.1f #pm %.1f",    nQCDInSignal_est[j],   nQCDInSignal_estErr[j]  ));
    res->AddText(Form("#Delta = %.2f", Deviation[j]         ));
    res->Draw();

  }

  printf("\n---------------------------------------------------------------------\n");
  cout << "  Range  (" << fit_from << "," << fit_upto <<")";
  printf("\n---------------------------------------------------------------------\n");
  printf("         njet        True         Est (+/-fit error)        Dev\n");
  for(int jj=0; jj<4; jj++){
    printf(" %12s  %10.2f  %10.2f (+/-%6.2f)  %10.2f\n", 
	   njlabel[jj], nQCDInSignal_true[jj], nQCDInSignal_est[jj], nQCDInSignal_estErr[jj], Deviation[jj] );   
  }
  printf("-----------------------------------------------------------------------\n");

  
  int ntemp = 1;
  if(use_one_template==false) ntemp = 4;

  // out1: "B3_e20_4T" or "B3_e20_1T_1j"
  string out = Form("%s_4T", control_sample.c_str());
  if(ntemp==1) out = Form("%s_1T_%s", control_sample.c_str(), templateBin.c_str());

  // for plots
  string outplot = Form("MLfit__%s__r%.1fto%.1f", 
			out.c_str(), fit_from, fit_upto);
			
  /*
  if(ntemp==1) {
    out = Form("MLfit__%s_1T_%s__r%.1fto%.1f", 
	       control_sample.c_str(),
	       templateBin.c_str(), 
	       fit_from, fit_upto);
  }
  */
  c1.SaveAs( Form("%s.gif", outplot.c_str() ) );

  //c1.SaveAs( Form("%s.pdf",out.c_str() ) );
  //gROOT->ProcessLine(Form(".!ps2pdf -dEPSCrop %s", out.c_str())); 


  // Write out results to a text file  
  //cout << "\n (c) Write out results text file"<< endl;
  
  ofstream myfile;
  myfile.open( Form( "est_template_%s.txt", out.c_str() ),
	       ios::app ); //out:(overwrite); ios:app (append)
  myfile.setf(ios::fixed,ios::floatfield);
  myfile << endl;
  for(int k=0; k<4; ++k){ // write out 1,2,3,4j only
    myfile << Deviation[k] << endl;
  }
  myfile.close();

}
	      
