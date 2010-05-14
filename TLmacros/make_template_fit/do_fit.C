//
// Try fit ranges of:
//
//  0.2  to  1.0, 1.2, 1.4
//  0.3  to  1.1, 1.3, 1.5
//  0.4  to  1.2, 1.4, 1.6
//---------------------------

double begin[9] = { 0.2, 0.2, 0.2,
		    0.3, 0.3, 0.3,
		    0.4, 0.4, 0.4 };
double end[9] = {  1.0, 1.2, 1.4,
		   1.1, 1.3, 1.5,
		   1.2, 1.4, 1.6 };

void do_fit(string sel="usual"){


  gROOT->ProcessLine(".!rm -f est_template*.txt");

  if(sel=="usual") do_fit_usual();
  if(sel=="dijet") do_fit_dijet();

}

// template from usual s+b (em-enrich + bce + others)
//--------------------------
void do_fit_usual(){

  cout << "---------------"<< endl;
  cout << " Template Fit"<< endl;
  cout << "---------------\n"<< endl;
  gROOT->LoadMacro("fit_template_new.C");


  // 1 template
  //-------------
  for(short i=0; i<9; ++i)  fit_template_1T( "B3_e20", "1j",  begin[i], end[i] );
  for(short i=0; i<9; ++i)  fit_template_1T( "B3_e20", "2j",  begin[i], end[i] );
  for(short i=0; i<9; ++i)  fit_template_1T( "B3_e20", "3j",  begin[i], end[i] );
  for(short i=0; i<9; ++i)  fit_template_1T( "B3_e20", "4mj", begin[i], end[i] );

  // 4 templates
  //-------------
  for(short i=0; i<9; ++i)  fit_template_4T( "B3_e20", begin[i], end[i] );


  cout << "\n\n";
  cout << "---------------------"<< endl;
  cout << " QCD estimates plots" << endl;
  cout << "---------------------\n"<< endl;
  gROOT->LoadMacro("plot_qcd_estimate.C");

  plot_qcd_estimate("template","B3_e20",1,"1j");
  plot_qcd_estimate("template","B3_e20",1,"2j");
  plot_qcd_estimate("template","B3_e20",1,"3j");
  plot_qcd_estimate("template","B3_e20",1,"4mj");
  plot_qcd_estimate("template","B3_e20",4,"");

}

// Use template from dijet pt15 only
//-------------------------
void do_fit_dijet(){


  cout << "---------------"<< endl;
  cout << " Template Fit"<< endl;
  cout << "---------------\n"<< endl;
  gROOT->LoadMacro("fit_template_new.C");

  // 1 template (1j) dijet
  //-------------
  for(short i=0; i<9; ++i)  fit_template_1T_dijet( "B3_e20", "1j",  begin[i], end[i] );

  cout << "\n\n";
  cout << "---------------------"<< endl;
  cout << " QCD estimates plots" << endl;
  cout << "---------------------\n"<< endl;
  gROOT->LoadMacro("plot_qcd_estimate.C");

  plot_qcd_estimate("template","B3_e20",1,"1j");

  
}
