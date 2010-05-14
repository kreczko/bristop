// 24 Mar 2010
{

  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);

  TCanvas c1("c1","event table",1000,600);

  TFile f("../test_mc_mixture.root");

  TH1D *allO = f.Get("Data_njetsVcuts");
  TH1D *sigO = f.Get("Signal_njetsVcuts");
  TH1D *wjO = f.Get("Wjets_njetsVcuts");
  TH1D *zjO = f.Get("Zjets_njetsVcuts");
  TH1D *qcdO = f.Get("QCD_njetsVcuts");
  TH1D *stopO = f.Get("SingleTop_njetsVcuts");

  double origlumi = 10.0; //original int lumi assumed
  //  double scaletolumi[4] ={ 20, 10, 1, 0.1 };//scale from 20/pb to 10,1,0.1
  //  string lumi[4]={"20pb","10pb","1pb","0.1pb"};
  double scaletolumi[3] ={ 10, 1, 0.1 };//scale from 20/pb to 10,1,0.1
  string lumi[3] = {"10pb","1pb","0.1pb"};

  for(int i=0; i<=2; i++){

    TH1D *all  = allO->Clone();
    TH1D *sig  = sigO->Clone();
    TH1D *wj   = wjO->Clone();
    TH1D *zj   = zjO->Clone();
    TH1D *qcd  = qcdO->Clone();
    TH1D *stop = stopO->Clone();

    all->Scale( scaletolumi[i]/origlumi );
    sig->Scale( scaletolumi[i]/origlumi );
    wj->Scale(  scaletolumi[i]/origlumi );
    zj->Scale(  scaletolumi[i]/origlumi );
    qcd->Scale( scaletolumi[i]/origlumi );
    stop->Scale(scaletolumi[i]/origlumi );

    all ->SetTitle( Form("%s for %.1f pb^{-1}",  all->GetTitle(), scaletolumi[i]) );
    sig ->SetTitle( Form("%s for %.1f pb^{-1}",  sig->GetTitle(), scaletolumi[i]) );
    wj  ->SetTitle( Form("%s for %.1f pb^{-1}",   wj->GetTitle(), scaletolumi[i]) );
    zj  ->SetTitle( Form("%s for %.1f pb^{-1}",   zj->GetTitle(), scaletolumi[i]) );
    qcd ->SetTitle( Form("%s for %.1f pb^{-1}",  qcd->GetTitle(), scaletolumi[i]) );
    stop->SetTitle( Form("%s for %.1f pb^{-1}", stop->GetTitle(), scaletolumi[i]) );


    all->Draw();
    c1.SaveAs(Form("all_%s.gif",lumi[i]));

    sig->Draw();
    c1.SaveAs(Form("sig_%s.gif",lumi[i]));
    
    wj->Draw();
    c1.SaveAs(Form("wj_%s.gif",lumi[i]));
    
    zj->Draw();
    c1.SaveAs(Form("zj_%s.gif",lumi[i]));
   
    qcd->Draw();
    c1.SaveAs(Form("qcd_%s.gif",lumi[i]));
    
    stop->Draw();
    c1.SaveAs(Form("stop_%s.gif",lumi[i]));

  }
}
