double fit_function(double *x, double *par)
{
  //double par1 = 330350.4; // halflife in seconds
  double par1 = 91.764; // halflife in hours

  double val = par[0]*TMath::Exp(-1*TMath::Log(2)/par1*x[0]) + par[1];

  return val;
}

void RadonPlot()
{
  gStyle->SetOptTitle(0);

  TTree* tree = new TTree;
  //tree->ReadFile("RadonData.dat");
  //tree->ReadFile("RadonData23.dat");
  //tree->ReadFile("RadonData_000.dat");
  //tree->ReadFile("RadonDataReprocessed.dat");
  tree->ReadFile("RadonDataReprocessed111116.dat");

  double time;
  double timeErr;
  double A;
  double AErr;
  double N;
  double NErr;

  tree->SetBranchAddress("time",&time);
  tree->SetBranchAddress("timeErr",&timeErr);
  tree->SetBranchAddress("A",&A);
  tree->SetBranchAddress("AErr",&AErr);
  tree->SetBranchAddress("N",&N);
  tree->SetBranchAddress("NErr",&NErr);

  int nData = tree->GetEntries();

  double *t = new double[nData];
  double *t_err = new double[nData];
  double *A_Bi = new double[nData];
  double *A_Bi_err = new double[nData];
  double *N_Rn = new double[nData];
  double *N_Rn_err = new double[nData];

  double *N_Rn_copy = new double[nData];
  double *N_Rn_err_copy = new double[nData];
  double *t_copy = new double[nData];

  double t_Unix = new double[nData];

  double integral = 0.0;
  double integralErr = 0.0;
  double RunTime = 0.0;
  int k = 0;
  for (const int i = 0; i < nData; i++) {
     tree->GetEntry(i);

     t[i] = time;
     t_err[i] = timeErr;
     A_Bi[i] = A;
     A_Bi_err[i] = AErr;
     N_Rn[i] = N;
     N_Rn_err[i] = NErr;
     t_Unix[i] = time * 3600.0 + 1304146800.0;

     if (t_Unix[i] >= 1305529200.0) {N_Rn_copy[k] = N; N_Rn_err_copy[k] = N_Rn_err[i]; t_copy[k] = t_Unix[i]; k++;}

     //if (i > 11 && i < 43) {integral += N*(1 - TMath::Exp(-2.0*timeErr/91.764)); integralErr = NErr*NErr*(1 - TMath::Exp(-2.0*timeErr/91.764))*(1 - TMath::Exp(-2.0*timeErr/91.764));}
     integral += N*(1 - TMath::Exp(-2.0*timeErr/132.387));
     integralErr += NErr*NErr*(1 - TMath::Exp(-2.0*timeErr/132.387))*(1 - TMath::Exp(-2.0*timeErr/132.387));
     RunTime += 2.0*timeErr;
  }

  cout << "Mean = " << TMath::Mean(k,N_Rn_copy) << " +- " << TMath::RMS(k,N_Rn_copy) << endl;

  TGraphErrors *grA = new TGraphErrors(nData, t, A_Bi, 0, A_Bi_err);
  TGraphErrors *grN = new TGraphErrors(nData, t, N_Rn, 0, N_Rn_err);
  TGraphErrors *grN_copy = new TGraphErrors(k, t_copy, N_Rn_copy, 0, N_Rn_err_copy);

  grA->SetTitle("^{214}Bi activity");
  grN->SetTitle("^{222}Rn atoms in xenon");
  grN_copy->SetTitle("^{222}Rn atoms in xenon");

  grA->GetXaxis()->SetTitle("time [h]");
  grN->GetXaxis()->SetTitle("time [h]");
  grN_copy->GetXaxis()->SetTitle("unix time");

  grA->GetYaxis()->SetTitle("^{214}Bi decay rate [mBq]");
  grN->GetYaxis()->SetTitle("^{222}Rn atoms in Xenon");
  grN_copy->GetYaxis()->SetTitle("^{222}Rn atoms in Xenon");

  grN->GetYaxis()->SetTitleOffset(1.3);
  grN_copy->GetYaxis()->SetTitleOffset(1.3);

  grA->SetMarkerStyle(20);
  grN->SetMarkerStyle(20);
  grN_copy->SetMarkerStyle(20);

  grA->SetMarkerSize(0.8);
  grN->SetMarkerSize(0.8);
  grN_copy->SetMarkerSize(0.8);

  //grA->SetMarkerColor(kRed);
  //grN->SetMarkerColor(kRed);
  //grN_copy->SetMarkerColor(kRed);

  TF1 *fit1 = new TF1("fit1",fit_function,0,3100,2);
  TF1 *fit2 = new TF1("fit2",fit_function,0,3100,2);
  TF1 *fit3 = new TF1("fit3",fit_function,0,200,2);

  grA->Fit("fit1","rn");
  grN->Fit("fit2","r");
  grN_copy->Fit("fit3","rn");


  double *par1 = fit1->GetParameters();
  double *par2 = fit2->GetParameters();
  double *par3 = fit3->GetParameters();

  fit1->SetParameters(par1);
  fit2->SetParameters(par2);
  fit3->SetParameters(par3);

  fit1->SetLineColor(kBlack);
  fit2->SetLineColor(kBlack);
  fit3->SetLineColor(kBlack);

  fit1->SetLineWidth(1);
  fit2->SetLineWidth(1);
  fit3->SetLineWidth(1);

  // lines indicating pump and purifier issues
  TLine *l1 = new TLine(2445,0,2445,1000);
  TLine *l2 = new TLine(2487,0,2487,1000);
  TLine *l3 = new TLine(2631,0,2631,1000);
  TLine *l4 = new TLine(2678,0,2678,1000);

  l1->SetLineWidth(1);
  l1->SetLineStyle(2);

  l2->SetLineWidth(1);
  l2->SetLineStyle(2);

  l3->SetLineWidth(1);
  l3->SetLineStyle(2);

  l4->SetLineWidth(1);
  l4->SetLineStyle(2);

  TText *txt1 = new TText(2400,1020,"pump off (1 SLPM)");
  TText *txt2 = new TText(2440,1070,"purifier 1 isolated");
  TText *txt3 = new TText(2590,1020,"pump on (8 SLPM)");
  TText *txt4 = new TText(2640,1070,"purifier 2 isolated");

  txt1->SetTextFont(2);
  txt1->SetTextSize(0.02);

  txt2->SetTextFont(2);
  txt2->SetTextSize(0.02);

  txt3->SetTextFont(2);
  txt3->SetTextSize(0.02);

  txt4->SetTextFont(2);
  txt4->SetTextSize(0.02);

  cout << fit2->Integral(528,1318) / 25.5 << " +- " << fit2->IntegralError(528,1318,fit2->GetParameters()) / 25.5 << endl;
  cout << "Integral: " << integral << " +- " << TMath::Sqrt(integralErr) << endl;
  cout << "Run time: " << RunTime << " h" << endl;

  
  grN->GetYaxis()->SetRangeUser(10,8000);
  
  //TGaxis *axis = new TGaxis(grN->GetXaxis()->GetXmax(),grN->GetYaxis()->GetXmin(),grN->GetXaxis()->GetXmax(),grN->GetYaxis()->GetXmax(),0,grA->GetYaxis()->GetXmax(),510,"+L");
  //TGaxis *axis = new TGaxis(grN->GetXaxis()->GetXmax(),grN->GetYaxis()->GetXmin(),grN->GetXaxis()->GetXmax(),grN->GetYaxis()->GetXmax(),0.08393,12.5893,510,"+G");
  TGaxis *axis = new TGaxis(grN->GetXaxis()->GetXmax(),10,grN->GetXaxis()->GetXmax(),8000,0.17115,136.9207,510,"+G");

  axis->SetTitle("^{222}Rn activity (#muBq/kg)");
  //axis->SetTitleFont(82);
  axis->SetTitleSize(0.04);
  //axis->SetLabelFont(82);
  axis->SetLabelSize(0.04);
  axis->SetLabelOffset(0.04);
  axis->CenterTitle(true);

  //TFile *fOut = new TFile("RadonPlot111116.root","RECREATE");
  //grN->Write();
  //fOut->Close();

  TCanvas *c1 = new TCanvas();
  grA->Draw("AZP");
  fit1->Draw("same");

  TCanvas *c2 = new TCanvas();
  grN->Draw("AZP");
  fit2->Draw("same");
  axis->Draw("same");
  c2->SetLogy(1);
/*  l1->Draw("same");
  l2->Draw("same");
  l3->Draw("same");
  l4->Draw("same");
  txt1->Draw("same");
  txt2->Draw("same");
  txt3->Draw("same");
  txt4->Draw("same");
*/
  TCanvas *c3 = new TCanvas();
  grN_copy->Draw("AZP");

  grN_copy->GetXaxis()->SetLimits(1305529200.0, 1317106800.0);
  c3->Update();

  return;
}
