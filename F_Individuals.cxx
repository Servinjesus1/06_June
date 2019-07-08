//Macro 5: Fit all four peaks individually in their regions
//DSAM Method Peak, ExModGaus peaks too
//Date: 06/07/2019
//Author: Spencer Fretwell

{
  // #include "A1_Inclusions.cxx+"
  #include "A1_Inclusions.cxx+"

  //ROOT Variables
  // gROOT->SetBatch(kTRUE);
  gStyle->SetOptFit(1);
  gStyle->SetStatW(0.08);
  gStyle->SetStatX(1);
  gStyle->SetStatY(1);
  gStyle->SetOptStat(0);

  //Initial Variables___________________________________________________________

  //Fit Range
  r1 = 10;
  r2 = 500;

  //Debug
  double p[10];

  //Files
  const char* f_in = "Run3.root";
  const char* f_out = "./Output/F_Individuals/Run3_Fit.root";

  //Setup TObjects______________________________________________________________
  TF1* FT; //Total Function, sum of Fits[]
  TFitResultPtr result; //Fit Results
  TCanvas* CF = new TCanvas("Fit","Fit Result",1920,1080); //Fit Canvas
  map<int,vector<int>> parVec = FTParVec();

  //Files
  TFile* FI = new TFile(f_in);
  TFile* FO = new TFile(f_out,"RECREATE"); FO->cd();

  //Folders
  TDirectory* Rts = FO->mkdir("Results");
  TDirectory* Hgs = FO->mkdir("Histograms");
  TDirectory* Fns = FO->mkdir("Functions");
  TDirectory* Cvs = FO->mkdir("Canvases");

  //Histograms
  TH1F* H1 = (TH1F*)FI->Get(FI->GetListOfKeys()->At(1)->GetName()); //Original H
  TH1F* R1; TH1F* R2;

  //Rebin H1____________________________________________________________________
  //First Guess: 0.2eV bins
  WidthRebin(H1,0.2);

  //Find Peaks
  H1->SetAxisRange(25,160);
  TSpectrum* TS1 = new TSpectrum(4);
  TS1->Search(H1,20,"nodraw",0.0005);
  H1->SetAxisRange(r1,r2);

  //Get Peak Centroids
  Double_t peaksX[4],peaksY[4];
  TPolyMarker* PKs = (TPolyMarker*)H1->GetListOfFunctions()->FindObject("TPolyMarker");
  for(int i=0;i<4;i++){
    peaksX[i]=PKs->GetX()[i]; //Peak X values
    peaksY[i]=H1->GetBinContent(peaksX[i]); //Peak Y values
  }
  H1->GetListOfFunctions()->Remove(PKs); //Remove Search Carats

  //Second Guess: ((Peak 1 -> 111.58eV),(Peak 3 -> 56.83eV))
  LinRebin(H1,{{peaksX[0],111.58},{peaksX[2],56.83}});

  //Final Tweaks________________________________________________________________
  //Set H1 Errors
  H1->Sumw2();
  H1->Scale(1., "width"); //Scale? From fitNormSum tutorial

  //Specify Minimizer
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1e-4);
  ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-4);
  ROOT::Math::IntegratorOneDimOptions::SetDefaultNPoints(10000);
  ROOT::Math::IntegratorOneDimOptions::SetDefaultWKSize(10000);

  //Set H1 Name/Title
  H1->SetNameTitle("H1","Consecutive Fits, Full");

  //Find Peaks
  H1->SetAxisRange(r1,160);
  TS1 = new TSpectrum(4);
  TS1->Search(H1,20,"nodraw",0.0005);
  H1->SetAxisRange(r1,r2);

  //Get Peak Centroids
  PKs = (TPolyMarker*)H1->GetListOfFunctions()->FindObject("TPolyMarker");
  for(int i=0;i<4;i++){
    peaksX[i]=PKs->GetX()[i]; //Peak X values
    peaksY[i]=H1->GetBinContent(peaksX[i]); //Peak Y values
  }
  H1->GetListOfFunctions()->Remove(PKs); //Remove Search Carats
  Hgs->cd(); H1->Write("H1");

//_______________________________FITTING________________________________________
//______________________________________________________________________________
  //(1) Background______________________________________________________________
  Fits.push_back(new TF1(FBg(r1,r2)));  //BKG (Const1 Slope1 Const2 Slope2)
  int curfit=Fits.size()-1;

  //Exclusion regions
  FitExcl.push_back({15,40});   //BG Excl Regions
  FitExcl.push_back({42,165});  //BG Excl Regions
  FitExcl.push_back({174,350}); //BG Excl Regions

  //Fit
  result = H1->Fit(Fits[curfit],"SQ0M");
  result->Print();

  //Write to File
  H1->SetTitle("(1) Background");
  CF->SetTitle("(1) Background");
  Fits[curfit]->SetTitle("(1) Background");
  R2 = (TH1F*)H1->Clone("R1");
  CanvF(CF,Fits[curfit],H1,R2,"1_Bkg");
  Rts->cd(); result->Write("R1");
  Hgs->cd(); R2->Write("R1_Bkg");
  Fns->cd(); Fits[curfit]->Write("F1");
  Cvs->cd(); CF->Write("C1");
  curfit++;

  //(2) Fit First Peak: 7LiK Voigt______________________________________________
  r2=200;
  Fits.push_back(new TF1(FVoigt(r1,r2))); //Voigt (Constant Centroid Sigma Gamma)

  //Exclusion regions
  FitExcl.push_back({r1,peaksX[0]-1});
  FitExcl.push_back({120,r2});

  //Parameter Guesses
  Fits[curfit]->SetParameter("Centroid",peaksX[0]);

  Fits[curfit]->SetParameter("Constant",peaksY[0]);

  //Fit
  result = R2->Fit(Fits[curfit],"SQ0M");
  result->Print();

  //Write to File
  R2->SetTitle("(2) Voigt");
  CF->SetTitle("(2) Voigt");
  Fits[curfit]->SetTitle("(2) Voigt");
  R1 = (TH1F*)R2->Clone("R1");
  CanvF(CF,Fits[curfit],R1,R2,"2_Vgt");
  Rts->cd(); result->Write("R2");
  Hgs->cd(); R2->Write("R2_Vgt");
  Fns->cd(); Fits[curfit]->Write("F2");
  Cvs->cd(); CF->Write("C2");
  curfit++;

  //(3) Fit First Peak: 7LiK ExMoVoigt__________________________________________
  TF1* FEMV1; int aEMV1 = 2;  int cEMV1 = 2;
  Fits.push_back(new TF1(FExMoVoigt(aEMV1,cEMV1,FEMV1))); //EMV1 (Const R_Escape Centroid Offset Sigma Tau)
  Fits[curfit]->SetName("EMV1");

  //Exclusion regions
  FitExcl.push_back({r1-50,99});
  FitExcl.push_back({113,r2+300});

  //Parameter Guesses
  FixSame(Fits[curfit],"Constant",Fits[FitNum("Voigt")]);

  FixSame(Fits[curfit],"Centroid",Fits[FitNum("Voigt")]);

  FixSame(Fits[curfit],"Gamma",Fits[FitNum("Voigt")]);

  //Fit
  result->Print();
  result = R2->Fit(Fits[curfit],"SQ0");
  Fits[2]->SetParameter(1,Fits[2]->GetParameter(1));

  //Write to File
  R2->SetTitle("(3) ExMoVoigt");
  CF->SetTitle("(3) ExMoVoigt");
  Fits[curfit]->SetTitle("(3) ExMoVoigt");
  R1 = (TH1F*)R2->Clone("R2");
  CanvF(CF,Fits[curfit],R1,R2,"3_EMV");
  Rts->cd(); result->Write("R3");
  Hgs->cd(); R2->Write("R3_EMV");
  Fns->cd(); Fits[curfit]->Write("F3");
  Cvs->cd(); CF->Write("C3");
  curfit++;

  TVirtualPad* fpad = CF->GetPad(1);
  parVec = FTParVec();
  FT = new TF1(FTotal(r1,300,parVec));
  fpad->cd(); FT->Draw("Same");
  CF->SaveAs("Output/F_Individuals/TS3.pdf");

  //(4) Fit Extra Peak: Gamma EMG_______________________________________________
  TF1* FEMG0; int aEMG0 = 1; int cEMG0 = 1;
  Fits.push_back(new TF1(FExMoGaus(aEMG0,cEMG0,FEMG0))); //EMG0 (Constant Centroid Offset Sigma Tau)

  //Exclusion regions
  FitExcl.push_back({r1,118});
  FitExcl.push_back({155,r2});

  //Parameter Guesses
  Fits[curfit]->SetParameter("Centroid",125); //Centroid Guess

  Fits[curfit]->SetParameter("Sigma",1.2); //Sigma Guess

  Fits[curfit]->SetParameter("Tau",-5); //Tau Guess
  Fits[curfit]->SetParLimits(Fits[curfit]->GetParNumber("Tau"),-10,0); //Tau Is Negative

  //Fit
  result = R2->Fit(Fits[curfit],"SQ0M");
  result->Print();

  //Write to File
  R2->SetTitle("(4) EMG");
  CF->SetTitle("(4) EMG");
  Fits[curfit]->SetTitle("(4) EMG");
  R1 = (TH1F*)R2->Clone("R1");
  CanvF(CF,Fits[curfit],R1,R2,"4_EMG");
  Rts->cd(); result->Write("48");
  Hgs->cd(); R2->Write("R4_EMG");
  Fns->cd(); Fits[curfit]->Write("F4");
  Cvs->cd(); CF->Write("C4");
  curfit++;

  fpad = CF->GetPad(1);
  parVec = FTParVec();
  FT = new TF1(FTotal(r1,300,parVec));
  fpad->cd(); FT->Draw("Same");
  CF->SaveAs("Output/F_Individuals/TS4.pdf");

  //(5) Fit Third Peak: L-G Gaus________________________________________________
  Fits.push_back(new TF1(FGaus(r1,r2))); //Create Gauss Fn (Constant R_LKCapt Centroid Sigma)

  //Exclusion regions
  FitExcl.push_back({r1,45});
  FitExcl.push_back({62,r2});

  //Parameter Guesses
  FixSame(Fits[curfit],"Constant",Fits[FitNum("Voigt")]);

  Fits[curfit]->SetParameter("Centroid",peaksX[2]);

  // FixSame(Fits[curfit],"Sigma",Fits[FitNum("Voigt")]);

  //Fit
  result = R2->Fit(Fits[curfit],"SQ0M");
  result->Print();

  //Write to File
  R2->SetTitle("(5) Gaus");
  CF->SetTitle("(5) Gaus");
  Fits[curfit]->SetTitle("(5) Gaus");
  R1 = (TH1F*)R2->Clone("R1");
  CanvF(CF,Fits[curfit],R1,R2,"5_Gaus");
  Rts->cd(); result->Write("R5");
  Hgs->cd(); R2->Write("R5_Gaus");
  Fns->cd(); Fits[curfit]->Write("F5");
  Cvs->cd(); CF->Write("C5");
  curfit++;

  fpad = CF->GetPad(1);
  parVec = FTParVec();
  FT = new TF1(FTotal(r1,300,parVec));
  fpad->cd(); FT->Draw("Same");
  CF->SaveAs("Output/F_Individuals/TS5.pdf");

  //(6) Fit Fourth Peak: L-X DSAMG______________________________________________
  TF1* DSAMG; int aDSAMG = 3;
  Fits.push_back(new TF1(FDSAMG(aDSAMG,DSAMG))); //DSAMG (Constant B_Ratio R_LKCapt Decay Energy AOffset Sigma)

  //Exclusion regions
  // FitExcl.push_back({r1,15});
  FitExcl.push_back({45,r2});

  //Parameter Guesses
  FixSame(Fits[curfit],"Constant",Fits[FitNum("Voigt")]);

  Fits[curfit]->SetParName(2,"R_LKCapt");
  Fits[curfit]->SetParLimits(2,0,0.06);
  // FixSame(Fits[curfit],"R_LKCapt",Fits[FitNum("LK_Gaussian")]);

  FixName(Fits[curfit],"AOffset",0); //Offset = 0

  // FixSame(Fits[curfit],"Sigma",Fits[FitNum("Voigt")]);

  //Fit
  result = R2->Fit(Fits[curfit],"SQ0M");
  result->Print();

  //Write to File
  R2->SetTitle("(6) DSAMG");
  CF->SetTitle("(6) DSAMG");
  Fits[curfit]->SetTitle("(6) DSAMG");
  R1 = (TH1F*)R2->Clone("R1");
  CanvF(CF,Fits[curfit],R1,R2,"6_DSAMG");
  Rts->cd(); result->Write("R6");
  Hgs->cd(); R2->Write("R6_DSAMG");
  Fns->cd(); Fits[curfit]->Write("F6");
  Cvs->cd(); CF->Write("C6");
  curfit++;

  fpad = CF->GetPad(1);
  parVec = FTParVec();
  FT = new TF1(FTotal(r1,300,parVec));
  fpad->cd(); FT->Draw("Same");
  CF->SaveAs("Output/F_Individuals/TS6.pdf");

  //(7) Fit Second Peak: 7Li*K DSAMV____________________________________________
  TF1* DSAMV;  int aDSAMV = 2;
  Fits.push_back(new TF1(FDSAMV(aDSAMV,DSAMV))); //DSAMV (Constant B_Ratio Decay Energy AOffset Sigma Gamma)

  //Exclusion regions
  FitExcl.push_back({r1,peaksX[1]-5});
  FitExcl.push_back({99,r2});

  //Parameter Guesses
  FixSame(Fits[curfit],"Constant",Fits[FitNum("Voigt")]);

  FixSame(Fits[curfit],"B_Ratio",Fits[FitNum("DSAMG")]);

  FixSame(Fits[curfit],"Decay",Fits[FitNum("DSAMG")]);

  FixSame(Fits[curfit],"Energy",Fits[FitNum("DSAMG")]);

  FixName(Fits[curfit],"AOffset",54.75); //Offset = Auger

  // FixSame(Fits[curfit],"Sigma",Fits[FitNum("Voigt")]);

  FixSame(Fits[curfit],"Gamma",Fits[FitNum("Voigt")]);

  //Fit
  result = R2->Fit(Fits[curfit],"SQ0M");
  result->Print();

  // //Forced Solution
  // Fits[curfit]->SetParameter("B_Ratio",0.1048);
  // Fits[curfit]->SetParameter("Decay",1.6);
  // Fits[curfit]->SetParameter("Energy",29);

  //Write to File
  R2->SetTitle("(7) DSAMV");
  CF->SetTitle("(7) DSAMV");
  Fits[curfit]->SetTitle("(7) DSAMV");
  R1 = (TH1F*)R2->Clone("R1");
  CanvF(CF,Fits[curfit],R1,R2,"7_DSAMV");
  Rts->cd(); result->Write("R7");
  Hgs->cd(); R2->Write("R7_DSAMV");
  Fns->cd(); Fits[curfit]->Write("F7");
  Cvs->cd(); CF->Write("C7");
  curfit++;

  fpad = CF->GetPad(1);
  parVec = FTParVec();
  FT = new TF1(FTotal(r1,300,parVec));
  fpad->cd(); FT->Draw("Same");
  CF->SaveAs("Output/F_Individuals/TS7.pdf");

  //(8) Fit Second Peak: 7Li*K ExMoVoigt________________________________________
  TF1* FEMV2;  int aEMV2 = 3;  int cEMV2 = 3;
  Fits.push_back(new TF1(FExMoVoigt(aEMV2,cEMV2,FEMV2))); //EMV2 (Constant R_Escape B_Ratio Energy AOffset Offset Sigma Tau Gamma)
  Fits[curfit]->SetName("EMV2");

  //Exclusion regions
  FitExcl.push_back({r1-50,74});
  FitExcl.push_back({90,r2+50});

  //Parameter Guesses
  FixSame(Fits[curfit],"Constant",Fits[FitNum("Voigt")]);

  SetSame(Fits[curfit],"R_Escape",Fits[FitNum("EMV1")]);
  ParLimit(Fits[curfit],1,0.1);

  Fits[curfit]->SetParName(2,"B_Ratio");
  FixSame(Fits[curfit],"B_Ratio",Fits[FitNum("DSAMV")]);

  Fits[curfit]->SetParName(3,"Energy");
  FixSame(Fits[curfit],"Energy",Fits[FitNum("DSAMV")]);

  Fits[curfit]->SetParName(4,"AOffset");
  FixSame(Fits[curfit],"AOffset",Fits[FitNum("DSAMV")]);

  SetSame(Fits[curfit],"Offset",Fits[FitNum("EMV1")]);
  ParLimit(Fits[curfit],5,0.2);

  // FixSame(Fits[curfit],"Sigma",Fits[FitNum("EMV1")]);

  FixSame(Fits[curfit],"Tau",Fits[FitNum("EMV1")]);

  FixSame(Fits[curfit],"Gamma",Fits[FitNum("Voigt")]);

  // Fit
  result = R2->Fit(Fits[curfit],"SQ0M");
  result->Print();

  //Write to File
  R2->SetTitle("(8) ExMoVoigt");
  CF->SetTitle("(8) ExMoVoigt");
  Fits[curfit]->SetTitle("(8) ExMoVoigt");
  R1 = (TH1F*)R2->Clone("R1");
  CanvF(CF,Fits[curfit],R1,R2,"8_EMV");
  Rts->cd(); result->Write("R8");
  Hgs->cd(); R2->Write("R8_EMV");
  Fns->cd(); Fits[curfit]->Write("F8");
  Cvs->cd(); CF->Write("C8");
  curfit++;

  fpad = CF->GetPad(1);
  parVec = FTParVec();
  FT = new TF1(FTotal(r1,300,parVec));
  fpad->cd(); FT->Draw("Same");
  CF->SaveAs("Output/F_Individuals/TS8.pdf");

  //FT__________________________________________________________________________
  //Recreate FT
  parVec = FTParVec();
  FT = new TF1(FTotal(r1,300,parVec));

  //Write to File
  H1->SetTitle("Total");
  H1->SetStats(0);
  CF->SetTitle("Total Function");
  R2 = (TH1F*)H1->Clone("R2");
  FT->SetTitle("Total Function");
  CanvF(CF,FT,H1,R2,"Total");
  Fns->cd(); FT->Write("Tot");
  Cvs->cd(); CF->Write("Tot");
  FO->Close();
  /*
*/
}
