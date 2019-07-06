//Macro 5: Using previous fits
//Date: 06/07/2019
//Author: Spencer Fretwell

{
  #include <A1_Inclusions.cxx+>
  //ROOT Variables
  // gROOT->SetBatch(kTRUE);

  //Initial Variables___________________________________________________________

  //Fit Range
  r1 = 10;
  r2 = 900;

  //Files
  const char* f_in = "Run3.root";
  const char* f_prev = "./Output/F_Individuals/Run3_Fit.root";
  const char* f_out = "./Output/G_Final/Run3_Fit.root";

  //Setup TObjects______________________________________________________________
  TF1* FT = new TF1(); //Total Function, sum of Fits[]
  TFitResultPtr result; //Fit Results
  TCanvas* CF = new TCanvas("Fit","Fit Result",1920,1080); //Fit Canvas
  CF->SetLogy();

  //Files
  TFile* FI = new TFile(f_in);
  TFile* FO = new TFile(f_out,"RECREATE");
  TFile* FP = new TFile(f_prev,"READ");

  //Folders
  TDirectory* Rts = FO->mkdir("Results");
  TDirectory* Hgs = FO->mkdir("Histograms");
  TDirectory* Fns = FO->mkdir("Functions");
  TDirectory* Cvs = FO->mkdir("Canvases");

  //Histograms
  TH1F* H1 = (TH1F*)FP->GetDirectory("Histograms")->Get("H1");
  TH1F* R1;
  TH1F* R2;

  //Past Fits
  vector<TF1*> IndFits;
  TDirectory* PFns = FP->GetDirectory("Functions");
  IndFits.resize(8);
  for(int i=0;i<8;i++){
    IndFits[i]=(TF1*)PFns->Get(("F"+to_string(i+1)).c_str());
  }

  //Final Tweaks________________________________________________________________
  //Specify Minimizer
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1e-4);
  ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-4);
  ROOT::Math::IntegratorOneDimOptions::SetDefaultNPoints(10000);
  ROOT::Math::IntegratorOneDimOptions::SetDefaultWKSize(10000);

  //Set H1 Name/Title
  H1->SetNameTitle("H1","Final Fitting"); H1->SetAxisRange(r1,r2,"X");
  Hgs->cd(); H1->Write("H1");

//_______________________________FITTING________________________________________
//______________________________________________________________________________
  //(1) First Peak, EMG0, and Background________________________________________
  Fits.push_back(new TF1(FBg(r1,r2)));  //Create Background Fn (4 parameters)
  Fits.push_back(new TF1(FVoigt(r1,r2))); //Create Voigt Fn (4 parameters)
  TF1* FEMV; int aEMV1 = 2;  int cEMV1 = 2;
  Fits.push_back(new TF1(FExMoVoigt(aEMV1,cEMV1,FEMV))); //Create EMV Fn (7 parameters)
  TF1* FEMG0; int aEMG0 = 1; int cEMG0 = 1;
  Fits.push_back(new TF1(FExMoGaus(aEMG0,cEMG0,FEMG0))); //Create EMG Fn (4 parameters)

  //Exclusion regions
  FitExcl.push_back({r1,40});
  FitExcl.push_back({42,99});
  // FitExcl.push_back({115,168}); //IF EMG0 makes things bad...
  FitExcl.push_back({177,500});

  //Previous Parameters
  FTInitPars(IndFits, 0);

  //FT
  parx[0] = {};
  parx[1] = {};
  parx[2] = {{"Constant",1},{"Centroid",1},{"Gamma",1}};
  parx[3] = {};
  delete FT; FT = new TF1(FTotal(r1,r2));

  //Fit
  result = H1->Fit(FT,"SQ0R");
  ParCopy(FT);
  double p1[11]; FT->GetParameters(p1);

  //Write to File
  H1->SetTitle("(1) Peak 1 + Background + High EMG");
  CF->SetTitle("(1) Peak 1 + Background + High EMG");
  FT->SetTitle("(1) Peak 1 + Background + High EMG");
  R1 = (TH1F*)H1->Clone("R1");
  FitExcl2 = FitExcl;
  r2=220; CanvG(CF,FT,H1,R1,"1_Pk1"); r2=900;
  FitExcl = FitExcl2;
  R1 = (TH1F*)H1->Clone("R1");
  CanvG(CF,FT,H1,R1,"1_Pk1_900");
  Rts->cd(); result->Write("R1");
  Hgs->cd(); R1->Write("R1_Pk1");
  Fns->cd(); FT->Write("F1");
  Cvs->cd(); CF->Write("C1");

  //Reset Objects
  parx.clear();
  Fits.clear();
/*
  //(2) Second/Third Peak w/o Background________________________________________
  TF1* FDSAMV1;  int aDSAMV1 = 2;
  Fits.push_back(new TF1(FDSAMV(aDSAMV1,FDSAMV1))); //Create DSAMV Fn (5 parameters)
  TF1* FEMV2;  int aEMV2 = 3;  int cEMV2 = 3;
  Fits.push_back(new TF1(FExMoVoigt(aEMV2,cEMV2,FEMV2))); //Create EMV Fn (9 parameters)
  Fits.push_back(new TF1(FGaus(r1,r2))); //Create Gauss Fn (4 parameters)

  //Exclusion regions
  FitExcl.push_back({15,40});
  FitExcl.push_back({100,168});
  FitExcl.push_back({177,500});

  //Previous Parameters
  FTInitPars(IndFits, 4);

  //Parameter Guesses
  FixName(Fits[0],"Constant",FT->GetParameter("(1) Constant")); //Voigt
  FixName(Fits[0],"Sigma",   FT->GetParameter("(1) Sigma"));    //Voigt
  FixName(Fits[0],"Gamma",   FT->GetParameter("(1) Gamma"));    //Voigt
  FixName(Fits[0],"Offset",  54.75); //Offset = Auger
  FixName(Fits[1],"Constant",FT->GetParameter("(1) Constant")); //Voigt
  FixName(Fits[1],"Sigma",   FT->GetParameter("(2) Sigma"));    //EMV
  FixName(Fits[1],"Tau",     FT->GetParameter("(2) Tau"));      //EMV
  FixName(Fits[1],"Gamma",   FT->GetParameter("(1) Gamma"));    //Voigt
  FixName(Fits[2],"Constant",FT->GetParameter("(1) Constant")); //Voigt

  //FT
  parx[0] = {};
  parx[1] = {{"B_Ratio",0},{"Energy",0},{"AOffset",0}};
  parx[2] = {{"Sigma",0}};
  delete FT; FT = new TF1(FTotal(r1,r2));
  FT->Draw();

  //Fit
  result = R1->Fit(FT,"SQ0R");
  ParCopy(FT);
  double p2[18]; FT->GetParameters(p2);
  Fits.push_back(new TF1(FDSAM(aDSAMV1,FDSAMV1))); //DSAM before gaus for comparison
  {double p[10]; Fits[0]->GetParameters(p); Fits[3]->SetParameters(p);}

  //Write to File
  R1->SetTitle("(2) Peak 2/3");
  CF->SetTitle("(2) Peak 2/3");
  FT->SetTitle("(2) Peak 2/3");
  R2 = (TH1F*)R1->Clone("R1");
  FitExcl2 = FitExcl;
  r2=220; CanvG(CF,FT,R1,R2,"2_Pk2+3"); r2=900;
  FitExcl = FitExcl2;
  CanvG(CF,FT,R1,R2,"2_Pk2+3_900");
  Rts->cd(); result->Write("R2");
  Hgs->cd(); R2->Write("R2_Pk2+3");
  Fns->cd(); FT->Write("F2");
  Cvs->cd(); CF->Write("C2");

  //Reset Objects
  parx.clear();
  Fits.clear();

  //(3) Everything______________________________________________________________
  Fits.push_back(new TF1(FBg(r1,r2)));  //Create Background Fn (4 parameters)
  Fits.push_back(new TF1(FVoigt(r1,r2))); //Create Voigt Fn (4 parameters)
  TF1* FEMV; int aEMV1 = 2;  int cEMV1 = 2;
  Fits.push_back(new TF1(FExMoVoigt(aEMV1,cEMV1,FEMV))); //Create EMV Fn (6 parameters)
  TF1* FEMG0; int aEMG0 = 1; int cEMG0 = 1;
  Fits.push_back(new TF1(FExMoGaus(aEMG0,cEMG0,FEMG0))); //Create EMG Fn (4 parameters)
  TF1* FDSAMV1;  int aDSAMV1 = 2;
  Fits.push_back(new TF1(FDSAMV(aDSAMV1,FDSAMV1))); //Create DSAMV Fn (5 parameters)
  TF1* FEMV2;  int aEMV2 = 3;  int cEMV2 = 3;
  Fits.push_back(new TF1(FExMoVoigt(aEMV2,cEMV2,FEMV2))); //Create EMV Fn (9 parameters)
  Fits.push_back(new TF1(FGaus(r1,r2))); //Create Gauss Fn (4 parameters)
  TF1* FDSAMG1; int aDSAMG1 = 3;
  Fits.push_back(new TF1(FDSAMG(aDSAMG1,FDSAMG1))); //Create DSAMG Fn (6 parameters)

  //Exclusion regions
  FitExcl.push_back({177,500});

  //Previous Parameters
  vector<double> pall;
  pall.insert(pall.end(),&p1[0],&p1[13]);
  pall.insert(pall.end(),&p2[0],&p2[18]);
  {
    int j=0;
    for(int k=0;k<6;k++){
      for(int i=0;i<Fits[k]->GetNpar();i++){
        Fits[k]->SetParameter(i,pall[j]);
        j++;
      }
    }
    for(int i=0; i<Fits.size(); i++){
      for(int j=0; j<IndFits[i]->GetNpar(); j++){
        double a,b; IndFits[i]->GetParLimits(j,a,b);
        Fits[i]->SetParLimits(j,a,b);
      }
    }
  }

  //Parameter Guesses

  //FT
  delete FT; FT = new TF1(FTotal(r1,r2));
  ParCopy(FT);

  //Prelim Write to File
  H1->SetTitle("(3) Complete Fit");
  CF->SetTitle("(3) Complete Fit");
  FT->SetTitle("(3) Complete Fit");
  R1 = (TH1F*)H1->Clone("R1");
  FitExcl2 = FitExcl;
  r2=160; CanvG(CF,FT,H1,R1,"3_All_Prelim"); r2=900;
  FitExcl = FitExcl2;
  CanvG(CF,FT,H1,R1,"3_All_Prelim_900");
  FitExcl = FitExcl2;
  /*

  //Fit
  result = H1->Fit(FT,"S0REM");
  ParCopy(FT);

  //Write to File
  H1->SetTitle("(3) Complete Fit");
  CF->SetTitle("(3) Complete Fit");
  FT->SetTitle("(3) Complete Fit");
  R1 = (TH1F*)H1->Clone("R1");
  FitExcl2 = FitExcl;
  r2=160; CanvG(CF,FT,H1,R1,"3_All"); r2=900;
  FitExcl = FitExcl2;
  CanvG(CF,FT,H1,R1,"3_All_900");
  Rts->cd(); result->Write("R3");
  Hgs->cd(); R1->Write("R3_All");
  Fns->cd(); FT->Write("F3");
  Cvs->cd(); CF->Write("C3");

  FO->Close();
  */
}
