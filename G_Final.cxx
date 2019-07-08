//Macro 5: Using previous fits
//Date: 06/07/2019
//Author: Spencer Fretwell

{
  #include <A1_Inclusions.cxx+>
  //ROOT Variables
  gROOT->SetBatch(kTRUE);

  //Initial Variables___________________________________________________________

  //Fit Range
  r1 = 10;
  r2 = 900;

  //Files
  const char* f_in = "Run3.root";
  const char* f_prev = "./Output/F_Individuals/Run3_Fit.root";
  const char* f_out = "./Output/G_Final/Run3_Fit.root";

  //Setup TObjects______________________________________________________________
  TF1* FT; //Total Function, sum of Fits[]
  TFitResultPtr result; //Fit Results
  TCanvas* CF = new TCanvas("Fit","Fit Result",1920,1080); //Fit Canvas
  CF->SetLogy();
  map<int,vector<int>> parVec;

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
  Fits.push_back(new TF1(FBg(r1,r2)));  //BKG (Const1 Slope1 Const2 Slope2)
  Fits.push_back(new TF1(FVoigt(r1,r2))); //Voigt (Constant Centroid Sigma Gamma)
  TF1* FEMV1; int aEMV1 = 2;  int cEMV1 = 2;
  Fits.push_back(new TF1(FExMoVoigt(aEMV1,cEMV1,FEMV1))); //EMV1 (Const R_Escape Centroid Offset Sigma Tau)
  Fits[2]->SetName("EMV1");
  TF1* FEMG0; int aEMG0 = 1; int cEMG0 = 1;
  Fits.push_back(new TF1(FExMoGaus(aEMG0,cEMG0,FEMG0))); //EMG0 (Constant Centroid Offset Sigma Tau)

  //Exclusion regions
  FitExcl.push_back({15,40});
  FitExcl.push_back({42,98});
  FitExcl.push_back({160,500});

  //Previous Parameters
  FTInitPars(IndFits, 0);

  //FT
  parx[0] = {};
  parx[1] = {};
  parx[2] = {{"Constant",1},{"Centroid",1},{"Gamma",1}};
  parx[3] = {};
  parVec = FTParVec();
  delete FT; FT = new TF1(  FTotal(r1,r2,parVec));
  PrintParameters(FT);

  //Fit
  result = H1->Fit(FT,"SQ0");
  vector<double> p; ParSave(p);

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
  // Rts->cd(); result->Write("R1");
  // Hgs->cd(); R1->Write("R1_Pk1");
  // Fns->cd(); FT->Write("F1");
  // Cvs->cd(); CF->Write("C1");

  //Reset Objects
  parx.clear();
  Fits.clear();

  //(2) Second/Third Peak w/o Background________________________________________
  Fits.push_back(new TF1(FGaus(r1,r2))); //Create Gauss Fn (Constant R_LKCapt Centroid Sigma)
  TF1* DSAMG; int aDSAMG = 3;
  Fits.push_back(new TF1(FDSAMG(aDSAMG,DSAMG))); //DSAMG (Constant B_Ratio R_LKCapt Decay Energy AOffset Sigma)
  TF1* DSAMV;  int aDSAMV = 2;
  Fits.push_back(new TF1(FDSAMV(aDSAMV,DSAMV))); //DSAMV (Constant B_Ratio Decay Energy AOffset Sigma Gamma)
  TF1* FEMV2;  int aEMV2 = 3;  int cEMV2 = 3;
  Fits.push_back(new TF1(FExMoVoigt(aEMV2,cEMV2,FEMV2))); //EMV2 (Constant R_Escape B_Ratio Energy AOffset Offset Sigma Tau Gamma)
  Fits[3]->SetName("EMV2");

  //Exclusion regions
  FitExcl.push_back({100,r2});

  //Previous Parameters
  FTInitPars(IndFits, 4);

  //Parameter Guesses
  FixName(Fits[0],"Constant",FT->GetParameter("(1) Constant")); //Voigt
  FixName(Fits[0],"Sigma",   FT->GetParameter("(1) Sigma"));    //Voigt

  FixName(Fits[1],"AOffset",  0); //Offset = 0
  FixName(Fits[1],"Sigma",   FT->GetParameter("(1) Sigma"));    //Voigt

  FixName(Fits[2],"Constant",FT->GetParameter("(1) Constant")); //Voigt
  FixName(Fits[2],"Gamma",   FT->GetParameter("(1) Gamma"));    //Voigt
  FixName(Fits[2],"AOffset",  54.75); //Offset = Auger

  FixName(Fits[3],"Constant",FT->GetParameter("(1) Constant")); //Voigt
  FixName(Fits[3],"Offset",  FT->GetParameter("(2) Offset"));    //EMV
  FixName(Fits[3],"Sigma",   FT->GetParameter("(2) Sigma"));    //EMV
  FixName(Fits[3],"Tau",     FT->GetParameter("(2) Tau"));      //EMV
  FixName(Fits[3],"Gamma",   FT->GetParameter("(1) Gamma"));    //Voigt

  //FT
  parx[0] = {};
  // parx[1] = {{"R_LKCapt",0}};
  parx[1] = {};
  parx[2] = {{"B_Ratio",1},{"Decay",1},{"Energy",1}};
  parx[3] = {{"B_Ratio",1},{"Energy",1},{"AOffset",2}};
  parVec = FTParVec();
  delete FT; FT = new TF1(FTotal(r1,r2,parVec));
  PrintParameters(FT);

  //Fit
  result = R1->Fit(FT,"SQN0");
  ParSave(p); //Appends new parameters to old parameters
  TF1* DSAM;
  Fits.push_back(new TF1(FDSAM(aDSAMV,DSAM))); //DSAM before gaus for comparison
  {double p[10]; Fits[2]->GetParameters(p); Fits[4]->SetParameters(p);}

  //Write to File
  R1->SetTitle("(2) Peak 2/3");
  CF->SetTitle("(2) Peak 2/3");
  FT->SetTitle("(2) Peak 2/3");
  R2 = (TH1F*)R1->Clone("R1");
  FitExcl2 = FitExcl;
  r2=220; CanvG(CF,FT,R1,R2,"2_Pk2+3"); r2=900;
  FitExcl = FitExcl2;
  CanvG(CF,FT,R1,R2,"2_Pk2+3_900");
  // Rts->cd(); result->Write("R2");
  // Hgs->cd(); R2->Write("R2_Pk2+3");
  // Fns->cd(); FT->Write("F2");
  // Cvs->cd(); CF->Write("C2");

  //Reset Objects
  parx.clear();
  Fits.clear();

  //(3) Everything______________________________________________________________
  Fits.push_back(new TF1(FBg(r1,r2)));                    //BKG (Const1 Slope1 Const2 Slope2)
  Fits.push_back(new TF1(FVoigt(r1,r2)));                 //Voigt (Constant Centroid Sigma Gamma)
  FEMV1; aEMV1 = 2; cEMV1 = 2;
  Fits.push_back(new TF1(FExMoVoigt(aEMV1,cEMV1,FEMV1))); //EMV1 (Const R_Escape Centroid Offset Sigma Tau Gamma)
  Fits[2]->SetName("EMV1");
  FEMG0; aEMG0 = 1; cEMG0 = 1;
  Fits.push_back(new TF1(FExMoGaus(aEMG0,cEMG0,FEMG0)));  //EMG0 (Constant Centroid Offset Sigma Tau)
  Fits.push_back(new TF1(FGaus(r1,r2)));                  //Gaus (Constant R_LKCapt Centroid Sigma)
  DSAMG; aDSAMG = 3;
  Fits.push_back(new TF1(FDSAMG(aDSAMG,DSAMG)));          //DSAMG (Constant B_Ratio R_LKCapt Decay Energy AOffset Sigma)
  DSAMV; aDSAMV = 2;
  Fits.push_back(new TF1(FDSAMV(aDSAMV,DSAMV)));          //DSAMV (Constant B_Ratio Decay Energy AOffset Sigma Gamma)
  FEMV2; aEMV2 = 3; cEMV2 = 3;
  Fits.push_back(new TF1(FExMoVoigt(aEMV2,cEMV2,FEMV2))); //EMV2 (Constant R_Escape B_Ratio Energy AOffset Offset Sigma Tau Gamma)
  Fits[7]->SetName("EMV2");

  //Exclusion regions
  FitExcl.push_back({177,500});

  //Previous Parameters
  FTInitPars(IndFits, 0);
    int j=0;
    for(unsigned k=0;k<Fits.size();k++){
      for(int i=0;i<Fits[k]->GetNpar();i++){
        Fits[k]->SetParameter(i,p[j]);
        j++;
      }
    }
    for(unsigned i=0; i<Fits.size(); i++){
      const char* fname = IndFits[i]->GetTitle();
      Fits[i]->SetTitle(fname);
      for(int j=0; j<IndFits[i]->GetNpar(); j++){
        // double a,b; IndFits[i]->GetParLimits(j,a,b);
        // Fits[i]->SetParLimits(j,a,b);
        const char* pname=IndFits[i]->GetParName(j);
        Fits[i]->SetParName(j,pname);
      }
    }

  //Parameter Guesses
  FixName(Fits[6],"AOffset",  54.75); //Offset = Auger
  FixName(Fits[5],"AOffset",  0); //Offset = 0

  //FT
  parx[0] = {};                                                                                 //BKG (Const1 Slope1 Const2 Slope2)
  parx[1] = {};                                                                                 //Voigt (Constant Centroid Sigma Gamma)
  parx[2] = {{"Constant",1},{"Centroid",1},{"Gamma",1}};                                        //EMV1 (Const R_Escape Centroid Offset Sigma Tau Gamma)
  parx[3] = {};                                                                                 //EMG0 (Constant Centroid Offset Sigma Tau)
  parx[4] = {{"Constant",1}};                                                                   //Gaus (Constant R_LKCapt Centroid Sigma)
  // parx[5] = {{"Constant",1},{"R_LKCapt",4}};                                                    //DSAMG (Constant B_Ratio R_LKCapt Decay Energy AOffset Sigma)
  parx[5] = {{"Constant",1},{"Sigma",1},{"R_LKCapt",4}};                                     //DSAMG (Constant B_Ratio R_LKCapt Decay Energy AOffset Sigma)
  // parx[6] = {{"Constant",1},{"B_Ratio",5},{"Decay",5},{"Energy",5},{"Gamma",1}};                //DSAMV (Constant B_Ratio Decay Energy AOffset Sigma Gamma)
  parx[6] = {{"Constant",1},{"B_Ratio",5},{"Decay",5},{"Energy",5},{"Sigma",1},{"Gamma",1}}; //DSAMV (Constant B_Ratio Decay Energy AOffset Sigma Gamma)
  // parx[7] = {{"Constant",1},{"B_Ratio",5},{"Energy",5},{"AOffset",6},{"Tau",2},{"Gamma",1}};    //EMV2 (Constant R_Escape B_Ratio Energy AOffset Offset Sigma Tau Gamma)
  parx[7] = {{"Sigma",2},{"Tau",2},{"Gamma",1}};                                             //EMV2 (Constant R_Escape B_Ratio Energy AOffset Offset Sigma Tau Gamma)
  parVec = FTParVec();
  delete FT; FT = new TF1(FTotal(r1,r2,parVec));
  PrintParameters(FT);

  //Prelim Write to File
  H1->SetTitle("(3) Complete Fit");
  CF->SetTitle("(3) Complete Fit");
  FT->SetTitle("(3) Complete Fit");
  R1 = (TH1F*)H1->Clone("R1");
  FitExcl2 = FitExcl;
  H1->SetStats(0);
  r2=160; CanvG(CF,FT,H1,R1,"3_All_Prelim"); r2=900;
  H1->SetStats(1);
  FitExcl = FitExcl2;
  CanvG(CF,FT,H1,R1,"3_All_Prelim_900");
  FitExcl = FitExcl2;

  //Fit
  result = H1->Fit(FT,"SQN0");

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
  /*
*/
}
