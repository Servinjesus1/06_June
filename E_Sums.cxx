//Macro 4: Fit all four peaks additively together
//DSAM Method Peak, ExModGaus peaks too
//Date: 06/07/2019
//Author: Spencer Fretwell

{
  #include <A1_Inclusions.cxx+>
  //ROOT Variables
  //gROOT->SetBatch(kTRUE);

  //Initial Variables___________________________________________________________

  //Fit Range
  r1 = 10;
  r2 = 900;

  //Files
  const char* f_in = "Run3.root";
  const char* f_out = "./Output/E_Sums/Run3_Fit.root";

  //Setup TObjects______________________________________________________________
  TF1* FT; //Total Function, sum of Fits[]
  TFitResultPtr result; //Fit Results
  TMatrix Rmx;
  TCanvas* CF = new TCanvas("Fit","Fit Result",1920,1080); //Fit Canvas

  //Files
  TFile* FI = new TFile(f_in);
  TFile* FO = new TFile(f_out,"RECREATE");
  // TFile* FP = new TFile("./Output/F_Individuals/Run3_Fit.root","READ");

  //Folders
  TDirectory* Rts = FO->mkdir("Results");
  TDirectory* Hgs = FO->mkdir("Histograms");
  TDirectory* Fns = FO->mkdir("Functions");
  TDirectory* Cvs = FO->mkdir("Canvases");

  //Histograms
  TH1F* H1 = (TH1F*)FI->Get(FI->GetListOfKeys()->At(1)->GetName()); //Original H
  TH1F* R1;

  //Past Fits
  // vector<TF1*> IndFits;
  // IndFits.resize(FP->cd("Functions")->GetListOfObjects()->GetSize());
  // for(int i=0;i<IndFits.size();i++){
  //   IndFits[i]=(TF1*)FP->cd("Functions")->GetListOfObjects()[i];
  // }

  //Rebin H1____________________________________________________________________
  //First Guess: 0.2eV bins
  WidthRebin(H1,0.2);

  //Find Peaks
  H1->SetAxisRange(25,160);
  TSpectrum* TS1 = new TSpectrum(4);
  TS1->Search(H1,20,"",0.0005);
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

  //Set H1 Name/Title
  H1->SetNameTitle("H1","Consecutive Fits, Full");

  //Find Peaks
  H1->SetAxisRange(r1,160);
  TS1 = new TSpectrum(4);
  TS1->Search(H1,20,"",0.0005);
  H1->SetAxisRange(r1,r2);

  //Get Peak Centroids
  PKs = (TPolyMarker*)H1->GetListOfFunctions()->FindObject("TPolyMarker");
  for(int i=0;i<4;i++){
    peaksX[i]=PKs->GetX()[i]; //Peak X values
    peaksY[i]=H1->GetBinContent(peaksX[i]); //Peak Y values
  }
  H1->GetListOfFunctions()->Remove(PKs); //Remove Search Carats
  CF->Clear();
  Hgs->cd(); H1->Write("H1");

//_______________________________FITTING________________________________________
//______________________________________________________________________________
  //(1) Background______________________________________________________________
  Fits.push_back(new TF1(FBg(r1,r2)));  //Create Background Fn (4 parameters)
  int curfit=Fits.size()-1;

  //Exclusion regions
  FitExcl.push_back({20,40});   //BG Excl Regions
  FitExcl.push_back({42,168});  //BG Excl Regions
  FitExcl.push_back({177,350}); //BG Excl Regions

  //FT
  FT = new TF1(FTotal(r1,r2)); //Create Total Function

  //Fit
  result = H1->Fit(FT,"SQ0R");
  ParCopy(FT); //Copy fit parameters to Fits

  //Write to File
  H1->SetTitle("(1) Background");
  CF->SetTitle("(1) Background");
  FT->SetTitle("(1) Background");
  R1 = (TH1F*)H1->Clone("R1");
  CanvE(CF,FT,H1,R1,"1_Bkg");
  Rts->cd(); result->Write("R1");
  Hgs->cd(); R1->Write("R1_Bkg");
  Fns->cd(); FT->Write("F1");
  Cvs->cd(); CF->Write("C1");
  /*

  //(2) Fit First Peak: 7LiK Voigt______________________________________________
  Fits.push_back(new TF1(FVoigt(peaksX[0]-1,120))); //Create Voigt Fn (4 parameters)
  curfit=Fits.size()-1;

  //Parameter Guesses
  Fits[curfit]->SetParameter("Centroid",peaksX[0]);
  Fits[curfit]->SetParameter("Constant",peaksY[0]);
  Fits[curfit]->FixParameter(3,0.03);

  //Exclusion regions
  FitExcl.push_back({20,40});
  FitExcl.push_back({42,peaksX[0]-1});
  FitExcl.push_back({120,168});
  FitExcl.push_back({177,500});

  //FT      {A,C,S,G}
  parx[1] = {4,5,6,7};
  delete FT;
  FT = new TF1(FTotal(r1,r2)); //Create Total Function

  //Fit
  result = H1->Fit(FT,"SQ0R");
  ParCopy(FT);

  //Write to File
  r2=160;
  H1->SetTitle("(2) Voigt");
  CF->SetTitle("(2) Voigt");
  FT->SetTitle("(2) Voigt");
  R1 = (TH1F*)H1->Clone("R1");
  CanvE(CF,FT,H1,R1,"2_Vgt");
  Rts->cd(); result->Write("R2");
  Hgs->cd(); R1->Write("R2_Vgt");
  Fns->cd(); FT->Write("F2");
  Cvs->cd(); CF->Write("C2");
  r2=900;

  //(3) Fit First Peak: 7LiK ExModGaus__________________________________________
  Fits.push_back(new TF1(FExMoGaus(99,peaksX[0]+5))); //Create EMG Fn (4 parameters)
  curfit=Fits.size()-1;

  //Exclusion regions
  FitExcl.push_back({20,40});
  FitExcl.push_back({42,99});
  FitExcl.push_back({120,168});
  FitExcl.push_back({177,500});

  //FT      {A,C,!S,T} //!:bound parameter
  // parx[2] = {8,9,6,10};
  parx[2] = {8,9,10,11};
  delete FT;
  FT = new TF1(FTotal(r1,r2)); //Create Total Function

  //Fit
  result = H1->Fit(FT,"SQ0R");
  ParCopy(FT);

  //Write to File
  r2=160;
  H1->SetTitle("(3) ExModGaus");
  CF->SetTitle("(3) ExModGaus");
  FT->SetTitle("(3) ExModGaus");
  R1 = (TH1F*)H1->Clone("R1");
  CanvE(CF,FT,H1,R1,"3_EMG");
  Rts->cd(); result->Write("R3");
  Hgs->cd(); R1->Write("R3_EMG");
  Fns->cd(); FT->Write("F3");
  Cvs->cd(); CF->Write("C3");
  r2=900;

  //(4) Fit Second Peak: 7Li*K DSAMG____________________________________________
  // TF1* FC; TF1* FD; TF1* FG; TF1Convolution* f_conv; //Convolved Functions floating around
  Fits.push_back(new TF1(FDSAMG(peaksX[1]-2,99.5))); //Create DSAMG Fn (5 parameters)
  curfit=Fits.size()-1;

  //Parameter Guesses
  Fits[curfit]->SetParameter("Constant",peaksY[1]);
  // Fits[curfit]->FixParameter(3,54.75);

  //Fit Attempt
  result = R1->Fit(Fits[curfit],"SQ0RM");
  result->Print();
  r2=160;
  R1->SetTitle("(4) DSAMG - Initial Try");
  CanvE(CF,Fits[curfit],R1,(TH1F*)R1->Clone(),"4_DSAMG_0");
  r2=900;

  //Exclusion regions
  FitExcl.push_back({20,40});
  FitExcl.push_back({42,peaksX[1]-2});
  FitExcl.push_back({120,168});
  FitExcl.push_back({177,500});

  //FT      {A, D, E, O,!S}
  // parx[3] = {11,12,13,14,6};
  parx[3] = {12,13,14,15,16};
  delete FT;
  FT = new TF1(FTotal(r1,r2)); //Create Total Function

  //Fit
  result = H1->Fit(FT,"SQ0R");
  ParCopy(FT);

  //Write to File
  r2=160;
  H1->SetTitle("(4) DSAMG");
  CF->SetTitle("(4) DSAMG");
  FT->SetTitle("(4) DSAMG");
  R1 = (TH1F*)H1->Clone("R1");
  CanvE(CF,FT,H1,R1,"4_DSAMG");
  Rts->cd(); result->Write("R4");
  Hgs->cd(); R1->Write("R4_DSAMG");
  Fns->cd(); FT->Write("F4");
  Cvs->cd(); CF->Write("C4");
  r2=900;

  //(5) Fit Second Peak: 7Li*K ExModGaus________________________________________
  Fits.push_back(new TF1(FExMoGaus(69,peaksX[1]))); //Create EMG Fn (4 parameters)
  curfit=Fits.size()-1;

  //Exclusion regions
  FitExcl.push_back({20,40});
  FitExcl.push_back({42,69});
  FitExcl.push_back({120,168});
  FitExcl.push_back({177,500});

  //FT      {A, C,!S, T}
  // parx[4] = {15,16,6,17};
  parx[4] = {17,18,19,20};
  delete FT;
  FT = new TF1(FTotal(r1,r2)); //Create Total Function

  //Fit
  result = H1->Fit(FT,"SQ0RM");
  ParCopy(FT);

  //Write to File
  r2=160;
  H1->SetTitle("(5) ExModGaus");
  CF->SetTitle("(5) ExModGaus");
  FT->SetTitle("(5) ExModGaus");
  R1 = (TH1F*)H1->Clone("R1");
  CanvE(CF,FT,H1,R1,"5_EMG");
  Rts->cd(); result->Write("R5");
  Hgs->cd(); R1->Write("R4_DSAMG");
  Fns->cd(); FT->Write("F5");
  Cvs->cd(); CF->Write("C5");
  r2=900;

  //(6) Fit Third Peak: 7LiL Gaus_______________________________________________
  Fits.push_back(new TF1(FGaus(43,67))); //Create Gauss Fn (3 parameters)
  curfit=Fits.size()-1;

  //Parameter Guesses
  Fits[curfit]->SetParameter("Centroid",peaksX[2]);
  Fits[curfit]->SetParameter("Constant",peaksY[2]);
  Fits[curfit]->SetParameter("Sigma",Fits[1]->GetParameter("Sigma"));

  //Exclusion regions
  FitExcl.push_back({20,40});
  FitExcl.push_back({42,50});
  FitExcl.push_back({120,168});
  FitExcl.push_back({177,500});

  //FT      {A, C,!S}
  // parx[5] = {18,19,6};
  parx[5] = {21,22,23};
  delete FT;
  FT = new TF1(FTotal(r1,r2)); //Create Total Function

  //Fit
  result = H1->Fit(FT,"SQ0RM");
  ParCopy(FT);

  //Write to File
  r2=160;
  H1->SetTitle("(6) Gaus");
  CF->SetTitle("(6) Gaus");
  FT->SetTitle("(6) Gaus");
  R1 = (TH1F*)H1->Clone("R1");
  CanvE(CF,FT,H1,R1,"6_Gaus");
  Rts->cd(); result->Write("R6");
  Fns->cd(); FT->Write("F6");
  Cvs->cd(); CF->Write("C6");
  r2=900;

  //(7) Fit Fourth Peak: 7Li*L DSAMG____________________________________________
  // TF1* FC2; TF1* FD2; TF1* FG2; TF1Convolution* f_conv2; //Convolved Functions floating around
  Fits.push_back(new TF1(FDSAMG(23,42))); //Create DSAMG Fn (5 parameters)
  curfit=Fits.size()-1;

  //Parameter Guesses
  Fits[curfit]->SetParameter("Centroid",peaksX[3]);
  // Fits[curfit]->FixParameter(3,0);

  //Exclusion regions
  FitExcl.push_back({119.5,168});
  FitExcl.push_back({177,500});

  //FT      {A,!D,!E,O,!S} //!:bound parameter
  // parx[6] = {20,12,13,21,6};
  parx[6] = {24,25,26,27,28};
  delete FT;
  FT = new TF1(FTotal(r1,r2)); //Create Total Function

  //Fit
  result = H1->Fit(FT,"SQ0RM");
  ParCopy(FT);

  //Write to File
  r2=160;
  H1->SetTitle("(7) DSAMG");
  CF->SetTitle("(7) DSAMG");
  FT->SetTitle("(7) DSAMG");
  R1 = (TH1F*)H1->Clone("R1");
  CanvE(CF,FT,H1,R1,"7_DSAMG");
  Rts->cd(); result->Write("R7");
  Hgs->cd(); R1->Write("R4_DSAMG");
  Fns->cd(); FT->Write("F7");
  Cvs->cd(); CF->Write("C7");
  FO->Close();
  */
}
