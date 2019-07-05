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

  //Rebin H1____________________________________________________________________
  //First Guess: 0.2eV bins
  WidthRebin(H1,0.2);
  TH1F* R1 = (TH1F*)H1->Clone("Prebin"); TH1F* R2;

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
  R1->SetLineColor(kRed);
  CF->SetLogy(); H1->SetAxisRange(10,300,"X");
  H1->SetAxisRange(5e2,1e7,"Y"); H1->SetTitle("Rebin Comparison");
  H1->Draw();
  R1->Draw("Same");
  CF->SaveAs("Output/LinRebinComp.pdf")
}
