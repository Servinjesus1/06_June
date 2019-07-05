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
  r2 = 300;

  //Files
  const char* fi1 = "Run1.root";
  const char* fi2 = "Run2.root";
  const char* fi3 = "Run3.root";
  const char* f_out = "./Output/G_Final/Run3_Fit.root";

  //Setup TObjects______________________________________________________________
  TCanvas* CF = new TCanvas("Fit","Fit Result",1920*2,2*1080); //Fit Canvas
  CF->cd(); CF->SetLogy();

  //Files
  TFile* F1 = new TFile(fi1);
  TFile* F2 = new TFile(fi2);
  TFile* F3 = new TFile(fi3);

  //Histograms
  TH1F* H1 = (TH1F*)F1->Get("Sum_Ch12+14");
  TH1F* H2 = (TH1F*)F2->Get("Sum Ch12+14");
  TH1F* H3 = (TH1F*)F3->Get("Ch_12+14");

  //Rebin H1____________________________________________________________________
  //First Guess: 0.2eV bins
  WidthRebin(H2,0.2);
  WidthRebin(H3,0.2);

  //Find Peaks
  H1->SetAxisRange(25,160);
  H2->SetAxisRange(25,160);
  H3->SetAxisRange(25,160);
  TSpectrum* TS1 = new TSpectrum(4);
  TSpectrum* TS2 = new TSpectrum(4);
  TSpectrum* TS3 = new TSpectrum(4);
  TS1->Search(H1,3,"",0.0005);
  TS2->Search(H2,5,"",0.0005);
  TS3->Search(H3,5,"",0.0005);
  H1->SetAxisRange(r1,r2);
  H2->SetAxisRange(r1,r2);
  H3->SetAxisRange(r1,r2);

  //Get Peak Centroids
  Double_t peaksX1[4],peaksY1[4];
  Double_t peaksX2[4],peaksY2[4];
  Double_t peaksX3[4],peaksY3[4];
  TPolyMarker* PKs1 = (TPolyMarker*)H1->GetListOfFunctions()->FindObject("TPolyMarker");
  TPolyMarker* PKs2 = (TPolyMarker*)H2->GetListOfFunctions()->FindObject("TPolyMarker");
  TPolyMarker* PKs3 = (TPolyMarker*)H3->GetListOfFunctions()->FindObject("TPolyMarker");
  for(int i=0;i<4;i++){
    peaksX1[i]=PKs1->GetX()[i]; //Peak X values
    peaksX2[i]=PKs2->GetX()[i]; //Peak X values
    peaksX3[i]=PKs3->GetX()[i]; //Peak X values
    peaksY1[i]=H1->GetBinContent(H1->GetXaxis()->FindBin(peaksX1[i])); //Peak Y values
    peaksY2[i]=H2->GetBinContent(H2->GetXaxis()->FindBin(peaksX2[i])); //Peak Y values
    peaksY3[i]=H3->GetBinContent(H3->GetXaxis()->FindBin(peaksX3[i])); //Peak Y values
  }
  H1->GetListOfFunctions()->Remove(PKs1); //Remove Search Carats
  H2->GetListOfFunctions()->Remove(PKs2); //Remove Search Carats
  H3->GetListOfFunctions()->Remove(PKs3); //Remove Search Carats

  //Second Guess: ((Peak 1 -> 111.58eV),(Peak 3 -> 56.83eV))
  LinRebin(H1,{{peaksX1[0],111.58},{peaksX1[1],56.83}});
  LinRebin(H2,{{peaksX2[0],111.58},{peaksX2[1],56.83}});
  LinRebin(H3,{{peaksX3[0],111.58},{peaksX3[1],56.83}});

  //Find Peaks Again
  H1->SetAxisRange(25,160);
  H2->SetAxisRange(25,160);
  H3->SetAxisRange(25,160);
  TS1->Search(H1,3,"",0.0005);
  TS2->Search(H2,5,"",0.0005);
  TS3->Search(H3,5,"",0.0005);
  H1->SetAxisRange(r1,r2);
  H2->SetAxisRange(r1,r2);
  H3->SetAxisRange(r1,r2);
  PKs1 = (TPolyMarker*)H1->GetListOfFunctions()->FindObject("TPolyMarker");
  PKs2 = (TPolyMarker*)H2->GetListOfFunctions()->FindObject("TPolyMarker");
  PKs3 = (TPolyMarker*)H3->GetListOfFunctions()->FindObject("TPolyMarker");
  for(int i=0;i<4;i++){
    peaksX1[i]=PKs1->GetX()[i]; //Peak X values
    peaksX2[i]=PKs2->GetX()[i]; //Peak X values
    peaksX3[i]=PKs3->GetX()[i]; //Peak X values
    peaksY1[i]=H1->GetBinContent(H1->GetXaxis()->FindBin(peaksX1[i])); //Peak Y values
    peaksY2[i]=H2->GetBinContent(H2->GetXaxis()->FindBin(peaksX2[i])); //Peak Y values
    peaksY3[i]=H3->GetBinContent(H3->GetXaxis()->FindBin(peaksX3[i])); //Peak Y values
  }
  H1->GetListOfFunctions()->Remove(PKs1); //Remove Search Carats
  H2->GetListOfFunctions()->Remove(PKs2); //Remove Search Carats
  H3->GetListOfFunctions()->Remove(PKs3); //Remove Search Carats

  //Scaling
  H1->Sumw2();
  H2->Sumw2();
  H3->Sumw2();
  // H1->Scale(1., "width"); //Scale? From fitNormSum tutorial
  // H2->Scale(1., "width"); //Scale? From fitNormSum tutorial
  // H3->Scale(1., "width"); //Scale? From fitNormSum tutorial
  H2->Scale(peaksY1[0]/peaksY2[0]); //Scale? From fitNormSum tutorial
  H3->Scale(peaksY1[0]/peaksY3[0]); //Scale? From fitNormSum tutorial

  //Final Tweaks________________________________________________________________
  //Set H1 Name/Title
  H1->SetNameTitle("H1","Run Comparison"); H1->SetAxisRange(r1,r2,"X");
  H2->SetNameTitle("H2","Run Comparison"); H2->SetAxisRange(r1,r2,"X");
  H3->SetNameTitle("H3","Run Comparison"); H3->SetAxisRange(r1,r2,"X");

  gStyle->SetOptStat(0);
  H1->SetLineWidth(2);
  H1->Draw("hist");
  H2->SetLineColor(kGreen);
  H2->Draw("hist same");
  H3->SetLineColor(kRed);
  H3->Draw("hist same");
  CF->SaveAs("Output/Run_Comparison.pdf")
  /*
  */
}
