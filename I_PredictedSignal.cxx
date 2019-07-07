//Plots of expected signal (ideal and STJ)
//Date: 06/28/2019
//Author: Spencer Fretwell

{
  #include "A1_Inclusions.cxx+"
  //ROOT Variables
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  //Initial Variables___________________________________________________________

  //Fit Range
  r1 = 10;
  r2 = 160;

  //Setup TObjects______________________________________________________________
  TCanvas* CF = new TCanvas("Fit","Fit Result",1920*2,2*1080); //Fit Canvas
  CF->cd(); CF->SetLogy();

  //Histograms
  TH1F* H1 = new TH1F("H1","Ideal Spectrum",(r2-r1)*10,r1,r2);
  H1->SetAxisRange(1,3e6,"Y");
  TH1F* H2=(TH1F*)H1->Clone();;
  H2->SetLineColor(kGreen);

  //Functions
  TF1* DSAM; int aDSAM = 1;
  TF1* fDSAM = new TF1(FDSAM(aDSAM,DSAM)); // C D E AO
  TF1* fVoigt = new TF1(FVoigt(r1,r2));
  TF1* fGaus = new TF1(FGaus(r1,r2));
  TF1* DSAMV; int aDSAMV = 1;
  TF1* fDSAMV = new TF1(FDSAMV(aDSAMV,DSAMV));
  TF1* DSAMG; int aDSAMG = 1;
  TF1* fDSAMG = new TF1(FDSAMG(aDSAMG,DSAMG));

  //Stack
  THStack* hs = new THStack("hs","Ideal Spectrum");

  //Ideal Spectrum______________________________________________________________
  fDSAM->SetParameters(1e6*0.1,1.908,28,56.8263);
  for(int i=0;i<H1->GetXaxis()->GetNbins();i++){//FDSAM: K-X
    H1->SetBinContent(i,H1->GetBinContent(i)+fDSAM->Eval(H1->GetBinLowEdge(i)));
  }

  fDSAM->SetParameters(1e6*0.1*0.0406,1.908,28,0);
  for(int i=0;i<H1->GetXaxis()->GetNbins();i++){//FSAM: L-X
    H1->SetBinContent(i,H1->GetBinContent(i)+fDSAM->Eval(H1->GetBinLowEdge(i)));
  }

  {
    int i=H1->GetXaxis()->FindBin(111.5263);
    H1->SetBinContent(i,H1->GetBinContent(i)+1e6); //Delta K-G

    i=H1->GetXaxis()->FindBin(56.8263);
    H1->SetBinContent(i,H1->GetBinContent(i)+1e6*0.0406); //Delta L-G
  }

  H1->Draw();
  CF->SaveAs("Output/Ideal_Spectrum.pdf");
  H1->Reset();
  delete hs; hs = new THStack("hs","Expected STJ Signal");

  //Expected STJ Signal_________________________________________________________
  fVoigt->SetParameters(1e6,111.5263,0.2,0.04);
  for(int i=0;i<H1->GetXaxis()->GetNbins();i++){//FVoigt: K-G
    H1->SetBinContent(i,H1->GetBinContent(i)+fVoigt->Eval(H1->GetBinLowEdge(i)));
    H2->SetBinContent(i,fVoigt->Eval(H2->GetBinLowEdge(i)));
  }
  hs->Add((TH1F*)H2->Clone());
  H2->Reset();

  fGaus->SetParameters(1e6,0.0406,56.8263,0.2);
  for(int i=0;i<H1->GetXaxis()->GetNbins();i++){//FGaus: L-G
    H1->SetBinContent(i,H1->GetBinContent(i)+fGaus->Eval(H1->GetBinLowEdge(i)));
    H2->SetBinContent(i,fGaus->Eval(H2->GetBinLowEdge(i)));
  }
  hs->Add((TH1F*)H2->Clone());
  H2->Reset();

  fDSAMV->SetParameters(1e6*0.1,1.908,28,56.8263,0.2,0.04);
  for(int i=0;i<H1->GetXaxis()->GetNbins();i++){//FDSAMV: K-X
    H1->SetBinContent(i,H1->GetBinContent(i)+fDSAMV->Eval(H1->GetBinLowEdge(i)));
    H2->SetBinContent(i,fDSAMV->Eval(H2->GetBinLowEdge(i)));
  }
  hs->Add((TH1F*)H2->Clone());
  H2->Reset();

  fDSAMG->SetParameters(1e6*0.1*0.0406,1.908,28,0,0.2);
  for(int i=0;i<H1->GetXaxis()->GetNbins();i++){//FDSAMG: L-X
    H1->SetBinContent(i,H1->GetBinContent(i)+fDSAMG->Eval(H1->GetBinLowEdge(i)));
    H2->SetBinContent(i,fDSAMG->Eval(H2->GetBinLowEdge(i)));
  }
  hs->Add((TH1F*)H2->Clone());
  H2->Reset();

  H1->SetNameTitle("H1","Expected STJ Signal");
  hs->Add(H1); hs->Draw("nostack");
  CF->SaveAs("Output/Expected_STJ_Signal.pdf");
  H1->Reset();
  delete hs; hs = new THStack("hs","Damaged STJ Signal");

  //Expected STJ Signal_________________________________________________________
  fVoigt->SetParameters(1e6,111.5263,3,0.04);
  for(int i=0;i<H1->GetXaxis()->GetNbins();i++){//FVoigt: K-G
    H1->SetBinContent(i,H1->GetBinContent(i)+fVoigt->Eval(H1->GetBinLowEdge(i)));
    H2->SetBinContent(i,fVoigt->Eval(H2->GetBinLowEdge(i)));
  }
  hs->Add((TH1F*)H2->Clone());
  H2->Reset();

  fGaus->SetParameters(1e6,0.0406,56.8263,3);
  for(int i=0;i<H1->GetXaxis()->GetNbins();i++){//FGaus: L-G
    H1->SetBinContent(i,H1->GetBinContent(i)+fGaus->Eval(H1->GetBinLowEdge(i)));
    H2->SetBinContent(i,fGaus->Eval(H2->GetBinLowEdge(i)));
  }
  hs->Add((TH1F*)H2->Clone());
  H2->Reset();

  fDSAMV->SetParameters(1e6*0.1,1.908,28,56.8263,3,0.04);
  for(int i=0;i<H1->GetXaxis()->GetNbins();i++){//FDSAMV: K-X
    H1->SetBinContent(i,H1->GetBinContent(i)+fDSAMV->Eval(H1->GetBinLowEdge(i)));
    H2->SetBinContent(i,fDSAMV->Eval(H2->GetBinLowEdge(i)));
  }
  hs->Add((TH1F*)H2->Clone());
  H2->Reset();

  fDSAMG->SetParameters(1e6*0.1*0.0406,1.908,28,0,3);
  for(int i=0;i<H1->GetXaxis()->GetNbins();i++){//FDSAMG: L-X
    H1->SetBinContent(i,H1->GetBinContent(i)+fDSAMG->Eval(H1->GetBinLowEdge(i)));
    H2->SetBinContent(i,fDSAMG->Eval(H2->GetBinLowEdge(i)));
  }
  hs->Add((TH1F*)H2->Clone());
  H2->Reset();

  H1->SetNameTitle("H1","Damaged STJ Signal");
  hs->Add(H1); hs->Draw("nostack");
  CF->SaveAs("Output/Damaged_STJ_Signal.pdf");
  /*
  */
}
