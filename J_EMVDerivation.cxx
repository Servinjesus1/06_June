{
  #include "A1_Inclusions.cxx+"
  TF1* FEMG = new TF1(FExMoGaus(-500,500));
  TF1* FL = new TF1("Lorentz","TMath::CauchyDist(x,0,[0])",-500,500);
  FL->SetNpx(2000);
  TF1Convolution f_conv(FEMG,FL,-500,500,kTRUE);
  f_conv.SetNofPointsFFT(100000);
  TF1* FC = new TF1("Convolved",f_conv,-500,500,f_conv.GetNpar());
  FC->SetNpx(2000); FC->SetLineColor(kBlue);
  FC->SetParName(0,"Constant");
  FC->SetParName(1,"R_Escape");
  FC->SetParName(2,"Centroid");
  FC->SetParName(3,"Offset");
  FC->SetParName(4,"Sigma");
  FC->SetParName(5,"Tau");
  FC->SetParName(6,"Gamma");
  FC->SetParameters(50,0.5,30,0,2,1,0.02);
  FEMG->SetParameters(50,0.5,30,0,2,1);
  FL->SetParameter(0,0.02);
  TCanvas* C1 = new TCanvas("C1","Convolved",1920,1080);
  C1->cd(); C1->SetLogy();
  FC->Draw(); FL->Draw("same");
  FEMG->Draw("same");
}
