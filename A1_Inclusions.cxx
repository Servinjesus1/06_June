//STJ Peak Fit Functions
//March 30, 2019
//Spencer Fretwell

#include <TStyle.h>
#include <TROOT.h>
#include <iostream>
#include <string>
#include "TPad.h"
#include "TPolyMarker.h"
#include "TSpectrum.h"
#include "TH1.h"
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TPaveStats.h"
#include "TVirtualPad.h"
#include "TF1Convolution.h"
#include "TF1NormSum.h"
#include "TExec.h"
#include "TGaxis.h"
#include "TComplex.h"

//Global Scope Variables________________________________________________________
vector<TF1*> Fits;              //Background plus npeaks
vector<vector<double>> FitExcl; //Exclusion Ranges
vector<vector<double>> FitExcl2; //Copy of Exclusion Ranges
map<int,vector<pair<const char*,int>>> parx; //Parameter Merge Map: {"parname",Source Fn}
Double_t r1,r2;                 //Fit Range

//Important Macros______________________________________________________________
//Get Max in a map of arrays (total # of parameters)
template<typename KeyType, typename ValueType> ValueType get_max( const std::map<KeyType,vector<ValueType>>& x ) {
  using pairtype=std::pair<KeyType,vector<ValueType>>;
  auto p = (*std::max_element(x.begin(), x.end(), [] (const pairtype & p1, const pairtype & p2) {
    return *max_element(p1.second.begin(),p1.second.end()) < *max_element(p2.second.begin(),p2.second.end());
  })).second;
  return *max_element(p.begin(),p.end());
}

//Check if parameters are out of bounds
void ParCheck(TF1* F){
  for(int i=0;i<F->GetNpar();i++){
    double a,b;
    F->GetParLimits(i,a,b);
    if((a!=0&&b!=0)&&(F->GetParameter(i)>b||F->GetParameter(i)<a)){
      cout << a << " < " << F->GetParName(i) << " = " << F->GetParameter(i) << " < " << b << endl;
    }
  }
}

//Signum
template <typename T> inline constexpr
int signum(T x, std::false_type is_signed) {
    return T(0) < x;
}

template <typename T> inline constexpr
int signum(T x, std::true_type is_signed) {
    return (T(0) < x) - (x < T(0));
}

template <typename T> inline constexpr
int signum(T x) {
    return signum(x, std::is_signed<T>());
}

//ERFCX Function
Double_t Erfcx(Double_t A){
  return TMath::Exp(A*A)*TMath::Erfc(A);
}

//Pretty parameter names for full function
string FTParName(int i,const char* pr){
  string nm = ((string)"(" + to_string(i) + (string)") " + (string)pr);
  return nm;
}

//Range Check (FitExcl)
Bool_t ExclRange(Double_t x){
  int InRng = 0;
  for(auto q:FitExcl){
    if(q[0]<=x&&x<=q[1]){InRng++;}
  }
  if(InRng>0){return kTRUE;} else{return kFALSE;}
}

//Merge Map to FT Parameter List
map<int,vector<int>> FTParVec(){          //bindings stored as {"par",4} {"par",fi}
  int totPars=0;                          //translate to map of p[i] for FT
  map<int,vector<int>> v;
  for(unsigned i=0;i<Fits.size();i++){ //For each Fit Function...
    v[i].resize(Fits[i]->GetNpar());    //Resize V by # of parameters...
    if(!parx[i].empty()){                //if parx has corresponding bindings..
      for(auto & mp : parx[i]){            //for each binding in parx...
        int pnumo = Fits[i]->GetParNumber(mp.first); //output (bound) fn
        int pnumi = Fits[mp.second]->GetParNumber(mp.first); //input (binding) fn
        v[i][pnumo]=v[mp.second][pnumi];     //assign numbers accordingly
      }
    }
    for(auto & p : v[i]) p=(p==0)?totPars++:p; //iterate the leftover zeros.
  }
  return v;
}

//SetParLimitses
void SetParLims(TF1 &F, vector<vector<double>> v){
  for(int i=0; i<F.GetNpar(); i++){
    F.SetParLimits(i,v[i][1],v[i][2]);
  }
}

//Functional Forms______________________________________________________________
Double_t FFBg(Double_t* x, Double_t* p){ //Functional Form of Background
  //Null Condition
  if(p[0]==0&&p[1]==0){return 0;}
  //Exclusion Condition
  if(ExclRange(x[0])){TF1::RejectPoint(); return 0;}
  //background
  return p[0]*exp(-p[1]*x[0])+p[2]*exp(-p[3]*x[0]);
}

Double_t FFGaus0(Double_t* x, Double_t* p){ //Functional Form of Gaussian
  //Exclusion Condition
  // if(ExclRange(x[0])){TF1::RejectPoint();}
  //Gaussian
  return TMath::Gaus(x[0],0,p[0],1);
}

Double_t FFGaus(Double_t* x, Double_t* p){ //Functional Form of Gaussian
  //Exclusion Condition
  if(ExclRange(x[0])){TF1::RejectPoint(); return 0;}
  //Gaussian
  return p[0]*p[1]*TMath::Gaus(x[0],p[2],p[3],1);
}

Double_t FFVoigt0(Double_t* x, Double_t* p){ //Functional Form of Voigt
  //Exclusion Condition
  // if(ExclRange(x[0])){TF1::RejectPoint();}
  //Voigt
  return TMath::Voigt(x[0],p[0],p[1]);
}

Double_t FFVoigt(Double_t* x, Double_t* p){ //Functional Form of Gaussian
  //Null Condition
  if(p[0]==0&&p[1]==0){return 0;}
  //Exclusion Condition
  if(ExclRange(x[0])){TF1::RejectPoint(); return 0;}
  //Voigt
  return p[0]*TMath::Voigt(x[0]-p[1],p[2],p[3]);
}

Double_t FFCauchy0(Double_t* x, Double_t* p){ //Functional Form of Lorentzian
  if(p[0]==0){return 0;}
  //Exclusion Condition
  // if(ExclRange(x[0])){TF1::RejectPoint();}
  //Voigt
  return TMath::CauchyDist(x[0],0,p[0]);
}

Double_t FFExMoGaus0(Double_t* x, Double_t* p){ //Constant-Free EMG to be normalized
  if(p[0]==0){return 0;}
  //First two parameters reserved for Amplitude
  Double_t S = p[0]; //Sigma
  Double_t T = p[1]; //Tau
  Double_t X = -x[0];
  if(T<0){X=-X; T=-T;} //Tau<0 Shorthand for switching exponential tail
  Double_t X1 = X/S;
  Double_t X2 = X/T;
  Double_t Z = 1/TMath::Sqrt(2)*(S/T-X1);
  Double_t Z1 = 0.5*TMath::Power((S/T),2)-X2;
  if(Z<0){return S/T*TMath::Sqrt(TMath::Pi()/2)*TMath::Exp(Z1)*TMath::Erfc(Z);}
  if(TMath::Exp(Z*Z)==HUGE_VAL){return TMath::Exp(-0.5*X1*X1)/(1+T/S*X1);}//TEMPORARY TEST
  if(0<=Z<=6.71E10){return S/T*TMath::Sqrt(TMath::Pi()/2)*TMath::Exp(-0.5*X1*X1)*Erfcx(Z);}
  if(Z>6.71E10){return TMath::Exp(-0.5*X1*X1)/(1+T/S*X1);}
  return 0;
}

Double_t FFNormAmpN(int &a, TF1* & F, Double_t* x, Double_t* p){
  //For normalized functions to exclude and apply amplitude and centroid parameters
  if(ExclRange(x[0])){TF1::RejectPoint(); return 0;}

  int ALT = 0;
  for(int i=0; i<F->GetNpar(); i++){
    if(F->GetParameter(i)!=p[i+a]){
      ALT++;
      F->SetParameter(i,p[i+a]);
      // cout << F->GetName() << "(" << i << ") is now " << F->GetParameter(i) << endl;
    }
  }
  Double_t A = 1; for(int i=0;i<a;i++) A*=p[i]; //a constants multiplied together
  if(ALT>0){
    F->Update(); //For Normalization Integral
  }

  return A*F->Eval(x[0]); //Does Eval work?
}

Double_t FFNormAmpNC(int &a, int &c, TF1* & F, Double_t* x, Double_t* p){
  //For normalized functions to exclude and apply amplitude and centroid parameters
  if(ExclRange(x[0])){TF1::RejectPoint(); return 0;}

  int ALT = 0;
  for(int i=0; i<F->GetNpar(); i++){
    if(F->GetParameter(i)!=p[i+a+c]){
      ALT++;
      F->SetParameter(i,p[i+a+c]);
    }
  }
  Double_t A = 1; for(int i=0;i<a;i++) A*=p[i]; //a constants multiplied together
  Double_t C = 0; for(int i=0;i<c;i++) C+=p[i+a]; //c offsets added together
  if(ALT!=0){
    F->Update(); //For Normalization Integral
  }

  return A*F->Eval(x[0]-C);
}

Double_t FFDSAM0(Double_t* x, Double_t* p){ //Constant-Free DSAM Function
  if(ExclRange(x[0])){TF1::RejectPoint();}
  //Does not exclude points
  if(p[0]==0&&p[1]==0){return 0;}
  Double_t E = x[0];
  Double_t S = p[0]; //Decay parameter
  Double_t EE = p[1]; //Energy Paramter
  Double_t E0 = p[2]; //AOffset Parameter
  if(E0<1){E0=0;} //Since I can't seem to fix a parameter at zero...
  Double_t F = (2/EE)*(S/(S-1))*(1-TMath::Power(2*TMath::Abs(E-E0-EE)/EE,S-1));
  if(F<0) return 0;
  return F;
}

// Double_t FFDSAM(Double_t* x, Double_t* p){//Normalized DSAM * Amplitude
//   TF1* F = (TF1*)gROOT->GetFunction("DSAMfn");
//
//   for(int i=0; i<F->GetNpar(); i++){
//       F->SetParameter(i,p[i+2]); //Take passed parameters and set to DSAM
//   }
//
//   F->Update(); //Normalizes DSAMNorm
//   Double_t A = p[0]*p[1];
//   return A*F->Eval(x[0]);
// }
//
// Double_t FFDSAMC(Double_t* x, Double_t* p){//Normalized DSAMConvolved * Amplitude
//   TF1* F = (TF1*)gROOT->GetFunction("DSAMconv");
//
//   for(int i=0; i<F->GetNpar(); i++){
//       F->SetParameter(i,p[i+2]); //Take passed parameters and set to DSAM
//   }
//
//   F->Update(); //Normalizes DSAMNorm
//   Double_t A = p[0]*p[1];
//   return A*F->Eval(x[0]);
// }

Double_t FFTotal(Double_t* x, Double_t* p, map<int,vector<int>> & parVec) {
  if(ExclRange(x[0])){TF1::RejectPoint(); return 0;}

  //Check for Parameter differences
  int ALT = 0;
  for(unsigned i=0; i<parVec.size(); i++){//Update fits in parVec if changed
    for(int j=0;j<Fits[i]->GetNpar();j++){
      if(Fits[i]->GetParameter(j)!=p[parVec[i][j]]){
        ALT++;
        Fits[i]->SetParameter(j,p[parVec[i][j]]);
      }
    }
    if(ALT>0){Fits[i]->Update(); ALT=0;}//Call Update if any parameters change
  }

  //Evaluate
  Double_t FT=0;
  for(auto y : parVec){
    int i = y.first;
    FT+=Fits[i]->Eval(x[0]);
  }

  return FT;
}

// Double_t FFTotSimple(Double_t* x, Double_t* p) {
//   //EvalPar
//   Double_t FT=0;
//   for(auto y : Fits){
//     FT+=y->Eval(*x);
//   }
//   return FT;
// }

// Double_t FReject(TF1* F,Double_t* x, Double_t* p) {//Functional to exclude
//   cout << p[0] << endl;
//   cout << "FReject Called" << endl;
//   //Exclude Points Now (e.g. post convolution)
//   if(ExclRange(x[0])){TF1::RejectPoint(); return 0;}
//   //EvalPar
//   return F->EvalPar(x,p);
// }

//TF1 Forms_____________________________________________________________________
TF1 FBg(Double_t r1, Double_t r2){
  TF1 F("Background",FFBg,r1,r2,4);
  F.SetNpx(2000); F.SetFillColorAlpha(kAzure-9,0.5); F.SetFillStyle(1001);
  F.SetLineStyle(2); F.SetLineWidth(1);
  int i=-1;
  i++; F.SetParName(i,"Const1"); F.SetParameter(i,1e4);  F.SetParLimits(i,0,1e7);
  i++; F.SetParName(i,"Slope1"); F.SetParameter(i,0.1);  F.SetParLimits(i,0,1);
  i++; F.SetParName(i,"Const2"); F.SetParameter(i,1e3);  F.SetParLimits(i,0,1e7);
  i++; F.SetParName(i,"Slope2"); F.SetParameter(i,0.01); F.SetParLimits(i,0,1);
  return F;
}

TF1 FVoigt(Double_t r1, Double_t r2){
  TF1 F("Voigt",FFVoigt,r1,r2,4);
  F.SetNpx(2000); F.SetFillColorAlpha(kAzure-9,0.5); F.SetFillStyle(1001);
  F.SetLineColor(kGreen); F.SetLineWidth(1);
  int i=-1;
  i++; F.SetParName(i,"Constant"); F.SetParameter(i,0);    F.SetParLimits(i,10,1e8);
  i++; F.SetParName(i,"Centroid"); F.SetParameter(i,0);    F.SetParLimits(i,r1-10,r2+10);
  i++; F.SetParName(i,"Sigma");    F.SetParameter(i,2);    F.SetParLimits(i,0,15);
  i++; F.SetParName(i,"Gamma");    F.SetParameter(i,0.04); F.SetParLimits(i,0.02,0.06);
  return F;
}

TF1 FGaus(Double_t r1, Double_t r2){
  TF1 F("LK_Gaussian",FFGaus,r1,r2,4);
 	F.SetNpx(2000); F.SetFillColorAlpha(kAzure-9,0.5); F.SetFillStyle(1001);
  F.SetLineStyle(1); F.SetLineColor(kGreen); F.SetLineWidth(1);
  int i=-1;
  i++; F.SetParName(i,"Constant"); F.SetParameter(i,1); F.SetParLimits(i,10,1e8);
  i++; F.SetParName(i,"R_LKCapt"); F.SetParameter(i,0.04); F.SetParLimits(i,0,0.1);
  i++; F.SetParName(i,"Centroid"); F.SetParameter(i,1); F.SetParLimits(i,r1,r2);
  i++; F.SetParName(i,"Sigma");    F.SetParameter(i,1); F.SetParLimits(i,0,15);
  return F;
}

TF1 FExMoGaus(int &a, int &c, TF1* & FN){
  FN = new TF1("EMGNorm",FFExMoGaus0,-300,300,2);
  FN->SetNormalized(kTRUE);
  TF1 F("EMG",[&](Double_t* x, Double_t* p){return FFNormAmpNC(a,c,FN,x,p);},r1,r2,2+a+c);
  F.SetNpx(2000); F.SetFillColorAlpha(kAzure-9,0.5); F.SetFillStyle(1001);
  F.SetLineColor(kGreen); F.SetLineWidth(1);
  int i=-1;
  i++; F.SetParName(i,"Constant"); F.SetParameter(i,1E5);       F.SetParLimits(i,10,1e8);
  while(i<a-1){
  i++; F.SetParName(i,"R_Escape"); F.SetParameter(i,0.1);       F.SetParLimits(i,0,1);
  }
  i++; F.SetParName(i,"Centroid"); F.SetParameter(i,(r1+r2)/2); F.SetParLimits(i,r1,r2);
  while(i<a+c-1){
  i++; F.SetParName(i,"Offset");      F.SetParameter(i,-5);       F.SetParLimits(i,-10,0);
  }
  i++; F.SetParName(i,"Sigma");    F.SetParameter(i,1);         F.SetParLimits(i,0,15);
  i++; F.SetParName(i,"Tau");      F.SetParameter(i,1);         F.SetParLimits(i,0,15);

  return F;
}

TF1 FExMoVoigt(int &a, int &c, TF1*& FC){
  //Convolve
  TF1* FEMV = new TF1("EMG",FFExMoGaus0,-300,r2+500,2);
  //FEMV->SetNpx(500); FEMV->SetParameters(1,5);
  TF1* FL = new TF1("Lorentzian Function",FFCauchy0,-300,r2+500,1);
  //FL->SetNpx(500); FL->SetParameters(1,1);
  TF1Convolution* ConvEMG = new TF1Convolution(FEMV,FL,true);
  ConvEMG->SetRange(-300,r2+500); //Implements AbsComposition?
  ConvEMG->SetNofPointsFFT(1e6); //FFT Points?
  FC = new TF1("EMVconv",*ConvEMG,-300,r2+500,ConvEMG->GetNpar());
  FC->SetParameters(1,1,1);
  FC->SetNormalized(kTRUE);

  //Normalize & Exclude
  TF1 F("EMV",[&](Double_t* x, Double_t* p){return FFNormAmpNC(a,c,FC,x,p);},r1,r2,ConvEMG->GetNpar()+a+c);
  F.SetNpx(2000); F.SetFillColorAlpha(kAzure-9,0.5); F.SetFillStyle(1001);
  F.SetLineColor(kGreen); F.SetLineWidth(1);
  int i=-1;
  i++; F.SetParName(i,"Constant"); F.SetParameter(i,1E5);
  while(i<a-1){
  i++; F.SetParName(i,"R_Escape"); F.SetParameter(i,0.1);       F.SetParLimits(i,0,1);
  }
  i++; F.SetParName(i,"Centroid"); F.SetParameter(i,(r1+r2)/2);
  while(i<a+c-1){
  i++; F.SetParName(i,"Offset");   F.SetParameter(i,-5);        F.SetParLimits(i,-10,0);
  }
  i++; F.SetParName(i,"Sigma");    F.SetParameter(i,3);         F.SetParLimits(i,0,4);
  i++; F.SetParName(i,"Tau");      F.SetParameter(i,8);         F.SetParLimits(i,0,15);
  i++; F.SetParName(i,"Gamma");    F.SetParameter(i,0.04);      F.SetParLimits(i,0.02,0.06);
  return F;
}

TF1 FDSAM(int &n, TF1* & FN){ //DSAM Function
  FN = new TF1("DSAMfn",FFDSAM0,-300,r2+250,3); //Apply Norm and Amp
  FN->SetNpx(500); FN->SetParameters(1,1);
  FN->SetNormalized(kTRUE);
  TF1 F("DSAM",[&](Double_t* x, Double_t* p){return FFNormAmpN(n,FN,x,p);},r1,r2,3+n);
	F.SetNpx(2000); F.SetFillColorAlpha(kAzure-9,0.5); F.SetFillStyle(1001);
  F.SetLineColor(kGreen); F.SetLineWidth(1); F.SetLineStyle(3);
  int i=-1;
  i++; F.SetParName(i,"Constant"); F.SetParameter(i,5e5); F.SetParLimits(i,1e5,5e6);
  while(i<n-1){
  i++; F.SetParName(i,"B_Ratio");  F.SetParameter(i,0.1); F.SetParLimits(i,0,1);
  }
  i++; F.SetParName(i,"Decay");    F.SetParameter(i,2);   F.SetParLimits(i,1.7,2.7);
  i++; F.SetParName(i,"Energy");    F.SetParameter(i,29);   F.SetParLimits(i,10,32);
  i++; F.SetParName(i,"AOffset");   F.SetParameter(i,59);  F.SetParLimits(i,50,100);

  return F;
}

TF1 FDSAMConvolved(int &n, TF1* & FC, TF1* &FX){
  //Convolve
  TF1* FD = new TF1("DSAM0",FFDSAM0,-300,r2+250,3); //Apply Norm and Amp
  FD->SetNpx(500);
  TF1Convolution ConvDSAM(FD,FX,-300,r2+250,true);
  ConvDSAM.SetRange(-300,r2+250); //Implements AbsComposition?
  ConvDSAM.SetNofPointsFFT(1e5); //FFT Points?
  FC = new TF1("DSAMconv",ConvDSAM,-300,r2+250,ConvDSAM.GetNpar());
  FC->SetParameters(1,1,1,1,1);
  FC->SetNormalized(kTRUE); //Needs a pass to initialize parameters properly

  //Normalize & Exclude
  TF1 F("DSAMC",[&](Double_t* x, Double_t* p){return FFNormAmpN(n,FC,x,p);},r1,r2,ConvDSAM.GetNpar()+n);
  F.SetNpx(500); F.SetFillColorAlpha(kAzure-9,0.5); F.SetFillStyle(1001);
  F.SetLineColor(kGreen); F.SetLineWidth(1);
  int i=-1;
  i++; F.SetParName(i,"Constant"); F.SetParameter(i,5.1e5);  F.SetParLimits(i,10,5e8);
  while(i<n-1){
  i++; F.SetParName(i,"B_Ratio");  F.SetParameter(i,0.1049); F.SetParLimits(i,0.1040,0.1056); //Branching Ratio: 10.44(4)%
  }
  i++; F.SetParName(i,"Decay");    F.SetParameter(i,1.6);    F.SetParLimits(i,1,2.5);
  i++; F.SetParName(i,"Energy");   F.SetParameter(i,29);     F.SetParLimits(i,10,32);
  i++; F.SetParName(i,"AOffset");   F.SetParameter(i,59);     F.SetParLimits(i,0,100);
  i++; F.SetParName(i,"Sigma");    F.SetParameter(i,3);      F.SetParLimits(i,0,10);
  return F;
}

TF1 FDSAMG(int &n,TF1* & FN){ //Convolved with Gaussian
  TF1* FG = new TF1("G0",FFGaus0,-300,r2+250,1);
  FG->SetNpx(500);
  TF1 F = FDSAMConvolved(n,FN,FG);
  F.SetNameTitle("DSAMG","DSAMG");
  return F;
}

TF1 FDSAMV(int &n,TF1* & FC){ //Convolved with Voigt
  TF1* FV = new TF1("Voigt0",FFVoigt0,-300,r2+250,2);
  FV->SetNpx(500);
  TF1 F = FDSAMConvolved(n,FC,FV);
  F.SetNameTitle("DSAMV","DSAMV");

  F.SetParName(6,"Gamma");    F.SetParameter(6,0.04);   F.SetParLimits(6,0.02,0.06);
  return F;
}

TF1 FTotal(Double_t r1, Double_t r2, map<int,vector<int>> & parVec){ //Combines pre-defined Peak Functions
  int ptotal = get_max(parVec)+1;
  vector<TF1*> tempFits = Fits;
  TF1 FT("Total_Function",[&](Double_t* x, Double_t* p){return FFTotal(x,p,parVec);},r1,r2,ptotal);
  FT.SetNpx(750); FT.SetLineWidth(2); //FT.SetFillStyle(1001);
  //FT.SetFillColorAlpha(kRed,0.4);

  //Convert parVec to names/parameters
  vector<int> psofar;
    for(auto x : parVec){
      int i = x.first;
      for(size_t j=0;j<parVec[i].size();j++){
        if(psofar.end()==find(psofar.begin(),psofar.end(),parVec[i][j])){
          const char* nm = FTParName(i,tempFits[i]->GetParName(j)).c_str();
          FT.SetParName(parVec[i][j],nm);
          double a,b;
          tempFits[i]->GetParLimits(j,a,b);
          FT.SetParLimits(parVec[i][j],a,b);
          FT.SetParameter(parVec[i][j],tempFits[i]->GetParameter(j));
          psofar.push_back(parVec[i][j]);
        }
      }
    }

  return FT;
}

//Useful Macros_________________________________________________________________
//Rebinning_____________________________________________________________________
//Rebin by width
void WidthRebin(TH1* H1, double bwdth){
  auto f = [&](double a){return bwdth*a;};
  const static Int_t nbins = H1->GetXaxis()->GetNbins();
  double new_bins[nbins+1];
  for(int i=0; i <= nbins; i++){
    new_bins[i] = f(H1->GetBinLowEdge(i+1));
  }
  H1->SetBins(nbins, new_bins);
}

//Linear Rebin
void LinRebin(TH1* H1, vector<vector<double>> pts){
  double x1=pts[0][0]; double y1 = pts[0][1];
  double x2=pts[1][0]; double y2 = pts[1][1];
    double m =(y1-y2)/(x1-x2);
  auto f = [&](double x){return m*(x-x1)+y1;};
  const static Int_t nbins = H1->GetXaxis()->GetNbins();
  double new_bins[nbins+1];
    for(int i=0; i <= nbins; i++){
      new_bins[i] = f(H1->GetBinLowEdge(i+1));
    }
  H1->SetBins(nbins, new_bins);
}

//Plot Modules__________________________________________________________________
void RegLinePlot(TF1* F){//Draw fit regions as well as function
  //Draw just the filling
  F->SetLineWidth(0); F->SetFillColorAlpha(kAzure-9,0.5);
  F->SetFillStyle(1001); F->DrawClone("Same");

  //Draw the line while clearing FitExcl & setting range
  F->SetFillStyle(0); FitExcl.clear(); F->SetRange(r1,r2);
  F->SetLineWidth(1); F->DrawClone("Same");
}

//Output Modules________________________________________________________________
void CanvBasic(TCanvas* C,TF1* F,TH1F* H,TH1F* R2){
  //Setup
  TGaxis::SetExponentOffset(0.02,-0.075,"y");
  C->Clear(); C->cd();
  TPad *fpad = new TPad("fpad","Fit",0,0.33,1,1);
  TPad *rpad = new TPad("rpad","Residuals",0,0,1,0.33);
  rpad->SetFillStyle(4000);
  fpad->Draw(); rpad->Draw();
  fpad->SetNumber(1); rpad->SetNumber(2);
  fpad->SetLogy(); //rpad->SetLogy();
  F->SetRange(r1,r2); //May screw up Fits. Need for plotting change of r2
  H->SetAxisRange(r1,r2,"X");

  //Fit Region Residuals
  TH1F* R = (TH1F*)H->Clone("Residuals"); R->GetListOfFunctions()->Clear();
  for(int i=0;i<R->GetNbinsX();i++){
    double x = R->GetBinLowEdge(i);
    double a,b; F->GetRange(a,b);
    R->SetBinError(i,R->GetBinError(i)/R->GetBinContent(i)); //Normalized
    R->SetBinContent(i,(R->GetBinContent(i)-F->Eval(x))/R->GetBinContent(i)); //Normalized
    if(ExclRange(x)||x<a||x>b){R->SetBinContent(i,0); R->SetBinError(i,0);}
  }
  double ymin = R->GetBinContent(R->GetMinimumBin());
  double ymax = R->GetBinContent(R->GetMaximumBin());
  R->SetAxisRange(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin),"Y");
  double boff = 30*(r2-r1)/(R->GetXaxis()->GetLast()-R->GetXaxis()->GetFirst());

  //Plot Fit graphic
  fpad->cd();
  fpad->SetRightMargin(0.05); fpad->SetLeftMargin(0.05);
  fpad->SetBottomMargin(0); fpad->SetTopMargin(0.1);
  H->SetAxisRange(r1-boff,r2,"X"); H->SetAxisRange(10,1E7,"Y");
  H->SetLineColor(kBlack); H->SetName("Histo1");
  H->GetYaxis()->SetTitle("Counts"); H->GetYaxis()->SetTitleOffset(0.7);
  H->GetYaxis()->SetTickSize(0.015);
  H->GetXaxis()->SetLabelSize(0); H->GetYaxis()->SetLabelSize(0.04);
  H->Draw();
  RegLinePlot(F); //Draw and reset FitExcl
  gStyle->SetStatW(0.13);
  gStyle->SetStatH(0.1);
  gStyle->SetStatX(1);
  gStyle->SetStatY(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  for(auto i:Fits){//Draw each fit function too
    double a,b; i->GetRange(a,b); i->SetRange(r1,r2);
    if(i->GetLineStyle()==1){i->SetLineStyle(2);} i->SetLineColor(kBlack); i->SetFillStyle(0);
    i->DrawClone("Same"); i->SetRange(a,b);
  }
  //Pass Subtracted Plot
  R2->Add(F,-1); //Return next function for fitting
  //Plot Resiudal graphic
  rpad->cd(); R->SetTitle("");
  rpad->SetRightMargin(0.05);
  rpad->SetLeftMargin(0.05);
  rpad->SetBottomMargin(0.2);
  rpad->SetTopMargin(0);
  rpad->SetTickx(1);
  R->SetAxisRange(r1-boff,r2,"X");
  R->GetXaxis()->SetTitle("Energy (eV)");
  R->GetXaxis()->SetTitleOffset(0.9);
  R->GetXaxis()->SetTitleSize(0.08);
  R->GetYaxis()->SetTitle("");
  R->GetYaxis()->SetTickSize(0.015);
  R->SetLineColor(kBlack);
  R->GetXaxis()->SetLabelSize(0.08);
  R->GetYaxis()->SetLabelSize(0.06);
  R->Draw();
  TPaveStats* STATS = (TPaveStats*)H->FindObject("stats");
  C->GetListOfPrimitives()->AddLast(STATS); C->Modified();
}

void CanvE(TCanvas* &C,TF1* &F,TH1F* &H,TH1F* &R2, TString nm){
  CanvBasic(C,F,H,R2);
  C->SaveAs("./Output/E_Sums/" + nm + ".pdf");
  H->GetListOfFunctions()->Clear();
}

void CanvG(TCanvas* C,TF1* F,TH1F* H,TH1F* R2, TString nm){
  CanvBasic(C,F,H,R2);
  C->SaveAs("./Output/G_Final/" + nm + ".pdf");
  H->GetListOfFunctions()->Clear();
}

void CanvF(TCanvas* &C,TF1* &F,TH1F* &H,TH1F* &R2, TString nm){
  CanvBasic(C,F,H,R2);
  C->SaveAs("./Output/F_Individuals/" + nm + ".pdf");
  H->GetListOfFunctions()->Clear();
}

//Miscellaneous_________________________________________________________________
void ParCopy(TF1* & F){//Copy parameters from FT back to Fits
  map<int,vector<int>> parVec = FTParVec();
  for(auto i : parVec){
    int ii=i.first;
    int jj=0;
    for(auto j : i.second){
      double a,b; F->GetParLimits(j,a,b);
      double pj = F->GetParameter(j);
      Fits[ii]->SetParameter(jj,pj);
      if(a!=b){
        if(TMath::Abs(1-pj/a)<=0.001||TMath::Abs(1-pj/b)<=0.001){
          cout<<"Parameter "<<j<<"("<<F->GetParName(j)<<")"<<
          " is likely bounded."<<endl;
          cout<<a<<"< P("<<j<<") = "<<pj<<" < "<<b<<endl;
        }
      // Fits[ii]->SetParLimits(jj,pj-0.5*pj,pj+0.5*pj);
      }
      jj++;
    }
  }
}

void ParLimit(TF1* &F, int i, double l){//Limit parameter to fraction l
    double a = F->GetParameter(i);
    F->SetParLimits(i,a-abs(l*a),a+abs(l*a));
}

void FixSame(TF1* fin, const char* pname, TF1* ffix){//Fix parameter from ffix to fin
    int pnum = fin->GetParNumber(pname);
    fin->FixParameter(pnum,ffix->GetParameter(pname));
}

void SetSame(TF1* fin, const char* pname, TF1* ffix){//Fix parameter from ffix to fin
    int pnum = fin->GetParNumber(pname);
    fin->SetParameter(pnum,ffix->GetParameter(pname));
}

void ParLimits(TF1* & F, double l){//Limit parameters to fraction l
  for(int i=0; i<F->GetNpar(); i++){
    double a = F->GetParameter(i);
    F->SetParLimits(i,a-abs(l*a),a+abs(l*a));
  }
}

void FixName(TF1* fin, const char* pname, Double_t val){//Fix parameter by NAME
  fin->FixParameter(fin->GetParNumber(pname),val);
  if(val==0){
    fin->SetParameter(fin->GetParNumber(pname),-1);
    fin->SetParLimits(fin->GetParNumber(pname),-1,-1);
  }
}

void FTInitPars(vector<TF1*> & IndFits, int fskip){
  for(unsigned i=0;i<Fits.size();i++){
    double p[10]; IndFits[i+fskip]->GetParameters(p);
    Fits[i]->SetParameters(p);
    const char* fname = IndFits[i+fskip]->GetTitle();
    Fits[i]->SetTitle(fname);
    // cout << fname << ": Fits[" << i << "] - " << p << endl;
    for(int j=0; j<IndFits[i+fskip]->GetNpar(); j++){
      double a,b; IndFits[i+fskip]->GetParLimits(j,a,b);
      Fits[i]->SetParLimits(j,a,b);
      const char* pname=IndFits[i+fskip]->GetParName(j);
      Fits[i]->SetParName(j,pname);
      // cout << "(" << j << ") - [" << a << "," << b << "] : " << pname << endl;
    }
  }
}

unsigned FitNum(const char* pname){
  unsigned i=0;
  while(i<Fits.size()){
    if(strcmp(pname,Fits[i]->GetName())==0){return i;}
    i++;
  }
  cout << "ERROR: Function (" << pname << ") NOT FOUND" << endl;
  return Fits.size();
}

void ParSave(vector<double> & p){//Save Fits parameters to a growing list for total function building
  for(unsigned i=0;i<Fits.size();i++){
    vector<double> ptemp; ptemp.resize(Fits[i]->GetNpar());
    Fits[i]->GetParameters(&*ptemp.begin());
    p.insert(p.end(),ptemp.begin(),ptemp.end());
  }
}

void UpdateParLims(TF1* FT){ //Update parameter limits after parameter merging, shouldn't matter...
  map<int,vector<int>> parVec = FTParVec();
  int ALT = 0;
  double a1,b1;
  double a2,b2;
  for(unsigned i=0; i<parVec.size(); i++){//Update fits in parVec if changed
    for(int j=0;j<Fits[i]->GetNpar();j++){
      FT->GetParLimits(parVec[i][j],a1,b1); //Merged FT parameter limits
      Fits[i]->SetParLimits(j,a1,b1); //Set limits of original fit in case merged
    }
  }
}

void FTCheck(TF1* & FT, TCanvas* CF){ //Half each parameter and redraw FT
  string name = "START";
  for(int i=0; i<FT->GetNpar(); i++){
    CF->Clear(); CF->cd(); FT->SetTitle(name.c_str());
    FT->Draw(); CF->SetLogy();
    for(unsigned j=0; j<Fits.size(); j++){
      Fits[j]->Draw("Same");
    }
    CF->SaveAs(((string)"Output/G_Final/"+name+(string)".pdf").c_str());
    name = to_string(i);
    FT->SetParameter(i,FT->GetParameter(i)*0.5);
  }
}

void PrintParameters(TF1* & F){
  cout << F->GetName() << endl;
  for(int i=0;i<F->GetNpar();i++){
    cout << "  (" << i << ") " << F->GetParName(i);
    double a,b; F->GetParLimits(i,a,b); if(b<=F->GetParameter(i)||F->GetParameter(i)<=a){cout << " - [ " << a << " < " << F->GetParameter(i) << " < " << b << " ]" << endl;}
    else{ cout << " - [ " << F->GetParameter(i) << " ]" << endl;}
  }
}
