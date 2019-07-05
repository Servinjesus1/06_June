// A simple import of Laser Data to ROOT branches and Histograms
// Author: Spencer Fretwell
// Date: 1/14/2019

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TTree.h"

//_____________________________
//Variables

const char *f_in = "Run3.txt";
const char *f_out = "Run3.root";

//_____________________________
//Importing Data

void RTImport(){
  //Import Data from txt file to TFile
  TFile* LDFile = new TFile(f_out,"RECREATE","Laser Test Data from Experiment");
  TTree* tree = new TTree("Data","Main Data Tree");
    //Silence Read Errors (empty bins)
    Int_t currentIgnoreLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError;
    Long64_t nlines = tree->ReadFile(f_in,"",'\t');
    //Return to Default
    gErrorIgnoreLevel= currentIgnoreLevel;
    LDFile->Write();

//___________________________
//Making the histograms

  TTreeReader reader("Data", LDFile); //Tree Reader Object

  //Histograms
  Int_t NumBranches = tree->GetListOfBranches()->GetEntries(); //Number of Branches

  for(int it=1;it<NumBranches;it++){ //Branch For Loop
    //Filling Histograms (2)
    int jt=1;
    LDFile->GetObject("Data",tree); //Retrieves tree from File
    const char* BrName = tree->GetListOfBranches()->At(it)->GetName(); //Branch Name String
    TH1D h = TH1D(Form("%s",BrName),Form("%s",BrName),nlines,0,nlines); //Recreate Histogram h#
    TTreeReaderValue<Float_t> Entries(reader,BrName);
    while(reader.Next()){ //Event For Loop
      for(int i=0; i<*Entries; i++) h.Fill(jt);
      jt++;
    }
    delete tree; //Deletes tree from memory so it doesn't duplicate
    LDFile->Write();
    reader.Restart();
  }
}
