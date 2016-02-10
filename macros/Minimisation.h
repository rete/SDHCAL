//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat May 31 15:48:59 2014 by ROOT version 5.34/05
// from TTree CaloHitAnalysis/CaloHitAnalysis
// found on file: ../output/ROOTOutputFile_single_pi-_10GeV.root
//////////////////////////////////////////////////////////

#ifndef Minimisation_h
#define Minimisation_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector>
#include <cmath>

class MinimizerData
{
public:
	std::vector<int>     m_nHit;
	std::vector<int>     m_nHit1;
	std::vector<int>     m_nHit2;
	std::vector<int>     m_nHit3;
	std::vector<float>   m_energy;
	std::vector<double>  m_bestParameters;
};

MinimizerData *pMinimizerData;

std::string getFileName(const std::string &fileLocation, int energy);

double minimisationFunction(double x1, double y1, double z1,
	      double x2, double y2, double z2,
	      double x3, double y3, double z3);

//double minimisationFunction(double w0, double w1, double w2);

void minuitFunction(int& nDim, double *gout, double& result, double par[], int flg);

double resolutionFunction(double *x, double *par);

void processMinimizer(int nEpoch = 1, const std::string &fileLocation = "");

void minimize(int nEpoch);

class Minimisation;

double getReconstructedEnergy(Minimisation *pMinimisation, MinimizerData *pMinimizerData);


class Minimisation
{

public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           NHit;
   Int_t           NHit1;
   Int_t           NHit2;
   Int_t           NHit3;

   // List of branches
   TBranch        *b_NHit;   //!
   TBranch        *b_NHit1;   //!
   TBranch        *b_NHit2;   //!
   TBranch        *b_NHit3;   //!

   Minimisation(TTree *tree=0);
   virtual ~Minimisation();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   TH1 *DrawEnergyHistogram();

   MinimizerData     *m_pMinimizerData;
   float             m_currentEnergy;
   float             m_meanEnergyFit;
   float             m_meanEnergyErrorFit;
   float             m_energyResolutionFit;
   float             m_energyResolutionErrorFit;
   TH1               *m_pEnergyHistogram;
};

#endif

#ifdef Minimisation_cxx
Minimisation::Minimisation(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0)
   {
   	throw;
   }
   Init(tree);
}

Minimisation::~Minimisation()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Minimisation::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Minimisation::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Minimisation::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("NHit", &NHit, &b_NHit);
   fChain->SetBranchAddress("NHit1", &NHit1, &b_NHit1);
   fChain->SetBranchAddress("NHit2", &NHit2, &b_NHit2);
   fChain->SetBranchAddress("NHit3", &NHit3, &b_NHit3);

   Notify();
}

Bool_t Minimisation::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Minimisation::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Minimisation::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Minimisation_cxx
