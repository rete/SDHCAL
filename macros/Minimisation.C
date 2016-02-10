#define Minimisation_cxx
#include "Minimisation.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TF1.h>
#include <TText.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TFitter.h>
#include <TH2D.h>

#include <sstream>
#include <exception>
#include <stdexcept>
#include <cmath>

enum DataType
{
	SIMULATION,
	TEST_BEAM
};

std::string GetStringEnergy(DataType dataType, int energy)
{
	if(SIMULATION == dataType)
	{
		std::stringstream ss;
		ss << energy;
		return ss.str();
	}
	else
	{
		switch (energy)
		{
		 case 10:
		 	return "715693";
		 case 20:
		 	return "715675";
		 case 30:
		 	return "715747";
		 case 40:
		 	return "715748";
		 case 50:
		 	return "715751";
		 case 60:
		 	return "715753";
		 case 70:
		 	return "715754";
		 case 80:
		 	return "715756";
		 default:
		 	throw std::invalid_argument("Invalid energy");
		}
	}
}


//--------------------------------------------------------------------------------

void Minimisation::Loop()
{
 if (fChain == 0)
 	return;

	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	for (Long64_t jentry=0; jentry<nentries;jentry++)
	{
		Long64_t ientry = LoadTree(jentry);

		if (ientry < 0)
			break;

		nb = fChain->GetEntry(jentry);
		nbytes += nb;

		m_pMinimizerData->m_nHit.push_back(NHit);
		m_pMinimizerData->m_nHit1.push_back(NHit1);
		m_pMinimizerData->m_nHit2.push_back(NHit2);
		m_pMinimizerData->m_nHit3.push_back(NHit3);
		m_pMinimizerData->m_energy.push_back(m_currentEnergy);
	}
}



double CrystalBall(double* x, double* par)
{
 //http://en.wikipedia.org/wiki/Crystal_Ball_function
 double xcur = x[0];
 double alpha = par[0];
 double n = par[1];
 double mu = par[2];
 double sigma = par[3];
 double N = par[4];
 TF1* exp = new TF1("exp","exp(x)",1e-20,1e20);
 double A; double B;
 if (alpha < 0){
 A = pow((n/(-1*alpha)),n)*exp->Eval((-1)*alpha*alpha/2);
 B = n/(-1*alpha) + alpha;}
 else {
 A = pow((n/alpha),n)*exp->Eval((-1)*alpha*alpha/2);
 B = n/alpha - alpha;}

 double f;
	if ((xcur-mu)/sigma > (-1)*alpha)
	f = N*exp->Eval((-1)*(xcur-mu)*(xcur-mu)/
(2*sigma*sigma));
	else
	f = N*A*pow((B- (xcur-mu)/sigma),(-1*n));
	delete exp;

  return f;
}

//--------------------------------------------------------------------------------

TH1 *Minimisation::DrawEnergyHistogram()
{
 if (fChain == 0)
 	return NULL;

	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	std::stringstream histoStrStream;
	histoStrStream << "h_" << m_currentEnergy << "GeV";

	m_pEnergyHistogram = new TH1D(histoStrStream.str().c_str(), histoStrStream.str().c_str(), 100, 0, 2*m_currentEnergy);

	for (Long64_t jentry=0; jentry<nentries;jentry++)
	{
		Long64_t ientry = LoadTree(jentry);

		if (ientry < 0)
			break;

		nb = fChain->GetEntry(jentry);
		nbytes += nb;

		double reconstructedEnergy = getReconstructedEnergy(this, m_pMinimizerData);

		m_pEnergyHistogram->Fill(reconstructedEnergy);
	}

 TF1 *pGaussianFunction = new TF1("gausFunc", "gaus", 0, 120);
 gStyle->SetOptFit(1);
 m_pEnergyHistogram->Fit(pGaussianFunction, "NQ", "");
 m_pEnergyHistogram->Fit(pGaussianFunction, "Q", "",
			 pGaussianFunction->GetParameter(1)-1.5*pGaussianFunction->GetParameter(2),
			 pGaussianFunction->GetParameter(1)+1.5*pGaussianFunction->GetParameter(2) );
 m_pEnergyHistogram->DrawCopy();

 double eRec = pGaussianFunction->GetParameter(1);
 double eRecError = pGaussianFunction->GetParError(1);
 double eResol = pGaussianFunction->GetParameter(2)/pGaussianFunction->GetParameter(1);
 double eResolError = std::sqrt(std::pow(pGaussianFunction->GetParError(2)/eRec, 2) + std::pow((eResol*eRecError)/(eRec*eRec), 2));

 m_meanEnergyFit            = eRec;
 m_energyResolutionFit      = eResol;
 m_meanEnergyErrorFit       = eRecError;
 m_energyResolutionErrorFit = eResolError;

 std::cout << "Beam energy = " << m_currentEnergy << ", "
	   << "ERec = " << eRec << " +- " << eRecError << " GeV , "
	   << "EResol = "  << eResol << " +- " << eResolError
	   << std::endl;

 TText *pText = new TText();
 std::stringstream ss;
 ss << "E_{beam} = " << m_currentEnergy << " GeV";
 pText->DrawTextNDC(0.2, 0.8, ss.str().c_str());

 return m_pEnergyHistogram;
}

//--------------------------------------------------------------------------------

double getReconstructedEnergy(Minimisation *pMinimisation, MinimizerData *pMinimizerData)
{
	return 	(pMinimizerData->m_bestParameters.at(0) + pMinimizerData->m_bestParameters.at(1)*pMinimisation->NHit + pMinimizerData->m_bestParameters.at(2)*pMinimisation->NHit*pMinimisation->NHit)*pMinimisation->NHit1
	    		+ (pMinimizerData->m_bestParameters.at(3) + pMinimizerData->m_bestParameters.at(4)*pMinimisation->NHit + pMinimizerData->m_bestParameters.at(5)*pMinimisation->NHit*pMinimisation->NHit)*pMinimisation->NHit2
	    		+ (pMinimizerData->m_bestParameters.at(6) + pMinimizerData->m_bestParameters.at(7)*pMinimisation->NHit + pMinimizerData->m_bestParameters.at(8)*pMinimisation->NHit*pMinimisation->NHit)*pMinimisation->NHit3;
}

//--------------------------------------------------------------------------------

double minimisationFunction(double w0, double w1, double w2,
	      double w3, double w4, double w5,
	      double w6, double w7, double w8)
{
  unsigned int size = pMinimizerData->m_nHit.size();

  double test = 0;
  double energyReco = 0.f;

  for(unsigned int i=0; i<size; i++)
  {
  	energyReco =
				(w0 + w1*pMinimizerData->m_nHit.at(i) + w2*pMinimizerData->m_nHit.at(i)*pMinimizerData->m_nHit.at(i))*pMinimizerData->m_nHit1.at(i)+
			 (w3 + w4*pMinimizerData->m_nHit.at(i) + w5*pMinimizerData->m_nHit.at(i)*pMinimizerData->m_nHit.at(i))*pMinimizerData->m_nHit2.at(i)+
			 (w6 + w7*pMinimizerData->m_nHit.at(i) + w8*pMinimizerData->m_nHit.at(i)*pMinimizerData->m_nHit.at(i))*pMinimizerData->m_nHit3.at(i);

			test +=	std::pow(pMinimizerData->m_energy.at(i) - energyReco, 2) / pMinimizerData->m_energy.at(i);
  }
  return test / (size-1);
}

//--------------------------------------------------------------------------------

void minuitFunction(int& nDim, double* gout, double& result, double *par, int flg)
{
	result = minimisationFunction(par[0], par[1], par[2],
			par[3], par[4], par[5],
			par[6], par[7], par[8]);
}

//--------------------------------------------------------------------------------

std::string getFileName(const std::string &fileLocation, int energy)
{
	std::stringstream ss;
	ss << fileLocation << "ROOTOutputFile_single_pi-_" << energy << "GeV.root";
//	ss << fileLocation << "/TDHCAL_minimisation_" << GetStringEnergy(TEST_BEAM, energy) << "_cut.root";
	return ss.str();
}

//--------------------------------------------------------------------------------

double resolutionFunction(double *x, double *par)
{
  return std::sqrt(par[0]*par[0]/x[0] + par[1]*par[1]);
}

//--------------------------------------------------------------------------------

void processMinimizer(int nEpoch, const std::string &fileLocation)
{
	if(nEpoch <= 0)
	{
		std::cerr << "Invalid number of epoch for minimization. Given : " << nEpoch << ". Expected > 0" << std::endl;
		return;
	}

	std::cout << "Processing minimizer with nEpoch = " << nEpoch << std::endl;

	pMinimizerData = new MinimizerData();
	// please feel free to modify the directory
	//	const std::string fileLocation = "/home/rete/soft/SDHCAL/output/minimisation/FTF_BIC/";

	// /home/rete/soft/SDHCAL/output/minimisation/FTFP_BERT_HP/
	// /home/rete/soft/SDHCAL/output/minimisation/FTF_BIC/
	// /home/rete/soft/SDHCAL/output/minimisation/TB_SeptAug2012/

	
	const std::string treeName = "CaloHitAnalysis";

	std::vector<int> energies;
	energies.push_back(10);
	energies.push_back(20);
	energies.push_back(30);
	energies.push_back(40);
	energies.push_back(50);
	energies.push_back(60);
	energies.push_back(70);
	energies.push_back(80);

	//===== grab all the data in one vector =====

	for(unsigned int e=0 ; e<energies.size() ; e++)
	{
		TFile *pFile = TFile::Open(getFileName(fileLocation, energies.at(e)).c_str());

		if(!pFile)
			throw runtime_error("Invalid file name !");

		pFile->ls();

		TTree *pTree = (TTree *) pFile->Get(treeName.c_str());

		std::cout << "TTree address : " << pTree << std::endl;

		if(!pTree)
			throw runtime_error("Invalid tree name !");

		std::cout << "Taking data from " << energies.at(e) << " GeV file ..." << std::endl;

		Minimisation *pMinimisation = new Minimisation(pTree);
		pMinimisation->m_currentEnergy = energies.at(e);
		pMinimisation->m_pMinimizerData = pMinimizerData;

		pMinimisation->Loop();

		delete pMinimisation;
	}

	std::cout << "Data size to minimize : " << pMinimizerData->m_nHit.size() << std::endl;
	minimize(nEpoch);

	//============

 TGraph *pReconstructedEnergyGraph = new TGraph();
 pReconstructedEnergyGraph->SetMarkerStyle(20);
 pReconstructedEnergyGraph->SetMarkerSize(0.8);

 TGraphErrors *pEnergyResolutionGraph = new TGraphErrors();
 pEnergyResolutionGraph->SetMarkerStyle(20);
 pEnergyResolutionGraph->SetMarkerSize(0.8);

 TCanvas *cc;

	for(unsigned int e=0 ; e<energies.size() ; e++)
 {
		TFile *pFile = TFile::Open(getFileName(fileLocation, energies.at(e)).c_str());
		TTree *pTree = (TTree *) pFile->Get(treeName.c_str());

		Minimisation *pMinimisation = new Minimisation(pTree);
		pMinimisation->m_currentEnergy = energies.at(e);
		pMinimisation->m_pMinimizerData = pMinimizerData;

	 gStyle->SetOptFit(1);
	 std::stringstream canvasStrStream;
	 canvasStrStream << energies.at(e) << "GeVCanvas";

	 std::cout << "Drawing histogram" << std::endl;
	 cc = new TCanvas(canvasStrStream.str().c_str(), canvasStrStream.str().c_str());
	 cc->SetWindowSize(600, 600);
	 cc->cd();
	 pMinimisation->DrawEnergyHistogram();//->Write();
	 cc->Update();
	 std::cout << "Drawing histogram OK" << std::endl;

	 double eRec = pMinimisation->m_meanEnergyFit;
	 double eRecError = pMinimisation->m_meanEnergyErrorFit;
	 double eResol = pMinimisation->m_energyResolutionFit;
	 double eResolError = pMinimisation->m_energyResolutionErrorFit;

	 pReconstructedEnergyGraph->SetPoint(e, energies.at(e), eRec);
	 pEnergyResolutionGraph->SetPoint(e, energies.at(e), eResol);
	 pEnergyResolutionGraph->SetPointError(e, 0, eResolError);
 }

	// draw the results !
	TCanvas *pCanvas = new TCanvas(fileLocation.c_str(), fileLocation.c_str());
 pCanvas->SetWindowSize(1000, 600);
 pCanvas->Divide(2, 1);

 pCanvas->cd(1);
 TF1 *pLinearFunction = new TF1("LinearFunction", "[0]*x", 0, 100);
 pReconstructedEnergyGraph->Draw("ap");
 pReconstructedEnergyGraph->Fit(pLinearFunction, "", "", 0, 75);

 pCanvas->cd(2);
 TF1 *pEnergyResolutionFunction = new TF1("ResolutionFunction", resolutionFunction, 0, 50, 3);
 pEnergyResolutionFunction->SetLineColor(1);
 pEnergyResolutionFunction->SetLineWidth(2);
 pEnergyResolutionFunction->SetParameter(0,0.6);
 pEnergyResolutionFunction->SetParameter(1,0.06);
 pEnergyResolutionFunction->SetParameter(2,0.);
 pEnergyResolutionGraph->Draw("ap");
 pEnergyResolutionGraph->Fit(pEnergyResolutionFunction);

	delete pMinimizerData;
}

//--------------------------------------------------------------------------------

void minimize(int nEpoch)
{
	std::vector<double> bestParameters;

	// these are by current parameters
	// bestParameters.push_back(0.0442674);
	// bestParameters.push_back(3.95986e-05);
	// bestParameters.push_back(-2.0133e-10);
	// bestParameters.push_back(0.0781586);
	// bestParameters.push_back(-5.9623e-05);
	// bestParameters.push_back(-3.06952e-08);
	// bestParameters.push_back(0.115295);
	// bestParameters.push_back(1.57508e-05);
	// bestParameters.push_back(2.28041e-09);

	// from Imad FTFP_BERT_HP
	bestParameters.push_back( 0.0293659   );
	bestParameters.push_back( 2.96392e-05 );
	bestParameters.push_back(-2.051233e-08);
	bestParameters.push_back( 0.0925993   );
	bestParameters.push_back( 1.12284e-05 );
	bestParameters.push_back(-2.016464e-09);	
	bestParameters.push_back( 0.166332    );
	bestParameters.push_back( 3.08655e-05 );
	bestParameters.push_back( 2.97939e-08 );
	

	std::vector<double> tryParameters(bestParameters);

	double minimum = minimisationFunction(bestParameters[0], bestParameters[1], bestParameters[2],
 		bestParameters[3], bestParameters[4], bestParameters[5],
 		bestParameters[6], bestParameters[7], bestParameters[8]);

 std::cout << "Minimum initialized to : " << minimum << std::endl;

 for(int epoch=0 ; epoch < nEpoch ; epoch++)
 {
		TFitter* pMinimizer = new TFitter(9);

		double p1 = -1;
		pMinimizer->ExecuteCommand("SET PRINTOUT",&p1,1);
		pMinimizer->SetFCN(minuitFunction);
		pMinimizer->SetParameter(0,"W0",tryParameters[0],.0001,0,0);
		pMinimizer->SetParameter(1,"W1",tryParameters[1],.0001,0,0);
		pMinimizer->SetParameter(2,"W2",tryParameters[2],.0001,0,0);
		pMinimizer->SetParameter(3,"W3",tryParameters[3],.0001,0,0);
		pMinimizer->SetParameter(4,"W4",tryParameters[4],.00001,0,0);
		pMinimizer->SetParameter(5,"W5",tryParameters[5],.000001,0,0);
		pMinimizer->SetParameter(6,"W6",tryParameters[6],.0001,0,0);
		pMinimizer->SetParameter(7,"W7",tryParameters[7],.00001,0,0);
		pMinimizer->SetParameter(8,"W8",tryParameters[8],.00000001,0,0);


		pMinimizer->ExecuteCommand("SIMPLEX",0,0);
//		minimizer->ExecuteCommand("MIGRAD",0,0);
//		minimizer->ExecuteCommand("MINOS",0,0);

		std::vector<double> newParameters(9, 0.);
		newParameters[0] = pMinimizer->GetParameter(0);
		newParameters[1] = pMinimizer->GetParameter(1);
		newParameters[2] = pMinimizer->GetParameter(2);
		newParameters[3] = pMinimizer->GetParameter(3);
		newParameters[4] = pMinimizer->GetParameter(4);
		newParameters[5] = pMinimizer->GetParameter(5);
		newParameters[6] = pMinimizer->GetParameter(6);
		newParameters[7] = pMinimizer->GetParameter(7);
		newParameters[8] = pMinimizer->GetParameter(8);

		double newMinimum = minimisationFunction(
				newParameters[0], newParameters[1], newParameters[2],
				newParameters[3], newParameters[4], newParameters[5],
				newParameters[6], newParameters[7], newParameters[8] );

		std::cout << "Epoch : " << epoch << ", new minimum : " << newMinimum << std::endl;

		// if the minimum is better than the previous one,
		// set the parameters as the best current ones
		if(newMinimum < minimum)
		{
			minimum = newMinimum;
			bestParameters = newParameters;
		}

		// copy the try parameters
		tryParameters = newParameters;

		delete pMinimizer;
 }

 std::cout << "Minimum reached at : " << minimum << std::endl;

 std::vector<double> hardcodedParameters;
 // from Imad TB
 hardcodedParameters.push_back(0.039081);
 hardcodedParameters.push_back(3.23877e-05);
 hardcodedParameters.push_back(-2.04628e-08); 
 hardcodedParameters.push_back(0.077336);
 hardcodedParameters.push_back(7.71267e-06);
 hardcodedParameters.push_back(-3.49567e-09);
 hardcodedParameters.push_back(0.119574);
 hardcodedParameters.push_back(2.40615e-05);
 hardcodedParameters.push_back(2.80403e-08);


 // from Imad FTFP_BERT_HP
 // hardcodedParameters.push_back( 0.0293659   );
 // hardcodedParameters.push_back( 2.96392e-05 );
 // hardcodedParameters.push_back(-2.051233e-08);
 // hardcodedParameters.push_back( 0.0925993   );
 // hardcodedParameters.push_back( 1.12284e-05 );
 // hardcodedParameters.push_back(-2.016464e-09);	
 // hardcodedParameters.push_back( 0.166332    );
 // hardcodedParameters.push_back( 3.08655e-05 );
 // hardcodedParameters.push_back( 2.97939e-08 );
 
 // uncomment if you want to use these parameters
 bestParameters = hardcodedParameters;

 for(unsigned int b=0 ; b<bestParameters.size() ; b++)
 	std::cout << "Best parameter " << b << " is : " << bestParameters.at(b) << std::endl;

 pMinimizerData->m_bestParameters = bestParameters;
}

