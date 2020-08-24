#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TH1D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TFormula.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TF1.h>
#include <string>
#include <TObject.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStopwatch.h>
#include "RandomWalker.h"
#include "Data.h"
#include "TFactor.h"
#endif

using namespace std;


void simulation(int searches, unsigned int seed = 123) {

	Bool_t debug = kFALSE;
	Bool_t debugg= kFALSE;

	
	const double s2 = 700;
	const double rmin = 5.5; 
    const double rstart = 6.5;
    const double rmax = 20;
    
	int l; 

	
	
	gRandom->SetSeed(seed); 

    TRandom* punta = gRandom;

    ///////////////////Instances////////////////////////////

	RandomWalker rwalker(0, 0, rmin, rstart, rmax);
    Data dati (rmin, rstart, rmax);
    TFactor TF (rmin, rstart, rmax, 20000, punta, s2); 
 	
 	///////////////////////////////////////////////////////

	const int walks=100000; 
	const long int maxsearches = 100000000000000;


	rwalker.PerformWalks(dati,walks,punta);

	if (debugg) dati.ShowHistos();	
	if (debugg) dati.ShowHistosAngles();	

	double DistTraOp[] = {0,10,20,30,40,50,60,80,100,120,140,160,180,200}; 
	
	int n = sizeof(DistTraOp)/sizeof(DistTraOp[0]);
	double errors[n];
	double medie[n];


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//One Operator
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	double *poi;
	poi = TF.Searches(0, searches, maxsearches, punta, dati);
	medie[0]  = *(poi);
	errors[0] = *(poi+1);


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Two Operators
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	for (l = 1; l < n; l ++){
		
		poi = TF.Searches(DistTraOp[l], searches, maxsearches, punta, dati);
		medie[l] = *(poi);
		errors[l] = *(poi+1);
		
	} 

	
	double rates[n];
	double errorirates[n];

	for (l = 0; l < n; l ++){
		rates[l] = medie[0]/medie[l];
		printf("The ratio of association rates (k/k0) with rmax =  %2.f nm and distance between operators = %3.f bp is: ", rmax, DistTraOp[l]);
		cout << rates[l] << endl;
		errorirates[l] = TMath::Sqrt(pow((1/medie[l])* errors[0],2)  + pow((medie[0]/pow(medie[l],2))*errors[l],2));
	}


	
   	TCanvas *c1 = new TCanvas("c1", "c1", 700, 500);
   	TGraphErrors *gr  = new TGraphErrors(n,DistTraOp,rates, 0, errorirates);
   	gr->SetTitle("");
   	gr->GetXaxis()->SetTitle("distance between operators [bp]");
   	gr->GetXaxis()->SetLimits(-5,DistTraOp[n-1]+10);    
   	gr->GetYaxis()->SetTitle("ratio of association rates (k/k0)");
   	gr->Draw("A*");
   
	
}


    
	


	




