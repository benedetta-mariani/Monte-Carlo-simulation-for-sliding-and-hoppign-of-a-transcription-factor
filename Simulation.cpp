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
	
	const double rmin = 5.5; 
    const double rmax = 11;
    const double rstart = 6.5;
    const double rstep = rstart - rmin;
   

    vector <double> distanze; 
	int j, l; 
	long int i;
	
	
	gRandom->SetSeed(seed); 
    TRandom* punta = gRandom;
    
	RandomWalker rwalker;
    Data dati (rmin, rstep, rmax);
	const int walks=100000; 
	const long int maxsearches = 100000000000000;
 
	for (l=0; l<walks; l++){
		int nstep = 0;
		distanze.clear();
		rwalker.Getstart(rstart);
		double r = rstart;
		int m = 0;
		distanze.push_back(r);
		
		while (  r > rmin ) {  
			rwalker.RandomWalk(rstep, punta);
			r = rwalker.Distanza();
			if (r > rmax){
				nstep = 0;
				if (debug) cout << "The "<< l << "-th Random Walk has finished because the Random Walker has overcome rmax" << endl; 
				break;
				}
			distanze.push_back(r);
			nstep = nstep + 1;

		}
		
		double angle = rwalker.Angolo();
		if (nstep != 0) { 
			for(i = 0; i <=nstep; i++) {
				dati.Add(distanze[i], nstep);
				dati.AddAngles(distanze[i], angle, nstep);
			}
		}
		
		
	}

	if (debugg) dati.ShowHistos();	
	if (debugg) dati.ShowHistosAngles();	

	//One Operator
	
	TFactor TF (rmin, rmax, rstart, 20000, punta); 
	
	TF.ChangeDistanza(0);
	long double sum = 0;
	double MOneOp[searches];
	for (i = 0, j = 0; i < maxsearches && j < searches; i++){
		
		TF.SetPositionOperator(punta);
		
		int c = TF.Search(punta, dati); 


		if(c!=0) {
			sum += c;
			MOneOp[j] = c;
			j++;

		}

		
	}
	
	cout << "Number of searches with binding to the site: " << j << endl;
	if (j < searches) cout << "maxsearches is too small" << endl;
	double MeanOneOp = sum/j;
	cout << "Average value of macroscopic dissociations with one operator: " << MeanOneOp << endl; 

	double varianceOne = 0;
	
	for ( i = 0; i < j; i ++){
		varianceOne += pow((MOneOp[i] - MeanOneOp),2)/(j -1); 
	}
	varianceOne = TMath::Sqrt(varianceOne/j); 
	
	//Two Operator
	
	double DistTraOp[] = {0, 10,20,30,40,50,60,80,100,120,140,160,180,200}; 
	
	int n = sizeof(DistTraOp)/sizeof(DistTraOp[0]);
	double errors[n];
	for (i=0; i<n; i++) errors[i]=0;

	double medie[n];
	medie[0] = MeanOneOp;
	double MTwoOp[searches];

	for (l = 1; l < n; l ++){
		
		TF.ChangeDistanza(DistTraOp[l]);
		
		long double sum2 = 0;

		for (i=0, j=0; i < maxsearches && j < searches; i++){
					
			TF.SetPositionOperator(punta);			
			int c = TF.Search(punta, dati); 
			
			
			if(c!=0) {
				sum2 += c;
				MTwoOp[j] = c;
				j++;
			}
		} 

		cout << "Number of searches with binding to the site: " << j << endl;
		if (j < searches) cout << "maxsearches is too small" << endl;
		double MeanTwoOp = sum2/j;
		printf("Average value of macroscopic dissociations with a distance of %3.f bp between operators: ", DistTraOp[l]);
		cout << MeanTwoOp << endl; 
		medie[l] = MeanTwoOp;
		
		for ( i = 0; i < j; i ++){
			errors[l] += pow((MTwoOp[i] - MeanTwoOp),2)/(j -1);
			
		}
		errors[l] = TMath::Sqrt(errors[l]/j); 
	
		
	} 
	errors[0] = varianceOne; 

	double rates[n];
	double errorirates[n];

	for (i = 0; i < n; i ++){
		rates[i] = MeanOneOp/medie[i];
		printf("The ratio of association rates (k/k0) with rmax =  %2.f nm and distance between operators = %3.f bp is: ", rmax, DistTraOp[i]);
		cout << rates[i] << endl;
		errorirates[i] = TMath::Sqrt(pow((1/medie[i])* varianceOne,2)  + pow((MeanOneOp/pow(medie[i],2))*errors[i],2));
	}

	
   	TCanvas *c1 = new TCanvas("c1", "c1", 700, 500);
   	TGraphErrors *gr  = new TGraphErrors(n,DistTraOp,rates, 0, errorirates);
   	gr->SetTitle("");
   	gr->GetXaxis()->SetTitle("distance between operators [bp]");
   	gr->GetXaxis()->SetLimits(-5,DistTraOp[n-1]+10);    
   	gr->GetYaxis()->SetTitle("ratio of association rates (k/k0)");
   	gr->Draw("A*");
	
}


    
	


	




