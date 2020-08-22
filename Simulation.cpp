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
    	//const int searches =300;


    	vector <double> distanze; //vector in cui salvare provvisoriamente le distanze raggiunte dal Random Walker durante un RW    
	int j, l; 
	long int i;
	
	
	gRandom->SetSeed(seed); 
    	TRandom* punta = gRandom;
    
	RandomWalker rwalker;
    	Data dati (rmin, rstep, rmax);
	const int walks=100000; //numero di random walk che si fanno per ottenere la distribuzione di probabilità degli n
	const long int maxsearches = 100000000000000;
 
	for (l=0; l<walks; l++){
		int nstep = 0;
		distanze.clear();
		rwalker.Getstart(rstart);
		double r = rstart;
		int m = 0;
		distanze.push_back(r);
		//if (debug) cout<< "il Random Walker è a distanza "<< r << endl;
		while (  r > rmin ) {  
			rwalker.RandomWalk(rstep, punta);
			r = rwalker.Distanza();
			if (r > rmax){
				nstep = 0;
				if (debug) cout << "il Random Walk " << l << "-esimo è terminato anzitempo perché il Random Walker ha superato rmax" << endl; 
				break;
				}
			distanze.push_back(r);
			//if (debug) cout<< "il Random Walker è a distanza "<< r << endl;
			nstep = nstep + 1;

		}
		
		double angle = rwalker.Angolo();
		if (nstep != 0) { 
			for(i = 0; i <=nstep; i++) {
				//if (debug) cout << distanze[i] <<","<<rmin <<"," <<rstep<<endl;
				dati.Add(distanze[i], nstep);
				dati.AddAngles(distanze[i], angle, nstep);
			}
		}
		
		
	}

	if (debugg) dati.ShowHistos();	
	if (debugg) dati.ShowHistosAngles();	

	//caso di un singolo operatore. 
	
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
	
	cout << "numero di ricerche con legame al sito " << j << endl;
	if (j < searches) cout << "maxsearches è troppo piccolo" << endl;
	double MeanOneOp = sum/j;
	cout << "numero medio di dissociazioni macroscopiche nel caso di un singolo operatore= " << MeanOneOp << endl; 

	double varianceOne = 0;
	
	for ( i = 0; i < j; i ++){
		varianceOne += pow((MOneOp[i] - MeanOneOp),2)/(j -1); 
	}
	varianceOne = TMath::Sqrt(varianceOne/j); 
	
	//caso di due operatori
	
	double DistTraOp[] = {10,20,30,40,50,60,80,100,120,140,160,180,200}; 
	
	int n = sizeof(DistTraOp)/sizeof(DistTraOp[0]);
	double errors[n];
	for (i=0; i<n; i++) errors[i]=0;

	double medie[n];
	double MTwoOp[searches];

	for (l = 0; l < n; l ++){
		
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
		cout << "numero di ricerche con legame al sito: " << j << endl;
		if (j < searches) cout << "maxsearches è troppo piccolo" << endl;
		double MeanTwoOp = sum2/j;
		printf("numero medio di dissociazioni macroscopiche nel caso di una distanza di %3.f bp tra i due operatori= ", DistTraOp[l]);
		cout << MeanTwoOp << endl; 
		medie[l] = MeanTwoOp;
		
		for ( i = 0; i < j; i ++){
			errors[l] += pow((MTwoOp[i] - MeanTwoOp),2)/(j -1);
			
		}
		errors[l] = TMath::Sqrt(errors[l]/j); 
	
		
	} 

	double rates[n];
	double errorirates[n];

	for (i = 0; i < n; i ++){
		rates[i] = MeanOneOp/medie[i];
		printf("i valori dei rate con rmax =  %2.f sono: ", rmax);
		cout << rates[i] << endl;
		errorirates[i] = TMath::Sqrt(pow((1/medie[i])* varianceOne,2)  + pow((MeanOneOp/pow(medie[i],2))*errors[i],2));
	}

   	TCanvas *c1 = new TCanvas("c1", "c1", 700, 500);
   	TGraphErrors *gr  = new TGraphErrors(n,DistTraOp,rates, 0, errorirates);
   	gr->SetTitle("");
   	gr->GetXaxis()->SetTitle("distance between operators [bp]");
   	gr->GetYaxis()->SetTitle("ratio of association rates (k/k0)");
   	gr->Draw("A*");
	
}


    
	


	




