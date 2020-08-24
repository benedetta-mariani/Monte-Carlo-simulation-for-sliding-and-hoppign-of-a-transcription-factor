#include <Riostream.h>
#include <TFormula.h>
#include <TRandom3.h>
#include <TRandom.h>
#include <TMath.h>
#include <vector>
#include <TF1.h>
#include <string>
#include <TFile.h>
#include <TCanvas.h>
#include "TFactor.h"
#include "Data.h"

using namespace std;

ClassImp(TFactor)


TFactor::TFactor(TRandom* poi): 
	fN(100),
	fs2(700), 
	fM(0), 
	fPosition(0),
	fp(1),
	frstart(6.5),
	frmin(5.5),
	frmax(30),
	fDebug(0),
	fTolleranza(36),
	fIterazioni(100000000000),
	fDistanza(0)
	
{
	fs2p = fs2/(1 + 2*fs2);
	frstep = frstart - frmin;
	int csi = poi -> Integer(fN);
	fPositionOperator = csi+1;

	if (fPositionOperator + fDistanza > fN){
		if (fPositionOperator - fDistanza > 0) fPositionOperator2 = fPositionOperator - fDistanza;
		else cout<<"The distance is too large"<<endl;
	}
	else if (fPositionOperator + fDistanza < fN && fPositionOperator - fDistanza < 0) fPositionOperator2 = fPositionOperator + fDistanza;
	
	else {
		double ran = poi -> Rndm();
		if (ran >= 0.5) fPositionOperator2 = fPositionOperator + fDistanza;
		else fPositionOperator2 = fPositionOperator - fDistanza;
	}
	if (fDebug) cout << "The positions of the operators are " << fPositionOperator << ", " << fPositionOperator2 << endl;

	
	char funzione[1000];
	sprintf(funzione, "1/TMath::Sqrt(TMath::Pi()*[0]*pow(%f,2))*TMath::Exp(-pow(x,2)/([1]*pow(%f,2)))",frstep,frstep);
			
	f1 = new TF1("myfunc",funzione,-400,400);
	jumps =  new TH1D("s","jumps",2000,-1000,-1000);	
}


TFactor::TFactor(double rmin, double rstart, double rmax, int N, TRandom* poi, double s2): 
	fN(N), 
	fs2(s2), 
	fM(0), 
	fPosition(0),
	fp(1), 
	frstart(rstart),
	frmin(rmin),
	frmax(rmax),
	fDebug(0),
	fIterazioni(100000000000),
	fTolleranza(36),
	fDistanza(0)

{
	fs2p = fs2/(1 + 2*fs2);
	frstep = frstart - frmin;
	int csi = poi -> Integer(fN);
	fPositionOperator = csi+1;

	if (fPositionOperator + fDistanza > fN){
		if (fPositionOperator - fDistanza > 0) fPositionOperator2 = fPositionOperator - fDistanza;
		else cout<<"The distance is too large"<<endl;
	}
	else if (fPositionOperator + fDistanza < fN && fPositionOperator - fDistanza < 0) fPositionOperator2 = fPositionOperator + fDistanza;
	
	else {
		double ran = poi -> Rndm();
		if (ran >= 0.5) fPositionOperator2 = fPositionOperator + fDistanza;
		else fPositionOperator2 = fPositionOperator - fDistanza;
	}
	if (fDebug) cout << "The positions of the operators are " << fPositionOperator << ", " << fPositionOperator2 << endl;

	char funzione[1000];
	sprintf(funzione, "1/TMath::Sqrt(TMath::Pi()*[0]*pow(%f,2))*TMath::Exp(-pow(x,2)/([1]*pow(%f,2)))",frstep,frstep);
			
	f1 = new TF1("myfunc",funzione,-400,400);
	jumps =  new TH1D("s","jumps",400,-200,200);	

}


TFactor::TFactor (const TFactor &source): 
	fN(source.fN), 
	fs2(source.fs2), 
	fs2p(source.fs2p),
	fM(source.fM), 
	fPosition(source.fPosition),
	fPositionOperator(source.fPositionOperator),
	fp(source.fp), 
	frstart(source.frstart),
	frmin(source.frmin),
	frmax(source.frmax),
	frstep(source.frstep),
	fDebug(source.fDebug),
	fIterazioni(source.fIterazioni),
	fTolleranza(source.fTolleranza),
	fDistanza(source.fDistanza),
	fPositionOperator2(source.fPositionOperator2),
	f1(source.f1),
	jumps(source.jumps)

{
}


TFactor::~TFactor(){
	delete f1; 
	delete jumps;
	cout<<"DESTRUCTOR OF CLASS TFACTOR - this = "<<this<<endl;
}


int TFactor::SetPosition(TRandom* poi) {
	
	int csi = poi -> Integer(fN);
	fPosition = csi+1;
	if (fPosition == fPositionOperator || fPosition == fPositionOperator2){
		if (fDebug) cout << "specific binding of the TF on its binding site"<< endl;
		if (fDebug) cout << "The value of M is " << fM << endl;
		return 2;

	}
	else {
		if (fDebug) cout << "The Position is " << fPosition << endl; 
		return 3;
	}
	
}


double TFactor::Selectr (TRandom* poi) { 
	double x = poi -> Rndm();
	double r = pow(frstart, 1/x)* pow(frmin, (1-x)/x);
	if (fDebug) cout <<"r vale " << r  << endl;
	return r;
	
}

	
int TFactor::Move(TRandom* poi) {

	if (fDebug) cout << "x is now extracted" << endl; 
	double x = poi -> Rndm();
	if( x <= fs2p){
		if (fDebug) cout << "Hopping has not happened, the TF slides" << endl;
		if (fPosition != fN) {fPosition = fPosition + 1;}
		if (fDebug) cout << "The current Position is " << fPosition << endl;
		
		if (fPosition == fPositionOperator || fPosition == fPositionOperator2){
			if (fDebug) cout << "specific binding of the TF on its binding site"<< endl;
			if (fDebug) cout << "The value of M is " << fM << endl;
			return 2;
			
		}
		else return 3;
		
	}
	
	else if (x > fs2p && x <= 2 * fs2p) {
		if (fDebug) cout << "Hopping has not happened, the TF slides" << endl;
		if (fPosition != 1) {fPosition = fPosition - 1;} 
		if (fDebug) cout << "The current Position is " << fPosition << endl;

		if (fPosition == fPositionOperator || fPosition == fPositionOperator2) {
			if (fDebug) cout << "specific binding of the TF on its binding site"<< endl;
			if (fDebug) cout << "The value of M is " << fM << endl;
			return 2;
			
		}
		else return 3;
		
	}
		
		
	else if( x > 2* fs2p && x <= 1) {

		if (fDebug) cout << "Now the TF tries to hop" << endl;
		return 4;
	}
	return 5;
}
	

int TFactor::Hopping (double r, int n , TRandom* poi, double angle) { 
	if (r >= frmax) { 
		fM = fM + 1;
		if (fDebug) cout << "Macroscopic dissociation" << endl;
		return 1;
	}
	else {
	
		if (n != 0) {
		    if (fDebug) cout << "Now the TF hops" << endl;
			f1->SetParameters(n,n);
			double z = f1 -> GetRandom();

			if (TMath::Abs(z) >= 23){ //23 nm = zmax
				fM = fM + 1;
				return 1;
			} 
	
			int jump = int(z /0.33); // conversion from nm to base pairs units
			double angolovecchio = fPosition*36%360; // 36 degrees: angle between two consecutive base pairs
			fPosition = fPosition + jump;
			if (fDebug) cout << "The position after hopping is " << fPosition << endl;
			if (fPosition > fN) fPosition = fN;
			else if (fPosition < 1) fPosition = 1;
			int angoloDna = fPosition*36%360; 
			int angolonuovo = (int)(angolovecchio + angle)%360;
			double x = poi -> Rndm();
			bool ang = kFALSE;
			if(angolonuovo <= angoloDna + fTolleranza/2 && angolonuovo >= angoloDna - fTolleranza/2) ang = kTRUE; //check if the angle follow the helix
			if (x <= fp && ang ) { 
				jumps-> Fill(jump);
				if (fDebug) cout <<"The TF binds with the DNA" << endl;
				if (fPosition == fPositionOperator || fPosition == fPositionOperator2) {
					if (fDebug) cout << "specific binding of the TF on its binding site after hopping"<< endl;
					if (fDebug) cout << "The value of M is " << fM << endl;
					return 2;
					
					}
				else return 3;
			}
			else {
				if (fDebug) cout <<"new hopping" << endl;
				return 4;
			}
		}

		else { 
			if (fDebug) cout << "This value of r is not valid" << endl;
			return 4;
		}
	}
}

int TFactor::Search(TRandom* poi, Data& p){
	fM = 0; 
	int b= SetPosition(poi);
	unsigned int i=0;
	double r;
	double c;

	while (b != 2 && b!=5 && i < fIterazioni){
		switch(b){

			case 1:{
				b = SetPosition(poi);
				if (b != 2)	b = Move(poi);
			}
			break;			

			case 3:{
				b = Move(poi);
			}
			break;

			case 4:{ 
				r = Selectr(poi);		
				c = (p.CasualAngles(r) + 360)%360;
				b = Hopping(r, p.Casual(r), poi, c);
			}
			break;

			case 5: cout << "No previous condition is valid " << endl;
			break;
		}
		i++;
	}
	if (b == 2) {
		return fM; 	


		}
	else {
		if (fDebug) cout << "The loop ended before the TFactor reached its binding site " << endl;
		return 0;
		} 

	

}	


double * TFactor::Searches(double distanza, int searches, long int maxsearches, TRandom* punta, Data& p){

	ChangeDistanza(distanza);
	int j, l; 
	long int i;
	long double sum = 0;
	double M[searches];
	for (i = 0, j = 0; i < maxsearches && j < searches; i++){
		
		SetPositionOperator(punta);
		
		int c = Search(punta, p); 


		if(c!=0) {
			sum += c;
			M[j] = c;
			j++;

		}

		
	}

	ReturnJumpHist(distanza); //uncomment these lines to show and then reset the histograms with hopping jumps length that lead to associations to DNA
	ResetJumpHist();
	
	cout << "Number of searches with binding to the site: " << j << endl;
	if (j < searches) cout << "maxsearches is too small" << endl;
	double Mean = sum/j;
	printf( "Average value of macroscopic dissociations with distance %3.f between operators: ", distanza);
	cout << Mean << endl; 

	double variance = 0;
	
	for ( i = 0; i < j; i ++){
		variance += pow((M[i] - Mean),2)/(j -1); 
	}
	variance = TMath::Sqrt(variance/j); 


	double *arr = new double[2]; 
	arr[0] =  Mean;
	arr[1] = variance;
	return arr;
	delete []arr;


}



void TFactor::ResetJumpHist(){

	jumps-> Reset();
}


void TFactor::ReturnJumpHist(double distance){

	char *tito = new char[1000];
	sprintf(tito,"Total hopping jump lengths for an operator distance of %3.f", distance);
	new TCanvas(tito,tito);
	jumps -> DrawCopy();
	delete []tito;
}

void TFactor::ChangeTolleranza(double tolleranza) {
	fTolleranza = tolleranza;
} 


void TFactor::ChangeDistanza(double distanza) {
	fDistanza = distanza;
	if (fDebug) cout << "The distance now is " << fDistanza << endl;
} 


void TFactor::SetPositionOperator(TRandom* poi){
	int csi = poi -> Integer(fN);
	fPositionOperator = csi+1;

	if (fPositionOperator + fDistanza > fN){
		if (fPositionOperator - fDistanza > 0) fPositionOperator2 = fPositionOperator - fDistanza;
		else cout<<"The distance is too large"<<endl;
	}
	else if (fPositionOperator + fDistanza < fN && fPositionOperator - fDistanza < 0) fPositionOperator2 = fPositionOperator + fDistanza;
	
	else {
		double ran = poi -> Rndm();
		if (ran >= 0.5) fPositionOperator2 = fPositionOperator + fDistanza;
		else fPositionOperator2 = fPositionOperator - fDistanza;
	}
	if (fDebug) cout << "The positions of the operators are " << fPositionOperator << ", " << fPositionOperator2 << endl;

}
