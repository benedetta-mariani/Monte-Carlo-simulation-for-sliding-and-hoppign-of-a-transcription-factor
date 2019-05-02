#include <Riostream.h>
#include <TFormula.h>
#include <TRandom3.h>
#include <TRandom.h>
#include <TMath.h>
#include <TF1.h>
#include <string>
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
		else cout<<"Distanza troppo grande"<<endl;
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
}


TFactor::TFactor(double rmin, double rmax, double rstart, int N, TRandom* poi): 
	fN(N), 
	fs2(700), 
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
		else cout<<"Distanza troppo grande"<<endl;
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
	f1(source.f1)
{
}


TFactor::~TFactor(){
	delete f1; 
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
		double r = Selectr(poi);
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
			if (TMath::Abs(z) >= 23){
				fM = fM + 1;
				return 1;
			}
			int jump = int(z /0.33); // conversione da z a numero di basi azotate
			double angolovecchio = fPosition*36%360;
			fPosition = fPosition + jump;
			if (fDebug) cout << "The position after hopping is " << fPosition << endl;
			if (fPosition > fN) fPosition = fN;
			else if (fPosition < 1) fPosition = 1;
			int angoloDna = fPosition*36%360; 
			int angolonuovo = (int)(angolovecchio + angle)%360;
			double x = poi -> Rndm();
			bool ang = kFALSE;
			if(angolonuovo <= angoloDna + fTolleranza/2 && angolonuovo >= angoloDna - fTolleranza/2) ang = kTRUE;
			if (x <= fp && ang ) { 
				if (fDebug) cout <<"The TF binds with the DNA" << endl;
				if (fPosition == fPositionOperator || fPosition == fPositionOperator2) {
					if (fDebug) cout << "specific binding of the TF on its binding site"<< endl;
					if (fDebug) cout << "The value of M is " << fM << endl;
					return 2;
					
					}
				else return 3;
			}
			else {
				if (fDebug) cout <<"new hopping" << endl;
				double r = Selectr(poi);
				return 4;
			}
		}

		else { 
			if (fDebug) cout << "This value of r is not valid" << endl;
			double r = Selectr(poi);
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

			case 5: cout << "Nessuna delle condizioni poste è valida " << endl;
			break;
		}
		i++;
	}
	if (b == 2)	return fM; 		
	else {
		if (fDebug) cout << "Il loop è terminato senza che il TFactor abbia raggiunto il sito di legame. " << endl;
		return 0;
		} 
}	


void TFactor::ChangeTolleranza(double tolleranza) {
	fTolleranza = tolleranza;
} 


void TFactor::ChangeDistanza(double distanza) {
	fDistanza = distanza;
	if (fDebug) cout << "la distanza ora vale " << fDistanza << endl;
} 


void TFactor::SetPositionOperator(TRandom* poi){
	int csi = poi -> Integer(fN);
	fPositionOperator = csi+1;

	if (fPositionOperator + fDistanza > fN){
		if (fPositionOperator - fDistanza > 0) fPositionOperator2 = fPositionOperator - fDistanza;
		else cout<<"Distanza troppo grande"<<endl;
	}
	else if (fPositionOperator + fDistanza < fN && fPositionOperator - fDistanza < 0) fPositionOperator2 = fPositionOperator + fDistanza;
	
	else {
		double ran = poi -> Rndm();
		if (ran >= 0.5) fPositionOperator2 = fPositionOperator + fDistanza;
		else fPositionOperator2 = fPositionOperator - fDistanza;
	}
	if (fDebug) cout << "The positions of the operators are " << fPositionOperator << ", " << fPositionOperator2 << endl;

}
