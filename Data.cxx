#include <Riostream.h>
#include "Data.h"
#include <vector>
#include <TCanvas.h>
#include <TH1D.h>
#include "TObject.h"
#include <TFile.h>

ClassImp(Data)

Data::Data(double rmin, double rstep, double rmax): TObject(),
	fflag(0),
	fflag2(0),
	fDebug(0.),
	fnamelength(10),
	ftitlelength(200),
	fnamelength2(10),
	ftitlelength2(200),
	frmin(rmin),
	frstep(rstep),
	frmax(rmax)	
{
	fnrange=int((frmax - frmin)/frstep) + 1;
	fcontainer = new vector <int> [fnrange];
	fnome = new char[fnamelength];
	ftitolo = new char[ftitlelength];
	fhist = new TH1D* [fnrange];
	fcontainer2 = new vector <int> [fnrange];
	fnome2 = new char[fnamelength2];
	ftitolo2 = new char[ftitlelength2];
	fhist2 = new TH1D* [fnrange];
}


Data::Data(const Data &source): TObject(source),
	fflag(0),
	fflag2(0),
	fDebug(source.fDebug),
	fnamelength(source.fnamelength),
	ftitlelength(source.ftitlelength),
  	frmin(source.frmin),
  	frstep(source.frstep),
	frmax(source.frmax),
	fnrange(source.fnrange),
	fnamelength2(source.fnamelength2),
	ftitlelength2(source.ftitlelength2)
{
	if(fnrange>0 && fnamelength>0 && ftitlelength>0){
		int i;
		fnome = new char[fnamelength];
		for(i=0;i<fnamelength;i++) fnome[i]=source.fnome[i];
		ftitolo = new char[ftitlelength];
		for(i=0;i<ftitlelength;i++) ftitolo[i]=source.ftitolo[i];
	  	fcontainer = new vector <int>[fnrange];
	  	fhist = new TH1D* [fnrange];
	  	fnome2 = new char[fnamelength2];
		for(i=0;i<fnamelength2;i++) fnome2[i]=source.fnome2[i];
		ftitolo2 = new char[ftitlelength2];
		for(i=0;i<ftitlelength2;i++) ftitolo2[i]=source.ftitolo2[i];
	  	fcontainer2 = new vector <int>[fnrange];
	  	fhist2 = new TH1D* [fnrange];

	  	for(i=0;i<fnrange;i++){
	  		fcontainer2[i]=source.fcontainer2[i];
	  		fhist2[i]=source.fhist2[i];
	  		
		}
  }
	
}


Data::~Data() {
	cout<<"DESTRUCTOR OF CLASS DATA - this = "<<this<<endl;
	if(fnrange>0) {
		for(int i=0; i<fnrange; i++){
			
			if (fcontainer[i].size()>0){    
				if(fflag){
				//if (fDebug) cout << "sto distruggendo l'istogramma numero "<<i<<endl; 
				delete fhist[i];
				}

			}
			if (fcontainer2[i].size()>0){   
				if(fflag2){
				//if (fDebug) cout << "sto distruggendo l'istogramma degli angoli numero "<<i<<endl; 
				delete fhist2[i];
				}
			}
		}
		
	  	delete []fcontainer;
	  	delete []fcontainer2;
		delete []fhist;
		delete []fhist2;
		
	}

	if(fnamelength>0)delete[]fnome;
	if(ftitlelength>0)delete[]ftitolo;
	if(fnamelength2>0)delete[]fnome2;
	if(ftitlelength2>0)delete[]ftitolo2;
}


void Data::Add (double dist, int nstep){
	if (nstep != 0 && dist < frmax) {
		int m;
		if ((dist-frmin) <0) m=0;
		else m = (int)((dist - frmin)/frstep);	
		fcontainer[m].push_back(nstep);
		
	}
	else if (nstep != 0 && dist > frmax) cout << "distanza troppo grande per essere raggiunta dal Random Walker"<< endl;
}


void Data::AddAngles (double dist, double angle, int nstep){
	if (nstep != 0 && dist < frmax) {
		int m;
		if ((dist-frmin) <0) m=0;
		else m = (int)((dist - frmin)/frstep);	
		fcontainer2[m].push_back(angle);
		
		
	}
	else if (nstep != 0 && dist > frmax) cout << "distanza troppo grande per essere raggiunta dal Random Walker"<< endl;
}


void Data::Stampa(){
	for (int i=0;i<fnrange;i++){
	cout<<"sto stampando il vettore numero"<<i<<endl;
		for (unsigned int j=0; j<fcontainer[i].size(); j++){
			cout << fcontainer[i][j]<< ","<< endl;
			
		}
		if (fcontainer[i].size()>0)cout <<"il massimo è: " << *max_element(fcontainer[i].begin(), fcontainer[i].end())<<endl;
	}
}

void Data::Histos(){
	fflag=1;
	for (int i=0; i<fnrange; i++) {
		sprintf(fnome,"h%d",i);
		sprintf(ftitolo,"numero di passi delle traiettorie che passano per le distanze del contenitore %d-esimo",i);
		if (fcontainer[i].size()>0){
			int x1 = *max_element(fcontainer[i].begin(), fcontainer[i].end());			
			fhist[i] = new TH1D(fnome,ftitolo, x1 + 1, 1, x1 + 2);	
			for(unsigned int j = 0; j < fcontainer[i].size(); j++) {
				fhist[i] -> Fill(fcontainer[i][j]);
			}
		}
	}		
}	


void Data::ShowHistos() {
	if(fflag){
		for (int i=0; i<fnrange; i++) {
			sprintf(fnome,"h%d",i);
			if (fcontainer[i].size()>0){
				new TCanvas(fnome,fnome);
				fhist[i] -> DrawCopy();
			}
		}
	}
	else{
		Data::Histos();
		for (int i=0; i<fnrange; i++) {
			sprintf(fnome,"h%d",i);
			if (fcontainer[i].size()>0){
				new TCanvas(fnome,fnome);
				fhist[i] -> DrawCopy();
			}
		}
	}
}


void Data::HistosAngles(){
	fflag2=1;
	for (int i=0; i<fnrange; i++) {
		sprintf(fnome2,"ha%d",i);
		sprintf(ftitolo2,"angoli di rebinding delle traiettorie che passano per le distanze del contenitore %d-esimo",i);
		if (fcontainer2[i].size()>0){
			int x1 = *max_element(fcontainer2[i].begin(), fcontainer2[i].end());
			int x2 = *min_element(fcontainer2[i].begin(), fcontainer2[i].end());			
			fhist2[i] = new TH1D(fnome2,ftitolo2, x1 - x2 + 2, x2 -1, x1 + 1);	
			for(unsigned int j = 0; j < fcontainer2[i].size(); j++) {
				fhist2[i] -> Fill(fcontainer2[i][j]);
			}
		}
	}		
}	


void Data::ShowHistosAngles() {
	if(fflag2){
		for (int i=0; i<fnrange; i++) {
			sprintf(fnome2,"ha%d",i);
			if (fcontainer2[i].size()>0){
				new TCanvas(fnome2,fnome2);
				fhist2[i] -> DrawCopy();
			}
		}
	}
	else{
		Data::HistosAngles();
		for (int i=0; i<fnrange; i++) {
			sprintf(fnome2,"ha%d",i);
			if (fcontainer2[i].size()>0){
				new TCanvas(fnome2,fnome2);
				fhist2[i] -> DrawCopy();
			}
		}
	}
}


int Data::Casual(double r){
	int m=0; 
	if (r >= frmin && r <frmax) m = (int)((r - frmin)/frstep);
	if (!fflag) Data::Histos();
	if (fcontainer[m].size()>0){
		return (int) fhist[m] -> GetRandom();
	} 
	else {
		cout << "valore di r che non ha un istogramma relativo al numero di passi" << endl; 
		return 0;
	}
	

}

int Data::CasualAngles(double r){
	int m=0; 
	if (r >= frmin && r <frmax) m = (int)((r - frmin)/frstep);
	if (!fflag2) Data::HistosAngles();
	if (fcontainer2[m].size()>0){
		return (int) fhist2[m] -> GetRandom();
	} 
	else {
		cout << "valore di r che non ha un istogramma relativo agli angoli" << endl;
		return 0;
	}
	

}