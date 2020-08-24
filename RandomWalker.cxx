#include <Riostream.h>
#include "RandomWalker.h"
#include "TMath.h"
#include "TRandom3.h"
#include "Data.h"

ClassImp(RandomWalker)

RandomWalker::RandomWalker():
fx(0.),
fy(0.),
fxx(0.),
fyy(0.),
frmin(5.5),
frstart(6.5), 
frmax(11)
{
  frstep = frstart - frmin;
}


RandomWalker::RandomWalker(double x, double y, double rmin, double rstart, double rmax): 
fx(x),
fy(y),
frmin(rmin),
frstart(rstart), 
frmax(rmax)
{
  fxx=fx;
  fyy=fy;
  frstep = frstart - frmin;
}

RandomWalker::RandomWalker(const RandomWalker &source): 
fx(source.fx),
fy(source.fy),
fxx(source.fxx),
fyy(source.fyy),
frstep(source.frstep),
frstart(source.frstart),
frmin(source.frmin),
frmax(source.frmax)
{
}


RandomWalker::~RandomWalker() {
	cout<<"DESTRUCTOR OF CLASS RANDOMWALKER - this = "<<this<<endl;
}


void RandomWalker::Getstart() {

	fx = frstart;
	fy = 0;
}


void RandomWalker::RandomStep(TRandom* a) {

	double t = a -> Rndm();
	t = 2*TMath::Pi()*t;
	fxx = fx;
	fyy = fy;
	fx += frstep*TMath::Cos(t);
	fy += frstep*TMath::Sin(t);
}


void RandomWalker::RandomWalk(TRandom* a) {
	int t = (int) 4*(a -> Rndm());
	switch(t){
	case 0: {
		fxx = fx;
		fx+= frstep;
	}
	break;
	
	case 1: {
		fyy = fy;
		fy+= frstep;
	}
	break;
	
	case 2:{
		fxx = fx;
		fx-= frstep;
	}
	break;
	
	case 3: {
		fyy = fy;
		fy-= frstep;
	} 
	break;
	}

}


double RandomWalker::Distanza() {
	double distanza = TMath::Sqrt(pow(fx,2) + pow(fy,2));
	return distanza;
}


double RandomWalker::Angolo() {
	
	double angolo1 = TMath::ATan2(fy,fx);
	double angolo2 = TMath::ATan2(fyy,fxx); // approximation
	double angolo = (angolo1 + angolo2)/2;
	angolo = angolo*180/TMath::Pi();
	return angolo; // angle in degrees! 
}



void RandomWalker::Stampaposizione(){
	cout << "The position of the Random Walker is " << fx <<","<< fy << " e " << fxx << "," << fyy << endl;
}

void RandomWalker::PerformWalks(Data& p, int walks, TRandom* punta){
	vector <double> distanze; 
	int l; 
	long int i;
	for (l=0; l<walks; l++){
		int nstep = 0;
		distanze.clear();
		Getstart();
		double r = frstart;
		int m = 0;
		distanze.push_back(r);
		
		while (  r > frmin ) {  
			RandomWalk(punta);
			r = Distanza();
			if (r > frmax){
				nstep = 0;
				//if (debug) cout << "The "<< l << "-th Random Walk has finished because the Random Walker has overcome rmax" << endl; 
				break;
				}
			distanze.push_back(r);
			nstep = nstep + 1;

		}
		
		double angle = Angolo();
		if (nstep != 0) { 
			for(i = 0; i <=nstep; i++) {
				p.Add(distanze[i], nstep);
				p.AddAngles(distanze[i], angle, nstep);
			}
		}
		
	}

}
