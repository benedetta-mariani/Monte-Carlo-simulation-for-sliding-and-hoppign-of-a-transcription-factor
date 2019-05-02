#include <Riostream.h>
#include "RandomWalker.h"
#include "TMath.h"
#include "TRandom3.h"

ClassImp(RandomWalker)

RandomWalker::RandomWalker(): 
fx(0.),
fy(0.),
fxx(0.),
fyy(0.)
{
}


RandomWalker::RandomWalker(double x, double y): 
fx(x),
fy(y)
{
  fxx=fx;
  fyy=fy;
}

RandomWalker::RandomWalker(const RandomWalker &source): 
fx(source.fx),
fy(source.fy),
fxx(source.fxx),
fyy(source.fyy)
{
}


RandomWalker::~RandomWalker() {
	cout<<"DESTRUCTOR OF CLASS RANDOMWALKER - this = "<<this<<endl;
}


void RandomWalker::Getstart(double frstart) {

	fx = frstart;
	fy = 0;
}


void RandomWalker::RandomStep(double rstep, TRandom* a) {

	double t = a -> Rndm();
	t = 2*TMath::Pi()*t;
	fxx = fx;
	fyy = fy;
	fx += rstep*TMath::Cos(t);
	fy += rstep*TMath::Sin(t);
}


void RandomWalker::RandomWalk(double rstep, TRandom* a) {
	int t = (int) 4*(a -> Rndm());
	switch(t){
	case 0: {
		fxx = fx;
		fx+= rstep;
	}
	break;
	
	case 1: {
		fyy = fy;
		fy+= rstep;
	}
	break;
	
	case 2:{
		fxx = fx;
		fx-= rstep;
	}
	break;
	
	case 3: {
		fyy = fy;
		fy-= rstep;
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
	double angolo2 = TMath::ATan2(fyy,fxx);
	double angolo = (angolo1 + angolo2)/2;
	angolo = angolo*180/TMath::Pi();
	return angolo;
}



void RandomWalker::Stampaposizione(){
	cout << "la posizione del RandomWalker e' " << fx <<","<< fy << " e " << fxx << "," << fyy << endl;
}

