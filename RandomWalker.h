#ifndef RandomWalker_H
#define RandomWalker_H
#include "TRandom3.h"

class RandomWalker  {

public:
  
	RandomWalker();
	RandomWalker(double x, double y);  
	RandomWalker (const RandomWalker &source);
	virtual ~RandomWalker();
	void Getstart(double frstart);
	void RandomStep(double rstep, TRandom* a);
	void RandomWalk(double rstep, TRandom* a);
	double Angolo();
	double Distanza();
	void Stampaposizione(); 

private:
	double fx;
	double fy;
	double fxx;
	double fyy;


ClassDef(RandomWalker,1) 
};

#endif
