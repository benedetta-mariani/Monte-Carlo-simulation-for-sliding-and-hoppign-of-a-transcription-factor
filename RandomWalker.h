#ifndef RandomWalker_H
#define RandomWalker_H
#include "TRandom3.h"
#include "Data.h"

class RandomWalker  {

public:
  
	RandomWalker();
	RandomWalker(double x, double y, double rmin, double rstart, double rmax); 
	RandomWalker(const RandomWalker &source);
	virtual ~RandomWalker();
	void Getstart();
	void RandomStep(TRandom* a);
	void RandomWalk(TRandom* a);
	double Angolo();
	double Distanza();
	void Stampaposizione();
	void PerformWalks(Data& p, int walks, TRandom* punta);

private:
	double fx;
	double fy;
	double fxx;
	double fyy;
	double frmin;
	double frstart;
	double frmax;
	double frstep;


ClassDef(RandomWalker,1) 
};

#endif
