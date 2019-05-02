#ifndef TFactor_H
#define TFactor_H
#include "Data.h"
#include "TRandom3.h"
#include "string"
#include "TF1.h"

class TFactor {


public:

	TFactor(TRandom* poi);
	TFactor(double rmin, double rmax, double frstart, int N, TRandom* poi);
	TFactor (const TFactor &source);
	virtual ~TFactor();	
	int Hopping (double r, int n, TRandom* poi, double angle);
	int SetPosition(TRandom* poi);
	double Selectr(TRandom* poi);
	int Move(TRandom* poi);
	void ChangeTolleranza(double tolleranza);
	int Search(TRandom* poi, Data& p);
	void ChangeDistanza(double distanza);
	void SetPositionOperator(TRandom* poi);


private:
	int fN;
	double frstart;
	double frmin;
	double frmax;
	double frstep;	
	double fs2;
	double fs2p;
	int fM;
	int fPosition;
	int fPositionOperator;
	int fPositionOperator2;
	double fp;
	bool fDebug;
	unsigned long int fIterazioni; 
	double fTolleranza;
	double fDistanza;
	TF1 *f1;
	
	
 ClassDef(TFactor,1) 
};

#endif
