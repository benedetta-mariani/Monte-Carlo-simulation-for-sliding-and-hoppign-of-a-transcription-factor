#include <vector>
#include <TCanvas.h>
#include <TH1D.h>
#include "TObject.h"
#ifndef DATA_H
#define DATA_H

class Data : public TObject{

	public:
  		Data(double rmin, double rstart, double rmax);
  		Data(const Data& source);
  		virtual ~Data();
  		void Add (double dist, int nstep);
  		void AddAngles(double dist, double angle, int nstep);
  		void Stampa(); 
		void ShowHistos();
		void ShowHistosAngles();
		int Casual(double r);
		int CasualAngles(double r);

	private:
		int fflag;
		int fflag2;
		bool fDebug;
		double frmax;
		double frmin;
		double frstart;
		double frstep;
		int fnrange;
		vector <int> *fcontainer;//[fnrange]
		vector <int> *fcontainer2;//[fnrange]
		TH1D **fhist;
		TH1D **fhist2;
		int fnamelength;
		int ftitlelength;
		char *fnome;//[fnamelength]
    	char *ftitolo;//[ftitlelength]
    	int fnamelength2;
		int ftitlelength2;
		char *fnome2;//[fnamelength]
    	char *ftitolo2;//[ftitlelength]	
    	void Histos();
    	void HistosAngles();

ClassDef(Data,1)
};

#endif
