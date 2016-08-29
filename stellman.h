#ifndef STELLMANN
#define STELLMANN

#include "base.h"

class CStellman
{
private:
	Sphere* steps;
	int nsteps; // number of steps in the walk  
	int generation;
//	int npivot; // number of accepted pivots not yet applied to walk
//	int* ptime; // array of pivot times 
//	Matrix* igroup; // array of group elements  
//	Point* shift; // array of shifts

	void line_initialize(int direction);

	void scan(FILE *fptr);

	inline bool Intersec_Stellman(Matrix& qt,int nt);

	double GetRg2();

	void deallocate();
public:
	CStellman(){}

	// if name is "0", initialize SAW with step n, if not, read n from file
	void ImportSAW(const char* name,int n=100);

	bool Attempt_pivot(Matrix qt,int nt);

	Sphere GetStepi(int i){return(steps[i]);}
	int GetNSteps(){return nsteps;}
	inline void IncreaseGeneration();

	void Writedown(const char* name);

	void Record(const char* name);
};

#endif