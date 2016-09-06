#pragma once

#include "base.h"

class CKennedy
{
private:
	Sphere* steps;
	int nsteps; // number of steps in the walk  
	int generation;
	int max_npivot;
	//double old_energy;
	int npivot; // number of accepted pivots not yet applied to walk
	int* ptime; // array of pivot times 
	Matrix* igroup; // array of group elements  
	Point* shift; // array of shifts
	int nsimplify;
	bool needsimplify;

	int find_segment(int itime, int npivot, int* ptime);

	void clean_pivot();

	void line_initialize(int direction);

	void scan(FILE *fptr);

	inline int pivot_strictly_saw(Matrix& qt,int nt);
		// This is the version that changes i and j "simultaneously"
		// If the pivot is accepted, this routine carries it out. Otherwise it is
		// leaves the walk unchanged. 
		// Routine returns count if the new walk  is self-avoiding, 
		// returns -count if the pivot produces a self-intersection or 
		// intersects the excluded region, 
		// where count is the number of distance computations done.
		// This return value is only used to study how long the routine takes. 

	void add_pivot(int pivot_loc, Matrix* poper, Point trans);

	void simplify();
		// carry out the pivot operations implicit in the walk, so npivot -> 0 

	//need simplified
	double GetRg2();

	void deallocate();
public:
	CKennedy(){}

	// if name is "0", initialize SAW with step n, if not, read n from file
	void ImportSAW(const char* name,int n=100);

	inline bool Attempt_pivot(Matrix qt,int nt);

	Sphere GetStepi(int i);

	int GetNSteps(){return nsteps;}
	inline void IncreaseGeneration();

	void Writedown(const char* name);

	void Record(const char* name);
};