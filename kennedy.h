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

	int find_segment(int itime, int npivot, int* ptime)
	{
		int iseg, isegl, isegu;

		if (itime >= ptime[npivot]) return(npivot);

		isegl = 0; isegu = npivot; iseg = 0;

		while (isegu > isegl + 1)
		{
			iseg = (isegl + isegu) / 2;
			if (itime < ptime[iseg]) isegu = iseg;
			else isegl = iseg;
		}
		return(isegl);
	} // find_segment()

	void clean_pivot()
	{
		npivot = 0;
		igroup[0].identity();
		shift[0].zero();
		// segment number iseg corresponds to times [ptime[iseg],ptime[iseg+1])
		// So when npivot=0 we should have : 
		ptime[0] = 0;
		ptime[1] = nsteps + 1;
		needsimplify = false;
	}

	void line_initialize(int direction)
	{ 
		steps = new Sphere[nsteps + 1];
		generation = 0;
		int i;
		double i1, i2, i3;
		i1 = 0; i2 = 0; i3 = 0;
		for (i = 0; i <= nsteps; i++)
		{
			steps[i].assign(i1, i2, i3, RADIUS, RIGID);
			//if(i == 5) steps[i].assign(i1, i2, i3, 0.2, 20);
			switch (direction) {
			case 1: i1++; break;
			case 2: if (i % 2 == 0) i2++; else i1++; break;
			case 3: i2++; break;
			default: printf("bad case in line_initialize() d\n"); exit(0); break;
			}
		}
	} // end line_initialize()
	void scan(FILE *fptr)
	{
		int i;
		fscanf_s(fptr, "N=%d ind=%d ", &nsteps, &generation);
		printf("Start SAW is: N=%d ind=%d \n", nsteps, generation);

		steps = new Sphere[nsteps + 1];
		for (i = 0; i <= nsteps; i++) 
		{
			steps[i].scan(fptr);
		}
	}

	int pivot_strictly_saw(Matrix& qt,int nt)
		// This is the version that changes i and j "simultaneously"
		// If the pivot is accepted, this routine carries it out. Otherwise it is
		// leaves the walk unchanged. 
		// Routine returns count if the new walk  is self-avoiding, 
		// returns -count if the pivot produces a self-intersection or 
		// intersects the excluded region, 
		// where count is the number of distance computations done.
		// This return value is only used to study how long the routine takes. 
	{
		Matrix *pgroup_jseg, *pgroup_iseg;
		int i, j, ip, jp, iseg, jseg, imin, imax, jmin, jmax;
		int separation, min_separation, sep_mod;
		Sphere  pp, stepsp, stepsi, stepsj;
		Point origin, transi, transj, shift_jseg, shift_iseg;
		int count, changei_flag;
		Matrix* poper = &qt;
		Matrix invoper = qt.inv();
		Matrix* pinvoper = &invoper;
		// pivot is given by w[t] -> g (w[t]-w[pivot_loc])+w[pivot_loc]
		// this is equivalent to w[t] -> g w[t] + trans 
		// where trans= w[pivot_loc] - g w[pivot_loc]
		// stepsp=w[pivot_loc]

		origin.zero();
		iseg = find_segment(nt, npivot, ptime);
		stepsp = igroup[iseg].dot(steps[nt])+shift[iseg];
		transi= poper->dot(Point(stepsp.x,stepsp.y,stepsp.z));
		transi.x = stepsp.x - transi.x;
		transi.y = stepsp.y - transi.y;
		transi.z = stepsp.z - transi.z;
		transj= pinvoper->dot(Point(stepsp.x,stepsp.y,stepsp.z));
		transj.x = stepsp.x - transj.x;
		transj.y = stepsp.y - transj.y;
		transj.z = stepsp.z - transj.z;
		count = 0;

		/////////////////////////////////////////////////////////////////////////////
		// 
		// "Simultaneous" changing of i and j. 
		// Given j<pivot_loc<i, we assume omega(jp)!=omega(ip) for j<jp<pivot_loc<ip<i
		// We then either decrease j or increase i.
		// 
		/////////////////////////////////////////////////////////////////////////////

		jmin = 0; jmax = nt - 1;
		imin = nt + 1; imax = nsteps;

		j = nt - 2;
		i = nt + 1;
		// if pivot_loc==0, j=-2 which would allow jp=-1
		if (j<jmin - 1) j = jmin - 1;

		if (nt>nsteps) i = imax + 1;

		if (true) while (i <= imax || j >= jmin)
		{
			// changei_flag=1 means we will increase i, =0 means we will increase j
			// We change the index that is closer to pivot_loc (if allowed).
			if (i - nt>nt - j) changei_flag = 0;
			else changei_flag = 1;
			if (i>imax) changei_flag = 0;
			if (j<jmin) changei_flag = 1;

			if (changei_flag)
			{// increase i. Need lower bound on distance from pivoted omega[i] to 
				// {omega[jp]: j<jp<pivot_loc}
				// This lower bound will be min_separation
				// stepsi is w[i]   
				iseg = find_segment(i, npivot, ptime);
				stepsi = igroup[iseg].dot(steps[i])+shift[iseg];
				// pp is w[i] after pivot
				pp=poper->dot(stepsi)+transi;
				min_separation = nsteps;
				jseg = npivot;   // can change to jseg=iseg; ?
				shift_jseg = shift[jseg];
				pgroup_jseg = &igroup[jseg];
				for (jp = jmax; jp>j;) // note that j is decreased
				{
					if (ptime[jseg]>jp)
					{
						//jseg = find_segment(jp, npivot, ptime);
						while (ptime[jseg]>jp) jseg--;
						shift_jseg = shift[jseg];
						pgroup_jseg = &igroup[jseg];
					}
					// stepsj is w[jp] - w[i]   
					stepsj=pgroup_jseg->dot(steps[jp])+shift_jseg;
					separation = stepsj.WellSeparate(pp);
					count++;
					if (separation == 0) return(-count);
					if (separation >= min_separation)
					{
						jp -= separation - min_separation + 1;
					}
					else
					{
						//sep_mod = separation % 3;
						//min_separation = (2 * separation + sep_mod) / 3;
						//jp -= 1 + (separation - sep_mod) / 3; 
						sep_mod = separation % 2;
						min_separation = (separation + sep_mod) / 2;
						jp -= 1 + (separation - sep_mod) / 2;
					}
				} // end loop on jp
				i += min_separation;
				if (i>imax + 1) i = imax + 1;
			} // end increase i 

			else
			{// decrease j. Need lower bound on distance from omega[j] to 
				// pivoted {omega[ip]: pivot_loc<ip<i}
				// Equivalently we can use a lower bound on distance from inverse
				// pivoted omega[j] and {omega[ip]: pivot_loc<ip<i}
				// This lower bound will be min_separation
				// stepsj is w[j] before this pivot
				jseg = find_segment(j, npivot, ptime);
				stepsj=igroup[jseg].dot(steps[j])+shift[jseg];
				pp = pinvoper->dot(stepsj)+transj;
				min_separation = nsteps;
				iseg = 0;
				shift_iseg = shift[iseg];
				pgroup_iseg = &igroup[iseg];
				for (ip = imin; ip<i;) // note that i is increased
				{
					if (ptime[iseg + 1] <= ip) // check this
					{
						//iseg = find_segment(ip, npivot, ptime);
						while (ptime[iseg + 1] <= ip) iseg++;
						shift_iseg = shift[iseg];
						pgroup_iseg = &igroup[iseg];
					}
					stepsi = pgroup_iseg->dot(steps[ip])+shift_iseg;
					//separation = int(stepsi.WellSeparate());
					separation = stepsi.WellSeparate(pp);
					count++;
					if (separation == 0) return(-count);
					if (separation >= min_separation)
					{
						ip += separation - min_separation + 1;
					}
					else
					{
						//sep_mod = separation % 3;
						//min_separation = (2 * separation + sep_mod) / 3;
						//ip += 1 + (separation - sep_mod) / 3;
						sep_mod = separation % 2;
						min_separation = (separation + sep_mod) / 2;
						ip += 1 + (separation - sep_mod) / 2;
					}
				} // end loop on ip
				j -= min_separation;
				if (j<jmin - 1) j = jmin - 1;
			} // end decrease j 
		} // end while 

		// If we reach this point the walk is self-avoiding. 
		add_pivot(nt, poper, transi);
		return(count);
	} // end pivot_strictly_saw()

	void add_pivot(int pivot_loc, Matrix* poper, Point trans)
	{
		int iseg, ipivot;
		Point pp;

		needsimplify=true;
		
		if (npivot>max_npivot - 1)
		{
			printf("number of implicit pivots in walk exceeds MAX_NPIVOT \n");
			exit(1);
		}

		// first, we add the pivot time to the list of pivot times. 
		iseg = find_segment(pivot_loc, npivot, ptime);
		if (pivot_loc != ptime[iseg])
		{
			ptime[npivot + 2] = ptime[npivot + 1];
			for (ipivot = npivot; ipivot >= iseg; ipivot--)
			{
				ptime[ipivot + 1] = ptime[ipivot];
				igroup[ipivot + 1] = igroup[ipivot];
				shift[ipivot + 1] = shift[ipivot];
			} // end loop on ipivot
			ptime[iseg + 1] = pivot_loc;
			npivot++;
			iseg++; // pivot will be applied to segments iseg to npivot
		}

		// second, we update igroup and shift
		for (ipivot = iseg; ipivot <= npivot; ipivot++)
		{
			pp = shift[ipivot];
			shift[ipivot] = poper->dot(pp)+trans;
			igroup[ipivot] = poper->dot(igroup[ipivot]);
		} // end loop on ipivot

	} // add_pivot()

	void simplify()
		// carry out the pivot operations implicit in the walk, so npivot -> 0 
	{	
		if(!needsimplify){return;}

		Matrix*  pOper_ipivot;
		int ipivot, itime;
		Sphere pp;
		Point shift_ipivot;

		// even on the 0th segment there may be something to do 
		for (ipivot = 0; ipivot<npivot; ipivot++)
		{
			shift_ipivot = shift[ipivot];
			pOper_ipivot = &igroup[ipivot];
			for (itime = ptime[ipivot]; itime<ptime[ipivot + 1]; itime++)
			{
				pp = steps[itime];
				steps[itime] = pOper_ipivot->dot(pp)+shift_ipivot;
			}
		} // end loop on ipivot

		// npivot segment is different
		for (itime = ptime[npivot]; itime <= nsteps; itime++)
		{
			pp = steps[itime];
			steps[itime] = igroup[npivot].dot(pp)+shift[npivot];
		}

		clean_pivot();
	} // end walk::simplify()

	//need simplified
	double GetRg2()
	{
		Point rc = steps[0].center();
		double dis = 0;
		double rg2 = 0;

		for (int i=1;i<=nsteps;i++)
		{
			rc = rc + steps[i].center();
		} 
		rc /= (double)(nsteps+1);

		for (int i=0;i<=nsteps;i++)
		{
			dis = (steps[i].center()-rc).norm();
			rg2 += dis*dis;
		} 
		rg2 /= (nsteps+1);
		return(rg2);
	}
	void deallocate()
	{
		delete[] steps;
		delete[] ptime;
		delete[] igroup;
		delete[] shift;
	};//weird
public:
	CKennedy(){}

	// if name is "0", initialize SAW with step n, if not, read n from file
	void ImportSAW(const char* name,int n=100)
	{
		nsteps = n;
		//needsimplify = false;
		nsimplify = int(sqrt(double(nsteps / 40)));
		max_npivot = nsimplify+20;
		ptime = new int[max_npivot + 2];// this may be one larger than needed 
		igroup = new Matrix[max_npivot + 1];
		shift = new Point[max_npivot + 1];
		clean_pivot();
		//if (name[0] == '0') line_initialize(2);
		if (name == "0") line_initialize(2);
		else
		{
			FILE *fptr = fopen(name, "r");
			if (fptr == NULL)
			{
				printf("FAILURE to open walk file %s, exiting \n", name);
				exit(0);
			}
			this->scan(fptr);
			fclose(fptr);
		}
	}

	bool Attempt_pivot(Matrix qt,int nt)
	{
		int accept_flag = pivot_strictly_saw(qt,nt);
		if (accept_flag>=0)
		{
			if (npivot >= nsimplify)simplify();
			return true;
		}
		//simplify();
		return false;
	}

	Sphere GetStepi(int i)
	{
		int iseg;
		iseg = find_segment(i, npivot, ptime);
		return igroup[iseg].dot(steps[i])+shift[iseg];
	}
	int GetNSteps(){return nsteps;}
	void IncreaseGeneration(){generation++;}

	void Writedown(const char* name)
	{
		simplify();
		FILE *fptr;
		char buffer[50]; // <- danger, only storage for 256 characters.
		sprintf_s(buffer, "%s%d", name,generation);

		fptr = fopen(buffer, "w");
		fprintf(fptr, "N=%d ind=%d \n", nsteps, generation);
		for (int i = 0; i <= nsteps; i++) steps[i].print(fptr);
// 		{
// 			int iseg;
// 			iseg = find_segment(i, npivot, ptime);
// 			(igroup[iseg].dot(steps[i])+shift[iseg]).print(fptr);
// 		}
		fclose(fptr);
	}
	void Record(const char* name)
	{
		simplify();
		FILE *fptr;
		char buffer[50]; // <- danger, only storage for 256 characters.
		sprintf_s(buffer, "%s_%d", name,nsteps);
		// record the "data"
		double endnorm = steps[nsteps].center().norm();
		fptr = fopen(buffer, "a");
		fprintf(fptr,"%14.10f    %14.10f \n",GetRg2(),endnorm*endnorm);
		fclose(fptr);
	}
};