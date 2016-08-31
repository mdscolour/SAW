/**************************************************************************
To start the program using simply in cmd:

$ ./start

This will write down a SAW in the file "FinalWalk".
Maximum 4 input parameters are possible:

$ ./start 1000 "0" 10 1000

1st length of the SAW
2nd the name of the initial walk, "0" means a new walk
3rd the number of pivot per step in MC chain, which will carry out at the same time at the end
4st number of step to go in MC chain
**************************************************************************/
#include "stellman.h"
#include "kennedy.h"


////////   ./start length 0or"FinalWalk"(name) discarded outerloop innerloop (in MCSs)
/////// when outerloop is 0, print autocorrelation time for 10 (MCSs) each record and 1000 records total.
int main(int argc,char *argv[])
{
/**************  random seed  *******************************/
SeedByTime();
/**************  main start here *******************************/
	//run_detect();
	CKennedy walk;
	int length=100;
	const char* init_name="0";
	int inner_loop = 1000;
	int outer_loop = 2;
	int discard = 100;

	if(argc>=2) length=atoi(argv[1]);
	if(argc>=3) init_name=argv[2];

	if(argc>=4) discard=atoi(argv[3]);
	if(argc>=5) outer_loop=atoi(argv[4]); // the forth is the outer loop and fifth is the inner
	if(argc>=6) inner_loop=atoi(argv[5]);

	walk.ImportSAW(init_name,length);
	for (int k=0;k<discard;k++)
	{
		SeedByTime();
		walk.Attempt_pivot(Random_symmetry(),Random_integer_uniform(0,length));
	}
	//walk.IncreaseGeneration();	
	for (int i=0;i<outer_loop;i++)
	{
		for (int j=0;j<inner_loop;j++)
		{
			SeedByTime();
			walk.Attempt_pivot(Random_symmetry(),Random_integer_uniform(0,length));
		}
		walk.IncreaseGeneration();
		//walk.Writedown("FinalWalk");
		walk.Record("data");
 	}

/**************  main end here *******************************/
#ifdef _WIN32
	printf("\n");
	system("pause");
#endif
	return 0;
}

void CStellman::line_initialize( int direction )
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

void CStellman::scan( FILE *fptr )
{
	int i;
	fscanf(fptr, "N=%d ind=%d ", &nsteps, &generation);
	printf("Start SAW is: N=%d ind=%d \n", nsteps, generation);

	steps = new Sphere[nsteps + 1];
	for (i = 0; i <= nsteps; i++) 
	{
		steps[i].scan(fptr);
	}
}

inline bool CStellman::Intersec_Stellman( Matrix& qt,int nt )
{
	Point pshift = steps[nt].center();
	for(int i=nt+1;i<=nsteps;i++)
	{
		Sphere pp = qt.dot(steps[i]-pshift)+pshift;
		int j = nt-1;

		while (j>=0)
		{
			double dis = steps[j].distance(pp);
			if (dis <= steps[i].r+steps[j].r)
			{
				//printf("try:%d, found %d and %d overlap\n",nt,i,j);
				return true;
			}
			else 
			{
				if(dis<=1) j--;
				else j -= (int)dis;
			}
		}
	}
	return false;
}

double CStellman::GetRg2()
{
	Point rc = GetStepi(0).center();
	double dis = 0;
	double rg2 = 0;

	for (int i=1;i<=nsteps;i++)
	{
		rc = rc + GetStepi(i).center();
	} 
	rc /= (double)(nsteps+1);

	for (int i=0;i<=nsteps;i++)
	{
		dis = (GetStepi(i).center()-rc).norm();
		rg2 += dis*dis;
	} 
	rg2 /= (nsteps+1);
	return(rg2);
}

void CStellman::deallocate()
{
	delete[] steps;
}

void CStellman::ImportSAW( const char* name,int n/*=100*/ )
{
	nsteps = n;
	if (name[0] == '0') line_initialize(2);
	//if (name == "0") line_initialize(2);
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

bool CStellman::Attempt_pivot( Matrix qt,int nt )
{
	bool intersec = Intersec_Stellman(qt,nt);
	if(!intersec)
	{
		Point pshift = steps[nt].center();
		for(int i=nt+1;i<=nsteps;i++)
		{
			Sphere ptemp = steps[i];
			steps[i] = qt.dot(ptemp-pshift)+pshift;
		}
		return true;
	}
	return false;
}

void CStellman::Writedown( const char* name )
{
	FILE *fptr;
	char buffer[50]; // <- danger, only storage for 256 characters.
	sprintf(buffer, "%s%d", name,generation);

	fptr = fopen(buffer, "w");
	fprintf(fptr, "N=%d ind=%d \n", nsteps, generation);
	for (int i = 0; i <= nsteps; i++) steps[i].print(fptr);
	fclose(fptr);
}

void CStellman::Record( const char* name )
{
	FILE *fptr;
	char buffer[50]; // <- danger, only storage for 256 characters.
	sprintf(buffer, "%s_%d", name,nsteps);
	// record the "data"
	double endnorm = GetStepi(nsteps).center().norm();
	fptr = fopen(buffer, "a");
	fprintf(fptr,"%14.10f    %14.10f \n",GetRg2(),endnorm*endnorm);
	fclose(fptr);
}

inline void CStellman::IncreaseGeneration()
{generation++;}

int CKennedy::find_segment( int itime, int npivot, int* ptime )
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

void CKennedy::clean_pivot()
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

void CKennedy::line_initialize( int direction )
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

void CKennedy::scan( FILE *fptr )
{
	int i;
	fscanf(fptr, "N=%d ind=%d ", &nsteps, &generation);
	printf("Start SAW is: N=%d ind=%d \n", nsteps, generation);

	steps = new Sphere[nsteps + 1];
	for (i = 0; i <= nsteps; i++) 
	{
		steps[i].scan(fptr);
	}
}

inline int CKennedy::pivot_strictly_saw( Matrix& qt,int nt )
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

void CKennedy::add_pivot( int pivot_loc, Matrix* poper, Point trans )
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

void CKennedy::simplify()
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

double CKennedy::GetRg2()
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

void CKennedy::deallocate()
{
	delete[] steps;
	delete[] ptime;
	delete[] igroup;
	delete[] shift;
}

void CKennedy::ImportSAW( const char* name,int n/*=100*/ )
{
	nsteps = n;
	//needsimplify = false;
	nsimplify = int(sqrt(double(nsteps / 40)));
	max_npivot = nsimplify+20;
	ptime = new int[max_npivot + 2];// this may be one larger than needed 
	igroup = new Matrix[max_npivot + 1];
	shift = new Point[max_npivot + 1];
	clean_pivot();
	if (name[0] == '0') line_initialize(2);
	//if (name == "0") line_initialize(2);
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

inline bool CKennedy::Attempt_pivot( Matrix qt,int nt )
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

Sphere CKennedy::GetStepi( int i )
{
	int iseg;
	iseg = find_segment(i, npivot, ptime);
	return igroup[iseg].dot(steps[i])+shift[iseg];
}

void CKennedy::IncreaseGeneration()
{generation++;}

void CKennedy::Writedown( const char* name )
{
	simplify();
	FILE *fptr;
	char buffer[50]; // <- danger, only storage for 256 characters.
	sprintf(buffer, "%s%d_%d", name,generation,nsteps);

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

void CKennedy::Record( const char* name )
{
	simplify();
	FILE *fptr;
	char buffer[50]; // <- danger, only storage for 256 characters.
	sprintf(buffer, "%s_%d", name,nsteps);
	// record the "data"
	double endnorm = steps[nsteps].center().norm();
	fptr = fopen(buffer, "a");
	fprintf(fptr,"%14.10f    %14.10f \n",GetRg2(),endnorm*endnorm);
	fclose(fptr);
}
