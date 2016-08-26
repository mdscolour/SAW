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

	bool Intersec_Stellman(Matrix& qt,int nt)
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
	double GetRg2()
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
	void deallocate()
	{
		delete[] steps;
	}
public:
	CStellman(){}

	// if name is "0", initialize SAW with step n, if not, read n from file
	void ImportSAW(const char* name,int n=100)
	{
		nsteps = n;
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

	Sphere GetStepi(int i){return(steps[i]);}
	int GetNSteps(){return nsteps;}
	void IncreaseGeneration(){generation++;}

	void Writedown(const char* name)
	{
		FILE *fptr;
		char buffer[50]; // <- danger, only storage for 256 characters.
		sprintf_s(buffer, "%s%d", name,generation);

		fptr = fopen(buffer, "w");
		fprintf(fptr, "N=%d ind=%d \n", nsteps, generation);
		for (int i = 0; i <= nsteps; i++) steps[i].print(fptr);
		fclose(fptr);
	}
	void Record(const char* name)
	{
		FILE *fptr;
		char buffer[50]; // <- danger, only storage for 256 characters.
		sprintf_s(buffer, "%s_%d", name,nsteps);
		// record the "data"
		double endnorm = GetStepi(nsteps).center().norm();
		fptr = fopen(buffer, "a");
		fprintf(fptr,"%14.10f    %14.10f \n",GetRg2(),endnorm*endnorm);
		fclose(fptr);
	}
};