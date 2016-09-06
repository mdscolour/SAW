// drand48 has 300 trillian period, not always good
#pragma once

#define CLISBY
#define KENNEDY
#define STELLMAN

#ifdef KENNEDY
#include "kennedy.h"
#endif
#ifdef STELLMAN
#include "stellman.h"
#endif
#ifdef CLISBY
#include "clisby.h"
#endif

bool thesameSAW(const char* saw1, const char* saw2)
{
	FILE *fptr1, *fptr2;
	int n1,n2,g1,g2;
	Sphere sp1,sp2;
	fptr1 = fopen(saw1, "r");
	fptr2 = fopen(saw2, "r");
	fscanf(fptr1, "N=%d ind=%d ", &n1, &g1);
	fscanf(fptr2, "N=%d ind=%d ", &n2, &g2);
	if (n1!=n2) {printf("not the same length\n");return false;}
	for (int i = 0; i <= n1; i++) 
	{
		sp1.scan(fptr1);
		sp2.scan(fptr2);
		if (sp1 == sp2) {}
		else
		{
			sp1.print(stdout);
			sp2.print(stdout);
			return false;
		}
	}
	fclose(fptr1);
	fclose(fptr2);
	return true;
}

// random number generator
double normalRandom(double mu = 0.0, double sigma = 1.0)
{
	double u1=RNG_NAME();
	double u2=RNG_NAME();
	return cos(8.*atan(1.)*u2)*sqrt(-2.*log(u1))*sigma+mu; 
}
void detect_rand(double true_mean=0, double true_stdev=1)
{
	std::cout << "Testing " << "normal" << " distribution\n";
	
	int sample_size = 10000;
	double sum = 0;
	for (int i = 0; i < sample_size; ++i)
		sum += normalRandom(true_mean,true_stdev);
	double mean = sum / sample_size;

	std::cout << "Computed mean: " << mean << "\n";
	double lower = true_mean - true_stdev/sqrt((double)sample_size);
	double upper = true_mean + true_stdev/sqrt((double)sample_size);
	std::cout << "Expected a value between " << lower << " and " << upper << "\n\n";
}

// point 
void detect_base()
{
	Point a(1,0,0),b(1,2,0);
	int arr[10000];
	for (int i=0;i<10000;i++)
	{
		Matrix m = Random_symmetry();
		arr[i]=(int)m.dot(a).coord_x()*100+100;
	}
	for (int j=0;j<=205;j++)
	{
		int i;
		for (i=0;i<10000;i++)
		{
			if(arr[i]==j)break;
		}
		if(i==10000)cout<<j<<endl;
	}
	
}
//stellman
#ifdef STELLMAN
void detect_stellman()
{
	CStellman walk;
	int n = 10;
	walk.ImportSAW("0",n);
	
// 	for(int i=0;i<walk.GetNSteps();i++)
// 	{
// 		printf("%.4lf \n",walk.GetStepi(i).distance(walk.GetStepi(i+1)));
// 	}

	Matrix qtarr[4];
	qtarr[0] = Matrix(1,0,0,0,-1,0,0,0,1);
	qtarr[1] = Matrix(-1,0,0,0,1,0,0,0,1);
	qtarr[2] = Matrix(1,0,0,0,-1,0,0,0,1);
	double theta = M_PI/4;
	qtarr[3] = Matrix(cos(theta),-sin(theta),0,sin(theta),cos(theta),0,0,0,1);

	int ntarr[] = {8,9,0,9};
	//walk.Writedown("FinalWalk");
	//walk.Record("data");
	for (int i=0;i<4;i++)
	{
		walk.Attempt_pivot(qtarr[i],ntarr[i]);
		walk.IncreaseGeneration();
		//walk.Writedown("FinalWalk");
		//walk.Record("data");
	}

// 	Matrix qt=Random_symmetry();
// 	int nt = 0;
// 	walk.Attempt_pivot(qt,nt);
// 	walk.Record("data");
}
#endif

//kennedy
#ifdef KENNEDY
void detect_kennedy()
{
	CKennedy walk;
	int n = 10;
	walk.ImportSAW("0",n);

	// 	for(int i=0;i<walk.GetNSteps();i++)
	// 	{
	// 		printf("%.4lf \n",walk.GetStepi(i).distance(walk.GetStepi(i+1)));
	// 	}

	Matrix qtarr[4];
	qtarr[0] = Matrix(1,0,0,0,-1,0,0,0,1);
	qtarr[1] = Matrix(-1,0,0,0,1,0,0,0,1);
	qtarr[2] = Matrix(1,0,0,0,-1,0,0,0,1);
	double theta = M_PI/4;
	qtarr[3] = Matrix(cos(theta),-sin(theta),0,sin(theta),cos(theta),0,0,0,1);

	int ntarr[] = {8,9,0,9};
	walk.Writedown("FinalWalk");
	walk.Record("data");
	for (int i=0;i<4;i++)
	{
		walk.Attempt_pivot(qtarr[i],ntarr[i]);
		walk.IncreaseGeneration();
		walk.Writedown("FinalWalk");
		walk.Record("data");
	}
}
#endif

#ifdef CLISBY
void detect_clisby()
{
	CClisby walk;
	int n = 10;
	walk.ImportSAW("0",n);

	// 	for(int i=0;i<walk.GetNSteps();i++)
	// 	{
	// 		printf("%.4lf \n",walk.GetStepi(i).distance(walk.GetStepi(i+1)));
	// 	}

	Matrix qtarr[4];
	qtarr[0] = Matrix(1,0,0,0,-1,0,0,0,1);
	qtarr[1] = Matrix(-1,0,0,0,1,0,0,0,1);
	qtarr[2] = Matrix(1,0,0,0,-1,0,0,0,1);
	double theta = M_PI/4;
	qtarr[3] = Matrix(cos(theta),-sin(theta),0,sin(theta),cos(theta),0,0,0,1);

	int ntarr[] = {8,9,0,9};
	walk.Writedown("FinalWalk");
	walk.Record("data");
	for (int i=0;i<4;i++)
	{
		walk.Attempt_pivot(qtarr[i],ntarr[i]);
		walk.IncreaseGeneration();
		walk.Writedown("FinalWalk");
		walk.Record("data");
	}
}
#endif

#ifdef KENNEDY
#ifdef STELLMAN
void detect_kenandstell()
{
	int n = 10000;
	Matrix qt;
	int nt;

	CStellman swalk;
	swalk.ImportSAW("0",n);
	CKennedy kwalk;
	kwalk.ImportSAW("0",n);

	for(int i=1;i<100000;i++)
	{
		qt = Random_symmetry();
		nt = Random_integer_uniform(0,n);
		swalk.Attempt_pivot(qt,nt);
		swalk.IncreaseGeneration();
		
		kwalk.Attempt_pivot(qt,nt);
		kwalk.IncreaseGeneration();
		if (i%1000 == 0)
		{
			swalk.Writedown("temp/stell");
			kwalk.Writedown("temp/ken");

			char buffer[50]; // <- danger, only storage for 256 characters.
			sprintf(buffer, "%s%d", "temp/stell",i);
			char buffer2[50]; // <- danger, only storage for 256 characters.
			sprintf(buffer2, "%s%d", "temp/ken",i);
			if (thesameSAW(buffer,buffer2)){printf("%d ",i);}
			else printf("\n%d not the same +++++++++++++++++\n",i);
		}
	}
}
#endif
#endif

void run_detect()
{
	//detect_rand();
	//detect_stellman();
	//detect_kennedy();
	//detect_clisby();  //not pass the below
	//detect_kenandstell();

	//thesameSAW("temp/stell59","temp/ken59");

	clock_t start, finish; 
	int n = 1024-1;

	start = clock();  
	CClisby walk;
	walk.ImportSAW("0",n);
	for (int i=0;i<1000;i++)
	{
		walk.Attempt_pivot(Random_symmetry(),Random_integer_uniform(0,n));
		walk.IncreaseGeneration();
	}
	for(int i=0;i<walk.GetNSteps();i++)
	{
		double t = walk.GetStepi(i).distance(walk.GetStepi(i+1));
		if(fabs(t-1)>1e-5) 
			printf("%.4lf \n",t);
	}
	finish = clock(); 
	printf( "%f seconds\n", (double)(finish - start) / CLOCKS_PER_SEC ); 
}