// drand48 has 300 trillian period, not always good

#define KENNEDY

#ifdef KENNEDY
#include "kennedy.h"
#endif
#ifdef STELLMAN
#include "stellman.h"
#endif
#ifdef CLISBY
#include "clisby.h"
#endif

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
}
#endif

void run_detect()
{
	//detect_rand();
	//detect_stellman();
	detect_kennedy();
	//detect_clisby();
	//check several step...to be done
}