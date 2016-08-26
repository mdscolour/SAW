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
#include "detect.h"


////////   ./start length 0or"FinalWalk"(name) discarded outerloop innerloop (in MCSs)
/////// when outerloop is 0, print autocorrelation time for 10 (MCSs) each record and 1000 records total.
int main(int argc,char *argv[])
{
/**************  random seed  *******************************/
#ifdef linux
	struct timeval tpstart;
	gettimeofday(&tpstart, NULL);
	srand48(tpstart.tv_usec);
	//srand48(floor(mytime()));
#endif
#ifdef _WIN32
	srand(unsigned int(time(NULL)));
#endif
/**************  main start here *******************************/
	run_detect();


/**************  main end here *******************************/
#ifdef _WIN32
	printf("\n");
	system("pause");
#endif
	return 0;
}