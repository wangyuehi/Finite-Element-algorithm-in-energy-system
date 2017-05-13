
#include "stdio.h"
#include <stdlib.h>
#include "femtestGPv1.h"
#include <windows.h>  
#include <time.h> //time_t time()  clock_t clock()  
#include <Mmsystem.h>             //timeGetTime()  
#pragma comment(lib, "Winmm.lib")   //timeGetTime()  

int main(int argc, int* argv[])
{
	int retval = 0;


	FILE *fp;
	int i, j;
	double n[1000] = { 0 };
	int m[1008][3] = { 0 };
	double x[600] = { 0 };
	double y[600] = { 0 };
	time_t timeBegin, timeEnd;
	DWORD  dwBegin, dwEnd;
	double F_thrust1test = 0;

	readFile(m, x, y);
	timeBegin = clock();



	TOP(n, m, x, y);


	
	int a = 0;
	scanf("%d", a);

	return 0;
}