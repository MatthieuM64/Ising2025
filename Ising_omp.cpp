/*C++ CODE - MANGEAT MATTHIEU - 2021*/
/*ACTIVE POTTS MODEL*/

//////////////////////
///// LIBRAIRIES /////
//////////////////////

//Public librairies.
#include <cstdlib>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <numeric>
#include <map>
#include <string>
#include <string.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <atomic>
#include <memory>
#include <omp.h>

using namespace std;

//Personal libraries.
#include "lib/random_OMP.cpp"
#include "lib/special_functions.cpp"

///////////////////////////
///// BASIC FUNCTIONS /////
///////////////////////////

//Index function
int index(const int &x0, const int &y0, const int &LX)
{
	return y0*LX+x0;
}

//Creation of the spin
void Ising_init(unique_ptr<atomic<int>[]> &SPIN, const int &Nsites, const int &init)
{
	for (int i=0; i<Nsites; i++)
	{
		if (init==0) //Disordered state.
		{
			SPIN[i].store(1-2*int(2*ran()),memory_order_relaxed);
		}
		else if (init==1) //Ordered state (positive mag).
		{
			SPIN[i].store(1,memory_order_relaxed);
		}
		else
		{
			cerr << "BAD INIT VALUE: " << init << endl;
			abort();
		}
	}
}

//Modulo function.
int modulo(const int &x, const int &LX)
{
	if (x<0)
	{
		return x+LX;
	}
	else if (x>=LX)
	{
		return x-LX;
	}
	else
	{
		return x;
	}
}

//Total average on all space.
double average(const unique_ptr<atomic<int>[]> &SPIN, const int &Nsites)
{
	int mag=0;
	for (int i=0; i<Nsites; i++)
	{
		mag+=SPIN[i].load(memory_order_relaxed);
	}
	return double(mag)/Nsites;
}

//Export state.
void exportState(const unique_ptr<atomic<int>[]> &SPIN, const double &beta, const double &h, const int &LX, const int &LY, const int &init, const int &RAN, const int &t)
{
	const int Nsites=LX*LY, nbytes=(Nsites+7)/8;
	vector<uint8_t> buf(nbytes,0);
	for (int i=0;i<Nsites;i++)
	{
		int s=SPIN[i].load(memory_order_relaxed);
		uint8_t bit = static_cast<uint8_t>((1 - s) >> 1);
		buf[i >> 3] |= static_cast<uint8_t>(bit << (7 - (i & 7)));
	}
	
	//Creation of the file.
	int returnSystem=system("mkdir -p data_Ising_dynamics2d/");
	stringstream ssState;
	ssState << "./data_Ising_dynamics2d/Ising_state_beta=" << beta << "_h=" << h << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << "_t=" << t << ".bin";
	string nameState = ssState.str();
	
	ofstream fileState(nameState.c_str(),ios::binary);
	fileState.write(reinterpret_cast<char*>(buf.data()), buf.size());
	fileState.close();
}

//Read parameters from command line.
void ReadCommandLine(int argc, char** argv, double &beta, double &h, int &LX, int &LY, int &tmax, int &init, int &RAN, int &THREAD_NUM)
{
 	for( int i = 1; i<argc; i++ )
	{
		if (strstr(argv[i], "-beta=" ))
		{
			beta=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-h=" ))
		{
			h=atof(argv[i]+3);
		}
		else if (strstr(argv[i], "-LX=" ))
		{
			LX=atoi(argv[i]+4);
		}
		else if (strstr(argv[i], "-LY=" ))
		{
			LY=atoi(argv[i]+4);
		}
		else if (strstr(argv[i], "-tmax=" ))
		{
			tmax=atoi(argv[i]+6);
		}
		else if (strstr(argv[i], "-init=" ))
		{
			init=atoi(argv[i]+6);
		}
		else if (strstr(argv[i], "-ran=" ))
		{
			RAN=atoi(argv[i]+5);
		}
		else if (strstr(argv[i], "-threads=" ))
		{
			THREAD_NUM=atoi(argv[i]+9);
		}
		else
		{
			cerr << "BAD ARGUMENT : " << argv[i] << endl;
			abort();
		}
	}
	cout << "-beta=" << beta << " -h=" << h << " -LX=" << LX << " -LY=" << LY << " -tmax=" << tmax << " -init=" << init << " -ran=" << RAN << " -threads=" << THREAD_NUM << endl;
}

/////////////////////
///// MAIN CODE /////
/////////////////////

int main(int argc, char *argv[])
{
	//Physical parameters: beta=inverse temperature, h=external magnetic field, LX*LY=size of the box.
	double beta=1., h=0.;
	int LX=512, LY=512;
	
	//Numerical parameters: init=initial condition, tmax=maximal time, RAN=index of RNG, THREAD_NUM=number of threads.
	int init=0, tmax=100000, RAN=0, THREAD_NUM=4;

	//Read imported parameters in command line.
	ReadCommandLine(argc,argv,beta,h,LX,LY,tmax,init,RAN,THREAD_NUM);

	//OpenMP.
	omp_set_dynamic(0);
	omp_set_num_threads(THREAD_NUM);
	cout << OMP_MAX_THREADS << " maximum threads on this node. " << THREAD_NUM << " threads will be used." << endl;

	//Start the random number generator.
	init_gsl_ran();
	for (int k=0; k<THREAD_NUM; k++)
	{
		gsl_rng_set(GSL_r[k],THREAD_NUM*RAN+k);
	}

	//Creation of the file to export global averages.
	const int dossier=system("mkdir -p ./data_Ising_averages/");
	stringstream ssAverages;
	ssAverages << "./data_Ising_averages/Ising_averages_beta=" << beta << "_h=" << h << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << ".txt";
	string nameAverages = ssAverages.str();
	ofstream fileAverages(nameAverages.c_str(),ios::trunc);
	fileAverages.precision(6);
	
	//Number of sites.
	int Nsites=LX*LY;
	
	//Ising spin in each site.
	unique_ptr<atomic<int>[]> SPIN(new atomic<int>[Nsites]);
	Ising_init(SPIN,Nsites,init);
	
	//Lock for the site.
	vector<omp_lock_t> lock_site(Nsites);
	for (int i=0; i<Nsites; i++)
	{
		omp_init_lock(&lock_site[i]);
	}
	
	//Time evolution.
	for(int t=0;t<=tmax;t++)
	{
		//Export data.
		const double mag=average(SPIN,Nsites);
		fileAverages <<  t << " " << mag << endl;
		
		if (t%25==0 or t==tmax)
		{
			exportState(SPIN,beta,h,LX,LY,init,RAN,t);
			cout << "time=" << t << " -mag=" << mag << running_time.TimeRun(" ") << endl;
		}
		
		//At each time-step update sites randomly.
		#pragma omp parallel for default(shared)
		for (int i=0; i<Nsites; i++)
		{
			//Choose a site randomly (j) unused by another thread (Hogwild-style asynchronous Metropolis algorithm).
			int j;
			do{
				j=int(Nsites*ran());
			}while(not omp_test_lock(&lock_site[j]));

			const int x0=j%LX, y0=j/LX, spin=SPIN[j].load(memory_order_relaxed);
			
			//Energy difference for a flip.
			const int xp=modulo(x0+1,LX), xm=modulo(x0-1,LX), yp=modulo(y0+1,LY), ym=modulo(y0-1,LY);
			const int MAG=SPIN[index(xp,y0,LX)].load(memory_order_relaxed) + SPIN[index(xm,y0,LX)].load(memory_order_relaxed) + SPIN[index(x0,yp,LX)].load(memory_order_relaxed) + SPIN[index(x0,ym,LX)].load(memory_order_relaxed);
			const double delH=2*spin*(MAG+h);
			
			//Metropolis algorithm.
			if (delH<0 or ran()<exp(-beta*delH))
			{
				SPIN[j].store(-spin,memory_order_relaxed);
			}
			omp_unset_lock(&lock_site[j]);
		}		
	}
	return 0;
}
