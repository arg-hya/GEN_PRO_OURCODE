#pragma once

#include <iostream>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <math.h>       /* pow , sqrt  */

class GenAlgo
{
private:

	/*=================
	TYPE DEFINTIONS :
	=================*/
	struct indiv
	{
		char *fChromo;				//Chromosome String of each individual.
		double fFitness;				//Fittness of each individual.
	};
	typedef struct indiv INDIVIDUAL;
	typedef INDIVIDUAL *POPULATION;        /* array of individuals */
	typedef INDIVIDUAL CHILDREN;        /* array of individuals */
	typedef INDIVIDUAL *TABLE;			/* tables used */

	typedef  double Fitness;   

	/*=================
	GLOBAL VARIABLES :
	=================*/
	unsigned int MaxFitness_Calculations;	//Stores the maximum number of fittness calculations permitted.

	int		nVAR,						//Number of variables in the obj func.
			nGEN,						//Number of Generations to be run.
			nPOPU,						//Population size.
		   *lenChromo_var,				//Chromosome length of each variable.
			lenChromo_tot,				//Chromosome length of each individual :: (nVAR * lenChromo_var).
		    func_num,					//CEC'13 Benchmark function number being evaluated
			func_num1,
			T,							//Length of Pb.
			D;							//Length of Pe.
			//indx_f,						//Index of one parent.
			//indx_m;					    //Index of another parent.


	float	 *LW_BND,				//Lower Bound and Higher Bound of each variable.
			 *HG_BND,				//Lower Bound and Higher Bound of each variable.				 
			  Pc ,					//Cross-Over Propability.
		      Pm ,					//Mutation Propability.
			  Ptb ,					//Crossover from Tables Propability.
			  BASIC_SEED;			//basic seed for rand().
			

	Fitness *best_fitness_array,	//Stores the best fittness level of each generation.
			 mean_best_fitness,		//Stores the mean value of the fitness levels reached in the entire run
			 std_best_fitness,		//Stores the standard deviation.
			 avgRunningFITNESS,		//average fittness of the current population.
		    *Xi;					//Variable importance factor table.			


	POPULATION	POPU = NULL;
	POPULATION	Temp_popu = NULL;
	CHILDREN child3;

	TABLE	Pb = NULL,		//Best Candidate table
			Pxb = NULL,		//Propable Candidate table
			Pe = NULL;		//Un-correlated candidate table

	/*=================
	FUNC PROTOTYPES :
	=================*/
	bool InputParams();
	bool Initialize();
	Fitness CalculateFitness(INDIVIDUAL const &individual);
	bool GetParents(int *indx_f, int *indx_m, const bool itsTime, bool *cross_frm_tables);
	bool CreateChildren(INDIVIDUAL & child1, INDIVIDUAL & child2, const int *indx_f, const int *indx_m, const bool itsTime, bool *cross_frm_tables);
	bool MutateChildren(INDIVIDUAL & child1);
	bool IdentifyChilds(INDIVIDUAL const & child1, INDIVIDUAL const & child2, const int indx_f, const int indx_m, int const indx);
	bool CopyPopulation(int gen_no);

	bool PopulateTables(int indx);
	int	 FindLeastCorrelatedCandidate(const int idx_Pb);
	float CalCorrelationCoff(INDIVIDUAL const & indv_main, INDIVIDUAL const & indv_obj);

	bool CalculateAvgFitness();
	double DecodeString(char * fChromo , const int i);
	int largestPowerOf2(const unsigned int n);
	bool BetterFit(const double fit_target, const double fit_base);
	int RouletteWheelSelection(INDIVIDUAL const * Popu, Fitness Start, Fitness End, int Popu_size);
	bool ShowPopu();
	bool ShowDude();
	bool ShowIndividual(INDIVIDUAL const &individual);
	bool ShowStatistics();

public:
	unsigned int Fitness_Calculations;				//Number of fittness calculations performed
	bool Minimizing_prob;							//Is this a minimizing problem
	Fitness optimum = 0.0;							//Optimum value of the function under consideration
	GenAlgo();
	~GenAlgo();
	bool Run(const int func_no);
	Fitness(*Objfunc)(const double *);

	void (*test_func)(double *, double *, int, int, int);
};