// Gen_Pro1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "GenAlgo.h"
#include "Header.h"


/*************************************************************************/
/// <b>Function: GenAlgo constructor</b>
///
/// \remarks	Constructor.
/*************************************************************************/
GenAlgo::GenAlgo()
{	
	srand((unsigned int)time(NULL));	// initialize random seed		

	Minimizing_prob = false;

}

/*************************************************************************/
/// <b>Function: GenAlgo destructor</b>
///
/// \remarks	destructor.
/*************************************************************************/
GenAlgo::~GenAlgo()
{
	if (LW_BND)
	{
		delete[] LW_BND;
		LW_BND = NULL;
	}

	if (HG_BND)
	{
		delete[] HG_BND;
		HG_BND = NULL;
	}

	if (lenChromo_var)
	{
		delete[] lenChromo_var;
		lenChromo_var = NULL;
	}

	for (int i = 0; i < nPOPU; i++)
	{
		if (POPU[i].fChromo)
		{
			delete[] POPU[i].fChromo;
			POPU[i].fChromo = NULL;
		}
		if (Temp_popu[i].fChromo)
		{
			delete[] Temp_popu[i].fChromo;
			Temp_popu[i].fChromo = NULL;
		}
	}
	for (int i = 0; i < T; i++)
	{
		if (Pb[i].fChromo)
		{
			delete[] Pb[i].fChromo;
			Pb[i].fChromo = NULL;
		}
		if (Pxb[i].fChromo)
		{
			delete[] Pxb[i].fChromo;
			Pxb[i].fChromo = NULL;
		}
	}
	for (int i = 0; i < D; i++)
	{
		if (Pe[i].fChromo)
		{
			delete[] Pe[i].fChromo;
			Pe[i].fChromo = NULL;
		}
	}
	
	
	if (POPU) 
	{
		delete[] POPU;
		POPU = NULL;
	}
	if (Temp_popu)
	{
		delete[] Temp_popu;
		Temp_popu = NULL;
	}

	if (best_fitness_array)
	{
		delete[] best_fitness_array;
		best_fitness_array = NULL;
	}
	if (Pb)
	{
		delete[] Pb;
		Pb = NULL;
	}
	if (Pxb)
	{
		delete[] Pxb;
		Pxb = NULL;
	} 
	if (Pe)
	{
		delete[] Pe;
		Pe = NULL;
	}
	if (Xi)
	{
		delete[] Xi;
		Xi = NULL;
	}
}

/*************************************************************************/
/// <b>Function: InputParams</b>
///
/// \param  
///
/// \return		Returns True\False
///
/// \remarks	This function takes the parameters required from the user
///             returns true if everything goes well, otherwise returns false.
/*************************************************************************/
bool GenAlgo::InputParams()
{
	/*std::cout << "\nHow many generations ? ------------- : ";*/
	nGEN = 51;
	//std::cin >> nGEN;

	
	/*std::cout << "\nPopulation Size ? ------------------ : ";	*/
	nPOPU = 100;
	//std::cin >> nPOPU;
	
	if ((nPOPU > MAXPOPSIZE) || (nPOPU < MINPOPSIZE))
	{
		std::cout << "\n Change the value of MAXPOPSIZE in program and re-run the program";
		exit(EXIT_FAILURE);
	}

	//std::cout << "\nCross Over Probability ? ( 0 to 1 )  : ";	
	////std::cin >> Pc;
	Pc = 0.9;
	//std::cout << "\nMutation Probability ? ( 0 to 1 ) -- : ";	
	////std::cin >> Pm;
	Pm = 0.01;

	Ptb_xb = 0.6;
	Ptb_e = 0.05;
	//std::cout << "\nNumber of variables (Maximum %d) ---- : " << MAXVECSIZE;
	nVAR = 10;
	//std::cin >> nVAR;

	LW_BND = new float[nVAR];
	HG_BND = new float[nVAR];

	/**NOTE: Bounds are considered Rigid.*/
	for (int k = 0; k < nVAR; k++)
	{
		/*std::cout << "\nLower and Upper bounds of x[%d] ----- : ", k + 1;*/
		LW_BND[k] = -100;
		HG_BND[k] =  100;
	//	std::cin >> LW_BND[k];
	//	std::cin >> HG_BND[k];
	}

	lenChromo_var = new int[nVAR];

	/*std::cout << "\n Calculating Total string length (each variable has equal string length)";	*/
	for (int k = 0; k < nVAR; k++)
	{		
		lenChromo_var[k] = largestPowerOf2(HG_BND[k] - LW_BND[k]);
	}		

	//std::cout << "\n Give random seed (0 to 1.0)";	
	//std::cin >> BASIC_SEED;

	return SUCCESS;
}

/*************************************************************************/
/// <b>Function: Initialize</b>
///
/// \param  
///
/// \return		Returns True\False
///
/// \remarks	This function Creates and Initialises 0th Gen Population and initializes
///             some global variables, returns true if everything goes well, otherwise returns false.
/*************************************************************************/
bool GenAlgo::Initialize()
{
	int i ,j;

	/*NOTE:: lenChromo_tot = Sum(lenChromo_var) . Thus all the variables are 
	*	     considered in a single choromosome.							   */
	lenChromo_tot = 0;
	for (int k = 0; k < nVAR; k++)	lenChromo_tot += lenChromo_var[k];
	
	MaxFitness_Calculations = 10000 * nVAR;		//According to CEC 2013 benchmark.
	
	POPU = new INDIVIDUAL[nPOPU];	
	Temp_popu = new INDIVIDUAL[nPOPU];

	best_fitness_array = new Fitness[nGEN];

	T = (int)(nGEN * 0.2);
	D = (int)(nGEN * 0.1);	//NOTE:: D should be always less than 2*T
	
	//Initializing tables
	Pb = new INDIVIDUAL[T];
	Pxb = new INDIVIDUAL[T];
	Pe = new INDIVIDUAL[D];

	Xi = new Fitness[nVAR];

	for (i = 0; i < nPOPU; i++)
	{
		POPU[i].fChromo = new char[lenChromo_tot];
		Temp_popu[i].fChromo = new char[lenChromo_tot];
		
		for (j = 0; j < lenChromo_tot; j++)
		{
			POPU[i].fChromo[j] = (rand() % 2) ;   // randomly initialises the chorosome string (mod2 random number generation)
		}	
		
		POPU[i].fFitness = CalculateFitness(POPU[i]);	//calculates the fittness level.
	}
	for (i = 0; i < T; i++)
	{
		Pb[i].fChromo = new char[lenChromo_tot];
		Pxb[i].fChromo = new char[lenChromo_tot];
	}
	for (i = 0; i < D; i++)
	{
		Pe[i].fChromo = new char[lenChromo_tot];
	}


	/*Populating Xi start*/
	for (int j = 0; j < nVAR; j++)
	{
		Xi[j] = 0.50;   // Initially Xi is set ro 0.5. i.e every variable is equally important.
	}
	/*Populating Xi end*/

	return SUCCESS;
	
}


/*************************************************************************/
/// <b>Function: GetParents</b>
///
/// \param  
///
/// \return		Returns True\False
///
/// \remarks	This function selects two different parents from the population
///             and stores their locations in indx_m and indx_f respectively.
/*************************************************************************/
bool GenAlgo::GetParents(int *indx_f, int *indx_m, const bool itsTime, byte *cross_frm_tables)
{
	double 	bCrossTables_Pxb = 0.0,		//decides whether to populate from the Pxb Table.
			bCrossTables_Pe = 0.0;		//decides whether to populate from the Pe Table.
	
	bCrossTables_Pxb = FLOAT_RANDOM(0, 1);
	bCrossTables_Pe = FLOAT_RANDOM(0, 1);

	*cross_frm_tables = PARENT_FROM_POPU;

	if (itsTime)
	{
		/*if (bCrossTables_Pxb <= Ptb_xb)
		{
			*cross_frm_tables = PARENT_FROM_PXB;
			*indx_m = RANDOM(0, T  - 1);
		}*/
		 if (bCrossTables_Pe <= Ptb_e)
		{
			*cross_frm_tables = PARENT_FROM_PE;
			*indx_m = RANDOM(0,  D - 1);
		}
		else if (bCrossTables_Pxb <= Ptb_xb)
		{
			*cross_frm_tables = PARENT_FROM_PXB;
			*indx_m = RANDOM(0, T - 1);
		}
	}
	if (PARENT_FROM_POPU == *cross_frm_tables)
	{		
		*indx_m = AboveAverageSelection(POPU, nPOPU);
				  //RouletteWheelSelection(POPU, optimum, nPOPU);
				 //RANDOM(0, nPOPU - 1);
	}	

	while (1)
	{
		*indx_f = AboveAverageSelection(POPU, nPOPU);
				 //RouletteWheelSelection(POPU, optimum, nPOPU);
				//RANDOM(0, nPOPU - 1);

		if ((!itsTime) ||
			(false == *cross_frm_tables)
			)
		{
			if (*indx_f != *indx_m)	break;
		}
		else break;
	}
	return SUCCESS;
}

/*************************************************************************/
/// <b>Function: CreateChildren</b>
///
/// \param  
///
/// \return		Returns True\False
///
/// \remarks	This function creates two children from the parents by randomly
///				selecting a crossover point and then merging the parents w.r.t the selected point.
/*************************************************************************/
bool GenAlgo::CreateChildren(INDIVIDUAL & child1, INDIVIDUAL & child2,
			  const int *indx_f, const int *indx_m, const bool itsTime, byte *cross_frm_tables)
{
	int crossPoint = 0,
		i;
	double bCrossover = 0.0;		//decides if crossover should take place	   
	
	bCrossover = FLOAT_RANDOM(0, 1);   	//generating random real number [0,1]	

	crossPoint = RANDOM(0, lenChromo_tot - 1);	
	 
	if (bCrossover <= Pc)
	{
		//CrossOver and create children
		if ((itsTime) && 
			((PARENT_FROM_PXB == *cross_frm_tables) || (PARENT_FROM_PE == *cross_frm_tables))	//One parent from tables
			)
		{			
			int indx_tables = *indx_m;	//indx_m contains the index from the table.//RANDOM(0, (T + D - 1));
			
			const char *fChromo_table;	//contains the chromo from the considered table
			fChromo_table = (PARENT_FROM_PXB == *cross_frm_tables) ? Pxb[indx_tables].fChromo : Pe[indx_tables].fChromo;

#if One_Point_Crossover 
				for (i = 0; i <= crossPoint; i++)
				{
					child1.fChromo[i] = POPU[*indx_f].fChromo[i];
					child2.fChromo[i] = fChromo_table[i];
				}
				for (i = crossPoint + 1; i < lenChromo_tot; i++)
				{
					child1.fChromo[i] = fChromo_table[i];
					child2.fChromo[i] = POPU[*indx_f].fChromo[i];
				}
				
#endif
		}
		else   //No parent from tables
		{
#if One_Point_Crossover 
			for (i = 0; i <= crossPoint; i++)
			{
				child1.fChromo[i] = POPU[*indx_f].fChromo[i];
				child2.fChromo[i] = POPU[*indx_m].fChromo[i];
			}
			for (i = crossPoint + 1; i < lenChromo_tot; i++)
			{
				child1.fChromo[i] = POPU[*indx_m].fChromo[i];
				child2.fChromo[i] = POPU[*indx_f].fChromo[i];
			}
		}
#endif
	}
	else   //Crossover is not taking place
	{
		for (i = 0; i < lenChromo_tot; i++)	//Copying parent chromo into child chromo
		{
			child1.fChromo[i] = POPU[*indx_f].fChromo[i];
			child2.fChromo[i] = POPU[*indx_m].fChromo[i];
		}
	}


	return SUCCESS;
	
}

/*************************************************************************/
/// <b>Function: MutateChildren</b>
///
/// \param  
///
/// \return		Returns True\False
///
/// \remarks	This function mutates the offspring at each locus with probability Pm (the mutation rate).
/*************************************************************************/
bool GenAlgo::MutateChildren(INDIVIDUAL & child)
{
	int mutationPoint = 0 , i;
	//	i;
	double bmutate = 0.0;			//decides if the child should be mutated

	//bmutate = FLOAT_RANDOM(0, 1); // rand() % 2;	//mod2 random number generation

	//if (bmutate <= Pm)
	//{
	//	mutationPoint = RANDOM(0, lenChromo_tot - 1);	//selects a random mutation point		

	//	child.fChromo[mutationPoint] = (child.fChromo[mutationPoint] + 1) % 2;	//fliping the value
	//}

	for (i = 0; i < lenChromo_tot; i++)
	{
		bmutate = FLOAT_RANDOM(0, 1); // rand() % 2;	//mod2 random number generation

		if (bmutate <= Pm)
		{
			//mutationPoint = RANDOM(0, lenChromo_tot - 1);	//selects a random mutation point		

			child.fChromo[i] = (child.fChromo[i] + 1) % 2;	//fliping the value
		}
	}
	return SUCCESS;

}

/*************************************************************************/
/// <b>Function: IdentifyChilds</b>
///
/// \param  
///
/// \return		Returns True\False
///
/// \remarks	This function  identifies the childs which are fitter than their parents and the stores them in a temporary population.
///				if a child is unfit then it stores the fit parents
/*************************************************************************/
bool GenAlgo::IdentifyChilds(INDIVIDUAL const & child1, INDIVIDUAL const & child2,
							 const int indx_f, const int indx_m, int const indx)
{	int i;	
	INDIVIDUAL Contestants[4] = { child1, child2, POPU[indx_f], POPU[indx_m] };		//Table of Contestants
	
	int Indi1_idx = 0, Indi2_idx = 0;
	for (i = 0; i < _countof(Contestants); i++)		//Finding the fittest individual indx
	{		
		if (BetterFit(Contestants[i].fFitness, Contestants[Indi1_idx].fFitness))	Indi1_idx = i;		
	}

	for (i = 0; i < _countof(Contestants); i++)		//Finding the second fittest individual indx
	{
		if (Indi1_idx != i)	//Both must not be the same individual
		{
			if (BetterFit(Contestants[i].fFitness, Contestants[Indi2_idx].fFitness))	Indi2_idx = i;
		}		
	}

	/*Copying First individual*/
	for (i = 0; i < lenChromo_tot; i++)
		Temp_popu[indx].fChromo[i] = Contestants[Indi1_idx].fChromo[i];

	Temp_popu[indx].fFitness = Contestants[Indi1_idx].fFitness;

	/*Copying Second individual*/
	for (i = 0; i < lenChromo_tot; i++)
		Temp_popu[indx + 1].fChromo[i] = Contestants[Indi2_idx].fChromo[i];

	Temp_popu[indx + 1].fFitness = Contestants[Indi2_idx].fFitness;

	return SUCCESS;

}

/*************************************************************************/
/// <b>Function: CopyPopulation</b>
///
/// \param  
///
/// \return		Returns True\False
///
/// \remarks	Copies population from temp population to new population,
///				and stores the best fitness level achieved in this generation.
/*************************************************************************/
bool GenAlgo::CopyPopulation(int gen_no)
{
	int max_fit_indx = 0;

	for (int i = 0; i < nPOPU; i++)
	{
		for (int j = 0; j < lenChromo_tot; j++)
		{
			POPU[i].fChromo[j] = Temp_popu[i].fChromo[j];   // copies the chromosome string from temp population to new population.
		}

		POPU[i].fFitness = Temp_popu[i].fFitness;	//copies the fittness level.
	}

	/*Stores the best fitness level achieved in this generation in best_fitness_array[]. */	  
	for (int i = 0; i < nPOPU; i++)
	{		
		if (BetterFit(POPU[i].fFitness , POPU[max_fit_indx].fFitness))	max_fit_indx = i;
	}
	best_fitness_array[gen_no] = POPU[max_fit_indx].fFitness;
	
	return SUCCESS;
}

/*************************************************************************/
/// <b>Function: PopulateTables</b>
///
/// \param  
///
/// \return		Returns True\False
///
/// \remarks	Populates the tables.
/*************************************************************************/
bool GenAlgo::PopulateTables(int indx)
{
	int idx_Pb = 0;				//Stores the current best candidate indx
	float maxFittness = POPU[0].fFitness;	

	/**********Populating Xi start**********/
	//NOTE:: Every generation starts with a fresh value of Xi[].
	for (int j = 0; j < nVAR; j++)
	{
		Xi[j] = 0.50;   // Initially Xi is set ro 0.5. i.e every variable is equally important.
	}
	/**********Populating Xi end**********/


	/**********Populating Pb start**********/
	for (int i = 0; i < nPOPU; i++)
	{
		if (BetterFit(POPU[i].fFitness , maxFittness))	idx_Pb = i;	//finding the current best candidate
	}

	for (int j = 0; j < lenChromo_tot; j++)
	{
		Pb[(indx % T)].fChromo[j] = POPU[idx_Pb].fChromo[j];   // copies the chromosome string to Pb and discards the oldest if indx > T-1 .
	}
	Pb[(indx % T)].fFitness = POPU[idx_Pb].fFitness;
	/**********Populating Pb end**********/

	/**********Populating Pxb start**********/
	//Right now Pxb is the average of Pb
	if (indx >= T)	//Pb is full
	{	
		//int indx_Pxb = RANDOM(0, T-1);
		int index = (indx % T);

		for (int j = 0; j < lenChromo_tot; j++)	//Setting all the genes to zero.
		{
			Pxb[index].fChromo[j] = 0;
		}		
		for (int i = 0; i < nVAR; i++)
		{
			double numerator_sum = 0.0, denom_sum = 0.0;
			int ran_indx = RANDOM(0, T - 2);

			numerator_sum = ((Pb[ran_indx].fFitness - optimum)*(DecodeString(Pb[ran_indx].fChromo, i))) + ((Pb[ran_indx + 1].fFitness - optimum)*(DecodeString(Pb[ran_indx + 1].fChromo, i)));
			denom_sum = (Pb[ran_indx].fFitness - optimum) + (Pb[ran_indx + 1].fFitness - optimum);
			numerator_sum /= denom_sum;			

			EncodeString(Pxb[index].fChromo, i, numerator_sum);
		}

		Pxb[index].fFitness = CalculateFitness(Pxb[index]);
	}
	/**********Populating Pxb end**********/

	/*********Populating Pe start**********/
	int indx_Pe = FindLeastCorrelatedCandidate(idx_Pb);
	for (int j = 0; j < lenChromo_tot; j++)
	{
		Pe[(indx % D)].fChromo[j] = POPU[indx_Pe].fChromo[j];   // copies the chromosome string to Pe and discards the oldest if indx > D-1 .
	}
	Pe[(indx % D)].fFitness = POPU[indx_Pe].fFitness;
	/**********Populating Pe end**********/	

	return SUCCESS;
}

/*************************************************************************/
/// <b>Function: FindLeastCorrelatedCandidate</b>
///
/// \param  
///
/// \return		Returns True\False
///
/// \remarks	Finds the least correlated candidate w.r.t the current best candidate, whose fitness is above avg.
/*************************************************************************/
int GenAlgo::FindLeastCorrelatedCandidate(const int idx_Pb)
{
	float corr_min = 1, corr;			//Stores the minimum corr. coeff. found
	int indx = idx_Pb;					//Stores the index of the minimum corr. coeff. found

	for (int i = 0; i < nPOPU; i++)
		{

		if (BetterFit(POPU[i].fFitness , avgRunningFITNESS) && (i != idx_Pb))
			{
				corr = CalCorrelationCoff(POPU[idx_Pb], POPU[i]);
				if (corr < corr_min)
				{
					corr_min = corr;
					indx = i;
				}
			}
		}

	return indx;
}

/*************************************************************************/
/// <b>Function: CalCorrelationCoff</b>
///
/// \param  
///
/// \return		Returns True\False
///
/// \remarks	Copies population from temp population to new population.
///				"indv_main" is the base and "indv_obj" is the comparable
/*************************************************************************/
float GenAlgo::CalCorrelationCoff(INDIVIDUAL const & indv_main, INDIVIDUAL const & indv_obj)
{
	double corr = 0.0, diff = 0, norm_diff, norm_diff_fitness, factor = -1, rms_square = 0.0, mod_highest_fitness = 0.0;
	
	/*for (int i = 0; i < nPOPU; i++)
	{
		if (mod_highest_fitness < POPU[i].fFitness)
		{
			mod_highest_fitness = POPU[i].fFitness;		
		}
	}	*/
	norm_diff_fitness = mod((indv_main.fFitness - indv_obj.fFitness));		
	//norm_diff_fitness /= mod((mod_highest_fitness - optimum));
	
	for (int i = 0; i < nVAR; i++)
	{
		
		/*Calcualting Correlation Coefficient*/
		diff = mod((DecodeString(indv_main.fChromo, i) - DecodeString(indv_obj.fChromo, i)));
		corr += Xi[i] * square(diff);  
		rms_square += square(diff);

		/*Updating Xi table*/
		if (0 != norm_diff_fitness)	//division by zero considered 2)
		{
			norm_diff = diff;// / (HG_BND[i] - LW_BND[i]);
			factor = norm_diff / norm_diff_fitness;
		}

		/*NOTE:: if factor is -1, then the fitness values are equal.*
		*In this case we are not changing the weight of the variable*/
		if (-1 != factor)	
		{
			if (factor < 1)	Xi[i] += ((1 - Xi[i]) / pow(2, factor ));
			else
			{				
				Xi[i] -= (Xi[i] / pow(2, (1 / factor)));
			}
		}		
	}

	corr /= rms_square;
	corr = sqrt(corr);

	return (1 - corr);
}

/*************************************************************************/
/// <b>Function: CopyPopulation</b>
///
/// \param  
///
/// \return		Returns fittness
///
/// \remarks	Calculates fittness of the individual.
/*************************************************************************/
GenAlgo::Fitness GenAlgo::CalculateFitness(INDIVIDUAL const &individual)
{
	double *x = new double[nVAR],
		result = 0;
	Fitness rslt;

	for (int i = 0; i < nVAR; i++)
	{			
		x[i] = DecodeString(individual.fChromo , i);
	}	
	
	//exit(EXIT_FAILURE);
		
	test_func(x, &rslt, nVAR, 1, func_num);

	//rslt =  Objfunc(x);
	if (x)
	{
		delete[] x;
		x = NULL;
	}

	return rslt;
}

/*************************************************************************/
/// <b>Function: DecodeString</b>
///
/// \param  
///
/// \return		Returns decoded string
///
/// \remarks	Binary to BCD conversion. [LSB...MSB]
/*************************************************************************/
double GenAlgo::DecodeString(char * fChromo, const int i)
{
	int temp = 0,
		start = 0,
		end = 0,
		var_highest = (int)pow(2, lenChromo_var[i]);

	double coff,
		   result = 0.0 ;

	for (int k = 0; k < i; k++)	start += lenChromo_var[k];
		
	end = start + lenChromo_var[i];
	

		for (int j = start; j < end; j++)
	{
			temp |= (fChromo[j] == 1) << (j - start);
	}
	
	///*NOTE:: Boundary conditions are taken care here and ONLY here.*/
		if (true)
		{
			coff = (double)(temp / (double)var_highest);
			result = LW_BND[i] + coff * (HG_BND[i] - LW_BND[i]);
		}

		if (false)
		{
			UNREFERENCED_PARAMETER(var_highest);
			result = (double)(temp % (unsigned int)(HG_BND[i] - LW_BND[i]));
			result = result - ((HG_BND[i] - LW_BND[i]) / 2);
		}

	return result;
		
}

/*************************************************************************/
/// <b>Function: DecodeString</b>
///
/// \param  
///
/// \return		Returns decoded string
///
/// \remarks	 BCD to Binary conversion. [LSB...MSB]
/*************************************************************************/
bool GenAlgo::EncodeString(char * fChromo, const int i , double value)
{
	int temp = 0,
		start = 0,
		end = 0,
		var_highest = (int)pow(2, lenChromo_var[i]);

	double coff;		

	for (int k = 0; k < i; k++)	start += lenChromo_var[k];

	end = start + lenChromo_var[i];
	
	coff = (value - LW_BND[i]) / (HG_BND[i] - LW_BND[i]);
	temp = (int)(coff *(double)var_highest);
	

	for (int j = start; j < end; j++)
	{
		if (temp & (1 << (j - start)))
		{
			fChromo[j] = 1;
		}
	}	


	return SUCCESS;

}

/*************************************************************************/
/// <b>Function: largestPowerOf2</b>
///
/// \param  
///
/// \return		Returns decoded string
///
/// \remarks	largest power of 2 in the given number.
/*************************************************************************/
int GenAlgo::largestPowerOf2(const unsigned int n)
{
	int res = 2, 
		count = 0;
	while (res < n) {
		res *= 2;
		count++;
	}

	return count + 1 ;

	//// given an int number (which would be the upper value of the range)
	//number |= number >> 1;
	//number |= number >> 2;
	//number |= number >> 4;
	//number |= number >> 8;
	//number |= number >> 16;
	//return (number >> 1) + 1;

}
/*************************************************************************/
/// <b>Function: BetterFit</b>
///
/// \param  
///
/// \return		Returns true/false
///
/// \remarks	Returns the better fitt candidate. fit_target is compared with fit_base
/*************************************************************************/
bool GenAlgo::BetterFit(const double fit_target, const double fit_base)
{
	if (Minimizing_prob)
	{
		if (fit_target <= fit_base)	return true;
		else return false;
	}
	else
	{
		if (fit_target >= fit_base)	return true;
		else return false;
	}	
}

/*************************************************************************/
/// <b>Function: CalculateAvgFitness</b>
///
/// \param  
///
/// \return		Returns decoded string
///
/// \remarks	Calculates the running average of the current population.
/*************************************************************************/
bool GenAlgo::CalculateAvgFitness()
{
	avgRunningFITNESS = 0;
	for (int i = 0; i < nPOPU; i++)
	{
		avgRunningFITNESS += POPU[i].fFitness;
	}
	avgRunningFITNESS /= nPOPU;

	return SUCCESS;
}

/*************************************************************************/
/// <b>Function: RouletteWheelSelection</b>
///
/// \param  
///
/// \return		Returns index of the selected individual
///
/// \remarks	Roulette Wheel Selection is done.
/*************************************************************************/
int GenAlgo::RouletteWheelSelection(INDIVIDUAL const * Popu, Fitness Start, int Popu_size)
{
	int i;
	Fitness sum = 0.0,
			fitness_indx = 0.0,
			mod_highest_fitness = 0.0,		//highest value of fittness encountered stored here [if (true == Minimizing_prob)]
			End = 0.0 ;						//Range of rand, i.e	0  to Rand_Max
	
	if (true == Minimizing_prob)
	{
		for (i = 0; i < Popu_size; i++)
			{		
				if (mod_highest_fitness < (POPU[i].fFitness - Start))
				{
					mod_highest_fitness = POPU[i].fFitness - Start;		//shifting origin to 0
				}					
			}		
	}

	for (i = 0; i < Popu_size; i++)
	{
		if (true == Minimizing_prob)
		{			
			End += (mod_highest_fitness - (POPU[i].fFitness - Start));		
		}
		else
		{
			End += POPU[i].fFitness - Start;		
		}
	}

	fitness_indx = FLOAT_RANDOM(0.0, End);		//Range [0  to  Rand_max]

	for (i = 0; i < Popu_size - 1; i++)
	{
		if (true == Minimizing_prob)
		{			
			sum += (mod_highest_fitness - (POPU[i].fFitness - Start));		//flipping fitness values in case of Minimizing_prob
			if (sum > fitness_indx) break;
		}
		else
		{
			sum += (POPU[i].fFitness - Start);
			if (sum > fitness_indx) break;
		}		
	}	

	return i ;	
}

/*************************************************************************/
/// <b>Function: AboveAverageSelection</b>
///
/// \param  
///
/// \return		Returns index of the selected individual
///
/// \remarks	Randomly selects a candidate whose fitness is above the average fitness, then finds its 
///				correlation coff. w.r.t the best candidate. Returns the highest correlated candidate's index.
///				No of iterations done is  Popu_size / 8
/*************************************************************************/
int GenAlgo::AboveAverageSelection(INDIVIDUAL const * Popu, int Popu_size)
{
	int i,
		j = 0,
		random_indx = 0,
		indx = -1,
		current_best_fit_indx = 0;
	float corr = 0.0,
		 tempcorr = 0.0;
	  
	
	/*Finds the best fitness level achieved in this generation. */
	//NOTE:: index can also be taken best_fitness_array[genNo], this loop is not required then
	for (i = 0; i < Popu_size; i++)
	{
		if (BetterFit(POPU[i].fFitness, POPU[current_best_fit_indx].fFitness))	current_best_fit_indx = i;
	}

	do
	{
		j++;

		random_indx = RANDOM(0, Popu_size - 1);
		
		if (BetterFit(POPU[random_indx].fFitness, avgRunningFITNESS))
		{
			tempcorr = CalCorrelationCoff(Popu[current_best_fit_indx], Popu[random_indx]);			
			if (tempcorr > corr)		//more correlation with the best candidate the better
			{
				corr = tempcorr;
				indx = random_indx;
			}

		}	
		
	} while (j < Popu_size/8);

	return (indx != -1) ? indx : RANDOM(0 , Popu_size - 1);
}

/*************************************************************************/
/// <b>Function: Run</b>
///
/// \param  
///
/// \return		Returns True\False
///
/// \remarks	This function takes the parameters required from the user
///             returns true if everything goes well, otherwise returns false.
/*************************************************************************/
bool GenAlgo::Run(const int func_no)
{

	int nPopu,
		totGen,
		indx_f = 0,				//Index of one parent.	
		indx_m = 0;				//Index of another parent.

	byte cross_frm_tables = PARENT_FROM_POPU;

	bool rslt,		 
		 itsTime = false;		//flag to indicate the start of our algo
	

	func_num = func_no;

	/*INPUT PARAMETERS FROM USER*/
	rslt = InputParams();

	totGen = nGEN;
	Fitness_Calculations = 0;	//Setting the fittness calculations count to 0.

	/*INITIALIZE THE 0TH GEN POPULATION*/
	Initialize(); 
	Fitness_Calculations += nPOPU;

	do
	{		
		CHILDREN child1, child2;	//This are used to store the children temporaryly
		nPopu = nPOPU;	

		//ShowPopu();

		if (Fitness_Calculations > MaxFitness_Calculations)	break;		//Maximum number of fittness calculations has occured.

		if ((nGEN - totGen) > 2 * T)	itsTime = true;	

		CalculateAvgFitness();

		while ((nPopu -= 2) >= 0)
		{
			/*allocation memory*/
			child1.fChromo = new char[lenChromo_tot];
			child2.fChromo = new char[lenChromo_tot];
			child1.fFitness = child2.fFitness = NULL;


			/*SELECT PARENTS*/
			GetParents(&indx_f, &indx_m, itsTime, &cross_frm_tables);

			/*DO CROSSOVER*/
			CreateChildren(child1, child2, &indx_f, &indx_m, itsTime, &cross_frm_tables);

			/*MUTATES THE CHILDREN*/
			MutateChildren(child1);
			MutateChildren(child2);

			/*CALCULATE FITNESS OF OFFSPRINGS*/
			child1.fFitness = CalculateFitness(child1);	
			child2.fFitness = CalculateFitness(child2);
			Fitness_Calculations += 2;

			/*IDENTIFY GOOD CHILDS AND STORE THEM*/
			IdentifyChilds(child1, child2, indx_f, indx_m, nPOPU - nPopu - 2);


			/*cleaning up*/
			indx_f = indx_m;
			delete[] child1.fChromo;
			delete[] child2.fChromo;
			child1.fFitness = child2.fFitness = NULL;
		}

		/*UPDATE NEW GEN POPULATION*/
		CopyPopulation(nGEN - totGen);

		/*POPULATE TABLES*/
		PopulateTables(nGEN - totGen);
		Fitness_Calculations++;

		/*std::cout << "\nGen no :: " << nGEN - totGen << std::endl;
		std::cout << "**************" << std::endl;*/
		//ShowDude();

	} while (--totGen);

	//ShowStatistics();

	ShowDude();		

return true;
}

bool GenAlgo::ShowStatistics()
{
	mean_best_fitness = std_best_fitness = 0;

	for (int j = 0; j < nGEN; j++) {
		mean_best_fitness += best_fitness_array[j];
	}

	mean_best_fitness /= nGEN;

	for (int j = 0; j < nGEN; j++) {
		std_best_fitness += square(mean_best_fitness - best_fitness_array[j]);
	}

	std_best_fitness /= nGEN;
	std_best_fitness = sqrt(std_best_fitness);

	std::cout << "\nMean = " << mean_best_fitness << std::endl;
	std::cout << "std = " << std_best_fitness << std::endl;

	return SUCCESS;
}

bool GenAlgo::ShowIndividual(INDIVIDUAL const &individual)
{
	/*for (int j = 0; j < nVAR; j++)
	{
		std::cout << "Var " << j << " : " << DecodeString(individual.fChromo, j) << " " << std::endl;
	}*/
	std::cout << "Fittness : " << individual.fFitness << std::endl;
	std::cout << "\n Optimum Fittness : " << optimum << std::endl;
	return SUCCESS;
}

bool GenAlgo::ShowDude()
{ 
	int max_fit_indx = 0;

	for (int i = 0; i < nPOPU; i++)
	{
		if (BetterFit(POPU[i].fFitness, POPU[max_fit_indx].fFitness))	max_fit_indx = i;
	}
	
	//std::cout << "\nDude is : " << std::endl;
	ShowIndividual(POPU[max_fit_indx]);
	
	return SUCCESS;
}

bool GenAlgo::ShowPopu()
{
	for (int i = 0; i < nPOPU; i++)
	{
		std::cout << "*************" << "Individual " << i << "*************" << std::endl;
		ShowIndividual(POPU[i]);		
	}

	return SUCCESS;
}