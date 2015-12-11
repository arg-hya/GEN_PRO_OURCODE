// Gen_Pro1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "GenAlgo.h"

#define square(x)  ((x)*(x))

void test_func(double *, double *, int, int, int);
double *OShift, *M, *y, *z, *x_bound;
int ini_flag = 0, n_flag, func_flag;

void OptimumValue(int function_number, double &optimum);

//#define prob2
//
//#ifdef prob2
//double objective(double const *x)
//{
//
//	double term1, term2, term3;
//	float g, penalty_coef;
//
//	//if (x == NULL) error_ptr_null("x in objective()");
//
//
//	term1 = (double)(x[0] + x[1]);
//	term2 = (double)((x[0] + x[1] * x[1] - 7.0)*(x[0] + x[1] * x[1] - 7.0));
//	term3 = term1 + term2;
//
//	penalty_coef = 0.0;
//	g = (float)((square(x[0] - 5.0) + square(x[1])) / 26.0 - 1.0);
//	if (g < 0.0) term3 = term3 + penalty_coef * g * g;
//	return (double)(term1);
//
//}
//#endif
//
//#ifdef prob1
//double objective(double const *x)
//{
//
//	double term1, term2, term3;
//	double g, penalty_coef;
//
//	//if (x == NULL) error_ptr_null("x in objective()");
//
//
//	term1 = (double)(square(x[0]) + square(x[1]));
//
//	term1 = (term1 == 0) ? term1 + 1 : term1;
//	
//	return (double)( 1/ term1);
//
//}
//#endif


int _tmain(int argc, _TCHAR* argv[])
{
	double optimum = 0.0;

	for (int i = 0; i < 28; i++)
	{
		std::cout << "************ Function number: " << i + 1 << " ************" << std::endl;

		OptimumValue(i + 1, optimum);
		//OptimumValue(1, optimum);

		GenAlgo *a = new GenAlgo();		

		a->test_func = test_func;		//Passes the objective function	

		a->Minimizing_prob = true;		//Mentioning this is a minimizing prob

		a->optimum = optimum;

		a->Run(i+1);			//Runs it
		//a->Run(1);			//Runs it

		a->~GenAlgo();		//Delete it
	}

	//getchar();

	return 0;
}

void OptimumValue(int function_number, double &optimum)
{
	//set optimal value
	switch (function_number) {
	case 1:
		optimum = -1400;
		break;
	case 2:
		optimum = -1300;
		break;
	case 3:
		optimum = -1200;
		break;
	case 4:
		optimum = -1100;
		break;
	case 5:
		optimum = -1000;
		break;
	case 6:
		optimum = -900;
		break;
	case 7:
		optimum = -800;
		break;
	case 8:
		optimum = -700;
		break;
	case 9:
		optimum = -600;
		break;
	case 10:
		optimum = -500;
		break;
	case 11:
		optimum = -400;
		break;
	case 12:
		optimum = -300;
		break;
	case 13:
		optimum = -200;
		break;
	case 14:
		optimum = -100;
		break;
	case 15:
		optimum = 100;
		break;
	case 16:
		optimum = 200;
		break;
	case 17:
		optimum = 300;
		break;
	case 18:
		optimum = 400;
		break;
	case 19:
		optimum = 500;
		break;
	case 20:
		optimum = 600;
		break;
	case 21:
		optimum = 700;
		break;
	case 22:
		optimum = 800;
		break;
	case 23:
		optimum = 900;
		break;
	case 24:
		optimum = 1000;
		break;
	case 25:
		optimum = 1100;
		break;
	case 26:
		optimum = 1200;
		break;
	case 27:
		optimum = 1300;
		break;
	case 28:
		optimum = 1400;
		break;
	}
}

