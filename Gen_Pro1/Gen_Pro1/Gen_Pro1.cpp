// Gen_Pro1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "GenAlgo.h"

#define square(x)  ((x)*(x))

void test_func(double *, double *, int, int, int);
double *OShift, *M, *y, *z, *x_bound;
int ini_flag = 0, n_flag, func_flag;

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


int main(int argc, _TCHAR* argv[])
{
	for (int i = 0; i < 28; i++)
	{
		std::cout << "************ Function number: " << i + 1 << " ************" << std::endl;
		GenAlgo *a = new GenAlgo();		

		a->test_func = test_func;		//Passes the objective function	

		a->Minimizing_prob = true;		//Mentioning this is a minimizing prob

		a->Run(i+1);			//Runs it

		a->~GenAlgo();		//Delete it
	}

	getchar();

	return 0;
}

