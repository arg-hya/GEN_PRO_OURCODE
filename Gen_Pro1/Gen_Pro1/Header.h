#pragma once



#define SUCCESS		true
#define FAIL		false

#define MAXVECSIZE	30
#define MAXPOPSIZE	1000
#define MINPOPSIZE	2

/*PARENT SELECTION FLAGS START*/
#define PARENT_FROM_POPU	0	//parents are from the main population
#define PARENT_FROM_PXB		1	//one parent is from the PXB population
#define PARENT_FROM_PE		2	//one parent is from the PE population
/*PARENT SELECTION FLAGS START*/


#define RANDOM(min, max)  min + rand() % (max - min + 1) //range : [min, max] . Seeding done in constructor. [INTEGER]

#define FLOAT_RANDOM(min, max)  min + (rand() / ( double(RAND_MAX) / (max - min )))  //range : [min, max] . Seeding done in constructor.[DOUBLE]
#define square(x)  ((x)*(x))
#define mod(x)  (x<0.0)? (-1 * x) : (x);

const double epsilon = pow(10.0, -8);

#define UNREFERENCED_PARAMETER(P) (P)


/*Selection , Crossover , Mutation Operators*/

//Crossover
#define One_Point_Crossover	1


