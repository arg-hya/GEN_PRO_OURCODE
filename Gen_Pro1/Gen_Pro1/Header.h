#pragma once

#define SUCCESS		true
#define FAIL		false

#define MAXVECSIZE	30
#define MAXPOPSIZE	1000
#define MINPOPSIZE	2

#define RANDOM(min, max)  min + rand() % (max - min + 1) //range : [min, max] . Seeding done in constructor.
#define UNREFERENCED_PARAMETER(P) (P)


