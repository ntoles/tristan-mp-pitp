/*
This module calls functions from par.c, which is the input file parser from athena
Its purpose is to interface the fortran calls on inputparserf.F90 with the par
module, which is also writen in c

Currently implemented fortran interface public functions: 

par_open_f()
par_gets_def_f()
par_getd_def_f()
par_geti_def_f()
par_close_f()

*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "par.h"

#define MACSTYLE

#define CONCAT(prefix, name) prefix ## name
#ifdef MACSTYLE
	#define F90(name) CONCAT(name,_)
#else
	#define F90(name) name
#endif

// par_open_f - opens and reads the parameter file

void F90(par_open_f)(char *filename)
{
	par_open(filename);
}

// par_gets_def_f - returns a string from input file, or the default value

void F90(par_gets_def_f)(char *block, char *name, char *def, char *value, int len1, int len2, int len3, int len4)
{
	char *returnvalue;
	int i;
	
	returnvalue=par_gets_def(block,name,def);
	
	strncpy(value,returnvalue,len4<strlen(returnvalue)?len4:strlen(returnvalue));
	
	for (i=strlen(returnvalue); i<len4; i++)
		value[i]=' ';
}

// par_getd_def_f - returns a double from input file, or the default value

void F90(par_getd_def_f)(char *block, char *name, double *def, double *value)
{
	*value = par_getd_def(block,name,*def);
}

// par_geti_def_f - returns an integer from input file, or the default value

void F90(par_geti_def_f)(char *block, char *name, int *def, int *value)
{
	*value = par_geti_def(block,name,*def);
}

// par_close_f - closes the module and frees memory

void F90(par_close_f)(void)
{
	par_close();
}
