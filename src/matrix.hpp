

#ifndef GUARD_MATRIX
#define GUARD_MATRIX


#include <cstdlib>
#include <cstring>
#include <string>
#include <stdexcept>
#include <iostream>
#include <map>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include "limits.h"
#include <limits>

#include "tools.hpp"

typedef unsigned long ulong;
typedef unsigned char uchar;


using namespace std;


class Class_Substitution_Matrix
{
	public:

	int ** matrix;


	//int * self_scores;
	//ulong * counts;
	//ulong * neighbors;
	//int * neighbors_scores;

	//IMAPS * index;
	int aminoacid_alphabetSize;

	int *  aminoacid_ascii2num;
	char * aminoacid_num2ascii;
	ulong gapChar;

	//Class_Substitution_Matrix();
	Class_Substitution_Matrix(string * filename);

	~Class_Substitution_Matrix();

	int ** readMatrix(string * filename);
	int getScore(char a, char b);



	void printMatrix(int ** substitution_matrix);
	void checkMatrix(int ** substitution_matrix);


};




#endif
