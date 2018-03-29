#pragma once
#include "poisson1d.h"
#include <fstream>
class poisson2d
{
	//Assume N_x=N_y, h_x=h_y
	int N;
	double h;
	double** stiffnessMatrix;
	double** sparseStiffnessMatrix;
	double* solutionVector;
	double* rhs;
	std::ofstream myoutput;

public:
	poisson2d(int _N);
	~poisson2d();

	void assembleStiffnessMatrix();
	void assembleSparseStiffnessMatrix();

	void assemblerhs();

	void printStiffnessMatrix();
	void printSparseStiffnessMatrix();
	void printrhs();
	void Jacobi();

	void calculateDefect();
};

