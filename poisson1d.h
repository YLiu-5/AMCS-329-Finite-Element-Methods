#pragma once
#include <iostream>
#include <fstream>
class poisson1d
{
protected:
	int N;
	double h;
	double tau;
	double** stiffnessMatrix;
	double** sparseStiffnessMatrix;
	double* solutionVector;
	double* rhs;
	double* Ulast;
	double bc;
	double norminf=0;
	double norm2=0;
	std::ofstream myoutput;

public:
	poisson1d();
	poisson1d(int _N,double _bc, double _tau);
	~poisson1d();
	void assembleStiffnessMatrix();
	void assembleSparseStiffnessMatrix();
	void assemblerhs();
	void solver();
	void printStiffnessMatrix();
	void printSparseStiffnessMatrix();
	void printrhs();
	void printsol();
	void printxOnehalf();
	void solcomp1plusx();
	void norms();
	double getNorm2();
	double getNormInf();
	void forwardEuler();
	void backwardEuler();

};

