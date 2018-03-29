#include "poisson1d.h"
#include "poisson2d.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <iomanip>
int main()
{
	std::vector<double> norm_inf(6);
	std::vector<double> norm_2(6);
	for(int k=4;k<=4;k++)
	{
		poisson1d P1d(pow(2,k),0,0.01);
		P1d.assembleStiffnessMatrix();
		P1d.assembleSparseStiffnessMatrix();
		P1d.printSparseStiffnessMatrix();
		P1d.assemblerhs();
		P1d.printrhs();
		P1d.solver();
		P1d.printsol();
		//P1d.norms();
		//norm_inf[k-2]=(P1d.getNormInf());
		//norm_2[k-2]=(P1d.getNorm2());
	}

/*	norm2Output.open("norm2.txt");
	for(int i=0;i<norm_2.size();i++)
	{
		norm2Output<<std::setw(15)<<std::left<<norm_2[i];
	}
	norm2Output.close();

	std::ofstream normInfOutput;
	normInfOutput.open("nornInf.txt");
	for (int i=0; i<norm_inf.size();i++)
	{
		normInfOutput<<std::setw(15)<<std::left<<norm_inf[i];
	}
	normInfOutput.close();*/
	
	//poisson1d P1d_2(8, 1);
	//P1d_2.assembleStiffnessMatrix();
	//P1d_2.assemblerhs();
	//P1d_2.solver();
	//P1d_2.printxOnehalf();
	//
	//poisson2d P2d(4);
	//P2d.assembleStiffnessMatrix();
	//P2d.printStiffnessMatrix();
	//P2d.assemblerhs();
	//P2d.printrhs();
	//P2d.calculateDefect();
	return 0;
}