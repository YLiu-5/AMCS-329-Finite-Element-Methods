#include "poisson1d.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <cmath>
#include <fstream>
#include <string>
poisson1d::poisson1d()
{
}

poisson1d::poisson1d(int _N,double _bc, double _tau)
{
	N = _N;
	tau=_tau;
	h = 1.0 / N;
	bc = _bc;
	//allocate memory
	//totally 0,1,2,...N grid points, but only N-1 unknowns. In the case N=128, we have 127 unknowns, the indices are 0,1,2,...,N-2
	//For Exercise 3.3 there are N unknowns because of Neumann BC.
	stiffnessMatrix = new double*[N-1];
	for (int i = 0; i < N-1; i++)
		stiffnessMatrix[i] = new double[N-1];
	sparseStiffnessMatrix = new double*[3];
	for (int i = 0; i < 3; i++)
		sparseStiffnessMatrix[i] = new double[N-1];
	
	solutionVector = new double[N-1];
	for (int i = 0; i < N-1; i++)
		solutionVector[i] = 0;
	//Initialize

	
	rhs = new double[N-1];	//P1d.printStiffnessMatrix();
	

	std::string filename;
	filename = "poisson1d_output.txt";
	myoutput.open(filename.c_str(),std::ios_base::out);
	myoutput<< "\n\n";
	myoutput<< "Case: N="<<N<<std::endl;

}


poisson1d::~poisson1d()
{
}

void poisson1d::assembleStiffnessMatrix()
{

	//initialization
	for(int i = 0;i < N-1; i++)
		for (int j = 0; j < N-1; j++)
		{
			stiffnessMatrix[i][j] = 0.0;
		}

	for (int i = 1;i < N - 2; i++)
		for (int j = 0; j < N - 1; j++)
		{
			if (j == i)
				stiffnessMatrix[i][j] = 2.0;
			else if (j - i == 1)
				stiffnessMatrix[i][j] = -1.0;
			else if (j - i == -1)
				stiffnessMatrix[i][j] = -1.0;
			else
				stiffnessMatrix[i][j] = 0.0;
		}
	//first row and last row
	stiffnessMatrix[0][0] = 2.0;
	stiffnessMatrix[0][1] = -2.0;
	stiffnessMatrix[N-2][N-2] = 2.0;
	stiffnessMatrix[N-2][N-3] = -1.0;
}

void poisson1d::assembleSparseStiffnessMatrix()
{
	//Initialization
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < N-1; j++)
		{
			sparseStiffnessMatrix[i][j] = 0.0;
		}

	for (int k = 0; k < N; k++)
		sparseStiffnessMatrix[1][k] = 2;
	for (int k = 0; k < N; k++)
	{
		if (k == 0)
			sparseStiffnessMatrix[2][k] = -1;
		else if (k == N - 2)
			sparseStiffnessMatrix[0][k] = -1;
		else
		{
			sparseStiffnessMatrix[0][k] = -1;
			sparseStiffnessMatrix[2][k] = -1;
		}
	}
}

void poisson1d::assemblerhs()
{

	for (int i = 0; i < N-1; i++)
	{
		rhs[i]=1*pow(h,2);
		//rhs[i] = sin(M_PI*(i+1)*h) * pow(h,2);
		//rhs[i] = 1* pow(h,2);
	}
	//boundary condition
	rhs[0] += 0;
	rhs[N - 2] += 0;

}


void poisson1d::solver()
{
	//Thomas Algorithm
	double ratio = 0;
	for (int k = 1; k < N - 1; k++)
	{
		//ratio = stiffnessMatrix[k][k - 1] / stiffnessMatrix[k - 1][k - 1];
		ratio = sparseStiffnessMatrix[0][k] / sparseStiffnessMatrix[1][k - 1];
		//stiffnessMatrix[k][k] -= stiffnessMatrix[k - 1][k] * ratio;
		sparseStiffnessMatrix[1][k] -= sparseStiffnessMatrix[2][k - 1] * ratio;
		sparseStiffnessMatrix[0][k] -= sparseStiffnessMatrix[1][k - 1] * ratio;
		rhs[k] -= rhs[k - 1] * ratio;
	}
	//solve by backward substitution

	//solutionVector[N - 2] = rhs[N - 2] / stiffnessMatrix[N - 2][N - 2];
	solutionVector[N - 2] = rhs[N - 2] / sparseStiffnessMatrix[1][N - 2];
	for (int k = N - 3; k >= 0; k--)
	{
		//solutionVector[k] = (rhs[k] - stiffnessMatrix[k][k + 1] * solutionVector[k + 1]) / stiffnessMatrix[k][k];
		solutionVector[k] = (rhs[k] - sparseStiffnessMatrix[2][k] * solutionVector[k + 1]) / sparseStiffnessMatrix[1][k];
	}

}

void poisson1d::printStiffnessMatrix()
{
	myoutput << "Stiffness Matrix: " << std::endl;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < N - 1; j++)
		{
			myoutput <<std::setw(5) << std::left << stiffnessMatrix[i][j];
		}
		myoutput << std::endl;
	}
}

void poisson1d::printSparseStiffnessMatrix()
{
	myoutput << "Sparse Stiffness Matrix: " << std::endl;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < N - 1; j++)
		{
			myoutput << std::setw(5) << std::left << sparseStiffnessMatrix[i][j];
		}
		myoutput << std::endl;
	}
}

void poisson1d::printrhs()
{

	myoutput << "RHS: " << std::endl;
	for (int i = 0; i < N - 1; i++)
	{
		myoutput << std::setw(10) << std::left << rhs[i];
	}
	myoutput << std::endl;

}

void poisson1d::printsol()
{
	myoutput << "Sol: " << std::endl;
	for (int i = 0; i < N - 1; i++)
	{
		myoutput << std::setw(10) << std::left << solutionVector[i];
	}
}

void poisson1d::printxOnehalf()
{
	myoutput << "value of u_h at x=1/2 : " << solutionVector[(N-2)/2] <<std::endl;

}

void poisson1d::solcomp1plusx()
{
	double* comp = new double[N - 1];
	for (int i = 0; i < N - 1; i++)
	{
		comp[i] = solutionVector[i] - (1 + (i + 1)*h);
	}
	myoutput << std::endl;
	myoutput << "solution - (1+x)" << std::endl;
	for (int i = 0; i < N - 1; i++)
	{
		myoutput << std::setw(15) << std::left << comp[i];
	}
}

void poisson1d::norms()
{
	double diff[N-1];
	for (int i=0;i<N-1;i++)
	{
		diff[i]=std::abs(solutionVector[i]-sin(M_PI*(i+1)*h)/(M_PI*M_PI));
		//std::cout << "diff[" << i << "]=" <<diff[i] <<std::endl;
	}

	//find maximum norm and 2 norm

	for (int i=0;i<N-1;i++)
	{
		if(diff[i]>norminf)
			norminf=diff[i];
		norm2 += pow(diff[i],2);
	}
	norm2 *= h;
	norm2 = sqrt(norm2);
}

double poisson1d::getNorm2()
{
	//std::cout <<"norm2="<<norm2<<'\n';
	return norm2;
}

double poisson1d::getNormInf()
{
	return norminf;
}

void poisson1d::forwardEuler()
{
	double* ulast= new double[N-1];
	double* Au = new double[N-1];
	//Initial Condition
	for (int k=0; k<N-1; k++)
		ulast[k]=4*(k+1)*h*(1-(k+1)*h);

	// myoutput << "Init: " << std::endl;
	// for (int i = 0; i < N - 1; i++)
	// {
	// 	myoutput << std::setw(10) << std::left << ulast[i];
	// }


	for (double time=tau; time<=0.4; time+=tau)
	{
		for (int k=0; k<N-1; k++)
			Au[k]=0;
		
		//Calculate Au (or u_xx)
		for (int k=0; k<N-1; k++)
		{
			for(int i=0; i<3; i++)
			{
				if (sparseStiffnessMatrix[i][k] - 0 < 1e-6) continue;
				switch(i){
					case 0 : Au[k] += tau*pow(h,-2)*(-sparseStiffnessMatrix[i][k])*ulast[k-1]  ; break;
					case 1 : Au[k] += tau*pow(h,-2)*(-sparseStiffnessMatrix[i][k])*ulast[k]  ; break;
					case 2 : Au[k] += tau*pow(h,-2)*(-sparseStiffnessMatrix[i][k])*ulast[k+1]  ; break;
					//Minus sign here because the sparseStiffnessMatrix is -u_(xx)
				}
			}
		}
		myoutput<<"\n"<<"Time = "<<time<<std::endl;
		printsol();

		for (int k=0; k<N-1; k++)
			solutionVector[k] = ulast[k]+ Au[k];

		//Store ulast
		for (int k=0; k<N-1; k++)
			ulast[k]=solutionVector[k];
		
		myoutput<<"\n"<<"Time = "<<time<<std::endl;
		printsol();
	}
}



void poisson1d::backwardEuler()
{
	
}