#include "poisson2d.h"
#include <iostream>
#include <iomanip>
#include <math.h>

poisson2d::poisson2d(int _N)
{
	N = _N;
	h = 1.0 / N;
	stiffnessMatrix = new double*[(N - 1) * (N - 1)];
	for (int i = 0; i < (N - 1)*(N-1); i++)
		stiffnessMatrix[i] = new double[(N - 1) * (N-1)];
	solutionVector = new double[(N - 1) * (N - 1)];
	rhs = new double[(N - 1) * (N - 1)];
	myoutput.open("poisson2d_output.txt",std::ios_base::out);

	sparseStiffnessMatrix = new double*[5];
	for (int k=0; k < 5; k++)
		sparseStiffnessMatrix[k]= new double [(N - 1)*(N - 1)];
}


poisson2d::~poisson2d()
{

}

void poisson2d::assembleStiffnessMatrix()
{
	//initialization
	for (int i = 0; i < (N - 1)*(N - 1); i++)
		for (int j = 0; j < (N - 1)*(N - 1); j++)
		{
			stiffnessMatrix[i][j] = 0.0;
		}

	for (int k = 0; k < N - 1; k++)
	{
		for (int i = 0; i < N-1; i++)
		{
			stiffnessMatrix[3 * k + i][3 * k + i] = 4;
			if (i-1>=0)
				stiffnessMatrix[3 * k + i][3 * k + i - 1] = -1;
			if (i+1<=N-2)
				stiffnessMatrix[3 * k + i][3 * k + i + 1] = -1;
		}
		if (k - 1 >= 0)
		{
			for (int j = 0; j < N - 1; j++)
			{
				stiffnessMatrix[3 * k + j][3 * (k - 1) + j] = -1;
			}
		}
		if (k + 1 <= N - 2)
		{
			for (int j = 0; j < N - 1; j++)
			{
				stiffnessMatrix[3 * k + j][3 * (k + 1) + j] = -1;
			}
		}
	}
}
void poisson2d::assembleSparseStiffnessMatrix()
{
	//initialization
	for (int i = 0; i < 5; i++)
		for (int j = 0; j < (N - 1)*(N - 1); j++)
			sparseStiffnessMatrix[i][j] = 0.0;

	for (int k=N-1; k<(N-1)*(N-1); k++)
		sparseStiffnessMatrix[0][k]=-1.0;

	for (int k=0; k<(N-1)*(N-1); k++)
	{
		if (k%(N-1)==0)
			continue;
		sparseStiffnessMatrix[1][k]=-1.0;
	}

	for (int k=0; k<(N-1)*(N-1); k++)
		sparseStiffnessMatrix[2][k]=4.0;

	for (int k=0; k<(N-1)*(N-1); k++)
	{
		if (k%(N-1)==N-2)
			continue;
		sparseStiffnessMatrix[3][k]=-1.0;
	}

	for (int k=0; k<(N-1)*(N-1)-(N-1); k++)
		sparseStiffnessMatrix[4][k]=-1.0;
		
}

void poisson2d::assemblerhs()
{
	//initialize
	for (int k=0; k<(N-1)*(N-1); k++)
		rhs[k]=0;
	//Boundary Condition x+y
	for (int k=0; k<N-1; k++)
		rhs[k]+=((k+1)*h);
	for (int k=(N-2)*(N-1); k<(N-1)*(N-1); k++)
		rhs[k]+=((k+1)*h+1);
	for (int k=0; k<(N-1)*(N-1); k+=N-1)
		rhs[k]+=(((k/(N-1))+1)*h);
	for (int k=N-2; k<(N-1)*(N-1); k+=N-1)
		rhs[k]+=(((k/(N-1))+1)*h+1);
}

void poisson2d::Jacobi()
{
	double* ulast = new double[(N-1)*(N-1)];
	double* r = new double[(N-1)*(N-1)];

	double rnorm=0, rnormlast=0;
	double rho=0;
	myoutput<<"rnorm"<<std::setw(10)<<"rho"<<std::endl;
	for (int j=0; j<(N-1)*(N-1); j++)
	{
		//initialize
		ulast[j]=0;
	}
	for (int count=0; count<=100; count++)
	{
		// myoutput << "loop: " <<count <<std::endl;
		// myoutput << "\n" << "ulast: " << std::endl;
		// for (int i = 0; i < (N - 1)*(N - 1); i++)
		// {
		// 	myoutput << std::setw(10) << std::left << ulast[i];
		// }		
		// myoutput<<std::endl;

		for (int j=0; j<(N-1)*(N-1); j++)
		{
			r[j]=rhs[j];
			for (int k=0; k<5; k++)
			{
				if(abs(sparseStiffnessMatrix[k][j]-0) < 1e-6) continue;
				if (j==4)
					myoutput<<k;
				switch(k)
				{				
					case 0: r[j] -= pow(h,-2)*(sparseStiffnessMatrix[k][j])*ulast[j-(N-1)]; 		break;
					case 1: r[j] -= pow(h,-2)*(sparseStiffnessMatrix[k][j])*ulast[j-1]; 			break;
					case 2: r[j] -= pow(h,-2)*(sparseStiffnessMatrix[k][j])*ulast[j];  			    break;
					case 3: r[j] -= pow(h,-2)*(sparseStiffnessMatrix[k][j])*ulast[j+1]; 			break;
					case 4: r[j] -= pow(h,-2)*(sparseStiffnessMatrix[k][j])*ulast[j+(N+1)]; 		break;
				}
			}
		}
		//calculate r norm
		for (int k=0; k<(N-1)*(N-1); k++)
			rnorm += pow(r[k],2);
		rnorm=sqrt(rnorm);

		if(count!=0)
			rho=rnorm/rnormlast;
		

		myoutput << std::setw(5) << std::left << rnorm;
		myoutput << std::setw(5) << std::left << rho;
		myoutput << std::endl;
		
		
		rnormlast=rnorm;
		rnorm=0;

		myoutput << "\n" << "r: " << std::endl;
		for (int i = 0; i < (N - 1)*(N - 1); i++)
		{
			myoutput << std::setw(10) << std::left << r[i];
		}		
		printsol();

		//update U^(k+1)		
		for (int k=0; k<(N-1)*(N-1); k++)
		{
			solutionVector[k]=ulast[k]+r[k]/(4*pow(h,-2));
			ulast[k]=solutionVector[k];
		}
		printsol();

	}
}

void poisson2d::printStiffnessMatrix()
{
	myoutput << "Stiffness Matrix: " << std::endl;
	for (int i = 0; i < (N - 1)*(N - 1); i++)
	{
		for (int j = 0; j < (N - 1)*(N - 1); j++)
		{
			myoutput << std::setw(5) << std::left << stiffnessMatrix[i][j];
		}
		myoutput << std::endl;
	}
}

void poisson2d::printSparseStiffnessMatrix()
{
	myoutput << "Sparse Stiffness Matrix:" <<std::endl;
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < (N - 1)*(N - 1); j++)
		{
			myoutput << std::setw(5) << std::left << sparseStiffnessMatrix[i][j];
		}
		myoutput << std::endl;
	}
}

void poisson2d::printrhs()
{
	myoutput << std::endl;
	myoutput << "RHS: " << std::endl;
	for (int i = 0; i < (N - 1)*(N-1); i++)
	{
		myoutput << std::setw(5) << std::left << rhs[i];
	}
	myoutput << std::endl;
}

void poisson2d::calculateDefect()
{
	myoutput << std::endl;
	myoutput << "Defect: " << std::endl;
	double *defect = new double[(N - 1)*(N - 1)];
	for (int k = 0; k < (N - 1)*(N - 1); k++)
	{
		defect[k] = -rhs[k];
		for (int m = 0; m < (N - 1)*(N - 1); m++)
		{
			defect[k] += stiffnessMatrix[k][m];
		}
	}
	
	for (int i = 0; i < (N - 1)*(N - 1); i++)
	{
		myoutput << std::setw(5) << std::left << defect[i];
	}
	myoutput << std::endl;

}
void poisson2d::printsol()
{
	
	myoutput << "\n" << "Sol: " << std::endl;
	for (int i = 0; i < (N - 1)*(N - 1); i++)
	{
		myoutput << std::setw(10) << std::left << solutionVector[i];
	}
}

// void poisson2d::printr()
// {
	
// 	myoutput << "\n" << "r: " << std::endl;
// 	for (int i = 0; i < (N - 1)*(N - 1); i++)
// 	{
// 		myoutput << std::setw(10) << std::left << r[i];
// 	}
// }