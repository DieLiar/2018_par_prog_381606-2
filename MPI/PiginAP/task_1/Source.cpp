#include <iostream>
#include <stdio.h>
#include <conio.h>
#include <mpi.h>
using namespace std;
void MatrixGenerate(int **&Matrix, int Size)
{
	for (int i = 0; i<Size; i++)
	{
		for (int j = 0; j<Size; j++)
		{
			Matrix[i][j] = rand() % 90 + 10;
		}
	}
}
void VectorSerialMin(int *&Vector, int *&MinVector, int Size)
{
	for (int i = 0; i<Size; i++)
	{
		MinVector[i] = INT_MAX;
		for (int j = 0; j<Size; j++)
		{
			if (Vector[i + j*Size] < MinVector[i])
				MinVector[i] = Vector[i + j*Size];
		}
	}
/*
	for (int i = 0; i < Size; i++)
	{
		cout << MinVector[i] << " " << endl;		//Min Vector COUT
	}
*/

}

void MatrixToVector(int **&Matrix, int *&Vector, int Size)
{
	for (int i = 0; i<Size; i++)
	{
		for (int j = 0; j<Size; j++)
		{
			Vector[i*Size + j] = Matrix[i][j];
		}
	}
}

int main(int argc, char* argv[])
{
	double starttime = 0;
	double endtime = 0;
	int TaskLength = 0;
	int Size = 4, Num, ID;
	int min = INT_MAX;
	int** Matrix;
	int *MatrixVector, *SerialMinVector, *ParallelMinVector;
	Matrix = new int*[Size];
	MatrixVector = new int[Size*Size];
	SerialMinVector = new int[Size];
	ParallelMinVector = new int[Size];
	for (int i = 0; i<Size; i++)
		Matrix[i] = new int[Size];

	MPI_Init(&argc, &argv);		// starting mpi
	MPI_Comm_size(MPI_COMM_WORLD, &Num);
	MPI_Comm_rank(MPI_COMM_WORLD, &ID);
	TaskLength = Size / Num;
	if (Size%Num != 0)
		TaskLength++;
	if (ID == 0)
	{
		MatrixGenerate(Matrix, Size);
		MatrixToVector(Matrix, MatrixVector, Size);
		starttime = MPI_Wtime();	             

			


                 //START TIME
	}
	MPI_Bcast(MatrixVector, Size*Size, MPI_INT, 0, MPI_COMM_WORLD);
	for (int i = TaskLength*ID; i<TaskLength*(ID + 1); i++)
	{
		min = INT_MAX;
		for (int j = 0; j<Size; j++)
		{
			if (MatrixVector[i + j*Size] < min)
			{
				min = MatrixVector[i + j*Size];
			}
			cout << " min = " << min << endl;
		}
		if (ID == 0)
		{
			ParallelMinVector[i] = min;
		}
		else
		{
			MPI_Send(&min, 1, MPI_INT, 0, 15, MPI_COMM_WORLD);
		}
		
	}
	if (ID == 0)
	{

		for (int i = 1; i<Num; i++)
		{
			for (int j = 0; j<TaskLength; j++)
			{
				MPI_Status status;
				MPI_Recv(&min, 1, MPI_INT, i, 15, MPI_COMM_WORLD, &status);
				ParallelMinVector[i*TaskLength + j] = min;
			}
		}
		endtime = MPI_Wtime();                                 //END TIME


		cout << "====================== MATRIX ======================" << endl;
	for (int i = 0; i<Size; i++)
	{
		for (int j = 0; j<Size; j++)
		{

			cout << Matrix[i][j]<< " ";
		}
	}
	cout << "  " << endl;


	cout << "======================" << endl;
	

	
		cout << "Time for parallel : " << endtime - starttime << endl;
		double s = MPI_Wtime();
		VectorSerialMin(MatrixVector, SerialMinVector, Size);
		double e = MPI_Wtime();
		cout << "Time for Serial : " << e - s;
		for (int i = 0; i<Size; i++)
			if (ParallelMinVector[i] != SerialMinVector[i])
			{
				cout << endl << i << " Wrong!!!" << endl;
				break;
			}
	}

	MPI_Finalize();    //end mpi


	return 0;
}