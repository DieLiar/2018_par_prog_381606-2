//Лабораторная работа 2
//7. Передача от всех одному 

#include <iostream>
#include <time.h>
#include <malloc.h>
#include <math.h>
#include <cassert>
#include<cstdlib>
#include<mpi.h>


int NEW_TREE_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm);
size_t stype(MPI_Datatype type);
void arrayinit(int *array, int n, int min, int max);
int NEW_Reduce(const void *a, void *a1, int n, MPI_Datatype t, MPI_Op o, int r, MPI_Comm c);

int main(int argc, char **argv)
{
	int ProcNum, ProcRank;
	int *array;
	int n, root, k;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	n = atoi(argv[1]);
	array = new int[n];
	arrayinit(array, n, atoi(argv[2]), atoi(argv[3]));
	root = atoi(argv[4]);
	int * res;
	if (ProcRank == root)
	{
		res = new int[n];
		/*for (int i = 0; i < n; i++)
		{
			std::cout << array[i] << " ";
		}*/
		std::cout << std::endl;
	}
	k = atoi(argv[5]);
	double tstart1;

	if (ProcRank == 0)
	{
		tstart1 = MPI_Wtime();
	}
	switch (k)
	{
	case 0:
		res = new int[n];
		MPI_Reduce(array, res, n, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		break;
	case 1:
		res = new int[n];
		NEW_TREE_Reduce(array, res, n, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		break;
	case 2:
		res = new int[n];
		NEW_Reduce(array, res, n, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
		break;
	default:
		MPI_Finalize();
		return 0;
	}

	if (ProcRank == 0)
	{
		std::cout << std::endl;
		if (k == 0)
			std::cout << "MPI_Reduce: " << std::endl;
		else if (k == 1)
			std::cout << "New_Reduce: " << std::endl;
		else if (k == 2)
			std::cout << "TREE_Reduce_Binary_Tree: " << std::endl;
		std::cout << "Elapsed time = " << MPI_Wtime() - tstart1 << std::endl;
		//std::cout << "You are in process " << 0 << std::endl;
		/*for (int i = 0; i < n; i++)
			std::cout << res[i] << " ";*/
	}


	MPI_Finalize();
	return 0;

}

void arrayinit(int *array, int n, int min, int max)
{

	for (int i = 0; i < n; i++)
	{
		array[i] = min + rand() % (max - min);
	}
}


size_t stype(MPI_Datatype type)
{
	size_t result;
	switch (type)
	{
	default:
		std::cout << "Incorrect datatype" << std::endl;
		exit(1);
	case MPI_INT:
		result = sizeof(int);
		break;
	case MPI_FLOAT:
		result = sizeof(float);
		break;
	case MPI_DOUBLE:
		result = sizeof(double);
		break;
	}
	return result;
}

void max(void * buf1, void * buf2, void * res, size_t size, MPI_Datatype type)
{
	switch (type)
	{
	case MPI_INT:
	{
		int * buffer = static_cast<int*>(res);
		int * buffer1 = static_cast<int*>(buf1);
		int * buffer2 = static_cast<int*>(buf2);
		for (size_t i = 0; i < size; i++)
		{
			buffer[i] = buffer1[i] > buffer2[i];
		}
		break;
	}
	case MPI_FLOAT:
	{
		float * buffer = static_cast<float*>(res);
		float * buffer1 = static_cast<float*>(buf1);
		float * buffer2 = static_cast<float*>(buf2);
		for (size_t i = 0; i < size; i++)
		{
			buffer[i] = buffer1[i] > buffer2[i];
		}
		break;
	}
	case MPI_DOUBLE:
	{
		double * buffer = static_cast<double*>(res);
		double * buffer1 = static_cast<double*>(buf1);
		double * buffer2 = static_cast<double*>(buf2);
		for (size_t i = 0; i < size; i++)
		{
			buffer[i] = buffer1[i] > buffer2[i];
		}
		break;
	}
	}
}

void min(void * buf1, void * buf2, void * res, size_t size, MPI_Datatype type)
{
	switch (type)
	{
	case MPI_INT:
	{
		int * buffer = static_cast<int*>(res);
		int * buffer1 = static_cast<int*>(buf1);
		int * buffer2 = static_cast<int*>(buf2);
		for (size_t i = 0; i < size; i++)
		{
			buffer[i] = buffer1[i] < buffer2[i];
		}
		break;
	}
	case MPI_FLOAT:
	{
		float * buffer = static_cast<float*>(res);;
		float * buffer1 = static_cast<float*>(buf1);
		float * buffer2 = static_cast<float*>(buf2);
		for (size_t i = 0; i < size; i++)
		{
			buffer[i] = buffer1[i] < buffer2[i];
		}
		break;
	}
	case MPI_DOUBLE:
	{
		double * buffer = static_cast<double*>(res);;
		double * buffer1 = static_cast<double*>(buf1);
		double * buffer2 = static_cast<double*>(buf2);
		for (size_t i = 0; i < size; i++)
		{
			buffer[i] = buffer1[i] < buffer2[i];
		}
		break;
	}
	}
}

void sum(void * buf1, void * buf2, void* res, size_t size, MPI_Datatype type)
{
	switch (type)
	{
	case MPI_INT:
	{
		int * buffer = static_cast<int*>(res);
		int * buffer1 = static_cast<int*>(buf1);
		int * buffer2 = static_cast<int*>(buf2);
		for (size_t i = 0; i < size; i++)
		{
			buffer[i] = buffer1[i] + buffer2[i];
		}
		break;
	}
	case MPI_FLOAT:
	{
		float * buffer = static_cast<float*>(res);
		float * buffer1 = static_cast<float*>(buf1);
		float * buffer2 = static_cast<float*>(buf2);
		for (size_t i = 0; i < size; i++)
		{
			buffer[i] = buffer1[i] + buffer2[i];
		}
		break;
	}
	case MPI_DOUBLE:
	{
		double * buffer = static_cast<double*>(res);
		double * buffer1 = static_cast<double*>(buf1);
		double * buffer2 = static_cast<double*>(buf2);
		for (size_t i = 0; i < size; i++)
		{
			buffer[i] = buffer1[i] + buffer2[i];
		}
		break;
	}
	}
}

void Operation(MPI_Datatype type, MPI_Op op, void* buf1, void* buf2, void* res, size_t size)
{
	switch (op)
	{
	case MPI_MAX:
		max(buf1, buf2, res, size, type);
		break;
	case MPI_MIN:
		min(buf1, buf2, res, size, type);
		break;
	case MPI_SUM:
		sum(buf1, buf2, res, size, type);
		break;
	default:
		break;
	}
}
void NEW_MIN(const void *a, void *a1, int n, MPI_Datatype t)
{
	int i;
	if (t == MPI_INT)
	{
		for (i = 0; i<n; i++)
		{
			if (((int *)a1)[i]>((int *)a)[i])
			{
				((int *)a1)[i] = ((int *)a)[i];
			}
		}
	}
	if (t == MPI_FLOAT)
	{
		for (i = 0; i<n; i++)
		{
			if (((float *)a1)[i]>((float *)a)[i])
			{
				((float *)a1)[i] = ((float *)a)[i];
			}
		}
	}
	if (t == MPI_DOUBLE)
	{
		for (i = 0; i<n; i++)
		{
			if (((double *)a1)[i]>((double *)a)[i])
			{
				((double *)a1)[i] = ((double *)a)[i];
			}
		}
	}
}

//операция суммы заданного типа, результат в "a1"
void NEW_SUM(const void *a, void *a1, int n, MPI_Datatype t)
{
	int i;
	if (t == MPI_INT)
	{
		for (i = 0; i<n; i++)
		{
			((int *)a1)[i] = ((int *)a)[i] + ((int *)a1)[i];
		}
	}
	if (t == MPI_FLOAT)
	{
		for (i = 0; i<n; i++)
		{
			((float *)a1)[i] = ((float *)a)[i] + ((float *)a1)[i];
		}
	}
	if (t == MPI_DOUBLE)
	{
		for (i = 0; i<n; i++)
		{
			((double *)a1)[i] = ((double *)a)[i] + ((double *)a1)[i];
		}
	}
}
void NEW_MAX(const void *a, void *a1, int n, MPI_Datatype t)
{
	int i;
	if (t == MPI_INT)
	{
		for (i = 0; i<n; i++)
		{
			if (((int *)a1)[i]<((int *)a)[i])
			{
				((int *)a1)[i] = ((int *)a)[i];
			}
		}
	}
	if (t == MPI_FLOAT)
	{
		for (i = 0; i<n; i++)
		{
			if (((float *)a1)[i]<((float *)a)[i])
			{
				((float *)a1)[i] = ((float *)a)[i];
			}
		}
	}
	if (t == MPI_DOUBLE)
	{
		for (i = 0; i<n; i++)
		{
			if (((double *)a1)[i]<((double *)a)[i])
			{
				((double *)a1)[i] = ((double *)a)[i];
			}
		}
	}
}
int NEW_Reduce( void *a, void *a1, int n, MPI_Datatype t, MPI_Op o, int r, MPI_Comm c)
{
	int rank1, size1;
	void *a2;
	int f = 0;
	int i;
	a2 = (void*)malloc((n) * sizeof(void*));
	MPI_Status Status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank1);
	MPI_Comm_size(MPI_COMM_WORLD, &size1);
	if (rank1 == r)
	{
		for (i = 0; i<size1; i++)
		{
			if (i != r)
			{
				if (f == 0)
				{
					MPI_Recv(a2, n, t, i, 0, c, &Status);
					f = 1;
				}
				else
				{
					MPI_Recv(a1, n, t, i, 0, c, &Status);
					if (o == MPI_MAX)
					{
						NEW_MAX(a1, a2, n, t);
					}
					if (o == MPI_MIN)
					{
						NEW_MIN(a1, a2, n, t);
					}
					if (o == MPI_SUM)
					{
						NEW_SUM(a1, a2, n, t);
					}
			

				}
			}
		}
		free(a2);
		return 0;
	}
	else
	{
		MPI_Send(a, n, t, r, 0, c);
	}
}


int NEW_TREE_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm)
{
	MPI_Status status;
	int ProcRank, ProcNum;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	int RankRecv, RankSend;
	void * res = malloc(count * stype(type));
	if (ProcRank != root)
		recvbuf = malloc(count * stype(type));
	memcpy(recvbuf, sendbuf, count*stype(type));
	int newProcRank = (ProcRank - root + ProcNum) % ProcNum;
	int mask = 1;
	while (mask < ProcNum)
	{
		if ((newProcRank & mask) == 0)
		{
			RankSend = newProcRank | mask;
			if (RankSend < ProcNum)
			{
				RankSend = (RankSend + root) % ProcNum;
				MPI_Recv(res, count, type, RankSend, 0, comm, &status);
				Operation(type, op, recvbuf, res, recvbuf, count);
			}
		}
		else
		{
			RankRecv = newProcRank&(~mask);
			RankRecv = (RankRecv + root) % ProcNum;
			MPI_Send(recvbuf, count, type, RankRecv, 0, comm);
			break;
		}
		mask = mask << 1;
	}
	if (ProcRank != root)
	{
		free(recvbuf);
	}

	if (ProcRank == root)
	{
		free(res);
	}
	return 0;
}