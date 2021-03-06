#include "CRS.h"

void initialize_matrix(const int N, const int NZ, crsMatrix *mtx)
{
	mtx->N = N, mtx->NZ = NZ;
	mtx->els = (double *)malloc(NZ * sizeof(double));
	mtx->columns = (int *)malloc(NZ * sizeof(int));
	mtx->rowIndexs = (int *)malloc((N + 1) * sizeof(int));
}
void free_matrix(crsMatrix *mtx)
{
	free(mtx->els), free(mtx->columns), free(mtx->rowIndexs);
}

void generateRegularCrs(const int N, const int cntInRow, crsMatrix *mtx)
{
	int notNull = cntInRow * N;
	int swap, currCol;
	bool exist;

	initialize_matrix(N, notNull, mtx);
	for (int row = 0; row < N; row++) {
		for (int col = 0; col < cntInRow; col++) {
			do {
				currCol = rand() % N;//�������� ��������� �������
				exist = false;
				for (int i = 0; i < col; i++)
					if (currCol == mtx->columns[row * cntInRow + i]) {
						exist = true; break;
					}
			} while (exist == true);
			mtx->columns[row * cntInRow + col] = currCol;
		}

		for (int col = 0; col < cntInRow - 1; col++)
			for (int val = 0; val < cntInRow - 1; val++)
				if (mtx->columns[row * cntInRow + val] >
					mtx->columns[row * cntInRow + val + 1]) {
					//������ �������
					swap = mtx->columns[row * cntInRow + val];
					mtx->columns[row * cntInRow + val] = mtx->columns[row * cntInRow + val + 1];
					mtx->columns[row * cntInRow + val + 1] = swap;
				}
	}

	for (int val = 0; val < notNull; val++) {
		mtx->els[val] = (double)rand() / RAND_MAX * MAX_VAL;
	}
	for (int row = 0, i = 0; row <= N; row++)
		mtx->rowIndexs[row] = i, i += cntInRow;
}

void generateSpecialCrs(const int N, const int cntInRow, crsMatrix *mtx)
{
	double end = pow((double)cntInRow, 1.0 / 3.0);
	double step = end / N;
	vector<int> *columns = new vector<int>[N];
	int rowNZ, NZ = 0;
	int num1, num2;

	for (int i = 0; i < N; i++) {
		rowNZ = int(pow((double(i + 1) * step), 3) + 1);
		NZ += rowNZ;
		num1 = (rowNZ - 1) / 2;
		num2 = rowNZ - 1 - num1;
		if (rowNZ != 0) {
			if (i < num1) {
				num2 += num1 - i;
				num1 = i;
				for (int j = 0; j < i; j++) columns[i].push_back(j);
				columns[i].push_back(i);
				for (int j = 0; j < num2; j++) columns[i].push_back(i + 1 + j);
			}
			else {
				if (N - i - 1 < num2) {
					num1 += num2 - (N - 1 - i);
					num2 = N - i - 1;
				}
				for (int j = 0; j < num1; j++) columns[i].push_back(i - num1 + j);
				columns[i].push_back(i);
				for (int j = 0; j < num2; j++) columns[i].push_back(i + j + 1);
			}
		}
	}
	initialize_matrix(N, NZ, mtx);

	int count = 0;
	int sum = 0;

	for (int i = 0; i < N; i++) {
		mtx->rowIndexs[i] = sum;
		sum += columns[i].size();
		for (unsigned int j = 0; j < columns[i].size(); j++) {
			mtx->columns[count] = columns[i][j];
			mtx->els[count] = (double)rand() / RAND_MAX * MAX_VAL;
			count++;
		}
	}
	mtx->rowIndexs[N] = sum;
	delete[] columns;
}

int symbolicMult(const int N, const crsMatrix *A, const crsMatrix *B, crsMatrix *C)
{
	int n = A->N, nz = 0, indexA, k, kEnd;
	vector<int> rowIndexs, columns;
	rowIndexs.push_back(0);

	int *temp = new int[N];

	for (int i = 0; i < n; i++) {
		memset(temp, EMPTY, N * sizeof(int));
		for (int j = A->rowIndexs[i]; j < A->rowIndexs[i + 1]; j++)
			temp[A->columns[j]] = j;
		for (int j = 0; j < N; j++) {
			indexA = EMPTY, k = B->rowIndexs[j], kEnd = B->rowIndexs[j + 1];
			while (indexA == EMPTY && k < kEnd) {
				indexA = temp[B->columns[k]];
				k++;
			}
			if (indexA != EMPTY) { columns.push_back(j); nz++; }
		}
		rowIndexs.push_back(nz);
	}

	initialize_matrix(n, nz, C);

	for (int j = 0; j < nz; ++j) C->columns[j] = columns[j];
	for (int i = 0; i <= n; ++i) C->rowIndexs[i] = rowIndexs[i];

	delete[] temp;
	return 0;
}
int numericMult(const int N, const crsMatrix *A, const crsMatrix *B, crsMatrix *C)
{
	int n = A->N, index = 0, indexA, startC, finishC, colC;

	int *temp = new int[N];
	double sum = 0, iSum = 0;

	for (int i = 0; i < n; ++i) {
		startC = C->rowIndexs[i], finishC = C->rowIndexs[i + 1];
		if (finishC > startC) {
			memset(temp, EMPTY, sizeof(int)* N);
			for (int j = A->rowIndexs[i]; j < A->rowIndexs[i + 1]; ++j) temp[A->columns[j]] = j;
			for (int j = startC; j < finishC; ++j, sum = 0, iSum = 0) {
				colC = C->columns[j];
				for (int k = B->rowIndexs[colC]; k < B->rowIndexs[colC + 1]; ++k) {
					indexA = temp[B->columns[k]];
					if (indexA != EMPTY) {
						sum += A->els[indexA] * B->els[k];
					}
				}
				C->els[index] = sum;
				index++;
			}
		}
	}

	delete[] temp;
	return 0;
}

void matr_transpose(const crsMatrix *A, crsMatrix *AT)
{
	int tmp, sum = 0;
	int index, rowIndex;

	memset(AT->rowIndexs, 0, (AT->N + 1) * sizeof(int));
	for (int i = 0; i < A->NZ; i++)
		AT->rowIndexs[A->columns[i] + 1]++;
	for (int i = 1; i <= A->N; i++) {
		tmp = AT->rowIndexs[i];
		AT->rowIndexs[i] = sum;
		sum += tmp;
	}
	for (int i = 0; i < A->N; i++)
		for (int j = A->rowIndexs[i]; j < A->rowIndexs[i + 1]; j++) {
			rowIndex = A->columns[j];
			index = AT->rowIndexs[rowIndex + 1];
			AT->els[index] = A->els[j];
			AT->columns[index] = i;
			AT->rowIndexs[rowIndex + 1]++;
		}
}
