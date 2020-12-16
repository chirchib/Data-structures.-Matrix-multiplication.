#include <iostream>
#include <ctime>
#include <mkl_cblas.h>
#include "Strassen's Algorithm.h"
using namespace std;

float Random(float min, float max)
{
	return  (float)(rand()) / RAND_MAX * (max - min) + min;
}

float** Trans(float** mat, int size)
{
	float** matTrans = new float*[size];
	for (int i = 0; i < size; ++i)
		matTrans[i] = new float[size];

	for (int i = 0; i < size; ++i)
		for (int j = 0; j < size; ++j)
			matTrans[i][j] = mat[j][i];

	return matTrans;
}

float dot_prod(float* a, float* b, int n)
{
	float s1 = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0;
	int n8 = (n / 8) * 8;
	for (int i = 0; i < n8; i += 8) 
	{
		s1 += a[i] * b[i] + a[i + 1] * b[i + 1];
		s2 += a[i + 2] * b[i + 2] + a[i + 3] * b[i + 3];
		s3 += a[i + 4] * b[i + 4] + a[i + 5] * b[i + 5];
		s4 += a[i + 6] * b[i + 6] + a[i + 7] * b[i + 7];
	}
	for (int i = n8; i < n; i++) 
		s4 += a[i] * b[i];
	
	return (s4 + s3) + (s2 + s1);
}

void DisplayMat(float** mat, int size)
{
	cout << "\n" ;
	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j < size; ++j)
		{
			cout << mat[i][j] << " ";
		}
		cout << endl;
	}
	cout << "\n";
}

int main()
{
	setlocale(LC_ALL, "Russian");

	const int N = 2;
	float min = 0.0, max = 1.0;

	float **A = new float*[N];
	float **B = new float*[N];

	float **C = new float*[N];
	float **D = new float*[N];
	float **E = new float*[N];
	for (int i = 0; i < N; ++i)
	{
		A[i] = new float[N];
		B[i] = new float[N];
		
		C[i] = new float[N];
		D[i] = new float[N];
		E[i] = new float[N];
	}

	srand(time(0));

	// initializing of the matrices
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			A[i][j] = Random(min, max);
			B[i][j] = Random(min, max);
		}
	}
	
	// linear algebra
	//////////////////////////////////////////////////////

	int start_time = clock();

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			C[i][j] = 0;
			for (int k = 0; k < N; ++k)
				C[i][j] += A[i][k] * B[k][j];
		}
	}

	int end_time = clock();

	for (int i = 0; i < N; ++i)
	{
		delete[] C[i];
	}
	delete[] C;

	double time1 = (double)(end_time - start_time) / CLOCKS_PER_SEC;
	cout << "1-ый вариант перемножения : Формула из линейной алгебры." << endl;
	cout << "Время работы алгоритма : " << time1 << endl;
	cout << "Сложность алгоритма : " << 2 * pow(N, 3) << endl;
	cout << "Производительность в MFlops : " << 2 * pow(N, 3) / time1 * pow(10, -6) << endl;

	// func cblas_sgemm from BLAS
	///////////////////////////////////////////////////////

	float *mat1 = new float[N * N];
	float *mat2 = new float[N * N];
	float *mat3 = new float[N * N];
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			mat1[i * N + j] = A[i][j];
			mat2[i * N + j] = B[i][j];
		}
	}

	start_time = clock();

	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, mat1, N, mat2, N, 0.0, mat3, N);

	end_time = clock();

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			D[i][j] = mat3[i * N + j];
		}
	}

	for (int i = 0; i < N; ++i)
	{
		delete[] D[i];
	}
	delete[] D;

	double time2 = (double)(end_time - start_time) / CLOCKS_PER_SEC;
	cout << "\n\n2-ой вариант перемножения : Результат работы функции cblas_sgemm \nиз библиотеки BLAS (рекомендуемая реализация из Intel MKL)." << endl;
	cout << "Время работы алгоритма : " << time2 << endl;
	cout << "Сложность алгоритма : " << 2 * pow(N, 3) << endl;
	cout << "Производительность в MFlops : " << 2 * pow(N, 3) / time2 * pow(10, -6) << endl;

	// my way
	////////////////////////////////////////////////////////////////////////

	start_time = clock();

	float** TransB = Trans(B, N);
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
			E[i][j] = dot_prod(A[i], TransB[j], N);

	end_time = clock();


	for (int i = 0; i < N; ++i)
	{
		delete[] E[i];
	}
	delete[] E;

	double time3 = (double)(end_time - start_time) / CLOCKS_PER_SEC;
	cout << "\n\n3-ий вариант перемножения : Оптимизированный алгоритм по моему выбору. Алгоритм Штрассена." << endl;
	cout << "Время работы алгоритма : " << time3 << endl;
	cout << "Сложность алгоритма : " << 2 * pow(N, 3) << endl;
	cout << "Производительность в MFlops : " << 2 * pow(N, 3) / time3 * pow(10, -6) << endl;


	
	
	cout << "\n\nФИТУ 2-5, Курочкин Дмитрий Сергеевич" << endl;
	return 0;
}
