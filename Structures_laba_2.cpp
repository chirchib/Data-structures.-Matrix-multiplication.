#include <iostream>
#include <ctime>
#include <mkl_cblas.h>
#include "Strassen's Algorithm.h"
using namespace std;

float Random(float min, float max)
{
	return  (float)(rand()) / RAND_MAX * (max - min) + min;
}

float* Linear_algebra(float* mat1, float* mat2, int size)
{
	float* mat = new float[size * size];

	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j < size; ++j)
		{
			mat[i * size + j] = 0;
			for (int k = 0; k < size; ++k)
				mat[i * size + j] += mat1[i * N + j] * mat2[i * size + j];
		}
	}

	return mat;
}

void DisplayMat(float* mat, int size)
{
	cout << "\n\n" ;
	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j < size; ++j)
		{
			cout << mat[size * i + j] << " ";
		}
		cout << endl;
	}
	cout << "\n\n";
}

int main()
{
	setlocale(LC_ALL, "Russian");

	const unsigned int N = 2048;
	float min = 0.0, max = 10.0;
	float *mat_A = new float[N * N];
	float *mat_B = new float[N * N];
	
	float *mat_C = new float[N * N];
	float *mat_D = new float[N * N];
	float *mat_E = new float[N * N];

	srand(time(0));

	// initializing of the matrices
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			mat_A[N*i + j] = Random(min, max);
			mat_B[N*i + j] = Random(min, max);
		}
	}
	
	// linear algebra

	int start_time = clock();
	mat_C = Linear_algebra(mat_A, mat_B, N);
	int end_time = clock();
	delete[] mat_C;

	double time1 = (double)(end_time - start_time) / CLOCKS_PER_SEC;
	cout << "1-ый вариант перемножения : Формула из линейной алгебры." << endl;
	cout << "Время работы алгоритма : " << time1 << endl;
	cout << "Сложность алгоритма : " << 2 * pow(N, 3) << endl;
	cout << "Производительность в MFlops : " << 2 * pow(N, 3) / time1 * pow(10, -6) << endl;

	// func cblas_sgemm from BLAS

	start_time = clock();
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, mat_A, N, mat_B, N, 0.0, mat_D, N);
	end_time = clock();
	delete[] mat_D;

	double time2 = (double)(end_time - start_time) / CLOCKS_PER_SEC;
	cout << "\n\n2-ой вариант перемножения : Результат работы функции cblas_sgemm \nиз библиотеки BLAS (рекомендуемая реализация из Intel MKL)." << endl;
	cout << "Время работы алгоритма : " << time2 << endl;
	cout << "Сложность алгоритма : " << 2 * pow(N, 3) << endl;
	cout << "Производительность в MFlops : " << 2 * pow(N, 3) / time2 * pow(10, -6) << endl;

	// my way
	// Strassen's algorithm

	start_time = clock();
	mat_E  = Strassen(mat_A, mat_B, N);
	end_time = clock();
	delete[] mat_E;

	double time3 = (double)(end_time - start_time) / CLOCKS_PER_SEC;
	cout << "\n\n3-ий вариант перемножения : Оптимизированный алгоритм по моему выбору, \nнаписанный мной, производительность должна быть не ниже 30% от 2-го варианта." << endl;
	cout << "Время работы алгоритма : " << time3 << endl;
	cout << "Сложность алгоритма : " << 2 * pow(N, 3) << endl;
	cout << "Производительность в MFlops : " << 2 * pow(N, 3) / time3 * pow(10, -6) << endl;


	cout << "\n\n\nФИТУ 2-5, Курочкин Дмитрий Сергеевич";

	return 0;
}
