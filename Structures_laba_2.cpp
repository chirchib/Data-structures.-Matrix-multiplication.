#include <iostream>
#include <ctime>
#include <mkl_cblas.h>
using namespace std;

float Random(float min, float max)
{
	return  (float)(rand()) / RAND_MAX * (max - min) + min;
}

int main()
{
	setlocale(LC_ALL, "Russian");

	const unsigned int N = 2048;
	float min = 0.0, max = 10.0;
	float *mat_A = new float[N * N];
	float *mat_B = new float[N * N];
	
	float *mat_C = new float[N * N];

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
	unsigned int start_time1 = clock();
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			mat_C[i * N + j] = 0;
			for (int k = 0; k < N; ++k)
				mat_C[i * N + j] += mat_A[i * N + j] * mat_B[i * N + j];
		}
	}
	unsigned int end_time1 = clock();
	delete[] mat_C;
	cout << "1-ый вариант перемножения: Формула из линейной алгебры." << endl;
	cout << "Сложность алгоритма: " << 2 * pow(N, 3) << endl;
	cout << "Производительность в MFlops: " << 2 * pow(N, 3) / end_time1 * pow(10, -6);

	// func cblas_sgemm from BLAS
	unsigned int start_time2 = clock();
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, mat_A, N, mat_B, N, 1, mat_C, N);
	unsigned int end_time2 = clock();
	delete[] mat_C;
	cout << "\n\n2-ой вариант перемножения: Результат работы функции cblas_sgemm \nиз библиотеки BLAS (рекомендуемая реализация из Intel MKL)." << endl;
	cout << "Сложность алгоритма: " << 2 * pow(N, 3) << endl;
	cout << "Производительность в MFlops: " << 2 * pow(N, 3) / end_time2 * pow(10, -6);

	// my way
	// Winograd's algorithm
	unsigned int start_time3 = clock();
	
	float rowFactor[N];
	float columnFactor[N];

	// вычисление rowFactors для MATRIX_1
	for (int i = 1; i < N; ++i)
	{
		rowFactor[N * i] = mat_A[N*i + 1] * mat_A[N * i + 2];
		for (int j = 2; j < N / 2; ++j)
			rowFactor[N * i] = rowFactor[N * i] + mat_A[N * i + 2 * j - 1] * mat_A[N * i + 2 * j];

	}

	//////////////////////////////////// ТОЧКА ОСТАНОВЕКИ
	// вычисление	columnFactors для MATRIX_2
	for (int i = 1; i < N; ++i)
	{
		columnFactor[i] = mat_B[1, i] * mat_B[2, i];
		for (int j = 2; j < N / 2; ++j)
			columnFactor[i] = columnFactor[i] + mat_B[2 * j - 1, i] * mat_B[2 * j, i];

	}


	// вычисление матрицы MATRIX_3
	for (int i = 1; i < N; ++i)
	{
		for (int j = 1; i < N; ++j)
		{
			mat_C[i, j] = -rowFactor[i] - columnFactor[j];
				for (int k = 1; k < N / 2; ++k)
					mat_C[i, j] = mat_C[i, j] + (mat_A[i, 2 * k - 1] + mat_B[2 * k, j]) * (mat_A[i, 2 * k] + mat_B[2 * k - 1, j]);
		}
	}

	unsigned int end_time3 = clock();
	cout << "\n\n3-ий вариант перемножения: Оптимизированный алгоритм по вашему выбору, \nнаписанный мной, производительность должна быть не ниже 30% от 2-го варианта." << endl;
	cout << "Сложность алгоритма: " << 2 * pow(N, 3) << endl;
	cout << "Производительность в MFlops: " << 2 * pow(N, 3) / end_time3 * pow(10, -6);


	cout << "\n\n\nФИТУ 2-5, Курочкин Дмитрий Сергеевич";

	return 0;
}
