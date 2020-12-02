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
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, mat_A, N, mat_B, N, 1, mat_D, N);
	unsigned int end_time2 = clock();
	delete[] mat_D;
	cout << "\n\n2-ой вариант перемножения: Результат работы функции cblas_sgemm \nиз библиотеки BLAS (рекомендуемая реализация из Intel MKL)." << endl;
	cout << "Сложность алгоритма: " << 2 * pow(N, 3) << endl;
	cout << "Производительность в MFlops: " << 2 * pow(N, 3) / end_time2 * pow(10, -6);

	// my way
	// Winograd's algorithm	
	unsigned int start_time3 = clock();
	float rowFactor[N];
	float columnFactor[N];
	for (int i = 0; i < N; ++i)
	{
		rowFactor[i] = mat_A[i * N + 1] * mat_A[i * N + 2];
		columnFactor[i] = mat_B[N + i] * mat_B[2 * N + i];
		for (int j = 1; j < N / 2; ++j)
		{
			rowFactor[i] += mat_A[i * N + 2 * j - 1] * mat_A[i * N + 2 * j];
			columnFactor[i] += mat_B[2 * j * N - 1 + i] * mat_B[2 * j * N + i];
		}
	}

	for (int i = 0; i < N; ++i)
	{
		columnFactor[i] = mat_B[N + i] * mat_B[2 * N + i];
		for (int j = 1; j < N / 2; ++j)
			columnFactor[i] += mat_B[2 * j * N - 1 + i] * mat_B[2 * j * N + i];
	}
	for (int i = 0; i < N; ++i)
	{
		for (int j = 1; j < N; ++j)
		{
			mat_E[i * N + j] = -rowFactor[i] - columnFactor[j];
			for (int k = 0; k < N / 2; ++k)
				mat_E[i * N + j] += (mat_A[i * N + 2 * k - 1] + mat_B[2 * k * N + j]) * (mat_A[i * N + 2 * k] + mat_B[2 * k * N + j]);
		}
	}
	unsigned int end_time3 = clock();
	cout << "\n\n3-ий вариант перемножения: Оптимизированный алгоритм по моему выбору, \nнаписанный мной, производительность должна быть не ниже 30% от 2-го варианта." << endl;
	cout << "Сложность алгоритма: " << 2 * pow(N, 3) << endl;
	cout << "Производительность в MFlops: " << 2 * pow(N, 3) / end_time3 * pow(10, -6);


	cout << "\n\n\nФИТУ 2-5, Курочкин Дмитрий Сергеевич";

	return 0;
}
