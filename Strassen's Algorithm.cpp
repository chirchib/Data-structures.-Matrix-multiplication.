#include "Strassen's Algorithm.h"
#include <iostream>
#include <ctime>

float* Strassen(float* matA, float* matB, int size)
{
	float** mat1 = create_mat(size);
	float** mat2 = create_mat(size);

	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j < size; ++j)
		{
			mat1[i][j] = matA[i * size + j];
			mat2[i][j] = matB[i * size + j];
		}
	}

	float** mat_res = mult_mat(mat1, mat2, size);

	float* matE = new float[size * size];
	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j < size; ++j)
		{
			matE[i * size + j] = mat_res[i][j];
		}
	}

	return matE;
}

float** create_mat(int size)
{
	float** mat = new float* [size];
	for (int i = 0; i < size; ++i)
		mat[i] = new float[size];
	
	for (int i = 0; i < size; ++i)
		for (int j = 0; j < size; ++j)
			mat[i][j] = 0;

	return mat;
}

float** mult(float** mat1, float** mat2, int size)
{
	float** mat = create_mat(size);
	for (int i = 0; i < size; ++i)
		for (int j = 0; j < size; ++j)
			for (int k = 0; k < size; ++k)
				mat[i][j] += mat1[i][k] * mat2[k][j];

	return mat;
}

void split_4mat(float** mat, float** mat11, float** mat12, float** mat21, float** mat22, int size)
{
	int new_size = size / 2;
	for (int i = 0; i < new_size; ++i)
		for (int j = 0; j < new_size; ++j)
			mat11[i][j] = mat[i][j];

	for (int i = 0; i < new_size; ++i)
		for (int j = new_size; j < size; ++j)
			mat12[i][j - new_size] = mat[i][j];

	for (int i = new_size; i < size; ++i)
		for (int j = 0; j < new_size; ++j)
			mat21[i - new_size][j] = mat[i][j];

	for (int i = new_size; i < size; ++i)
		for (int j = new_size; j < size; ++j)
			mat22[i - new_size][j - new_size] = mat[i][j];
}

float** collect_4mat(float** mat11, float** mat12, float** mat21, float** mat22, int size)
{
	int new_size = 2 * size;
	float** mat = create_mat(new_size);

	for (int i = 0; i < size; ++i)
		for (int j = 0; j < size; ++j)
			mat[i][j] = mat11[i][j];

	for (int i = 0; i < size; ++i)
		for (int j = size; j < new_size; ++j)
			mat[i][j] = mat12[i][j - size];

	for (int i = size; i < new_size; ++i)
		for (int j = 0; j < size; ++j)
			mat[i][j] = mat21[i - size][j];

	for (int i = size; i < new_size; ++i)
		for (int j = size; j < new_size; ++j)
			mat[i][j] = mat22[i - size][j - size];

	return mat;
}

float** sum_mat(float** mat1, float** mat2, int size)
{
	float** mat = create_mat(size);

	for (int i = 0; i < size; ++i)
		for (int j = 0; j < size; ++j)
			mat[i][j] = mat1[i][j] + mat2[i][j];

	return mat;
}

float** diff_mat(float** mat1, float** mat2, int size)
{
	float** mat = create_mat(size);

	for (int i = 0; i < size; ++i)
		for (int j = 0; j < size; ++j)
			mat[i][j] = mat1[i][j] - mat2[i][j];

	return mat;
}

float** mult_mat(float** mat1, float** mat2, int size)
{
	if (size <= 64)
		return mult(mat1, mat2, size);
	else
	{
		size /= 2;

		float** a11 = create_mat(size);
		float** a12 = create_mat(size);
		float** a21 = create_mat(size);
		float** a22 = create_mat(size);
		
		float** b11 = create_mat(size);
		float** b12 = create_mat(size);
		float** b21 = create_mat(size);
		float** b22 = create_mat(size);

		split_4mat(mat1, a11, a12, a21, a22, size * 2);
		split_4mat(mat2, b11, b12, b21, b22, size * 2);

		float** p1 = mult_mat(sum_mat(a11, a22, size), sum_mat(b11, b22, size), size);
		float** p2 = mult_mat(sum_mat(a21, a22, size), b11, size);
		float** p3 = mult_mat(a11, diff_mat(b12, b22, size), size);
		float** p4 = mult_mat(a22, diff_mat(b21, b11, size), size);
		float** p5 = mult_mat(sum_mat(a11, a12, size), b22, size);
		float** p6 = mult_mat(sum_mat(a21, a11, size), sum_mat(b11, b12, size), size);
		float** p7 = mult_mat(sum_mat(a12, a22, size), sum_mat(b21, b22, size), size);
		
		float** c11 = sum_mat(diff_mat(p1, p5, size), sum_mat(p4, p7, size), size);
		float** c12 = sum_mat(p3, p5, size);
		float** c21 = sum_mat(p2, p4, size);
		float** c22 = sum_mat(diff_mat(p1, p2, size), sum_mat(p3, p6, size), size);

		return collect_4mat(c11, c12, c21, c22, size);
	}
}
