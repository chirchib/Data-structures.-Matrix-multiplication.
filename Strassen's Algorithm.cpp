#include "Strassen's Algorithm.h"
#include <iostream>
#include <ctime>
using namespace std;

float** create_mat(int size_mat)
{
	float** mat = new float* [size_mat];
	for (int i = 0; i < size_mat; ++i)
	{
		mat[i] = new float[size_mat];
		for (int j = 0; j < size_mat; ++j)
			mat[i][j] = 0;
	}

	return mat;
}

float** mult(float** mat1, float** mat2, int size_mat)
{
	float** mat_res = create_mat(size_mat);
	for (int i = 0; i < size_mat; ++i)
		for (int j = 0; j < size_mat; ++j)
			for (int k = 0; k < size_mat; ++k)
				mat_res[i][j] += mat1[i][k] * mat2[k][j];

	return mat_res;
}

void split_4mat(float** mat, float** mat11, float** mat12, float** mat21, float** mat22, int size_mat)
{
	int new_size_mat = size_mat / 2;
	for (int i = 0; i < new_size_mat; ++i)
		for (int j = 0; j < new_size_mat; ++j)
			mat11[i][j] = mat[i][j];

	for (int i = 0; i < new_size_mat; ++i)
		for (int j = new_size_mat; j < size_mat; ++j)
			mat12[i][j - new_size_mat] = mat[i][j];

	for (int i = new_size_mat; i < size_mat; ++i)
		for (int j = 0; i < new_size_mat; ++j)
			mat21[i - new_size_mat][j] = mat[i][j];

	for (int i = new_size_mat; i < size_mat; ++i)
		for (int j = new_size_mat; i < size_mat; ++j)
			mat22[i - new_size_mat][j - new_size_mat] = mat[i][j];
}

float** collect_4mat(float** mat11, float** mat12, float** mat21, float** mat22, int size_mat)
{
	int new_size_mat = 2 * size_mat;
	float** mat = create_mat(new_size_mat);

	for (int i = 0; i < size_mat; ++i)
		for (int j = 0; j < size_mat; ++j)
			mat[i][j] = mat11[i][j];

	for (int i = 0; i < size_mat; ++i)
		for (int j = size_mat; j < new_size_mat; ++j)
			mat[i][j] = mat12[i][j - size_mat];

	for (int i = size_mat; i < new_size_mat; ++i)
		for (int j = 0; i < size_mat; ++j)
			mat[i][j] = mat21[i - size_mat][j];

	for (int i = size_mat; i < new_size_mat; ++i)
		for (int j = size_mat; i < new_size_mat; ++j)
			mat[i][j] = mat22[i - size_mat][j - size_mat];

	return mat;
}

float** sum_mat(float** mat1, float** mat2, int size_mat)
{
	float** mat = create_mat(size_mat);
	for (int i = 0; i < size_mat; ++i)
		for (int j = 0; j < size_mat; ++j)
			mat[i][j] = mat1[i][j] + mat2[i][j];

	return mat;
}

float** diff_mat(float** mat1, float** mat2, int size_mat)
{
	float** mat = create_mat(size_mat);
	for (int i = 0; i < size_mat; ++i)
		for (int j = 0; j < size_mat; ++j)
			mat[i][j] = mat1[i][j] - mat2[i][j];

	return mat;
}

float** mult_mat(float** mat1, float** mat2, int size_mat)
{
	if (size_mat <= 64)
		return mult(mat1, mat2, size_mat);
	else
	{
		size_mat /= 2;

		float** a11 = create_mat(size_mat);
		float** a12 = create_mat(size_mat);
		float** a21 = create_mat(size_mat);
		float** a22 = create_mat(size_mat);
		
		float** b11 = create_mat(size_mat);
		float** b12 = create_mat(size_mat);
		float** b21 = create_mat(size_mat);
		float** b22 = create_mat(size_mat);

		split_4mat(mat1, a11, a12, a21, a22, size_mat * 2);
		split_4mat(mat2, b11, b12, b21, b22, size_mat * 2);

		float** p1 = mult_mat(sum_mat(a11, a22, size_mat), sum_mat(b11, b22, size_mat), size_mat);
		float** p2 = mult_mat(sum_mat(a21, a22, size_mat), b11, size_mat);
		float** p3 = mult_mat(a11, diff_mat(b12, b22, size_mat), size_mat);
		float** p4 = mult_mat(a22, diff_mat(b21, b11, size_mat), size_mat);
		float** p5 = mult_mat(sum_mat(a11, a12, size_mat), b22, size_mat);
		float** p6 = mult_mat(sum_mat(a21, a11, size_mat), sum_mat(b11, b12, size_mat), size_mat);
		float** p7 = mult_mat(sum_mat(a12, a22, size_mat), sum_mat(b21, b22, size_mat), size_mat);
		
		float** c11 = sum_mat(diff_mat(p1, p5, size_mat), sum_mat(p4, p7, size_mat), size_mat);
		float** c12 = sum_mat(p3, p5, size_mat);
		float** c21 = sum_mat(p2, p4, size_mat);
		float** c22 = sum_mat(diff_mat(p1, p2, size_mat), sum_mat(p3, p6, size_mat), size_mat);

		return collect_4mat(c11, c12, c21, c22, size_mat);
	}
}
