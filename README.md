# Data structures. Matrix multiplication.

I have implemented three ways of matrix multiplication:

The first way is the simplest. It is a common formula from linear algebra.

``` 
for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j < size; ++j)
		{
			mat[i * size + j] = 0;
			for (int k = 0; k < size; ++k)
				mat[i * size + j] += mat1[i * size + k] * mat2[k * size + j];
		}
	}
```

The second way is the result of funcion from BLAS(Basic Linear Algebra Subprograms) libraries. I used cblas_?gemm. In my version, 
there was exactly the cblas_sgemm, because the matrix was of the float type.
```
Example from my code:

  cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, mat_A, N, mat_B, N, 0.0, mat_D, N);
  
```
And finaly the third way. I was trying realized the Strassen's algorithm. More information [here](https://en.wikipedia.org/wiki/Strassen_algorithm) (wiki).
The algorithm was perfect for my quare matrixes of size 2048x2048. It is the one of the most fastest matrix multiplication algorithms.
