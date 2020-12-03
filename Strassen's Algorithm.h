class Strassen
{
public:
	float** mult(float**, float**, int);

	float** mult_mat(float**, float**, int);
private:
	float** create_mat(int);

	float** create_mat_rand(int);
	
	void split_4mat(float**, float**, float*, float**, float**, int);

	float** collect_4mat(float**, float**, float**, float**, int);

	float** sum_mat(float**, float**, int);

	float** diff_mat(float**, float**, int);
};
