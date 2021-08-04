#include<iostream>
#include"matrix.h"
#include<cassert>
#include<cmath>
#include <algorithm>
#include "test.h"

void FillFunction_first(int32_t n_rows_, int32_t n_columns_, int32_t **a)
{
	for(int i = 0; i < n_rows_; i++)
		for(int j = 0; j < n_columns_; j++)
		{
			a[i][j] = std::max(i,j) + 1;
		}

}

void FillFunction_second(int32_t n_rows_, int32_t n_columns_, int32_t **a)
{
	for(int i = 0; i < n_rows_; i++)
		for(int j = 0; j < n_columns_; j++)
		{
			a[i][j] = std::max(n_rows_ - i, n_columns_ - j) + 1;
		}
}

void ArrayMultiplication(int32_t num, int32_t n_rows_, int32_t n_columns_, int32_t **arr_)
{
	for(int i = 0; i < n_rows_; i++)
		for(int j = 0; j < n_columns_; j++)
		{
			arr_[i][j] *= num;
		}
}


void CreationTest()
{
	int32_t n_rows_ = 5, n_columns_ = 7;
	int32_t **arr_ = new int32_t* [n_rows_];
	for(int i = 0; i< n_rows_;i++)
	{
		arr_[i] = new int32_t [n_columns_];
	}

	FillFunction_first(n_rows_, n_columns_, arr_);

	Matrix A(n_rows_, n_columns_, arr_); 

	assert(A[0][0] == 1 && A[1][1] == 2 && A[2][2] == 3 && A[4][6] == 7);

	assert(A.GetRowNum() == 5 && A.GetColumnNum() == 7); 

	for(int i = 0; i < n_rows_; i++)
	{
		delete [] arr_[i];
	}
	delete [] arr_;
}

void EqualityTest()
{
	int32_t n_rows_ = 4, n_columns_ = 2;
	int32_t **arr_ = new int32_t* [n_rows_];
	for(int i = 0; i< n_rows_;i++)
	{
		arr_[i] = new int32_t [n_columns_];
	}

	FillFunction_first(n_rows_, n_columns_, arr_);

	Matrix A(n_rows_, n_columns_, arr_);
	Matrix B(n_rows_, n_columns_, arr_);

	assert(A == B);

	FillFunction_second(n_rows_, n_columns_, arr_);
	Matrix C(n_rows_, n_columns_, arr_);

	assert(!(A == C));

	for(int i = 0; i<n_rows_; i++)
	{
		delete [] arr_[i];
	}
	delete [] arr_;

	int32_t **a_ = new int32_t* [n_rows_ - 1];
	for(int i = 0; i < n_rows_ - 1; i++)
	{
		a_[i] = new int32_t [n_columns_];
	}

	FillFunction_first(n_rows_ - 1, n_columns_, a_);
	Matrix D(n_rows_ - 1, n_columns_, a_);

	assert(!(A == D));

	for(int i = 0; i<n_rows_ - 1; i++)
	{
		delete [] a_[i];
	}
	delete [] a_;

}

void AssignmentTest()
{
	int32_t n_rows_ = 4, n_columns_ = 2;
	int32_t **arr_ = new int32_t* [n_rows_];
	for(int i = 0; i < n_rows_;i++)
	{
		arr_[i] = new int32_t [n_columns_];
	}
	FillFunction_first(n_rows_, n_columns_, arr_);
	Matrix A(n_rows_, n_columns_, arr_);

	for(int i = 0; i < n_rows_; i++)
	{
		delete [] arr_[i];
	}
	delete [] arr_;

	n_rows_ = 1; 
	n_columns_ = 10;

	int32_t **a_ = new int32_t* [n_rows_];
	for(int i = 0; i < n_rows_; i++)
	{
		a_[i] = new int32_t [n_columns_];
	}

	FillFunction_second(n_rows_, n_columns_, a_);

	Matrix B(n_rows_, n_columns_, a_);

	A = B;

	assert(A == B);
}

void MultiplicationTest()
{
	int32_t n_rows_ = 10, n_columns_ = 20;
	int32_t **arr_ = new int32_t* [n_rows_];
	for(int i = 0; i < n_rows_; i++)
	{
		arr_[i] = new int32_t [n_columns_];
	}
	FillFunction_first(n_rows_, n_columns_, arr_);
	Matrix A(n_rows_, n_columns_, arr_);

	A *= 10;

	ArrayMultiplication(10, n_rows_, n_columns_, arr_);

	Matrix B(n_rows_, n_columns_, arr_);

	assert(A == B);

	for(int i = 0; i < n_rows_; i++)
	{
		delete [] arr_[i];
	}
	delete [] arr_;
}
