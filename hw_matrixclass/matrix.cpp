#include"matrix.h"

Matrix::Matrix()
{
	n_columns = 0;
	n_rows = 0;
	arr = nullptr;
}


Matrix::Matrix(const int32_t n_rows_, const int32_t n_columns_, int32_t **arr_)
{
	n_columns = n_columns_;
	n_rows = n_rows_;
	arr = new int32_t* [n_rows_];
	for(int i = 0; i < n_rows_; i++)
	{
		arr[i] = new int32_t [n_columns_];
	}

	for (int i = 0; i < n_rows; i++)
		for(int j = 0; j < n_columns; j++)
		{
			arr[i][j] = arr_[i][j];
		}

}

Matrix::Matrix(const Matrix &A)
{
	n_rows = A.GetRowNum(); 
	n_columns = A.GetColumnNum(); 

	arr = new int32_t* [n_rows];
	for(int i = 0; i < n_rows; i++)
	{
		arr[i] = new int32_t [n_columns];
	}
	for(int i = 0; i < n_rows; i++)
		for(int j = 0; j< n_columns; j++)
		{
			arr[i][j] = A[i][j];
		}
}



int32_t Matrix::GetRowNum() const
{
	return n_rows;
}

int32_t Matrix::GetColumnNum() const
{
	return n_columns;
}

Matrix &Matrix::operator *= (int32_t num)
{
	for(int i = 0; i < n_rows; i++)
		for (int j = 0; j < n_columns; j++)
		{
			arr[i][j] *= num;
		}
	return *this;
}

bool Matrix::operator == (Matrix &B) const
{
	bool res = true;
	if (n_rows != B.GetRowNum() || n_columns != B.GetColumnNum())
	{
		res = false;
	}
	else
	{
		for(int i = 0; i < n_rows; i++)
			for(int j = 0; j < n_columns; j++)
			{
				if(arr[i][j] != B[i][j])
				{
					res = false;
				}
			}
	}
	return res;
}

Matrix &Matrix::operator = (Matrix &B)
{
	if (&B == this)
	{
		return *this;
	}
	else
	{
		for(int i = 0; i < n_rows; i++)
		{
			delete [] arr[i];
		}
		delete [] arr;

		n_columns = B.GetColumnNum();
		n_rows = B.GetRowNum();

		arr = new int32_t* [n_rows];
		for(int i = 0; i < n_rows; i++)
		{
			arr[i] = new int32_t [n_columns];
		}

		for(int i = 0; i < n_rows; i++)
			for(int j = 0; j < n_columns; j++)
			{
				arr[i][j] = B[i][j];
			}
		return *this;
	}
}

Matrix::Row Matrix::operator [] (const int32_t i) 
{
	if(i >= n_rows)
	{
		throw std::out_of_range("");
	}

	Row row(i, n_columns, arr);
	return row;
}

const Matrix::Row  Matrix::operator [] (const int32_t i) const
{
	if(i >= n_rows)
	{
		throw std::out_of_range("");
	}

	Row row(i, n_columns, arr);
	return row;
}

Matrix::~Matrix()
{
	for(int i = 0; i < n_rows; i++)
	{
		delete [] arr[i];
	}
	delete [] arr;
}


Matrix::Row::Row()
{
	length = 0;
	a = nullptr;
}

Matrix::Row::Row(const int32_t i, const int32_t n_columns_, int32_t **arr_)
{
	length = n_columns_;
	a = arr_[i];
}

Matrix::Row::Row(const Matrix::Row &r)
{
	length = r.length;
	a = r.a;
}

int32_t &Matrix::Row::operator [] (const int32_t j)
{
	if(j >= length)
	{
		throw std::out_of_range("");
	}
	return a[j];
}

const int32_t &Matrix::Row::operator [] (const int32_t j) const
{
	if(j >= length)
	{
		throw std::out_of_range("");
	}
	return a[j];
}


