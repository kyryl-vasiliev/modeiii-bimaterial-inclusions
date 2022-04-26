#ifndef MATRIX_H
#define MATRIX_H
#include <vector>
#include <utility>
#include <initializer_list>
#include "Row.h"

template<typename T>
class Matrix 
{
public:
	Matrix() = default;
	Matrix(size_t rowCount_, size_t colCount_);
	Matrix(const Matrix &lhs);
	Matrix(Matrix &&lhs) noexcept;
	Matrix &operator=(const Matrix &lhs);
	Matrix &operator=(Matrix &&lhs) noexcept; 

	~Matrix();

	Row<T> &operator[](size_t col_);
	Row<T> operator[](size_t col_) const;
	Matrix &operator=(std::initializer_list<std::vector<T> > m_);
	std::pair<size_t, size_t> findMaxElementByModule(std::pair<size_t, size_t> p1_, std::pair<size_t, size_t> p2_);
	size_t size() const;
private:
	std::vector<Row<T> > matrix;
	void check(size_t index) const;
};

template<typename T>
void Matrix<T>::check(size_t index) const
{
	if (!(index < size())) 
	{
		throw std::invalid_argument("index out of range");
	}
}

template<typename T>
Matrix<T>::Matrix(size_t rowCount_, size_t colCount_): matrix(rowCount_, 
													Row<T>(colCount_) )
{
}

template<typename T>
Matrix<T>::Matrix(const Matrix &lhs) 
{
	matrix.resize(lhs.size());
	for (size_t i = 0; i < lhs.size(); ++i) 
	{
		matrix[i] = lhs[i];
	}
}

template<typename T>
Matrix<T>::Matrix(Matrix &&lhs) noexcept
{
	matrix = std::move(lhs.matrix);
}

template<typename T>
Matrix<T> &Matrix<T>::operator=(Matrix &&lhs) noexcept 
{
	matrix =std::move(lhs.matrix);
	return *this;
}
template<typename T>
Matrix<T> &Matrix<T>::operator=(const Matrix &lhs) 
{
	matrix.resize(lhs.size());
	for (size_t i = 0; i < lhs.size(); ++i) 
	{
		matrix[i] = lhs[i];
	}
	return *this;
}

template<typename T>
Matrix<T>::~Matrix() 
{
	matrix.clear();
}

template<typename T>
std::pair<size_t, size_t> Matrix<T>::findMaxElementByModule(std::pair<size_t, size_t> p1, std::pair<size_t, size_t> p2) 
{
	T max_value = T();
	std::pair<size_t, size_t> max_pos = { p1.first, p1.second };
	for (size_t i = p1.first; i < p2.first; ++i) 
	{
		for (size_t j = p1.second; j < p2.second; ++j) 
		{
			if ((matrix[i][j] > T() ? T(1) : T(-1))*(matrix[i][j]) > max_value)
			{
				max_value = (matrix[i][j] > T() ? T(1) : T(-1)) * (matrix[i][j]);
				max_pos = { i, j };
			}
		}
	}
	return max_pos;
}

template<typename T>
Matrix<T> &Matrix<T>::operator=(std::initializer_list<std::vector<T> > m_) 
{
	if (m_.size() == 0) 
	{
		matrix.clear();
		return *this;
	}

	matrix.resize(m_.size(), Row<T>{ (*m_.begin()).size() }); ///??????
	int j = 0;
	for (auto i = m_.begin(); i != m_.end(); ++i) 
	{
		matrix[j] = *i;
		++j;
	}
	return *this;
}

template<typename T>
Row<T> &Matrix<T>::operator[](size_t col) 
{
	check(col);
	return matrix[col];
}

template<typename T>
Row<T> Matrix<T>::operator[](size_t col) const 
{
	check(col);
	return matrix[col];
}

template<typename T>
size_t Matrix<T>::size() const
{
	return matrix.size();
}

template<typename T>
std::ostream &operator<<(std::ostream &out, const Matrix<T> &m_) 
{
	for (size_t i = 0; i < m_.size(); ++i) 
	{
		out << m_[i] << '\n';
	}
	return out;
}

#endif