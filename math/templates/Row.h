#ifndef ROW_H
#define ROW_H
#include <ostream>
#include <vector>
#include <initializer_list>

#include <cassert>
#include <iostream>

template<typename T>
class Row 
{
public:
	Row() = default;
	Row(size_t count_);
	Row(const Row<T> & lhs);
	Row &operator=(const Row<T> & lhs);

	Row &operator=(Row<T> && lhs) noexcept;
	Row(Row<T> && lhs) noexcept;

	~Row();

	T &operator[](size_t col_);
	T operator[](size_t col_) const;
	
	Row &operator=(std::vector<T> lhs);
	Row &operator=(std::initializer_list<T> lhs);
	Row operator+(const Row<T> &lhs);
	Row operator-(const Row<T> &lhs);
	Row operator*(const T &lhs);
	Row operator/(const T &lhs);
	operator std::vector<T>() const { return row; } ;

	void swap(Row<T> & lhs) noexcept;
	size_t size() const;
private:
	std::vector<T> row;
	void check(size_t index) const;
	void check(const Row<T> &lhs) const;
};

template<typename T> 
void Row<T>::check(size_t index) const 
{
	if (!(index < size())) 
	{
		throw std::invalid_argument("index out of range");
	}
}

template<typename T> 
void Row<T>::check(const Row<T> &lhs) const 
{
	if (size()!= lhs.size()) 
	{
		throw std::invalid_argument("different row ranges");
	}
}

template<typename T>
Row<T>::Row(size_t count_) : 
	row(count_, T())
{
	
}

template<typename T>
Row<T>::Row(const Row<T> & lhs) : 
	row(lhs.row) 
{
	
}

template<typename T>
Row<T> &Row<T>::operator=(const Row<T> & lhs) 
{
	row = lhs.row;
	return *this;
}

template<typename T>
Row<T>& Row<T>::operator=(Row<T> && lhs) noexcept
{
	row = std::move(lhs.row);
	return *this;	
}

template<typename T>
Row<T>::Row(Row<T> && lhs) noexcept
{
	row = std::move(lhs.row);
}

template<typename T>
Row<T>::~Row()
{
	row.clear();
}

template<typename T>
T& Row<T>::operator[](size_t col_) 
{
	check(col_);
	return row[col_];
}

template<typename T>
T Row<T>::operator[](size_t col_) const 
{
	check(col_);
	return row[col_];
}

template<typename T>
Row<T> &Row<T>::operator=(std::initializer_list<T> lhs)
{
	row = lhs;
	return *this;
}

template<typename T>
Row<T>& Row<T>::operator=(std::vector<T> lhs) 
{
	row = lhs;
	return *this;
}

template<typename T>
void Row<T>::swap(Row<T> & lhs) noexcept 
{
	row.swap(lhs.row);
}

template<typename T>
size_t Row<T>::size() const
{
	return row.size();
}

template<typename T>
Row<T> Row<T>::operator+(const Row<T> &lhs) 
{
	check(lhs);
	Row<T> tmp(*this);
	for (size_t i = 0; i < lhs.row.size(); ++i) 
	{
		tmp.row[i] += lhs.row[i];
	}
	return tmp;
}

template<typename T>
Row<T> Row<T>::operator-(const Row<T> &lhs) 
{
	check(lhs);
	Row<T> tmp(*this);
	for (size_t i = 0; i < lhs.row.size(); ++i) 
	{
		tmp.row[i] -= lhs.row[i];
	}
	return tmp;
}

template<typename T>
Row<T> Row<T>::operator*(const T &lhs) 
{
	Row<T> tmp(*this);
	for (size_t i = 0; i < row.size(); ++i) 
	{
		tmp.row[i] *= lhs;
	}
	return tmp;
}

template<typename T>
Row<T> Row<T>::operator/(const T &lhs) 
{
	if (lhs == T()) 
	{
		throw std::invalid_argument("division by 0");
	}
	Row<T> tmp(*this);
	for (size_t i = 0; i < row.size(); ++i) 
	{
		tmp.row[i] /= lhs;
	}
	return tmp;
}

template<typename T>
std::ostream &operator<<(std::ostream &out, const Row<T> &row_) 
{
	for (size_t i = 0; i < row_.size(); ++i ) 
	{
		out << row_[i] << ' ';
	}
	return out;
}

template<typename T>
Row<T> swapByPatternBackward(const Row<T> &vctr, const std::vector<std::pair<size_t, size_t> > &pattern_) 
{
	Row<T> tmp(vctr);
	for (auto p = crbegin(pattern_); p != crend(pattern_); ++p) 
	{
		std::swap(tmp[(*p).first], tmp[(*p).second]);
	}
	return tmp;
}

#endif