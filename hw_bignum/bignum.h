#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include "vector.h"
#include <iomanip>
#include <sstream>

const uint16_t POWER = 4;
const size_t BASE = pow(10, POWER);

bool long_isdigit(const std::string &num) 
{
	if (num.empty() || (num.size() == 1 && num[0] == '-')) 
	{
		return false;
	}
	auto it_begin = (num[0] == '-') ? next(begin(num)) : begin(num);
	auto it = std::find_if(it_begin, end(num), [](const char &c) {
			return !isdigit(c);
			});
	return it == num.end();
}

int64_t my_stoi(const char *first, const char *last) 
{
	int64_t resnum = 0;
	while (first != last) 
	{
		resnum = resnum * 10 + *first - '0';
		++first;
	}
	return resnum;
}

class BigNum 
{
	private:
		Vector<size_t> _data_;
		bool _sign_; // true == +

	public:

		BigNum() = default;
		BigNum(const BigNum &num) = default;
		~BigNum() = default;
		BigNum &operator=(const BigNum &num) = default;

		BigNum(BigNum &&other) : _sign_(std::move(other._sign_))
		{
			_data_ = other._data_;
		}

		BigNum &operator=(BigNum &&moved) 
		{
			if (this == &moved)
			{
				return *this;
			}
			_sign_ = std::move(moved._sign_);
			_data_ = moved._data_;
			return *this;
		}

		BigNum(int num)
		{
			_sign_ = num > 0;
			num = _sign_ ? num : -num;
			size_t length = (num == 0) ? 1u : static_cast<size_t>(log10(num) + 1);
			_data_.resize(length % POWER == 0 ? length / POWER : length / POWER + 1);

			for (size_t &digit : _data_) 
			{
				digit = num % BASE;
				num /= BASE;
			}
		}


		BigNum(const std::string &data_str)
		{
			if (!long_isdigit(data_str)) 
			{
				throw std::invalid_argument(data_str);
			} 
			else 
			{
				_sign_ = data_str[0] != '-';
				auto starts = (_sign_) ? data_str.c_str() : data_str.c_str() + 1;
				auto ends = data_str.c_str() + data_str.size();

				if ((ends - starts) % POWER == 0) 
				{
					_data_.reserve((ends - starts) / POWER);
				} 
				else 
				{
					_data_.reserve((ends - starts) / POWER + 1);
				}

				auto pointer = ends;
				for (pointer; pointer > starts + POWER; pointer -= POWER) 
				{
					_data_.push_back(my_stoi(pointer - POWER, pointer));
				}
				if (starts != pointer) 
				{
					_data_.push_back(my_stoi(starts, pointer));
				}
			}
		}


		friend bool operator==(const BigNum &num1, const BigNum &num2)
		{
			if (num1._data_.size() == 1 && num1._data_.size() == 1 && num1._data_[0] == 0 && num2._data_[0] == 0) 
			{
				return true;
			}

			if (num1._sign_ != num2._sign_ ||
					num1._data_.size() != num2._data_.size()) 
			{
				return false;
			}

			for (size_t i = 0; i < num1._data_.size(); ++i) 
			{
				if (num1._data_[i] != num2._data_[i]) 
				{
					return false;
				}
			}
			return true;
		}


		friend bool operator!=(const BigNum &num1, const BigNum &num2)
		{
			return !(num1 == num2);
		}


		friend bool operator>(const BigNum &num1, const BigNum &num2)
		{
			return !(num1 <= num2);
		}


		friend bool operator<(const BigNum &num1, const BigNum &num2)
		{
			if (num1._sign_ != num2._sign_) 
			{
				return num2._sign_;
			} 
			else if (num1._data_.size() != num2._data_.size()) 
			{
				return !num1._sign_ ^ (num1._data_.size() < num2._data_.size());
			} 
			else 
			{
				return !(num1 == num2) && (!num1._sign_) ^ compare_mod(num1, num2);
			}
		}


		friend bool operator>=(const BigNum &num1, const BigNum &num2)
		{
			return !(num1 < num2);
		}


		friend bool operator<=(const BigNum &num1, const BigNum &num2)
		{
			return (num1 < num2 || num1 == num2);
		}


		friend BigNum abs_plus(const BigNum &num1, const BigNum &num2, bool sign)
		{
			bool res = 0;
			BigNum resnum;
			resnum._sign_ = sign;
			resnum._data_.reserve(num1._data_.size() + 1);

			for (size_t i = 0; i < num2._data_.size(); ++i) 
			{
				if (BASE - res - num1._data_[i] > num2._data_[i]) 
				{
					resnum._data_.push_back(num1._data_[i] + num2._data_[i] + res);
					res = 0;
				} 
				else 
				{
					resnum._data_.push_back(num1._data_[i] - (BASE - res - num2._data_[i]));
					res = 1;
				}
			}
			for (size_t i = num2._data_.size(); i < num1._data_.size(); ++i) 
			{
				if (BASE - res > num1._data_[i]) 
				{
					resnum._data_.push_back(num1._data_[i] + res);
					res = 0;
				} 
				else 
				{
					resnum._data_.push_back(num1._data_[i] - (BASE - res));
				}
			}
			if (res) 
			{
				resnum._data_.push_back(res);
			}

			return resnum;
		}


		friend BigNum abs_minus(const BigNum &num1, const BigNum &num2, bool sign)
		{
			bool res = 0;
			BigNum resnum;
			resnum._sign_ = sign;
			resnum._data_.reserve(num1._data_.size());
			for (size_t i = 0; i < num2._data_.size(); ++i) 
			{
				if (num1._data_[i] < num2._data_[i] + res) 
				{
					resnum._data_.push_back(num1._data_[i] + (BASE - num2._data_[i] - res));
					res = 1;
				} 
				else 
				{
					resnum._data_.push_back(num1._data_[i] - res - num2._data_[i]);
					res = 0;
				}
			}

			for (size_t i = num2._data_.size(); i < num1._data_.size(); ++i) 
			{
				if (num1._data_[i] < res) 
				{
					resnum._data_.push_back(BASE - res);
				} 
				else 
				{
					resnum._data_.push_back(num1._data_[i] - res);
					res = 0;
				}
			}
			resnum.del_beg_zeroes();
			return resnum;
		}


		const BigNum operator-() const
		{
			BigNum resnum(*this);
			resnum._sign_ = !resnum._sign_;
			return resnum;
		}


		friend BigNum operator+(const BigNum &num1, const BigNum &num2)
		{
			if (num1._sign_ == num2._sign_) 
			{
				if (num1._data_.size() < num2._data_.size()) 
				{
					return abs_plus(num2, num1, num1._sign_);
				} 
				else 
				{
					return abs_plus(num1, num2, num1._sign_);
				}
			} 
			else if (compare_mod(num1, num2)) 
			{
				return abs_minus(num2, num1, num2._sign_);
			} 
			else 
			{
				return abs_minus(num1, num2, num1._sign_);
			}
		}


		friend BigNum operator-(const BigNum &num1, const BigNum &num2)
		{
			if (num1._sign_ != num2._sign_) 
			{
				if (num1._data_.size() < num2._data_.size()) 
				{
					return abs_plus(num2, num1, num1._sign_);
				} 
				else 
				{
					return abs_plus(num1, num2, num1._sign_);
				}
			} 
			else if (compare_mod(num1, num2)) 
			{
				return abs_minus(num2, num1, !num2._sign_);
			} 
			else
			{
				return abs_minus(num1, num2, num1._sign_);
			}
		}


		friend BigNum operator*(const BigNum &num1, const BigNum &num2)
		{
			if(num1._data_.size() <  num2._data_.size()) 
			{
				return num2 * num1;
			}
			BigNum resnum;
			resnum._data_.resize(num1._data_.size() + num2._data_.size() + 1);
			resnum._sign_ = num1._sign_ == num2._sign_;

			for (size_t i = 0; i < num1._data_.size(); ++i) 
			{
				size_t res = 0;
				for (size_t j = 0; j < num2._data_.size(); ++j) 
				{
					int64_t prod = res + resnum._data_[i + j] + num1._data_[i] * num2._data_[j];
					resnum._data_[i + j] = prod % BASE;
					res = prod / BASE;
				}
				size_t j = num2._data_.size();
				while(res != 0) 
				{
					resnum._data_[i + j++] = res % BASE;
					res /= BASE;
				}
			}
			resnum.del_beg_zeroes();
			return resnum;
		}


		friend BigNum operator^(BigNum num1, size_t num2)
		{
			BigNum resnum(1);
			resnum._sign_ |= num2 % 2 == 0; 
			if (num2 == 0) 
			{
				return resnum;
			}
			while (num2) 
			{
				if (num2 % 2 == 0) 
				{
					num2 /= 2;
					num1 = num1 * num1;
				} 
				else 
				{
					--num2;
					resnum = resnum * num1;
				}
			}
			return resnum;
		}


		int rev(int num, int lg_n) 
		{
			int res = 0;
			for (int i = 0; i < lg_n; ++i)
				if (num & (1 << i))
					res |= 1 << (lg_n - 1 - i);
			return res;
		}




		friend bool compare_mod(const BigNum &num1, const BigNum &num2)
		{
			if (num1._data_.size() != num2._data_.size()) 
			{
				return num1._data_.size() < num2._data_.size();
			} 
			else 
			{
				for (size_t i = num1._data_.size(); i > 0; --i) {
					if (num1._data_[i - 1] < num2._data_[i - 1]) {
						return true;
					}
					if (num1._data_[i - 1] > num2._data_[i - 1])
						return false;
				}
				return false;
			}
		}


		friend std::ostream &operator<<(std::ostream &os, const BigNum &num2)
		{
			if (!num2._sign_) 
			{
				os << "-";
			}
			os << num2._data_[num2._data_.size() - 1];
			for (size_t i = num2._data_.size() - 1; i > 0; --i) 
			{
				os << std::setfill('0') << std::setw(POWER);
				os << num2._data_[i-1];
			}
			return os;
		}

		operator std::string() const
		{
			std::stringstream ss;
			ss << *this;
			return ss.str();
		}

		void del_beg_zeroes()
		{
			while (this->_data_.size() > 1 && this->_data_.back() == 0) 
			{
				this->_data_.pop_back();
			}
		}

};
