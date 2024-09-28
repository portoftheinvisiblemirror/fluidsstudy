#ifndef MY_TENSOR_HPP
#define MY_TENSOR_HPP
#include "vector.hpp"
#include <iostream>
class tensor {
private:
	vector column1, column2, column3;
public:
	tensor()
	{
		vector v1(-1,-1,-1), v2(-1,-1,-1), v3(-1,-1,-1);
		column1 = v1;
		column2 = v2;
		column3 = v3;
	}
	tensor(vector v1, vector v2, vector v3)
	{
		column1 = v1;
		column2 = v2;
		column3 = v3;
	}
	tensor(double inputs[9])
	{
		vector v1(inputs[0], inputs[1], inputs[2]), v2(inputs[3], inputs[4], inputs[5]), v3(inputs[6], inputs[7], inputs[8]);
		column1 = v1;
		column2 = v2;
		column3 = v3;
	}
	tensor(double inputs[3][3])
	{
		vector v1(inputs[0][0], inputs[0][1], inputs[0][2]), v2(inputs[1][0], inputs[1][1], inputs[1][2]),v3(inputs[2][0], inputs[2][1], inputs[2][2]);
		column1 = v1;
		column2 = v2;
		column3 = v3;
	}
	tensor& operator=(const tensor& T)
	{
		if (this != &T)
		{
			column1 = T.column1;
			column2 = T.column2;
			column3 = T.column3;
		}
		return *this;
	}
	tensor transpose()
	{
		vector A(column1.X(), column2.X(), column3.X());
		vector B(column1.Y(), column2.Y(), column3.Y());
		vector C(column1.Z(), column2.Z(), column3.Z());
		tensor transposed(A, B, C);
		return transposed;
	}
	tensor lowertri()
	{
		vector c2 = column2;
		c2.set(0, column2.Y(), column2.Z());
		vector c3 = column3;
		c2.set(0, 0, column3.Z());
		tensor result(column1, c2, c3);
		return result;
	}
	tensor suppertri()
	{
		tensor result = *this;
		result = result - result.lowertri();
		return result;
	}
	vector getColumn(int c)
	{
		switch (c)
		{
		case 1:
			return column1;
		case 2:
			return column2;
		case 3:
			return column3;
		default:
			std::cout << "fail";
			return column1;
		}
	}
	double getValue(int c, int d)
	{
		vector column = (*this).getColumn(c);
		switch (d)
		{
		case 1:
			return column.X();
		case 2:
			return column.Y();
		case 3:
			return column.Z();
		default:
			std::cout << "fail";
			return column.X();
		}
		
	}
	vector operator*(vector v)
	{
		tensor trans = (*this).transpose();
		vector result(trans.getColumn(1)*v, trans.getColumn(2) * v, trans.getColumn(3) * v);
		return result;
	}
	tensor operator*(tensor t)
	{
		vector result1 = (*this) * t.getColumn(1);
		vector result2 = (*this) * t.getColumn(2);
		vector result3 = (*this) * t.getColumn(3);
		tensor result(result1, result2, result3);
		return result;
	}
	tensor operator*(double d)
	{
		tensor result(column1 * d, column2 * d, column2 * d);
		return result;
	}
	tensor operator/(double d)
	{
		tensor result(column1 / d, column2 / d, column2 / d);
		return result;
	}
	tensor operator+(tensor t)
	{
		tensor result(column1 + t.getColumn(1), column2 + t.getColumn(2), column3 + t.getColumn(3));
		return result;
	}
	tensor operator-(tensor t)
	{
		tensor result(column1 - t.getColumn(1), column2 - t.getColumn(2), column3 - t.getColumn(3));
		return result;
	}
};
#endif
