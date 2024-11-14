#ifndef MY_VECTOR_HPP
#define MY_VECTOR_HPP
#include <cmath>

class vector {
private:
	double x, y, z; //coordinates
public:
	vector() {
		x = 0;
		y = 0;
		z = 0;
	}
	vector(double a, double b, double c) //constructor
	{
		x = a;
		y = b;
		z = c;
	}

	void set(double a, double b, double c)
	{
		x = a;
		y = b;
		z = c;
	}

	//accessor methods
	double X()
	{
		return x;
	}
	double Y()
	{
		return y;
	}
	double Z()
	{
		return z;
	}

	//dot product
	double operator*(vector const& v)
	{
		double result;
		result = x * v.x + y * v.y + z * v.z;
		return result;
	}
	//cross product
	vector operator%(vector const& v)
	{
		vector result((y * v.z - z * v.y), - (x * v.z - z * v.x), (x * v.y - y * v.x));
		return result;
	}
	//addition and subtraction
	vector operator+(vector const& v)
	{
		vector result(v.x + x, v.y + y, v.z + z);
		return result;
	}
	vector operator-(vector const& v)
	{
		vector result(x - v.x, y - v.y, z - v.z);
		return result;
	}
	//scalar multiplication
	vector operator*(double c)
	{
		vector result(x * c, y * c, z * c);
		return result;
	}
	//vector operator*(tensor t)
	//{
	//	vector result(*this * t.getColumn(1), *this * t.getColumn(2), *this * t.getColumn(3));
	//	return result;
	//}
	//scalar division
	vector operator/(double d)
	{
		vector result(x / d, y / d, z / d);
		return result;
	}
	//copying
	vector &operator=(const vector &V)
	{
		if (this != &V)
		{
			x = V.x;
			y = V.y;
			z = V.z;
		}
		return *this;
	}
	double length()
	{
		return sqrt(x * x + y * y + z * z);
	}
};
#endif
