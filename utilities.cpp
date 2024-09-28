#include "utilities.hpp"
#include "vector.hpp"
tensor shear(double x, double y, double z)
{
	vector c1 = diff(x, y, z, 1);
	vector c2 = diff(x, y, z, 2);
	vector c3 = diff(x, y, z, 3);
	tensor grad(c1, c2, c3);
	return grad + grad.transpose();
}
tensor stress(double p, double x, double y, double z)
{
	double c[] = { p, 0, 0, 0, p, 0, 0, 0, p };
	tensor pressure(c);
	tensor output = pressure + shear(x, y, z);
	return output;
}
vector diff(double x, double y, double z, int index)
{
	double f = 0.00001;
	switch (index)
	{
	case 1:
		return (velocity(x + f, y, z) - velocity(x, y, z)) / f;
	case 2:
		return (velocity(x, y + f, z) - velocity(x, y, z)) / f;
	case 3:
		return (velocity(x, y, z + f) - velocity(x, y, z)) / f;
	}
}
vector velocity(double x, double y, double z)
{
	vector output(x + y + z, x * y * z, z);
	return output;
}
vector divvadvecv(double x, double y, double z)
{
	vector v = velocity(x, y, z);
	vector diff1 = diff(x, y, z, 1);
	vector diff2 = diff(x, y, z, 2);
	vector diff3 = diff(x, y, z, 3);
	vector result(v.X() * diff1.X() + v.Y() * diff2.X() + v.Z() * diff3.X(), v.X() * diff1.Y() + v.Y() * diff2.Y() + v.Z() * diff3.Y(), v.X() * diff1.Z() + v.Y() * diff2.Z() + v.Z() * diff3.Z());
	return result;
}