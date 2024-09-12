#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <chrono>
long double** allocate2d(const unsigned int latitudes, const unsigned int longitudes, const unsigned int narea, long double*** spheremesh);
long double*** allocate3d(const unsigned int latitudes, const unsigned int longitudes, const unsigned int narea, long double*** spheremesh);
long double*** compute_spheremesh(const double radius, const unsigned int latitudes, const unsigned int longitudes, long double*** spheremesh);
long double** compute_cubes(const double radius, const unsigned int cube, long double** cubes);

long double** allocate2d(const unsigned int x, const unsigned int y, long double** array2d) {
	// Allocate memory blocks of size
	// x i.e., no of 2D Arrays
	array2d = new long double* [x];

	for (unsigned int i = 0; i < x; i++) {

		// Allocate memory blocks for
		// rows of each 2D array
		array2d[i] = new long double[y];
blink
	}

	return array2d;
}

long double*** allocate3d(const unsigned int latitudes, const unsigned int longitudes, const unsigned int narea, long double*** spheremesh) {
	// Dimensions of the 3D array
	int x = latitudes, y = longitudes, z = narea;
	int count = 0;

	// Allocate memory blocks of size
	// x i.e., no of 2D Arrays
	spheremesh = new long double** [x];

	for (int i = 0; i < x; i++) {

		// Allocate memory blocks for
		// rows of each 2D array
		spheremesh[i] = new long double* [y];

		for (int j = 0; j < y; j++) {

			// Allocate memory blocks for
			// columns of each 2D array
			spheremesh[i][j] = new long double[z];
		}
	}

	//for (int i = 0; i < x; i++) {
	//    for (int j = 0; j < y; j++) {
	//        for (int k = 0; k < z; k++) {

	//            // Assign values to the
	//            // memory blocks created
	//            spheremesh[i][j][k] = ++count;
	//            std::cout << spheremesh[i][j][k] << std::endl;
	//        }
	//    }
	//}
	std::cout << latitudes << " " << longitudes << " " << narea << std::endl;
	//int *** ispheremesh = new int **[latitudes];
	//for(unsigned int i = 0; i < latitudes; i ++){
	//  ispheremesh[i] = new int *[longitudes];
	//  for(unsigned int j = 0; j < narea; j ++){
	//    ispheremesh[i][j] = new int [narea];
	//  }
	//}
	printf("Address: %X \n", spheremesh);
	printf("Check value: %d \n", spheremesh[10][10][3]);
	return spheremesh;
}
long double angleof(double x, double y)
{
	long double angle;
	if (x != 0)
	{
		angle = atan(y / x);
		if (y < 0)
		{
			if (x > 0)
				angle += 2 * M_PI;
			else
				angle += M_PI;
		}
		else if (y > 0)
		{
			if (x < 0)
				angle = M_PI + angle;
		}
		else
			if (x > 0)
				angle = 0;
			else
				angle = M_PI;
	}
	else
	{
		if (y > 0)
			angle = M_PI / 2;
		else if (y < 0)
			angle = 3 * M_PI / 2;
		else
			angle = 0;
	}
	return angle;
}
long double*** compute_spheremesh(const double radius, const unsigned int latitudes, const unsigned int longitudes, long double*** spheremesh) {
	//generate sphere mesh
		//const int latitudes=150;//change this number for more precision
		//const int longitudes=150; //change this number for more precision
		//long double spheremesh[latitudes][longitudes][4]; //an array that contains all the coordinates in the mesh. The coordinates will be (x,y,z,area)
	spheremesh = allocate3d(latitudes, longitudes, 4, spheremesh);
	printf("Address 2: %X \n", spheremesh);
	printf("Check value 2: %f \n", spheremesh[10][10][3]);
	long double latitude = M_PI / (2 * latitudes); //note: 90 degrees north will be called 0 radians "latitude" here, and 90 degrees south will be called pi radians "latitude"
	for (unsigned int i = 0; i < latitudes; i++)
	{
		long double z = radius * cos(latitude); //under the normal definition of latitude, this would be sine. But our "latitude" in the program switchees cosines and sines
		double angle = 0;
		long double a = (radius * radius) * (cos(latitude - 0.5 * M_PI / (latitudes)) - cos(latitude + 0.5 * M_PI / (latitudes))) * (2 * M_PI / longitudes);
		for (unsigned int j = 0; j < longitudes; j++)
		{
			long double x = radius * sin(latitude) * cos(angle), y = radius * sin(latitude) * sin(angle);
			spheremesh[i][j][0] = x;
			spheremesh[i][j][1] = y;
			spheremesh[i][j][2] = z;
			spheremesh[i][j][3] = a;
			angle += 2 * M_PI / longitudes;
		}
		latitude += M_PI / (latitudes);
		angle = 0;
	}
	return spheremesh;
}
long double** compute_cubes(const double radius, const unsigned int cube, long double** cubes) {
	//generate sphere mesh
		//const int latitudes=150;//change this number for more precision
		//const int longitudes=150; //change this number for more precision
		//long double spheremesh[latitudes][longitudes][4]; //an array that contains all the coordinates in the mesh. The coordinates will be (x,y,z,area)
	unsigned long int ncubes = (int)8 * cube * cube * cube;
	cubes = allocate2d(ncubes, 3, cubes);
	long double s = radius / cube;// cube side length
	// double cubes[8*(80*80*80)][3]; // center coordinates of an array of cubes, radius/s is size per dim
	//unsigned long int ncubes = 8*(8*8*8);
	long double h = s / 2;
	int i = 0;
	long double z = s / 2;
	long double y = s / 2;
	long double x = s / 2;
	for (int j = 0; j < (int)(2 * cube); j++)
	{
		for (int k = 0; k < (int)(2 * cube); k++)
		{
			for (int l = 0; l < (int)(2 * cube); l++)
			{
				cubes[i][0] = x - radius;
				cubes[i][1] = y - radius;
				cubes[i][2] = z - radius;
				++i;
				x += s;
			}
			y += s;
			x = s / 2;
		}
		z += s;
		y = s / 2;
	}
	return cubes;
}
// Compute the area of the intersection between a sphere and one rectangular box
static long double areal(const double radius, const double xlow, const double xup, const double ylow, const double yup, const double zlow, const double zup, unsigned int latitudes, unsigned int longitudes, long double*** spheremesh)
{
	//pt. 2 find the area
	long double area = 0;
	unsigned int lowlat = latitudes;
	if (zlow < -radius)
		lowlat = (int)ceil((acos(zlow / radius) - M_PI / (2 * latitudes)) / (M_PI / (latitudes)));
	unsigned int uplat = 0;
	if (uplat > radius)
		uplat = (int)floor((acos(zup / radius) - M_PI / (2 * latitudes)) / (M_PI / (latitudes)));
	unsigned int lowlong = 0;
	unsigned int uplong = longitudes;
	long double angles[4] = { angleof(xlow, ylow), angleof(xlow, yup), angleof(xup, ylow), angleof(xup, yup) };
	std::sort(angles, angles + 4);
	if (!(xlow <= 0 && ylow <= 0 && xup >= 0 && yup >= 0)) //find the longitude limits, assuming the projection of the box to the x-y plane does not contain the origin
	{
		lowlong = (int)floor(angles[0] / (2 * M_PI / longitudes));
		uplong = (int)ceil(angles[3] / (2 * M_PI / longitudes));
	}
	if (!(xlow > 0 && ylow < 0 && yup >= 0)) //case 1: the polar axis does not intersect the projection to the x-y plane, or the lower y-bound is the polar axis
	{
		for (unsigned int i = uplat; i < lowlat; i++)   //for loop starts at latitude corresponding to zup -> latitude corresponding to zlow, since the latitudes start at the top of the sphere
		{
			for (unsigned int j = lowlong; j < uplong; j++) //start at the lowest longitude and go to the highest longitude
				if (spheremesh[i][j][0] >= xlow && spheremesh[i][j][0] <= xup && spheremesh[i][j][1] >= ylow && spheremesh[i][j][1] <= yup && spheremesh[i][j][2] >= zlow && spheremesh[i][j][2] <= zup) //check if the point is inside the box. If yes, add the area
					area += spheremesh[i][j][3];
		}
	}
	else //case 2: the polar axis is inside the projection to the x-y plane, or the upper y-bound is the polar axis
	{
		lowlong = (int)ceil(angles[1] / (2 * M_PI / longitudes));
		if (yup == 0)
		{
			lowlong = 1;
		}
		uplong = (int)floor(angles[2] / (2 * M_PI / longitudes));
		for (unsigned int i = uplat; i < lowlat; i++)   //for loop starts at latitude corresponding to zup -> latitude corresponding to zlow, since the latitudes start at the top of the sphere
		{
			for (unsigned int j = 0; j < lowlong; j++) //start from 0 degrees, and end one above the lowest longitude
				if (spheremesh[i][j][0] >= xlow && spheremesh[i][j][0] <= xup && spheremesh[i][j][1] >= ylow && spheremesh[i][j][1] <= yup && spheremesh[i][j][2] >= zlow && spheremesh[i][j][2] <= zup) //check if the point is inside the box. If yes, add the area
					area += spheremesh[i][j][3];
			for (unsigned int j = uplong; j < latitudes; j++)
				if (spheremesh[i][j][0] >= xlow && spheremesh[i][j][0] <= xup && spheremesh[i][j][1] >= ylow && spheremesh[i][j][1] <= yup && spheremesh[i][j][2] >= zlow && spheremesh[i][j][2] <= zup) //check if the point is inside the box. If yes, add the area
					area += spheremesh[i][j][3];
		}
	}
	return area;
}

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
class tensor {
private:
	vector column1, column2, column3;
public:
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
vector velocity(double x, double y, double z)
{
	vector output(x+y+z,x*y*z,z);
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
		return (velocity(x, y+f, z) - velocity(x, y, z)) / f;
	case 3:
		return (velocity(x, y, z+f) - velocity(x, y, z)) / f;
	}
}
double pressure(double x, double y, double z)
{
	return x + y + z;
}
vector pgrad(double x, double y, double z)
{
	double f = 0.00001;
	vector result(pressure(x + f, y, z) - pressure(x, y, z), pressure(x, y+f, z) - pressure(x, y, z), pressure(x , y, z+f) - pressure(x, y, z));
	result = result / f;
	return result;
}
double plaplace(double x, double y, double z,double h)
{
	
	return (pressure(x+h,y,z)+ pressure(x - h, y, z)-2*pressure(x, y, z)+ pressure(x , y+h, z) + pressure(x, y-h, z) - 2 * pressure(x, y, z)+ pressure(x, y, z+h) + pressure(x, y, z-h) - 2 * pressure(x, y, z))/h/h;
}
vector vadvecv(double x, double y, double z)
{
	vector v = velocity(x, y, z);
	vector diff1 = diff(x, y, z, 1);
	vector diff2 = diff(x, y, z, 2);
	vector diff3 = diff(x, y, z, 3);
	vector result(v.X()*diff1.X() +v.Y()*diff2.X()+v.Z()*diff3.X(), v.X() * diff1.Y() + v.Y() * diff2.Y() + v.Z() * diff3.Y(), v.X() * diff1.Z() + v.Y() * diff2.Z() + v.Z() * diff3.Z());
	return result;
}
//vector vfirst(vector* velocities)
//{
//	vector fvelocities[sizeof(velocities)]; //0, 0, 0 -> 2*cube-1, 2*cube-1, 2*cube-1 base 2*cube, 4*cube*cube, 2*cube, 1
//	for (int i = 0; i <= sizeof(fvelocities); ++i)
//	{
//		fvelocities[i] =
//	}
//	
//}
tensor shear(double x, double y, double z)
{
	vector c1 = diff(x, y,z,1);
	vector c2 = diff(x, y, z, 2);
	vector c3 = diff(x, y, z, 3);
	tensor grad(c1, c2, c3);
	return grad + grad.transpose();
}
tensor stress(double p, double x, double y, double z)
{
	double c[] = {p, 0, 0, 0, p, 0, 0, 0, p};
	tensor pressure(c);
	tensor output = pressure + shear(x,y,z);
	return output;
}
static vector forcel(const double radius, const double xlow, const double xup, const double ylow, const double yup, const double zlow, const double zup, unsigned int latitudes, unsigned int longitudes, long double*** spheremesh,double p, double sx, double sy, double sz)
{
	//pt. 3 find the force
	// maybe restrict the range of indices that can potentially intersect with the cube to speed up area calculation.
	vector force(0,0,0);
	unsigned int lowlat = latitudes;
	if (zlow < -radius)
		lowlat = (int)ceil((acos(zlow / radius) - M_PI / (2 * latitudes)) / (M_PI / (latitudes)));
	unsigned int uplat = 0;
	if (uplat > radius)
		uplat = (int)floor((acos(zup / radius) - M_PI / (2 * latitudes)) / (M_PI / (latitudes)));
	unsigned int lowlong = 0;
	unsigned int uplong = longitudes;
	long double angles[4] = { angleof(xlow, ylow), angleof(xlow, yup), angleof(xup, ylow), angleof(xup, yup) };
	std::sort(angles, angles + 4);
	if (!(xlow <= 0 && ylow <= 0 && xup >= 0 && yup >= 0)) //find the longitude limits, assuming the projection of the box to the x-y plane does not contain the origin
	{
		lowlong = (int)floor(angles[0] / (2 * M_PI / longitudes));
		uplong = (int)ceil(angles[3] / (2 * M_PI / longitudes));
	}
	if (!(xlow > 0 && ylow < 0 && yup >= 0)) //case 1: the polar axis does not intersect the projection to the x-y plane, or the lower y-bound is the polar axis
	{
		for (unsigned int i = uplat; i < lowlat; i++)   //for loop starts at latitude corresponding to zup -> latitude corresponding to zlow, since the latitudes start at the top of the sphere
		{
			for (unsigned int j = lowlong; j < uplong; j++) //start at the lowest longitude and go to the highest longitude
				if (spheremesh[i][j][0] >= xlow && spheremesh[i][j][0] <= xup && spheremesh[i][j][1] >= ylow && spheremesh[i][j][1] <= yup && spheremesh[i][j][2] >= zlow && spheremesh[i][j][2] <= zup) //check if the point is inside the box. If yes, add the area
				{
					vector r(spheremesh[i][j][0], spheremesh[i][j][1], spheremesh[i][j][2]);
					force =force+ stress(p, spheremesh[i][j][0] + sx, spheremesh[i][j][1] + sy, spheremesh[i][j][2] + sz)*r*spheremesh[i][j][3]/radius;
				}
		}
	}
	else //case 2: the polar axis is inside the projection to the x-y plane, or the upper y-bound is the polar axis
	{
		lowlong = (int)ceil(angles[1] / (2 * M_PI / longitudes));
		if (yup == 0)
		{
			lowlong = 1;
		}
		uplong = (int)floor(angles[2] / (2 * M_PI / longitudes));
		for (unsigned int i = uplat; i < lowlat; i++)   //for loop starts at latitude corresponding to zup -> latitude corresponding to zlow, since the latitudes start at the top of the sphere
		{
			for (unsigned int j = 0; j < lowlong; j++) //start from 0 degrees, and end one above the lowest longitude
				if (spheremesh[i][j][0] >= xlow && spheremesh[i][j][0] <= xup && spheremesh[i][j][1] >= ylow && spheremesh[i][j][1] <= yup && spheremesh[i][j][2] >= zlow && spheremesh[i][j][2] <= zup) //check if the point is inside the box. If yes, add the area
				{
					vector r(spheremesh[i][j][0], spheremesh[i][j][1], spheremesh[i][j][2]);
					force = force + stress(p, spheremesh[i][j][0]+sx, spheremesh[i][j][1]+sy, spheremesh[i][j][2]+sz) * r * spheremesh[i][j][3] / radius;
				}
			for (unsigned int j = uplong; j < latitudes; j++)
				if (spheremesh[i][j][0] >= xlow && spheremesh[i][j][0] <= xup && spheremesh[i][j][1] >= ylow && spheremesh[i][j][1] <= yup && spheremesh[i][j][2] >= zlow && spheremesh[i][j][2] <= zup) //check if the point is inside the box. If yes, add the area
				{
					vector r(spheremesh[i][j][0], spheremesh[i][j][1], spheremesh[i][j][2]);
					force = force + stress(p, spheremesh[i][j][0] + sx, spheremesh[i][j][1] + sy, spheremesh[i][j][2] + sz) * r * spheremesh[i][j][3] / radius;
				}
		}
	}
	return force;
}
static vector torquel(const double radius, const double xlow, const double xup, const double ylow, const double yup, const double zlow, const double zup, unsigned int latitudes, unsigned int longitudes, long double*** spheremesh, double p,double sx, double sy, double sz)
{
	//pt. 4 find the torque
	// maybe restrict the range of indices that can potentially intersect with the cube to speed up area calculation.
	vector torque(0, 0, 0);
	unsigned int lowlat = latitudes;
	if (zlow < -radius)
		lowlat = (int)ceil((acos(zlow / radius) - M_PI / (2 * latitudes)) / (M_PI / (latitudes)));
	unsigned int uplat = 0;
	if (uplat > radius)
		uplat = (int)floor((acos(zup / radius) - M_PI / (2 * latitudes)) / (M_PI / (latitudes)));
	unsigned int lowlong = 0;
	unsigned int uplong = longitudes;
	long double angles[4] = { angleof(xlow, ylow), angleof(xlow, yup), angleof(xup, ylow), angleof(xup, yup) };
	std::sort(angles, angles + 4);
	if (!(xlow <= 0 && ylow <= 0 && xup >= 0 && yup >= 0)) //find the longitude limits, assuming the projection of the box to the x-y plane does not contain the origin
	{
		lowlong = (int)floor(angles[0] / (2 * M_PI / longitudes));
		uplong = (int)ceil(angles[3] / (2 * M_PI / longitudes));
	}
	if (!(xlow > 0 && ylow < 0 && yup >= 0)) //case 1: the polar axis does not intersect the projection to the x-y plane, or the lower y-bound is the polar axis
	{
		for (unsigned int i = uplat; i < lowlat; i++)   //for loop starts at latitude corresponding to zup -> latitude corresponding to zlow, since the latitudes start at the top of the sphere
		{
			for (unsigned int j = lowlong; j < uplong; j++) //start at the lowest longitude and go to the highest longitude
				if (spheremesh[i][j][0] >= xlow && spheremesh[i][j][0] <= xup && spheremesh[i][j][1] >= ylow && spheremesh[i][j][1] <= yup && spheremesh[i][j][2] >= zlow && spheremesh[i][j][2] <= zup) //check if the point is inside the box. If yes, add the area
				{
					vector r(spheremesh[i][j][0] + sx, spheremesh[i][j][1] + sy, spheremesh[i][j][2] + sz);
					vector n(spheremesh[i][j][0], spheremesh[i][j][1], spheremesh[i][j][2]);
					n = n / radius;
					torque = torque + r%(stress(p, spheremesh[i][j][0] + sx, spheremesh[i][j][1] + sy, spheremesh[i][j][2] + sz)*n)* spheremesh[i][j][3];
				}
		}
	}
	else //case 2: the polar axis is inside the projection to the x-y plane, or the upper y-bound is the polar axis
	{
		lowlong = (int)ceil(angles[1] / (2 * M_PI / longitudes));
		if (yup == 0)
		{
			lowlong = 1;
		}
		uplong = (int)floor(angles[2] / (2 * M_PI / longitudes));
		for (unsigned int i = uplat; i < lowlat; i++)   //for loop starts at latitude corresponding to zup -> latitude corresponding to zlow, since the latitudes start at the top of the sphere
		{
			for (unsigned int j = 0; j < lowlong; j++) //start from 0 degrees, and end one above the lowest longitude
				if (spheremesh[i][j][0] >= xlow && spheremesh[i][j][0] <= xup && spheremesh[i][j][1] >= ylow && spheremesh[i][j][1] <= yup && spheremesh[i][j][2] >= zlow && spheremesh[i][j][2] <= zup) //check if the point is inside the box. If yes, add the area
				{
					vector r(spheremesh[i][j][0] + sx, spheremesh[i][j][1] + sy, spheremesh[i][j][2] + sz);
					vector n(spheremesh[i][j][0], spheremesh[i][j][1], spheremesh[i][j][2]);
					n = n / radius;
					torque = torque + r % (stress(p, spheremesh[i][j][0] + sx, spheremesh[i][j][1] + sy, spheremesh[i][j][2] + sz) * n) * spheremesh[i][j][3];
				}
			for (unsigned int j = uplong; j < latitudes; j++)
				if (spheremesh[i][j][0] >= xlow && spheremesh[i][j][0] <= xup && spheremesh[i][j][1] >= ylow && spheremesh[i][j][1] <= yup && spheremesh[i][j][2] >= zlow && spheremesh[i][j][2] <= zup) //check if the point is inside the box. If yes, add the area
				{
					vector r(spheremesh[i][j][0] + sx, spheremesh[i][j][1] + sy, spheremesh[i][j][2] + sz);
					vector n(spheremesh[i][j][0], spheremesh[i][j][1], spheremesh[i][j][2]);
					n = n / radius;
					torque = torque + r % (stress(p, spheremesh[i][j][0] + sx, spheremesh[i][j][1] + sy, spheremesh[i][j][2] + sz) * n) * spheremesh[i][j][3];
				}
		}
	}
	return torque;
}

int checkcube(double radius, double s, double x, double y, double z)
{
	//center of cube
	vector center(x, y, z);
	if (center.length() + s<radius || center.length() - s >radius)
	{
		return -1;
	}
	////create array of cube vertices as vectors
	//vector cube[8];
	//int index = 0;
	//double xin = x - s / 2;
	//double yin = y - s / 2;
	//double zin = z - s / 2;
 // // Maybe as you construct your box coordinates, return immediately if the box is too far away.
	//for (int i = 0; i <= 1; i++)
	//{
	//	yin = y - s / 2;
	//	for (int j = 0; j <= 1; j++)
	//	{
	//		xin = x - s / 2;
	//		for (int k = 0; k <= 1; k++)
	//		{
	//			cube[index].set(xin,yin,zin);
	//			xin += s;
	//		}
	//		yin += s;
	//	}
	//	zin += s;
	//}
	////array of minkowski differences
	//vector m[4];
	////start with a given support point
	//vector support1 = cube[1];
	////compute 1st minkowski difference
	//m[0] = cube[1] + (cube[1]-center) * (radius / (cube[1]-center).length());
	//
	////generate simplex
	//int k = 1;
	//for (int i = 1; i < 4; ++i)
	//{
	//	m[3] = m[2];
	//	m[2] = m[1];
	//	m[1] = m[0];
	//	m[0] = cube[support(cube, center, (m[0]+m[1]+ m[2] + m[3])/k)] + m[0] * (radius / m[0].length());
	//}

	//if (tetrahedral(m[0], m[1], m[2], m[3]) == -1)
	//	return -1;
	//else
	return 0;
}
vector force(double x0, double y0, double z0, unsigned int longitudes, unsigned int latitudes, long double radius, int cube, double p, long double **cubes)
{
	long double s = radius / cube;
	long double h = s / 2;
	vector sum(0,0,0);
	auto start = std::chrono::high_resolution_clock::now();
	long double*** spheremesh = nullptr;
	spheremesh = compute_spheremesh(radius, latitudes, longitudes, spheremesh);
	auto finish = std::chrono::high_resolution_clock::now();
	auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
	std::cout << "time of sphere: " << milliseconds.count() / 1000.0 << "s\n";
	start = std::chrono::high_resolution_clock::now();
	int l = 0;
	assert(spheremesh != nullptr);
	for (unsigned long int kk = 0; kk < 8*cube*cube*cube; kk++)
	{
		if (checkcube(radius, s, cubes[kk][0], cubes[kk][1], cubes[kk][2]) == 0)
		{
			double a = cubes[kk][0] - h, b = cubes[kk][0] + h, c = cubes[kk][1] - h, d = cubes[kk][1] + h, e = cubes[kk][2] - h, f = cubes[kk][2] + h;
			sum = sum+ forcel(radius, a, b, c, d, e, f, latitudes, longitudes, spheremesh,p,x0,y0,z0);
		}
	}
	finish = std::chrono::high_resolution_clock::now();
	milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
	std::cout << "time of area: " << milliseconds.count() / 1000.0 << "s\n";
	std::cout << "Expected area: " << 4 * M_PI * radius * radius << std::endl;
	std::cout << "Calculated force: " << "(" << sum.X() << ", " << sum.Y() << ", " << sum.Z() << ")" << std::endl;
	return sum;
}
vector torque(double x0, double y0, double z0, unsigned int longitudes, unsigned int latitudes, long double radius, int cube, double p,long double **cubes)
{
	long double s = radius / cube;
	long double h = s / 2;
	vector sum(0, 0, 0);
	auto start = std::chrono::high_resolution_clock::now();
	long double*** spheremesh = nullptr;
	spheremesh = compute_spheremesh(radius, latitudes, longitudes, spheremesh);
	auto finish = std::chrono::high_resolution_clock::now();
	auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
	std::cout << "time of sphere: " << milliseconds.count() / 1000.0 << "s\n";
	start = std::chrono::high_resolution_clock::now();
	int l = 0;
	assert(spheremesh != nullptr);
	for (unsigned long int kk = 0; kk < 8*cube*cube*cube; kk++)
	{
		if (checkcube(radius, s, cubes[kk][0], cubes[kk][1], cubes[kk][2]) == 0)
		{
			double a = cubes[kk][0] - h, b = cubes[kk][0] + h, c = cubes[kk][1] - h, d = cubes[kk][1] + h, e = cubes[kk][2] - h, f = cubes[kk][2] + h;
			sum = sum + torquel(radius, a, b, c, d, e, f, latitudes, longitudes, spheremesh, p, x0, y0, z0);
		}
	}
	finish = std::chrono::high_resolution_clock::now();
	milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
	std::cout << "time of area: " << milliseconds.count() / 1000.0 << "s\n";
	std::cout << "Expected area: " << 4 * M_PI * radius * radius << std::endl;
	std::cout << "Calculated torque: " << "(" << sum.X() << ", " << sum.Y() << ", " << sum.Z() << ")" << std::endl;
	return sum;
}
int main()
{
	std::cout << "Input xyz coordinates of the sphere in this format (separated by spaces no commas): x y z\n";
	double spherex, spherey, spherez;
	std::cin >> spherex >> spherey >> spherez;
	std::cout << "Input the number of longitudes used in generating the spheremesh.";
	int longitudes, latitudes;
	std::cin >> longitudes;
	std::cout << "Input the number of latitudes used in generating the spheremesh.";
	std::cin >> latitudes;
	std::cout << "What is the radius of the sphere?";
	long double radius;
	std::cin >> radius;
	std::cout << "For the cartesian grid, by what factor should the step distance be smaller than the radius? (If you want 8*8 cubes, input 2, if you want 8*27 cubes, input 3 and so on)";
	int cube;
	std::cin >> cube;
	std::cout << "What is the pressure?";
	double pressure;
	std::cin >> pressure;
	std::cout << "What is the mass?";
	double mass;
	std::cin >> mass;
	std::cout << "WHat is the time period?";
	double time;
	std::cin >> time;
	double increment = time / 100;
	long double** cubes = nullptr;
	cubes = compute_cubes(radius, cube, cubes);
	vector* velocities;
	int l = sizeof(cubes);
	velocities = new vector  [l];
	assert(cubes != nullptr);
	vector ucm(0, 0, 0);
	vector xyz(spherex, spherey, spherez);
	vector angv(0, 0, 0);
	for (int i = 0; i < 100; i++)
	{
		vector F = force(xyz.X(), xyz.Y(), xyz.Z(), longitudes, latitudes, radius, cube, pressure,cubes);
		angv = angv+torque(xyz.X(), xyz.Y(), xyz.Z(), longitudes, latitudes, radius, cube, pressure,cubes)*increment/(2*mass*radius*radius/5);
		xyz = xyz + ucm * increment + F * (increment) * (increment) / (2 * mass);
	}
	std::cout << xyz.X() <<" " << xyz.Y()<<" " << xyz.Z();
	
	return 0;
}
