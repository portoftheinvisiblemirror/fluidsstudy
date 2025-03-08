#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <vector>

#include <stdlib.h>

#include "utilities.hpp"
#include "mem.hpp"
#include "vector.hpp"
#include "tensor.hpp"
#include "vtk.hpp"
#include <Dense>

double*** compute_spheremesh(const double radius, const unsigned int latitudes, const unsigned int longitudes, double*** spheremesh);
double** compute_cubes(const double radius, const unsigned int cube, double** cubes);

double angleof(double x, double y)
{
	double angle;
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
double*** compute_spheremesh(const double radius, const unsigned int latitudes, const unsigned int longitudes) {
	//generate sphere mesh
	double *** spheremesh = allocate3d(latitudes, longitudes, 4);
	double latitude = M_PI / (2 * latitudes); //note: 90 degrees north will be called 0 radians "latitude" here, and 90 degrees south will be called pi radians "latitude"
	for (unsigned int i = 0; i < latitudes; i++)
	{
		double z = radius * cos(latitude); //under the normal definition of latitude, this would be sine. But our "latitude" in the program switchees cosines and sines
		double angle = 0;
		double a = (radius * radius) * (cos(latitude - 0.5 * M_PI / (latitudes)) - cos(latitude + 0.5 * M_PI / (latitudes))) * (2 * M_PI / longitudes);
		for (unsigned int j = 0; j < longitudes; j++)
		{
			double x = radius * sin(latitude) * cos(angle), y = radius * sin(latitude) * sin(angle);
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
double** compute_cubes(const double radius, const unsigned int cube, double** cubes) {
	//generate sphere mesh
		//const int latitudes=150;//change this number for more precision
		//const int longitudes=150; //change this number for more precision
		//double spheremesh[latitudes][longitudes][4]; //an array that contains all the coordinates in the mesh. The coordinates will be (x,y,z,area)
	unsigned long int ncubes = (int)8 * cube * cube * cube;
	cubes = allocate2d(ncubes, 3, cubes);
	double s = radius / cube;// cube side length
	// double cubes[8*(80*80*80)][3]; // center coordinates of an array of cubes, radius/s is size per dim
	//unsigned long int ncubes = 8*(8*8*8);
	double h = s / 2;
	int i = 0;
	double z = s / 2;
	double y = s / 2;
	double x = s / 2;
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
static double areal(const double radius, const double xlow, const double xup, const double ylow, const double yup, const double zlow, const double zup, unsigned int latitudes, unsigned int longitudes, double*** spheremesh)
{
	//pt. 2 find the area
	double area = 0;
	int lowlat = latitudes;
	if (zlow < -radius)
		lowlat = (int)ceil((acos(zlow / radius) - M_PI / (2 * latitudes)) / (M_PI / (latitudes)));
	int uplat = 0;
	if (uplat > radius)
		uplat = (int)floor((acos(zup / radius) - M_PI / (2 * latitudes)) / (M_PI / (latitudes)));
	int lowlong = 0;
	int uplong = longitudes;
	double angles[4] = { angleof(xlow, ylow), angleof(xlow, yup), angleof(xup, ylow), angleof(xup, yup) };
	std::sort(angles, angles + 4);
	if (!(xlow <= 0 && ylow <= 0 && xup >= 0 && yup >= 0)) //find the longitude limits, assuming the projection of the box to the x-y plane does not contain the origin
	{
		lowlong = (int)floor(angles[0] / (2 * M_PI / longitudes));
		uplong = (int)ceil(angles[3] / (2 * M_PI / longitudes));
	}
	if (!(xlow > 0 && ylow < 0 && yup >= 0)) //case 1: the polar axis does not intersect the projection to the x-y plane, or the lower y-bound is the polar axis
	{
		for (int i = uplat; i < lowlat; i++)   //for loop starts at latitude corresponding to zup -> latitude corresponding to zlow, since the latitudes start at the top of the sphere
		{
			for (int j = lowlong; j < uplong; j++) //start at the lowest longitude and go to the highest longitude
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
		for (int i = uplat; i < lowlat; i++)   //for loop starts at latitude corresponding to zup -> latitude corresponding to zlow, since the latitudes start at the top of the sphere
		{
			for (int j = 0; j < lowlong; j++) //start from 0 degrees, and end one above the lowest longitude
				if (spheremesh[i][j][0] >= xlow && spheremesh[i][j][0] <= xup && spheremesh[i][j][1] >= ylow && spheremesh[i][j][1] <= yup && spheremesh[i][j][2] >= zlow && spheremesh[i][j][2] <= zup) //check if the point is inside the box. If yes, add the area
					area += spheremesh[i][j][3];
			for (int j = uplong; j < latitudes; j++)
				if (spheremesh[i][j][0] >= xlow && spheremesh[i][j][0] <= xup && spheremesh[i][j][1] >= ylow && spheremesh[i][j][1] <= yup && spheremesh[i][j][2] >= zlow && spheremesh[i][j][2] <= zup) //check if the point is inside the box. If yes, add the area
					area += spheremesh[i][j][3];
		}
	}
	return area;
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

//vector vfirst(vector* velocities)
//{
//	vector fvelocities[sizeof(velocities)]; //0, 0, 0 -> 2*cube-1, 2*cube-1, 2*cube-1 base 2*cube, 4*cube*cube, 2*cube, 1
//	for (int i = 0; i <= sizeof(fvelocities); ++i)
//	{
//		fvelocities[i] =
//	}
//	
//}

static vector forcel(const double radius, double xlow, double xup, double ylow, double yup, double zlow, double zup, unsigned int latitudes, unsigned int longitudes, double*** spheremesh,double p, double sx, double sy, double sz,tensor***st, int cube)
{
	xlow -= sx;
	xup -= sx;
	ylow -= sy;
	yup -= sy;
	zlow -= sz;
	zup -= sz;
	//pt. 3 find the force
	// maybe restrict the range of indices that can potentially intersect with the cube to speed up area calculation.
	vector force(0,0,0);
	unsigned int lowlat = latitudes;
	if (zlow > -radius)
		lowlat = (int)ceil((acos(zlow / radius) - M_PI / (2 * latitudes)) / (M_PI / (latitudes)));
	unsigned int uplat = 0;
	if (zup < radius)
		uplat = (int)floor((acos(zup / radius) - M_PI / (2 * latitudes)) / (M_PI / (latitudes)));
	if (zup <-radius || zlow >radius)
	{
		return force;
	}
	unsigned int lowlong = 0;
	unsigned int uplong = longitudes;
	double angles[4] = { angleof(xlow, ylow), angleof(xlow, yup), angleof(xup, ylow), angleof(xup, yup) };
	std::sort(angles, angles + 4);
	if (!(xlow <= 0 && ylow <= 0 && xup >= 0 && yup >= 0)) //find the longitude limits, assuming the projection of the box to the x-y plane does not contain the origin
	{
		lowlong = (int)floor(angles[0] / (2 * M_PI / longitudes));
		uplong = (int)ceil(angles[3] / (2 * M_PI / longitudes));
		if (uplong > longitudes)
			uplong = longitudes;
	}
	if (!(xlow > 0 && ylow < 0 && yup >= 0)) //case 1: the polar axis does not intersect the projection to the x-y plane, or the lower y-bound is the polar axis
	{
		for (unsigned int i = uplat; i < lowlat; i++)   //for loop starts at latitude corresponding to zup -> latitude corresponding to zlow, since the latitudes start at the top of the sphere
		{
			for (unsigned int j = lowlong; j < uplong; j++) //start at the lowest longitude and go to the highest longitude
				if (spheremesh[i][j][0] >= xlow && spheremesh[i][j][0] <= xup && spheremesh[i][j][1] >= ylow && spheremesh[i][j][1] <= yup && spheremesh[i][j][2] >= zlow && spheremesh[i][j][2] <= zup) //check if the point is inside the box. If yes, add the area
				{
					vector r(spheremesh[i][j][0], spheremesh[i][j][1], spheremesh[i][j][2]);
					force =force+ (st[(int)round((xlow+sx) / (xup-xlow) + cube)][(int)round((ylow+sy) / (xup - xlow) + cube)][(int)round((zlow+sz) / (xup - xlow) + cube)]) * r * spheremesh[i][j][3] / radius;
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
					force = force + (st[(int)round((xlow + sx) / (xup - xlow) + cube)][(int)round((ylow + sy) / (xup - xlow) + cube)][(int)round((zlow + sz) / (xup - xlow) + cube)]) * r * spheremesh[i][j][3] / radius;
				}
			for (unsigned int j = uplong; j < latitudes; j++)
				if (spheremesh[i][j][0] >= xlow && spheremesh[i][j][0] <= xup && spheremesh[i][j][1] >= ylow && spheremesh[i][j][1] <= yup && spheremesh[i][j][2] >= zlow && spheremesh[i][j][2] <= zup) //check if the point is inside the box. If yes, add the area
				{
					vector r(spheremesh[i][j][0], spheremesh[i][j][1], spheremesh[i][j][2]);
					force = force + (st[(int)round((xlow + sx) / (xup - xlow) + cube)][(int)round((ylow + sy) / (xup - xlow) + cube)][(int)round((zlow + sz) / (xup - xlow) + cube)]) * r * spheremesh[i][j][3] / radius;
				}
		}
	}
	/*std::cout << force.X() << " " << force.Y() << " " << force.Z() << "\n";*/
	//std::cout << "Velocity:" << velocities[0][x][y - 1][z] << " " << velocities[1][x][y - 1][z] << " " << velocities[2][x][y - 1][z] << " " << velocities[0][x][y][z] << " " << velocities[1][x][y][z] << " " << velocities[2][x][y][z] << " " << velocities[0][x][y + 1][z] << " " << velocities[1][x][y + 1][z] << " " << velocities[2][x][y + 1][z] << "\n";
	/*std::cout << "Stress:" << st[(int)round((xlow+sx) / (xup-xlow) + cube)][(int)round((ylow+sy) / (xup - xlow) + cube)][(int)round((zlow+sz) / (xup - xlow) + cube)].getValue(1, 1) << " " << st[(int)round((xlow+sx) / (xup-xlow) + cube)][(int)round((ylow+sy) / (xup - xlow) + cube)][(int)round((zlow+sz) / (xup - xlow) + cube)].getValue(1, 2) << " " << st[(int)round((xlow+sx) / (xup-xlow) + cube)][(int)round((ylow+sy) / (xup - xlow) + cube)][(int)round((zlow+sz) / (xup - xlow) + cube)].getValue(1, 3) << "\n";
	std::cout << "Stress:" << st[(int)round((xlow+sx) / (xup-xlow) + cube)][(int)round((ylow+sy) / (xup - xlow) + cube)][(int)round((zlow+sz) / (xup - xlow) + cube)].getValue(2, 1) << " " << st[(int)round((xlow+sx) / (xup-xlow) + cube)][(int)round((ylow+sy) / (xup - xlow) + cube)][(int)round((zlow+sz) / (xup - xlow) + cube)].getValue(2, 2) << " " << st[(int)round((xlow+sx) / (xup-xlow) + cube)][(int)round((ylow+sy) / (xup - xlow) + cube)][(int)round((zlow+sz) / (xup - xlow) + cube)].getValue(2, 3) << "\n";
	std::cout << "Stress:" << st[(int)round((xlow+sx) / (xup-xlow) + cube)][(int)round((ylow+sy) / (xup - xlow) + cube)][(int)round((zlow+sz) / (xup - xlow) + cube)].getValue(3, 1) << " " << st[(int)round((xlow+sx) / (xup-xlow) + cube)][(int)round((ylow+sy) / (xup - xlow) + cube)][(int)round((zlow+sz) / (xup - xlow) + cube)].getValue(3, 2) << " " << st[(int)round((xlow+sx) / (xup-xlow) + cube)][(int)round((ylow+sy) / (xup - xlow) + cube)][(int)round((zlow+sz) / (xup - xlow) + cube)].getValue(3, 3) << "\n"; */
	std::cout << (xlow) << " " << xup << "\n";
	//std::cout << stresses[x][y][z].getValue(2, 1) << " " << stresses[x][y][z].getValue(2, 2) << " " << stresses[x][y][z].getValue(2, 3) << "\n";
	//std::cout << stresses[x][y][z].getValue(3, 1) << " " << stresses[x][y][z].getValue(3, 2) << " " << stresses[x][y][z].getValue(3, 3) << "\n";
	//std::cout << "Position:" << cubes[i][0] << " " << cubes[i][1] << " " << cubes[i][2] << "\n";
	//std::cout << "Increment:" << s / cube << "\n";

	return force;
}
static vector torquel(const double radius, const double xlow, const double xup, const double ylow, const double yup, const double zlow, const double zup, unsigned int latitudes, unsigned int longitudes, double*** spheremesh, double p,double sx, double sy, double sz, tensor*** st, int cube)
{
	//pt. 4 find the torque
	// maybe restrict the range of indices that can potentially intersect with the cube to speed up area calculation.
	vector torque(0, 0, 0);
	unsigned int lowlat = latitudes;
	if (zlow > -radius)
		lowlat = (int)ceil((acos(zlow / radius) - M_PI / (2 * latitudes)) / (M_PI / (latitudes)));
	unsigned int uplat = 0;
	if (zup < radius)
		uplat = (int)floor((acos(zup / radius) - M_PI / (2 * latitudes)) / (M_PI / (latitudes)));
	if (zup <-radius || zlow >radius)
	{
		return torque;
	}
	unsigned int lowlong = 0;
	unsigned int uplong = longitudes;
	double angles[4] = { angleof(xlow, ylow), angleof(xlow, yup), angleof(xup, ylow), angleof(xup, yup) };
	std::sort(angles, angles + 4);
	if (!(xlow <= 0 && ylow <= 0 && xup >= 0 && yup >= 0)) //find the longitude limits, assuming the projection of the box to the x-y plane does not contain the origin
	{
		lowlong = (int)floor(angles[0] / (2 * M_PI / longitudes));
		uplong = (int)ceil(angles[3] / (2 * M_PI / longitudes));
		if (uplong > longitudes)
			uplong = longitudes;
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
					torque = torque + r%((st[(int)round((xlow) / (xup - xlow) + cube)][(int)round((ylow) / (xup - xlow) + cube)][(int)round((zlow) / (xup - xlow) + cube)]) *n)* spheremesh[i][j][3];
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
					torque = torque + r % ((st[(int)round((xlow) / (xup - xlow) + cube)][(int)round((ylow) / (xup - xlow) + cube)][(int)round((zlow) / (xup - xlow) + cube)]) * n) * spheremesh[i][j][3];
				}
			for (unsigned int j = uplong; j < latitudes; j++)
				if (spheremesh[i][j][0] >= xlow && spheremesh[i][j][0] <= xup && spheremesh[i][j][1] >= ylow && spheremesh[i][j][1] <= yup && spheremesh[i][j][2] >= zlow && spheremesh[i][j][2] <= zup) //check if the point is inside the box. If yes, add the area
				{
					vector r(spheremesh[i][j][0] + sx, spheremesh[i][j][1] + sy, spheremesh[i][j][2] + sz);
					vector n(spheremesh[i][j][0], spheremesh[i][j][1], spheremesh[i][j][2]);
					n = n / radius;
					torque = torque + r % ((st[(int)round((xlow) / (xup - xlow) + cube)][(int)round((ylow) / (xup - xlow) + cube)][(int)round((zlow) / (xup - xlow) + cube)]) * n) * spheremesh[i][j][3];
				}
		}
	}
	return torque;
}

int checkcube(double radius, double s, double x, double y, double z, double sx, double sy, double sz)
{
	//center of cube
	vector center(x-sx, y-sy, z-sz);
	if (center.length() + s<radius || center.length() - s >radius)
	{
		return -1;
	}
	return 0;
}
tensor gradv(double x, double y, double z, int n, double h, double**** a)
{

	double arr[9] = { partdif(a, x, y, z, n, h, 0, 0), partdif(a, x, y, z, n, h, 1, 0), partdif(a, x, y, z, n, h, 2, 0),
	partdif(a, x, y, z, n, h, 0, 1), partdif(a, x, y, z, n, h, 1, 1), partdif(a, x, y, z, n, h, 2, 1),
	partdif(a, x, y, z, n, h, 0, 2), partdif(a, x, y, z, n, h, 1, 2), partdif(a, x, y, z, n, h, 2, 2) };
	tensor grad(arr);
	/*std::cout <<"\n" << grad.getValue(1, 1) << grad.getValue(1, 2) << grad.getValue(1, 3) << "\n";
	std::cout << grad.getValue(2, 1) << grad.getValue(2, 2) << grad.getValue(2, 3) << "\n";
	std::cout << grad.getValue(3, 1) << grad.getValue(3, 2) << grad.getValue(3, 3) << "\n";*/
	return grad;
}



vector force(double x0, double y0, double z0, unsigned int longitudes, unsigned int latitudes, double radius, double cl,int cube, double p, double **cubes, const std::vector<int>& vect, double *** spheremesh,tensor***st)
{
	double s = cl / cube;
	double h = s / 2;
	vector sum(0,0,0);
	/*auto start = std::chrono::high_resolution_clock::now();
	auto finish = std::chrono::high_resolution_clock::now();
	auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
	std::cout << "time of sphere: " << milliseconds.count() / 1000.0 << "s\n";
	start = std::chrono::high_resolution_clock::now();*/
	int l = 0;
	assert(spheremesh != nullptr);
	for (auto &kk :vect)
	{
			double a = cubes[kk][0] - h, b = cubes[kk][0] + h, c = cubes[kk][1] - h, d = cubes[kk][1] + h, e = cubes[kk][2] - h, f = cubes[kk][2] + h;
			sum = sum+ forcel(radius, a, b, c, d, e, f, latitudes, longitudes, spheremesh,p,x0,y0,z0,st,cube);
	}
	/*finish = std::chrono::high_resolution_clock::now();
	milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
	std::cout << "time of area: " << milliseconds.count() / 1000.0 << "s\n";
	std::cout << "Expected area: " << 4 * M_PI * radius * radius << std::endl;
	std::cout << "Calculated force: " << "(" << sum.X() << ", " << sum.Y() << ", " << sum.Z() << ")" << std::endl;*/
	return sum;
}
vector torque(double x0, double y0, double z0, unsigned int longitudes, unsigned int latitudes, double radius, double cl, int cube, double p,double **cubes, const std::vector<int>& vect, double *** spheremesh,tensor***st)
{
	double s = cl / cube;
	double h = s / 2;
	vector sum(0, 0, 0);
	/*auto start = std::chrono::high_resolution_clock::now();
	auto finish = std::chrono::high_resolution_clock::now();
	auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
	std::cout << "time of sphere: " << milliseconds.count() / 1000.0 << "s\n";
	start = std::chrono::high_resolution_clock::now();*/
	int l = 0;
	assert(spheremesh != nullptr);
	for (auto &kk :vect)
	{
			double a = cubes[kk][0] - h, b = cubes[kk][0] + h, c = cubes[kk][1] - h, d = cubes[kk][1] + h, e = cubes[kk][2] - h, f = cubes[kk][2] + h;
			sum = sum + torquel(radius, a, b, c, d, e, f, latitudes, longitudes, spheremesh, p, x0, y0, z0,st,cube);
	}
	/*finish = std::chrono::high_resolution_clock::now();
	milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
	std::cout << "time of area: " << milliseconds.count() / 1000.0 << "s\n";
	std::cout << "Expected area: " << 4 * M_PI * radius * radius << std::endl;
	std::cout << "Calculated torque: " << "(" << sum.X() << ", " << sum.Y() << ", " << sum.Z() << ")" << std::endl;
	*/return sum;
}
int bleargh()
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
	double radius,hlength;
	std::cin >> radius;
	std::cout << "What is the half the entire length of the grid?";
	std::cin >> hlength;
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
	double** cubes = nullptr;
	cubes = compute_cubes(hlength, cube, cubes);
  double*** spheremesh = compute_spheremesh(radius, latitudes, longitudes);

  double **** velocities = allocate4d(4, 2*cube, 2*cube, 2*cube);//and pressure
	int l = sizeof(cubes);
	for (int x = 0; x < 2*cube; x++)
		for (int y = 0; y < 2*cube; y++)
			for (int z = 0; z < 2*cube; z++)
			{
        // calculate i from x,y,z
        int i = x+ y * (2 * cube) + z* (2 * cube)* (2 * cube);
				velocities[0][x][y][z] = velocity(cubes[i][0], cubes[i][1], cubes[i][2]).X();
				velocities[1][x][y][z] = velocity(cubes[i][0], cubes[i][1], cubes[i][2]).Y();
				velocities[2][x][y][z] = velocity(cubes[i][0], cubes[i][1], cubes[i][2]).Z();
				velocities[3][x][y][z] = pressure;
			}
  writeVTKFile(0, velocities, 2*cube, 2*cube, 2*cube, hlength/cube, hlength/cube, hlength/cube);
	tensor *** stresses = allocate3dt(2*cube, 2*cube, 2*cube);
	for (int x = 0; x < 2*cube; x++)
		for (int y = 0; y < 2*cube; y++)
			for (int z = 0; z < 2*cube; z++)
			{
				// calculate i from x,y,z
				int i = x + y * (2 * cube) + z * (2 * cube) * (2 * cube);
				stresses[x][y][z] = stress(pressure,cubes[i][0], cubes[i][1], cubes[i][2]);
			}
	assert(cubes != nullptr);
	vector ucm(0, 0, 0);
	vector rcm(spherex, spherey, spherez);
	vector angv(0, 0, 0);
	std::vector <int> correctindexes;
	for (unsigned long int kk = 0; kk < 8 * cube * cube * cube; kk++)
	{
		if (checkcube(radius, hlength/cube, cubes[kk][0], cubes[kk][1], cubes[kk][2],0,0,0) == 0)
		{
			correctindexes.push_back(kk);
		}
	}
	for (int i = 0; i < 100; i++)
	{
		vector F = force(rcm.X(), rcm.Y(), rcm.Z(), longitudes, latitudes, radius, hlength ,cube, pressure, cubes, correctindexes, spheremesh,stresses);
		angv = angv+torque(rcm.X(), rcm.Y(), rcm.Z(), longitudes, latitudes, radius, hlength,cube, pressure,cubes, correctindexes, spheremesh,stresses)*increment/(2*mass*radius*radius/5);
		rcm = rcm + ucm * increment + F * (increment) * (increment) / (2 * mass);
		Eigen::VectorXd b(8*cube*cube*cube);
		Eigen::MatrixXd A(8*cube*cube*cube,8*cube*cube*cube);
		int c=0;
		for (int x = 0; x < 2 * cube; x++)
			for (int y = 0; y < 2 * cube; y++)
				for (int z = 0; z < 2 * cube; z++)
				{
					b(c)=divvadvecv(velocities,x,y,z,8*cube*cube*cube, hlength/cube);
					++c;
				}
		vector *** stressdiv= divtenall(stresses, 2 * cube, hlength / cube); //cjamge radois to cube lengths
		for (int x = 0; x < 2 * cube; x++)
			for (int y = 0; y < 2 * cube; y++)
				for (int z = 0; z < 2 * cube; z++)
				{
					// calculate i from x,y,z
					int i = x + y * (2 * cube) + z * (2 * cube) * (2 * cube);
					vector vel(velocities[0][x][y][z], velocities[1][x][y][z], velocities[2][x][y][z]);
					tensor grad = gradv(x, y, z, 2 * cube,hlength / cube, velocities);
					vector secterm = grad.transpose() * vel;
					vector change = stressdiv[x][y][z] - secterm;
					velocities[0][x][y][z] += increment * change.X();
					velocities[1][x][y][z] += increment * change.Y();
					velocities[2][x][y][z] += increment * change.Z();
					double a[] = { pressure, 0, 0, 0, pressure, 0, 0, 0, pressure };
					tensor p(a);
					stresses[x][y][z] = grad.transpose()+grad;
				}
		writeVTKFile(i+1, velocities, 2*cube, 2*cube, 2*cube, hlength/cube, hlength/cube, hlength/cube);
	}
	std::cout << rcm.X() <<" " << rcm.Y()<<" " << rcm.Z();
	//vodt=divstress-gradv*v
	
	
	return 0;
}
int main()
{
  //bleargh();
  //exit(0);
	int latitudes=70, longitudes=70;
	double cube=97;
	double** cubes = nullptr;
	double s = 0.5;
	print_rss_memory("start program");
	cubes = compute_cubes(s, cube, cubes);
	print_rss_memory("after compute cubes");
	double radius = 0.01;
	double p =0.1,pressure=1;
	
	double**** velocities = allocate4d(4, 2 * cube, 2 * cube, 2 * cube);//and pressure
	print_rss_memory("after velocities");
	int l = sizeof(cubes);
	for (int x = 0; x < 2 * cube; x++)
		for (int y = 0; y < 2 * cube; y++)
			for (int z = 0; z < 2 * cube; z++)
			{
				// calculate i from x,y,z
				int i = x + y * (2 * cube) + z * (2 * cube) * (2 * cube);
				velocities[0][x][y][z] = velocity(cubes[i][0], cubes[i][1], cubes[i][2]).X();
				velocities[1][x][y][z] = velocity(cubes[i][0], cubes[i][1], cubes[i][2]).Y();
				velocities[2][x][y][z] = velocity(cubes[i][0], cubes[i][1], cubes[i][2]).Z();
				velocities[3][x][y][z] = pressure;
			}
	double*** spheremesh = compute_spheremesh(radius, latitudes, longitudes);
	print_rss_memory("after spheremesh");
	tensor*** stresses = allocate3dt(2 * cube, 2 * cube, 2 * cube);
	print_rss_memory("after stresses");
	for (int x = 0; x < 2 * cube; x++)
		for (int y = 0; y < 2 * cube; y++)
			for (int z = 0; z < 2 * cube; z++)
			{
				//// calculate i from x,y,z
				int i = x + y * (2 * cube) + z * (2 * cube) * (2 * cube);
				
				tensor grad = gradv(x, y, z, 2 * cube, s / cube, velocities);
				double a[] = { pressure, 0, 0, 0, pressure, 0, 0, 0, pressure };
				tensor p(a);
				stresses[x][y][z] = grad.transpose() + grad - p;
				if (x == cube +50&& y == cube +27&& z == cube+28)
				{
					std::cout << "Velocity:" << velocities[0][x][y - 1][z] << " " << velocities[1][x][y - 1][z] << " " << velocities[2][x][y - 1][z] << " " << velocities[0][x][y][z] << " " << velocities[1][x][y][z] << " " << velocities[2][x][y][z] << " " << velocities[0][x][y + 1][z] << " " << velocities[1][x][y + 1][z] << " " << velocities[2][x][y + 1][z] << "\n";
					std::cout <<"Stress:" << stresses[x][y][z].getValue(1, 1) << " " << stresses[x][y][z].getValue(1, 2) <<" " << stresses[x][y][z].getValue(1, 3) << "\n";
					std::cout << stresses[x][y][z].getValue(2, 1) <<" " << stresses[x][y][z].getValue(2, 2) <<" " << stresses[x][y][z].getValue(2, 3) << "\n";
					std::cout << stresses[x][y][z].getValue(3, 1) << " " << stresses[x][y][z].getValue(3, 2) <<" " << stresses[x][y][z].getValue(3, 3) << "\n";
					std::cout << "Position:" << cubes[i][0] <<" " << cubes[i][1] << " " << cubes[i][2]<< "\n";
					std::cout << "Increment:" << s / cube << "\n";
				}
				
				
			}
	std::cout << "radius forcex forcey forcez torquex torquey torquez error1 error2\n";
	//
	/*std::cout << "N area error\n";*/
	for (int i = 1; i <= 1; i++)
	{
		std::vector <int> correctindexes;
		for (unsigned long int kk = 0; kk < 8 * cube * cube * cube; kk++)
		{
			if (checkcube(radius, s / cube, cubes[kk][0], cubes[kk][1], cubes[kk][2],0,p,0) == 0)
			{
				correctindexes.push_back(kk);
			}
		}
		vector sum(0,0,0);
		vector sumt(0, 0, 0);
		/*for (unsigned long int kk = 0; kk < 8 * cube * cube * cube; kk++)
		{
			if (checkcube(radius, 50 / cube, cubes[kk][0], cubes[kk][1], cubes[kk][2]) == 0)
			{
				correctindexes.push_back(kk);
			}
		}*/


			sumt= torque(0, p, 0, longitudes, latitudes, radius,s, cube, pressure, cubes, correctindexes, spheremesh,stresses);
			sum = force(0, p, 0, longitudes, latitudes, radius, s,cube, pressure, cubes, correctindexes, spheremesh,stresses);
			correctindexes.clear();


		std::cout << radius<<" " << sum.X() << " "<< sum.Y() << " " <<sum.Z() << "  " << sumt.X()<< " " << sumt.Y() << " " << sumt.Z() << " "
			<< abs(1 - sum.X() / (8 * M_PI * radius*radius*radius * 4 / 3)) * 100 << " " << abs(1 - sumt.Z() / (-8 * M_PI*p * radius*radius*radius * 4 / 3)) * 100 << "\n";
		radius += 0.01;
		
		/*double sum = 0;


		for (int kk:correctindexes)
		{
			sum+=areal(radius, cubes[kk][0]-s/cube/2, cubes[kk][0]+ s / cube / 2, cubes[kk][1]- s / cube / 2, cubes[kk][1]+ s / cube / 2, cubes[kk][2]- s / cube / 2, cubes[kk][2]+ s / cube / 2, latitudes, longitudes, spheremesh);
		}
		std::cout << radius << " " << sum << " " << 100 * abs(4 * M_PI*radius*radius - sum) / 4 / M_PI/radius/radius << "\n";*/
		free(spheremesh);
		/*radius += 0.01;*/
		spheremesh = compute_spheremesh(radius, latitudes, longitudes);
		print_rss_memory("after recompute sphere");
	}

}
