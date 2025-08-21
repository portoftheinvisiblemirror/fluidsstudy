#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <vector>
#include <fstream>
#include <stdlib.h>

#include "utilities.hpp"
#include "mem.hpp"
#include "vector.hpp"
#include "tensor.hpp"
#include "vtk.hpp"
#include <tuple>

double*** compute_spheremesh(const double radius, const unsigned int latitudes, const unsigned int longitudes, double*** spheremesh);

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
	double*** spheremesh = allocate3d(latitudes, longitudes, 4);
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
static double areal(double radius, double xlow, double xup, double ylow, double yup, double zlow, double zup, unsigned int latitudes, unsigned int longitudes, double*** spheremesh, double sx, double sy, double sz, tensor*** st, const double dx, const double dy, const double dz, int x, int y, int z)
{
	xlow -= sx;
	xup -= sx;
	ylow -= sy;
	yup -= sy;
	zlow -= sz;
	zup -= sz;
	//pt. 4 find the torque
	// maybe restrict the range of indices that can potentially intersect with the cube to speed up area calculation.
	double area=0;
	unsigned int lowlat = latitudes;
	if (zlow > -radius)
		lowlat = (int)ceil((acos(zlow / radius) - M_PI / (2 * latitudes)) / (M_PI / (latitudes)));
	unsigned int uplat = 0;
	if (zup < radius)
		uplat = (int)floor((acos(zup / radius) - M_PI / (2 * latitudes)) / (M_PI / (latitudes)));
	if (zup <-radius || zlow >radius)
	{
		return area;
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
					area +=spheremesh[i][j][3];
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
					area += spheremesh[i][j][3];
				}
			for (unsigned int j = uplong; j < latitudes; j++)
				if (spheremesh[i][j][0] >= xlow && spheremesh[i][j][0] <= xup && spheremesh[i][j][1] >= ylow && spheremesh[i][j][1] <= yup && spheremesh[i][j][2] >= zlow && spheremesh[i][j][2] <= zup) //check if the point is inside the box. If yes, add the area
				{
					area += spheremesh[i][j][3];
				}
		}
	}
	return area;
}

static vector forcel(double radius, double xlow, double xup, double ylow, double yup, double zlow, double zup, unsigned int latitudes, unsigned int longitudes, double*** spheremesh, double sx, double sy, double sz, tensor*** st, const double dx, const double dy, const double dz, int x, int y, int z)
{
	xlow -= sx;
	xup -= sx;
	ylow -= sy;
	yup -= sy;
	zlow -= sz;
	zup -= sz;
	//pt. 3 find the force
	// maybe restrict the range of indices that can potentially intersect with the cube to speed up area calculation.
	vector force(0, 0, 0);
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
					force = force + (st[x][y][z]) * r * spheremesh[i][j][3] / radius;
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
					force = force + (st[x][y][z]) * r * spheremesh[i][j][3] / radius;
				}
			for (unsigned int j = uplong; j < latitudes; j++)
				if (spheremesh[i][j][0] >= xlow && spheremesh[i][j][0] <= xup && spheremesh[i][j][1] >= ylow && spheremesh[i][j][1] <= yup && spheremesh[i][j][2] >= zlow && spheremesh[i][j][2] <= zup) //check if the point is inside the box. If yes, add the area
				{
					vector r(spheremesh[i][j][0], spheremesh[i][j][1], spheremesh[i][j][2]);
					force = force + (st[x][y][z]) * r * spheremesh[i][j][3] / radius;
				}
		}
	}
	/*std::cout << force.X() << " " << force.Y() << " " << force.Z() << "\n";*/
	//std::cout << "Velocity:" << velocities[0][x][y - 1][z] << " " << velocities[1][x][y - 1][z] << " " << velocities[2][x][y - 1][z] << " " << velocities[0][x][y][z] << " " << velocities[1][x][y][z] << " " << velocities[2][x][y][z] << " " << velocities[0][x][y + 1][z] << " " << velocities[1][x][y + 1][z] << " " << velocities[2][x][y + 1][z] << "\n";
	/*std::cout << "Stress:" << st[(int)round((xlow+sx) / (xup-xlow) + cube)][(int)round((ylow+sy) / (xup - xlow) + cube)][(int)round((zlow+sz) / (xup - xlow) + cube)].getValue(1, 1) << " " << st[(int)round((xlow+sx) / (xup-xlow) + cube)][(int)round((ylow+sy) / (xup - xlow) + cube)][(int)round((zlow+sz) / (xup - xlow) + cube)].getValue(1, 2) << " " << st[(int)round((xlow+sx) / (xup-xlow) + cube)][(int)round((ylow+sy) / (xup - xlow) + cube)][(int)round((zlow+sz) / (xup - xlow) + cube)].getValue(1, 3) << "\n";
	std::cout << "Stress:" << st[(int)round((xlow+sx) / (xup-xlow) + cube)][(int)round((ylow+sy) / (xup - xlow) + cube)][(int)round((zlow+sz) / (xup - xlow) + cube)].getValue(2, 1) << " " << st[(int)round((xlow+sx) / (xup-xlow) + cube)][(int)round((ylow+sy) / (xup - xlow) + cube)][(int)round((zlow+sz) / (xup - xlow) + cube)].getValue(2, 2) << " " << st[(int)round((xlow+sx) / (xup-xlow) + cube)][(int)round((ylow+sy) / (xup - xlow) + cube)][(int)round((zlow+sz) / (xup - xlow) + cube)].getValue(2, 3) << "\n";
	std::cout << "Stress:" << st[(int)round((xlow+sx) / (xup-xlow) + cube)][(int)round((ylow+sy) / (xup - xlow) + cube)][(int)round((zlow+sz) / (xup - xlow) + cube)].getValue(3, 1) << " " << st[(int)round((xlow+sx) / (xup-xlow) + cube)][(int)round((ylow+sy) / (xup - xlow) + cube)][(int)round((zlow+sz) / (xup - xlow) + cube)].getValue(3, 2) << " " << st[(int)round((xlow+sx) / (xup-xlow) + cube)][(int)round((ylow+sy) / (xup - xlow) + cube)][(int)round((zlow+sz) / (xup - xlow) + cube)].getValue(3, 3) << "\n"; */
	//std::cout << stresses[x][y][z].getValue(2, 1) << " " << stresses[x][y][z].getValue(2, 2) << " " << stresses[x][y][z].getValue(2, 3) << "\n";
	//std::cout << stresses[x][y][z].getValue(3, 1) << " " << stresses[x][y][z].getValue(3, 2) << " " << stresses[x][y][z].getValue(3, 3) << "\n";
	//std::cout << "Position:" << cubes[i][0] << " " << cubes[i][1] << " " << cubes[i][2] << "\n";
	//std::cout << "Increment:" << s / cube << "\n";

	return force;
}
static vector torquel(double radius, double xlow, double xup, double ylow, double yup, double zlow, double zup, unsigned int latitudes, unsigned int longitudes, double*** spheremesh, double sx, double sy, double sz, tensor*** st, const double dx, const double dy, const double dz, int x, int y, int z)
{
	xlow -= sx;
	xup -= sx;
	ylow -= sy;
	yup -= sy;
	zlow -= sz;
	zup -= sz;
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
					torque = torque + r % ((st[x][y][z]) * n) * spheremesh[i][j][3];
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
					torque = torque + r % ((st[x][y][z]) * n) * spheremesh[i][j][3];
				}
			for (unsigned int j = uplong; j < latitudes; j++)
				if (spheremesh[i][j][0] >= xlow && spheremesh[i][j][0] <= xup && spheremesh[i][j][1] >= ylow && spheremesh[i][j][1] <= yup && spheremesh[i][j][2] >= zlow && spheremesh[i][j][2] <= zup) //check if the point is inside the box. If yes, add the area
				{
					vector r(spheremesh[i][j][0] + sx, spheremesh[i][j][1] + sy, spheremesh[i][j][2] + sz);
					vector n(spheremesh[i][j][0], spheremesh[i][j][1], spheremesh[i][j][2]);
					n = n / radius;
					torque = torque + r % ((st[x][y][z]) * n) * spheremesh[i][j][3];
				}
		}
	}
	return torque;
}

tensor gradv(double x, double y, double z, int nx,int ny, int nz, double dx, double dy, double dz, std::vector<std::vector<std::vector<std::vector<double>>>> a)
{

	double arr[9] = { partdif(a, x, y, z, nx, dx, 0, 0), partdif(a, x, y, z, nx, dx, 1, 0), partdif(a, x, y, z, nx, dx, 2, 0),
	partdif(a, x, y, z, ny, dy, 0, 1), partdif(a, x, y, z, ny, dy, 1, 1), partdif(a, x, y, z, ny, dy, 2, 1),
	partdif(a, x, y, z, nz, dz, 0, 2), partdif(a, x, y, z, nz, dz, 1, 2), partdif(a, x, y, z, nz, dz, 2, 2) };
	tensor grad(arr);
	/*std::cout <<"\n" << grad.getValue(1, 1) << grad.getValue(1, 2) << grad.getValue(1, 3) << "\n";
	std::cout << grad.getValue(2, 1) << grad.getValue(2, 2) << grad.getValue(2, 3) << "\n";
	std::cout << grad.getValue(3, 1) << grad.getValue(3, 2) << grad.getValue(3, 3) << "\n";*/
	return grad;
}



vector force(const double x0, const double y0, const double z0, const double dx, const double dy, const double dz, const double nx, const int ny, const int nz, unsigned int longitudes, unsigned int latitudes, double radius, double*** spheremesh, std::vector<std::vector<std::vector<bool>>> sphere, tensor*** st, const double R)
{
	vector sum(0, 0, 0);
	/*auto start = std::chrono::high_resolution_clock::now();
	auto finish = std::chrono::high_resolution_clock::now();
	auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
	std::cout << "time of sphere: " << milliseconds.count() / 1000.0 << "s\n";
	start = std::chrono::high_resolution_clock::now();*/
	int l = 0;
	assert(spheremesh != nullptr);
	for (int i = 0; i < nx; ++i)
		for (int j = 0; j < ny; ++j)
			for (int k = 0; k < nx; ++k)
			{
				if (sphere[i][j][k])
				{
					double a = i * dx, b = i * (1 + dx), c = j * dy - R, d = (j + 1) * dy - R, e = k * dz - R, f = (k + 1) * dz - R;
					sum = sum + forcel(radius, a, b, c, d, e, f, latitudes, longitudes, spheremesh, x0, y0, z0, st, dx, dy, dz, i, j, k);
				}
			}
	/*finish = std::chrono::high_resolution_clock::now();
	milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
	std::cout << "time of area: " << milliseconds.count() / 1000.0 << "s\n";
	std::cout << "Expected area: " << 4 * M_PI * radius * radius << std::endl;
	std::cout << "Calculated force: " << "(" << sum.X() << ", " << sum.Y() << ", " << sum.Z() << ")" << std::endl;*/
	return sum;
}
vector torque(const double x0, const double y0, const double z0, const double dx, const double dy, const double dz, const double nx, const int ny, const int nz, unsigned int longitudes, unsigned int latitudes, double radius, double*** spheremesh, std::vector<std::vector<std::vector<bool>>> sphere, tensor*** st, const double R)
{
	vector sum(0, 0, 0);
	/*auto start = std::chrono::high_resolution_clock::now();
	auto finish = std::chrono::high_resolution_clock::now();
	auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
	std::cout << "time of sphere: " << milliseconds.count() / 1000.0 << "s\n";
	start = std::chrono::high_resolution_clock::now();*/
	int l = 0;
	assert(spheremesh != nullptr);
	for (int i = 0; i<nx;++i)
		for (int j = 0; j < ny; ++j)
			for (int k = 0; k < nx; ++k)
	{
				if (sphere[i][j][k])
				{
					double a = i * dx, b = i * (1 + dx), c = j * dy - R, d = (j + 1) * dy - R, e = k * dz - R, f = (k + 1) * dz - R;
					sum = sum + torquel(radius, a, b, c, d, e, f, latitudes, longitudes, spheremesh, x0, y0, z0, st, dx, dy, dz, i, j, k);
				}
	}
	/*finish = std::chrono::high_resolution_clock::now();
	milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
	std::cout << "time of area: " << milliseconds.count() / 1000.0 << "s\n";
	std::cout << "Expected area: " << 4 * M_PI * radius * radius << std::endl;
	std::cout << "Calculated torque: " << "(" << sum.X() << ", " << sum.Y() << ", " << sum.Z() << ")" << std::endl;
	*/return sum;
}
double area(const double x0, const double y0, const double z0, const double dx, const double dy, const double dz, const double nx, const int ny, const int nz, unsigned int longitudes, unsigned int latitudes, double radius, double*** spheremesh, std::vector<std::vector<std::vector<bool>>> sphere, tensor*** st, const double R)
{
	double sum=0;
	/*auto start = std::chrono::high_resolution_clock::now();
	auto finish = std::chrono::high_resolution_clock::now();
	auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
	std::cout << "time of sphere: " << milliseconds.count() / 1000.0 << "s\n";
	start = std::chrono::high_resolution_clock::now();*/
	int l = 0;
	assert(spheremesh != nullptr);
	for (int i = 0; i < nx; ++i)
		for (int j = 0; j < ny; ++j)
			for (int k = 0; k < nx; ++k)
			{
				if (sphere[i][j][k])
				{
					double a = i * dx, b = i * (1 + dx), c = j * dy - R, d = (j + 1) * dy - R, e = k * dz - R, f = (k + 1) * dz - R;
					sum = sum + areal(radius, a, b, c, d, e, f, latitudes, longitudes, spheremesh, x0, y0, z0, st, dx, dy, dz, i, j, k);
					std::cout << "Sum: " << sum << "\n";
				}
			}
	/*finish = std::chrono::high_resolution_clock::now();
	milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
	std::cout << "time of area: " << milliseconds.count() / 1000.0 << "s\n";
	std::cout << "Expected area: " << 4 * M_PI * radius * radius << std::endl;
	std::cout << "Calculated force: " << "(" << sum.X() << ", " << sum.Y() << ", " << sum.Z() << ")" << std::endl;*/
	return sum;
}
vector *forceandtorque(std::vector<std::vector<std::vector<bool>>> sphere, std::vector<std::vector<std::vector<double>>> u, std::vector<std::vector<std::vector<double>>> v, std::vector<std::vector<std::vector<double>>> w, std::vector<std::vector<std::vector<double>>> P, int nx, int ny, int nz, double dx, double dy, double dz, double x0, double y0, double z0, double R, double radius)
{
	int latitudes = 70, longitudes = 70;
	double*** spheremesh = compute_spheremesh(radius, latitudes, longitudes);
	print_rss_memory("after spheremesh");
	tensor*** stresses = allocate3dt(nx, ny, nz);
	print_rss_memory("after stresses");
	std::vector<std::vector<std::vector<std::vector<double>>>> a = { u,v,w };
	for (int x = 0; x < nx; x++)
		for (int y = 0; y < ny; y++)
			for (int z = 0; z < nz; z++)
			{ 
				tensor grad = gradv(x, y, z, nx,ny,nz, dx,dy,dz,a); //uvw
				tensor p(P[x][y][z]);
				stresses[x][y][z] = grad.transpose() + grad - p;
			}
	//determine the force and torque 
	vector answers[2]={ force(x0, y0, z0, dx, dy, dz, nx, ny, nz, longitudes, latitudes, radius, spheremesh, sphere, stresses, R), torque(x0, y0, z0, dx, dy, dz, nx, ny, nz, longitudes, latitudes, radius, spheremesh, sphere, stresses, R) };
	return answers;
}
std::tuple<vector, vector> forceandtorquestag(std::vector<std::vector<std::vector<bool>>> sphere, std::vector<std::vector<std::vector<double>>> u, std::vector<std::vector<std::vector<double>>> v, std::vector<std::vector<std::vector<double>>> w, std::vector<std::vector<std::vector<double>>> P, int nx, int ny, int nz, double dx, double dy, double dz, double x0, double y0, double z0, double R, double radius)
{
	int latitudes = 70, longitudes = 70;
	double*** spheremesh = compute_spheremesh(radius, latitudes, longitudes);
	print_rss_memory("after spheremesh");
	tensor*** stresses = allocate3dt(nx, ny, nz);
	print_rss_memory("after stresses");
	std::vector<std::vector<std::vector<double>>> us(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz)));
	std::vector<std::vector<std::vector<double>>> vs(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz)));
	std::vector<std::vector<std::vector<double>>> ws(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz)));
	for (int x = 0; x < nx; x++)
		for (int y = 0; y < ny; y++)
			for (int z = 0; z < nz; z++)
			{
				us[x][y][z] = (u[x][y][z] + u[x + 1][y][z]) / 2;
				vs[x][y][z] = (v[x][y][z] + v[x][y+1][z]) / 2;
				ws[x][y][z] = (w[x][y][z] + w[x][y][z+1]) / 2;
			}
	std::vector<std::vector<std::vector<std::vector<double>>>> a = { us,vs,ws};
	for (int x = floor((x0-radius)/dx); x <= floor((x0 + radius) / dx); x++)
		for (int y = floor((R+y0 - radius) / dy); y <= floor((R + y0 + radius) / dy); y++)
			for (int z = floor((R + z0 - radius) / dz); z < floor((R + z0 + radius) / dz); z++)
			{
				tensor grad = gradv(x, y, z, nx, ny, nz, dx, dy, dz, a); //uvw
				tensor p(P[x][y][z]);
				stresses[x][y][z] = grad.transpose() + grad - p;
			}
	std::cout << "Area: " << area(x0, y0, z0, dx, dy, dz, nx, ny, nz, longitudes, latitudes, radius, spheremesh, sphere, stresses, R);
	//determine the force and torque
	vector Force = force(x0, y0, z0, dx, dy, dz, nx, ny, nz, longitudes, latitudes, radius, spheremesh, sphere, stresses, R);
	vector Torque = torque(x0, y0, z0, dx, dy, dz, nx, ny, nz, longitudes, latitudes, radius, spheremesh, sphere, stresses, R);
	return std::make_tuple(Force, Torque);
}