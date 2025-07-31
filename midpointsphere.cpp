#include <vector>
#include <cmath>
#include <algorithm>
void fillin(std::vector<bool> row, int i1, int i2)
{
	for (int k = i1; k <= i2; ++k)
	{
		row[k] = true;
	}
}
//index 0, 0, 0 <=> cartesian 0,-R,-R
void filledmidpointcircle(const double y0, const double z0, const double dy, const double dz, const int ny, const int nz, const double radius, std::vector<std::vector<bool>> circle, const double R)
{
	//Step 1: find the diameter at y = ycenter
	
	//Part a: find the indices of the center
	int ycenter = floor((y0 + R) / dy), zcenter = floor((z0 + R) / dz);
	//Part b: fill in the z values across the y=ycenter diameter.
	fillin(circle[ycenter], floor((z0 + R - radius) / dz), floor((z0 + R + radius) / dz));
	
	//Step 2: repeat filling in y rows for each y index in the circle
	//Part a: determine the limits
	int yhigh = floor((y0 + R + radius) / dy), ylow = floor((y0 - radius + R) / dy);

	//Part b: start looping
	for (int i = ycenter + 1; i <= yhigh; ++i) //bound this to make sure yhigh is not greater than ny-1?
	{
		//determine half length
		double coord = -R + i * dy - y0;
		double halfz = sqrt(radius * radius - coord*coord); 
		//fill in
		fillin(circle[i], floor((z0 + R - halfz) / dz), floor((z0 + R + halfz) / dz));
	}
	for (int i = ycenter - 1; i >= ylow; --i)
	{
		//determine half length
		double coord = -R + (i+1) * dy - y0; //include i+1 to maximize half z
		double halfz = sqrt(radius * radius - coord * coord);
		//fill in
		fillin(circle[i], floor((z0 + R - halfz) / dz), floor((z0 + R + halfz) / dz));
	}
}

std::vector<std::vector<std::vector<bool>>> filledmidpointsphere(const double x0, const double y0, const double z0, const double dx, const double dy, const double dz, const double nx, const int ny, const int nz, const double radius, const double R)
{
	std::vector<std::vector<std::vector<bool>>> sphere(nx, std::vector<std::vector<bool>>(ny, std::vector<bool>(nz)));
	//Step 1: find the great circle at x = xcenter

	//Part a: find the indices of the center
	int xcenter = floor(x0/dx), ycenter = floor((y0 + R) / dy), zcenter = floor((z0 + R) / dz);
	//Part b: fill in the y, z values x=xcenter great circle.
	filledmidpointcircle(y0, z0, dy, dz, ny, nz, radius, sphere[xcenter], R);

	//Step 2: repeat filling in y,z circles for each x index in the sphere
	//Part a: determine the limits
	int xhigh = floor((x0+ radius) / dx), xlow = floor((x0 - radius) / dx);

	//Part b: start looping
	for (int i = xcenter + 1; i <= xhigh; ++i)
	{
		//determine half radius
		double coord = i * dx - x0;
		double halfr = sqrt(radius * radius - coord*coord); 
		//fill in
		filledmidpointcircle(y0, z0, dy, dz, ny, nz, halfr, sphere[i], R);
	}

	for (int i = xcenter - 1; i >= xlow; --i)
	{
		//determine half radius
		double coord = (i+1) * dx - x0;
		double halfr = sqrt(radius * radius - coord * coord); 
		//fill in
		filledmidpointcircle(y0, z0, dy, dz, ny, nz, halfr, sphere[i], R);
	}
	return sphere;
}