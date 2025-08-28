#include "utilities.hpp"
#include "vector.hpp"
#include "mem.hpp"
#include <vector>
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
	tensor output =shear(x, y, z)-pressure;
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
tensor sheara(double**** a, int x, int y, int z)
{
	vector c1 = diff(x, y, z, 1);
	vector c2 = diff(x, y, z, 2);
	vector c3 = diff(x, y, z, 3);
	tensor grad(c1, c2, c3);
	return grad + grad.transpose();
}
inline __attribute__((always_inline)) double partdif(double **** a,int x, int y, int z, int n, double h, int d,int e)//array,x,y,z ,number of spacings, spacing, component of vector, direction of differentiation: 0=x,1=y,2=z
{
	double answer = 0;
	switch (e)
	{
	case 0://x diff
		if (x >= 1) 
		{
			if (x < n - 1)//center
			{
				answer = (a[d][x + 1][y][z] -a[d][x-1][y][z])/2/h;
			}
			else//left side
			{
				answer = (3 * a[d][x][y][z] - 4 * a[d][x - 1][y][z] + 1 * a[d][x - 2][y][z]) / 2 / h;

			}
		}
		else//right side
		{
			answer = (-3 * a[d][x][y][z] + 4 * a[d][x + 1][y][z] - 1 * a[d][x + 2][y][z]) / 2 / h;
		}
		break;
	case 1: //y diff
		if (y >= 1)
		{
			if (y < n - 1)
			{
				answer = (a[d][x][y+1][z] - a[d][x][y-1][z]) / 2 / h;
			}
			else
			{
				answer = (3 * a[d][x][y][z] - 4 * a[d][x][y - 1][z] + 1 * a[d][x][y - 2][z]) / 2 / h;
			}
		}
		else
		{
			answer = (-3 * a[d][x][y][z] + 4 * a[d][x][y + 1][z] - 1 * a[d][x][y + 2][z]) / 2 / h;
		}
		break;
	case 2:
		//diffz
		if (z >= 1)
		{
			if (z < n - 1)
			{
				answer = (a[d][x][y][z+1] - a[d][x][y][z-1]) / 2 / h;
			}
			else
			{
				answer = (3 * a[d][x][y][z] - 4 * a[d][x][y][z - 1] + 1 * a[d][x][y][z - 2]) / 2 / h;
			}
		}
		else
		{
			answer = (-3 * a[d][x][y][z] + 4 * a[d][x][y][z + 1] - 1 * a[d][x][y][z + 2]) / 2 / h;
		}
		break;
	default:
		std::cout << "error in partdif: " << x << " " << y  << " " << z << " " << n << " " << h << " " << d << " "<< e;
		exit(1);
	}
	return answer;
}
double partdif(std::vector< std::vector<std::vector<std::vector<double>>>>& a, int x, int y, int z, int n, double h, int d, int e)//array,x,y,z ,number of spacings, spacing, component of vector, direction of differentiation: 0=x,1=y,2=z
{
	double answer = 0;
	switch (e)
	{
	case 0://x diff
		if (x >= 1)
		{
			if (x < n - 1)//center
			{
				answer = (a[d][x + 1][y][z] - a[d][x - 1][y][z]) / 2 / h;
			}
			else//left side
			{
				answer = (3 * a[d][x][y][z] - 4 * a[d][x - 1][y][z] + 1 * a[d][x - 2][y][z]) / 2 / h;

			}
		}
		else//right side
		{
			answer = (-3 * a[d][x][y][z] + 4 * a[d][x + 1][y][z] - 1 * a[d][x + 2][y][z]) / 2 / h;
		}
		break;
	case 1: //y diff
		if (y >= 1)
		{
			if (y < n - 1)
			{
				answer = (a[d][x][y + 1][z] - a[d][x][y - 1][z]) / 2 / h;
			}
			else//n-1 only
			{
				answer = (3 * a[d][x][y][z] - 4 * a[d][x][y - 1][z] + 1 * a[d][x][y - 2][z]) / 2 / h;
			}
		}
		else //0 only
		{
			answer = (-3 * a[d][x][y][z] + 4 * a[d][x][y + 1][z] - 1 * a[d][x][y + 2][z]) / 2 / h;
		}
		break;
	case 2:
		//diffz
		if (z >= 1)
		{
			if (z < n - 1)
			{
				answer = (a[d][x][y][z + 1] - a[d][x][y][z - 1]) / 2 / h;
			}
			else
			{
				answer = (3 * a[d][x][y][z] - 4 * a[d][x][y][z - 1] + 1 * a[d][x][y][z - 2]) / 2 / h;
			}
		}
		else
		{
			answer = (-3 * a[d][x][y][z] + 4 * a[d][x][y][z + 1] - 1 * a[d][x][y][z + 2]) / 2 / h;
		}
		break;
	default:
		std::cout << "error in partdif: " << x << " " << y << " " << z << " " << n << " " << h << " " << d << " " << e;
		exit(1);
	}
	return answer;
}
double partdif(std::vector<std::vector<std::vector<double>>>& a, int x, int y, int z, int n, double h, int e)//array,x,y,z ,number of spacings, spacing, direction of differentiation: 0=x,1=y,2=z
{
	double answer = 0;
	switch (e)
	{
	case 0://x diff
		if (x >= 1)
		{
			if (x < n - 1)//center
			{
				answer = (a[x + 1][y][z] - a[x - 1][y][z]) / 2 / h;
			}
			else//left side
			{
				answer = (3 * a[x][y][z] - 4 * a[x - 1][y][z] + 1 * a[x - 2][y][z]) / 2 / h;

			}
		}
		else//right side
		{
			answer = (-3 * a[x][y][z] + 4 * a[x + 1][y][z] - 1 * a[x + 2][y][z]) / 2 / h;
		}
		break;
	case 1: //y diff
		if (y >= 1)
		{
			if (y < n - 1)
			{
				answer = (a[x][y + 1][z] - a[x][y - 1][z]) / 2 / h;
			}
			else//n-1 only
			{
				answer = (3 * a[x][y][z] - 4 * a[x][y - 1][z] + 1 * a[x][y - 2][z]) / 2 / h;
			}
		}
		else //0 only
		{
			answer = (-3 * a[x][y][z] + 4 * a[x][y + 1][z] - 1 * a[x][y + 2][z]) / 2 / h;
		}
		break;
	case 2:
		//diffz
		if (z >= 1)
		{
			if (z < n - 1)
			{
				answer = (a[x][y][z + 1] - a[x][y][z - 1]) / 2 / h;
			}
			else
			{
				answer = (3 * a[x][y][z] - 4 * a[x][y][z - 1] + 1 * a[x][y][z - 2]) / 2 / h;
			}
		}
		else
		{
			answer = (-3 * a[x][y][z] + 4 * a[x][y][z + 1] - 1 * a[x][y][z + 2]) / 2 / h;
		}
		break;
	default:
		std::cout << "error in partdif: " << x << " " << y << " " << z << " " << n << " " << h  << " " << e;
		exit(1);
	}
	return answer;
}
double partdif2(std::vector< std::vector<std::vector<std::vector<double>>>>& a, int x, int y, int z, int n, double h, int d, int e)//array,x,y,z ,number of spacings, spacing, component of vector, direction of differentiation: 0=x,1=y,2=z
{
	double answer = 0;
	switch (e)
	{
	case 0://x diff
		if (x >= 1)
		{
			if (x < n - 1)//center
			{
				answer = (a[d][x + 1][y][z] -2* a[d][x][y][z] + a[d][x - 1][y][z]) / h/h;
			}
			else//left side 2 -5 4 -1
			{
				answer = (2 * a[d][x][y][z] - 5 * a[d][x - 1][y][z] + 4 * a[d][x - 2][y][z]- a[d][x - 3][y][z]) / h/h;

			}
		}
		else//right side
		{
			answer = (2 * a[d][x][y][z] - 5 * a[d][x + 1][y][z] + 4 * a[d][x + 2][y][z] - a[d][x + 3][y][z]) / h / h;
		}
		break;
	case 1: //y diff
		if (y >= 1)
		{
			if (y < n - 1)//center
			{
				answer = (a[d][x][y+1][z] - 2 * a[d][x][y][z] + a[d][x][y-1][z]) / h / h;
			}
			else//left side 2 -5 4 -1
			{
				answer = (2 * a[d][x][y][z] - 5 * a[d][x][y-1][z] + 4 * a[d][x][y-2][z] - a[d][x][y-3][z]) / h / h;

			}
		}
		else//right side
		{
			answer = (2 * a[d][x][y][z] - 5 * a[d][x][y+1][z] + 4 * a[d][x][y+2][z] - a[d][x][y+3][z]) / h / h;
		}
		break;
	case 2:
		//diffz
		if (z >= 1)
		{
			if(z < n - 1)//center
			{
				answer = (a[d][x][y][z+1] - 2 * a[d][x][y][z] + a[d][x][y][z-1]) / h / h;
			}
			else//left side 2 -5 4 -1
			{
				answer = (2 * a[d][x][y][z] - 5 * a[d][x][y][z-1] + 4 * a[d][x][y][z-2] - a[d][x][y][z-3]) / h / h;

			}
		}
		else//right side
		{
			answer = (2 * a[d][x][y][z] - 5 * a[d][x][y][z+1] + 4 * a[d][x][y][z+2] - a[d][x][y][z+3]) / h / h;
		}
		break;
	default:
		std::cout << "error in partdif: " << x << " " << y << " " << z << " " << n << " " << h << " " << d << " " << e;
		exit(1);
	}
	return answer;
}

vector *** divtenall(tensor*** a, const int n, double h)
{
	double **** array= allocate4d(9,n,n,n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			for (int k = 0; k < n; k++)
			{
				for (int l = 1; l <= 3; l++)//column
				{
					for (int m = 1; m <= 3; m++)//row
					{
						array[3 * l + m - 4][i][j][k] = a[i][j][k].getValue(l, m); //0,3,6  1,4,7   2,4,8
					}
				}
			}
		}
	}
	vector *** answer = allocate3dv(n,n,n);
	for (int x = 0; x < n; x++)
	{
		for (int y = 0; y < n; y++)
		{
			for (int z = 0; z < n; z++)
			{
				double temp[3] = { 0,0,0 };
				for (int l = 1; l <= 3; l++)//column
				{
					for (int m = 1; m <= 3; m++)//row
					{
						temp[l-1] += partdif(array, x, y, z, n, h, 3 * l + m - 4, m-1);
					}
				}
				answer[x][y][z].set(temp[0],temp[1],temp[2]);
			}
		}
	}
	return answer;
}
vector velocity(double x, double y, double z)
{
	vector output(4*y*y,0,0);
	return output;
}
double divvadvecv(double **** a, double x, double y, double z,int n, double h)
{
	double sum=0;
	for (int d=0 ; d<=2; d++)
	{
		for (int e=0 ; e<=2; e++)
			{
				sum+=partdif(a,x, y, z, n, h, d,e)*partdif(a,x, y, z, n, h, e,d);
			}
	}
	return sum;
}
