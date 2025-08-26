#include "mem.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

void print_rss_memory(const std::string &message) {
	// Path to the statm file for the current process
	std::cout << message << "\n";
	const char* statm_path = "/proc/self/statm";
	FILE* file = fopen(statm_path, "r");

	if (!file) {
		perror("Failed to open statm file");
		return;
	}

	// Variables to hold memory usage information
	unsigned long size, rss;

	// Read the first two values from the statm file
	if (fscanf(file, "%lu %lu", &size, &rss) != 2) {
		perror("Failed to read statm file");
		fclose(file);
		return;
	}

	fclose(file);

	// Page size in bytes
	long page_size = sysconf(_SC_PAGESIZE);
	if (page_size == -1) {
		perror("Failed to get page size");
		return;
	}

	// Calculate RSS in bytes
	unsigned long rss_in_bytes = rss * page_size;

	// Print the results
	printf("Resident Set Size (RSS): %lu bytes\n", rss_in_bytes);
}

// Implementation file for memory management related utilities
double** allocate2d(const unsigned int x, const unsigned int y, double** array2d) {
	// Allocate memory blocks of size
	// x i.e., no of 2D Arrays
	array2d = new double* [x];

	for (unsigned int i = 0; i < x; i++) {

		// Allocate memory blocks for
		// rows of each 2D array
		array2d[i] = new double[y];
	}
	return array2d;
}

double*** allocate3d(const unsigned int latitudes, const unsigned int longitudes, const unsigned int narea) {
	// Dimensions of the 3D array
	int x = latitudes, y = longitudes, z = narea;
	int count = 0;

	// Allocate memory blocks of size
	// x i.e., no of 2D Arrays
	double *** spheremesh = new double** [x];

	for (int i = 0; i < x; i++) {

		// Allocate memory blocks for
		// rows of each 2D array
		spheremesh[i] = new double* [y];

		for (int j = 0; j < y; j++) {

			// Allocate memory blocks for
			// columns of each 2D array
			spheremesh[i][j] = new double[z];
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
	//std::cout << latitudes << " " << longitudes << " " << narea << std::endl;
	//int *** ispheremesh = new int **[latitudes];
	//for(unsigned int i = 0; i < latitudes; i ++){
	//  ispheremesh[i] = new int *[longitudes];
	//  for(unsigned int j = 0; j < narea; j ++){
	//    ispheremesh[i][j] = new int [narea];
	//  }
	//}
#ifdef DEBUG
	printf("Address: %X \n", spheremesh);
#endif
	//printf("Check value: %d \n", spheremesh[10][10][3]);
	return spheremesh;
}
tensor*** allocate3dt(const unsigned int latitudes, const unsigned int longitudes, const unsigned int narea) {
	// Dimensions of the 3D array
	int x = latitudes, y = longitudes, z = narea;
	int count = 0;

	// Allocate memory blocks of size
	// x i.e., no of 2D Arrays
	tensor*** spheremesh = new tensor** [x];

	for (int i = 0; i < x; i++) {

		// Allocate memory blocks for
		// rows of each 2D array
		spheremesh[i] = new tensor* [y];

		for (int j = 0; j < y; j++) {

			// Allocate memory blocks for
			// columns of each 2D array
			spheremesh[i][j] = new tensor[z];
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
#ifdef DEBUG
	printf("Address: %X \n", spheremesh);
#endif
	//printf("Check value: %d \n", spheremesh[10][10][3]);
	return spheremesh;
}
double**** allocate4d(const unsigned int nc, const unsigned int nx, const unsigned int ny, const unsigned int nz) {
	// Allocate memory blocks of size
	// nc, nx, ny, nz
	double **** ptr = new double*** [nc]; // ptr is a pointer to an array of 3d matrices

	for (int i = 0; i < nc; i++) {
		ptr[i] = new double** [nx]; // ptr[i] is a pointer to an array 2d matrices
		for (int j = 0; j < nx; j++) {
			ptr[i][j] = new double* [ny];
			for (int k = 0; k < nz; k++) {
				ptr[i][j][k] = new double [nz];
			} 
		}
	}
#ifdef DEBUG
	printf("Address: %X \n", ptr);
#endif
	return ptr;
}
vector*** allocate3dv(const unsigned int nc, const unsigned int nx, const unsigned int ny) {
	// Allocate memory blocks of size
	// nc, nx, ny, nz
	vector*** ptr = new vector** [nc]; // ptr is a pointer to an array of 3d matrices

	for (int i = 0; i < nc; i++) {
		ptr[i] = new vector* [nx]; // ptr[i] is a pointer to an array 2d matrices
		for (int j = 0; j < nx; j++) {
			ptr[i][j] = new vector [ny];
		}
	}
#ifdef DEBUG
	printf("Address: %X \n", ptr);
#endif
	return ptr;
}
void deletia(tensor***& spheremesh, const unsigned int latitudes, const unsigned int longitudes, const unsigned int narea)
{
	int x = latitudes, y = longitudes, z = narea;
	for (int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++) {
			delete [] spheremesh[i][j];
		}
		delete [] spheremesh[i];
	}
	delete [] spheremesh;
}
void deletia(double***& spheremesh, const unsigned int latitudes, const unsigned int longitudes, const unsigned int narea)
{
	int x = latitudes, y = longitudes, z = narea;
	for (int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++) {
			delete[] spheremesh[i][j];
		}
		delete[] spheremesh[i];
	}
	delete[] spheremesh;
}