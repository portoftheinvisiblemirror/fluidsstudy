#include "mem.hpp"
// Implementation file for memory management related utilities
long double** allocate2d(const unsigned int x, const unsigned int y, long double** array2d) {
	// Allocate memory blocks of size
	// x i.e., no of 2D Arrays
	array2d = new long double* [x];

	for (unsigned int i = 0; i < x; i++) {

		// Allocate memory blocks for
		// rows of each 2D array
		array2d[i] = new long double[y];
	}
	return array2d;
}

long double*** allocate3d(const unsigned int latitudes, const unsigned int longitudes, const unsigned int narea) {
	// Dimensions of the 3D array
	int x = latitudes, y = longitudes, z = narea;
	int count = 0;

	// Allocate memory blocks of size
	// x i.e., no of 2D Arrays
	long double *** spheremesh = new long double** [x];

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
long double**** allocate4d(const unsigned int nc, const unsigned int nx, const unsigned int ny, const unsigned int nz) {
	// Allocate memory blocks of size
	// nc, nx, ny, nz
	long double **** ptr = new long double*** [nc]; // ptr is a pointer to an array of 3d matrices

	for (int i = 0; i < nc; i++) {
		ptr[i] = new long double** [nx]; // ptr[i] is a pointer to an array 2d matrices
		for (int j = 0; j < nx; j++) {
			ptr[i][j] = new long double* [ny];
			for (int k = 0; k < nz; k++) {
				ptr[i][j][k] = new long double [nz];
			} 
		}
	}
#ifdef DEBUG
	printf("Address: %X \n", ptr);
#endif
	return ptr;
}
