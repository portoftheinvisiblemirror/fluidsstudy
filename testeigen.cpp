#include <Dense>
#include <Sparse>
#include <iostream>
#include <vector>
int main()
{
	int cube = 1;
	Eigen::VectorXd x(8*cube*cube*cube), b(8*cube*cube*cube);
	std::vector<Eigen::Triplet<double>> tripletlist;
	tripletlist.reserve(7*8 * cube * cube * cube);
	for (int i = 0; i < 8 * cube * cube * cube; ++i)
	{
		tripletlist.push_back(Eigen::Triplet<double>(i, i, 1));
		b(i)=i+1;
	}
	tripletlist.push_back(Eigen::Triplet <double>(1, 0, 1));
	Eigen::SparseMatrix<double> A(8 * cube * cube * cube, 8 * cube * cube * cube);
	A.setFromTriplets(tripletlist.begin(), tripletlist.end());
	std::cout << A << std::endl;
	
	
	{
		Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver;
		// fill A and b;
		// Compute the ordering permutation vector from the structural pattern of A
		solver.analyzePattern(A);
		// Compute the numerical factorization 
		solver.factorize(A);
		//Use the factors to solve the linear system 
		x = solver.solve(b);
		std::cout << x << std::endl;
	}
	{
		Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>   solver;
		// fill A and b;
		// Compute the ordering permutation vector from the structural pattern of A
		solver.analyzePattern(A);
		// Compute the numerical factorization 
		solver.factorize(A);
		//Use the factors to solve the linear system 
		x = solver.solve(b);
		std::cout << x << std::endl;
	}
}