#include <iostream>
#include <Eigen/Sparse>
#include <bench/BenchTimer.h>
#include <vector>

bool debug_flag;

using namespace Eigen;
int main() {
	int n =186000;
	int nnz = 24;

	typedef Triplet<double> T;
	std::vector<T> elements;
	elements.reserve(n*nnz);
	for (int j = 0; j < n; ++j)
		for (int k = 0; k < nnz; ++k)
		{
			int i;
			do {
				i = internal::random<int>(0, n - 1);
			} while (D(i, j) != 0);
			D(i, j) = 1;
			elements.push_back(T(i, j, 1));
		}


	BenchTimer t1, t2, t3, t4, t5, t6;

	t1.start();
	{
		std::vector<T> el;
		el.reserve(n*nnz);
		for (int k = 0; k < elements.size(); ++k)
			el.push_back(T(elements[k].row(), elements[k].col(), elements[k].value()));
		SparseMatrix<double> S(n, n);
		S.setFromTriplets(el.begin(), el.end());
	}
	t1.stop();
	std::cout << "setFromTriplets:     " << t1.best() * 1000 << "ms\n";

	t3.start();
	{
		SparseMatrix<double> S(n, n);
		S.reserve(VectorXi::Constant(n, nnz));
		for (int k = 0; k < elements.size(); ++k)
			S.insert(elements[k].row(), elements[k].col()) = elements[k].value();
		S.makeCompressed();
	}
	t3.stop();
	std::cout << "reserve+insert:      " << t3.best() * 1000 << "ms\n";

	t4.start();
	{
		SparseMatrix<double> S(n, n);
		S.reserve(VectorXi::Constant(n, nnz));
		for (int k = 0; k < elements.size(); ++k)
			S.coeffRef(elements[k].row(), elements[k].col()) = elements[k].value();
		S.makeCompressed();
	}
	t4.stop();
	std::cout << "reserve+coeffRef:    " << t4.best() * 1000 << "ms\n";

	t5.start();
	{
		SparseMatrix<double> S(n, n);
		for (int k = 0; k < elements.size(); ++k)
			S.coeffRef(elements[k].row(), elements[k].col()) = elements[k].value();
		S.makeCompressed();
	}
	t5.stop();
	std::cout << "naive coeffRef:      " << t5.best() * 1000 << "ms\n";

}