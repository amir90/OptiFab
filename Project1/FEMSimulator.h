#pragma once
#include "HelperTop.h" //need to generalize this
#include <set>

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

Eigen::MatrixXd B_mat(float x, float y, float z);
Eigen::MatrixXd calc_N(float x, float y, float z);
Eigen::MatrixXd K_mat();
Eigen::MatrixXd M_mat();
Eigen::SparseMatrix<double> build_global_stiffness_matrix(std::vector<designVariable> const &x, int numNodes, Eigen::MatrixXd KE, int power, std::vector<int> freeDOF, int numOfFree);
Eigen::SparseMatrix<double> build_global_stiffness_matrix_full(std::vector<designVariable> const &x, Eigen::MatrixXd KE, int power, std::vector<std::vector<double>> &nodes);
Eigen::VectorXd do_FEM(SpMat A, Eigen::VectorXd f);
struct comp {
	bool operator() (const std::pair<int, Eigen::Vector3f>& a, const  std::pair<int, Eigen::Vector3f>& b) const {

		return a.first > b.first;
	}
};
Eigen::VectorXd doFEM(std::vector<std::vector<double>> const nodes, std::vector<designVariable>  const x, std::set < std::pair<int, Eigen::Vector3f>, comp> forceVertexIdSet, std::set<int> constraintVertexIdSet, float dx, double E);
Eigen::SparseMatrix<double> build_global_mass_matrix(std::vector<designVariable> const &x, int numNodes, Eigen::MatrixXd KE, int power, std::vector<int> freeDOF, int numOfFree);
