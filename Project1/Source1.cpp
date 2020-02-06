#pragma once
#include<iostream>
#include<fstream>
#include <MMASolver.h>
#include <Eigen\Core>
#include <Eigen\Dense>
#include<igl\cotmatrix.h>
#include <igl\readOBJ.h>
#include<igl\opengl\glfw\Viewer.h>
#include <igl/unproject_onto_mesh.h>
#include "HelperTop.h"
#include "FEMSimulator.h"
#include "StressTuner.h"
#include "ComplianceTuner.h"
#include "fatigueOptimizer.h"
#include <igl/jet.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/png/texture_from_file.h>
#include <glut.h>
#include "NSGA3/problem_base.h"
#include "NSGA3/alg_nsgaiii.h"
#include "NSGA3/alg_population.h"
#include "NSGA3/exp_experiment.h"
#include <fstream>
#include <iostream>
#include "NSGAIIIComp.h"
#include <igl/edge_topology.h>
#include <igl/per_face_normals.h>


template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

	// initialize original index locations
	std::vector<size_t> idx(v.size());
	std::iota(idx.begin(), idx.end(), 0);
	// sort indexes based on comparing values in v
	// using std::stable_sort instead of std::sort
	// to avoid unnecessary index re-orderings
	// when v contains elements of equal values 
	std::stable_sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

	return idx;
}

bool debug_flag = false;

template <typename T>
std::vector<T> linspace(T a, T b, size_t N) {
	T h = (b - a) / static_cast<T>(N - 1);
	std::vector<T> xs(N);
	typename std::vector<T>::iterator x;
	T val;
	for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
		*x = val;
	return xs;
}

int main() {

	bool voxelTestFlag = false;

	std::vector<int> mask;

	int argc = 1;
	char *argv[1] = { (char*)"Something" };

	//initialize voxel parameters
	int numOfVoxelsX =60;
	int numOfVoxelsY = 30;
	int numOfVoxelsZ = 1;
	double dx = 0.001; // size
	int NeighbourLayers = 1;
	float r0 = 1.5*dx*NeighbourLayers;

	std::set < std::pair<int, Eigen::Vector3f>, comp> forceVertexIdSet;

	std::set < std::pair<int, Eigen::Vector3f>, comp> forceVertexIdSet2; //for NSGAIII test

	std::set<int> constraintVertexIdSet;


	//declare geometry parameters
	std::vector<designVariable> x;
	std::vector<std::vector<double>> nodes; //create array that connects node indices to coordinates.
	nodes.resize((numOfVoxelsX + 1) *(numOfVoxelsY + 1)*(numOfVoxelsZ + 1));

	create_mesh(nodes, x,  numOfVoxelsX,  numOfVoxelsY,  numOfVoxelsZ,dx);

//	getBCforVoxelTest(x, nodes, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ,forceVertexIdSet, constraintVertexIdSet);


	double change = 1000;

	int power = 3;

	int C = 1;

	double E = 70000000000;

	std::vector<int> freeDOF;
	freeDOF.resize(nodes.size());
	int ind = 0; //number of free nodes
	for (int i = 0; i < freeDOF.size(); i++) {

		if (constraintVertexIdSet.find(i) == constraintVertexIdSet.end()) {
			//if (constraints.find(i) == constraints.end()) {

			freeDOF[i] = ind;
			ind++;
		}
		else {

			freeDOF[i] = -1;
		}

	}

	//int numOfFixed = nodes.size() - ind; //number of fixed nodes


	Eigen::VectorXd f(3 * ind);

	for (int j = 0; j < 3 * ind; j++) {
		f(j) = 0;
	}

	for (auto j = forceVertexIdSet.begin(); j != forceVertexIdSet.end(); j++) {
		int tempForceid = freeDOF[j->first] * 3;
		f[tempForceid + 0] = j->second[0];// force[temp_f].second[0];
		f[tempForceid + 1] = j->second[1];//force[temp_f].second[1];
		f[tempForceid + 2] = j->second[2]; //force[temp_f].second[2];//

	}

//	SpMat K_global= build_global_stiffness_matrix(x, nodes.size(), K_mat(), power, freeDOF,ind);

	double delta = 0;

	Eigen::VectorXd u_dot_0, u0;

	u_dot_0 = Eigen::VectorXd::Zero(nodes.size());

	u0 = Eigen::VectorXd::Zero(nodes.size());

	std::vector<int> targets = linspace<int>(0, x.size() - 1, x.size());

	std::cout << "targets size: " << targets.size() << std::endl;

	std::vector<std::vector<double>> sigma_m;

	std::vector<std::vector<Stress>> sigma_a;

	std::vector<std::vector<Eigen::VectorXd>> u_max;

	std::vector<std::vector<Eigen::VectorXd>> u_min;

	int steps = 15;

	Eigen::MatrixXd disp_test;

	igl::opengl::glfw::Viewer viewer; //sets opengl context (important for images of buttons)

		//initialize matrices for faces and vertices - for libigl display
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	Eigen::MatrixXd N;

	Eigen::MatrixXd Color;

	Eigen::VectorXd test_u = Eigen::VectorXd::Zero(nodes.size()*3);

	std::vector<double> ck_iter(C, 1);

	double DPN_k = 0;

	double max_Dk = -1;

	double tau=0.5;

	for (int iter = 0; iter < 0; iter++) {

		sigma_a.clear();

		sigma_m.clear();

		u_max.clear();

		u_min.clear();

		u_dot_0 = Eigen::VectorXd::Zero(nodes.size());

		u0 = Eigen::VectorXd::Zero(nodes.size());

	//	u_max.clear(); u_min.clear();

		calc_rho(nodes, x, -1, dx, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ, r0);

		Eigen::VectorXd u_temp = doFEM(nodes, x, forceVertexIdSet, constraintVertexIdSet, dx, E,power);

		std::cout << "max FEM displacement: " << u_temp.maxCoeff() << std::endl;

	//	ClacStressTimeHistory(targets, x, 0.001, 0.001, power, forceVertexIdSet, dx, E, nodes, constraintVertexIdSet, u_dot_0, u0, 0.0001, steps, sigma_m, sigma_a, disp_test, u_max, u_min);

	//	std::cout << "number of u_max: " << u_max[3].size() << "/" << sigma_a[3].size() << std::endl;

		//	std::cout << "calculated stress histories" << std::endl;
		if (debug_flag) {
			for (int i = 0; i < targets.size(); i++) {

				for (int j = 0; j < sigma_m[i].size(); j++) {

					std::cout << "i: " << i << " j: " << j << " sigma_m: " << sigma_m[i][j] << " sigma_a: " << sigma_a[i][j].vonMises << std::endl;

				}
			}
		}

		double bf = -0.075;
		double sigma_f = 650000000;
		double ni = 1000000;
		double Sut = 380000000;
		double Sf = 3;
/*
		std::vector<double> damage(x.size(), 0);

		for (int i = 0; i < targets.size(); i++) {

			std::vector<double> N(sigma_m[i].size());
			//calculate damage for target
			double D = 0;
			for (int j = 0; j < N.size(); j++) {
				if (((sigma_f*(Sut - std::max(sigma_m[i][j], 0.0)))) < 0) {
					D = 30;
					break;
				}
				D += ni * 2 * std::pow(Sut*sigma_a[i][j].vonMises / (sigma_f*(Sut - std::max(sigma_m[i][j], 0.0))), -1.0 / bf);
				if (D > 1) { //is clipping the correct approach?
					D = 30;
					break;
				}
			}

			damage[i] = D;

			assert(D >= 0);
		}


		//split the damages in to groups

		auto sorted = sort_indexes(damage);

		std::vector<std::vector<int>> Omega_k;

		for (int i = 0; i < C; i++) {

			g[i] = 0;

			Omega_k.push_back(std::vector<int>(sorted.begin() + std::ceil((float)i / C * sorted.size()), sorted.begin() + std::ceil((float)(i + 1) / C * sorted.size())));

			std::cout << Omega_k[0].size() << std::endl;

			for (int id : Omega_k[i]) {

				g[i] += std::pow(damage[id],power)*x[id].rho;

				if (max_Dk < damage[id]) {

					max_Dk = damage[id];

				}

			}

	//		std::cout << i << ": " << max_Dk << std::endl;

			g[i] = ck_iter[i]*std::pow(g[i], 1.0 / power)-1;

			ck_iter[i] = ck_iter[i] *(tau * max_Dk/(g[i]+1) + (1 - tau));

			assert(ck_iter[i] >= 0);

		}

		std::vector<int> n(sigma_a.size(), ni);

		full_dD_PN_k_dgamma_e(x, Omega_k, n, power, u_max, u_min, sigma_a, damage, Sut, bf, Sf, sigma_m, E, nodes, r0,dg);

		change = optimizeFatigue(xnew, dg, g, df, xmin, xmax, mma, nodes, x, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ, r0, dx);

		if (iter == -1) {

			makeSurfaceMesh3(x, nodes);
			igl::readOBJ("mesh.obj", V, F);
			viewer.data().set_mesh(V, F);

			Eigen::VectorXd rhoFaces(F.rows());

			for (int i = 0; i < x.size(); i++) {
				std::cout << dg[i] << std::endl;
				//		std::cout << "element: " << i << ", Damage: " << x[i].rho << std::endl;
				for (int j = 0; j < 12; j++) {
					rhoFaces(i * 12 + j) = x[i].rho;//damage[i];

				}
			}

			igl::jet(rhoFaces, true, Color);

			viewer.data().set_colors(Color);

			viewer.launch();

		}
*/
		/*
		//for debug
		for (int i = 0; i < nodes.size(); i++) {

			nodes[i][0] += u_temp[3*i+0]*50; nodes[i][1] += u_temp[3 * i + 1]*50; nodes[i][2] += u_temp[3 * i + 2]*50;

		}

		makeMesh(x, nodes, mask);

		return 1;
		*/
		//std::cout << "g[0]: " << g[0] << std::endl;
		/*

		Eigen::VectorXd displacement(24);
		
		for (int i = 0; i < x.size(); i++) {


			displacement(0) = u_temp[x[i].nodeIdx[0] * 3]; displacement(1) = u_temp[x[i].nodeIdx[0] * 3 + 1]; displacement(2) = u_temp[x[i].nodeIdx[0] * 3 + 2];
			displacement(3) = u_temp[x[i].nodeIdx[1] * 3]; displacement(4) = u_temp[x[i].nodeIdx[1] * 3 + 1]; displacement(5) = u_temp[x[i].nodeIdx[1] * 3 + 2];
			displacement(6) = u_temp[x[i].nodeIdx[2] * 3]; displacement(7) = u_temp[x[i].nodeIdx[2] * 3 + 1]; displacement(8) = u_temp[x[i].nodeIdx[2] * 3 + 2];
			displacement(9) = u_temp[x[i].nodeIdx[4] * 3]; displacement(10) = u_temp[x[i].nodeIdx[4] * 3 + 1]; displacement(11) = u_temp[x[i].nodeIdx[4] * 3 + 2];
			displacement(12) = u_temp[x[i].nodeIdx[3] * 3]; displacement(13) = u_temp[x[i].nodeIdx[3] * 3 + 1]; displacement(14) = u_temp[x[i].nodeIdx[3] * 3 + 2];
			displacement(15) = u_temp[x[i].nodeIdx[6] * 3]; displacement(16) = u_temp[x[i].nodeIdx[6] * 3 + 1]; displacement(17) = u_temp[x[i].nodeIdx[6] * 3 + 2];
			displacement(18) = u_temp[x[i].nodeIdx[5] * 3]; displacement(19) = u_temp[x[i].nodeIdx[5] * 3 + 1]; displacement(20) = u_temp[x[i].nodeIdx[5] * 3 + 2];
			displacement(21) = u_temp[x[i].nodeIdx[7] * 3]; displacement(22) = u_temp[x[i].nodeIdx[7] * 3 + 1]; displacement(23) = u_temp[x[i].nodeIdx[7] * 3 + 2];

			dc[i] = -E*power * std::pow(x[i].value, power - 1)* displacement.transpose()*K_mat()*displacement;
		//	std::cout << dc[i] << std::endl;

		}
		*/

		//calculate discrete derivative constraints due to changing x[0].value, and compare to the analytical value.

		std::vector<Stress> Stresses = calcStresses(ni, E, nodes, x, u_temp, dx);

		for (int i = 0; i < Stresses.size(); i++) {

			std::cout << "stresses" << Stresses[i].vonMises << std::endl;
		}

//		change = optimizeStress(Stresses, u_temp, C, power, dx, change, nodes, x, xnew, dg, g, df, xmin, xmax, density, E, ni, mma, -1, r0);

	//	if (iter == 0) {
	//		deriv = dg[0];
	//	}

		std::cout << change << std::endl;

	//	delta -= (g[0] * (0.5 - iter) * 2)/ 0.0000000001;

		//test_u += u_temp * (0.5 - iter) * 2;

	//	change = optimizeCompliance(xnew, dg, g, dc, xmin, xmax, dc_filtered, V_allowed, mma, nodes, x, dx,
	//		numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ, r0);

	//	std::cout << "iter: " << iter <<" change: "<<change<< std::endl;

	}

	//std::cout <<"u_norm: "<< test_u.norm() <<" "<< test_u.minCoeff() <<" "<< test_u.maxCoeff()<< std::endl;

	//std::cout << "discrete derivative: : "<< delta << std::endl;
	//std::cout << "analytical derivative: : " << deriv << std::endl;

	for (int i = 0; i < x.size(); i++) {

		std::cout << x[i].value << std::endl;

	}

//	std::cout << nodes.size() << std::endl;
//	std::cout << "test forces and constraints" << std::endl;
//	std::cout << forceVertexIdSet.size() << std::endl;
//	std::cout << forceVertexIdSet.begin()->first << " " << forceVertexIdSet.begin()->second << std::endl; 
//	std::cout << std::next(forceVertexIdSet.begin(),1)->first <<" "<< std::next(forceVertexIdSet.begin(), 1)->second<< std::endl;
//	std::cout << constraintVertexIdSet.size() << std::endl;



//	for (int i = 0; i < nodes.size(); i++) {

//		nodes[i][0] += u_temp(3 * i + 0) * 100; nodes[i][1] += u_temp(3 * i + 1) * 100; nodes[i][2] += u_temp(3 * i + 2) * 100;

//	}

	makeMesh(x, nodes,mask);

	std::cout << "done making surface mesh" << std::endl;

	//return 1;

	//voxelize("assets/beam.obj", x, 5, 1, nodes);

	std::cout << "done creating voxel mesh" << std::endl;

	bool renderFlag = false; //render or don't render with libigl

	//base material properties (Aluminum)
//	double ni = 0.3;
	//double E = 70000000000;
//	double density = 2800;
	double allowableStress = 320000000;

	//Define and Initialize local stiffness matrix
	Eigen::MatrixXd KE(24, 24); //assumes x[0] value of node is 1 at the start	(?)
	KE = K_mat();

	//makeSurfaceMesh3(x, nodes); //writes mesh to file "mesh.obj"

	//igl::readOBJ("mesh.obj", V, F);

	//constraint and force selection on voxel grid by mouse input

	static bool boolVariable = true;
	static bool planeVariable = true;
	static float forceX = 0;
	static float  forceY = 0;
	static float forceZ = 0;
	
	int optimizerFlag = -1;

	std::set < std::pair<int, Eigen::Vector3f>, comp> forceFaceList;
	std::set<int> constraintfaceList;
	
	viewer.callback_mouse_down =
		[&V, &F, &N, &forceVertexIdSet, &constraintVertexIdSet, &numOfVoxelsX,&numOfVoxelsY,&numOfVoxelsZ,dx,&nodes, &voxelTestFlag, &constraintfaceList, &forceFaceList](igl::opengl::glfw::Viewer& viewer, int, int)->bool
	{
		std::cout <<"test voxelization flag: "<< voxelTestFlag << std::endl;
		if (voxelTestFlag) return false;
		int fid;
		Eigen::Vector3f bc; //holds the hit position in barycentric coordinates
		// Cast a ray in the view direction starting from the mouse position
		double x_Ray = viewer.current_mouse_x;
		double y = viewer.core.viewport(3) - viewer.current_mouse_y;
		if (igl::unproject_onto_mesh(Eigen::Vector2f(x_Ray, y), viewer.core.view,
			viewer.core.proj, viewer.core.viewport, V, F, fid, bc))
		{
			Eigen::MatrixXd V_point(1, 3), endPoint(1, 3), colorPointConstraints(1, 3), colorPointForces(1, 3);
			colorPointConstraints.row(0) = Eigen::Vector3d(1, 0, 0);
			colorPointForces.row(0) = Eigen::Vector3d(0, 1, 0);
			std::pair<int, Eigen::Vector3f> tempForce;

			int selectedV = 0;
			for (int j = 1; j < 3; j++) {
				if (bc(j) > bc(selectedV)) {
					selectedV = j;
				}
			}

			if (bc(selectedV) > 0.6) {

				if (!planeVariable) {

					if (boolVariable) { //make constraints, remove forces or existing constraint
						tempForce.first = F(fid, selectedV);
						tempForce.second = Eigen::Vector3f(forceX, forceY, forceZ);

						if (constraintVertexIdSet.find(tempForce.first) != constraintVertexIdSet.end()) {
							constraintVertexIdSet.erase(tempForce.first);
						}
						else {
							constraintVertexIdSet.insert(tempForce.first);
						}

						if (forceVertexIdSet.find(tempForce) != forceVertexIdSet.end()) {
							forceVertexIdSet.erase(tempForce);
						}

					}
					else { //make forces, remove constraints or existing force

						if (forceVertexIdSet.find(tempForce) != forceVertexIdSet.end()) {
							forceVertexIdSet.erase(tempForce);
						}
						else {
							forceVertexIdSet.insert(tempForce);
						}

						if (constraintVertexIdSet.find(tempForce.first) != constraintVertexIdSet.end()) {
							constraintVertexIdSet.erase(tempForce.first);
						}
					}

				} else {	//select all plane

					Eigen::MatrixXi EV, FE, EF;
					igl::edge_topology(V, F, EV, FE, EF);

					//perform BFS from initial face. define set of visited faces, and a list of faces to visit

					std::list<int> faces_to_visit; std::set<int> visited_faces;

					faces_to_visit.push_back(fid);
					
					while (!faces_to_visit.empty()) {

						int currFace = faces_to_visit.front();

						visited_faces.insert(currFace);

						//add force/constraint to face vertices
				//		for (int v = 0; v < 3; v++) {
				//			tempForce.first = F(currFace, v);
							tempForce.first = currFace;
							tempForce.second = Eigen::Vector3f(forceX, forceY, forceZ);

							if (boolVariable) { //make constraints, remove forces or existing constraint

								constraintfaceList.insert(tempForce.first);

								if (forceFaceList.find(tempForce) != forceFaceList.end()) {
									forceFaceList.erase(tempForce);
								}

							}
							else {

									forceFaceList.insert(tempForce);
					
								if (constraintfaceList.find(tempForce.first) != constraintfaceList.end()) {
									constraintfaceList.erase(tempForce.first);
								}

							}
				//		}
						//add unvisited neighbours to list it change in normal is at most tol

						double tol = 0.8;

						for (int e = 0; e < 3; e++) {

							int nb_face = EF(FE(currFace, e), 0) == currFace ? EF(FE(currFace, e), 1) : EF(FE(currFace, e), 0);

							double dCos_angle = N.row(currFace).dot(N.row(nb_face));

							if (dCos_angle > tol) {

								if (visited_faces.find(nb_face) == visited_faces.end()) {
									std::cout << nb_face << std::endl;
									faces_to_visit.push_back(nb_face);
									visited_faces.insert(nb_face);
								}
							}

						}
						
						faces_to_visit.pop_front();
						
					}
					
				}
			}

			viewer.data().clear();
			viewer.data().set_mesh(V, F);
			viewer.data().show_lines = true;

		}


		Eigen::RowVector3d V_point(3), endPoint(3), colorPointConstraints(3), colorPointForces(3);
		colorPointConstraints.row(0) = Eigen::Vector3d(1, 0, 0);
		colorPointForces.row(0) = Eigen::Vector3d(0, 1, 0);

		std::cout<<"drawing constraints and forces"<<std::endl;
		//draw points and lines (for forces)
		for (auto j = forceVertexIdSet.begin(); j != forceVertexIdSet.end(); j++) {
			viewer.data().add_points(V.row(j->first), colorPointForces);
			double x = j->second(0); double y = j->second(1); double z = j->second(2);
			viewer.data().add_edges(V.row(j->first), V.row(j->first) + Eigen::RowVector3d(x, y, z)*dx, colorPointForces);

		}
		for (auto j = constraintVertexIdSet.begin(); j != constraintVertexIdSet.end(); j++) {
			viewer.data().add_points(V.row(*j), colorPointConstraints);
		//	viewer.data().add_points(Eigen::RowVector3d(nodes[*j][0], nodes[*j][1], nodes[*j][2]), colorPointConstraints);
		}

		for (auto j = constraintfaceList.begin(); j != constraintfaceList.end(); j++) {
			for (int k = 0; k < 3; k++) {
				viewer.data().add_points(V.row(F(*j,k)), colorPointConstraints);
				//	viewer.data().add_points(Eigen::RowVector3d(nodes[*j][0], nodes[*j][1], nodes[*j][2]), colorPointConstraints);
			}
		}

		for (auto j = forceFaceList.begin(); j != forceFaceList.end(); j++) {
			for (int k = 0; k < 3; k++) {
				viewer.data().add_points(V.row(F(j->first, k)), colorPointForces);
				double x = j->second(0); double y = j->second(1); double z = j->second(2);
				viewer.data().add_edges(V.row(F(j->first, k)), V.row(F(j->first, k)) + Eigen::RowVector3d(x, y, z)*dx, colorPointForces);
			}

		}
		// Show mesh
		return false;
		};
	std::cout << R"(Usage:
  [click]  Pick face on shape
)";

	// Attach a menu plugin
	igl::opengl::glfw::imgui::ImGuiMenu menu;

	GLuint tex_2d = 0;
	GLuint tex_2d_1 = 1;
	GLuint tex_2d_2 = 2;

	menu.callback_draw_custom_window = [&]()
	{
		// Define next window position + size
		ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(250, 650), ImGuiSetCond_FirstUseEver);
		ImGui::Begin(
			"Preprocessing", nullptr,
			ImGuiWindowFlags_NoSavedSettings
		);

		double doubleVariable = 0.1f; // Shared between two menus

		// Expose the same variable directly ...
		ImGui::PushItemWidth(-80);
		//ImGui::DragScalar("double", ImGuiDataType_Double, &doubleVariable, 0.1, 0, 0, "%.4f");
		ImGui::DragInt("Nx", &numOfVoxelsX);
		ImGui::DragInt("Ny", &numOfVoxelsY);
		ImGui::DragInt("Nz", &numOfVoxelsZ);
		if (ImGui::Button("Create", ImVec2(100, 50)))
		{

			viewer.data().clear();
			create_mesh(nodes, x, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ, dx);
		//	forceVertexIdSet.clear();
		//	constraintVertexIdSet.clear();
			remove("C:\\Users\\barda\\source\\repos\\Project1\\x64\\Release\\mesh.obj");
			makeSurfaceMesh3(x, nodes); //writes mesh to file "finalObject.obj"
			igl::readOBJ("mesh.obj", V, F);
			viewer.data().set_mesh(V, F);
			igl::per_face_normals(V, F, N);
			std::cout << "create mesh" << std::endl;

		}

		if (ImGui::Button("Open", ImVec2(100, 50))) {

			viewer.data().clear();
			voxelTestFlag = false;
			char filename[MAX_PATH];

			OPENFILENAME ofn;
			ZeroMemory(&filename, sizeof(filename));
			ZeroMemory(&ofn, sizeof(ofn));
			ofn.lStructSize = sizeof(ofn);
			ofn.hwndOwner = NULL;  // If you have a window to center over, put its HANDLE here
			ofn.lpstrFilter = "Obj Files\0*.obj";
			ofn.lpstrFile = filename;
			ofn.nMaxFile = MAX_PATH;
			ofn.lpstrTitle = "Select a File";
			ofn.Flags = OFN_DONTADDTORECENT | OFN_FILEMUSTEXIST;

			if (GetOpenFileNameA(&ofn))
			{
				std::cout << "You chose the file \"" << filename << "\"\n";
			}
			
			//voxelize(filename, x, 5, dx, nodes, V, F);
			forceVertexIdSet.clear();
			constraintVertexIdSet.clear();
			forceFaceList.clear();
			constraintfaceList.clear();
			constraintVertexIdSet.clear();
			//makeMesh(x, nodes, mask); //writes mesh to file "finalObject.obj"
			igl::readOBJ(filename, V, F);
			viewer.data().set_mesh(V, F);
			igl::per_face_normals(V, F, N);
			std::cout << "create mesh" << std::endl;
			viewer.data().set_mesh(V, F);
			
		}

		ImGui::Checkbox("constraint\force", &boolVariable);
		ImGui::Checkbox("select plane", &planeVariable);

		ImGui::PopItemWidth();
		ImGui::DragFloat("x", &forceX);
		ImGui::DragFloat("y", &forceY);
		ImGui::DragFloat("z", &forceZ);
		static std::string str = "bunny";
		//	ImGui::InputText("Name", str);

		igl::png::texture_from_file("fist.png", tex_2d); //works becuase of the opengl context initialized above
		igl::png::texture_from_file("beam.png", tex_2d_1); 
		igl::png::texture_from_file("FEM.png", tex_2d_2);

		ImGui::BulletText("Choose Simulator");
		ImGui::Indent();
		ImGui::ImageButton((void*)tex_2d_2, ImVec2(50, 50), ImVec2(0, 0), ImVec2(1, 1), -1, ImVec4(1, 1, 1, 0), ImVec4(1, 1, 1, 1));
		ImGui::Unindent();
		//for test - can remove
		{
			ImGui::BulletText("Choose Optimizer");
			ImGui::Indent();
			ImGuiIO& io = ImGui::GetIO();
			ImTextureID my_tex_id = io.Fonts->TexID;

			if (ImGui::Button("Compliance MMA Stress Optimizer", ImVec2(200, 25))) {

				optimizerFlag = 2;

				std::cout << "chose stress optimizer" << std::endl;

			}

			if (ImGui::Button("Compliance MMA Fatigue Optimizer", ImVec2(200, 25))) {

				optimizerFlag = 3;

				std::cout << "choose fatigue optimizer" << std::endl;
			}


				ImGui::PushID(0);
			//	if (ImGui::ImageButton((void*)0, ImVec2(50, 50), ImVec2(0, 0), ImVec2(1, 1), -1, ImVec4(1, 1, 1, 0), ImVec4(1, 1, 1, 1))) {
					if (ImGui::Button("NSGAIII Test", ImVec2(200, 25))) {
					optimizerFlag = 0;

					std::vector<std::pair<int, std::vector<double>>> force;

					force.resize((forceVertexIdSet.size()));
					int it = 0;
					for (auto i = forceVertexIdSet.begin(); i != forceVertexIdSet.end(); i++) {
						force[it].second.resize(3);
						force[it].first = i->first;
						force[it].second[0] = i->second(0);
						force[it].second[1] = i->second(1);
						force[it].second[2] = i->second(2);
						it++;
					}

				//	int numofuntunable = turn_untunable(nodes, x, NeighbourLayers, dx, force);

				//	calc_rho(nodes, x, NeighbourLayers, dx, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ,r0);

					std::cout << "pressed stress optimizer" << std::endl;
				}
				ImGui::PopID();
				ImGui::PushID(1);
				//	if (ImGui::ImageButton((void*)1, ImVec2(50, 50), ImVec2(0, 0), ImVec2(1, 1), -1, ImVec4(1, 1, 1, 0), ImVec4(1, 1, 1, 1))) {
				if (ImGui::Button("Compliance MMA optimizer", ImVec2(200, 25))) {

						optimizerFlag = 1;
						std::cout << "pressed compliance optimizer" << std::endl;
					}
				ImGui::PopID();

			ImGui::Unindent();
		}

		if (ImGui::Button("Test Voxelization", ImVec2(100, 50))) {

			voxelTestFlag = true;
			viewer.data().clear();
			viewer.data_list.resize(2);
			viewer.data_list[0].set_mesh(V, F);
			voxelize(x, 10, dx, nodes, V, F, constraintVertexIdSet, forceVertexIdSet, constraintfaceList, forceFaceList, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ);
		//	getBCforVoxelTest(x, nodes, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ, forceVertexIdSet, constraintVertexIdSet);
			makeMesh(x, nodes, mask);
			igl::readOBJ("finalObject.obj", V, F);
			viewer.data().set_mesh(V, F);
			viewer.data_list[1].set_mesh(V, F);
			Eigen::RowVector3d colorPointConstraints(1, 0, 0);
			Eigen::RowVector3d colorPointForces(0, 1, 0);
			//	std::cout << "number of nodes: " << nodes.size() << " " << nodes[1][0] << nodes[1][1] << nodes[1][2] <<std::endl;
			for (auto j = constraintVertexIdSet.begin(); j != constraintVertexIdSet.end(); j++) {
				std::cout << *j << std::endl;
				viewer.data().add_points(Eigen::RowVector3d(nodes[*j][0], nodes[*j][1], nodes[*j][2]), colorPointConstraints);
			}
			std::cout << "number of forces: " << forceVertexIdSet.size() << std::endl;
			for (auto j = forceVertexIdSet.begin(); j != forceVertexIdSet.end(); j++) {
				double x = j->second(0); double y = j->second(1); double z = j->second(2);
				viewer.data().add_points(Eigen::RowVector3d(nodes[j->first][0], nodes[j->first][1], nodes[j->first][2]), colorPointForces);
				viewer.data().add_edges(Eigen::RowVector3d(nodes[j->first][0], nodes[j->first][1], nodes[j->first][2]), Eigen::RowVector3d(nodes[j->first][0], nodes[j->first][1], nodes[j->first][2]) + Eigen::RowVector3d(x, y, z)*dx, colorPointForces);

			}
		}

		if (ImGui::Button("Go", ImVec2(50, 50))) {

			if (viewer.core.is_animating == true) {
				viewer.core.is_animating = false;
			}
			else {
				viewer.core.is_animating = true;
			}
		}

		if (ImGui::Button("Clear", ImVec2(50, 50))) {
			viewer.data().clear();
			x.clear();
			nodes.clear();
			//create_mesh(nodes, x, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ,dx);
			forceVertexIdSet.clear();
			constraintVertexIdSet.clear();
			remove("C:\\Users\\barda\\source\\repos\\Project1\\x64\\Release\\mesh.obj");
			//makeSurfaceMesh3(x, nodes); //writes mesh to file "mesh.obj"
			//igl::readOBJ("mesh.obj", V, F);
			//viewer.data().set_mesh(V, F);
		}
		ImGui::End();
	};
	
	viewer.plugins.push_back(&menu);

	//double change = 1000;

	//int power = 3;

	//initialize for compliance optimizer
	
	//double* dc = new double[x.size()];
	//double V_allowed = 0.3*numOfVoxelsX*numOfVoxelsY*numOfVoxelsZ;
	/*
	double* xnew = new double[x.size()];
	double* dg = new double[x.size()];
	double* g = new double[1];
	double *xmin = new double[x.size()];
	double *xmax = new double[x.size()];
	double* dc_filtered = new double[x.size()];
	MMASolver* mma = new MMASolver(x.size(), 1);
	


	g[0] = -V_allowed;
	for (int i = 0; i < x.size(); i++) {
		dg[i] = 1;
		xmin[i] = 0.001;
		xmax[i] = 1;
		g[0] += x[i].value;
		xnew[i] = x[i].value;
	}
	*/
	bool animation_flag = false;

	//initialize for stress optimizer
	//int C = 2; //number of stress groups
	//double* df = new double[x.size()];
	//double* xnew = new double[x.size()];
	//double* dg = new double[x.size()*C];
	//double* g = new double[C];
	//double *xmin = new double[x.size()];
	//double *xmax = new double[x.size()];
	//double* dc_filtered = new double[x.size()];
	//MMASolver* mma = new MMASolver(x.size(), 1);

	//for (int i = 0; i < x.size(); i++) {
	//	xmin[i] = 0.001;
	//	xmax[i] = 1;
	//	xnew[i] = x[i].value;
	//	for (int j = 0; j < C; j++) {
	//		dg[i*C + j]=0;
	//	}
	//}
	
	viewer.core.is_animating = false;

	int curr_time_step = 0;

	viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer &)->bool
	{

		if (!viewer.core.is_animating) {
			return false;
		}
		
		if (optimizerFlag == 0) { //test NSGA3

			// To use genetic NSGA III genetic algorithm:
			// 1. specify INI file
			// 2. Specify porblem to solve

			std::cout << "testing NSGAIII" << std::endl;
			std::pair<int, Eigen::Vector3f> tempForce;

			tempForce.first = numOfVoxelsX * numOfVoxelsY * numOfVoxelsZ - 1;
			tempForce.second = Eigen::Vector3f(20, 0, 0);
			forceVertexIdSet2.insert(tempForce);
			tempForce.first = numOfVoxelsX * numOfVoxelsY * numOfVoxelsZ - 2;
			tempForce.second = Eigen::Vector3f(0, -20, 0);
		//	forceVertexIdSet2.insert(tempForce);

			 
			if (debug_flag) {
				std::cout << "testing test forces " << forceVertexIdSet2.size() << std::endl;

				for (auto i = forceVertexIdSet2.begin(); i != forceVertexIdSet2.end(); i++) { //works

					int id = i->first;
					std::cout << "testing test forces" << std::endl;
					viewer.data().add_points(Eigen::RowVector3d(nodes[id][0], nodes[id][1], nodes[id][2]), Eigen::RowVector3d(1, 0.5, 1)); //test passed
				}
			}

			//define NSGAIII
			CNSGAIII nsgaiii;
			ComplianceProblem *problem = 0;
			std::ifstream exp_ini("ComplianceTest.ini");
			nsgaiii.Setup(exp_ini);
			problem = new ComplianceProblem(2,x.size(),"ComplianceTest",x,nodes, forceVertexIdSet
			,dx, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ, 0.3, constraintVertexIdSet); //define the problem to be solved (variables, num of objectives, evaluation function)
			std::ofstream output("output.txt");
			for (size_t r = 0; r < 100; r += 1)
			{
				std::cout << "run: " << r << std::endl;
				CPopulation solutions;
				nsgaiii.Solve(&solutions, *problem);

				for (int i=0; i<solutions.size(); i++) {
				
					output <<r<<" "<< solutions[i].objs()[0] << " " << solutions[i].objs()[1]<<std::endl;

				}

			}

			output.close();
			std::cout << "passed NSGAIII" << std::endl;

			Eigen::MatrixXd color;
		//	for (int i = 0; i < solutions.size(); i++) {

			int best_fit_ind = 0;
			/*
			std::cout << solutions.size() << std::endl;

			for (int i = 0; i < solutions.size(); i++) {

				if (solutions[i].objs()[0] < solutions[best_fit_ind].objs()[0]) {

				//	double sum = 0;

				//	for (int j = 0; j < solutions[i].vars().size(); j++) {
				//		sum += solutions[i].vars()[j];
				//	}
				//	std::cout << sum << std::endl;
					best_fit_ind = i;

				}

			}
			
			for (int i = 0; i < x.size(); i++) {

				if (debug_flag) {
					std::cout << "value in cell " << i << ": " << solutions[best_fit_ind].vars()[i] << std::endl;
				}

				for (int j = 0; j < 12; j++) {
					rhoFaces(i * 12 + j) = solutions[best_fit_ind].vars()[i];
				}
			}*/
		//	igl::jet(rhoFaces, true, color);

		//	std::cout << "fitness: " << solutions[best_fit_ind].objs()[0];

			delete problem;

			viewer.data().set_colors(color);

		//	}

			viewer.core.is_animating = false;
			return false;

		}
		
		if (optimizerFlag == 1) {

		//	voxelize(x, 10, dx, nodes, V/10, F, constraintVertexIdSet, forceVertexIdSet, constraintfaceList, forceFaceList, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ);

			create_mesh(nodes, x, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ, dx);

			getBCforVoxelTest(x, nodes, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ,forceVertexIdSet, constraintVertexIdSet);

			std::cout << "number of voxels: " << x.size() << std::endl;

			double start_weight=0;
			
			int num_untunable=0;

			//initialize for compliance optimizer


			
			for (int i = 0; i < x.size(); i++) {
				
				if (x[i].tunable) {
					start_weight += x[i].value;
				}
				else {
					num_untunable++;
				}

			}

			std::cout << "number of untunable: " << num_untunable << std::endl;

			int numOfParameters = x.size() - num_untunable;
			double* dc = new double[x.size()- num_untunable];
			double V_allowed = 0.5*start_weight;
			double* xnew = new double[x.size()- num_untunable];
			double* dg = new double[(x.size()- num_untunable )*C];
			double* df = new double[x.size()- num_untunable];
			double* g = new double[C];
			double *xmin = new double[x.size()- num_untunable];
			double *xmax = new double[x.size()- num_untunable];
			double* dc_filtered = new double[x.size()- num_untunable];
			MMASolver* mma = new MMASolver(x.size()- num_untunable, 1);

			double density = 2800;

			double ni = 0.3;

			double check_weight=0;


			g[0] = -V_allowed+ start_weight;
			for (int i = 0; i < numOfParameters; i++) {

				for (int j = 0; j < C; j++) {
			//		g[j] = 0;
			//		dg[i*C + j] = 1;
				}
				xmin[i] = 0.001;
				xmax[i] = 1;
				dg[i] = 1;
				xnew[i] = 1;
				df[i] = dx * dx*dx*density;
			}

		//	std::cout <<"V_allowed: "<< V_allowed<< " start_weight: " << start_weight << " check weight: " <<  check_weight << std::endl;

			std::cout << "dx: " <<dx<< std::endl;

	//		std::cout << "constraint size: " << constraintVertexIdSet.size() << std::endl;

	//		std::cout << "force size: " << forceVertexIdSet.size() << std::endl;

			for (int iter = 0; iter < 100; iter++) {

				double goal_test = 0;

				calc_rho(nodes, x, NeighbourLayers, dx, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ, r0);

		//		for (int i = 0; i < x.size(); i++) {

		//			if (x[i].tunable) {
		//				std::cout << "x[" << i << "]: " << x[i].rho << std::endl;
			//		}

			//	}

		//		exit(-1);

				Eigen::VectorXd u = doFEM(nodes, x, forceVertexIdSet, constraintVertexIdSet, dx, E, power);

				for (auto f = forceVertexIdSet.begin(); f != forceVertexIdSet.end(); f++) {

					int id = f->first;
					goal_test += u[3 * id + 0] * f->second(0)+ u[3 * id + 1] * f->second(1)+ u[3 * id + 2] * f->second(2);

				}

				std::cout << "iter: " << iter << " constraint: " << g[0] << " goal: " <<goal_test<< std::endl;

				std::cout << "calculated u, min u: " << u.minCoeff() << std::endl;



				//	igl::jet(rhoFaces, true, Color);
			//	viewer.data().set_mesh(V, F);
			//	viewer.data().set_colors(Color);

				Eigen::VectorXd displacement(24);
				int ind=0;
				for (int i = 0; i < x.size(); i++) {

					if (x[i].tunable == false ) {
						continue;
					}

					displacement(0) = u[x[i].nodeIdx[0] * 3]; displacement(1) = u[x[i].nodeIdx[0] * 3 + 1]; displacement(2) = u[x[i].nodeIdx[0] * 3 + 2];
					displacement(3) = u[x[i].nodeIdx[1] * 3]; displacement(4) = u[x[i].nodeIdx[1] * 3 + 1]; displacement(5) = u[x[i].nodeIdx[1] * 3 + 2];
					displacement(6) = u[x[i].nodeIdx[2] * 3]; displacement(7) = u[x[i].nodeIdx[2] * 3 + 1]; displacement(8) = u[x[i].nodeIdx[2] * 3 + 2];
					displacement(9) = u[x[i].nodeIdx[4] * 3]; displacement(10) = u[x[i].nodeIdx[4] * 3 + 1]; displacement(11) = u[x[i].nodeIdx[4] * 3 + 2];
					displacement(12) = u[x[i].nodeIdx[3] * 3]; displacement(13) = u[x[i].nodeIdx[3] * 3 + 1]; displacement(14) = u[x[i].nodeIdx[3] * 3 + 2];
					displacement(15) = u[x[i].nodeIdx[6] * 3]; displacement(16) = u[x[i].nodeIdx[6] * 3 + 1]; displacement(17) = u[x[i].nodeIdx[6] * 3 + 2];
					displacement(18) = u[x[i].nodeIdx[5] * 3]; displacement(19) = u[x[i].nodeIdx[5] * 3 + 1]; displacement(20) = u[x[i].nodeIdx[5] * 3 + 2];
					displacement(21) = u[x[i].nodeIdx[7] * 3]; displacement(22) = u[x[i].nodeIdx[7] * 3 + 1]; displacement(23) = u[x[i].nodeIdx[7] * 3 + 2];

					dc[ind] = -dx*power * std::pow(x[i].value, power - 1)* E*displacement.transpose()*KE*displacement;
					ind++;

				}

				change = optimizeCompliance(xnew, dg, g, dc, xmin, xmax, dc_filtered, V_allowed, mma, nodes, x, dx, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ, 1);

	//			if (iter == 49) {
	//				for (int i = 0; i < nodes.size(); i++) {
	//					nodes[i][0] += u(3 * i + 0) * 0.000001; nodes[i][1] += u(3 * i + 1) * 0.000001; nodes[i][2] += u(3 * i + 2) * 0.000001;
	//				}
	//			}

				std::cout << "change: " << change << std::endl;
			}
			
			std::cout << "x.size(): " <<x.size()<< std::endl;

			for (int i = 0; i < x.size(); i++) {

				if (x[i].tunable) {
					std::cout << "x[" << i << "]: " << x[i].value << std::endl;
				}

			}

				makeMesh(x, nodes,mask);

				std::cout << "done" << std::endl;

				viewer.core.is_animating = false;
				return false;


		}

		if (optimizerFlag == 2) { //test stress optimizer

		//	calc_rho(nodes, x, -1, dx, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ, r0);

			std::cout << "begin stress" << std::endl;

				auto constraintVertexIdSet2 = constraintVertexIdSet;
				forceVertexIdSet2 = forceVertexIdSet;

				if (!mask.empty()) {
				constraintVertexIdSet2.clear();
				forceVertexIdSet2.clear();

				for (auto i = forceVertexIdSet.begin(); i != forceVertexIdSet.end(); i++) {
					forceVertexIdSet2.insert(std::pair<int, Eigen::Vector3f>(mask[i->first], i->second));
				}
				for (auto i = constraintVertexIdSet.begin(); i != constraintVertexIdSet.end(); i++) {
					constraintVertexIdSet2.insert(mask[*i]);
				}
			}

			Eigen::VectorXd u = doFEM(nodes, x, forceVertexIdSet2, constraintVertexIdSet2, dx, E,power);

			for (int i = 0; i < nodes.size(); i++) {
				nodes[i][0] += u(3 * i + 0)*100000; nodes[i][1] += u(3 * i + 1)*100000; nodes[i][2] += u(3 * i + 2)*100000;
			}

			if (!mask.empty()) {

				makeMesh(x, nodes, mask);

				igl::readOBJ("finalObject.obj", V, F);

			}
			
			std::cout << "max displacement: " << u.maxCoeff() << std::endl;

			//auto StressVec = calcStresses(ni, E, nodes, x, u, dx);


			//optimizeStress(StressVec, u, C, power, dx, change, nodes, x, xnew, dg, g, df, xmin, xmax, density, E, ni, mma, 1, r0);

//			rhoFaces.resize(F.rows());
			std::cout << F.rows() << std::endl;
			int ind = 0;
			for (int i = 0; i < x.size(); i++) {
				if (x[i].value>0.001) {
//				std::cout <<  "stress at " <<i<<" : "<< calcVonMises(StressVec[i].Stresses) << std::endl;
				for (int j = 0; j < 12; j++) {
	//				rhoFaces(ind * 12 + j) = x[i].value;//calcVonMises(StressVec[i].Stresses);
					}
				ind++;
				}
			}

		//	igl::jet(rhoFaces, true, Color);
			viewer.data().set_mesh(V, F);
			viewer.data().set_colors(Color);

			Eigen::VectorXd displacement(24);

		}


		if (optimizerFlag == 3) { //test fatigue optimizer

			/*
			int steps = 60;

			if (!animation_flag) {

				//	calc_rho(nodes, x, -1, dx, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ, r0);

			//		Eigen::VectorXd u = doFEM(nodes, x, forceVertexIdSet, constraintVertexIdSet, dx, 1);

			//		std::cout << t << std::endl;

				Eigen::VectorXd u_dot_0, u0;

				u_dot_0 = Eigen::VectorXd::Zero(nodes.size());

				u0 = Eigen::VectorXd::Zero(nodes.size());

				std::vector<int> targets = linspace<int>(0, x.size() - 1, x.size());

				std::cout << "targets size: " << targets.size() << std::endl;

				std::vector<std::vector<double>> sigma_m;

				std::vector<std::vector<Stress>> sigma_a;

				std::vector<std::vector<Eigen::VectorXd>> u_max(targets.size());

				std::vector<std::vector<Eigen::VectorXd>> u_min(targets.size());

				ClacStressTimeHistory(targets, x, 0.001, 0.001, power, forceVertexIdSet, dx, E, nodes, constraintVertexIdSet, u_dot_0, u0, 0.0001, steps, sigma_m, sigma_a, disp_test,u_max,u_min);

				std::cout << "number of u_max: " << u_max[3].size() <<"/"<< sigma_a[3].size() << std::endl;

			//	std::cout << "calculated stress histories" << std::endl;
				if (debug_flag) {
					for (int i = 0; i < targets.size(); i++) {

						for (int j = 0; j < sigma_m[i].size(); j++) {

							std::cout << "i: " << i << " j: " << j << " sigma_m: " << sigma_m[i][j] << " sigma_a: " << sigma_a[i][j].vonMises << std::endl;

						}
					}
				}

				double bf = -0.075;
				double sigma_f = 650000000;
				double ni = 1000000;
				double Sut = 380000000;
				double Sf = 3;

				std::vector<double> damage(x.size(),0);

				for (int i = 0; i < targets.size(); i++) {

					std::vector<double> N(sigma_m[i].size());
					//calculate damage for target
					double D = 0;
					for (int j = 0; j < N.size(); j++) {
						if (((sigma_f*(Sut - std::max(sigma_m[i][j], 0.0)))) < 0) {
							D = 1;
							break;
						}
						D += ni*2*std::pow(Sut*sigma_a[i][j].vonMises / (sigma_f*(Sut - std::max(sigma_m[i][j], 0.0))), -1.0 / bf);
						if (D >= 1) {
							D = 1;
							break;
						}
						//	std::cout << D << std::endl;
					}

					x[targets[i]].rho = D;//sigma_a[i][0].vonMises;
					damage[i] = D;

	
				}


				//need to split the damages in to groups

				int K = 2; //number of damage groups

				auto sorted = sort_indexes(damage);

				std::vector<std::vector<int>> Omega_k;

				for (int i = 0; i < K; i++) {

					Omega_k.push_back(std::vector<int>(sorted.begin() + std::ceil((float)i / K * sorted.size()), sorted.begin() + std::ceil((float)(i+1)/K * sorted.size())));

				//	std::cout<<Omega_k[i].size() << std::endl;
				}

				
				std::vector<int> n(sigma_a.size(), ni);

				std::vector<double> derivs = full_dD_PN_k_dgamma_e(x, Omega_k, n, power, u_max, u_min, sigma_a, damage, Sut, bf, Sf,sigma_m, E, nodes, r0);


				for (int i = 0; i < x.size(); i++) {
					std::cout << derivs[i] << std::endl;
			//		std::cout << "element: " << i << ", Damage: " << x[i].rho << std::endl;
					for (int j = 0; j < 12; j++) {
						rhoFaces(i * 12 + j) = derivs[i];//x[i].rho;

					}
				}

				igl::jet(rhoFaces, true, Color);
	//			viewer.data().set_mesh(V, F);
				viewer.data().set_colors(Color);

				animation_flag = false; //change to true for animation

			}

			std::ofstream writer("matrix.txt");

			writer << disp_test;

			if (animation_flag) {

				// do animation
				auto tempNodes = nodes;
				double factor = 1;
				for (int i = 0; i < tempNodes.size(); i++) {
					tempNodes[i][0] = nodes[i][0] + disp_test(curr_time_step, 3 * i)*factor;
					tempNodes[i][1] = nodes[i][1] + disp_test(curr_time_step, 3 * i + 1)*factor;
					tempNodes[i][2] = nodes[i][2] + disp_test(curr_time_step, 3 * i + 2)*factor;
				}

			//	std::cout << "max movement: " << disp_test.row(curr_time_step).cwiseAbs().maxCoeff() << std::endl;
				makeSurfaceMesh3(x, tempNodes); //writes mesh to file "mesh.obj"
				igl::readOBJ("mesh.obj", V, F);
				viewer.data().set_mesh(V, F);
				curr_time_step++;
				if (curr_time_step < steps) {
					viewer.core.is_animating = true;
				}
				else {
					curr_time_step = 0;
					animation_flag = false;
					viewer.core.is_animating = false;
				}
				return false;

			}
			
			else {



			}
			*/
		}
			
		viewer.core.is_animating = false;
		return false;

		};

	viewer.launch();
	
			return 0;
		}

