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

	std::vector<int> mask;

	int argc = 1;
	char *argv[1] = { (char*)"Something" };

	//initialize voxel parameters
	int numOfVoxelsX =30;
	int numOfVoxelsY = 10;
	int numOfVoxelsZ = 1;
	float dx =1;//0.01;// 0.001; // size
	int NeighbourLayers = 1;
	float r0 = 1.5*dx*NeighbourLayers;

	//declare geometry parameters
	std::vector<designVariable> x;
	std::vector<std::vector<double>> nodes; //create array that connects node indices to coordinates.
	nodes.resize((numOfVoxelsX + 1) *(numOfVoxelsY + 1)*(numOfVoxelsZ + 1));

	create_mesh(nodes, x,  numOfVoxelsX,  numOfVoxelsY,  numOfVoxelsZ,dx);

	//voxelize("assets/beam.obj", x, 5, 1, nodes);

	std::cout << "done creating voxel mesh" << std::endl;

	

	//base material properties (Aluminum)
	double ni = 0.3;
	double E = 70000000000;
	double density = 2800;
	double allowableStress = 320000000;

	//Define and Initialize local stiffness matrix
	Eigen::MatrixXd KE(24, 24); //assumes x[0] value of node is 1 at the start	(?)
	KE = K_mat();

	makeSurfaceMesh3(x, nodes); //writes mesh to file "mesh.obj"

	//initialize matrices for faces and vertices - for libigl display
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	//igl::readOBJ("mesh.obj", V, F);

	//constraint and force selection on voxel grid by mouse input

	std::set < std::pair<int, Eigen::Vector3f>, comp> forceVertexIdSet;

	std::set < std::pair<int, Eigen::Vector3f>, comp> forceVertexIdSet2; //for NSGAIII test

	std::set<int> constraintVertexIdSet;
	
	static bool boolVariable = true;
	static bool planeVariable = true;
	static float forceX = 0;
	static float  forceY = 0;
	static float forceZ = 0;
	
	igl::opengl::glfw::Viewer viewer; //sets opengl context (important for images of buttons)

	int optimizerFlag = -1;
	
	viewer.callback_mouse_down =
		[&V, &F, &forceVertexIdSet, &constraintVertexIdSet, &numOfVoxelsX,&numOfVoxelsY,&numOfVoxelsZ, &x](igl::opengl::glfw::Viewer& viewer, int, int)->bool
	{
		int fid;
		Eigen::Vector3f bc; //holds the hit position in barycentric coordinates
		// Cast a ray in the view direction starting from the mouse position
		double x_Ray = viewer.current_mouse_x;
		double y = viewer.core.viewport(3) - viewer.current_mouse_y;
		if (igl::unproject_onto_mesh(Eigen::Vector2f(x_Ray, y), viewer.core.view,
			viewer.core.proj, viewer.core.viewport, V, F, fid, bc))
		{
			std::cout << V.row(F(fid, 1)) << std::endl;
			Eigen::MatrixXd V_point(1, 3), endPoint(1, 3), colorPointConstraints(1, 3), colorPointForces(1, 3);
			colorPointConstraints.row(0) = Eigen::Vector3d(1, 0, 0);
			colorPointForces.row(0) = Eigen::Vector3d(0, 1, 0);

			//if vertex select mode - which of the three triangle nodes is closest to click site?
			int selectedV = 0;
			for (int j = 1; j < 3; j++) {
				if (bc(j) > bc(selectedV)) {
					selectedV = j;
				}
			}

			//if plane select mode - check which plane
			int planeInd=-1;
			if (planeVariable) {
				int v1 = F(fid, 0); int v2 = F(fid, 1); int v3 = F(fid, 2); //negative yz plane
				if ((V(v1, 0) == V(v2, 0)) && (V(v1, 0) == V(v3, 0)) && V(v1, 0) == 0.0) {
					planeInd = 0;
				}
				if ((V(v1, 0) == V(v2, 0)) && (V(v1, 0) == V(v3, 0)) && V(v1, 0) > 0.0) { //positive yz plane
					planeInd = 1;
				}
			}
			if (bc(selectedV) > 0.6 || planeInd > -1) {
				viewer.data().clear();
				std::pair<int, Eigen::Vector3f> tempForce;
				int numOfVoxelsinPlane=0;
				int startVoxel = -1;
				int nodecounter = 0;
				std::vector<int> vertexList;
				if (planeVariable) {

					if (planeInd == 0) {
						startVoxel = 0;
						numOfVoxelsinPlane = numOfVoxelsY * numOfVoxelsZ;
						nodecounter = 4;
						vertexList = {0,2,3,5};
						
					}
					if (planeInd == 1) {
						startVoxel = numOfVoxelsX - 1;
						numOfVoxelsinPlane = numOfVoxelsY * numOfVoxelsZ;
						nodecounter = 4;
						vertexList = { 1,4,6,7 };
					}
				}
				else {
					std::cout <<"testing select: " <<fid<<" "<<selectedV<<" "<< F(fid, selectedV) <<" "<<F.rows()<< " "<<V.rows()<< std::endl;
					tempForce.first = F(fid, selectedV);
					tempForce.second = Eigen::Vector3f(forceX, forceY, forceZ);
					numOfVoxelsinPlane = 1;
					nodecounter = 1;
					std::cout << "testing end:"  << std::endl;
				}

				int currVoxel = startVoxel;

				while (numOfVoxelsinPlane) { //hack for selecting plane or individual vertices

					std::cout << "current voxel is: " << currVoxel << std::endl;

					while (nodecounter>0) {

						if (planeVariable) {
							tempForce.first = x[currVoxel].nodeIdx[vertexList[nodecounter - 1]];
							tempForce.second = Eigen::Vector3f(forceX, forceY, forceZ);
						}
						nodecounter--;

						if (boolVariable) { //constraints are on - delete forces
							std::cout << "got here" << std::endl;
							if (constraintVertexIdSet.find(tempForce.first) == constraintVertexIdSet.end()) {
								constraintVertexIdSet.insert(tempForce.first);
								forceVertexIdSet.erase(tempForce);
								std::cout << "got here1" << std::endl;
							}
							else {
								if (!planeVariable) {
									constraintVertexIdSet.erase(tempForce.first);
								}
							}

						}
						else {
							if (forceVertexIdSet.find(tempForce) == forceVertexIdSet.end()) {
								std::cout << "got here2" << std::endl;
								constraintVertexIdSet.erase(tempForce.first);
								forceVertexIdSet.insert(tempForce);
								std::cout << "got here3" << std::endl;
							}
							else {
								if (!planeVariable) {
									forceVertexIdSet.erase(tempForce);
								}
							}

						}

						std::cout << nodecounter << std::endl;

					}

					numOfVoxelsinPlane--;
					if (planeVariable) {
						currVoxel = moveY(x, currVoxel, 1);
						if (currVoxel == -1) {
							startVoxel = moveZ(x, startVoxel, 1);
							currVoxel = startVoxel;

						}
					}
					nodecounter = 4;

					std::cout << "got to end of loop" << std::endl;

				}

				//draw points and lines (for forces)
				for (auto j = forceVertexIdSet.begin(); j != forceVertexIdSet.end(); j++) {
					V_point.row(0) = V.row(j->first);
					endPoint.row(0) = V.row(j->first) + j->second.cast<double>().transpose().normalized();
					viewer.data().add_points(V_point, colorPointForces);
					viewer.data().add_edges(V_point, endPoint, colorPointForces);
				}
				for (auto j = constraintVertexIdSet.begin(); j != constraintVertexIdSet.end(); j++) {
					V_point.row(0) = V.row(*j);
					viewer.data().add_points(V_point, colorPointConstraints);
				}
			}
			// Show mesh
			viewer.data().set_mesh(V, F);
			viewer.data().show_lines = true;


			return true;
		}
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
			forceVertexIdSet.clear();
			constraintVertexIdSet.clear();
			remove("C:\\Users\\barda\\source\\repos\\Project1\\x64\\Release\\mesh.obj");
			makeSurfaceMesh3(x, nodes); //writes mesh to file "finalObject.obj"
			igl::readOBJ("mesh.obj", V, F);
			viewer.data().set_mesh(V, F);
			std::cout << "create mesh" << std::endl;

		}

		if (ImGui::Button("Open", ImVec2(100, 50))) {

			viewer.data().clear();

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
			
			voxelize(filename, x, 10, dx, nodes, V, F);
			forceVertexIdSet.clear();
			constraintVertexIdSet.clear();
			makeMesh(x, nodes, mask); //writes mesh to file "finalObject.obj"
			igl::readOBJ("finalObject.obj", V, F);
			viewer.data().set_mesh(V, F);
			std::cout << "create mesh" << std::endl;

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

				std::cout << "chose fatigue optimizer" << std::endl;
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

	double change = 1000;

	int power = 3;

	//initialize for compliance optimizer
	
	double* dc = new double[x.size()];
	
	double* xnew = new double[x.size()];
	double* dg = new double[x.size()];
	double* g = new double[1];
	double *xmin = new double[x.size()];
	double *xmax = new double[x.size()];
	double* dc_filtered = new double[x.size()];
	MMASolver* mma = new MMASolver(x.size(), 1);
	
	double V_allowed = 0.3*numOfVoxelsX*numOfVoxelsY*numOfVoxelsZ;

	g[0] = -V_allowed;
	for (int i = 0; i < x.size(); i++) {
		dg[i] = 1;
		xmin[i] = 0.001;
		xmax[i] = 1;
		g[0] += x[i].value;
		xnew[i] = x[i].value;
	}

	bool animation_flag = false;

	Eigen::MatrixXd disp_test;
	/*
	//initialize for stress optimizer
	int C = 2; //number of stress groups
	double* df = new double[x.size()];
	double* xnew = new double[x.size()];
	double* dg = new double[x.size()*C];
	double* g = new double[C];
	double *xmin = new double[x.size()];
	double *xmax = new double[x.size()];
	double* dc_filtered = new double[x.size()];
	MMASolver* mma = new MMASolver(x.size(), 1);

	for (int i = 0; i < x.size(); i++) {
		xmin[i] = 0.001;
		xmax[i] = 1;
		xnew[i] = x[i].value;
		for (int j = 0; j < C; j++) {
			dg[i*C + j]=0;
		}
	}
	*/
	viewer.core.is_animating = false;

	int iter2 = 0;

	int curr_time_step = 0;

	viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer &)->bool
	{

		if (!viewer.core.is_animating) {
			return false;
		}

		Eigen::MatrixXd Color;

		Eigen::VectorXd rhoFaces(F.rows());
		
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
			igl::jet(rhoFaces, true, color);

		//	std::cout << "fitness: " << solutions[best_fit_ind].objs()[0];

			delete problem;

			viewer.data().set_colors(color);

		//	}

			viewer.core.is_animating = false;
			return false;

		}
		
		if (optimizerFlag == 1) {

			std::cout << "optimizing compliance" << std::endl;

			calc_rho(nodes, x, -1, dx, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ,r0);

			Eigen::VectorXd u = doFEM(nodes, x, forceVertexIdSet, constraintVertexIdSet,dx,1);

			for (int i = 0; i < x.size(); i++) {
				for (int j = 0; j < 12; j++) {
					rhoFaces(i * 12 + j) = x[i].rho;
				}
			}

			igl::jet(rhoFaces, true, Color);
			viewer.data().set_mesh(V, F);
			viewer.data().set_colors(Color);

			Eigen::VectorXd displacement(24);
			for (int i = 0; i < x.size(); i++) {

				displacement(0) = u[x[i].nodeIdx[0] * 3]; displacement(1) = u[x[i].nodeIdx[0] * 3 + 1]; displacement(2) = u[x[i].nodeIdx[0] * 3 + 2];
				displacement(3) = u[x[i].nodeIdx[1] * 3]; displacement(4) = u[x[i].nodeIdx[1] * 3 + 1]; displacement(5) = u[x[i].nodeIdx[1] * 3 + 2];
				displacement(6) = u[x[i].nodeIdx[2] * 3]; displacement(7) = u[x[i].nodeIdx[2] * 3 + 1]; displacement(8) = u[x[i].nodeIdx[2] * 3 + 2];
				displacement(9) = u[x[i].nodeIdx[4] * 3]; displacement(10) = u[x[i].nodeIdx[4] * 3 + 1]; displacement(11) = u[x[i].nodeIdx[4] * 3 + 2];
				displacement(12) = u[x[i].nodeIdx[3] * 3]; displacement(13) = u[x[i].nodeIdx[3] * 3 + 1]; displacement(14) = u[x[i].nodeIdx[3] * 3 + 2];
				displacement(15) = u[x[i].nodeIdx[6] * 3]; displacement(16) = u[x[i].nodeIdx[6] * 3 + 1]; displacement(17) = u[x[i].nodeIdx[6] * 3 + 2];
				displacement(18) = u[x[i].nodeIdx[5] * 3]; displacement(19) = u[x[i].nodeIdx[5] * 3 + 1]; displacement(20) = u[x[i].nodeIdx[5] * 3 + 2];
				displacement(21) = u[x[i].nodeIdx[7] * 3]; displacement(22) = u[x[i].nodeIdx[7] * 3 + 1]; displacement(23) = u[x[i].nodeIdx[7] * 3 + 2];

				dc[i] = -power * std::pow(x[i].value, power - 1)* displacement.transpose()*KE*displacement;

			}

			change = optimizeCompliance(xnew, dg, g, dc, xmin, xmax, dc_filtered, V_allowed, mma, nodes, x, dx, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ,1);

			if (iter2 > 30) {
				makeMesh(x, nodes,mask);
				viewer.data().clear();
				igl::readOBJ("finalObject.obj", V, F);
				viewer.data().set_mesh(V, F);
				viewer.core.is_animating = false;
				iter2 = 0;
				return false;
			}
			iter2++;
		}

		if (optimizerFlag == 2) { //test stress optimizer

		//	calc_rho(nodes, x, -1, dx, numOfVoxelsX, numOfVoxelsY, numOfVoxelsZ, r0);

			std::cout << "begin stress" << std::endl;

			auto constraintVertexIdSet2 = constraintVertexIdSet;
			constraintVertexIdSet2.clear();
			forceVertexIdSet2.clear();

			for (auto i = forceVertexIdSet.begin(); i != forceVertexIdSet.end(); i++) {
				forceVertexIdSet2.insert(std::pair<int,Eigen::Vector3f>(mask[i->first],i->second));
			}
			for (auto i = constraintVertexIdSet.begin(); i != constraintVertexIdSet.end(); i++) {
				constraintVertexIdSet2.insert(mask[*i]);
			}

			Eigen::VectorXd u = doFEM(nodes, x, forceVertexIdSet2, constraintVertexIdSet2, dx, E);

			for (int i = 0; i < nodes.size(); i++) {
				nodes[i][0] += u(3 * i + 0)*100000; nodes[i][1] += u(3 * i + 1)*100000; nodes[i][2] += u(3 * i + 2)*100000;
			}

			makeMesh(x, nodes, mask);

			igl::readOBJ("finalObject.obj",V,F);
			
			std::cout << "max displacement: " << u.maxCoeff() << std::endl;

			auto StressVec = calcStresses(ni, E, nodes, x, u, dx);
			rhoFaces.resize(F.rows());
			std::cout << F.rows() << std::endl;
			int ind = 0;
			for (int i = 0; i < x.size(); i++) {
				if (x[i].value>0.001) {
				std::cout <<  "stress at " <<i<<" : "<< calcVonMises(StressVec[i].Stresses) << std::endl;
				for (int j = 0; j < 12; j++) {
					rhoFaces(ind * 12 + j) = calcVonMises(StressVec[i].Stresses);
					}
				ind++;
				}
			}

			igl::jet(rhoFaces, true, Color);
			viewer.data().set_mesh(V, F);
			viewer.data().set_colors(Color);

			Eigen::VectorXd displacement(24);

		}


		if (optimizerFlag == 3) { //test fatigue optimizer

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

				std::vector<std::vector<double>> sigma_a;

				ClacStressTimeHistory(targets, x, 0.001, 0.001, power, forceVertexIdSet, dx, E, nodes, constraintVertexIdSet, u_dot_0, u0, 0.0001, steps, sigma_m, sigma_a, disp_test);

			//	std::cout << "calculated stress histories" << std::endl;
				if (debug_flag) {
					for (int i = 0; i < targets.size(); i++) {

						for (int j = 0; j < sigma_m[i].size(); j++) {

							std::cout << "i: " << i << " j: " << j << " sigma_m: " << sigma_m[i][j] << " sigma_a: " << sigma_a[i][j] << std::endl;

						}
					}
				}

				double bf = -0.075;
				double sigma_f = 650000000;
				double ni = 1000000;
				double Sut = 380000000;

				for (int i = 0; i < targets.size(); i++) {

					std::vector<double> N(sigma_m[i].size());
					//calculate damage for target
					double D = 0;
					for (int j = 0; j < N.size(); j++) {
						if (((sigma_f*(Sut - std::max(sigma_m[i][j], 0.0)))) < 0) {
							D = 1;
							break;
						}
						D += ni / std::pow(0.5*(Sut*sigma_a[i][j] / (sigma_f*(Sut - std::max(sigma_m[i][j], 0.0)))), 1.0 / bf);
						if (D >= 1) {
							D = 1;
							break;
						}
						//	std::cout << D << std::endl;
					}

					x[targets[i]].rho = D;
				}
				for (int i = 0; i < x.size(); i++) {

			//		std::cout << "element: " << i << ", Damage: " << x[i].rho << std::endl;
					for (int j = 0; j < 12; j++) {
						rhoFaces(i * 12 + j) = x[i].rho;

					}
				}

				igl::jet(rhoFaces, true, Color);
	//			viewer.data().set_mesh(V, F);
				viewer.data().set_colors(Color);

				animation_flag = true;

			}

			std::ofstream writer("matrix.txt");

			writer << disp_test;

			if (animation_flag) {

				// do animation
				auto tempNodes = nodes;
				double factor = 100000000;
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

		}
			
		viewer.core.is_animating = false;
		return false;

		};

	viewer.launch();
	
			return 0;
		}

