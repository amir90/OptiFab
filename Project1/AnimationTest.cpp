#include <igl/opengl/glfw/Viewer.h>
#include <igl/jet.h>
int main(int argc, char *argv[])
{
	// Inline mesh of a cube
	const Eigen::MatrixXd V = (Eigen::MatrixXd(8, 3) << -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.0, -1.0, 1.0, 1.0, 1.0).finished();
	const Eigen::MatrixXi F = (Eigen::MatrixXi(12, 3) << 1, 7, 5, 1, 3, 7, 1, 4, 3, 1, 2, 4, 3, 8, 7, 3, 4, 8, 5, 7, 8, 5, 8, 6, 1, 5, 6, 1, 6, 2, 2, 6, 8, 2, 8, 4).finished().array() - 1;
	igl::opengl::glfw::Viewer viewer;
	// Set mesh
	viewer.data().set_mesh(V, F);

	viewer.data().set_face_based(true);
	viewer.core.is_animating = true;
	// Initialize point
	Eigen::MatrixXd P = (Eigen::MatrixXd(1, 3) << 1.42, 0, 0).finished();
	// function will be  called before every draw
	Eigen::MatrixXd C;
	viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer &)->bool
	{
		// Create orbiting animation
		// update point location. no .clear() necessary
		Eigen::VectorXd Vals = Eigen::VectorXd::Random(12);
	//	std::cout << C << " " << V.size() << " " << C.size() << std::endl;
		igl::jet(Vals,true, C);
		viewer.data().set_colors(C);
		//viewer.data().set_points(P, Eigen::RowVector3d(1, 1, 1));
		return false;
	};
	viewer.launch();
}