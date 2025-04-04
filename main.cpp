/*
 * @Descripttion: SurfGenDemo
 * @version: 1.0
 * @Author: YuKaixin
 * @Date: 2023-12-06 09:43:37
 * @LastEditors: YuKaixin
 * @LastEditTime: 2023-12-06 09:43:38
 */

#include <iostream>
#include <queue>
#include <common_type.h>
#include <fstream>
#include "CLI11.hpp"
#include <string>


int main(int argc, char** argv)
{
	std::string linear_wed_file;
	std::string output_file;
	std::string savepath;
	double scale_factor = 20;

	CLI::App app{ "High Order Mesh Generation Preparation" };
	//*************************************************************************
	// Input file name
	//*************************************************************************
	app.add_option("-w", savepath,
		"Specify the work path(for saving the output file)")->required();
	app.add_option("-b", linear_wed_file,
		"Specify the linear boundary layer mesh flie path(*.vtk)")->required();
	app.add_option("-o", output_file,
		"Specify the output mesh flie path(*.vtk)")->required();
	app.add_option("--scale", scale_factor,
		"Specify the scale factor for box");

	double size = 1.0;
	app.add_option("--size", size,
		"Specify the size factor for tet mesh");
	
	try
	{
		app.parse(argc, argv);

	}
	catch (const CLI::ParseError& e)
	{
		return app.exit(e);
	}

	std::vector<std::vector<int>> LinearWedList;
	std::vector<vec3d>v_3d_coord;

	std::vector<std::vector<int>> BoxTri;

	ReadLinearWedgeFile(linear_wed_file, v_3d_coord, LinearWedList);

	RemoveNoReferPoint(v_3d_coord, LinearWedList);

	AddBoundingBox(v_3d_coord, LinearWedList, BoxTri);

	ScaleBoundingBox(scale_factor, v_3d_coord, LinearWedList, BoxTri);

	dt::Mesh Tetmesh;
	Tetrahedralize(size, v_3d_coord, LinearWedList, BoxTri, Tetmesh);

	OutPutMesh(output_file, v_3d_coord, LinearWedList, BoxTri, Tetmesh);

	return 0;
}

