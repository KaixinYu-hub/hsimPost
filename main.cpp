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
	
	try
	{
		app.parse(argc, argv);

	}
	catch (const CLI::ParseError& e)
	{
		return app.exit(e);
	}

	std::vector<std::vector<int>> LinearTetList;
	std::vector<std::vector<int>> LinearWedList;
	std::vector<vec3d>v_3d_coord;

	std::vector<double> pressureList;
	ReadLinearHybridMeshFile_vtk5(
		linear_wed_file,
		v_3d_coord,
		LinearTetList,
		LinearWedList, 
		pressureList);

	std::vector<Tet> MyTetList;
	std::vector<Wed> MyWedList;
	std::map<int, double>MyPressureList;
	MeshAssemble(
		v_3d_coord,
		LinearTetList, LinearWedList, pressureList,
		MyTetList, MyWedList, MyPressureList);

	WriteMeshTopologyVTK(output_file, MyTetList, MyWedList, v_3d_coord, MyPressureList);

	return 0;

}
