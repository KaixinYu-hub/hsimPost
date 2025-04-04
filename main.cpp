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


#pragma optimize("",off)
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

	//数据结构转换
	std::vector<Tet> MyTetList;
	for (int i = 0; i < LinearTetList.size(); i++)
	{
		Tet temp;
		temp.topo = LinearTetList[i];

		for (auto &index : temp.topo)
		{
			temp.pressure.push_back(pressureList[index]);
		}

		MyTetList.push_back(temp);
	}

	std::vector<Wed> MyWedList;
	for (int i = 0; i < LinearWedList.size(); i++)
	{
		Wed temp;
		temp.topo = LinearWedList[i];

		for (auto& index : temp.topo)
		{
			temp.pressure.push_back(pressureList[index]);
		}

		MyWedList.push_back(temp);
	}

	//去重
	std::map<int, int> map_assemble;
	PointAssemble(1e-6, 0, v_3d_coord, map_assemble);


	for (auto& Tet : MyTetList)
	{
		for (auto& id : Tet.topo)
		{
			id = map_assemble[id];
		}
	}

	for (auto& Wed : MyWedList)
	{
		for (auto& id : Wed.topo)
		{
			id = map_assemble[id];
		}
	}

	std::map<int, double> MyPressureList;

	for (auto& Tet : MyTetList)
	{
		for (int i = 0; i < Tet.topo.size(); i++)
		{
			MyPressureList[Tet.topo[i]] = Tet.pressure[i];
		}
	}
	for (auto& Wed : MyWedList)
	{
		for (int i = 0; i < Wed.topo.size(); i++)
		{
			MyPressureList[Wed.topo[i]] = Wed.pressure[i];
		}
	}

	WriteMeshTopologyVTK(output_file, MyTetList, MyWedList, v_3d_coord, MyPressureList);

	return 0;

}
#pragma optimize("",on)
