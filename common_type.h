/*
 * @Descripttion: SurfGenDemo
 * @version: 1.0
 * @Author: YuKaixin
 * @Date: 2022-07-20 21:11:51
 * @LastEditors: YuKaixin
 * @LastEditTime: 2023-12-06 17:07:35
 */
#pragma once

#include <Eigen/Dense>
#include <functional>
#include <vector>
#include <fstream>
#include <string>
#include <climits>
#include <set>
#include <iostream>
#include <array>
#include <map>


typedef Eigen::Vector3d vec3d;
typedef Eigen::Vector2d vec2d;
typedef Eigen::Vector3i vec3i;
typedef Eigen::Vector4i vec4i;
typedef Eigen::Vector2i vec2i;
typedef Eigen::MatrixXd mat;

#include "kdtree_m.h"

struct Tet
{
public:
	std::vector<int>topo;
	std::vector<double>pressure;
};

struct Wed
{
public:
	std::vector<int>topo;
	std::vector<double>pressure;
};

void ReadLinearHybridMeshFile_vtk5(
	std::string volumefilename,
	std::vector<vec3d>& v_3d_coord,
	std::vector<std::vector<int>>& LinearTetList,
	std::vector<std::vector<int>>& LinearWedList,
	std::vector<double>& pressureList
);

void PointAssemble(
	double eps,
	int start_id,
	std::vector<vec3d>& pointList_need_assembled,
	std::map<int, int>& map_assemble
);

void MeshAssemble(
	std::vector<vec3d>& v_3d_coord,
	std::vector<std::vector<int>>& LinearTetList,
	std::vector<std::vector<int>>& LinearWedList,
	std::vector<double>& pressureList,
	std::vector<Tet>& MyTetList,
	std::vector<Wed>& MyWedList,
	std::map<int, double>& MyPressureList
);

void WriteMeshTopologyVTK(const std::string& filename,
	const std::vector<Tet>& tets,
	const std::vector<Wed>& wedges,
	const std::vector<vec3d>& v_3d_coord,
	std::map<int, double>& MyPressureList);





