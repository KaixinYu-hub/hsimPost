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

typedef Eigen::Vector3d vec3d;
typedef Eigen::Vector2d vec2d;
typedef Eigen::Vector3i vec3i;
typedef Eigen::Vector4i vec4i;
typedef Eigen::Vector2i vec2i;
typedef Eigen::MatrixXd mat;


void ReadLinearHybridMeshFile_vtk5(
	std::string volumefilename,
	std::vector<vec3d>& v_3d_coord,
	std::vector<std::vector<int>>& LinearTetList,
	std::vector<std::vector<int>>& LinearWedList,
	std::vector<double>& pressureList
);




