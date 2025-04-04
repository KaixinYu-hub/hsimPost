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
#include "boundingbox.h"
#include <climits>
#include <set>
#include <array>
#include <dt_define.h>
#include <dt_API.h>

typedef Eigen::Vector3d vec3d;
typedef Eigen::Vector2d vec2d;
typedef Eigen::Vector3i vec3i;
typedef Eigen::Vector4i vec4i;
typedef Eigen::Vector2i vec2i;
typedef Eigen::MatrixXd mat;

void RemoveNoReferPoint(
	std::vector<vec3d>& v_3d_coord,
	std::vector<std::vector<int>>& LinearWedList)
{
	//checkForClosePoints(BoundingBox::points_array, 1e-5);

	int max_val = INT_MIN;

	// ����������������Ԫ
	for (const auto& wedge : LinearWedList) {

		// ������Ԫ�еĶ�������
		for (int idx : wedge) {
			if (idx > max_val) {
				max_val = idx;
			}
		}
	}

	// ���� v_3d_coord ��ǰ max_val + 1 ��Ԫ��
	if (v_3d_coord.size() > static_cast<size_t>(max_val + 1)) {
		v_3d_coord.resize(max_val + 1);
	}

	return;
}

void AddBoundingBox(
	std::vector<vec3d>& v_3d_coord,
	std::vector<std::vector<int>>& LinearWedList,
	std::vector<std::vector<int>>& BoxTri
)
{
	int origin_size = v_3d_coord.size();
	// Add all points from points_array to v_3d_coord
	for (int i = 0; i < 2288; ++i) {
		vec3d point;
		point << 
			BoundingBox::points_array[i][0],
			BoundingBox::points_array[i][1],
			BoundingBox::points_array[i][2];
		v_3d_coord.push_back(point);
	}

	// Loop through each triangle in cells_array
	for (int i = 0; i < 4572; ++i) {
		std::vector<int> new_triangle;
		// For each vertex in the cell, adjust the index by adding current_size
		new_triangle.push_back(BoundingBox::cells_array[i][0] + origin_size);
		new_triangle.push_back(BoundingBox::cells_array[i][1] + origin_size);
		new_triangle.push_back(BoundingBox::cells_array[i][2] + origin_size);

		BoxTri.push_back(new_triangle); // Add the adjusted triangle to BoxTri
	}

}

void ScaleBoundingBox(
	double scale_factor,
	std::vector<vec3d>& v_3d_coord,
	std::vector<std::vector<int>>& LinearWedList,
	std::vector<std::vector<int>>& BoxTri
)
{
	// ���� LinearWedList �İ�Χ��
	vec3d linear_min = vec3d::Ones() * std::numeric_limits<double>::infinity();
	vec3d linear_max = vec3d::Ones() * -std::numeric_limits<double>::infinity();

	for (const auto& wedge : LinearWedList) {
		for (int idx : wedge) {
			const vec3d& point = v_3d_coord[idx];
			linear_min = linear_min.cwiseMin(point);
			linear_max = linear_max.cwiseMax(point);
		}
	}
	vec3d linear_size = linear_max - linear_min;

	//�������ű���
	double factor = (linear_size.array() * scale_factor).mean();

	// ���ڼ����Χ����صĵ�
	std::set<int> unique_points;

	// ����BoxTri�е�ÿ�������Σ��ҳ�������صĶ���
	for (const auto& triangle : BoxTri) {
		for (int i = 0; i < 3; ++i) {
			unique_points.insert(triangle[i]);
		}
	}

	// �����Χ����ص������
	vec3d center(0.0, 0.0, 0.0);
	for (int index : unique_points) {
		const vec3d& point = v_3d_coord[index];
		center[0] += point[0];
		center[1] += point[1];
		center[2] += point[2];
	}
	center[0] /= unique_points.size();  // ����X�����ƽ��ֵ
	center[1] /= unique_points.size();  // ����Y�����ƽ��ֵ
	center[2] /= unique_points.size();  // ����Z�����ƽ��ֵ

	for (int index : unique_points) {
		vec3d& point = v_3d_coord[index];

		//// ��ƽ�ƶ��㵽ԭ��
		//point[0] -= center[0];
		//point[1] -= center[1];
		//point[2] -= center[2];

		// ���Ÿö��������
		point[0] *= factor;
		point[1] *= factor;
		point[2] *= factor;

		//// �ٽ������ƻذ�Χ������
		//point[0] += center[0];
		//point[1] += center[1];
		//point[2] += center[2];
	}
}

void WriteTriMeshToVTK_fordebug(
	const std::string& filename,
	const std::vector<std::array<double, 3>>& vertices,
	const std::vector<std::array<int, 4>>& smesh
) {
	std::ofstream vtk_file(filename);

	// д��VTKͷ��Ϣ
	vtk_file << "# vtk DataFile Version 2.0\n";
	vtk_file << "Triangular Mesh\n";
	vtk_file << "ASCII\n";
	vtk_file << "DATASET UNSTRUCTURED_GRID\n\n";

	// д�������
	vtk_file << "POINTS " << vertices.size() << " double\n";
	for (const auto& point : vertices) {
		vtk_file << point[0] << " " << point[1] << " " << point[2] << "\n";
	}

	// ���㵥Ԫ������������
	const int total_cells = smesh.size();
	const int total_conn = total_cells * 4;  // ÿ����������Ҫ 3+1=4 ���ֶ�

	// д�뵥Ԫ�����ԣ���ȡÿ��smesh��ǰ3��������Ϊ�����Σ�
	vtk_file << "\nCELLS " << total_cells << " " << total_conn << "\n";
	for (const auto& tri : smesh) {
		vtk_file << "3 " << tri[0] << " " << tri[1] << " " << tri[2] << "\n";
	}

	// д�뵥Ԫ���ͣ�����������Ϊ5��
	vtk_file << "\nCELL_TYPES " << total_cells << "\n";
	for (int i = 0; i < total_cells; ++i) {
		vtk_file << "5\n";
	}

	vtk_file.close();
}

void Tetrahedralize(
	double _size,
	std::vector<vec3d>& v_3d_coord,
	std::vector<std::vector<int>>& LinearWedList,
	std::vector<std::vector<int>>& BoxTri,
	dt::Mesh& Tetmesh
	)
{
	// get smesh
	std::vector<std::array<double, 3>> vertices;
	std::vector<std::array<int, 4>> smesh;

	for (auto& v : v_3d_coord)
	{
		vertices.push_back(std::array<double, 3>{v[0], v[1], v[2]});
	}

	for (auto& wed : LinearWedList)
	{
		smesh.push_back(std::array<int, 4>{wed[3], wed[4], wed[5], 0});
	}

	for (auto& tri : BoxTri)
	{
		smesh.push_back(std::array<int, 4>{tri[0], tri[1], tri[2], 0});
	}

	dt::Args args;
	args.layer.resize(1);
	args.layer[0] = 2;
	args.size = _size;

	Tetmesh.V = vertices;
	Tetmesh.F = smesh;

	//std::string filename = "C:\\Users\\R9000P\\Desktop\\hBLgenPre\\test\\horten\\horten_surface.vtk";
	//API_WriteMesh(filename, Tetmesh,true);

	API_Tetrahedralize(Tetmesh, args);

	//std::string filename = "C:\\Users\\R9000P\\Desktop\\hBLgenPre\\test\\horten\\horten_tet.vtk";
	//API_WriteMesh(filename, Tetmesh);

}


void OutPutMesh(
	const std::string volumefilename,
	std::vector<vec3d>& v_3d_coord,
	std::vector<std::vector<int>>& LinearWedList,
	std::vector<std::vector<int>>& BoxTri,
	dt::Mesh& Tetmesh)
{
	//v_3d_coord.clear();
	//LinearWedList.clear();
	//BoxTri.clear();

	std::ofstream vtk_file(volumefilename);

	// д��VTKͷ��Ϣ
	vtk_file << "# vtk DataFile Version 2.0\n";
	vtk_file << "3D Mesh with Prisms and Triangles\n";
	vtk_file << "ASCII\n";
	vtk_file << "DATASET UNSTRUCTURED_GRID\n\n";

	// д�������
	vtk_file << "POINTS " << v_3d_coord.size() + Tetmesh.V.size() << " double\n";
	for (const auto& point : v_3d_coord) {
		vtk_file << point.x() << " " << point.y() << " " << point.z() << "\n";
	}

	for (const auto& point : Tetmesh.V) {
		vtk_file << point[0] << " " << point[1] << " " << point[2] << "\n";
	}


	// ���㵥Ԫ������������
	const int prism_conn_count = 7 * LinearWedList.size(); // 6+1 per prism
	const int tri_conn_count = 4 * BoxTri.size();          // 3+1 per triangle
	const int tet_conn_count = 5 * Tetmesh.T.size();          // 4+1 per triangle
	const int total_cells = LinearWedList.size() + BoxTri.size() + Tetmesh.T.size();
	const int total_conn = prism_conn_count + tri_conn_count + tet_conn_count;

	// д�뵥Ԫ������
	vtk_file << "CELLS " << total_cells << " " << total_conn << "\n";

	// д�������ε�Ԫ��VTK����5��
	for (const auto& tri : BoxTri) {
		vtk_file << "3";  // ÿ����������3������
		for (int idx : tri) {
			vtk_file << " " << idx;
		}
		vtk_file << "\n";
	}

	// д����������Ԫ��VTK����13��
	for (const auto& prism : LinearWedList) {
		vtk_file << "6";  // ÿ����������6������
		for (int idx : prism) {
			vtk_file << " " << idx;
		}
		vtk_file << "\n";
	}

	for (auto& tet : Tetmesh.T)
	{
		vtk_file << "4";  
		for (int i = 0; i < 4; i++)
		{
			vtk_file << " " << tet[i] + v_3d_coord.size();
		}
		vtk_file << "\n";
	}

	// д�뵥Ԫ����
	vtk_file << "CELL_TYPES " << total_cells << "\n";

	// �����ε�Ԫ���ͣ�5��
	for (size_t i = 0; i < BoxTri.size(); ++i) {
		vtk_file << "5\n";
	}

	// ��������Ԫ���ͣ�13��
	for (size_t i = 0; i < LinearWedList.size(); ++i) {
		vtk_file << "13\n";
	}

	// �����嵥Ԫ���ͣ�10��
	for (size_t i = 0; i < Tetmesh.T.size(); ++i) {
		vtk_file << "10\n";
	}
	
	vtk_file.close();
}


void ReadLinearWedgeFile(
	const std::string volumefilename,
	std::vector<vec3d>& v_3d_coord,
	std::vector<std::vector<int>>& LinearWedList)
{
	std::ifstream file(volumefilename);
	if (!file.is_open()) {
		std::cerr << "Error: Unable to open file " << volumefilename << std::endl;
		return;
	}

	std::string token;
	// ���ΰ� token ��ȡ�ļ�����
	while (file >> token) {
		if (token == "POINTS") {
			int numPoints;
			file >> numPoints;  // ��ȡ����
			std::string dataType;
			file >> dataType;   // ��ȡ�������ͣ����� "float"���������в�������

			// ���ζ�ȡ numPoints �� 3D ����
			for (int i = 0; i < numPoints; ++i) {
				double x, y, z;
				file >> x >> y >> z;
				v_3d_coord.push_back(vec3d(x, y, z));
			}
		}
		else if (token == "CELLS") {
			int numCells;
			file >> numCells;  // ��ȡ��Ԫ����
			int totalInts;
			file >> totalInts; // ���� CELLS ���������������ݵ�������������У�飬���˴����ԣ�

			// ���ζ�ȡÿ����Ԫ����������
			for (int i = 0; i < numCells; ++i) {
				int numVertices;
				file >> numVertices;  // ÿ����Ԫ�Ľڵ����
				if (numVertices == 6) {  // Wedge ��ԪӦ�� 6 ���ڵ�
					int v0, v1, v2, v3, v4, v5;
					file >> v0 >> v1 >> v2 >> v3 >> v4 >> v5;

					std::vector<int>temp{ v0,v1,v2,v3,v4,v5 };

					LinearWedList.push_back(temp);
				}
				else {
					// �����Ԫ�ڵ�����Ϊ6����������Щ����
					for (int j = 0; j < numVertices; ++j) {
						int dummy;
						file >> dummy;
					}
					std::cerr << "Warning: Encountered a cell with " << numVertices
						<< " vertices. Skipping this cell." << std::endl;
				}
			}
		}
		else if (token == "CELL_TYPES") {
			// ��ȡ CELL_TYPES ���֡�ͨ����������֤��Ԫ�����Ƿ���ȷ�����磬VTK_WEDGE ��Ӧ13��
			int numCellTypes;
			file >> numCellTypes;
			for (int i = 0; i < numCellTypes; ++i) {
				int cellType;
				file >> cellType;
				// ��ѡ����� cellType �Ƿ�Ϊ 13��Wedge��
				if (cellType != 13) {
					std::cerr << "Warning: Encountered an unexpected cell type: " << cellType << std::endl;
				}
			}
		}
		else {
			// �������� token����ѡ�����
			// std::cout << "Ignoring token: " << token << std::endl;
		}
	}

	file.close();

}



