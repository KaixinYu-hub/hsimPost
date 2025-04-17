#include "common_type.h"

void ReadLinearHybridMeshFile_vtk5(
	std::string volumefilename,
	std::vector<vec3d>& v_3d_coord,
	std::vector<std::vector<int>>& LinearTetList,
	std::vector<std::vector<int>>& LinearWedList,
	std::vector<double>& pressureList,
	std::vector<vec3d>&velocity

)
{
	std::ifstream fin(volumefilename);
	if (!fin)
	{
		std::cerr << "Cannot open file: " << volumefilename << std::endl;
		return;
	}

	std::string line;

	// 一直读取直到找到包含 "POINTS" 的那一行
	while (std::getline(fin, line))
	{
		if (line.find("POINTS") != std::string::npos)
			break;
	}

	if (line.empty() || line.find("POINTS") == std::string::npos)
	{
		std::cerr << "Error: 'POINTS' section not found in file: " << volumefilename << std::endl;
		return;
	}

	std::istringstream issPoints(line);
	std::string pointsTag;
	int numPoints = 0;
	std::string dataType;
	issPoints >> pointsTag >> numPoints >> dataType;
	if (pointsTag != "POINTS")
	{
		std::cerr << "Error: expected 'POINTS', but got: " << pointsTag << std::endl;
		return;
	}

	// 读取所有点的坐标，每个点有3个浮点数
	v_3d_coord.resize(numPoints);
	for (int i = 0; i < numPoints; i++)
	{
		double x, y, z;
		fin >> x >> y >> z;
		v_3d_coord[i] = vec3d(x, y, z);
	}
	// 读取完点坐标后的换行
	std::getline(fin, line);


	// 3. 读取 CELLS 部分（例如：CELLS 4 16）
	while (std::getline(fin, line)) {
		if (line.find("CELLS") != std::string::npos)
			break;
	}
	if (line.empty() || line.find("CELLS") == std::string::npos) {
		std::cerr << "Error: 'CELLS' section not found in file: " << volumefilename << std::endl;
		return;
	}

	std::istringstream issCells(line);
	std::string cellsTag;
	int numOffsets = 0, totalConnectivity = 0;
	issCells >> cellsTag >> numOffsets >> totalConnectivity;
	if (cellsTag != "CELLS") {
		std::cerr << "Error: expected 'CELLS', but got: " << cellsTag << std::endl;
		return;
	}

	// 4. 读取 OFFSETS 部分（例如：OFFSETS vtktypeint64）
	while (std::getline(fin, line)) {
		if (line.find("OFFSETS") != std::string::npos)
			break;
	}
	if (line.empty() || line.find("OFFSETS") == std::string::npos) {
		std::cerr << "Error: 'OFFSETS' section not found in file: " << volumefilename << std::endl;
		return;
	}

	// 跳过 OFFSETS 的标签行，开始读 offset 数据
	std::vector<int> offsets(numOffsets);
	for (int i = 0; i < numOffsets; i++) {
		fin >> offsets[i];
	}
	// 清理掉末尾的换行符
	std::getline(fin, line);

	// 5. 读取 CONNECTIVITY 部分
	//    例如：CONNECTIVITY vtktypeint64
	std::getline(fin, line);
	while (line.empty()) { std::getline(fin, line); }
	std::istringstream issConnTag(line);
	std::string connTag;
	issConnTag >> connTag;  // 应该为 "CONNECTIVITY"
	// 读取 connectivity 数组，总数为 offsets.back()
	int totalConnCount = offsets.back();
	std::vector<int> connectivity(totalConnCount);
	for (int i = 0; i < totalConnCount; i++)
	{
		fin >> connectivity[i];
	}
	std::getline(fin, line);

	// 6. 读取 CELL_TYPES 部分
	//    例如：CELL_TYPES 3
	std::getline(fin, line);
	while (line.empty()) { std::getline(fin, line); }
	std::istringstream issCT(line);
	std::string cellTypesTag;
	int numCellTypes = 0;
	issCT >> cellTypesTag >> numCellTypes;
	if (cellTypesTag != "CELL_TYPES")
	{
		std::cerr << "Error: expected 'CELL_TYPES', but got: " << cellTypesTag << std::endl;
		return;
	}
	std::vector<int> cellTypes(numCellTypes);
	for (int i = 0; i < numCellTypes; i++)
	{
		fin >> cellTypes[i];
	}
	std::getline(fin, line);

	// 7. 按照 offsets 数组和 cellTypes 对单元进行划分
	//     offsets 中有 numOffsets 个数字，因此单元数为 numOffsets - 1
	int numCells = numOffsets - 1;
	for (int c = 0; c < numCells; c++)
	{
		int start = offsets[c];
		int end = offsets[c + 1];
		int nVertices = end - start;
		// 将该单元的所有顶点编号读出
		std::vector<int> cellConnectivity(connectivity.begin() + start, connectivity.begin() + end);
		int cellType = cellTypes[c];
		if (cellType == 10 || cellType == 71)  // 四面体（tetrahedron）通常有4个顶点
		{
			LinearTetList.push_back(cellConnectivity);
		}
		else if (cellType == 13) // 三棱柱（三角柱）通常有6个顶点
		{
			//std::swap(cellConnectivity[1], cellConnectivity[2]);
			//std::swap(cellConnectivity[4], cellConnectivity[5]);
			LinearWedList.push_back(cellConnectivity);
		}
		else if (cellType == 73)
		{
			if (cellConnectivity.size() == 18)//p2
			{
				std::swap(cellConnectivity[1], cellConnectivity[2]);
				std::swap(cellConnectivity[4], cellConnectivity[5]);
				std::swap(cellConnectivity[13], cellConnectivity[14]);
				std::swap(cellConnectivity[9], cellConnectivity[11]);
				std::swap(cellConnectivity[15], cellConnectivity[17]);
				std::swap(cellConnectivity[6], cellConnectivity[8]);
			}
			else if (cellConnectivity.size() == 40)//p3
			{
				std::swap(cellConnectivity[1], cellConnectivity[2]);
				std::swap(cellConnectivity[20], cellConnectivity[22]);
				std::swap(cellConnectivity[21], cellConnectivity[23]);
				std::swap(cellConnectivity[4], cellConnectivity[5]);

				std::swap(cellConnectivity[7], cellConnectivity[10]);
				std::swap(cellConnectivity[27], cellConnectivity[34]);
				std::swap(cellConnectivity[29], cellConnectivity[36]);
				std::swap(cellConnectivity[13], cellConnectivity[16]);

				std::swap(cellConnectivity[6], cellConnectivity[11]);
				std::swap(cellConnectivity[26], cellConnectivity[35]);
				std::swap(cellConnectivity[28], cellConnectivity[37]);
				std::swap(cellConnectivity[12], cellConnectivity[17]);
			}
			LinearWedList.push_back(cellConnectivity);
		}
		else
		{
			std::cerr << "Warning: Unknown cell type " << cellType << " encountered. Skipping." << std::endl;
		}
	}


	// 9. 读取 Velocity 部分（向量数据）
	while (std::getline(fin, line)) {
		if (!line.empty() && line.find("Velocity") != std::string::npos)
			break;
	}

	if (line.empty() || line.find("Velocity") == std::string::npos) {
		std::cerr << "Warning: 'Velocity' section not found in file: " << volumefilename << std::endl;
		// 不是致命错误
	}
	else {
		std::istringstream issVelocity(line);
		std::string velocityLabel;
		int components = 0;
		int velocityCount = 0;
		std::string dataType;
		issVelocity >> velocityLabel >> components >> velocityCount >> dataType;

		if (components != 3) {
			std::cerr << "Error: Velocity should have 3 components per vector!" << std::endl;
			return;
		}

		if (velocityCount != numPoints) {
			std::cerr << "Warning: Velocity count does not match number of points ("
				<< velocityCount << " vs " << numPoints << ")" << std::endl;
		}

		// 每个点读3个分量
		velocity.resize(velocityCount);
		for (int i = 0; i < velocityCount; i++) {
			double vx, vy, vz;
			fin >> vx >> vy >> vz;
			velocity[i] = vec3d(vx, vy, vz);
		}

		std::getline(fin, line); // 读完清空一行
	}

	// 8. 读取 Pressure 部分

	// 一直读取直到找到包含 "Pressure" 的那一行
	while (std::getline(fin, line)) {
		if (!line.empty() && line.find("Pressure") != std::string::npos)
			break;
	}

	if (line.empty() || line.find("Pressure") == std::string::npos) {
		std::cerr << "Warning: 'Pressure' section not found in file: " << volumefilename << std::endl;
		// 不是致命错误，不中断流程
	}
	else {
		std::istringstream issPressure(line);
		std::string pressureLabel;
		int unknownFlag = 0;  // Pressure 后可能是 1（标志），再是数量
		int pressureCount = 0;
		std::string dataType;
		issPressure >> pressureLabel >> unknownFlag >> pressureCount >> dataType;

		if (pressureCount != numPoints) {
			std::cerr << "Warning: Pressure count does not match number of points ("
				<< pressureCount << " vs " << numPoints << ")" << std::endl;
		}

		// 读取 pressureCount 个压力值
		pressureList.resize(pressureCount);
		for (int i = 0; i < pressureCount; i++) {
			fin >> pressureList[i];
		}

	}

}

void PointAssemble(
	double eps,
	int start_id,
	std::vector<vec3d>& pointList_need_assembled,
	std::map<int, int>& map_assemble
)
{
	std::vector<vec3d>pointList_assembled;

	/*
	Remove duplicate points at the boundaries of higher-order elements
	*/
	std::map<int, int>map_duplicate;//Record duplicate point
	std::map<int, int>map_nonduplicate;//Record non-duplicate point


	KdTreeM kdtree(3);
	std::vector<int> result;

	for (int i = 0; i < pointList_need_assembled.size(); i++)
	{
		result = kdtree.Query3DNodeByDistance(pointList_need_assembled[i], eps);

		int index = i + start_id;

		//duplicate points exist
		if (result.size() > 0)
		{
			map_duplicate[index] = result[0];
		}
		else
		{
			kdtree.Insert3DNode(pointList_need_assembled[i], index);

			pointList_assembled.push_back(pointList_need_assembled[i]);

			map_nonduplicate[index] = -1;

		}

	}

	pointList_need_assembled = pointList_assembled;

	//map the non-duplicate point to new topology 
	int offset = 0;
	for (auto& map_pair : map_nonduplicate)
	{
		if (map_pair.second == -1)
		{
			map_pair.second = start_id + offset;
			offset++;
		}
	}

	//map the duplicate point to non-duplicate point
	for (auto& map_pair : map_duplicate)
	{
		map_pair.second = map_nonduplicate[map_pair.second];
	}

	//assemble the above two hashtable 
	for (auto& map_pair : map_duplicate)
	{
		map_assemble.insert(map_pair);
	}
	for (auto& map_pair : map_nonduplicate)
	{
		map_assemble.insert(map_pair);
	}

}


void MeshAssemble(
	std::vector<vec3d>&v_3d_coord,
	std::vector<std::vector<int>>& LinearTetList,
	std::vector<std::vector<int>>& LinearWedList,
	std::vector<double>& pressureList,
	std::vector<vec3d>&velocityList,
	std::vector<Tet>& MyTetList,
	std::vector<Wed>& MyWedList,
	std::map<int, double>& MyPressureList,
	std::map<int, vec3d>& MyVelocityList
	)
{
	//数据结构转换
	for (int i = 0; i < LinearTetList.size(); i++)
	{
		Tet temp;
		temp.topo = LinearTetList[i];

		for (auto& index : temp.topo)
		{
			temp.pressure.push_back(pressureList[index]);
			temp.velocity.push_back(velocityList[index]);
		}

		MyTetList.push_back(temp);
	}

	for (int i = 0; i < LinearWedList.size(); i++)
	{
		Wed temp;
		temp.topo = LinearWedList[i];

		for (auto& index : temp.topo)
		{
			temp.pressure.push_back(pressureList[index]);
			temp.velocity.push_back(velocityList[index]);
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

	std::map<int, double> pressureSum;
	std::map<int, vec3d> velocitySum;
	std::map<int, int> pointTimes;

	for (auto& Tet : MyTetList)
	{
		for (int i = 0; i < Tet.topo.size(); i++)
		{
			int pid = Tet.topo[i];
			pressureSum[pid] += Tet.pressure[i];
			if (velocitySum.find(pid) == velocitySum.end())
				velocitySum[pid] = vec3d::Zero();
			velocitySum[pid] += Tet.velocity[i];
			pointTimes[pid]++;
		}
	}

	for (auto& Wed : MyWedList)
	{
		for (int i = 0; i < Wed.topo.size(); i++)
		{
			int pid = Wed.topo[i];
			pressureSum[pid] += Wed.pressure[i];
			if (velocitySum.find(pid) == velocitySum.end())
				velocitySum[pid] = vec3d::Zero();
			velocitySum[pid] += Wed.velocity[i];
			pointTimes[pid]++;
		}
	}

	MyPressureList.clear();
	MyVelocityList.clear();

	for (const auto& [pid, count] : pointTimes)
	{
		MyPressureList[pid] = pressureSum[pid] / count;
		MyVelocityList[pid] = velocitySum[pid] / count;
	}

}

void WriteMeshTopologyVTK(const std::string& filename,
	const std::vector<Tet>& tets,
	const std::vector<Wed>& wedges,
	const std::vector<vec3d>& v_3d_coord,
	std::map<int, double>& MyPressureList,
	std::map<int, vec3d>& MyVelocityList)
{
	std::ofstream fout(filename);
	if (!fout)
	{
		std::cerr << "Cannot open file for writing: " << filename << std::endl;
		return;
	}

	fout << "# vtk DataFile Version 3.0\n";
	fout << "Hybrid Mesh Output\n";
	fout << "ASCII\n";
	fout << "DATASET UNSTRUCTURED_GRID\n";

	// 1. 统计所有出现过的点编号（去重）
	std::set<int> unique_ids;
	for (const auto& t : tets) unique_ids.insert(t.topo.begin(), t.topo.end());
	for (const auto& w : wedges) unique_ids.insert(w.topo.begin(), w.topo.end());

	// 2. 建立旧ID -> 新ID 映射
	std::map<int, int> global_id_map;
	std::vector<vec3d> used_coords;
	int new_id = 0;
	for (int old_id : unique_ids)
	{
		global_id_map[old_id] = new_id++;
		used_coords.push_back(v_3d_coord[old_id]);
	}

	// 3. 输出 POINTS 坐标
	fout << "POINTS " << used_coords.size() << " float\n";
	for (const auto& pt : used_coords)
	{
		fout << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
	}

	// 4. 输出 CELLS
	int total_cells = static_cast<int>(tets.size() + wedges.size());
	int total_entries = static_cast<int>(tets.size()) * (tets[0].topo.size() + 1)
		+ static_cast<int>(wedges.size()) * (wedges[0].topo.size() + 1);

	fout << "CELLS " << total_cells << " " << total_entries << "\n";

	for (const auto& t : tets)
	{
		fout << t.topo.size() << " ";
		for (int id : t.topo) fout << global_id_map[id] << " ";
		fout << "\n";
	}

	for (const auto& w : wedges)
	{
		fout << w.topo.size() << " ";
		for (int id : w.topo) fout << global_id_map[id] << " ";
		fout << "\n";
	}

	// 5. 输出 CELL_TYPES
	fout << "CELL_TYPES " << total_cells << "\n";
	for (size_t i = 0; i < tets.size(); ++i) fout << "71\n";
	for (size_t i = 0; i < wedges.size(); ++i) fout << "73\n";

	// 4. 输出点数据：压力字段（直接遍历 map）
	fout << "POINT_DATA " << MyPressureList.size() << "\n";
	fout << "SCALARS Pressure float 1\n";
	fout << "LOOKUP_TABLE default\n";
	for (const auto& [id, pressure] : MyPressureList)
	{
		fout << pressure << "\n";
	}

	// 6.2 输出 Velocity（矢量）
	fout << "NORMALS Velocity float\n";
	for (const auto& [id, velocity] : MyVelocityList)
	{
		fout << velocity[0] << " " << velocity[1] << " " << velocity[2] << "\n";
	}

	fout.close();
	std::cout << "Mesh topology written to: " << filename << std::endl;
}
