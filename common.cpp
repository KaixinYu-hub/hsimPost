#include "common_type.h"

void ReadLinearHybridMeshFile_vtk5(
	std::string volumefilename,
	std::vector<vec3d>& v_3d_coord,
	std::vector<std::vector<int>>& LinearTetList,
	std::vector<std::vector<int>>& LinearWedList,
	std::vector<double>& pressureList

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
		if (cellType == 10)  // 四面体（tetrahedron）通常有4个顶点
		{
			if (nVertices != 4)
				std::cerr << "Warning: Tetrahedron cell does not have 4 vertices!" << std::endl;
			LinearTetList.push_back(cellConnectivity);
		}
		else if (cellType == 13) // 三棱柱（三角柱）通常有6个顶点
		{
			if (nVertices != 6)
				std::cerr << "Warning: Wedge cell does not have 6 vertices!" << std::endl;

			//only for zgridgen output
			std::swap(cellConnectivity[1], cellConnectivity[2]);
			std::swap(cellConnectivity[4], cellConnectivity[5]);

			LinearWedList.push_back(cellConnectivity);
		}
		else
		{
			std::cerr << "Warning: Unknown cell type " << cellType << " encountered. Skipping." << std::endl;
		}
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

		// 读取完后清理末尾空行
		std::getline(fin, line);
	}




}