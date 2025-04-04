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

	// һֱ��ȡֱ���ҵ����� "POINTS" ����һ��
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

	// ��ȡ���е�����꣬ÿ������3��������
	v_3d_coord.resize(numPoints);
	for (int i = 0; i < numPoints; i++)
	{
		double x, y, z;
		fin >> x >> y >> z;
		v_3d_coord[i] = vec3d(x, y, z);
	}
	// ��ȡ��������Ļ���
	std::getline(fin, line);


	// 3. ��ȡ CELLS ���֣����磺CELLS 4 16��
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

	// 4. ��ȡ OFFSETS ���֣����磺OFFSETS vtktypeint64��
	while (std::getline(fin, line)) {
		if (line.find("OFFSETS") != std::string::npos)
			break;
	}
	if (line.empty() || line.find("OFFSETS") == std::string::npos) {
		std::cerr << "Error: 'OFFSETS' section not found in file: " << volumefilename << std::endl;
		return;
	}

	// ���� OFFSETS �ı�ǩ�У���ʼ�� offset ����
	std::vector<int> offsets(numOffsets);
	for (int i = 0; i < numOffsets; i++) {
		fin >> offsets[i];
	}
	// �����ĩβ�Ļ��з�
	std::getline(fin, line);

	// 5. ��ȡ CONNECTIVITY ����
	//    ���磺CONNECTIVITY vtktypeint64
	std::getline(fin, line);
	while (line.empty()) { std::getline(fin, line); }
	std::istringstream issConnTag(line);
	std::string connTag;
	issConnTag >> connTag;  // Ӧ��Ϊ "CONNECTIVITY"
	// ��ȡ connectivity ���飬����Ϊ offsets.back()
	int totalConnCount = offsets.back();
	std::vector<int> connectivity(totalConnCount);
	for (int i = 0; i < totalConnCount; i++)
	{
		fin >> connectivity[i];
	}
	std::getline(fin, line);

	// 6. ��ȡ CELL_TYPES ����
	//    ���磺CELL_TYPES 3
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

	// 7. ���� offsets ����� cellTypes �Ե�Ԫ���л���
	//     offsets ���� numOffsets �����֣���˵�Ԫ��Ϊ numOffsets - 1
	int numCells = numOffsets - 1;
	for (int c = 0; c < numCells; c++)
	{
		int start = offsets[c];
		int end = offsets[c + 1];
		int nVertices = end - start;
		// ���õ�Ԫ�����ж����Ŷ���
		std::vector<int> cellConnectivity(connectivity.begin() + start, connectivity.begin() + end);
		int cellType = cellTypes[c];
		if (cellType == 10)  // �����壨tetrahedron��ͨ����4������
		{
			if (nVertices != 4)
				std::cerr << "Warning: Tetrahedron cell does not have 4 vertices!" << std::endl;
			LinearTetList.push_back(cellConnectivity);
		}
		else if (cellType == 13) // ����������������ͨ����6������
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


	// 8. ��ȡ Pressure ����

	// һֱ��ȡֱ���ҵ����� "Pressure" ����һ��
	while (std::getline(fin, line)) {
		if (!line.empty() && line.find("Pressure") != std::string::npos)
			break;
	}

	if (line.empty() || line.find("Pressure") == std::string::npos) {
		std::cerr << "Warning: 'Pressure' section not found in file: " << volumefilename << std::endl;
		// �����������󣬲��ж�����
	}
	else {
		std::istringstream issPressure(line);
		std::string pressureLabel;
		int unknownFlag = 0;  // Pressure ������� 1����־������������
		int pressureCount = 0;
		std::string dataType;
		issPressure >> pressureLabel >> unknownFlag >> pressureCount >> dataType;

		if (pressureCount != numPoints) {
			std::cerr << "Warning: Pressure count does not match number of points ("
				<< pressureCount << " vs " << numPoints << ")" << std::endl;
		}

		// ��ȡ pressureCount ��ѹ��ֵ
		pressureList.resize(pressureCount);
		for (int i = 0; i < pressureCount; i++) {
			fin >> pressureList[i];
		}

		// ��ȡ�������ĩβ����
		std::getline(fin, line);
	}




}