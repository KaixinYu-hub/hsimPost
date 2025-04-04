/*
 * @Author: liutaoran
 * @Date: 2022-06-12 19:48:27
 * @Last Modified by:   liutaoran
 * @Last Modified time: 2022-06-12 19:48:27
 */

#include "kdtree_m.h"
KdTreeM::KdTreeM(int dimension)
{
	kd = kd_create(dimension);
}

KdTreeM::~KdTreeM()
{
	kd_free(kd);
}

bool KdTreeM::Insert3DNode(vec3d& point, int id)
{
	int* id_ptr = new int(1);
	*id_ptr = id;
	kd_insert3(this->kd, point.x(), point.y(), point.z(), id_ptr);
	return true;
}

std::vector<int> KdTreeM::Query3DNodeByDistance(vec3d& point, double dis)
{
	res = kd_nearest_range3(kd, point.x(), point.y(), point.z(), dis);
	std::vector<int> vec;
	for (int i = 0; i < kd_res_size(res); i++)
	{
		int* id = static_cast<int*>(kd_res_item_data(res));
		vec.push_back(id[0]);
		if (kd_res_next(res) == 0) break;
	}
	if (res != nullptr) kd_res_free(res);
	return vec;
}


bool KdTreeM::Insert2DNode(vec2d& point, int id)
{
	int* id_ptr = new int[1];
	*id_ptr = id;
	double* point_ptr = new double[2];
	point_ptr[0] = point.x();
	point_ptr[1] = point.y();
	kd_insert(this->kd, point_ptr, id_ptr);
	return true;
}

std::vector<int> KdTreeM::Query2DNodeByDistance(vec2d& point, double dis)
{

	double* point_ptr = new double[2];
	point_ptr[0] = point.x();
	point_ptr[1] = point.y();
	res = kd_nearest_range(kd, point_ptr, dis);
	int number = kd_res_size(res);
	std::vector<int> vec;
	for (int i = 0; i < number; i++)
	{
		int* id = static_cast<int*>(kd_res_item(res, 0));
		vec.push_back(id[0]);
		if (kd_res_next(res) == 0) break;
	}
	if (res != nullptr) kd_res_free(res);
	return vec;
}
