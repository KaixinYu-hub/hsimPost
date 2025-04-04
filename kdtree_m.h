/*
 * @Author: liutaoran 
 * @Date: 2022-06-12 19:48:32 
 * @Last Modified by:   liutaoran 
 * @Last Modified time: 2022-06-12 19:48:32 
 */

#ifndef _KDTREE_M_H_
#define _KDTREE_M_H_

#include "kdtree.h"
#include "Eigen/Dense"
#include <vector>
#include "common_type.h"

namespace HmeshGen
{
	class KdTreeM {

	public:

		kdtree* kd;
		kdres* res;


		KdTreeM(int dimension);
		~KdTreeM();

		bool Insert3DNode(vec3d& point, int id);
		std::vector<int> Query3DNodeByDistance(vec3d& point, double dis);

		bool Insert2DNode(vec2d& point, int id);
		std::vector<int> Query2DNodeByDistance(vec2d& point, double dis);
	};
}
#endif
