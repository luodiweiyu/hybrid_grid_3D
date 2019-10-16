#include"include/mesh_class.h"
template <class T>
double Coordinate::distance(const T& c)
{
	double dx = x - c.x;
	double dy = y - c.y;
	double dz = z - c.z;
	return sqrt(dx * dx + dy * dy + dz * dz);
}
template <class T>
Mesh& Mesh::operator=(const T& U)
{
	x = U.x;
	y = U.y;
	z = U.z;
	id = U.id;
	connectId = U.connectId;
	rho = U.rho;
	u = U.u;
	p = U.p;
	section = U.section;
	sec_num = U.sec_num;
	step = U.step;
	neighborsec = U.neighborsec;
	neighborsec_ad = U.neighborsec_ad;
	neighbor = U.neighbor;
	moveConnct = U.moveConnct;
	type = U.type;
	return *this;
}
Mesh& Mesh::operator =(const Coordinate& S)
{
	x = S.x;
	y = S.y;
	z = S.z;
	return *this;
}

