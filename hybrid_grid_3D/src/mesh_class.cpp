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
Mesh_S& Mesh_S::operator=(const T& U)
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
	neiborsec = U.neiborsec;
	neiborsec_ad = U.neiborsec_ad;
	neibor = U.neibor;
	moveConnct = U.moveConnct;
	type = U.type;
	return *this;
}
Mesh_U& Mesh_U::operator =(const Mesh_S& S)
{
	x = S.x;
	y = S.y;
	z = S.z;
	id = S.id;
	connectId = S.connectId;
	rho = S.rho;
	u = S.u;
	p = S.p;
	section = S.section;
	sec_num = S.sec_num;
	step = S.step;
	neiborsec = S.neiborsec;
	neiborsec_ad = S.neiborsec_ad;
	neibor = S.neibor;
	moveConnct = S.moveConnct;
	type = S.type;
	return *this;
}

