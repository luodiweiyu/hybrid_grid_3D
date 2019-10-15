#include"include/const.h"
#include"include/functions.h"
#include"include/shockwave.h"
#include<iostream>
#include<ctime>
#include<stdlib.h>
#include<vector>
#include"include/Prandtl-Meyer.h"
#include"include/init.h"
using std::vector;
extern string flowType;
void init_mesh()
{
	extern vector <Mesh_S> ms;
	extern vector<vector<int>> ad;
	extern string flowType;
	extern double meshType;
	extern int xnum, ynum, znum, pnum;
	extern double dx, dy,dz;
	extern double xl, yl, zl;
	Mesh_S t;
	vector<int> a;
	int i, j, k;
	int size;
	for (i = 0; i < znum; i++)
	{
		for (j = 0; j < ynum; j++)
		{
			for (k = 0; k < xnum; k++)
			{
				ms.push_back(t);
				size = ms.size() - 1;
				if (k == 0)
				{
					ms[size].x = xl;
					if (j == 0)
					{
						ms[size].y = yl;
						if (i == 0)
							ms[size].z = zl;
						else
							ms[size].z = ms[size - xnum * ynum].z + dz;
					}
					else
					{
						ms[size].y = ms[size - xnum].y + dy;
						if (i == 0)
							ms[size].z = zl;
						else
							ms[size].z = ms[size - xnum * ynum].z + dz;
					}
				}
				else
				{
					ms[size].x = ms[size - 1].x + dx;
					if (j == 0)
					{
						ms[size].y = yl;
						if (i == 0)
							ms[size].z = zl;
						else
							ms[size].z = ms[size - xnum * ynum].z + dz;
					}
					else
					{
						ms[size].y = ms[size - xnum].y + dy;
						if (i == 0)
							ms[size].z = zl;
						else
							ms[size].z = ms[size - xnum * ynum].z + dz;
					}

				}
			}
		}
	}
}

void getType()
{
	extern string flowType;
	extern vector <Mesh_S> ms;
	extern int pnum;
	extern double xl, xr, yl, yr, zl, zr;
	int i;
	float DELTA = 1e-10;
	for (i = 0; i < pnum; i++)
	{
		if (abs(ms[i].x - xl) < DELTA)
			ms[i].type = "X-";
		else if (abs(ms[i].x - xr) < DELTA)
			ms[i].type = "X+";
		else if (abs(ms[i].y - yl) < DELTA)
			ms[i].type = "Y-";
		else if (abs(ms[i].y - yr) < DELTA)
			ms[i].type = "Y+";
		else if (abs(ms[i].z - zl) < DELTA)
			ms[i].type = "Z-";
		else if (abs(ms[i].z - zr) < DELTA)
			ms[i].type = "Z+";
		else
			ms[i].type = "IN";
	}
}

void initFlow()
{
	init_flow_cylinder();//Ô²ÖùÈÆÁ÷
}
void init_flow_cylinder()//Ô²ÖùÈÆÁ÷
{
	extern vector <Mesh_S> ms;
	double Ma1;
	int i;
	for (i = 0; i < ms.size(); i++)
	{
		ms[i].rho = 10;
		ms[i].u.x = 0;
		ms[i].u.y = 0;
		ms[i].u.z = 0;
		ms[i].p = 10;
	}
}

