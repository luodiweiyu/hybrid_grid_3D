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
	extern int xnum, ynum, znum;
	extern vector <Mesh> ap;
	extern string flowType;
	extern double meshType;
	extern int pnum;
	extern double dx, dy,dz;
	extern double xl, yl, zl;
	Mesh t;
	vector<int> a;
	int i, j, k;
	int size;
	for (i = 0; i < znum; i++)
	{
		for (j = 0; j < ynum; j++)
		{
			for (k = 0; k < xnum; k++)
			{
				ap.push_back(t);
				size = ap.size() - 1;
				if (k == 0)
				{
					ap[size].x = xl;
					if (j == 0)
					{
						ap[size].y = yl;
						if (i == 0)
							ap[size].z = zl;
						else
							ap[size].z = ap[size - xnum * ynum].z + dz;
					}
					else
					{
						ap[size].y = ap[size - xnum].y + dy;
						if (i == 0)
							ap[size].z = zl;
						else
							ap[size].z = ap[size - xnum * ynum].z + dz;
					}
				}
				else
				{
					ap[size].x = ap[size - 1].x + dx;
					if (j == 0)
					{
						ap[size].y = yl;
						if (i == 0)
							ap[size].z = zl;
						else
							ap[size].z = ap[size - xnum * ynum].z + dz;
					}
					else
					{
						ap[size].y = ap[size - xnum].y + dy;
						if (i == 0)
							ap[size].z = zl;
						else
							ap[size].z = ap[size - xnum * ynum].z + dz;
					}

				}
			}
		}
	}
}

void getType()
{
	extern string flowType;
	extern vector <Mesh> ap;
	extern int pnum;
	extern double xl, xr, yl, yr, zl, zr;
	int i;
	float DELTA = 1e-10;
	for (i = 0; i < pnum; i++)
	{
		if (abs(ap[i].x - xl) < DELTA)
			ap[i].type = "X-";
		else if (abs(ap[i].x - xr) < DELTA)
			ap[i].type = "X+";
		else if (abs(ap[i].y - yl) < DELTA)
			ap[i].type = "Y-";
		else if (abs(ap[i].y - yr) < DELTA)
			ap[i].type = "Y+";
		else if (abs(ap[i].z - zl) < DELTA)
			ap[i].type = "Z-";
		else if (abs(ap[i].z - zr) < DELTA)
			ap[i].type = "Z+";
		else
			ap[i].type = "IN";
	}
}

void initFlow()
{
	init_flow_cylinder();//Ô²ÖùÈÆÁ÷
}
void init_flow_cylinder()//Ô²ÖùÈÆÁ÷
{
	extern vector <Mesh> ap;
	double Ma1;
	int i;
	for (i = 0; i < ap.size(); i++)
	{
		ap[i].rho = 10;
		ap[i].u.x = 0;
		ap[i].u.y = 0;
		ap[i].u.z = 0;
		ap[i].p = 10;
	}
}

