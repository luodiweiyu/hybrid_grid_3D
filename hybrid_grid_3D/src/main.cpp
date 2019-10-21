#include<iostream>
#include<fstream>
#include<cmath>
#include<omp.h>
#include"include/const.h"
#include"include/functions.h"
#include"include/init.h"
#include"include/output.h"
#include"include/patition.h"
#include<stdlib.h>
#include<ctime>
#include<iomanip>
#include <vector>
#include<ctime>
#include<algorithm>
#include"include/mesh_class.h"
using namespace std;
vector<Mesh> ap;//all point
vector <Mesh*> ms;//所有结构网格节点
vector <Mesh*> mu;//非结构网格节点
vector <Mesh*> bound;
vector<Mesh> apr;//记录上一时刻的物理量；
vector <polygon_mesh> PM;//多边形网格，用于输出，与计算无关
vector<Coordinate> poly;//物体表面，变成很多散点，用于判断结构网格节点处于网格内外
double dt;
double t_sim = 0;
int step = 0;
double res;
extern double t_end;
extern int xnum, ynum, znum;
int main()
{
	control();
	init_mesh();
	polygonPoint(poly);
	getType();
	partition_Point();
	ofstream fout("poly.dat");
	fout << "variables = x, y, z" << endl;
	for (int i = 0; i < poly.size(); i++)
		fout << poly[i].x << "  " << poly[i].y << "  " << poly[i].z << endl;
	fout.close();
	sortPoint();
	fout.open("ms.dat");
	fout << "variables = x, y, z" << endl;
	for (int i = 0; i < ms.size(); i++)
	{
		if (ms[i]->section != 0)
			fout << ms[i]->x << "  " << ms[i]->y << "  " << ms[i]->z << endl;
	}
	fout.close();
	fout.open("mu.dat");
	fout << "variables = x, y, z" << endl;
	for (int i = 0; i < mu.size(); i++)
		fout << mu[i]->x << "  " << mu[i]->y << "  " << mu[i]->z << endl;
	fout.close();
	fout.open("bound.dat");
	fout << "variables = x, y, z" << endl;
	for (int i = 0; i < bound.size(); i++)
	{
		if (bound[i]->type == "Body")
			fout << bound[i]->x << "  " << bound[i]->y << "  " << bound[i]->z << endl;
	}
	fout.close();
	//polymesh();
	//out_M("mesh/step = " + to_string(step));
	//out_neighbor();
	coordinate_trans();
	initFlow();
	int i;

	while (t_sim < t_end)
	{
		record();
		get_dt();
#pragma  omp parallel for
		for (i = 0; i < ms.size(); i++)
			update_p4_s(*ms[i]);
		for (i = 0; i < mu.size(); i++)
		{
			update_p4_s(*mu[i]);
			//mu[i]->rho = ap[mu[i]->connectId - 1].rho;
			//mu[i]->p = ap[mu[i]->connectId - 1].p;
			//mu[i]->u.x = ap[mu[i]->connectId - 1].u.x;
			//mu[i]->u.y = ap[mu[i]->connectId - 1].u.y;
			//mu[i]->u.z = ap[mu[i]->connectId - 1].u.z;
			//mu[i]->rho = 10;
			//mu[i]->p = 10;
			//mu[i]->u.x = 0;
			//mu[i]->u.y = 0;
			//mu[i]->u.z = 0;
			//double DELTA = 1e-10;
			//Coordinate n, en;
			//Coordinate u = mu[i]->u;
			//n.x = mu[i]->neighbor[0]->x - a;
			//n.y = mu[i]->neighbor[0]->y - b;
			//n.z = mu[i]->neighbor[0]->z -c;
			//en.x = n.x / sqrt(n.x * n.x + n.y * n.y + n.z * n.z);
			//en.y = n.y / sqrt(n.x * n.x + n.y * n.y + n.z * n.z);
			//en.z = n.z / sqrt(n.x * n.x + n.y * n.y + n.z * n.z);
			//if ((abs(u.x) < DELTA && abs(u.y) < DELTA) && abs(u.z) < DELTA/* || (abs(tx) < DELTA || abs(ty) < DELTA)*/)
			//{
			//	mu[i]->rho = mu[i]->neighbor[0]->rho;
			//	mu[i]->u = mu[i]->neighbor[0]->u;
			//	mu[i]->p = mu[i]->neighbor[0]->p;
			//}
			//else
			//{
			//	double n_projection = (u.x * n.x + u.y * n.y + u.z * n.z) / sqrt(n.x * n.x + n.y * n.y + n.z * n.z);
			//	mu[i]->rho = mu[i]->neighbor[0]->rho;
			//	mu[i]->u.x = u.x - en.x * n_projection;
			//	mu[i]->u.y = u.y - en.y * n_projection;
			//	mu[i]->u.z = u.z - en.z * n_projection;
			//	mu[i]->p = mu[i]->neighbor[0]->p;
			//}

			ap[mu[i]->connectId] = *mu[i];
		}
		update_bound();
		if (++step % 100 == 0)
		{
			//if (abs(res - compute_res()) < 1e-20)
			//	break;
			//else
			res = compute_res();
			cout << "step = " << step << "  t_sim = " << t_sim << "  dt = " << dt << "  res = " << res << endl;
			fout.open("mesh/ap step = " + to_string(step) + ".dat");
			fout << "variables = x, y, z,rho" << endl;
			for (int i = 0; i < ap.size(); i++)
			{
				if(ap[i].section==1||ap[i].section==0)
				fout << ap[i].x << "  " << ap[i].y << "  " << ap[i].z << "  " << ap[i].rho<< endl;
			}
			fout.close();
			//fout.open("mesh/ms step = " + to_string(step) + ".dat");
			//fout << "variables = x, y, z,rho" << endl;
			//for (int i = 0; i < ms.size(); i++)
			//{
			//	//if (ms[i]->section != 1|| ms[i]->section != 0)
			//	//	cout << ms[i]->section << endl;
			//		fout << ms[i]->x << "  " << ms[i]->y << "  " << ms[i]->z << "  " << ms[i]->rho << endl;
			//}
			//fout.close();
		}
		//out_res();
		t_sim = t_sim + dt;
	}
	//out_M("mesh/step = " + to_string(step));

}