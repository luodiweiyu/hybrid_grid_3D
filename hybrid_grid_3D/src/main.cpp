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
vector <Mesh_S> ms;//所有结构网格节点
vector <Mesh_S*> ps;//结构网格节点
vector <Mesh_U> mu;//非结构网格节点
vector <Mesh_U*> bound;
vector<Mesh_S> msr;//记录上一时刻的物理量；
vector<Mesh_U> mur;//记录上一时刻的物理量；
vector <polygon_mesh> PM;//多边形网格，用于输出，与计算无关
vector<Coordinate> poly;//物体表面，变成很多散点，用于判断结构网格节点处于网格内外
double dt;
double t_sim = 0;
int step = 0;
double res;
extern double t_end;
int main()
{
	control();
	init_mesh();
	polygonPoint(poly);
	ofstream fout("poly.dat");
	fout << "variables = x, y, z" << endl;
	for (int i=0; i < poly.size(); i++)
		fout << poly[i].x << "  " << poly[i].y << "  " << poly[i].z << endl;
	system("PAUSE");
	/*
	getType();
	partition_Point();
	sortPoint();
	polymesh();
	out_M("mesh/step = " + to_string(step));
	//out_neighbor();
	//coordinate_trans();
	initFlow();
	int i;
	while (t_sim < t_end)
	{
		record();
		get_dt();
		//for (i = 0; i < ps.size(); i++)
		//	update_p4_s(*ps[i]);
		//for (i = 0; i < pu.size(); i++)
		//{
		//	if (pu[i]->neibor.size() == 3)
		//		update_p3(*pu[i]);
		//	else
		//		update_p4_u(*pu[i]);
		//	AP[pu[i]->connectId] = *pu[i];//replace
		//}
		update_bound();
		if (++step % 100 == 0)
		{
			if (abs(res - compute_res()) < 1e-20)
				break;
			else
				res = compute_res();
			cout << "step = " << step << "  t_sim = " << t_sim << "  dt = " << dt << "  res = " << res << endl;
			out_M("mesh/step = " + to_string(step));
		}
		out_res();
		t_sim = t_sim + dt;
	}
	out_M("mesh/step = " + to_string(step));
	system("PAUSE");
	*/
}