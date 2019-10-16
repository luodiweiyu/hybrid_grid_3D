#pragma once
#include<cmath>
#include<string>
#include<vector>
#include"include/const.h"
using std::vector;
using std::string;
const double pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640628;
const double random = 0.495;
const double r = 0.175;
const double a = 1 + r;
const double b = 1.5;
const double c = b;
extern double dx, dy, dz;
extern string flowType;
const int method[12][4] = {
{1,2,0,0},//方法3
{0,1,2,2},
{2,0,1,1},

{0,1,1,2},//方法2
{1,2,2,0},
{2,0,0,1},

{0,0,1,2},//方法1
{1,1,2,0},
{2,2,0,1},

{0,1,2,0},//方法4
{2,0,1,2},
{1,2,0,1},
};
const double delta_r = 1.5;
struct Flux
{
	double f1;
	double f2;
	double f3;
	double f4;
	Flux& operator +(const Flux& f)
	{
		f1 += f.f1;
		f2 += f.f2;
		f3 += f.f3;
		f4 += f.f4;
		return *this;
	}
	Flux& operator *(double x)
	{
		f1 = f1 * x;
		f2 = f2 * x;
		f3 = f3 * x;
		f4 = f4 * x;
		return *this;
	}
	Flux& operator /(double x)
	{
		f1 = f1 / x;
		f2 = f2 / x;
		f3 = f3 / x;
		f4 = f4 / x;
		return *this;
	}
};
struct Line//表示直线Ax+Bx+C=0
{
	double A;
	double B;
	double C;
};
struct polygon_mesh
{
	//都只记录点的序号，用于输出多边形网格
	vector <int> node;
	vector <int>face_start;
	vector <int>face_end;
	int section;
};