#include<string>
#include<fstream>
#include<iostream>
#include<stdlib.h>
using namespace std;
double cfl, t_end, gama;
string flowType, methodType;
double meshType;
double xl, xr, yl, yr, zl, zr;
int xnum, ynum, znum,pnum;
double dx, dy, dz;
void control()
{
	cout << "�����ļ���ʼ��ȡ:" << endl;
	ifstream fin("IO_files/control.txt");
	string temp;
	fin >> temp, fin >> cfl;
	fin >> temp, fin >> t_end;
	fin >> temp, fin >> gama;
	fin >> temp, fin >> flowType;
	fin >> temp, fin >> meshType;
	fin >> temp, fin >> methodType;
	fin >> temp, fin >> xl;
	fin >> temp, fin >> xr;
	fin >> temp, fin >> yl;
	fin >> temp, fin >> yr;
	fin >> temp, fin >> zl;
	fin >> temp, fin >> zr;
	fin >> temp, fin >> xnum;
	fin >> temp, fin >> ynum;
	fin >> temp, fin >> znum;
	if (cfl <= 0 || cfl > 1)
	{
		cout << "CFL = " << cfl << "���󣡣�" << endl;
		cout << "������ֹ��" << endl;
		exit(0);
	}
	if (t_end <= 0)
	{
		cout << "t_end = " << t_end << "���󣡣�" << endl;
		cout << "������ֹ��" << endl;
		exit(0);
	}
	if (xl > xr)
	{
		cout << "xl = " << xl << " > xr = " << xr << endl;
		cout << "�����ϳ����߼���" << endl;
		cout << "������ֹ��" << endl;
		exit(0);
	}
	if (yl > yr)
	{
		cout << "yl = " << yl << " > yr = " << yr << endl;
		cout << "�����ϳ����߼���" << endl;
		cout << "������ֹ��" << endl;
		exit(0);
	}
	if (zl > zr)
	{
		cout << "zl = " << zl << " > zr = " << zr << endl;
		cout << "�����ϳ����߼���" << endl;
		cout << "������ֹ��" << endl;
		exit(0);
	}
	if (xnum < 0)
	{
		cout << "xnum = " << xnum << "���󣡣�" << endl;
		cout << "������ֹ��" << endl;
		exit(0);
	}
	if (ynum < 0)
	{
		cout << "ynum = " << ynum << "���󣡣�" << endl;
		cout << "������ֹ��" << endl;
		exit(0);
	}
	if (znum < 0)
	{
		cout << "znum = " << znum << "���󣡣�" << endl;
		cout << "������ֹ��" << endl;
		exit(0);
	}
	cout << "�����ļ���ȡ������" << endl;
	pnum = xnum * ynum * znum;
	dx = (xr - xl) / (xnum - 1);
	dy = (yr - yl) / (ynum - 1);
	dz = (zr - zl) / (znum - 1);
}