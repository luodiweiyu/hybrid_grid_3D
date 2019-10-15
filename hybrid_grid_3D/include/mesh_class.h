#pragma once
#include"include/const.h"
void control();
class Coordinate//坐标
{
public:
	double x;
	double y;
	double z;
	template <class T>
	double distance(const T& c);
};
class Mesh_S:public Coordinate//结构网格
{
public:
	int id;//the address of point
	int connectId;//connected to same loacation of point
	double rho;//密度
	Coordinate u;//速度向量
	double p;//压强
	int section = 1;//分区
	int sec_num;//分区后计算的参数1或0
	int step = 0;//表明该点在何步骤更新的
	int neiborsec = -1;//相邻分区，如果该值为负，则表示该点不是激波点，反之为相邻分区共用的激波点
	int neiborsec_ad = -1;//相邻分区的地址
	vector <Mesh_S*>neibor;//记录该格点的相邻格点位置信息
	vector <int>moveConnct;//运动关联点，即本点运动与该moveConnect点有关
	string type;//格点类型，分为上下左右边界以及内部，U，D,L,R,IN,激波点SHOCK，接触间断点DISCON（contact discontinuity）,激波相交点CENTER
	template <class T>
	Mesh_S& operator =(const T& U);//等于号重载
};
class Mesh_U :public Mesh_S//非结构网格继承结构网格
{
public:
	vector <double> xix;
	vector <double> xiy;
	vector <double> xiz;
	vector <double> xit;
	vector <double> etax;
	vector <double> etay;
	vector <double> etaz;
	vector <double> etat;
	vector <double> zetax;
	vector <double> zetay;
	vector <double> zetaz;
	vector <double> zetat;
	vector <double> xxi;
	vector <double> xeta;
	vector <double> xzeta;
	vector <double> xtau;
	vector <double> yxi;
	vector <double> yeta;
	vector <double> yzeta;
	vector <double> ytau;
	vector <double> zxi;
	vector <double> zeta;
	vector <double> zzeta;
	vector <double> ztau;
	vector <double> J;//雅可比行列式
	Mesh_U& operator =(const Mesh_S& S);
};
