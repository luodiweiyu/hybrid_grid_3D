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
class Mesh :public Coordinate//结构网格
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
	int neighborsec = -1;//相邻分区，如果该值为负，则表示该点不是激波点，反之为相邻分区共用的激波点
	int neighborsec_ad = -1;//相邻分区的地址
	vector <Mesh*>neighbor;//记录该格点的相邻格点位置信息
	vector <int>moveConnct;//运动关联点，即本点运动与该moveConnect点有关
	string type;//格点类型，分为上下左右边界以及内部，U，D,L,R,IN,激波点SHOCK，接触间断点DISCON（contact discontinuity）,激波相交点CENTER
	double xix;
	double xiy;
	double xiz;
	double xit;
	double etax;
	double etay;
	double etaz;
	double etat;
	double zetax;
	double zetay;
	double zetaz;
	double zetat;
	double xxi;
	double xeta;
	double xzeta;
	double xtau;
	double yxi;
	double yeta;
	double yzeta;
	double ytau;
	double zxi;
	double zeta;
	double zzeta;
	double ztau;
	double J;//雅可比行列式
	Mesh& operator =(const Coordinate& S);
	//template <class T>
	//Mesh& operator =(const T& U);//等于号重载

};
