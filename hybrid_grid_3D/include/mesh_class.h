#pragma once
#include"include/const.h"
void control();
class Coordinate//����
{
public:
	double x;
	double y;
	double z;
	template <class T>
	double distance(const T& c);
};
class Mesh_S:public Coordinate//�ṹ����
{
public:
	int id;//the address of point
	int connectId;//connected to same loacation of point
	double rho;//�ܶ�
	Coordinate u;//�ٶ�����
	double p;//ѹǿ
	int section = 1;//����
	int sec_num;//���������Ĳ���1��0
	int step = 0;//�����õ��ںβ�����µ�
	int neiborsec = -1;//���ڷ����������ֵΪ�������ʾ�õ㲻�Ǽ����㣬��֮Ϊ���ڷ������õļ�����
	int neiborsec_ad = -1;//���ڷ����ĵ�ַ
	vector <Mesh_S*>neibor;//��¼�ø������ڸ��λ����Ϣ
	vector <int>moveConnct;//�˶������㣬�������˶����moveConnect���й�
	string type;//������ͣ���Ϊ�������ұ߽��Լ��ڲ���U��D,L,R,IN,������SHOCK���Ӵ���ϵ�DISCON��contact discontinuity��,�����ཻ��CENTER
	template <class T>
	Mesh_S& operator =(const T& U);//���ں�����
};
class Mesh_U :public Mesh_S//�ǽṹ����̳нṹ����
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
	vector <double> J;//�ſɱ�����ʽ
	Mesh_U& operator =(const Mesh_S& S);
};
