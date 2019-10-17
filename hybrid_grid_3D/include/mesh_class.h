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
class Mesh :public Coordinate//�ṹ����
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
	int neighborsec = -1;//���ڷ����������ֵΪ�������ʾ�õ㲻�Ǽ����㣬��֮Ϊ���ڷ������õļ�����
	int neighborsec_ad = -1;//���ڷ����ĵ�ַ
	vector <Mesh*>neighbor;//��¼�ø������ڸ��λ����Ϣ
	vector <int>moveConnct;//�˶������㣬�������˶����moveConnect���й�
	string type;//������ͣ���Ϊ�������ұ߽��Լ��ڲ���U��D,L,R,IN,������SHOCK���Ӵ���ϵ�DISCON��contact discontinuity��,�����ཻ��CENTER
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
	double J;//�ſɱ�����ʽ
	Mesh& operator =(const Coordinate& S);
	//template <class T>
	//Mesh& operator =(const T& U);//���ں�����

};
