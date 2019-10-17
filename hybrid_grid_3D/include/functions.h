#pragma once
#include"const.h"
#include"mesh_class.h"
//double distance(mesh& a, mesh& b);
//double distance(double x1, double y1, double x2, double y2);
void get_dt();
//Flux HLLC_��(mesh& CL, mesh& CR, mesh& C, int method);
//Flux HLLC_��2(mesh& CL, mesh& CR, mesh& C, int method);//�������Ҳ��㣬����任�ο���
//
//Flux HLLC_��(mesh& CD, mesh& CU, mesh& C, int method);
//Flux HLLC_��2(mesh& CD, mesh& CU, mesh& C, int method);
void coordinate_trans();

Flux VanLeer(Mesh& N, string direction, string side);//������任
Flux VanLeer(Mesh& N, Mesh& C, string direction, string side);//������任

void record();
template<class T>
double get_beta(T A, T B);//��������������x��ļн�
void reorder_neighbor();
template <class T>
double area(T A, T B, T C, T D);//�������ĵ㹹���ı������ 
void movemesh();
void findNeiborSec();
void update_Vm();
void clear_Vm();
void update_bound();
double get_theta(double x1, double y1, double x2, double y2);//��ֱ����x��ļн�

double max(double a, double b);
double min(double a, double b);
double absmax(double a, double b);
double absmin(double a, double b);
Coordinate getCrossPoint(double theta, double a, double b, double r);
void polygonPoint(vector <Coordinate>& poly);

Coordinate getCrossPoint(Line L1, Line L2);
template<class T>
Coordinate getCrossPoint(T M, double a, double b, double r);//ĳ���Բ�ĵ�������Բ�Ľ���
template <class T>
Line getLine(T A, T B);
Line getLine(double x1, double y1, double x2, double y2);
//template <class T>
bool judgeFieldInOut(Mesh& A, vector<Coordinate>& poly);
//bool judgeFieldInOut(double x, double y, vector <mesh> &poly);
bool judgeFieldInOut(double x, double y);
template <class T>
Line getLine(double theta, T A);
double compute_res();//����в�
void sortPoint();
//void update_p3(mesh& p);
void update_p4_u(Mesh& p);
void update_p4_s(Mesh& p);
void update_bound();
void polymesh();
template<class T>
int findNearPoint(T A, vector <Coordinate>& poly);//find the closest point of given grid point
int findNearPoint(double x, double y,double z, vector<Coordinate>& poly);//find the closest point of given grid point
