#pragma once
double get_Ma(double u, double v, double rho, double p);//�������
double get_Ma2(double Ma1, double beta);//б���������������Ǯ��𢡶��������ѧ��p241,7-131
double get_delta(double Ma1, double beta);//���������
double get_ufromMa2(double Ma2, double rho2, double p2, double delta);
double get_vfromMa2(double Ma2, double rho2, double p2, double delta);
double get_p2(double Ma1, double beta, double p1);
double get_rho2(double Ma1, double beta, double rho1);
