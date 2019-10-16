#include"include/const.h"
#include"include/functions.h"
#include<omp.h>
#include"include/shockwave.h"
#include<iostream>
#include"include/Prandtl-Meyer.h"
using std::vector;
void get_dt()
{
	extern vector<Mesh>ap;
	double maxxi = 0, maxeta = 0, maxzeta = 0;
	extern double dt, t_sim;
	double t;
	int i, j, k;
	double max1, max2, max3;
	double Sxi, Seta, Szeta, c, uxi, ueta, uzeta;
	extern double t_end, gama, cfl;
	dt = t_end;
	max1 = max2 = max3 = 0;
	for (i = 0; i < ap.size(); i++)
	{
		maxxi = maxeta = 0;
		if (ap[i].type != "IN")
			continue;
		if (ap[i].neighbor.size() == 6)
		{
			Sxi = sqrt(ap[i].xix * ap[i].xix + ap[i].xiy * ap[i].xiy + ap[i].xiz * ap[i].xiz);
			Seta = sqrt(ap[i].etax * ap[i].etax + ap[i].etay * ap[i].etay + ap[i].etaz * ap[i].etaz);
			Szeta = sqrt(ap[i].zetax * ap[i].zetax + ap[i].zetay * ap[i].zetay + ap[i].zetaz * ap[i].zetaz);
			c = sqrt(gama * ap[i].p / ap[i].rho);
			uxi = ap[i].u.x * ap[i].xix + ap[i].u.y * ap[i].xiy + ap[i].u.z * ap[i].xiz;
			ueta = ap[i].u.x * ap[i].etax + ap[i].u.y * ap[i].etay + ap[i].u.z * ap[i].etaz;
			uzeta = ap[i].u.x * ap[i].zetax + ap[i].u.y * ap[i].zetay + ap[i].u.z * ap[i].zetaz;

			maxxi = abs(uxi) + c * Sxi;
			maxeta = abs(ueta) + c * Seta;
			maxzeta = abs(uzeta) + c * Szeta;
			max1 = max(max1, maxxi);
			max2 = max(max2, maxeta);
			max2 = max(max3, maxzeta);
		}
		else
			std::cout << "id = " << i << "不是6个相邻节点！neighborSize = "<<ap[i].neighbor.size() << std::endl;
	}

	t = cfl / (max1 + max2 + max3);
	//t = CFL / (maxxi + maxeta);
	dt = min(dt, t);
	if (t_sim + dt > t_end)
		dt = t_end - t_sim;
}
double compute_res()//计算残差
{
	extern vector <Mesh> ap;
	extern vector <Mesh> apr;
	extern double dt;
	int i, j;
	double res = 0;
	int n = 0;
	for (i = 0; i < ap.size(); i++)
	{
		res = max(res, abs(ap[i].rho - apr[ap[i].id].rho) / apr[ap[i].id].rho);
		n++;
	}
	return res / n;
}
void record()
{
	extern vector <Mesh> ap;
	extern vector <Mesh> apr;
	extern double t_sim;
	int i, j;
	if (t_sim == 0)
		for (i = 0; i < ap.size(); i++)
			apr.push_back(ap[i]);
	else
		for (i = 0; i < ap.size(); i++)
			apr[i] = ap[i];
}

//void update_p3(Mesh& p)
////unstructral grid point,3 neighbor points
//{
//	extern vector <Mesh> mur;
//	extern double dt;
//	int i, j;
//	int n1, n2, n3, n4;
//	double U[4], U1[4], U2[4];
//	extern double gama;
//	Flux Fll, Flr, Fcl, Fcr, Frl, Frr;
//	Flux Gdd, Gdu, Gcd, Gcu, Gud, Guu;
//
//	int id;
//	id = p.id;
//	U[0] = mur[id].rho;
//	U[1] = mur[id].rho * p.u.x;
//	U[2] = mur[id].rho * p.u.y;
//	U[3] = 0.5 * mur[id].rho * (mur[id].u.x * mur[id].u.x + mur[id].u.y * mur[id].u.y) + mur[id].p / (gama - 1);
//
//	for (j = 0; j < 1; j++)
//	{
//		n1 = method[j][0];
//		n2 = method[j][1];
//		n3 = method[j][2];
//		n4 = method[j][3];
//		n1 = p.neighbor[n1]->id;
//		n2 = p.neighbor[n2]->id;
//		n3 = p.neighbor[n3]->id;
//		n4 = p.neighbor[n4]->id;
//		Fll = Flr = Fcl = Fcr = Frl = Frr = { 0 };
//		Gdd = Gdu = Gcd = Gcu = Gud = Guu = { 0 };
//
//		Fll = Fll + VanLeerB(mur[n3], mur[id].xix[j], mur[id].xiy[j], mur[id].xit[j], mur[id].J[j]) * p.J[j];
//		Flr = Flr + VanLeerA(mur[n3], mur[id].xix[j], mur[id].xiy[j], mur[id].xit[j], mur[id].J[j]) * p.J[j];
//		Fcl = Fcl + VanLeerB(mur[id], mur[id].xix[j], mur[id].xiy[j], mur[id].xit[j], mur[id].J[j]) * p.J[j];
//		Fcr = Fcr + VanLeerA(mur[id], mur[id].xix[j], mur[id].xiy[j], mur[id].xit[j], mur[id].J[j]) * p.J[j];
//		Frl = Frl + VanLeerB(mur[n1], mur[id].xix[j], mur[id].xiy[j], mur[id].xit[j], mur[id].J[j]) * p.J[j];
//		Frr = Frr + VanLeerA(mur[n1], mur[id].xix[j], mur[id].xiy[j], mur[id].xit[j], mur[id].J[j]) * p.J[j];
//
//		Gdd = Gdd + VanLeerB(mur[n4], mur[id].etax[j], mur[id].etay[j], mur[id].etat[j], mur[id].J[j]) * p.J[j];
//		Gdu = Gdu + VanLeerA(mur[n4], mur[id].etax[j], mur[id].etay[j], mur[id].etat[j], mur[id].J[j]) * p.J[j];
//		Gcd = Gcd + VanLeerB(mur[id], mur[id].etax[j], mur[id].etay[j], mur[id].etat[j], mur[id].J[j]) * p.J[j];
//		Gcu = Gcu + VanLeerA(mur[id], mur[id].etax[j], mur[id].etay[j], mur[id].etat[j], mur[id].J[j]) * p.J[j];
//		Gud = Gud + VanLeerB(mur[n2], mur[id].etax[j], mur[id].etay[j], mur[id].etat[j], mur[id].J[j]) * p.J[j];
//		Guu = Guu + VanLeerA(mur[n2], mur[id].etax[j], mur[id].etay[j], mur[id].etat[j], mur[id].J[j]) * p.J[j];
//	}
//	Fll = Fll / j, Flr = Flr / j, Fcl = Fcl / j, Fcr = Fcr / j, Frl = Frl / j, Frr = Frr / j;
//	Gdd = Gdd / j, Gdu = Gdu / j, Gcd = Gcd / j, Gcu = Gcu / j, Gud = Gud / j, Guu = Guu / j;
//
//	U[0] = U[0] - dt * (Fcr.f1 - Flr.f1 + Frl.f1 - Fcl.f1 + Gcu.f1 - Gdu.f1 + Gud.f1 - Gcd.f1);
//	U[1] = U[1] - dt * (Fcr.f2 - Flr.f2 + Frl.f2 - Fcl.f2 + Gcu.f2 - Gdu.f2 + Gud.f2 - Gcd.f2);
//	U[2] = U[2] - dt * (Fcr.f3 - Flr.f3 + Frl.f3 - Fcl.f3 + Gcu.f3 - Gdu.f3 + Gud.f3 - Gcd.f3);
//	U[3] = U[3] - dt * (Fcr.f4 - Flr.f4 + Frl.f4 - Fcl.f4 + Gcu.f4 - Gdu.f4 + Gud.f4 - Gcd.f4);
//	p.rho = U[0];
//	p.u.x = U[1] / U[0];
//	p.u.y = U[2] / U[0];
//	p.p = (gama - 1) * (U[3] - 0.5 * p.rho * (p.u.x * p.u.x + p.u.y * p.u.y));
//}
//void update_p4_s(Mesh& p)
////structral grid point,4 neighbor points
////事实上不应该坐标变换
//{
//	extern vector <Mesh> msr;
//	extern double dt;
//	int i, j;
//	int n1, n2, n3, n4;
//	double U[4], U1[4], U2[4];
//	extern double gama;
//	Flux Fll, Flr, Fcl, Fcr, Frl, Frr;
//	Flux Gdd, Gdu, Gcd, Gcu, Gud, Guu;
//
//	int id;
//	id = p.id;
//	U[0] = msr[id].rho;
//	U[1] = msr[id].rho * p.u.x;
//	U[2] = msr[id].rho * p.u.y;
//	U[3] = 0.5 * msr[id].rho * (msr[id].u.x * msr[id].u.x + msr[id].u.y * msr[id].u.y) + msr[id].p / (gama - 1);
//
//	n1 = p.neighbor[0]->id;
//	n2 = p.neighbor[1]->id;
//	n3 = p.neighbor[2]->id;
//	n4 = p.neighbor[3]->id;
//
//
//	Fll = VanLeerB(msr[n3], msr[id].xix[0], msr[id].xiy[0], msr[id].xit[0], msr[id].J[0]);
//	Flr = VanLeerA(msr[n3], msr[id].xix[0], msr[id].xiy[0], msr[id].xit[0], msr[id].J[0]);
//	Fcl = VanLeerB(msr[id], msr[id].xix[0], msr[id].xiy[0], msr[id].xit[0], msr[id].J[0]);
//	Fcr = VanLeerA(msr[id], msr[id].xix[0], msr[id].xiy[0], msr[id].xit[0], msr[id].J[0]);
//	Frl = VanLeerB(msr[n1], msr[id].xix[0], msr[id].xiy[0], msr[id].xit[0], msr[id].J[0]);
//	Frr = VanLeerA(msr[n1], msr[id].xix[0], msr[id].xiy[0], msr[id].xit[0], msr[id].J[0]);
//
//	Gdd = VanLeerB(msr[n4], msr[id].etax[0], msr[id].etay[0], msr[id].etat[0], msr[id].J[0]);
//	Gdu = VanLeerA(msr[n4], msr[id].etax[0], msr[id].etay[0], msr[id].etat[0], msr[id].J[0]);
//	Gcd = VanLeerB(msr[id], msr[id].etax[0], msr[id].etay[0], msr[id].etat[0], msr[id].J[0]);
//	Gcu = VanLeerA(msr[id], msr[id].etax[0], msr[id].etay[0], msr[id].etat[0], msr[id].J[0]);
//	Gud = VanLeerB(msr[n2], msr[id].etax[0], msr[id].etay[0], msr[id].etat[0], msr[id].J[0]);
//	Guu = VanLeerA(msr[n2], msr[id].etax[0], msr[id].etay[0], msr[id].etat[0], msr[id].J[0]);
//
//	U[0] = U[0] - dt * p.J[0] * p.sec_num * (Fcr.f1 - Flr.f1 + Frl.f1 - Fcl.f1 + Gcu.f1 - Gdu.f1 + Gud.f1 - Gcd.f1);
//	U[1] = U[1] - dt * p.J[0] * p.sec_num * (Fcr.f2 - Flr.f2 + Frl.f2 - Fcl.f2 + Gcu.f2 - Gdu.f2 + Gud.f2 - Gcd.f2);
//	U[2] = U[2] - dt * p.J[0] * p.sec_num * (Fcr.f3 - Flr.f3 + Frl.f3 - Fcl.f3 + Gcu.f3 - Gdu.f3 + Gud.f3 - Gcd.f3);
//	U[3] = U[3] - dt * p.J[0] * p.sec_num * (Fcr.f4 - Flr.f4 + Frl.f4 - Fcl.f4 + Gcu.f4 - Gdu.f4 + Gud.f4 - Gcd.f4);
//	p.rho = U[0];
//	p.u = U[1] / U[0];
//	p.v = U[2] / U[0];
//	p.p = (gama - 1) * (U[3] - 0.5 * p.rho * (p.u * p.u + p.v * p.v));
//}
//void update_p4_u(Mesh& p)
////unstructral grid point,4 neighbor points
//{
//	extern vector <Mesh> mur;
//	extern double dt;
//	int i, j;
//	int n1, n2, n3, n4;
//	double U[4], U1[4], U2[4];
//	extern double gama;
//	Flux Fll, Flr, Fcl, Fcr, Frl, Frr;
//	Flux Gdd, Gdu, Gcd, Gcu, Gud, Guu;
//
//	int id;
//	id = p.id;
//	U[0] = mur[id].rho;
//	U[1] = mur[id].rho * p.u.x;
//	U[2] = mur[id].rho * p.u.y;
//	U[3] = 0.5 * mur[id].rho * (mur[id].u.x * mur[id].u.x + mur[id].u.y * mur[id].u.y) + mur[id].p / (gama - 1);
//
//	n1 = p.neighbor[0]->id;
//	n2 = p.neighbor[1]->id;
//	n3 = p.neighbor[2]->id;
//	n4 = p.neighbor[3]->id;
//	Fll = VanLeerB(mur[n3], mur[id].xix[0], mur[id].xiy[0], mur[id].xit[0], mur[id].J[0]);
//	Flr = VanLeerA(mur[n3], mur[id].xix[0], mur[id].xiy[0], mur[id].xit[0], mur[id].J[0]);
//	Fcl = VanLeerB(mur[id], mur[id].xix[0], mur[id].xiy[0], mur[id].xit[0], mur[id].J[0]);
//	Fcr = VanLeerA(mur[id], mur[id].xix[0], mur[id].xiy[0], mur[id].xit[0], mur[id].J[0]);
//	Frl = VanLeerB(mur[n1], mur[id].xix[0], mur[id].xiy[0], mur[id].xit[0], mur[id].J[0]);
//	Frr = VanLeerA(mur[n1], mur[id].xix[0], mur[id].xiy[0], mur[id].xit[0], mur[id].J[0]);
//
//	Gdd = VanLeerB(mur[n4], mur[id].etax[0], mur[id].etay[0], mur[id].etat[0], mur[id].J[0]);
//	Gdu = VanLeerA(mur[n4], mur[id].etax[0], mur[id].etay[0], mur[id].etat[0], mur[id].J[0]);
//	Gcd = VanLeerB(mur[id], mur[id].etax[0], mur[id].etay[0], mur[id].etat[0], mur[id].J[0]);
//	Gcu = VanLeerA(mur[id], mur[id].etax[0], mur[id].etay[0], mur[id].etat[0], mur[id].J[0]);
//	Gud = VanLeerB(mur[n2], mur[id].etax[0], mur[id].etay[0], mur[id].etat[0], mur[id].J[0]);
//	Guu = VanLeerA(mur[n2], mur[id].etax[0], mur[id].etay[0], mur[id].etat[0], mur[id].J[0]);
//
//	U[0] = U[0] - dt * p.J[0] * (Fcr.f1 - Flr.f1 + Frl.f1 - Fcl.f1 + Gcu.f1 - Gdu.f1 + Gud.f1 - Gcd.f1);
//	U[1] = U[1] - dt * p.J[0] * (Fcr.f2 - Flr.f2 + Frl.f2 - Fcl.f2 + Gcu.f2 - Gdu.f2 + Gud.f2 - Gcd.f2);
//	U[2] = U[2] - dt * p.J[0] * (Fcr.f3 - Flr.f3 + Frl.f3 - Fcl.f3 + Gcu.f3 - Gdu.f3 + Gud.f3 - Gcd.f3);
//	U[3] = U[3] - dt * p.J[0] * (Fcr.f4 - Flr.f4 + Frl.f4 - Fcl.f4 + Gcu.f4 - Gdu.f4 + Gud.f4 - Gcd.f4);
//	p.rho = U[0];
//	p.u.x = U[1] / U[0];
//	p.u.y = U[2] / U[0];
//	p.p = (gama - 1) * (U[3] - 0.5 * p.rho * (p.u.x * p.u.x + p.u.y * p.u.y));
//}
void update_bound()
{
	extern vector <Mesh*>bound;
	extern vector<Mesh> ap;
	extern int xnum, ynum, znum;
	extern double gama;
	int i;
	int id;

	for (i = 0; i < bound.size(); i++)
	{
		id = bound[i]->id;
		if (bound[i]->type == "X-")
		{
			bound[i]->rho = 10;
			bound[i]->p = 10;
			bound[i]->u.x = 13 * sqrt(gama * bound[i]->p / bound[i]->rho);
			bound[i]->u.y = 0;
			bound[i]->u.z= 0;
		}
		else if (bound[i]->type == "X+")
		{
			bound[i]->rho = ap[id - 1].rho;
			bound[i]->u.x = ap[id - 1].u.x;
			bound[i]->u.y = ap[id - 1].u.y;
			bound[i]->u.z = ap[id - 1].u.z;
			bound[i]->p = ap[id - 1].p;
		}
		else if (bound[i]->type == "Y-")
		{
			bound[i]->rho = ap[id +xnum].rho;
			bound[i]->u.x = ap[id + xnum].u.x;
			bound[i]->u.y = ap[id + xnum].u.y;
			bound[i]->u.z = ap[id + xnum].u.z;
			bound[i]->p = ap[id + xnum].p;

		}
		else if (bound[i]->type == "Y+")
		{
			bound[i]->rho = ap[id - xnum].rho;
			bound[i]->u.x = ap[id - xnum].u.x;
			bound[i]->u.y = ap[id - xnum].u.y;
			bound[i]->u.z = ap[id - xnum].u.z;
			bound[i]->p = ap[id - xnum].p;

		}
		else if (bound[i]->type == "Z-")
		{
			bound[i]->rho = ap[id +xnum*ynum].rho;
			bound[i]->u.x = ap[id +xnum*ynum].u.x;
			bound[i]->u.y = ap[id +xnum*ynum].u.y;
			bound[i]->u.z = ap[id +xnum*ynum].u.z;
			bound[i]->p = ap[id +xnum*ynum].p;

		}
		else if (bound[i]->type == "Z+")
		{
			bound[i]->rho = ap[id -xnum*ynum].rho;
			bound[i]->u.x = ap[id -xnum*ynum].u.x;
			bound[i]->u.y = ap[id -xnum*ynum].u.y;
			bound[i]->u.z = ap[id -xnum*ynum].u.z;
			bound[i]->p = ap[id -xnum*ynum].p;

		}
		else if (bound[i]->type == "Body")
		{
			bound[i]->rho = 10;
			bound[i]->p = 10;
			bound[i]->u.x = 0;
			bound[i]->u.y = 0;
			bound[i]->u.z = 0;

		}
	}

	//for (i = 0; i < bb.size(); i++)
	//{
	//	double DELTA = 1e-10;
	//	double nx, ny, tx, ty, n1, n2;
	//	nx = bb[i]->neighbor[0]->x - bb[i]->x;
	//	ny = bb[i]->neighbor[0]->y - bb[i]->y;
	//	tx = -ny, ty = nx;
	//	n1 = bb[i]->neighbor[0]->u;
	//	n2 = bb[i]->neighbor[0]->v;
	//	if ((abs(n1) < DELTA && abs(n2) < DELTA)/* || (abs(tx) < DELTA || abs(ty) < DELTA)*/)
	//	{
	//		bb[i]->rho = bb[i]->neighbor[0]->rho;
	//		bb[i]->u = bb[i]->neighbor[0]->u;
	//		bb[i]->v = bb[i]->neighbor[0]->v;
	//		bb[i]->p = bb[i]->neighbor[0]->p;
	//	}
	//	else
	//	{
	//		double costheta = (tx * n1 + ty * n2) / (sqrt(tx * tx + ty * ty) * sqrt(n1 * n1 + n2 * n2));
	//		double ex = tx / sqrt(tx * tx + ty * ty);
	//		double ey = ty / sqrt(tx * tx + ty * ty);
	//		double u = sqrt(n1 * n1 + n2 * n2) * costheta * ex;
	//		double v = sqrt(n1 * n1 + n2 * n2) * costheta * ey;
	//		//std::cout << costheta << "  " << ex << "  " << ey << std::endl;
	//		//std::cout << A0[i].neighbor[0]->u << "  " << A0[i].neighbor[0]->v << "  " << u << "  " << v << std::endl;
	//		bb[i]->rho = bb[i]->neighbor[0]->rho;
	//		bb[i]->u = u;
	//		bb[i]->v = v;
	//		bb[i]->p = bb[i]->neighbor[0]->p;
	//	}
	//}

}


void sortPoint()
//put mesh point into different arrays
{
	extern vector<Mesh> ap;
	extern vector<Mesh*> ms;
	extern vector<Mesh*> mu;
	extern vector<Mesh*> bound;

	int n;
	for (int i = 0; i < ap.size(); i++)
	{
		if (ap[i].type == "IN")
		{
			n = 0;
			for (int j = 0; j < ap[i].neighbor.size(); j++)
				if (ap[i].neighbor[j]->type == "Body")
					n++;
			if (n == 0)
				ms.push_back(&ap[i]);
			else
				mu.push_back(&ap[i]);
		}
		else
			bound.push_back(&ap[i]);
	}
}

//void polymesh()
////get polygon mesh from grid points
//{
//	extern vector<mesh> AP;
//	extern vector<polygon_mesh>PM;
//	polygon_mesh Ptemp;
//	int i, j;
//	double DELTA = 1e-10;
//	int n1, n2, n3, n4, n5;
//	int size;
//	int id, id1, id2;
//	for (i = 0; i < Xnum * Ynum; i++)
//	{
//		if (i + Xnum + 1 < Xnum * Ynum && abs(AP[i + 1].x - AP[i].x - dx) < DELTA)
//		{
//			n1 = i;
//			n2 = i + 1;
//			n3 = i + Xnum + 1;
//			n4 = i + Xnum;
//			//if only one point is out of body
//			if ((AP[n1].section != 0 && AP[n2].section == 0 && AP[n3].section == 0 && AP[n4].section == 0) ||
//				(AP[n1].section == 0 && AP[n2].section != 0 && AP[n3].section == 0 && AP[n4].section == 0) ||
//				(AP[n1].section == 0 && AP[n2].section == 0 && AP[n3].section != 0 && AP[n4].section == 0) ||
//				(AP[n1].section == 0 && AP[n2].section == 0 && AP[n3].section == 0 && AP[n4].section != 0) ||
//				(AP[n1].section == 0 && AP[n2].section == 0 && AP[n3].section == 0 && AP[n4].section == 0))
//				continue;
//			//if all four points are out of body
//			else if (AP[n1].section != 0 && AP[n2].section != 0 && AP[n3].section != 0 && AP[n4].section != 0)
//			{
//				PM.push_back(Ptemp);
//				size = PM.size() - 1;
//				PM[size].node.push_back(n1);
//				PM[size].node.push_back(n2);
//				PM[size].node.push_back(n3);
//				PM[size].node.push_back(n4);
//				PM[size].face_start.push_back(n1);
//				PM[size].face_start.push_back(n2);
//				PM[size].face_start.push_back(n3);
//				PM[size].face_start.push_back(n4);
//				PM[size].face_end.push_back(n2);
//				PM[size].face_end.push_back(n3);
//				PM[size].face_end.push_back(n4);
//				PM[size].face_end.push_back(n1);
//			}
//			else if (AP[n1].section == 0 && AP[n2].section != 0 && AP[n3].section != 0 && AP[n4].section != 0)
//			{
//				id = AP[n2].connectId;
//				for (j = 0; j < AP[id].neighbor.size(); j++)
//				{
//					if (AP[id].neighbor[j]->type == "Body")
//						id1 = AP[id].neighbor[j]->id;
//				}
//				id = AP[n4].connectId;
//				for (j = 0; j < AP[id].neighbor.size(); j++)
//				{
//					if (AP[id].neighbor[j]->type == "Body")
//						id2 = AP[id].neighbor[j]->id;
//				}
//				if (id1 == id2)
//				{
//					n1 = id1;
//					PM.push_back(Ptemp);
//					size = PM.size() - 1;
//					PM[size].node.push_back(n1);
//					PM[size].node.push_back(n2);
//					PM[size].node.push_back(n3);
//					PM[size].node.push_back(n4);
//					PM[size].face_start.push_back(n1);
//					PM[size].face_start.push_back(n2);
//					PM[size].face_start.push_back(n3);
//					PM[size].face_start.push_back(n4);
//					PM[size].face_end.push_back(n2);
//					PM[size].face_end.push_back(n3);
//					PM[size].face_end.push_back(n4);
//					PM[size].face_end.push_back(n1);
//				}
//				else
//				{
//					n1 = id1;
//					n5 = id2;
//					//judge whether n1 and n5 are already be marked as neighbors
//					int n = 0;
//					for (j = 0; j < AP[n1].neighbor.size(); j++)
//					{
//						if (AP[n5].id == AP[n1].neighbor[j]->id)
//						{
//							n++;
//							break;
//						}
//					}
//					if (n == 0)
//						AP[n1].neighbor.push_back(&AP[n5]);
//					n = 0;
//					for (j = 0; j < AP[n5].neighbor.size(); j++)
//					{
//						if (AP[n1].id == AP[n5].neighbor[j]->id)
//						{
//							n++;
//							break;
//						}
//					}
//					if (n == 0)
//						AP[n5].neighbor.push_back(&AP[n1]);
//					PM.push_back(Ptemp);
//					size = PM.size() - 1;
//					PM[size].node.push_back(n1);
//					PM[size].node.push_back(n2);
//					PM[size].node.push_back(n3);
//					PM[size].node.push_back(n4);
//					PM[size].node.push_back(n5);
//					PM[size].face_start.push_back(n1);
//					PM[size].face_start.push_back(n2);
//					PM[size].face_start.push_back(n3);
//					PM[size].face_start.push_back(n4);
//					PM[size].face_start.push_back(n5);
//					PM[size].face_end.push_back(n2);
//					PM[size].face_end.push_back(n3);
//					PM[size].face_end.push_back(n4);
//					PM[size].face_end.push_back(n5);
//					PM[size].face_end.push_back(n1);
//				}
//			}
//			else if (AP[n1].section != 0 && AP[n2].section == 0 && AP[n3].section != 0 && AP[n4].section != 0)
//			{
//				id = AP[n3].connectId;
//				for (j = 0; j < AP[id].neighbor.size(); j++)
//				{
//					if (AP[id].neighbor[j]->type == "Body")
//						id1 = AP[id].neighbor[j]->id;
//				}
//				id = AP[n1].connectId;
//				for (j = 0; j < AP[id].neighbor.size(); j++)
//				{
//					if (AP[id].neighbor[j]->type == "Body")
//						id2 = AP[id].neighbor[j]->id;
//				}
//				if (id1 == id2)
//				{
//					n2 = id1;
//					PM.push_back(Ptemp);
//					size = PM.size() - 1;
//					PM[size].node.push_back(n1);
//					PM[size].node.push_back(n2);
//					PM[size].node.push_back(n3);
//					PM[size].node.push_back(n4);
//					PM[size].face_start.push_back(n1);
//					PM[size].face_start.push_back(n2);
//					PM[size].face_start.push_back(n3);
//					PM[size].face_start.push_back(n4);
//					PM[size].face_end.push_back(n2);
//					PM[size].face_end.push_back(n3);
//					PM[size].face_end.push_back(n4);
//					PM[size].face_end.push_back(n1);
//				}
//				else
//				{
//					n2 = id1;
//					n5 = id2;
//					//judge whether n1 and n5 are already be marked as neighbors
//					int n = 0;
//					for (j = 0; j < AP[n2].neighbor.size(); j++)
//					{
//						if (AP[n5].id == AP[n2].neighbor[j]->id)
//						{
//							n++;
//							break;
//						}
//					}
//					if (n == 0)
//						AP[n2].neighbor.push_back(&AP[n5]);
//					n = 0;
//					for (j = 0; j < AP[n5].neighbor.size(); j++)
//					{
//						if (AP[n2].id == AP[n5].neighbor[j]->id)
//						{
//							n++;
//							break;
//						}
//					}
//					if (n == 0)
//						AP[n5].neighbor.push_back(&AP[n2]);
//					PM.push_back(Ptemp);
//					size = PM.size() - 1;
//					PM[size].node.push_back(n1);
//					PM[size].node.push_back(n5);
//					PM[size].node.push_back(n2);
//					PM[size].node.push_back(n3);
//					PM[size].node.push_back(n4);
//					PM[size].face_start.push_back(n1);
//					PM[size].face_start.push_back(n5);
//					PM[size].face_start.push_back(n2);
//					PM[size].face_start.push_back(n3);
//					PM[size].face_start.push_back(n4);
//					PM[size].face_end.push_back(n5);
//					PM[size].face_end.push_back(n2);
//					PM[size].face_end.push_back(n3);
//					PM[size].face_end.push_back(n4);
//					PM[size].face_end.push_back(n1);
//				}
//			}
//			else if (AP[n1].section != 0 && AP[n2].section != 0 && AP[n3].section == 0 && AP[n4].section != 0)
//			{
//				id = AP[n4].connectId;
//				for (j = 0; j < AP[id].neighbor.size(); j++)
//				{
//					if (AP[id].neighbor[j]->type == "Body")
//						id1 = AP[id].neighbor[j]->id;
//				}
//				id = AP[n2].connectId;
//				for (j = 0; j < AP[id].neighbor.size(); j++)
//				{
//					if (AP[id].neighbor[j]->type == "Body")
//						id2 = AP[id].neighbor[j]->id;
//				}
//				if (id1 == id2)
//				{
//					n3 = id1;
//					PM.push_back(Ptemp);
//					size = PM.size() - 1;
//					PM[size].node.push_back(n1);
//					PM[size].node.push_back(n2);
//					PM[size].node.push_back(n3);
//					PM[size].node.push_back(n4);
//					PM[size].face_start.push_back(n1);
//					PM[size].face_start.push_back(n2);
//					PM[size].face_start.push_back(n3);
//					PM[size].face_start.push_back(n4);
//					PM[size].face_end.push_back(n2);
//					PM[size].face_end.push_back(n3);
//					PM[size].face_end.push_back(n4);
//					PM[size].face_end.push_back(n1);
//				}
//				else
//				{
//					n3 = id1;
//					n5 = id2;
//					//judge whether n1 and n5 are already be marked as neighbors
//					int n = 0;
//					for (j = 0; j < AP[n3].neighbor.size(); j++)
//					{
//						if (AP[n5].id == AP[n3].neighbor[j]->id)
//						{
//							n++;
//							break;
//						}
//					}
//					if (n == 0)
//						AP[n3].neighbor.push_back(&AP[n5]);
//					n = 0;
//					for (j = 0; j < AP[n5].neighbor.size(); j++)
//					{
//						if (AP[n3].id == AP[n5].neighbor[j]->id)
//						{
//							n++;
//							break;
//						}
//					}
//					if (n == 0)
//						AP[n5].neighbor.push_back(&AP[n3]);
//					PM.push_back(Ptemp);
//					size = PM.size() - 1;
//					PM[size].node.push_back(n1);
//					PM[size].node.push_back(n2);
//					PM[size].node.push_back(n5);
//					PM[size].node.push_back(n3);
//					PM[size].node.push_back(n4);
//					PM[size].face_start.push_back(n1);
//					PM[size].face_start.push_back(n2);
//					PM[size].face_start.push_back(n5);
//					PM[size].face_start.push_back(n3);
//					PM[size].face_start.push_back(n4);
//					PM[size].face_end.push_back(n2);
//					PM[size].face_end.push_back(n5);
//					PM[size].face_end.push_back(n3);
//					PM[size].face_end.push_back(n4);
//					PM[size].face_end.push_back(n1);
//				}
//			}
//			else if (AP[n1].section != 0 && AP[n2].section != 0 && AP[n3].section != 0 && AP[n4].section == 0)
//			{
//				id = AP[n1].connectId;
//				for (j = 0; j < AP[id].neighbor.size(); j++)
//				{
//					if (AP[id].neighbor[j]->type == "Body")
//						id1 = AP[id].neighbor[j]->id;
//				}
//				id = AP[n3].connectId;
//				for (j = 0; j < AP[id].neighbor.size(); j++)
//				{
//					if (AP[id].neighbor[j]->type == "Body")
//						id2 = AP[id].neighbor[j]->id;
//				}
//				if (id1 == id2)
//				{
//					n4 = id1;
//					PM.push_back(Ptemp);
//					size = PM.size() - 1;
//					PM[size].node.push_back(n1);
//					PM[size].node.push_back(n2);
//					PM[size].node.push_back(n3);
//					PM[size].node.push_back(n4);
//					PM[size].face_start.push_back(n1);
//					PM[size].face_start.push_back(n2);
//					PM[size].face_start.push_back(n3);
//					PM[size].face_start.push_back(n4);
//					PM[size].face_end.push_back(n2);
//					PM[size].face_end.push_back(n3);
//					PM[size].face_end.push_back(n4);
//					PM[size].face_end.push_back(n1);
//				}
//				else
//				{
//					n4 = id1;
//					n5 = id2;
//					//judge whether n1 and n5 are already be marked as neighbors
//					int n = 0;
//					for (j = 0; j < AP[n4].neighbor.size(); j++)
//					{
//						if (AP[n5].id == AP[n4].neighbor[j]->id)
//						{
//							n++;
//							break;
//						}
//					}
//					if (n == 0)
//						AP[n4].neighbor.push_back(&AP[n5]);
//					n = 0;
//					for (j = 0; j < AP[n5].neighbor.size(); j++)
//					{
//						if (AP[n4].id == AP[n5].neighbor[j]->id)
//						{
//							n++;
//							break;
//						}
//					}
//					if (n == 0)
//						AP[n5].neighbor.push_back(&AP[n4]);
//					PM.push_back(Ptemp);
//					size = PM.size() - 1;
//					PM[size].node.push_back(n1);
//					PM[size].node.push_back(n2);
//					PM[size].node.push_back(n3);
//					PM[size].node.push_back(n5);
//					PM[size].node.push_back(n4);
//					PM[size].face_start.push_back(n1);
//					PM[size].face_start.push_back(n2);
//					PM[size].face_start.push_back(n3);
//					PM[size].face_start.push_back(n5);
//					PM[size].face_start.push_back(n4);
//					PM[size].face_end.push_back(n2);
//					PM[size].face_end.push_back(n3);
//					PM[size].face_end.push_back(n5);
//					PM[size].face_end.push_back(n4);
//					PM[size].face_end.push_back(n1);
//				}
//			}
//			else
//			{
//				if (AP[n1].section == 0 && AP[n2].section == 0 && AP[n3].section != 0 && AP[n4].section != 0)
//				{
//					id = AP[n3].connectId;
//					for (j = 0; j < AP[id].neighbor.size(); j++)
//					{
//						if (AP[id].neighbor[j]->type == "Body")
//							id1 = AP[id].neighbor[j]->id;
//					}
//					id = AP[n4].connectId;
//					for (j = 0; j < AP[id].neighbor.size(); j++)
//					{
//						if (AP[id].neighbor[j]->type == "Body")
//							id2 = AP[id].neighbor[j]->id;
//					}
//					n1 = id2;
//					n2 = id1;
//				}
//				else if (AP[n1].section != 0 && AP[n2].section == 0 && AP[n3].section == 0 && AP[n4].section != 0)
//				{
//					id = AP[n4].connectId;
//					for (j = 0; j < AP[id].neighbor.size(); j++)
//					{
//						if (AP[id].neighbor[j]->type == "Body")
//							id1 = AP[id].neighbor[j]->id;
//					}
//					id = AP[n1].connectId;
//					for (j = 0; j < AP[id].neighbor.size(); j++)
//					{
//						if (AP[id].neighbor[j]->type == "Body")
//							id2 = AP[id].neighbor[j]->id;
//					}
//					n2 = id2;
//					n3 = id1;
//				}
//				else if (AP[n1].section != 0 && AP[n2].section != 0 && AP[n3].section == 0 && AP[n4].section == 0)
//				{
//					id = AP[n1].connectId;
//					for (j = 0; j < AP[id].neighbor.size(); j++)
//					{
//						if (AP[id].neighbor[j]->type == "Body")
//							id1 = AP[id].neighbor[j]->id;
//					}
//					id = AP[n2].connectId;
//					for (j = 0; j < AP[id].neighbor.size(); j++)
//					{
//						if (AP[id].neighbor[j]->type == "Body")
//							id2 = (*AP[id].neighbor[j]).id;
//					}
//					n3 = id2;
//					n4 = id1;
//				}
//				else if (AP[n1].section == 0 && AP[n2].section != 0 && AP[n3].section != 0 && AP[n4].section == 0)
//				{
//					id = AP[n2].connectId;
//					for (j = 0; j < AP[id].neighbor.size(); j++)
//					{
//						if (AP[id].neighbor[j]->type == "Body")
//							id1 = AP[id].neighbor[j]->id;
//					}
//					id = AP[n3].connectId;
//					for (j = 0; j < AP[id].neighbor.size(); j++)
//					{
//						if (AP[id].neighbor[j]->type == "Body")
//							id2 = AP[id].neighbor[j]->id;
//					}
//					n4 = id2;
//					n1 = id1;
//				}
//				int n = 0;
//				for (j = 0; j < AP[id1].neighbor.size(); j++)
//				{
//					if (AP[id2].id == AP[id1].neighbor[j]->id)
//					{
//						n++;
//						break;
//					}
//				}
//				if (n == 0)
//					AP[id1].neighbor.push_back(&AP[id2]);
//				n = 0;
//				for (j = 0; j < AP[id2].neighbor.size(); j++)
//				{
//					if (AP[id1].id == AP[id2].neighbor[j]->id)
//					{
//						n++;
//						break;
//					}
//				}
//				if (n == 0)
//					AP[id2].neighbor.push_back(&AP[id1]);
//				PM.push_back(Ptemp);
//				size = PM.size() - 1;
//				PM[size].node.push_back(n1);
//				PM[size].node.push_back(n2);
//				PM[size].node.push_back(n3);
//				PM[size].node.push_back(n4);
//				PM[size].face_start.push_back(n1);
//				PM[size].face_start.push_back(n2);
//				PM[size].face_start.push_back(n3);
//				PM[size].face_start.push_back(n4);
//				PM[size].face_end.push_back(n2);
//				PM[size].face_end.push_back(n3);
//				PM[size].face_end.push_back(n4);
//				PM[size].face_end.push_back(n1);
//
//			}
//		}
//	}
//}