#include<cmath>
#include"include/const.h"
#include<iostream>
#include<omp.h>
#include"include/functions.h"
using namespace std;
extern double gama;
Flux VanLeer(Mesh& N, string direction, string side)//ÎÞ×ø±ê±ä»»
{
	double u = N.u.x;
	double v = N.u.y;
	double w = N.u.z;
	double p = N.p;
	double rho = N.rho;
	double c = sqrt(gama * p / rho);
	double U;
	if (direction == "X" || direction == "x")
		U = u;
	else if (direction == "Y" || direction == "y")
		U = v;
	else if (direction == "Z" || direction == "z")
		U = w;
	double M = abs(U / c);
	Flux F;
	if (side == "R")
	{
		if (U + c <= 0)
		{
			F = { 0 };
		}
		else if (U - c >= 0)
		{
			double E = 0.5 * rho * (u * u + v * v + w * w) + p / (gama - 1);
			F.f1 = rho * U;
			F.f2 = rho * U * u + p;
			F.f3 = rho * U * v + p;
			F.f4 = rho * U * w + p;
			F.f5 = U * (E + p) - p;
		}
		if (M <= 1)
		{
			double FM = 0.25 * rho * c * (M + 1) * (M + 1);
			F.f1 = FM;
			F.f2 = FM * ((-U + 2 * c) / gama + u);
			F.f3 = FM * ((-U + 2 * c) / gama + v);
			F.f4 = FM * ((-U + 2 * c) / gama + w);
			F.f5 = FM * (U * (-U + 2 * c) / (gama + 1) + 2 * c * c / (gama * gama - 1) + 0.5 * (u * u + v * v + w * w));
		}
	}
	else if (side == "L")
	{
		if (U - c >= 0)
		{
			F = { 0 };
		}
		else if (U + c <= 0)
		{
			double E = 0.5 * rho * (u * u + v * v + w * w) + p / (gama - 1);
			F.f1 = rho * U;
			F.f2 = rho * U * u + p;
			F.f3 = rho * U * v + p;
			F.f4 = rho * U * w + p;
			F.f5 = U * (E + p) - p;
		}
		else if (M <= 1)
		{
			double FM = -0.25 * rho * c * (M - 1) * (M - 1);
			F.f1 = FM;
			F.f2 = FM * ((-U - 2 * c) / gama + u);
			F.f3 = FM * ((-U - 2 * c) / gama + v);
			F.f4 = FM * ((-U - 2 * c) / gama + w);
			F.f5 = FM * (U * (-U - 2 * c) / (gama + 1) + 2 * c * c / (gama * gama - 1) + 0.5 * (u * u + v * v + w * w));
		}
	}
	return F;
}

Flux VanLeer(Mesh& N, Mesh& C, string direction, string side)
{
	double u = N.u.x;
	double v = N.u.y;
	double w = N.u.z;
	double p = N.p;
	double rho = N.rho;
	double c = sqrt(gama * p / rho);
	double phix, phiy, phiz, phit, J;
	if (direction == "xi")
		phix = C.xix, phiy = C.xiy, phiz = C.xiz, phit = C.xit, J = C.J;
	if (direction == "eta")
		phix = C.etax, phiy = C.etay, phiz = C.etaz, phit = C.etat, J = C.J;
	if (direction == "zeta")
		phix = C.zetax, phiy = C.zetay, phiz = C.zetaz, phit = C.zetat, J = C.J;
	double Dphi = sqrt(phix * phix + phiy * phiy + phiz * phiz);
	double ub = u * phix + v * phiy + w * phiz + phit;
	double uc = ub / Dphi;
	phix = phix / Dphi;
	phiy = phiy / Dphi;
	phiz = phiy / Dphi;
	phit = phit / Dphi;
	double M = abs(uc / c);
	Flux F;
	if (side == "R")
	{
		if (uc + c <= 0)
		{
			F = { 0 };
		}
		else if (uc - c >= 0)
		{
			double E = 0.5 * rho * (u * u + v * v + w * w) + p / (gama - 1);
			F.f1 = rho * uc;
			F.f2 = rho * uc * u + p * phix;
			F.f3 = rho * uc * v + p * phiy;
			F.f4 = rho * uc * w + p * phiz;
			F.f5 = uc * (E + p) - p * phit;
		}
		else if (M <= 1)
		{
			double FM = 0.25 * rho * c * (M + 1) * (M + 1);
			F.f1 = FM;
			F.f2 = FM * (phix * (-uc + 2 * c) / gama + u);
			F.f3 = FM * (phiy * (-uc + 2 * c) / gama + v);
			F.f4 = FM * (phiz * (-uc + 2 * c) / gama + w);
			F.f5 = FM * (uc * (-uc + 2 * c) / (gama + 1) + 2 * c * c / (gama * gama - 1) + 0.5 * (u * u + v * v + w * w) - phit * (-uc + 2 * c) / gama);
		}
	}
	else if (side == "L")
	{
		if (uc - c >= 0)
		{
			F = { 0 };
		}
		else if (uc + c <= 0)
		{
			double E = 0.5 * rho * (u * u + v * v + w * w) + p / (gama - 1);
			F.f1 = rho * uc;
			F.f2 = rho * uc * u + p * phix;
			F.f3 = rho * uc * v + p * phiy;
			F.f4 = rho * uc * w + p * phiz;
			F.f5 = uc * (E + p) - p * phit;
		}
		else if (M <= 1)
		{
			double FM = -0.25 * rho * c * (M - 1) * (M - 1);
			F.f1 = FM;
			F.f2 = FM * (phix * (-uc - 2 * c) / gama + u);
			F.f3 = FM * (phiy * (-uc - 2 * c) / gama + v);
			F.f4 = FM * (phiz * (-uc - 2 * c) / gama + w);
			F.f5 = FM * (uc * (-uc - 2 * c) / (gama + 1) + 2 * c * c / (gama * gama - 1) + 0.5 * (u * u + v * v + w * w) - phit * (-uc - 2 * c) / gama);
		}

	}
	F.f1 = F.f1 * Dphi;
	F.f2 = F.f2 * Dphi;
	F.f3 = F.f3 * Dphi;
	F.f4 = F.f4 * Dphi;
	F.f5 = F.f5 * Dphi;
	return F;
}
