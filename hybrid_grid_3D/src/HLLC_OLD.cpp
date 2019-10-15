//#include<cmath>
//#include"const.h"
//#include<iostream>
//#include<omp.h>
//
//using namespace std;
////2άŷ����������ɢ���ŵ���������������ѧ�̡̳�p433 13.3.7
//double * HLLC_��(mesh CL, mesh CR, mesh C, int method)//�������Ҳ��㣬����任�ο���
//{
//	double xix = C.xix[method];//�˴�����xix,��xix*J-1=yeta
//	double xiy = C.xiy[method];
//	double J = C.J[method];
//	double Dxi = sqrt(xix*xix + xiy * xiy);
//	double xicL = (xix * CL.u + xiy * CL.v) / Dxi;
//	double xicR = (xix * CR.u + xiy * CR.v) / Dxi;
//	double aL = sqrt(gama*CL.p / CL.rho);
//	double aR = sqrt(gama*CR.p / CR.rho);
//	double EL = CL.p / (gama - 1) + 0.5*CL.rho * CL.u * CL.u + 0.5*CL.rho * CL.v * CL.v;
//	double ER = CR.p / (gama - 1) + 0.5*CR.rho * CR.u * CR.u + 0.5*CR.rho * CR.v * CR.v;
//	double FL[4], FR[4];
//
//	FL[0] = (CL.rho*xicL);
//	FL[1] = (CL.rho*CL.u*xicL + xix * CL.p);
//	FL[2] = (CL.rho*CL.v*xicL + xiy * CL.p);
//	FL[3] = (xicL * (EL + CL.p));
//	FR[0] = (CR.rho*xicR);
//	FR[1] = (CR.rho*CR.u*xicR + xix * CR.p);
//	FR[2] = (CR.rho*CR.v*xicR + xiy * CR.p);
//	FR[3] = (xicR * (ER + CR.p));
//	double SL = xicL - aL;
//	double SR = xicR + aR;
//	double SM = (CL.p - CR.p + CL.rho * xicL * (xicL - SL) + CR.rho * xicR* (SR - xicR)) / (CR.rho * (SR - xicR) - CL.rho * (SL - xicL));
//	double pM = 0.5*(CL.rho * (xicL - SL)*(xicL - SM) + CR.rho * (xicR - SR)*(xicR - SM) + CR.p + CL.p);
//	//Average-State Jacobians and Implicit Methods for Compressible Viscous and Turbulent Flows,(7)(8)
//	double ��L = 1 / (SL - SM);
//	double ��R = 1 / (SR - SM);
//
//	double rhoLM = ��L * CL.rho*(SL - xicL);
//	double rhouLS = ��L * ((SL - xicL)*(CL.rho*CL.u) + (pM - CL.p)*xix);
//	double rhovLS = ��L * ((SL - xicL)*(CL.rho*CL.v) + (pM - CL.p)*xiy);
//	double eLS = ��L * ((SL - xicL)*EL - CL.p*xicL + pM * SM);
//
//	double rhoRM = ��R * CR.rho*(SR - xicR);
//	double rhouRS = ��R * ((SR - xicR)*(CR.rho*CR.u) + (pM - CR.p)*xix);
//	double rhovRS = ��R * ((SR - xicR)*(CR.rho*CR.v) + (pM - CR.p)*xiy);
//	double eRS = ��R * ((SR - xicR)*ER - CR.p*xicR + pM * SM);
//
//	double  FML[4], FMR[4];
//	FML[0] = SM * rhoLM;
//	FML[1] = rhouLS * SM + pM * xix;
//	FML[2] = rhovLS * SM + pM * xiy;
//	FML[3] = (eLS + pM)*SM;
//	FMR[0] = SM * rhoRM;
//	FMR[1] = rhouRS * SM + pM * xix;
//	FMR[2] = rhovRS * SM + pM * xiy;
//	FMR[3] = (eRS + pM)*SM;
//
//	double *F_HLLC = new double[4];
//	if (SL > 0)
//		for (int i = 0; i < 4; i++)
//			F_HLLC[i] = FL[i];
//	else if (SL <= 0 && SM > 0)
//		for (int i = 0; i < 4; i++)
//			F_HLLC[i] = FML[i];
//	else if (SM <= 0 && SR >= 0)
//		for (int i = 0; i < 4; i++)
//			F_HLLC[i] = FMR[i];
//	else
//		for (int i = 0; i < 4; i++)
//			F_HLLC[i] = FR[i];
//	for (int i = 0; i < 4; i++)
//		F_HLLC[i] = F_HLLC[i] *Dxi;
//	return F_HLLC;
//}
//
//double * HLLC_��(mesh CD, mesh CU, mesh C, int method)
//{
//	double etax = C.etax[method];
//	double etay = C.etay[method];
//	double J = C.J[method];
//	double Deta = sqrt(etax*etax + etay * etay);
//	double etacU = (C.etax[method] * CU.u + C.etay[method] * CD.v) / Deta;
//	double etacD = (C.etax[method] * CD.u + C.etay[method] * CU.v) / Deta;
//	double aU = sqrt(gama*CU.p / CU.rho);
//	double aD = sqrt(gama*CD.p / CD.rho);
//	double EU = CU.p / (gama - 1) + 0.5*CU.rho * CU.u * CU.u + 0.5*CU.rho * CU.v * CU.v;
//	double ED = CD.p / (gama - 1) + 0.5*CD.rho * CD.u * CD.u + 0.5*CD.rho * CD.v * CD.v;
//	double GD[4], GU[4];
//
//	GD[0] = (CD.rho*etacD);
//	GD[1] = (CD.rho*CD.u*etacD + etax * CD.p);
//	GD[2] = (CD.rho*CD.v*etacD + etay * CD.p);
//	GD[3] = (etacD * (ED + CD.p));
//
//	GU[0] = (CU.rho*etacU);
//	GU[1] = (CU.rho*CU.u*etacU + etax * CD.p);
//	GU[2] = (CU.rho*CU.v*etacU + etay * CD.p);
//	GU[3] = (etacU * (EU + CU.p));
//
//	double SD = etacD - aD;
//	double SU = etacU + aU;
//	double SM = (CD.p - CU.p - CD.rho * etax * (SD - etacD) + CU.rho * etax* (SU - etacU)) / (CU.rho * (SU - etacU) - CD.rho * (SD - etacD));
//	double pM = 0.5*(CD.rho * (etacD - SD)*(etacD - SM) + CU.rho * (etacU - SU)*(etacU - SM) + CU.p + CD.p);
//	//Average-State Jacobians and Implicit Methods for Compressible Viscous and Turbulent Flows,(7)(8)
//	double ��D = 1 / (SD - SM);
//	double ��U = 1 / (SU - SM);
//
//	double rhoDM = ��D * CD.rho*(SD - etacD);
//	double rhouDS = ��D * ((SD - etacD)*(CD.rho*CD.u) + (pM - CD.p)*etax);
//	double rhovDS = ��D * ((SD - etacD)*(CD.rho*CD.v) + (pM - CD.p)*etay);
//	double eDS = ��D * ((SD - etacD)*ED - CD.p*etacD + pM * SM);
//
//	double rhoUM = ��U * CU.rho*(SU - etacU);
//	double rhouUS = ��U * ((SU - etacU)*(CU.rho*CU.u) + (pM - CU.p)*etax);
//	double rhovUS = ��U * ((SU - etacU)*(CU.rho*CU.v) + (pM - CU.p)*etay);
//	double eUS = ��U * ((SU - etacU)*EU - CU.p*etacU + pM * SM);
//
//	double  FMD[4], FMU[4];
//
//	FMD[0] = SM * rhoDM;
//	FMD[1] = rhouDS * SM + pM * etax;
//	FMD[2] = rhovDS * SM + pM * etay;
//	FMD[3] = (eDS + pM)*SM;
//	FMU[0] = SM * rhoUM;
//	FMU[1] = rhouUS * SM + pM * etax;
//	FMU[2] = rhovUS * SM + pM * etay;
//	FMU[3] = (eUS + pM)*SM;
//
//	double *G_HLLC = new double[4];
//	if (SD >= 0)
//		for (int i = 0; i < 4; i++)
//			G_HLLC[i] = GD[i];
//	else if (SD <= 0 && SM >= 0)
//		for (int i = 0; i < 4; i++)
//			G_HLLC[i] = FMD[i];
//	else if (SM <= 0 && SU >= 0)
//		for (int i = 0; i < 4; i++)
//			G_HLLC[i] = FMU[i];
//	else
//		for (int i = 0; i < 4; i++)
//			G_HLLC[i] = GU[i];
//	for (int i = 0; i < 4; i++)
//		G_HLLC[i] = G_HLLC[i] *Deta;
//
//	return G_HLLC;
//}
