#include"include/const.h"
#include"include/functions.h"
#include<vector>
using std::vector;
void coordinate_trans()
{
	extern vector<Mesh*> mu;
	extern vector<Mesh*> ms;
	int i, j, k;
	int n1, n2, n3, n4;
	extern double t_sim;
	extern double dt;
	for (i = 0; i < mu.size(); i++)
	{
		if (mu[i]->neighbor.size() == 6)
		{
			mu[i]->xxi = 0.5 * (mu[i]->neighbor[0]->x - mu[i]->neighbor[1]->x);
			mu[i]->yxi = 0.5 * (mu[i]->neighbor[0]->y - mu[i]->neighbor[1]->y);
			mu[i]->zxi = 0.5 * (mu[i]->neighbor[0]->z - mu[i]->neighbor[1]->z);
			mu[i]->xeta = 0.5 * (mu[i]->neighbor[2]->x - mu[i]->neighbor[3]->x);
			mu[i]->yeta = 0.5 * (mu[i]->neighbor[2]->y - mu[i]->neighbor[3]->y);
			mu[i]->zeta = 0.5 * (mu[i]->neighbor[2]->z - mu[i]->neighbor[3]->z);
			mu[i]->xzeta = 0.5 * (mu[i]->neighbor[4]->x - mu[i]->neighbor[5]->x);
			mu[i]->yzeta = 0.5 * (mu[i]->neighbor[4]->y - mu[i]->neighbor[5]->y);
			mu[i]->zzeta = 0.5 * (mu[i]->neighbor[4]->z - mu[i]->neighbor[5]->z);
			mu[i]->xtau = 0;
			mu[i]->ytau = 0;
			mu[i]->ztau = 0;
			mu[i]->J = mu[i]->xxi * (mu[i]->yeta * mu[i]->zzeta - mu[i]->zeta * mu[i]->yzeta)
				- mu[i]->yxi * (mu[i]->xxi * mu[i]->zzeta - mu[i]->zxi * mu[i]->yzeta)
				+ mu[i]->zxi * (mu[i]->xeta * mu[i]->yzeta - mu[i]->yeta * mu[i]->xzeta);
			mu[i]->J = 1.0 / mu[i]->J;

			mu[i]->xix = mu[i]->J * (mu[i]->yeta * mu[i]->zzeta - mu[i]->yzeta * mu[i]->zeta);
			mu[i]->xiy = mu[i]->J * (mu[i]->xzeta * mu[i]->zeta - mu[i]->xeta * mu[i]->zzeta);
			mu[i]->xiz = mu[i]->J * (mu[i]->xeta * mu[i]->zzeta - mu[i]->xxi * mu[i]->zeta);
			mu[i]->etax = mu[i]->J * (mu[i]->yzeta * mu[i]->zxi - mu[i]->yxi * mu[i]->zzeta);
			mu[i]->etay = mu[i]->J * (mu[i]->xxi * mu[i]->zzeta - mu[i]->xzeta * mu[i]->zxi);
			mu[i]->etaz = mu[i]->J * (mu[i]->xzeta * mu[i]->yxi - mu[i]->xxi * mu[i]->yzeta);
			mu[i]->zetax = mu[i]->J * (mu[i]->yxi * mu[i]->zeta - mu[i]->yeta * mu[i]->zxi);
			mu[i]->zetay = mu[i]->J * (mu[i]->xeta * mu[i]->zxi - mu[i]->xxi * mu[i]->xeta);
			mu[i]->zetaz = mu[i]->J * (mu[i]->xxi * mu[i]->yeta - mu[i]->xeta * mu[i]->yxi);
			mu[i]->xit = -(mu[i]->xtau * mu[i]->xix + mu[i]->ytau * mu[i]->xiy + mu[i]->ztau * mu[i]->xiz);
			mu[i]->etat = -(mu[i]->xtau * mu[i]->etax + mu[i]->ytau * mu[i]->etay + mu[i]->ztau * mu[i]->etaz);
			mu[i]->zetat = -(mu[i]->xtau * mu[i]->zetax + mu[i]->ytau * mu[i]->zetay + mu[i]->ztau * mu[i]->zetaz);

		}
	}
	for (i = 0; i < ms.size(); i++)
	{
		if (ms[i]->neighbor.size() == 6)
		{
			ms[i]->xxi = 0.5 * (ms[i]->neighbor[0]->x - ms[i]->neighbor[1]->x);
			ms[i]->yxi = 0.5 * (ms[i]->neighbor[0]->y - ms[i]->neighbor[1]->y);
			ms[i]->zxi = 0.5 * (ms[i]->neighbor[0]->z - ms[i]->neighbor[1]->z);
			ms[i]->xeta = 0.5 * (ms[i]->neighbor[2]->x - ms[i]->neighbor[3]->x);
			ms[i]->yeta = 0.5 * (ms[i]->neighbor[2]->y - ms[i]->neighbor[3]->y);
			ms[i]->zeta = 0.5 * (ms[i]->neighbor[2]->z - ms[i]->neighbor[3]->z);
			ms[i]->xzeta = 0.5 * (ms[i]->neighbor[4]->x - ms[i]->neighbor[5]->x);
			ms[i]->yzeta = 0.5 * (ms[i]->neighbor[4]->y - ms[i]->neighbor[5]->y);
			ms[i]->zzeta = 0.5 * (ms[i]->neighbor[4]->z - ms[i]->neighbor[5]->z);
			ms[i]->xtau = 0;
			ms[i]->ytau = 0;
			ms[i]->ztau = 0;
			ms[i]->J = ms[i]->xxi * (ms[i]->yeta * ms[i]->zzeta - ms[i]->zeta * ms[i]->yzeta)
				- ms[i]->yxi * (ms[i]->xxi * ms[i]->zzeta - ms[i]->zxi * ms[i]->yzeta)
				+ ms[i]->zxi * (ms[i]->xeta * ms[i]->yzeta - ms[i]->yeta * ms[i]->xzeta);
			ms[i]->J = 1.0 / ms[i]->J;

			ms[i]->xix = ms[i]->J * (ms[i]->yeta * ms[i]->zzeta - ms[i]->yzeta * ms[i]->zeta);
			ms[i]->xiy = ms[i]->J * (ms[i]->xzeta * ms[i]->zeta - ms[i]->xeta * ms[i]->zzeta);
			ms[i]->xiz = ms[i]->J * (ms[i]->xeta * ms[i]->zzeta - ms[i]->xxi * ms[i]->zeta);
			ms[i]->etax = ms[i]->J * (ms[i]->yzeta * ms[i]->zxi - ms[i]->yxi * ms[i]->zzeta);
			ms[i]->etay = ms[i]->J * (ms[i]->xxi * ms[i]->zzeta - ms[i]->xzeta * ms[i]->zxi);
			ms[i]->etaz = ms[i]->J * (ms[i]->xzeta * ms[i]->yxi - ms[i]->xxi * ms[i]->yzeta);
			ms[i]->zetax = ms[i]->J * (ms[i]->yxi * ms[i]->zeta - ms[i]->yeta * ms[i]->zxi);
			ms[i]->zetay = ms[i]->J * (ms[i]->xeta * ms[i]->zxi - ms[i]->xxi * ms[i]->xeta);
			ms[i]->zetaz = ms[i]->J * (ms[i]->xxi * ms[i]->yeta - ms[i]->xeta * ms[i]->yxi);
			ms[i]->xit = -(ms[i]->xtau * ms[i]->xix + ms[i]->ytau * ms[i]->xiy + ms[i]->ztau * ms[i]->xiz);
			ms[i]->etat = -(ms[i]->xtau * ms[i]->etax + ms[i]->ytau * ms[i]->etay + ms[i]->ztau * ms[i]->etaz);
			ms[i]->zetat = -(ms[i]->xtau * ms[i]->zetax + ms[i]->ytau * ms[i]->zetay + ms[i]->ztau * ms[i]->zetaz);

		}
	}

}

