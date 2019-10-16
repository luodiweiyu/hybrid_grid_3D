#include"include/const.h"
#include"include/functions.h"
#include<vector>
using std::vector;
void coordinate_trans()
{
	extern vector<Mesh*> mu;
	int i, j, k;
	int n1, n2, n3, n4;
	extern double t_sim;
	extern double dt;
	for (i = 0; i < mu.size(); i++)
	{
		mu[i]->xxi = 0.5 * (mu[i]->neibor[0]->x - mu[i]->neibor[1]->x);
		mu[i]->yxi = 0.5 * (mu[i]->neibor[0]->y - mu[i]->neibor[1]->y);
		mu[i]->zxi = 0.5 * (mu[i]->neibor[0]->z - mu[i]->neibor[1]->z);
		mu[i]->xeta = 0.5 * (mu[i]->neibor[2]->x - mu[i]->neibor[3]->x);
		mu[i]->yeta = 0.5 * (mu[i]->neibor[2]->y - mu[i]->neibor[3]->y);
		mu[i]->zeta = 0.5 * (mu[i]->neibor[2]->z - mu[i]->neibor[3]->z);
		mu[i]->xzeta = 0.5 * (mu[i]->neibor[4]->x - mu[i]->neibor[5]->x);
		mu[i]->yzeta = 0.5 * (mu[i]->neibor[4]->y - mu[i]->neibor[5]->y);
		mu[i]->zzeta = 0.5 * (mu[i]->neibor[4]->z - mu[i]->neibor[5]->z);
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

