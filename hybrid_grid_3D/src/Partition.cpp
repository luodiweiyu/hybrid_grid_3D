//Partition file for computing , divid the original mesh into several regions
#include<iostream>
#include"include/const.h"
#include"include/functions.h"
#include<vector>
#include<string>
#include<fstream>
using std::vector;


void partition_Point()//Partition existing grid points
{
	extern int xnum, ynum, znum;
	extern string methodType;
	if (methodType == "C")
	{
		int i;
		if (flowType == "cylinder")
		{
			extern vector <Coordinate> poly;
			extern vector<Mesh>ap;
			double r1;
			Mesh tempu;
			int sizeA, size;
			sizeA = ap.size();
			for (i = 0; i < sizeA; i++)
			{
				ap[i].id = i;
				int n[6] = { 0 };
				int sum = 0;
				if (ap[i].type != "IN")
					continue;
				else if (ap[i].type == "IN")
				{
					ap[i].neighbor.push_back(&ap[i + 1]);
					ap[i].neighbor.push_back(&ap[i - 1]);
					ap[i].neighbor.push_back(&ap[i + xnum]);
					ap[i].neighbor.push_back(&ap[i - xnum]);
					ap[i].neighbor.push_back(&ap[i + xnum * ynum]);
					ap[i].neighbor.push_back(&ap[i - xnum * ynum]);
				}
				if (judgeFieldInOut(ap[i], poly))
					ap[i].section = 0, ap[i].sec_num = 0;
				else
				{
					for (int j = 0; j < 6; j++)
					{
						if (!judgeFieldInOut(*ap[i].neighbor[j], poly))
							//std::cout << ap[i].neighbor[0] << std::endl;
							n[j] = 1;
						sum += n[j];
					}
					if (sum == 6)
						ap[i].section = 1, ap[i].sec_num = 1;
					else
					{
						double r1 = r;
						int near_id = findNearPoint(ap[i], poly);
						ap[i].section = -1, ap[i].sec_num = 0;
						ap.push_back(ap[i]);
						size = ap.size() - 1;
						ap[size].id = size;
						ap[size].connectId = i;
						ap[i].connectId = size;
						ap.push_back(tempu);
						ap[ap.size() - 1] = poly[near_id];
						ap[ap.size() - 1].id = ap.size() - 1;
						ap[ap.size() - 1].type = "Body";
						ap[ap.size() - 1].neighbor.push_back(&ap[size]);
						ap[ap.size() - 1].section = -2;
						if (sum == 5)
						{
							for (int j = 0; j < 6; j++)
							{
								if (n[j] == 0)
									ap[size].neighbor[j] = &ap[ap.size() - 1];
							}
						}
						else if (sum == 4)
						{
							int n1, n2;
							if (n[0] == 0)
							{
								n1 = 0;
								if (n[4] == 0)
									n2 = 4;
								else if (n[2] == 0)
									n2 = 2;
								else if (n[5] == 0)
									n2 = 5;
								else if (n[3] == 0)
									n2 = 3;
								else if (n[1] == 0)
									std::cout << "对侧坐标0,1为零，不符逻辑！" << std::endl;
							}
							else if (n[1] == 0)
							{
								n1 = 1;
								if (n[4] == 0)
									n2 = 4;
								else if (n[2] == 0)
									n2 = 2;
								else if (n[5] == 0)
									n2 = 5;
								else if (n[3] == 0)
									n2 = 3;
							}
							else if (n[2] == 0)
							{
								n1 = 2;
								if (n[4] == 0)
									n2 = 4;
								else if (n[5] == 0)
									n2 = 5;
								else if (n[3] == 0)
									std::cout << "对侧坐标2,3为零，不符逻辑！" << std::endl;
							}
							else if (n[3] == 0)
							{
								n1 = 3;
								if (n[4] == 0)
									n2 = 4;
								else if (n[5] == 0)
									n2 = 5;
							}
							else if (n[4] == 0)
							{
								n1 = 4;
								if (n[5] == 0)
									std::cout << "对侧坐标4,5为零，不符逻辑！" << std::endl;
							}
							else
								std::cout << "sum != 0!!" << std::endl;
							ap[size].neighbor[n1] = &ap[ap.size() - 1];
							ap[size].neighbor[n2] = &ap[ap.size() - 1];
						}
						else if (sum == 3)
						{
							int n1, n2, n3;
							if (n[0] + n[2] + n[4] == 0)
								n1 = 0, n2 = 2, n3 = 4;
							else if (n[0] + n[2] + n[5] == 0)
								n1 = 0, n2 = 2, n3 = 5;
							else if (n[0] + n[5] + n[3] == 0)
								n1 = 0, n2 = 5, n3 = 3;
							else if (n[0] + n[3] + n[4] == 0)
								n1 = 0, n2 = 3, n3 = 4;
							else if (n[1] + n[2] + n[4] == 0)
								n1 = 1, n2 = 2, n3 = 4;
							else if (n[1] + n[2] + n[5] == 0)
								n1 = 1, n2 = 2, n3 = 5;
							else if (n[1] + n[5] + n[3] == 0)
								n1 = 1, n2 = 5, n3 = 3;
							else if (n[1] + n[3] + n[4] == 0)
								n1 = 1, n2 = 3, n3 = 4;
							else
								std::cout << "存在一些奇怪的组合！" << std::endl;
							ap[size].neighbor[n1] = &ap[ap.size() - 1];
							ap[size].neighbor[n2] = &ap[ap.size() - 1];
							ap[size].neighbor[n3] = &ap[ap.size() - 1];
						}
						else
							std::cout << "sum = " << sum << "无法构成非结构点！" << std::endl;
					}
				}
			}
		}
	}


}

