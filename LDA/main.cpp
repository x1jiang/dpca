#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#define DEBUG
double * c1m1, * c2m1, * m2;
int * nom;
int p2, p1, p, n;
double prob;
inline double change(double x)//the change of L1 norm if we change one element from 0 to 1
{
	return (abs(1 + x)-abs(x));
}
inline double change(double x, double _m)//the change of L1 norm if we change one element from 0 to 1
{
	return (abs(_m + x)-abs(x));
}
void readstat()
{
	FILE * f = fopen("cppin.txt", "r");
	fscanf(f, "%d", &p);
	int n1, n2;
	for (int i = 1; i < p; i++)
		fscanf(f, "%d", &p);
	for (int i = 0; i < p; i++)
		fscanf(f, "%d", &n1);
	for (int i = 0; i < p; i++)
		fscanf(f, "%d", &n2);
	n = n1 + n2;
	prob = n1 / (double) n;
	nom = new int[p + 1];
	c1m1 = new double[p];
	c2m1 = new double[p];
	m2 = new double[p * p];
	nom[0] = 0;
	for (int i = 0; i < p; i++)
	{
		int _d;
		fscanf(f, "%d", &_d);
		if (_d != 0)
			nom[i + 1] = _d + nom[i];
		else
		{
			p2 = i;
			p1 = p - nom[i];
			for (i++; i < p; i++)
				fscanf(f, "%d", &_d);
		}
	}
	for (int i = 0; i < p; i++)
	{
		fscanf(f, "%lf", c1m1 + i);
		if (c1m1[i] < 0)
			c1m1[i] = 0;
	}
	for (int i = 0; i < p; i++)
	{
		fscanf(f, "%lf", c2m1 + i);
		if (c2m1[i] < 0)
			c2m1[i] = 0;
	}
	for (int i = 0; i < p * p; i++)
		fscanf(f, "%lf", m2 + i);
	fclose(f);
}
void sample()
{
	FILE * f = fopen("cppout.txt", "w");
	int * bestdis = new int[p2];
	double * bestcont = new double[p1];
	double * cont = new double[p1];
	int * stack = new int[p2 + 1];
	double * m1, *em1;
	double * ec1m1 = new double[p];//emperical
	double * ec2m1 = new double[p];
	double * em2 = new double[p * p];
	double * tm1 = new double[p];//target
	double * tm2 = new double[p * p];
	double * sump2 = new double[p + 1];//sums in the path
	double sum = 0;
	stack[0] = 0;
	stack++;
	memset(ec1m1, 0, sizeof(double) * p);
	memset(ec2m1, 0, sizeof(double) * p);
	memset(em2, 0, sizeof(double) * p * p);
	int count1 = 0, count2 = 0;
	for (int count = 1; count < n + 1; count++)
	{
		if (count % 100 == 0)
			printf("Sample Count: %d\n", count);
		int cc;
		if (rand() < prob * 32768.0)
		{
			//class 1
			m1 = c1m1;
			em1 = ec1m1;
			count1++;
			cc = count1;
		}
		else
		{
			//class 2
			m1 = c2m1;
			em1 = ec2m1;
			count2++;
			cc = count2;
		}
		//compute targets
		for (int j = 0; j < p - p1; j++)
			tm1[j] = change(em1[j] - m1[j] * cc);
		for (int j = p - p1; j < p; j++)
			tm1[j] = em1[j] - m1[j] * cc;
		for (int j = 0; j < p * p; j++)
			tm2[j] = em2[j] - m2[j] * count;
		for (int j = 0; j < p - p1; j++)
			for (int k = 0; k < p - p1; k++)
				tm2[j * p + k] = change(tm2[j * p + k]);
		//printf("%lf\n", tm2[0]);
		//sampling
		double maxnum = 100000000.0;
/*		double * rank[3];
		rank[0] = new double [p * 3 + 3];
		rank[1] = rank[0] + p + 1;
		rank[2] = rank[1] + p + 1;
*/		double * rank = new double [p * 3 + 3];
		
		double * answer = new double [p + 3];
		int c = 1;
		stack[0] = 0;
		sum = 0;
		stack[1] = 0;
		while (c != 0)
		{
			//push in or go back
			if (c <= p2)
			{
				if (stack[c] < nom[c])
				{
					sump2[c] = tm2[stack[c] * p + stack[c]] + tm1[stack[c]];
					for (int i = 0; i < c; i++)
						sump2[c] += tm2[stack[c] * p + stack[i]];
					sum += sump2[c];
					c++;
					stack[c] = nom[c - 1];
				}
				else
				{
					c--;
					sum -= sump2[c];
					stack[c]++;
				}
			}
			else
			{
				for (int j = p - p1; j < p; j++)
				{
					int j0 = j - p + p1;
					int po = 0;
/*					//rank means ax+b as a penalty. the first comes from the tm2 with known continuous numbers
					for (int k = p - p1; k < j; k++)
					{
						int k0 = k - p + p1;
						rank[0][k0] = cont[k0];
						rank[1][k0] = tm2[j * p + k];
						if (rank[0][k0] != 0)
						{
							rank[2][k0] = -rank[1][k0] / rank[0][k0];
							if (rank[2][k0] < 1 && rank[2][k0] > 0)
								answer[po++] = rank[2][k0];
						}
					}
					//then those that are unknown. use their expectations
					for (int k = j; k < p; k++)
					{
						int k0 = k - p + p1;
						rank[0][k0] = m1[k];
						rank[1][k0] = tm2[j * p + k];
						if (rank[0][k0] != 0)
						{
							rank[2][k0] = -rank[1][k0] / rank[0][k0];
							if (rank[2][k0] < 1 && rank[2][k0] > 0)
								answer[po++] = rank[2][k0];
						}
					}
					//discrete ones
					for (int k = 0; k < p2; k++)
					{
						rank[1][k + p1] = tm2[j * p + stack[k]];
						rank[0][k] = 1;
						rank[2][k] = -rank[1][k];
						if (rank[2][k] < 1 && rank[2][k] > 0)
							answer[po++] = rank[2][k];
					}
					//mean
					rank[0][p1 + p2] = 1;
					rank[1][p1 + p2] = tm1[j];
					rank[2][p1 + p2] = -tm1[j];
					if (tm1[j] > -1 && tm1[j] < 0)
						answer[po++] = tm1[j];
					answer[po++] = 0;
					answer[po++] = 1;
					double _min = 100000;
					int _argmin = -1;
					for (int k = 0; k < p; k++)
					{
						double _value = 0;
						for (int l = 0; l < p1 + p2 + 1; l++)
							_value += abs(rank[0][l] * answer[k] + rank[1][l]);
						if (_value < _min)
						{
							_min = _value;
							_argmin = k;
						}
					}
*/
					//rank means ax+b as a penalty. the first comes from the tm2 with known continuous numbers
					//if (j == p - p1)
					//	printf("\n");
					for (int k = p - p1; k < j; k++)
					{
						int k0 = k - p + p1;
						rank[k0] = cont[k0];
						rank[p + 1 + k0] = tm2[j * p + k];
						if (rank[k0] != 0)
						{
							rank[2 * p + 2 + k0] = -rank[1 + p + k0] / rank[k0];
							if (rank[2 + p * 2 + k0] < 1 && rank[2 + p * 2 + k0] > 0)
								answer[po++] = rank[2 + p * 2 + k0];
						}
					}
					//then those that are unknown. use their expectations
					for (int k = j; k < p; k++)
					{
						int k0 = k - p + p1;
						rank[k0] = m1[k];
						rank[1 + p + k0] = tm2[j * p + k];
						if (rank[k0] != 0)
						{
							rank[2 + p * 2 + k0] = -rank[1 + p + k0] / rank[k0];
							if (rank[2 + p * 2 + k0] < 1 && rank[2 + p * 2 + k0] > 0)
								answer[po++] = rank[2 + p * 2 + k0];
						}
					}
					//discrete ones
					for (int k = 0; k < p2; k++)
					{
						rank[1 + p + k + p1] = tm2[j * p + stack[k]];
						rank[k + p1] = 1;
						rank[2 + p * 2 + k + p1] = -rank[1 + p + k + p1];
						if (rank[2 + p * 2 + k + p1] < 1 && rank[2 + p * 2 + k + p1] > 0)
							answer[po++] = rank[2 + p * 2 + k + p1];
					}
					//mean
					rank[p1 + p2] = 1;
					rank[1 + p + p1 + p2] = tm1[j];
					rank[2 + p * 2 + p1 + p2] = -tm1[j];
					if (tm1[j] > -1 && tm1[j] < 0)
						answer[po++] = -tm1[j];
					answer[po++] = 0;
					answer[po++] = 1;
					double _min = 100000;
					int _argmin = -1;
					for (int k = 0; k < po; k++)
					{
						//if (j == p-p1)
						//printf("%lf ", answer[k]);
						double _value = 0;
						for (int l = 0; l < p1 + p2 + 1; l++)
							_value += abs(rank[l] * answer[k] + rank[1 + p + l]);
						if (_value < _min)
						{
							_min = _value;
							_argmin = k;
						}
					}
					//if (j == p-p1)
					//printf("%d\n", _argmin);
					cont[j0] = answer[_argmin];
				}
				double _d = 0;
				for (int j = p - p1; j < p; j++)
				{
					double value = cont[j - p + p1];
					for (int k = 1; k < p2 + 1; k++)
						_d += change(tm2[j * p + stack[k]], value);
					for (int k = p - p1; k <= j; k++)
						_d += change(tm2[j * p + k], value * cont[k - p + p1]);
				}
				if (maxnum > _d + sum)
				{
					maxnum = _d + sum;
					memcpy(bestdis, stack + 1, sizeof(int) * p2);
					memcpy(bestcont, cont, sizeof(double) * p1);
				}
				c--;
				sum -= sump2[c];
				stack[c]++;
			}
		}
		//output
		for (int i = 0; i < p2; i++)
			fprintf(f, "%d ", bestdis[i]);
		for (int i = 0; i < p1; i++)
			fprintf(f, "%lf ", bestcont[i]);
		if (m1 == c1m1)
			fprintf(f, "1\n");
		else
			fprintf(f, "0\n");
		for (int i = 0; i < p2; i++)
		{
			em1[bestdis[i]] += 1;
			for (int j = 0; j < p2; j++)
				em2[bestdis[i] * p + bestdis[j]] += 1;
		}
		for (int i = 0; i < p1; i++)
			for (int j = 0; j < p2; j++)
			{
				em2[bestdis[j] * p + i + p - p1] += bestcont[i];
				em2[bestdis[j] + (i + p - p1) * p] += bestcont[i];
			}
		for (int i = 0; i < p1; i++)
		{
			int i0 = i + p - p1;
			em1[i0] += bestcont[i];
			for (int j = 0; j < p1; j++)
			{
				int j0 = j + p - p1;
				em2[i0 * p + j0] += bestcont[i] * bestcont[j];
			}
		}
	}
	fclose(f);
#ifdef DEBUG
	double _max = 0;
	for (int i = 0; i < p * p; i++)
	{
		if (abs(em2[i] - m2[i]) > _max)
			_max = abs(em2[i] - m2[i]);
	}
	printf("Max error: %lf\n", _max / n);
#endif
}

int main()
{
	//initialize all statistics and parameters. read from cppin.txt
	readstat();
	//sampling and write to cppout.txt
	sample();
	return 0;
}