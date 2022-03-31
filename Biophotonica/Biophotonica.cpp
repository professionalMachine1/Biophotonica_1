#include <iostream>
#include <vector>
#include <fstream>
#include <tuple>
#include <time.h>
#include <string>

struct SystemP
{
	double sigma = 1, a = 0, normalfact = 1;
    int N = 0;
};

SystemP SplitOpt;
const double PI = 2 * asin(1);

//const std::string setting = "set terminal pngcairo size 1080,720 enhanced font 'Times New Roman,12'\n";
void PLOT(const std::string& commands)
{
	FILE* gpipe = _popen("gnuplot.exe -persist", "w");
	fprintf(gpipe, commands.c_str());
	_pclose(gpipe);
}

void PointDistribution(const std::vector<double>& x, const std::vector<double>& y);
void HistogrammImage(const std::vector<double>& R);
double NormalDistribution(bool RorNot);
double _2DNormalDistributionPolar();
double AdaptiveQuadrature(double a_, double b_, double eps);
double SimpsonForm(double a_, double b_);
double fValue(double x);
std::vector<double>& shellasort(std::vector<double>& v);

#define PRINT(name) std::cout << #name << " = " << name << std::endl

int funcNum = 1;

int main()
{
	system("chcp 1251");
	system("cls");

	srand(time(0));

	std::cout << "Введите sigma\n";
	std::cin >> SplitOpt.sigma;
	std::cout << "Введите число испытаний\n";
	std::cin >> SplitOpt.N;
	
	/*
		
		f = exp(-r^2/(2*sigma^2))

	*/

	funcNum = 1;
	SplitOpt.normalfact = 1 / AdaptiveQuadrature(0, 1000, pow(10, -8));
	std::vector<double> r(SplitOpt.N), x(SplitOpt.N), y(SplitOpt.N);
	double fi;
	for (int i = 0; i < SplitOpt.N; i++)
	{
		r[i] = NormalDistribution(true);
		fi = rand(), fi *= 2 * PI / RAND_MAX;
		x[i] = r[i] * cos(fi), y[i] = r[i] * sin(fi);
	}
	shellasort(r);
	PointDistribution(x, y);
	HistogrammImage(r);
	
	/*
	
		f = r*exp(-r^2/(2*sigma^2))
	
	*/
	
	funcNum = 2;
	SplitOpt.normalfact = 1;
	SplitOpt.normalfact = 1 / (AdaptiveQuadrature(pow(10, -8), 1000, pow(10, -8)));
	for (int i = 0; i < SplitOpt.N; i++)
	{
		r[i] = _2DNormalDistributionPolar();
		fi = rand(), fi *= 2 * PI / RAND_MAX;
		x[i] = r[i] * cos(fi), y[i] = r[i] * sin(fi);
	}
	shellasort(r);
	PointDistribution(x, y);
	HistogrammImage(r);

	/*
		
		x = exp(-x^2/(2*sigma^2)), y= exp(-y^2/(2*sigma^2))
		
	*/

	funcNum = 2;
	SplitOpt.normalfact = 1;
	SplitOpt.normalfact = 1 / (AdaptiveQuadrature(pow(10, -8), 1000, pow(10, -8)));
	for (int i = 0; i < SplitOpt.N; i++)
	{
		x[i] = NormalDistribution(false);
		y[i] = NormalDistribution(false);
		r[i] = sqrt(pow(x[i], 2) + pow(y[i], 2));
	}
	shellasort(r);
	PointDistribution(x, y);
	HistogrammImage(r);

	return 0;
}

/*

	Вывод на экран

*/

void PointDistribution(const std::vector<double>& x, const std::vector<double>& y)
{
	std::ofstream points("points.txt");

	for (int i = 0; i < SplitOpt.N; i++)
		points << x[i] << " " << y[i] << std::endl;

	std::string s;
	s = "set terminal wxt\n";
	s += "set grid\n";
	s += "plot 'points.txt' with points lc 'black' title'XY'\n";
	PLOT(s);
}

void HistogrammImage(const std::vector<double>& R)
{
	std::ofstream RandomPoints("Rand.txt"), fReal("RealValues.txt");
	double spanlength = 0.5, span, alignment, height = 0;
	span = spanlength, alignment = spanlength / 2.0;

	for (int i = 0; i < SplitOpt.N; )
	{
		while ((i < SplitOpt.N) && (R[i] <= span))
		{
			height += 1;
			i++;
		}
		RandomPoints << span - alignment << " " << height / (SplitOpt.N * spanlength) << std::endl;
		span += spanlength;
		height = 0;
	}

	double r = floor(R[0]), rmax = ceil(R[R.size() - 1]), h = pow(10, -2);
	while (r <= rmax)
	{
		fReal << r << " " << fValue(r) << std::endl;
		r += h;
	}

	std::string s;
	s = "set terminal wxt\n";
	s += "set grid\n";
	s += "plot 'Rand.txt' with boxes lc 'black' title'Rand', ";
	s += "'RealValues.txt' with lines lc 'red' title 'DistributionDensity'\n";
	PLOT(s);

}

/*

	Генерируем случайные значения

*/

double NormalDistribution(bool RorNot)
{
	double S, u1, u2, v1, v2;

	do
	{
		do
		{
			u1 = rand(), u1 /= RAND_MAX;
			u2 = rand(), u2 /= RAND_MAX;
			v1 = 2 * u1 - 1, v2 = 2 * u2 - 1;
			S = pow(v1, 2) + pow(v2, 2);
		} while (S >= 1);
		v1 = SplitOpt.a + SplitOpt.sigma * v1 * sqrt(-2 * log(S) / S);
		v2 = SplitOpt.a + SplitOpt.sigma * v2 * sqrt(-2 * log(S) / S);

	} while ((v1 < 0) && (RorNot));

	return v1;
}

double _2DNormalDistributionPolar()
{
	double u;
	u = rand(), u /= RAND_MAX;
	return sqrt(log(pow(1 - u, -2 * pow(SplitOpt.sigma, 2))));
}

/*

	Интегрирование

*/

double AdaptiveQuadrature(double a_, double b_, double eps)
{
	double I, I1, I2;

	I = SimpsonForm(a_, b_);
	I1 = SimpsonForm(a_, (a_ + b_) / 2), I2 = SimpsonForm((a_ + b_) / 2, b_);

	if (fabs(I1 + I2 - I) > eps)
		I = AdaptiveQuadrature(a_, (a_ + b_) / 2, eps / 2) + AdaptiveQuadrature((a_ + b_) / 2, b_, eps / 2);

	return I;
}

double SimpsonForm(double a_, double b_)
{
	return (b_ - a_) * (fValue(a_) + 4 * fValue((a_ + b_) / 2) + fValue(b_)) / 6;
}

double fValue(double x)
{
	double res = 0;
	switch (funcNum)
	{
	case 1:
		res = SplitOpt.normalfact * exp(-pow((x - SplitOpt.a) / SplitOpt.sigma, 2) / 2);
		break;
	case 2:
		res = SplitOpt.normalfact * x * exp(-pow((x - SplitOpt.a) / SplitOpt.sigma, 2) / 2);
		break;
	default:
		break;
	}
	return res;
}

/*

	Сортировка Шелла

*/

std::vector<double>& shellasort(std::vector<double>& v)
{
	double temp = 0;
	for (int d = v.size() / 2; d >= 1; d /= 2)
		for (int i = 0; i < v.size() - d; i++)
			for (int j = i; j >= 0; j--)
				if (v[j] > v[j + d])
				{
					temp = v[j];
					v[j] = v[j + d];
					v[j + d] = temp;
				}
	return v;
}