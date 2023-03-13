#include <string>
#include <random>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include "otrpo_func.cpp"


using namespace std;

random_device rd;
default_random_engine e(rd());

class MyDistribution
{
    public:
        MyDistribution(string desc, double form)
        {
            _description = desc;
            _v = form;
            updateK();
            updateA();
            updateP();
        }

        void setDescription(string desc)
        {
            _description = desc;
        }

        string getDescription()
        {
            return _description;
        }

        void setForm(double form)
        {
            _v = form;
            updateK();
            updateA();
            updateP();
        }

        double getForm()
        {
            return _v;
        }

        double variance()
        {
            auto pow_7 = pow(_v, 7);
            return 2 / _k * (igamma(3./7., pow_7) / 7 + exp(-pow_7) * (2 / pow(_a, 3) 
                + 2 * _v / (_a * _a) + (_v * _v) / _a));
        }

        double density(double x)
        {
            auto f = 0.;
            if (abs(x) > _v)
            {
                f = exp(-_a * (abs(x) - _v) - pow(_v, 7));
            }
            else
            {
                f = exp(-pow(abs(x), 7));
            }
            return f / _k;
        }

        double excess()
        {
            auto pow_7 = pow(_v, 7);
            return 2 / pow(variance(), 2) / _k * (igamma(5./7., pow_7) / 7 
                + exp(-pow_7) * (24 / pow(_a, 5) + 24 * _v / pow(_a, 4) 
                + 12 * _v * _v / pow(_a, 3) + 4 * pow(_v, 3) / _a / _a + pow(_v, 4) / _a));
        }

        double get_pseudorandom_number()
        {
            discrete_distribution<> disc_dist({_p, (1 - _p) / 2, (1 - _p) / 2});
            auto z = disc_dist(rd);
            double x1 = _v + 1;
            if (!z)
            {
                while (abs(x1) > _v)
                {
                    gamma_distribution<double> g_dist(1./7., 1);
                    bernoulli_distribution b_dist(0.5);

                    auto b_rand = b_dist(rd);
                    auto g_rand = g_dist(rd);
                    auto s = b_rand ? 1 : -1;

                    x1 = pow(g_rand, 1./7.) * s;
                }
                return x1;
            }
            else
            {
                exponential_distribution<double> e_dist(1 / _a);
                auto x2 = e_dist(rd);
                return z == 1 ? _v + x2 : -_v - x2;
            }
        }

    private:
        string _description;
        double _v;
        double _k;
        double _a;
        double _p;

        void updateK()
        {
            _k = 2. / 7. * (exp(-pow(_v, 7)) / pow(_v, 6) + igamma(1. / 7., pow(_v, 7)));
        }

        void updateA()
        {
            _a = 7 * pow(_v, 6);
        }

        void updateP()
        {
            _p = 2. / (7 * _k) * igamma(1./7., pow(_v, 7));
        }
};

// Вычисление среднего арифметического значения выборки (выборочное мат.ожидание)
double calcMean(vector<double> sample)
{
	return accumulate(begin(sample), end(sample), 0.) / sample.size();
}

// Вычисление медианы, минимального, максимального значений выборки
void calcMedianMinMax(vector<double> sample, double& median, double& min, double& max)
{
	auto pos = begin(sample) + sample.size() / 2; // [N/2 + 1]-порядковаястатистика
	nth_element(begin(sample), pos, end(sample));

	// Для нечетного размера выборки
	if (sample.size() % 2)
	{
		median = *pos;
		min = *min_element(begin(sample), pos); // Ищем в первой половине выборки
	}
	// Для четного размера выборки
	else
	{
		auto minmax = minmax_element(begin(sample), pos); // Ищем в первой половине выборки
		median = (*pos + *minmax.second) / 2.;// second - [N/2]-порядковая статистика
		min = *minmax.first; // first - минимальный элемент выборки
	}
	max = *max_element(pos, end(sample)); // Ищем во второй половине выборки
}

// Вычисление выборочной дисперсии распределения
double calcVariance(vector<double> sample, double mean)
{
	vector<double> x2(sample.size());
	for (int i = 0; i < x2.size(); i++)
		x2[i] = pow(sample[i] - mean, 2);
	return calcMean(x2);
}

// Вычисление коэффициента эксцесса распределения
double calcExcess(vector<double> sample, double mean)
{
	vector<double> x4(sample.size());
	for (int i = 0; i < sample.size(); i++)
		x4[i] = pow(sample[i] - mean, 4);
	double mu4 = calcMean(x4);
	return mu4 / pow(calcVariance(sample, mean), 2);
}

// Вычисление коэффициента асимметрии распределения
double calcAsym(vector<double> sample, double mean)
{
	vector<double> x3(sample.size());
	for (int i = 0; i < sample.size(); i++)
		x3[i] = pow(sample[i] - mean, 3);
	double mu3 = calcMean(x3);
	return mu3 / pow(calcVariance(sample, mean), 3 / 2);
}

vector<double> calcDensityOfSamplePoints(vector<double> sample, MyDistribution distr)
{
    vector<double> density(sample.size());
    for (auto i = 0; i < sample.size(); i++)
    {
        density[i] = distr.density(sample[i]);
    }
    return density;
}

int inputSampleSizeFromConsole()
{
    int N = 0;
    
    cout << "Enter size of the sample (>0): ";
    cin >> N;

    while (cin.fail() || N <= 0)
    {
        cout << "Type integer value more than zero" << endl;
        cin.clear();
		cin.ignore(10000, '\n');
        cout << "Enter size of the sample (>0): ";
        cin >> N;
    }
    
    return N;
}

double inputFormParamFromConsole()
{
    double param = 0.;

    cout << "Enter v (>0): ";
    cin >> param;

    while (cin.fail() || param <= 0.)
    {
        cout << "Type numeric value more than zero" << endl;
        cin.clear();
		cin.ignore(10000, '\n');
        cout << "Enter v (>0): ";
        cin >> param;
    }

    return param;
}

vector<double> generateSample(MyDistribution distr)
{
    auto N = inputSampleSizeFromConsole();
    vector<double> sample(N);
    for (auto i = 0; i < N; i++)
    {
        sample[i] = distr.get_pseudorandom_number();
    }
    return sample;
}

void writeVectorToFile(string fname, vector<double> vectorToWrite)
{
    ofstream out(fname);
    for (auto i = 0; i < vectorToWrite.size(); i++)
        out << vectorToWrite[i] << endl;
    out.close();
}

int main()
{
    auto v = inputFormParamFromConsole();
    MyDistribution md("My destribution", v);

    auto sample = generateSample(md);
    auto density = calcDensityOfSamplePoints(sample, md);

    writeVectorToFile("sample.txt", sample);
    writeVectorToFile("density.txt", density);

    double median = 0., min = 0., max = 0.;
    calcMedianMinMax(sample, median, min, max);
    auto mean = calcMean(sample);
    auto variance = calcVariance(sample, mean);
    auto excess = calcExcess(sample, mean);
    auto asym = calcAsym(sample, mean);
    cout << "Median: " << median << endl
         << "Minimum value: " << min << endl
         << "Maximum value: " << max << endl
         << "Mean: " << mean << endl
         << "Variance: " << variance << endl
         << "Theoretical variance: " << md.variance() << endl
         << "Excess: " << excess << endl
         << "Theoretical excess: " << md.excess() << endl
         << "Asymmetry: " << asym << endl;
    return 0;
}