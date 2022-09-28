#include <iostream>
#include <fstream>
#include <string>
#include <iomanip> // для setw()
#include "Problems.h"

/*
В данной программе реализовано численное решение задачи Дирихле для уравнения теплопроводности в прямоугольнике. 

Краевые условия:
При x = xMin u(t, x, y) = mu1(t, y)
При x = xMax u(t, x, y) = mu2(t, y)
При y = yMin u(t, x, y) = mu3(t, x)
При y = yMax u(t, x, y) = mu4(t, x)

Начальное условие: u(0, x, y) = phi(x, y). 

Выбрана явная схема метода конечных разностей. 

Решение задачи на данном временном шаге представлено в виде одномерного вектора u.
*/

// Структура с параметрами дискретизации
struct Grid {
	double xMin, xMax, yMin, yMax; // Границы прямоугольника
	int Nx, Ny;                    // Число узлов сетки по оси X и оси Y соответственно
	int N;                         // Общее число узлов
	double hx, hy;                 // Шаг по оси X и оси Y соответственно
	double hx_sq, hy_sq;           // Шаги по оси X и оси Y, возведённые в квадрат
	double tau;                    // Шаг по времени
	double T;                      // Максимальное значение параметра t
};

// Указатели на функции реализованы как глобальные переменные, чтобы не писать огромную сигнатуру
// для подпрограмм, которые используют эти функции
double (*uExact)(double, double, double);
double (*f)(double, double, double);
double (*phi)(double, double);
double (*mu_xMin)(double, double, double); // Значение решения на границе x = xMin
double (*mu_xMax)(double, double, double);
double (*mu_yMin)(double, double, double);
double (*mu_yMax)(double, double, double);

void assignUOnBorders(double* u, Grid grid, double t);                 // Присвоить вектору u значения краевых функций
void initU(double* u, Grid grid);                                      // Инициализировать решение u
void gridInfoToFile(std::string filename, Grid grid);                  // Записать в заголовок файла информацию о сетке grid
void writeArrayToStream(double* ar, int N, std::ostream& streamOut);   // Дописать в поток вывода ostream массив в виде строки
void getExactSolution(double* ar, Grid grid, double t);                // Записать в массив ar значения точного решения
void getDiff(double* diffAr, double* a1, double* a2, int N);           // Записать в массив diffAr модуль разности векторов a1 и a2
void copyArray(double* dest, double* source, int N);                   // Скопировать массив source в массив dest

int main() {
	// Условия задачи
	// Вариант 1
	uExact = uExact1;
	f = f1;
	phi = phi1;
	mu_xMin = mu1;
	mu_xMax = mu1;
	mu_yMin = mu1;
	mu_yMax = mu1;
	Grid grid1; 
	grid1.xMin = 0;
	grid1.xMax = M_PI;
	grid1.yMin = 0;
	grid1.yMax = M_PI;

	// Вариант 2
	/*uExact = uExact2;
	f = f2;
	phi = phi2;
	mu_xMin = mu_xMin2;
	mu_xMax = mu_xMax2;
	mu_yMin = mu_yMin2;
	mu_yMax = mu_yMax2;
	Grid grid1;
	grid1.xMin = 3;
	grid1.xMax = 6;
	grid1.yMin = -2;
	grid1.yMax = 2;*/

	// Инициализация параметров сетки
	grid1.Nx = 20;
	grid1.Ny = 30;
	grid1.N = grid1.Nx * grid1.Ny;

	grid1.hx = (grid1.xMax - grid1.xMin) / (grid1.Nx - 1.0);
	grid1.hy = (grid1.yMax - grid1.yMin) / (grid1.Ny - 1.0);
	grid1.hx_sq = grid1.hx * grid1.hx;
	grid1.hy_sq = grid1.hy * grid1.hy;

	grid1.tau = 2e-3;
	grid1.T = 3;
	double tauThreshold = grid1.hx_sq * grid1.hy_sq / (grid1.hx_sq + grid1.hy_sq) / 2.0;
	if (grid1.tau > tauThreshold) {
		printf("Time step (tau = %6.4f) is too large!\n", grid1.tau);
		printf("For given hx and hy time step must be <= %6.2e, or the solution will be unstable.\n", tauThreshold);
		return 0;
	}

	// Вывод в файл шапку с необходимой информацией
	std::string filename = "Results.txt";
	gridInfoToFile(filename, grid1);

	// Инициализация u, u1 и diffArray в момент времени 0
	double* uPast = new double[grid1.N];      // Массив для приближенного решения
	double* u1 = new double[grid1.N];         // Массив для точного решения
	double* diffArray = new double[grid1.N];  // Массив для разницы точного и приближенного решений
	initU(uPast, grid1);
	getExactSolution(u1, grid1, 0);
	getDiff(diffArray, uPast, u1, grid1.N);

	// Вывод массивов в файл
	std::ofstream outstream(filename, std::ios::app);
	outstream << std::setw(14) << std::setprecision(4) << 0.0000;
	writeArrayToStream(uPast, grid1.N, outstream);
	outstream << std::setw(14) << "";
	writeArrayToStream(u1, grid1.N, outstream);
	outstream << std::setw(14) << "";
	writeArrayToStream(diffArray, grid1.N, outstream);
	outstream << "\n";

	double t = 0;
	double x, y;
	double Cx = grid1.tau / grid1.hx_sq;
	double Cy = grid1.tau / grid1.hy_sq;
	int s;
	double tau2 = 1.0 / 30.0;  // Период сохранения решения в файл
	double t2 = 0.0;         
	double* uNew = new double[grid1.N];
	while (t < grid1.T) {
		t += grid1.tau;
		printf("t = %8.4f of %8.4f\r", t, grid1.T);
		assignUOnBorders(uNew, grid1, t);

		//Вычисление решения на новом шаге
		for (int i = 1; i < grid1.Ny - 1; i++) {
			for (int j = 1; j < grid1.Nx - 1; j++) {
				x = grid1.xMin + grid1.hx * j;
				y = grid1.yMin + grid1.hy * i;
				s = j + i * grid1.Nx;

				uNew[s] = uPast[s] + Cx * (uPast[s + 1] - 2 * uPast[s] + uPast[s - 1]) + Cy * (uPast[s + grid1.Nx] - 2*uPast[s] + uPast[s - grid1.Nx]) + grid1.tau*f(t, x, y);
			}
		}

		if (t > t2) {
			t2 += tau2;

			//Вычисление точного решения и разницы между точным и приближённым решениями
			getExactSolution(u1, grid1, t);
			getDiff(diffArray, uNew, u1, grid1.N);

			//Вывод в файл
			outstream << std::setw(14) << std::setprecision(4) << t;
			writeArrayToStream(uNew, grid1.N, outstream);
			outstream << std::setw(14) << "";
			writeArrayToStream(u1, grid1.N, outstream);
			outstream << std::setw(14) << "";
			writeArrayToStream(diffArray, grid1.N, outstream);
			outstream << "\n";
		}
		copyArray(uPast, uNew, grid1.N);
	}
	printf("t = %8.4f of %8.4f\n", t, grid1.T);

	// Удаление сущностей
	delete[] uPast;
	delete[] uNew;
	delete[] u1;
	delete[] diffArray;

	outstream.close();
	return 0;
}


void assignUOnBorders(double* u, Grid grid, double t) {
	double y;
	for (int i = 0; i < grid.Ny; i++) {
		y = grid.yMin + i * grid.hy;
		u[i * grid.Nx] = mu_xMin(t, grid.xMin, y);                 // Решение вдоль границы x = xMin
		u[grid.Nx - 1 + i * grid.Nx] = mu_xMax(t, grid.xMax, y);   // Решение вдоль границы x = xMax
	}
	double x;
	for (int j = 0; j < grid.Nx; j++) {
		x = grid.xMin + j * grid.hx;
		u[j] = mu_yMin(t, x, grid.yMin);                           // Решение вдоль границы y = yMin
		u[j + (grid.Ny - 1) * grid.Nx] = mu_yMax(t, x, grid.yMax); // Решение вдоль границы y = yMax
	}
}

void initU(double* u, Grid grid) {
	assignUOnBorders(u, grid, 0);

	double x, y;
	int s;
	for (int i = 1; i < grid.Ny - 1; i++) {
		for (int j = 1; j < grid.Nx - 1; j++) {
			x = grid.xMin + grid.hx * j;
			y = grid.yMin + grid.hy * i;
			s = j + i * grid.Nx;
			u[s] = phi(x, y);
		}
	}
}

void gridInfoToFile(std::string filename, Grid grid) {
	std::ofstream outstream(filename);

	outstream << "// xMin = " << grid.xMin << "\n";
	outstream << "// xMax = " << grid.xMax << "\n";
	outstream << "// yMin = " << grid.yMin << "\n";
	outstream << "// yMax = " << grid.yMax << "\n";
	outstream << "//\n";
	outstream << "// Nx = " << grid.Nx << "\n";
	outstream << "// Ny = " << grid.Ny << "\n";
	outstream << "// hx = " << grid.hx << "\n";
	outstream << "// hy = " << grid.hy << "\n";
	outstream << "//\n";
	outstream << "// tau = " << grid.tau << "\n";
	outstream << "// T = " << grid.T << "\n";
	outstream << "//\n";
	outstream << "// x = xMin + hx*j, y = yMin + hy*i, s = j + i*Nx\n";
	outstream << "//\n";
	outstream << "//           t" << std::setw(13) << "" << "uApproximate, uExact, difference\n";

	outstream.close();
}

void writeArrayToStream(double* ar, int N, std::ostream& streamOut)
{
	for (int s = 0; s < N; s++)
		streamOut << std::setw(14) << std::setprecision(4) << ar[s];
	streamOut << "\n";
}

void getExactSolution(double* ar, Grid grid, double t) {
	double x, y;
	int s;
	for (int i = 0; i < grid.Ny; i++) {
		for (int j = 0; j < grid.Nx; j++) {
			x = grid.xMin + grid.hx * j;
			y = grid.yMin + grid.hy * i;
			s = j + i * grid.Nx;
			ar[s] = uExact(t, x, y);
		}
	}
}

void getDiff(double* diffAr, double* a1, double* a2, int N) {
	for (int s = 0; s < N; s++)
		diffAr[s] = abs(a1[s] - a2[s]);
}

void copyArray(double* dest, double* source, int N) {
	for (int s = 0; s < N; s++)
		dest[s] = source[s];
}