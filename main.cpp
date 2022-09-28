#include <iostream>
#include <fstream>
#include <string>
#include <iomanip> // ��� setw()
#include "Problems.h"

/*
� ������ ��������� ����������� ��������� ������� ������ ������� ��� ��������� ���������������� � ��������������. 

������� �������:
��� x = xMin u(t, x, y) = mu1(t, y)
��� x = xMax u(t, x, y) = mu2(t, y)
��� y = yMin u(t, x, y) = mu3(t, x)
��� y = yMax u(t, x, y) = mu4(t, x)

��������� �������: u(0, x, y) = phi(x, y). 

������� ����� ����� ������ �������� ���������. 

������� ������ �� ������ ��������� ���� ������������ � ���� ����������� ������� u.
*/

// ��������� � ����������� �������������
struct Grid {
	double xMin, xMax, yMin, yMax; // ������� ��������������
	int Nx, Ny;                    // ����� ����� ����� �� ��� X � ��� Y ��������������
	int N;                         // ����� ����� �����
	double hx, hy;                 // ��� �� ��� X � ��� Y ��������������
	double hx_sq, hy_sq;           // ���� �� ��� X � ��� Y, ���������� � �������
	double tau;                    // ��� �� �������
	double T;                      // ������������ �������� ��������� t
};

// ��������� �� ������� ����������� ��� ���������� ����������, ����� �� ������ �������� ���������
// ��� �����������, ������� ���������� ��� �������
double (*uExact)(double, double, double);
double (*f)(double, double, double);
double (*phi)(double, double);
double (*mu_xMin)(double, double, double); // �������� ������� �� ������� x = xMin
double (*mu_xMax)(double, double, double);
double (*mu_yMin)(double, double, double);
double (*mu_yMax)(double, double, double);

void assignUOnBorders(double* u, Grid grid, double t);                 // ��������� ������� u �������� ������� �������
void initU(double* u, Grid grid);                                      // ���������������� ������� u
void gridInfoToFile(std::string filename, Grid grid);                  // �������� � ��������� ����� ���������� � ����� grid
void writeArrayToStream(double* ar, int N, std::ostream& streamOut);   // �������� � ����� ������ ostream ������ � ���� ������
void getExactSolution(double* ar, Grid grid, double t);                // �������� � ������ ar �������� ������� �������
void getDiff(double* diffAr, double* a1, double* a2, int N);           // �������� � ������ diffAr ������ �������� �������� a1 � a2
void copyArray(double* dest, double* source, int N);                   // ����������� ������ source � ������ dest

int main() {
	// ������� ������
	// ������� 1
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

	// ������� 2
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

	// ������������� ���������� �����
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

	// ����� � ���� ����� � ����������� �����������
	std::string filename = "Results.txt";
	gridInfoToFile(filename, grid1);

	// ������������� u, u1 � diffArray � ������ ������� 0
	double* uPast = new double[grid1.N];      // ������ ��� ������������� �������
	double* u1 = new double[grid1.N];         // ������ ��� ������� �������
	double* diffArray = new double[grid1.N];  // ������ ��� ������� ������� � ������������� �������
	initU(uPast, grid1);
	getExactSolution(u1, grid1, 0);
	getDiff(diffArray, uPast, u1, grid1.N);

	// ����� �������� � ����
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
	double tau2 = 1.0 / 30.0;  // ������ ���������� ������� � ����
	double t2 = 0.0;         
	double* uNew = new double[grid1.N];
	while (t < grid1.T) {
		t += grid1.tau;
		printf("t = %8.4f of %8.4f\r", t, grid1.T);
		assignUOnBorders(uNew, grid1, t);

		//���������� ������� �� ����� ����
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

			//���������� ������� ������� � ������� ����� ������ � ����������� ���������
			getExactSolution(u1, grid1, t);
			getDiff(diffArray, uNew, u1, grid1.N);

			//����� � ����
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

	// �������� ���������
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
		u[i * grid.Nx] = mu_xMin(t, grid.xMin, y);                 // ������� ����� ������� x = xMin
		u[grid.Nx - 1 + i * grid.Nx] = mu_xMax(t, grid.xMax, y);   // ������� ����� ������� x = xMax
	}
	double x;
	for (int j = 0; j < grid.Nx; j++) {
		x = grid.xMin + j * grid.hx;
		u[j] = mu_yMin(t, x, grid.yMin);                           // ������� ����� ������� y = yMin
		u[j + (grid.Ny - 1) * grid.Nx] = mu_yMax(t, x, grid.yMax); // ������� ����� ������� y = yMax
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