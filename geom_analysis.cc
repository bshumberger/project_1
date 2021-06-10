//Project # 1
//Header.....................................................................
#include <iostream>
#include <fstream>
#include <cstdio>
#include "molecule_class.h"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

using namespace std;

double mass[] = {
        0.000000,
        1.007825,
        4.002603,
        7.016003,
	9.012183,
        11.009305,
        12.000000,
        14.003074,
        15.994914,
        18.998403,
        19.992440,
        22.989769,
        23.985041,
        26.981538
        };


int main()
{
Molecule input("geom.dat");

cout << "Number of Atoms: " << input.natom << endl;
cout << "Input Coordinates:" << endl;
input.print_geom();

cout << "Interatomic Distances (bohr):" << endl;
for (int i = 0; i < input.natom; i++)
	for (int j = 0; j < i; j++)
		printf("%d %d %8.5f\n", i, j, input.R(i, j));

//cout << "Unit Vectors:" << endl;
//for (int i = 0; i < input.natom; i++)
//	for (int j = 0; j < i; j++)
//		printf("%d %d %8.5f %8.5f %8.5f\n", i, j, input.ex(i, j), input.ey(i, j), input.ez(i, j));

cout << "Bond Angles:" << endl;
for (int i = 0; i < input.natom; i++)
	for (int j = 0; j < i; j++)
		for (int k = 0; k < j; k++)
			if (input.R(i, j) < 4.0 && input.R(j, k) < 4.0)
				printf("%d-%d-%d %10.6f\n", i, j, k, input.phi(i, j, k)*(180/acos(-1)));

cout << "Out-of-Plane Angles:" << endl;
for (int i = 0; i < input.natom; i++)
        for (int k = 0; k < input.natom; k++)
                for (int j = 0; j < input.natom; j++)
			for (int l = 0; l < j; l++)
				if (i!=j && i!=k && i!=l && j!=k && j!=l && k!=l && input.R(k, j) < 4.0 && input.R(k, l) < 4.0 && input.R(k, i) < 4.0)
					 printf("%d-%d-%d-%d %10.6f\n", i, j, k, l, input.theta(i, j, k, l)*(180/acos(-1)));

cout << "Dihedral Angles:" << endl;
for (int i = 0; i < input.natom; i++)
        for (int j = 0; j < i; j++)
                for (int k = 0; k < j; k++)
                        for (int l = 0; l < k; l++)
                                if (i!=j && i!=k && i!=l && j!=k && j!=l && k!=l && input.R(i, j) < 4.0 && input.R(j, k) < 4.0 && input.R(k, l) < 4.0)
                                         printf("%d-%d-%d-%d %10.6f\n", i, j, k, l, input.tau(i, j, k, l)*(180/acos(-1)));

//Center of Mass..................................................................................


double M = 0;
double xM = 0;
double yM = 0;
double zM = 0;

for (int i = 0; i < input.natom; i++){
        M += mass[input.zvals[i]];
        xM += input.x[i]*mass[input.zvals[i]];
	yM += input.y[i]*mass[input.zvals[i]];
	zM += input.z[i]*mass[input.zvals[i]];
}

cout << "Center of Mass:" << endl;
printf("%12.8f %12.8f %12.8f\n", xM/M, yM/M, zM/M);


//Principle Moments of Inertia.....................................................................

input.translate(-xM/M, -yM/M, -zM/M);

Matrix I(3,3);

for (int i = 0; i < input.natom; i++){
	I(0,0) += mass[input.zvals[i]]*(input.y[i]*input.y[i]+input.z[i]*input.z[i]);
	I(1,1) += mass[input.zvals[i]]*(input.x[i]*input.x[i]+input.z[i]*input.z[i]);
	I(2,2) += mass[input.zvals[i]]*(input.x[i]*input.x[i]+input.y[i]*input.y[i]);
	I(0,1) -= mass[input.zvals[i]]*input.x[i]*input.y[i];
	I(0,2) -= mass[input.zvals[i]]*input.x[i]*input.z[i];
	I(1,2) -= mass[input.zvals[i]]*input.y[i]*input.z[i];
}

I(1,0) = I(0,1);
I(2,0) = I(0,2);
I(2,1) = I(1,2);

cout << "Moment of Inertia Tensor (amu * bohr^2):" << endl;
cout << I << endl;

Eigen::SelfAdjointEigenSolver<Matrix> solver(I);
Matrix evecs = solver.eigenvectors();
Matrix evals = solver.eigenvalues();

cout << "Principle Moments of Inertia (amu * bohr^2):" << endl;
cout << evals << endl;

double conv = 0.529177249 * 0.529177249;
cout << "\nPrincipal moments of inertia (amu * AA^2):\n";
cout << evals * conv << endl;
 
conv = 1.6605402E-24 * 0.529177249E-8 * 0.529177249E-8;
cout << "\nPrincipal moments of inertia (g * cm^2):\n";
cout << evals * conv << endl;

if (evals(0) == evals(1) == evals(2))
	cout << "Molecule is a spherical rotor." << endl;
else if (evals(0) == evals(1) < evals(2))
	cout << "Molecule is an oblate symmetric rotor." << endl;
else if (evals(0) < evals(1) == evals(2))
        cout << "Molecule is a prolate symmetric rotor." << endl;
else if (evals(0) != evals(1) != evals(2))
	cout << "Molecule is an asymmetric rotor." << endl;


//Rotational Constants.............................................................................


double _pi = acos(-1.0);
conv = 6.6260755E-34/(8.0 * _pi * _pi);
conv /= 1.6605402E-27 * 0.529177249E-10 * 0.529177249E-10;
conv *= 1e-6;
cout << "\nRotational constants (MHz):\n";
cout << "\tA = " << conv/evals(0) << "\t B = " << conv/evals(1) << "\t C = " << conv/evals(2) << endl;

conv = 6.6260755E-34/(8.0 * _pi * _pi);
conv /= 1.6605402E-27 * 0.529177249E-10 * 0.529177249E-10;
conv /= 2.99792458E10;
cout << "\nRotational constants (cm-1):\n";
cout << "\tA = " << conv/evals(0) << "\t B = " << conv/evals(1) << "\t C = " << conv/evals(2) << endl;


return 0;
}
