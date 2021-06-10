//Molecule Functions File
#include "molecule_class.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
using namespace std;

void Molecule::print_geom()
{
for (int i = 0; i < natom; i++)
	printf("%d %8.5f %8.5f %8.5f\n", zvals[i], x[i], y[i], z[i]);
}


//Unit Vectors Function.................................................................................


double Molecule::ex(int i, int j)
{
return -(x[i]-x[j])/R(i, j);
}
double Molecule::ey(int i, int j)
{
return -(y[i]-y[j])/R(i, j);
}
double Molecule::ez(int i, int j)
{
return -(z[i]-z[j])/R(i, j);
}


//Dot Product Function.................................................................................


double Molecule::edot(int i, int j, int k, int l)
{
return (ex(i, j)*ex(k, l)+ey(i, j)*ey(k, l)+ez(i, j)*ez(k, l));
}


//Cross Product Function................................................................................


double Molecule::excross(int i, int j, int k, int l)
{
return ey(i, j)*ez(k, l)-ez(i, j)*ey(k, l);
}
double Molecule::eycross(int i, int j, int k, int l)
{
return -(ex(i, j)*ez(k, l)-ez(i, j)*ex(k, l));
}
double Molecule::ezcross(int i, int j, int k, int l)
{
return ex(i, j)*ey(k, l)-ey(i, j)*ex(k, l);
}


//Bond Angles Function..................................................................................


double Molecule::phi(int i, int j, int k)
{
return acos( edot(j, i, j, k) );
}


//Out-of-Plane Angles Function..........................................................................


double Molecule::theta(int i, int j, int k, int l)
{
return asin( (((excross(k, j, k, l)/sin(phi(j, k, l)))*ex(k, i))+((eycross(k, j, k, l)/sin(phi(j, k, l)))*ey(k, i))+((ezcross(k, j, k, l)/sin(phi(j, k, l)))*ez(k, i))) );
}


//Dihedral Angles Function..............................................................................


double Molecule::tau(int i, int j, int k, int l)
{
return acos( (excross(i, j, j, k)*excross(j, k, k, l)+eycross(i, j, j, k)*eycross(j, k, k, l)+ezcross(i, j, j, k)*ezcross(j, k, k, l))/(sin(phi(i, j, k))*sin(phi(j, k, l))) );
}


//Translation Function..................................................................................


void Molecule::translate(double a, double b, double c)
{
for (int i = 0; i < natom; i++){
	x[i] += a;
	y[i] += b;
	z[i] += c;
}
}


//Bond Length Function..................................................................................


double Molecule::R(int i, int j)
{
return sqrt( (x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]));
}


//Constructor...........................................................................................


Molecule::Molecule(const char *filename)
{
//Open File..............................................
ifstream input(filename);
assert(input.good());

//Allocate Memory........................................
input >> natom;
zvals = new int [natom];
x = new double [natom];
y = new double [natom];
z = new double [natom];

//Read Values to Memory Space............................
for (int i = 0; i < natom; i++)
	input >> zvals[i] >> x[i] >> y[i] >> z[i];

input.close();
}


//Destructor............................................................................................


Molecule::~Molecule()
{
delete[] zvals;
delete[] x;
delete[] y;
delete[] z;
}

