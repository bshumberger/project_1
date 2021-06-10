//Molecule Header File
#include <string>

using namespace std;

class Molecule
{
public:
	int natom;
	int *zvals;
	double *x;
	double *y;
	double *z;

	void print_geom();
	double R(int i, int j);
	double ex(int i, int j);
	double ey(int i, int j);
	double ez(int i, int j);
	double excross(int i, int j, int k, int l);
	double eycross(int i, int j, int k, int l);
	double ezcross(int i, int j, int k, int l);
	double edot(int i, int j, int k, int l);
	double phi(int i, int j, int k);
	double theta(int i, int j, int k, int l);
	double tau(int i, int j, int k, int l);
	void translate(double a, double b, double c);

	Molecule(const char *filename);
	~Molecule();
};
