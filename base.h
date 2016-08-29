#ifndef BASE
#define BASE
#include "framework.h"

#define PRECISION 1e-5

class Point
{
public:
	double x, y, z;
	//everything just inline
	Point(){}
	Point(const char* key){ if (key == "zero")this->zero(); }
	Point(double a, double b, double c) :x(a), y(b), z(c){}
	//GPoint(const Sphere& p) :x(p.x), y(p.y), z(p.z){printf("aaa done once\n");}
	//GPoint(const GPoint& p) :x(p.x), y(p.y), z(p.z){}
	Point& operator=(const Point& p){ x = p.x; y = p.y; z = p.z; return(*this); }

	double coord_x()const { return(x); }
	double coord_y()const { return(y); }
	double coord_z()const { return(z); }
	Point operator+(const Point& p){ return Point(x + p.x, y + p.y, z + p.z); }
	Point operator-(const Point& p){ return Point(x - p.x, y - p.y, z - p.z); }
	Point operator*(const double& c){ return Point(x*c, y*c, z*c); }
	Point operator/(const double& c){ return Point(x/c, y/c, z/c); }
	Point& operator+=(const Point& p){ x += p.x; y += p.y; z += p.z; return *this; }
	Point& operator-=(const Point& p){ x -= p.x; y -= p.y; z -= p.z; return *this; }
	Point& operator*=(const double& c){ x *= c; y *= c; z *= c; return *this; }
	Point& operator/=(const double& c){ x /= c; y /= c; z /= c; return *this; }

	bool operator==(const Point& p){ return((abs(p.x - x)<1e-8 && abs(p.y - y)<1e-8 && abs(p.z - z)<1e-8)); }
	//bool close(GPoint p, double preci = 1.e-5){ return((abs(p.x - x)<preci && abs(p.y - y)<preci && abs(p.z - z)<preci)); }
	void print(FILE *fptr){ fprintf(fptr, "%lf %lf %lf \n", x, y, z); }

	void point_assign(double xx, double yy, double zz) { x = xx; y = yy; z = zz; }
	double dot(const Point& p)const { return(x*p.x+y*p.y+z*p.z); }

	void zero()	{x = 0;y = 0;z = 0;}
	double norm()const {return(double(sqrt(x*x + y*y + z*z)));}
};

// The real model which is a hard sphere
class Sphere
{
public:
	double x,y,z,r,k;
	Sphere(){}
	Sphere(double a, double b, double c, double d, double e) :x(a),y(b),z(c),r(d),k(e){}
	//Sphere(const GPoint& p) :GPoint(p.x, p.y, p.z),r(RADIUS){}

	Sphere operator+(const Point& p){ return Sphere(x + p.x, y + p.y, z + p.z, r, k); }
	Sphere operator-(const Point& p){ return Sphere(x - p.x, y - p.y, z - p.z, r, k); }
	bool operator==(const Sphere& p){ return((abs(p.x - x)<PRECISION && abs(p.y - y)<PRECISION && abs(p.z - z)<PRECISION && abs(p.r - r)<PRECISION && abs(p.k - k)<PRECISION)); }
	//Point operator-(Sphere p){ return Point(x - p.x, y - p.y, z - p.z); }
	Point center()const {return Point(x,y,z);}

	void assign(double xx, double yy, double zz, double rr,double kk) { x = xx; y = yy; z = zz; r = rr; k = kk; }
	void print(FILE *fptr){ fprintf(fptr, "%lf %lf %lf %lf %lf \n", x, y, z, r, k); }
	void scan(FILE *fptr){ fscanf(fptr, "%lf %lf %lf %lf %lf \n", &x, &y, &z, &r, &k); }
	double distance(const Sphere& s){return(sqrt((x-s.x)*(x-s.x) + (y-s.y)*(y-s.y) + (z-s.z)*(z-s.z)));}

	//int WellSeparate(const Sphere& s){ double t = sqrt((x-s.x)*(x-s.x) + (y-s.y)*(y-s.y) + (z-s.z)*(z-s.z)); return((t > s.r+this->r) ? int(t+(1-s.r-this->r)) : 0); }
	int WellSeparate(const Sphere& s){ double t = sqrt((x-s.x)*(x-s.x) + (y-s.y)*(y-s.y) + (z-s.z)*(z-s.z)); return((t>=1)?int(t):((t>=(s.r+this->r))?1:0)); }
	//double norm(){return(sqrt(x*x + y*y + z*z));}
	//inline void euclidean_op(Sphere* p, Point* ref, OpMatrix* op);{(*this) = op->dot(*p) + (*ref);}
};

// 3x3 matrix
class Matrix
{
public:
	Point row1;
	Point row2;
	Point row3;
	Matrix(){}
	Matrix(const char* key){ if (key == "identity")this->identity(); }
	Matrix(Point a, Point b, Point c) :row1(a), row2(b), row3(c){}
	Matrix(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3) :row1(x1, y1, z1), row2(x2, y2, z2), row3(x3, y3, z3){}
	//Matrix(RMatrix x) :row1(x.row1), row2(x.row2), row3(x.row3){}

	Point col1()const { return Point(row1.coord_x(), row2.coord_x(), row3.coord_x()); }
	Point col2()const { return Point(row1.coord_y(), row2.coord_y(), row3.coord_y()); }
	Point col3()const { return Point(row1.coord_z(), row2.coord_z(), row3.coord_z()); }
	Matrix& operator=(const Matrix& x){ row1 = x.row1;  row2 = x.row2; row3 = x.row3; return(*this); }
	//bool operator==(Matrix m){ return((row1 == m.row1 && row2 == m.row2 && row3 == m.row3)); }
	//bool close(Matrix m, double preci = 0.01){ return((row1.close(m.row1, preci) && row2.close(m.row2, preci) && row3.close(m.row3, preci))); }
	Matrix inv()
	{
		double determinant = row1.x*(row2.y*row3.z - row3.y*row2.z)
			- row1.y*(row2.x*row3.z - row2.z*row3.x)
			+ row1.z*(row2.x*row3.y - row2.y*row3.x);
		double invdet = 1 / determinant;
		return Matrix(
			(row2.y*row3.z - row3.y*row2.z)*invdet,
			-(row1.y*row3.z - row1.z*row3.y)*invdet,
			(row1.y*row2.z - row1.z*row2.y)*invdet,
			-(row2.x*row3.z - row2.z*row3.x)*invdet,
			(row1.x*row3.z - row1.z*row3.x)*invdet,
			-(row1.x*row2.z - row2.x*row1.z)*invdet,
			(row2.x*row3.y - row3.x*row2.y)*invdet,
			-(row1.x*row3.y - row3.x*row1.y)*invdet,
			(row1.x*row2.y - row2.x*row1.y)*invdet);
	}

	void print(FILE *fptr){ row1.print(fptr); row2.print(fptr); row3.print(fptr); }
	Matrix& identity(){ row1 = Point(1, 0, 0); row2 = Point(0, 1, 0); row3 = Point(0, 0, 1); return(*this); }

	Point dot(const Point& v){ return Point(row1.dot(v), row2.dot(v), row3.dot(v)); }		
	Matrix dot(const Matrix& x){ return Matrix(row1.dot(x.col1()), row1.dot(x.col2()), row1.dot(x.col3()), row2.dot(x.col1()), row2.dot(x.col2()), row2.dot(x.col3()), row3.dot(x.col1()), row3.dot(x.col2()), row3.dot(x.col3())); }
	Sphere dot(const Sphere& v){Point c=v.center();return Sphere(row1.dot(c),row2.dot(c),row3.dot(c),v.r,v.k);}
};

int Random_integer_uniform(int a,int b){return a+int((b-a)*RNG_NAME());}
int Random_integer_log(int a,int b)//not efficiency
{
	if(b==a)return 0;
	int t,i;
	for (i=0;i<1000;i++)
	{
		t = a+int((b-a)*RNG_NAME());
		if(RNG_NAME()<=log(1+1.0/(t-a+1))/log(b-a+1.0)) break;
	}
	if(i==999)printf("error(Random_integer_log): 1000 trial not success but still in domain.\n");
	return t;
}

Matrix Random_symmetry()
{
	double ta,phi;
	ta = M_PI * RNG_NAME();
	phi = 2* M_PI * RNG_NAME();
	double l =sin(ta)*cos(phi),m = sin(ta)*sin(phi), n = cos(ta);

	double theta = 2* M_PI * RNG_NAME();
	Matrix rot(l*l*(1-cos(theta))+cos(theta), m*l*(1-cos(theta))-n*sin(theta), n*l*(1-cos(theta))+m*sin(theta), 
		l*m*(1-cos(theta))+n*sin(theta), m*m*(1-cos(theta))+cos(theta), n*m*(1-cos(theta))-l*sin(theta),
		l*n*(1-cos(theta))-m*sin(theta), m*n*(1-cos(theta))+l*sin(theta), n*n*(1-cos(theta))+cos(theta));

	ta = M_PI * RNG_NAME();
	phi = 2* M_PI * RNG_NAME();
	double a =sin(ta)*cos(phi),b = sin(ta)*sin(phi), c = cos(ta);
	Matrix ref(1-2*a*a, -2*a*b, -2*a*c, -2*a*b, 1-2*b*b, -2*b*c, -2*a*c, -2*b*c, 1-2*c*c);
	return ref.dot(rot);
}

int twopowof(int power){int temp = 1; for(int i=0;i<power;i++) temp*=2; return temp;}
int factorial(int n){if(n==1|| n==0)return 1;else return n*factorial(n-1);}

Point origin("zero");
Matrix I("identity");

#endif