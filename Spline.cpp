/*
Author: Justin Eugene
MATH 3180-D01
Project to create a spline function

Input: 
 None

Output: 
 print the data points of a function
 print a tridiagonal matrix before and after Forward Elimination
 print the solutions for the second derivative of S(x), z
 print a final evaluation using an array of x values
*/

#include <iostream>
#include <iomanip>

using namespace std;

//function prototypes
void datapoints(double *t, double* y, int size);
void trisys(double *u,double *h, double *volh, double *v, int size);
void createNewV(double *u,double *v,double *h, double *volh, double *y, double *yoh, int size);
void forwardelim(double *u,double *v, double*h, double *volh, int height);
void printdash(int amt);
void findz(double *u, double *h, double *volh, double *v, double *z, int size);
void findCoefficients(double *a, double *b, double *c, double *d, double *y, double *h, double *z, double *yoh, int size);
void printC(double *a, double *b, double *c, double *d, int size);
void printSpline(int size, double *a, double *b, double *c, double *d, double *t);
double spline(int size, double *a, double *b, double *c, double *d, double *t, double x);
void printE(int size, double *a, double *b, double *c, double *d, double *t, double *x);

int main() {
	double t[] = {-8,-6,-4, -2, 0, 2, 4,6,8};
	int size = sizeof(t)/sizeof(double);
	double y[size], z[size], h[size-1], u[size-2], v[size-2];

	//variables t, y, and h are now defined
	for (int i = 0; i < size; i++) {
		y[i] = 1/((t[i] * t[i]) + 1);

				if(i != size-1) {
						h[i] = t[i+1] - t[i];

				}
	}

	//print out the datapoints;
	datapoints(t, y, size);
	cout << endl;

	//define u and v next
	double yoh[size-1];
	double volh[size - 1];
	createNewV(u, v, h, volh, y, yoh, size);

	//print matrix before and after elimination
	cout << "Tri-diagonal system for z1 through z(n-1)" << endl;
	trisys(u, h, volh, v, size);
	//elimation
	forwardelim(u, v, h, volh, size-2);
	
	//printing after forward elimination
	cout << "Tri-diagonal system after Forward Elimination" << endl;
	trisys(u, h, volh, v, size);

	//printing the solutions
	cout << "Solutions for S''(ti) = zi" << endl;
	printdash(27);
	findz(u, h, volh, v, z, size);

	//create coefficient arrays
	double a[size-1], b[size-1], c[size-1], d[size-1];
	findCoefficients(a, b, c, d, y, h, z, yoh, size);
	//print coefficients
	printC(a,b, c, d, size);

	//print nested form
	printSpline(size, a, b, c, d, t);

	//final evaluation
    //an array for the x values we are testing is needed
	double x[] = {-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8};
	printE(size, a, b, c, d, t, x);
}

void datapoints(double *t, double* y, int size) {
	cout << "Data points (ti, yi)" << endl << endl;

	for (int i = 0; i < size; i++) {
		cout << "(t" << i << ",. y" << i << ") = " 
			<< "(" << t[i] << "," << y[i] << ')' << endl;
	}

}
void createNewV(double *u,double *v,double *h, double *volh, double *y, double *yoh, int size) {
	yoh[0] = (y[1] -y[0])/h[0];
	volh[0] = h[0];
	
	for (int i = 0; i < size-2; i++) {
		u[i] = 2*(h[i+1]+h[i]);
		yoh[i+1] = (y[i+2]-y[i+1])/h[i+1];
		v[i] = 6*(yoh[i+1]-yoh[i]);
		volh[i+1] = h[i+1]; 
	}
}
void trisys(double *u,double *h, double *volh, double *v, int size) {
    //height width of the matrix w/o augment
	int height = size-2;
    
	for(int i = 0; i < height; i++) {
		//print the leading 0s
		for(int p = height - i + 1; p < height; p++)
			cout << setw(10) << 0;
		
		if(i != 0)
			cout << setw(10) << volh[i];

		cout << setw(10) << u[i];

		if(i != height-1)
			cout << setw(10) << h[i+1];

		//now print the trailing 0's
		for(int p = height - (height-2) + i; p < height; p++)
			cout << setw(10) << 0;

		cout << setprecision(6) << setw(10) << v[i];
		cout << endl;
	}
	cout << endl;
}
//function used to help format the output
void printdash(int amt) {
	for(int i = 0; i < amt; i++)
		cout << '-';
	cout << endl;
}
//forward elimination function used in matrices
void forwardelim(double *u,double *v, double*h, double *volh, int height) {
	double m;
	for(int i = 0; i < height; i++) {
		//find multiplier
		m = volh[i+1]/u[i];

		//update row
		volh[i+1] -= m*u[i];
		u[i+1] -= m*h[i+1];
		v[i+1] -= m*v[i];

	}
}
void findz(double *u, double *h, double *volh, double *v, double *z, int size) {
	double m;
	
	z[0] = 0;
	z[size-1] = 0;
	
	for(int i = size-2; i > 0; i--) {
		m = 1/u[i-1];
		//substract
		if(i < size-2) {
			volh[i+1] = h[i+1]*z[i+1];
			v[i] -= volh[i]; 
		}
		//cross multiply
		z[i] = v[i-1]*m;
	}

	//print solutions
	for(int i = 0; i < size; i++) {
		cout << "z" << i << " : " << z[i] << endl;
	}
	
	cout << endl;
}
//function finds the coefficient needed to create the spline
void findCoefficients(double *a, double *b, double *c, double *d, double *y, double *h, double *z, double *yoh, int size) {
    //find number of coefficients
	int cnum = size - 1;
    
	for(int i = 0; i < cnum; i++) {
		a[i] = y[i];
		b[i] = (-h[i]/6)*z[i+1] - ((h[i]/3)*z[i]) + yoh[i];
		c[i] = z[i]/2;
		d[i] = (z[i+1] - z[i])/(6*h[i]);
	}
}

//function prints the coefficients of the Natural Cubic Spline
void printC(double *a, double *b, double *c, double *d, int size) {
	cout << "Coefficients for Natural Cubic Spline in nested form" << endl;
	printdash(67);
	
	for(int i = 0; i < size-1; i++) {
		cout << "A" << i << ":" << setw(11) << setprecision(6) << a[i];
		cout << setw(4) << "B" << i << ":" << setw(11) <<  setprecision(6) << b[i];
		cout << setw(4) << "C" << i << ":" << setw(11) <<  setprecision(6) << c[i];
		cout << setw(4) << "D" << i << ":" << setw(11) <<  setprecision(6) << d[i];
		cout << endl;
	}

	cout << endl;
}
//print spline function
void printSpline(int size, double *a, double *b, double *c, double *d, double *t) {
	cout << "Natural Cubic Spline in nested form" << endl;
	printdash(67);

	int snum = size - 1;
	char sym;

	for(int i = 0; i < snum; i++) {
		cout << 'S' << i << " :  " << right << a[i] << "+(x";

        //if t is positive use a '-' sym otherwise use '+' sign
		if(t[i] <= 0)
			sym = '+';
		else if(t[i] > 0)
			sym = '-';
		
		cout << sym << abs(t[i]) << ")*(" << b[i] << "+(x";
		cout << sym << abs(t[i]) << ")*(" << c[i] << "+(x";
		cout << sym << abs(t[i]) << ")*(" << d[i] << ")))";
		cout << endl;
	}

	cout << endl;
}

double spline(int size, double *a, double *b, double *c, double *d, double *t, double x) {
	double s = 0;
	double diff = 0;
	for(int i = 0; i < size-1; i++) {
		diff = x - t[i];
		if(x <= t[i+1]) {
			s = a[i] + diff*(b[i] + diff*(c[i] + diff*d[i]));
			return s;
		}
	}
	
	return s;
}
//print evaluation
void printE(int size, double *a, double *b, double *c, double *d, double *t, double *x) {
	cout << "Evaluation of Original and Spline functions and the absolute errors" << endl;
	printdash(67);
	cout << " i " << setw(9) << "xi" << setw(16) << "f(xi)" << setw(15) << "S(xi)" << setw(19) << "|f(xi) -S(xi)|" << endl;
	printdash(67);

	double func = 0;
	double s = 0;
	
	for(int i = 0; i < 17; i++) {

		func = 1/((x[i] * x[i]) + 1);
		s = spline(size, a, b, c, d, t, x[i]);
		
		cout << setw(2) << i << fixed << setw(12) 
			<< setprecision(6) << x[i] 
			<< setw(15) << func;

		cout << setw(15) << s;
		cout << setw(15) << func - s;

		cout << endl;
	}
}