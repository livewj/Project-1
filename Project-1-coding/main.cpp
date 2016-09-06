#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;
ofstream ofile; //for writing to file

//This function uses Gaussian elimination to solve a tridiagonal matrix.
//Our main problem is to solve the equation - u''(x) = f(x)

//the given function f(x):
double f(double x) {
    return 100*exp(-10*x);
}

//this is the ananytical solution u(x):
double analytical(double x) {
    return 1 - (1 - exp(-10))*x - exp(-10*x);
}


//this function gives us a numerical solution to the differential-equation
//by using Gaussian elimination.
//this program takes two arguments from the command-line: the filename and the size, n
int main(int argc, char* argv[]) {
    char *writetofile; //output filename
    int n;

    //command-line arguments
    writetofile = argv[1];
    n = atoi(argv[2]); //converting from ascii to int

    //declaring constants
    double h = 1./(n+1); //step size
    double *x = new double[n+2];
    double *b_tilde = new double[n+1];
    b_tilde[0] = 0;
    //vectors of the matrix A:
    int *a = new int[n+1];
    int *b = new int[n+1];
    int *c = new int[n+1];
    //the new elements on the diagonal
    //double *diagonal_new = new double[n+1];
    //solutions
    double *u = new double[n+2]; //analytical
    u[0] = 0;
    double *v = new double[n+2]; //numerical
    v[0] = 0;
    double *e_tilde = new double[n+2];
    e_tilde[0] = 0;
    double *d_tilde = new double[n+2];
    d_tilde[0] = 0;

    //finding the values of x
    for (int i=0; i<=n+1; i++) {
        x[i] = i*h;
    }

    //finding the values of b_tilde and u.
    //filling up the vectors a, b and c
    for (int i=1; i<=n; i++) {
        b_tilde[i] = h*h*f(x[i]);
        u[i] = analytical(x[i]);
        a[i] = -1;
        b[i] = 2;
        c[i] = -1;

    }
    a[1] = 0;
    c[n] = 0;

    //finding the numerical solution using Gauss elimination
    //forward substitution:
    //initial conditions
    d_tilde[1] = b[1];

    e_tilde[1] = b_tilde[1]/d_tilde[1];
    for (int i=2; i<=n; i++) {
        //diagonal_new[i] = c[i-1]/d_tilde;

        d_tilde[i] = b[i] - (c[i-1]/d_tilde[i-1])*a[i];
        e_tilde[i] =(b_tilde[i] - e_tilde[i-1]*a[i])/d_tilde[i];
    }

    //backward substitution:
    for (int i=n-1; i>=1; i--) {
        v[i] = e_tilde[i] - (c[i]/d_tilde[i])*v[i+1];

    }

    //finally, writing the results to a file:
    ofile.open(writetofile);
    ofile << "x:            u(x):           v(x):" << endl;
    for (int i=1; i <=n; i++) { //writing to file
        ofile << setw(10) << x[i]; //x-values
        ofile << setw(10) << u[i]; //analytical solution
        ofile << setw(10) << v[i] << endl;  //numerical solution

    }
    ofile.close();

    //free memory
    delete [] x;
    delete [] u;
    delete [] v;
    delete [] a;
    delete [] b;
    delete [] c;
    delete [] b_tilde;

    return 0;
}

















