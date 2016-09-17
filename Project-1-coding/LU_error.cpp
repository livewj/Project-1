/*
**     Project1: e) Comparing the relative error for the LU-decomposition and the
**     algorithm for solving the tridiagonal matrix
**     
*/
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <cstring>
# include <time.h>
# include <string>
# include <cstdio>
# include <cstdlib>
# include <new>
# include "lib.cpp"

using namespace std;
ofstream ofile;

//Declaring functions
double Solution(double);
double f(double);
double LU_error(int);
double tridiagonal_error(int);
double error(double, double);


// Functions that will be used:
double Solution(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}
double f(double x) {return 100*exp(-10*x);}
double error(double u, double v) {return log10(fabs(v-u));} //calculating the error

void initialize (int *initial_n, int *number_of_steps, int *step_increase){
    printf("Read in initial_n, number_of_steps and step_increase \n");
    scanf("%d %d %d", initial_n, number_of_steps, step_increase);
    return;
}

double LU_error(int n) {
    double *x = new double[n+2];
    double **A; //our matrix
    double max_error;
    double h = 1.0/(n+1.0);
    double *b_twidd = new double[n];
    double *u = new double [n+2]; //analytical solution
    double *v = new double [n+2]; //numerical solution
    u[0] = 0;
    v[0] = 0;
    double CPU;

    //Filling up x-array
    for (int i=0; i<=n+1; i++) {
    x[i] = i*h; }

    //Right hand side of equation and analytical solution
    for (int i=1; i<=n; i++) {
        b_twidd[i] = h*h*f(x[i]); 
        u[i] = Solution(x[i]);
    }


    //Constructing the matrix
    A = new double*[n]; //nxn -matrix
    for (int i=0; i<n; i++) {
        A[i] = new double[n];
    }
    int flag;
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            flag = 0;
            if (i==j) {  //diagonal
                A[i][j] = 2;
                flag = 1;
            }
            if (i==j+1 || j == i+1) { //upper or lower diagonal
                A[i][j] = -1;
                flag = 1;
            }
            if (flag==0) {
                A[i][j] = 0;
            }
        }
    //cout << A << endl; 
    }


    //LU-decomposition of matrix A
    int *index;
    double d;
    double *b_twidd_backup = new double [n];
    index = new int[n]; //keeps track of the number of interchanges of rows.

    //backup of b_twidd
    for (int i=1; i<=n; i++) {
        b_twidd_backup[i] = h*h*f(x[i]);
    }

    //start timer
    clock_t start, finish;
    start = clock();

    ludcmp(A,n,index,&d); //LU-decomposition
    lubksb(A,n,index,b_twidd); //solve. The solution is now in b_twidd

    //stop timer
    finish = clock();
    CPU = ((double) (finish - start)/CLOCKS_PER_SEC );

    //print solution in terminal
    //for (int i=0; i<n; i++) {
    //    cout << "LU solution:" << b_twidd[i] << endl;}

    //largest relative error:
    max_error = fabs(error(u[1], b_twidd[1]));
    for (int i=1; i<=n; i++) {
        if (fabs(error(u[i], b_twidd[i])) < fabs(max_error)) {
            max_error = error(u[i], b_twidd[i]); //update largest value
        }
    }
    

    //free memory
    delete [] x;
    delete [] b_twidd;
    delete [] b_twidd_backup;
    for (int i=0; i<n; i++) { //delete matrix elements
        delete [] A[i];
    }
    delete [] A; //delete matrix

    return max_error;
}

double tridiagonal_error(int n) {
    double max_error;
    double CPU;
    double h = 1.0/(n+1.0);
    double *x = new double[n+2];
    double *b_twidd = new double[n+1]; // construction with n+1 points to make
                                       // indexing close to mathematics.
    b_twidd[0] = 0;

    // The constituents of the tridiagonal matrix A:
    // Zeroth element not needed, but included to make indexing easy:
    int *a = new int[n+1];
    int *b = new int[n+1];
    int *c = new int[n+1];

    // Temporal variabel in Gaussian elimination:
    double *diag_temp = new double[n+1];

    // Real solution and approximated one:
    double *u = new double[n+2]; // Analytical solution
    double *v = new double[n+2]; // Numerical solution
    // Including extra points to make the indexing easy:
    u[0] = 0;
    v[0] = 0;

    // Filling up x-array, making x[0] = 0 and x[n+1] = 1:
    for (int i=0; i<=n+1; i++) {
        x[i] = i*h; }

    // Filling up b_twiddle array, i.e. right hand side of equation:
    for (int i=1; i<=n; i++) {
        b_twidd[i] = h*h*f(x[i]);
        u[i] = Solution(x[i]);
        b[i] = 2;
        a[i] = -1;
        c[i] = -1;
    }
    c[n] = 0;
    a[1] = 0;


    // Algorithm for finding v:
    // a(i)*v(i-1) + b(i)*v(i) + c(i)*v(i+1) = b_twidd(i)

    //Start timer
    clock_t start, finish;
    start = clock();

    // Row reduction; forward substitution:
    double b_temp = b[1];
    v[1] = b_twidd[1]/b_temp;
    for (int i=2;i<=n;i++) {
        // Temporary value needed also in next loop:
        diag_temp[i] = c[i-1]/b_temp;
        // Temporary diagonal element:
        b_temp = b[i] - a[i]*diag_temp[i];
        // Updating right hand side of matrix equation:
        v[i] = (b_twidd[i]-v[i-1]*a[i])/b_temp;
    }

    // Row reduction; backward substition:
    for (int i=n-1;i>=1;i--) {
        v[i] -= diag_temp[i+1]*v[i+1];
    }

    //stop timer
    finish = clock();
    CPU = ( (double) (finish - start)/CLOCKS_PER_SEC );

    //print result to terminal
    //for (int i=1; i<n+1; i++) {
    //    cout << "Tridiagonal solution:" << v[i] << endl;}
  
    //largest relative error:
    max_error = fabs(error(u[1], v[1]));
    for (int i=1; i<=n; i++) {
        if (fabs(error(u[i], v[i])) < fabs(max_error)) {
            max_error = error(u[i], v[i]); //update largest value
        }
    }

    //Free memory
    delete [] x;
    delete [] b_twidd;
    delete [] a;
    delete [] b;
    delete [] c;
    delete [] u;
    delete [] v;

return max_error;

}


// Main program reads filename and n from command line:
int main(int argc, char* argv[]) {

    // Declaration of initial variables:
    int initial_n, number_of_steps, step_increase;
    double max_error;
    string outfilename;
    outfilename = "compare_error.txt";

    initialize (&initial_n, &number_of_steps, &step_increase);

    //Write to file
    ofile.open(outfilename);
    ofile << setiosflags(     ios::showpoint | ios::uppercase);
    ofile << "            n:        LU (error):      Tridiagonal (error): " << endl;
    int n = initial_n;
    for (int i=0; i<number_of_steps; i++) {
       ofile << setw(16) << setprecision(8) << n;
       ofile << setw(16) << setprecision(8) << LU_error(n);
       ofile << setw(16) << setprecision(8) << tridiagonal_error(n) << endl;
       n *= step_increase;
    }
    ofile.close();

    return 0;


}

