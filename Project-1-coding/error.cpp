/*
**     Project1: d) computing the relative error, when increasing n to n=10^7
**     The algorithm for solving the tridiagonal matrix
**     equation is implemented (O(8n) FLOPS).
*/
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include "time.h"
# include "lib.h"

using namespace std;
ofstream ofile;

// Declaring functions that will be used:
double Solution(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}
double f(double x) {return 100*exp(-10*x);}
double error(double u, double v) {return log10(fabs(v-u));} //calculating the error

void initialize (int *initial_n, int *number_of_steps, int *step_increase){
    printf("Read in initial_n, number_of_steps and step_increase \n");
    scanf("%d %d %d", initial_n, number_of_steps, step_increase);
    return;
}

double max_rel_error(int n) {
    double max_error;
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
       // Could print results to check:
       //cout << "x = " << x[i] << " and " << "h^2*f(x) = " << h*h*f(x[i]) << endl;

    // Filling up b_twiddle array, i.e. right hand side of equation:
    for (int i=1; i<=n; i++) {
        b_twidd[i] = h*h*f(x[i]);
        // Could print here to check:
        //cout << "b_twidd = " << b_twidd[i] << "for x = " << x[i] << endl;
        u[i] = Solution(x[i]);
        //cout << "u = " << u[i] << " for x = " << x[i] <<  endl;
        b[i] = 2;
        a[i] = -1;
        c[i] = -1;
    }
    c[n] = 0;
    a[1] = 0;
    // Algorithm for finding v:
    // a(i)*v(i-1) + b(i)*v(i) + c(i)*v(i+1) = b_twidd(i)
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

    //largest relative error:
    max_error = fabs(error(u[1], v[1]));
    for (int i=1; i<=n; i++) {
        if (fabs(error(u[i], v[i])) < fabs(max_error)) {
            max_error = error(u[i], v[i]); //update largest value
        }
    }

    //Freeing memory
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
    double max;
    char *outfilename;

    // Read in output file from command-line,
    // abort if there are too few command-line arguments:
    if( argc <= 1 ){
      cout << "Bad Usage: " << argv[0] <<
          " PLIS read in output file on same line" << endl;
      exit(1);
    }
    else{
      outfilename = argv[1]; // first command line argument
    }

    initialize (&initial_n, &number_of_steps, &step_increase);
    double *h_log10 = new double[number_of_steps];
    double *max_error_log10 = new double[number_of_steps];


    //Write to file
    ofile.open(outfilename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "       n:             log10(h):          log10(max_error):  " << endl;
    int n = initial_n;
    for (int i=0; i<number_of_steps; i++) {
       h_log10[i] = log10(1./(n+1));
       max_error_log10[i] = max_rel_error(n);
       ofile << setw(15) << setprecision(8) << n;
       ofile << setw(15) << setprecision(8) << h_log10[i];
       ofile << setw(15) << setprecision(8) << max_error_log10[i] << endl;
       n *= step_increase;
    }
    ofile.close();

    return 0;


}