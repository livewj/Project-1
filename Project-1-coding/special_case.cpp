/*
**     Project1: a) and b)
**     The algorithm for solving the tridiagonal matrix
**     equation is implemented (O(8n) FLOPS).
**     Also finding the CPU time
*/
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>

using namespace std;
ofstream ofile;

// Declaring two functions that will be used:
double Solution(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}

double f(double x) {return 100*exp(-10*x);}

// Main program reads filename and n from command line:
int main(int argc, char* argv[]) {

    // Declaration of initial variables:
    char *outfilename;
    int n;

    // Read in output file and n,
    // abort if there are too few command-line arguments:
    if( argc <= 2 ){
      cout << "Bad Usage: " << argv[0] <<
          " read also output file and n (int) on same line" << endl;
      exit(1);
    }
    else{
      outfilename = argv[1]; // first command line argument.
      n = atoi(argv[2]);  // second command line argument.
    }

    // Constants of the problem:
    double h = 1.0/(n+1.0);
    double *x = new double[n+2];
    double *f_ = new double[n+1]; // construction with n+1 points to make
                                       // indexing close to mathematics.
    f_[0] = 0;

    // The constituents of the tridiagonal matrix A:
    // Zeroth element not needed, but included to make indexing easy:
    int *a = new int[n+1];
    int *b = new int[n+1];
    int *c = new int[n+1];

    // Temporal variabel in Gaussian elimination:
    double *diag_temp = new double[n+1];

    // Real solution and approximated one:
    double *u = new double[n+2]; // Analytical solution
    double *u_num= new double[n+2]; // Numerical solution
    // Including extra points to make the indexing easy:
    u[0] = 0;
    u_num[0] = 0;

    // Filling up x-array, making x[0] = 0 and x[n+1] = 1:
    for (int i=0; i<=n+1; i++) {
        x[i] = i*h;
        // Could print results to check:
        //cout << "x = " << x[i] << " and " << "h^2*f(x) = " << h*h*f(x[i]) << endl;
    }

    // Filling up b_twiddle array, i.e. right hand side of equation:
    for (int i=1; i<=n; i++) {
        f_[i] = h*h*f(x[i]);
        // Could print here to check:
        //cout << "f_ = " << f_[i] << "for x = " << x[i] << endl;
        u[i] = Solution(x[i]);
        //cout << "u = " << u[i] << " for x = " << x[i] <<  endl;
        b[i] = 2;
        a[i] = -1;
        c[i] = -1;
    }
    c[n] = 0;
    a[1] = 0;

   

    //special case
    double *b_tilde = new double [n+2];
    double *f_tilde = new double [n+2];
    b_tilde[1] = b[1]; //=2
    f_tilde[1] = f_[1];

    //precalculate updated diagonal elements
    b_tilde[0] = b_tilde[n] = 2;
    u_num[0] = u_num[n] = 0.0;
    for (int i = 1; i < n; i++){
        b_tilde[i] = (i+1.0)/((double) i);
    }

    double special_time;
    //Start timer
    clock_t start, finish;
    start = clock();
    
    //forward substitution
    for (int i=2; i<n; i++) {
        f_tilde[i] = f_[i] + f_tilde[i-1]/b_tilde[i-1];
    }

    //backward substitution
    u_num[n-1] = f_tilde[n-1]/b_tilde[n-1];
    for (int i = n-2; i > 0; i--) {
        u_num[i] = (f_tilde[i] + u_num[i+1])/b_tilde[i];
    }

    //stop timer
    finish = clock();
    special_time = ( (double) (finish - start)/CLOCKS_PER_SEC );



    // Open file and write results to file:
    ofile.open(outfilename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "       x:             u(x):          u_num(x):        time:" << endl;
    for (int i=1;i<=n;i++) {
       ofile << setw(15) << setprecision(8) << x[i];
       ofile << setw(15) << setprecision(8) << u[i];
       ofile << setw(15) << setprecision(8) << u_num[i];
       ofile << setw(15) << setprecision(8) << special_time << endl;
    }
    ofile.close();

    delete [] x;
    delete [] f_;
    delete [] a;
    delete [] b;
    delete [] c;
    delete [] u;
    delete [] u_num;
    delete [] f_tilde;
    delete [] b_tilde;

    cout << "CPU-time = " << special_time << "for n = " << n << endl;

    return 0;
}