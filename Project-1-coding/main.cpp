#include <iostream>

using namespace std;
//This function uses Gaussian elimination to solve a tridiagonal matrix.
//Our main problem is to solve the equation - u''(x) = f(x)

double * vectors (unsigned int size,
                          double *a, //lower diagonal is now a pointer
                          double *b, //diagonal
                          double *c, //upper diagonal
                          double *b_tilde, //the solution, f(x)*h^2
                          bool special = 0) //the function can use the general algorithm
                                            //or the algorithm developed for our
                                            // special case
{
    double *v = new double[size]; //the solution goes here
    //General case
    if (!special){
        for (unsigned int i=1, i<si)
    }
}
                  

int main()
{

}

