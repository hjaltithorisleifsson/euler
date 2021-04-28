#ifndef EULER_UTILS_H_
#define EULER_UTILS_H_

#include <stdbool.h>
#include <stdlib.h>

/*
############################################################
# Various utility functions for mathematical problems. 
# Author: Hjalti Thor Isleifsson
############################################################

############################################################
#################### Number theory #########################
############################################################
*/

int gcd(int a, int b);


int lcm(int a, int b);


void bezout(int m, int n, int* a, int* b, int* c);


int inv(int x, int p);


bool isprime(int n);

//All primes <= n
void allprimes(int n, int** primes, int* n_primes);

//Computes the Euler's totient function of n.
int phi(int n);

/*
############################################################
######################## Sequences #########################
############################################################
*/

//Returns the nth Fibonacci number in ln(n) time.
//f_0 = 1, f_1 = 1 and f_(n+2) = f_(n+1) + f_n
unsigned long fibonacci(unsigned int n);

/**
 *
 * Gives a list of the first n fibonacci numbers, defined as 
 * f_0 = 1, f_1 = 1, and f_(n+2) = f_(n+1) + f_n.
 * 
 */
unsigned long* fiblist11(int n);

/**
 *
 * Gives a list of the first n fibonacci numbers, defined as 
 * f_0 = 1, f_1 = 2, and f_(n+2) = f_(n+1) + f_n.
 *
 */
unsigned long* fiblist12(int n);

/**
 *
 * Returns the greatest n such that f_n <= N.
 * Here f_n are 1, 2, 3, 5, 8,...
 *
 */
unsigned int fibinv(unsigned long N);

/*
#########################################################
################## Numerical schemes ####################
#########################################################
*/

double trapezoidal(double (*f)(double), double a, double b, int n);

double rombergintegrate(double (*f)(double), double a, double b, int n, double tol);

/*
#######################################################
################# Basic operations ####################
#######################################################
*/

unsigned int isqrt(unsigned int n);

//Computes x^n where n in ln(n) time where n is a positve integer.
double dipow(double x, int n);

/**
 *
 * Takes in a sorted array of integers. Returns the smallest idx such that 
 * array[idx] = tofind, if such an index exists. Else it returns -(idx + 1) 
 * where idx is the least number such that  tofind < array[idx], 
 * or len if all the elements are smaller than tofind.
 *
**/
int binsearch(int tofind, int* array, size_t len);


#endif