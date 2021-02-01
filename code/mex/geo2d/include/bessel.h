/*
 * bessel.h
 */
#ifndef BESSEL_H
#define BESSEL_H

/* Set following makro to 
   BESSELF to get fast polynomial evaluation of Bessel functions and 
           accuracy of 1.0E-008
   BESSELD to get slower evaluation of Bessel functions but with 
           higher accuracy in the order of EPS (1.0E-016)
 */
#define BESSELF

#ifdef BESSELF
double bessk0(double x);
double bessk1(double x);
double bessi0(double x);
double bessi1(double x);
#endif

#ifdef BESSELD
double bessk0(double x);
double bessk1(double x);
void bessik(double x, double xnu, double *ri, double *rk, double *rip, double *rkp);
void beschb(double x, double *gam1, double *gam2, double *gampl, double *gammi);
double chebev(double a, double b, double c[], int m, double x);
#endif

#endif

