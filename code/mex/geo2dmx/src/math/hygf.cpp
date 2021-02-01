/********************************************************************
 *
 * hygf.cpp
 * 
 * Routine to compute the hypergeometric function f(a,b,c,x)
 * for real arguments a = 1, b, c = 1 + b, x
 * with b >= 0 and x <= 1
 *
 * Translation of a FORTRAN 77 code from
 * "Computation of Special Functions"
 * by Shanjie Zhang and Jianming Jin
 * Copyright 1996 by John Wiley & Sons, Inc.
 * 
 * January 2003, Christoph Schwarzbach,
 * cschwarz@student.tu-freiberg.de
 *
 ********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double hygf(double b, double x);
double psi(double x);

double hygf(double b, double x)
/* 
 * Compute the hypergeometric function for a real argument
 * f(a,b,c,x), with a = 1, b >= 0,  c = 1 + b, -1 <= x <= 1.
 */
{
  double eps, c0, c1, p, r, w, y;
  int k;
  
  c0 = fabs(x);
  if (c0 > 1.0 || b < 0)
  {
    fprintf(stderr,
      "Invalid arguments for hygf.\n");
    exit(0);
  }

  if (c0  > 0.95)
  {
    eps = 1.0e-8;
  }
  else
  {
    eps = 1.0e-15;
  }
  
  if (fabs(1.0 - x) < eps)
  {
//    fprintf(stderr,"The hypergeometric series is divergent.\n");
    return(HUGE_VAL);
  }
  if (x == 0.0 || b == 0.0)
  {
    return(1.0);
  }

  if (x < 0.0)
  {
    p = x / (x - 1.0);
    y = 1.0;
    w = 0.0;
    r = 1.0;
    if (b < 1.0)
    {
      c0 = b;
      for (k = 1; k <= 500; k++)
      {
        c1 = b + (double) k;
        r *= c0 * c0 / (c1 * (double) k)  * p;
        y += r;
        if(fabs(y-w) < eps * fabs(y)) break;
        w = y;
        c0 = c1;
      }
      y /= pow(1.0 - x, b);
    }
    else
    {
      for (k = 1; k <= 500; k++)
      {
        r *= ((double) k) / (b + (double) k) * p;
        y += r;
        if(fabs(y-w) < eps * fabs(y)) break;
        w = y;
      }
      y /= (1.0 - x);
    }
  }
  else
  {
    if (x < 0.9)
    {
      y = 1.0;
      r = 1.0;
      w = 0.0;
      c0 = b;
      for (k = 1; k <= 1500; k++)
      {
        c1 = b + (double) k;
        r *= c0 / c1 * x;
        y += r;
        if(fabs(y - w) < eps * fabs(y)) break;
        w = y;
        c0 = c1;
      }
    }
    else
    {
      p = 0.0;
      r = 1.0;
      c0 = b;
      c1 = 0.5772156649015329 + psi(b) + log(1.0 - x);
      y = c1;
      w = 0.0;
      for (k = 1; k <= 500; k++)
      {
        p += (1.0 - b) / (c0 * (double) k);
        r *= c0 / ((double) k) * (1.0 - x);
        y += r * (c1 + p);
        if(fabs(y - w) < eps * fabs(y)) break;
        w = y;
        c0 = b + (double) k;
      }
      y *= -b;
    }
  }
  return(y);
}

double psi(double x)
/* Compute Digamma function Psi(x)*/
{
  const double pi = 3.141592653589793,
    el = 0.5772156649015329;
  double xa = fabs(x),
    s = 0.0,
    psi, x2;
  int i, n;
  
  if (x == floor(x) && x <= 0.0)
  {
    psi = HUGE_VAL;
  }
  else
  {
    if (xa == floor(xa))
    {
      n = (int) xa;
      for (i = 1; i <= n-1; i++)
      {
        s += 1.0 / (double) i;
      }
      psi = -el + s;
    }
    else
    {
      if (xa + 0.5 == floor(xa + 0.5))
      {
        n = (int) (xa - 0.5);
        for (i = 1; i <= n; i++)
        {
          s += 1.0 / (double) (2 * i - 1);
        }
        psi = -el + 2.0 * s - 1.386294361119891;
      }
      else
      {
        if (xa < 10.0)
        {
          n = 10 - (int) xa;
          for (i = 0; i <= n-1; i++)
          {
            s += 1.0 / (xa + (double) i);
          }
          xa += (double) n;
        }
        x2 = 1.0 / (xa * xa);
        psi = log(xa) - 0.5 / xa + 
          (((((((0.4432598039215686 * x2
          - 0.83333333333333333e-01) * x2
          + 0.21092796092796093e-01) * x2
          - 0.75757575757575758e-02) * x2
          + 0.41666666666666667e-02) * x2
          - 0.39682539682539683e-02) * x2
          + 0.83333333333333333e-02) * x2
          - 0.8333333333333e-01) * x2;
        psi -= s;
      }
    }
    if (x < 0.0)
    {
      psi -= pi / tan(pi*x) + 1.0 / x;
    }
  }
  return(psi);
}

