/*
 * bessel.cpp
 */
#include "bessel.h"


#ifdef BESSELF

#include <math.h>

double bessi0(double xx)
{
  double          fax, y, ans;
  if (fabs(xx) < 3.75)
  {
    y = (xx / 3.75) * (xx / 3.75);
    ans = 1.0 + y * (3.5156229 + y *
		     (3.0899424 + y *
		      (1.2067492 + y *
		       (0.2659732 + y *
			(0.360768e-1 + y * 0.45813e-2)))));
  } else
  {
    fax = fabs(xx);
    y = 3.75 / fax;
    ans = (exp(fax) / sqrt(fax)) * (0.39894228 + y * (0.1328592e-1 + y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2 + y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1 + y * 0.392377e-2))))))));
  }
  return ans;
}

double bessi1(double x)
{
  double          ax, y, ans;

  if (fabs(x) < 3.75)
  {
    y = (x / 3.75) * (x / 3.75);
    ans = x * (0.5 + y *
	       (0.87890594 + y *
		(0.51498869 + y *
		 (0.15084934 + y *
		  (0.2658733e-1 + y *
		   (0.301532e-2 + y * 0.32411e-3))))));
  } else
  {
    ax = fabs(x);
    y = 3.75 / ax;
    ans = 0.2282967e-1 + y *
      (-0.2895312e-1 + y *
       (0.1787654e-1 - y * 0.420059e-2));
    ans = 0.39894228 + y *
      (-0.3988024e-1 + y *
       (-0.362018e-2 + y *
	(0.163801e-2 + y *
	 (-0.1031555e-1 + y * ans))));
    ans = (exp(ax) / sqrt(ax)) * ans;
  }
  return ans;
}

double bessk0(double xx)
{
  double          y, ans;
  if (xx <= 2.0)
  {
    y = xx * xx / 4.0;
    ans = (-log(xx / 2.0) * bessi0(xx)) +
      (-0.57721566 + y *
       (0.42278420 + y *
	(0.23069756 + y *
	 (0.3488590e-1 + y *
	  (0.262698e-2 + y *
	   (0.10750e-3 + y * 0.74e-5))))));
  } else
  {
    y = (2.0 / xx);
    ans = (exp(-xx) / sqrt(xx)) *
      (1.25331414 + y *
       (-0.7832358e-1 + y *
	(0.2189568e-1 + y *
	 (-0.1062446e-1 + y *
	  (0.587872e-2 + y *
	   (-0.251540e-2 + y * 0.53208e-3))))));
  }
  return ans;
}

double bessk1(double x)
{
  double          y, ans;

  if (x <= 2.0)
  {
    y = x * x / 4.0;
    ans = (log(x / 2.0) * bessi1(x)) + (1.0 / x) *
      (1.0 + y *
       (0.15443144 + y *
	(-0.67278579 + y *
	 (-0.18156897 + y *
	  (-0.1919402e-1 + y *
	   (-0.110404e-2 + y * (-0.4686e-4)))))));
  } else
  {
    y = 2.0 / x;
    ans = (exp(-x) / sqrt(x)) *
      (1.25331414 + y *
       (0.23498619 + y *
	(-0.3655620e-1 + y *
	 (0.1504268e-1 + y *
	  (-0.780353e-2 + y *
	   (0.325614e-2 + y * (-0.68245e-3)))))));
  }
  return ans;
}
#endif

#ifdef BESSELD

#include <iostream.h>
#include <math.h>
#define EPS 1.0e-16
#define FPMIN 1.0e-30
#define MAXIT 10000
#define XMIN 2.0
#ifndef M_PI
#define M_PI 3.1415926535897931
#endif

double bessk0(double x)
{
  double ri, rk, rip, rkp;
  bessik(x,0.0,&ri,&rk,&rip,&rkp);
  return(rk);
}

double bessk1(double x)
{
  double ri, rk, rip, rkp;
  bessik(x,1.0,&ri,&rk,&rip,&rkp);
  return(rk);
}

/*
 * For the following routines see
 * Press, W. H., Teukolsky, S. A., Vetterling, W. T. and Flannery, P. B.:
 * Numerical Recipes in C: The Art of Scientific Computing.
 * 
 */

void bessik(double x, double xnu, double *ri, double *rk, double *rip, double *rkp)
/*
 * Returns the modified Bessel functions ri = I_xnu, rk = K_xnu and their derivatives
 * rip = I'_xnu, rkp = K'_xnu, for positive x and for xnu >= 0. The relative accuracy
 * is within one or two significant digits of EPS. FPMIN is a number close to the
 * machine's smallest floating-point number. All internal arithmetic is in double
 * precision. To convert the entire routine to double precision, change the float
 * declarations above to double and decrease EPS to 10-16. Also convert the function
 * beschb. 
 */
{
	double a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,gam1,gam2,
		gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1,ripl,
		ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2;
	int i,l,nl;

	if (x <= 0.0 || xnu < 0.0) cerr << "bad arguments in bessik" << endl;
	nl=(int)(xnu+0.5);	// nl is the number of downward
	xmu=xnu-nl;					//	recurrences of the I's and upward
	xmu2=xmu*xmu;				//	recurrences of K's. xmu lies
	xi=1.0/x;						//	between -1/2 and 1/2.
	xi2=2.0*xi;
	h=xnu*xi;									// Evaluate CF1 by modified Lentz's
	if (h < FPMIN) h=FPMIN;		//	method (§5.2). 
	b=xi2*xnu;
	d=0.0;
	c=h;
	for (i=1;i<=MAXIT;i++) {
		b += xi2;
		d=1.0/(b+d);		// Denominators cannot be zero here,
		c=b+1.0/c;			//	so no need for special precautions.
		del=c*d;
		h=del*h;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > MAXIT) cerr << "x too large in bessik; try asymptotic expansion" << endl;
	ril=FPMIN;				// Initialize I_xnu and I'_xnu for downward
	ripl=h*ril;				//	recurrence.
	ril1=ril;					// Store values for later rescaling.
	rip1=ripl;
	fact=xnu*xi;
	for (l=nl;l>=1;l--) {
		ritemp=fact*ril+ripl;
		fact -= xi;
		ripl=fact*ritemp+ril;
		ril=ritemp;
	}
	f=ripl/ril;				// Now have unnormalized I_xnu and I'_xnu.
	if (x < XMIN) {		// Use series.
		x2=0.5*x;
		pimu=M_PI*xmu;
		fact = (fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
		d = -log(x2);
		e=xmu*d;
		fact2 = (fabs(e) < EPS ? 1.0 : sinh(e)/e);
		beschb(xmu,&gam1,&gam2,&gampl,&gammi);	// Chebyshev evaluation of Gamma_1 and Gamma_2
		ff=fact*(gam1*cosh(e)+gam2*fact2*d);		// f0.
		sum=ff;
		e=exp(e);
		p=0.5*e/gampl;		// p0.
		q=0.5/(e*gammi);	// q0.
		c=1.0;
		d=x2*x2;
		sum1=p;
		for (i=1;i<=MAXIT;i++) {
			ff=(i*ff+p+q)/(i*i-xmu2);
			c *= (d/i);
			p /= (i-xmu);
			q /= (i+xmu);
			del=c*ff;
			sum += del;
			del1=c*(p-i*ff);
			sum1 += del1;
			if (fabs(del) < fabs(sum)*EPS) break;
		}
		if (i > MAXIT) cerr << "bessk series failed to converge" << endl;
		rkmu=sum;
		rk1=sum1*xi2;
	} else {					// Evaluate CF2 by Steed's algorithm
		b=2.0*(1.0+x);	//	(§5.2), which is OK because there
		d=1.0/b;				//	can be no zero denominators.
		h=delh=d;
		q1=0.0;					// Initializations for recurrence (6.7.35).
		q2=1.0;
		a1=0.25-xmu2;
		q=c=a1;					// First term in equation (6.7.34).
		a = -a1;
		s=1.0+q*delh;
		for (i=2;i<=MAXIT;i++) {
			a -= 2*(i-1);
			c = -a*c/i;
			qnew=(q1-b*q2)/a;
			q1=q2;
			q2=qnew;
			q += c*qnew;
			b += 2.0;
			d=1.0/(b+a*d);
			delh=(b*d-1.0)*delh;
			h += delh;
			dels=q*delh;
			s += dels;
			if (fabs(dels/s) < EPS) break;
			// Need only test convergence of sum since CF2 itself converges more quickly.
		}
		if (i > MAXIT) cerr << "bessik: failure to converge in cf2" << endl;
			h=a1*h;
			rkmu=sqrt(M_PI/(2.0*x))*exp(-x)/s;	// Omit the factor exp(-x) to scale
			rk1=rkmu*(xmu+x+0.5-h)*xi;				//	all the returned functions by exp(x)
		}																		//	for x >= XMIN.											
	rkmup=xmu*xi*rkmu-rk1;
	rimu=xi/(f*rkmu-rkmup);								// Get I_xnu from Wronskian.
	*ri=(rimu*ril1)/ril;									// Scale original I_xnu and I'_xnu.
	*rip=(rimu*rip1)/ril;
	for (i=1;i<=nl;i++) {									// Upward recurrence of K_xnu.
		rktemp=(xmu+i)*xi2*rk1+rkmu;
		rkmu=rk1;
		rk1=rktemp;
	}
	*rk=rkmu;
	*rkp=xnu*xi*rkmu-rk1;
}

void beschb(double x, double *gam1, double *gam2, double *gampl, double *gammi)
/* 
 * Evaluates Gamma_1 and Gamma_2 by Chebyshev expansion for |x|<=1/2. Also returns
 * 1/Gamma_1(1+x) and 1/Gamma_1(1-x). If converting to double precision,
 * set NUSE1 = 7, NUSE2 = 8. */
{
	double xx;
	static double c1[] = {
		-1.142022680371168e0,6.5165112670737e-3,
		3.087090173086e-4,-3.4706269649e-6,6.9437664e-9,
		3.67795e-11,-1.356e-13};
	static double c2[] = {
		1.843740587300905e0,-7.68528408447867e-2,
		1.2719271366546e-3,-4.9717367042e-6,-3.31261198e-8,
		2.423096e-10,-1.702e-13,-1.49e-15};

	xx=8.0*x*x-1.0;											// Multiply x by 2 to make range be -1 to 1,
	*gam1=chebev(-1.0,1.0,c1,7,xx); //	and then apply transformation for evaluating
	*gam2=chebev(-1.0,1.0,c2,8,xx);	//	even Chebyshev series.
	*gampl= *gam2-x*(*gam1);
	*gammi= *gam2+x*(*gam1);
}


double chebev(double a, double b, double c[], int m, double x)
/* Chebyshev evaluation: All arguments are input. c[0..m-1] is an array of Chebyshev coef-
		cients, the first m elements of c output from chebft (which must have been called with
		the same a and b). The Chebyshev polynomial SUM_{k=0}^{m-1} c_k T_k(y) - c_0 = 2 is
		evaluated at a point y = [x-(b+a)/2]=[(b-a)/2], and the result is returned as the
		function value. */
{
	double d=0.0,dd=0.0,sv,y,y2;
	int j;

	if ((x-a)*(x-b) > 0.0) cerr << "x not in range in routine chebev" << endl;
	y2=2.0*(y=(2.0*x-a-b)/(b-a));	// Change of variable.
	for (j=m-1;j>=1;j--) {				// Clenshaw's recurrence.
		sv=d;
		d=y2*d-dd+c[j];
		dd=sv;
	}

	return(y*d-dd+0.5*c[0]);				// Last step is different.
}
#endif

