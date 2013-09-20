//This file is a function library amassed from the sw files

#include "model.h"
#include <string>
#include "fortran.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "gsl_sf_gamma.h"

static double fB(int i, int l1, int l2, double px, double ax, double bx, 
		 int r, double g){
  return binomial_prefactor(i,l1,l2,px-ax,px-bx)*Bfunc(i,r,g);
}

static double Bfunc(int i, int r, double g){
  return fact_ratio2(i,r)*pow(4*g,(double)(r-i));
}

struct Len{
	int		len, l, m, n, i;
	
	double	*exps, *coefs, *norms,
			x, y, z;	
};

static double contr_coulomb(Len& a, Len& b, Len& c, Len& d){

  int i,j,k,l;
  double ret = 0, incr = 0;

  for (i=0; i<a.len; i++)
	  for (j=0; j<b.len; j++)
		  for (k=0; k<c.len; k++)
			  for (l=0; l<d.len; l++){
				  a.i = i; b.i = j; c.i = k; d.i = l;
				  incr = coulomb_repulsion(a,b,c,d);
				  ret += a.coefs[i]*b.coefs[j]*c.coefs[k]*d.coefs[l]*incr;
	}
  return ret;
}

struct PQ{ double xp, yp, zp, xq, yq, zq; };

//coulumb_repulsion object rewrite
static double coulomb_repulsion(Len& a, Len& b, Len& c, Len& d){

  double rab2, rcd2,rpq2,xp,yp,zp,xq,yq,zq,gamma1,gamma2,delta,sum;
  double *Bx, *By, *Bz;
  int I,J,K;

  rab2 = dist2(a,b);
  rcd2 = dist2(c,d);

  PQ pqs = product_center_1D_all(a,b,c,d);

  //======================================================//

  rpq2 = dist2(xp,yp,zp,xq,yq,zq);
  gamma1 = alphaa+alphab;
  gamma2 = alphac+alphad;
  delta = (1./gamma1+1./gamma2)/4.;

  Bx = B_array(la,lb,lc,ld,xp,xa,xb,xq,xc,xd,gamma1,gamma2,delta);
  By = B_array(ma,mb,mc,md,yp,ya,yb,yq,yc,yd,gamma1,gamma2,delta);
  Bz = B_array(na,nb,nc,nd,zp,za,zb,zq,zc,zd,gamma1,gamma2,delta);

  sum = 0.;
  for (I=0; I<la+lb+lc+ld+1;I++)
    for (J=0; J<ma+mb+mc+md+1;J++)
      for (K=0; K<na+nb+nc+nd+1;K++)
	sum += Bx[I]*By[J]*Bz[K]*Fgamma(I+J+K,0.25*rpq2/delta);

  free(Bx);
  free(By);
  free(Bz);  
  
  return 2.*pow(M_PI,2.5)/(gamma1*gamma2*sqrt((double)(gamma1+gamma2)))
    *exp(-alphaa*alphab*rab2/gamma1) 
    *exp(-alphac*alphad*rcd2/gamma2)*sum*norma*normb*normc*normd;
}

static double *B_array(int l1, int l2, int l3, int l4, double p, double a,
		double b, double q, double c, double d,
		double g1, double g2, double delta){
  int Imax,i1,i2,r1,r2,u,I,i;
  double *B;
  Imax = l1+l2+l3+l4+1;
  B = (double *)malloc(Imax*sizeof(double));
  for (i=0; i<Imax; i++) B[i] = 0.;

  for (i1=0; i1<l1+l2+1; i1++)
    for (i2=0; i2<l3+l4+1; i2++)
      for (r1=0; r1<i1/2+1; r1++)
	for (r2=0; r2<i2/2+1; r2++)
	  for (u=0; u<(i1+i2)/2-r1-r2+1; u++){
	    I = i1+i2-2*(r1+r2)-u;
	    B[I] = B[I] + B_term(i1,i2,r1,r2,u,l1,l2,l3,l4,
				 p,a,b,q,c,d,g1,g2,delta);
	  }

  return B;
}

static double Model_pot(int l1, int l2, int m1, int m2, int n1, int n2, double alp1, double alp2){

  	double answer;
	double alp = alp1 + alp2;
        int l=0;

  	FORTRANNAME(qromb)(&l,&alp,&answer);                     
        return answer;  
}



static double Zeta(double x1, double x2, int l1, int l2, double alpha1, double alpha2, double c){

	int npower=2;
       	double alpsum=alpha1+alpha2;

// Calculate integrals over CAP per Cederbaum & Santra
	double Total = 0.0;
	double t1,t2,t3,t4,t5,t6,t7,t8;
	for(int rho=0; rho <= l1; rho++){
		for (int sigma=0; sigma<=l2; sigma++){
              //     cout << rho << "  rho  "  <<sigma << "  sigma"  << endl;
		 	t1 = binomial(l1,rho)*binomial(l2,sigma);
   		        t2 = product_center_1D(alpha1,x1,alpha2,x2);
   		        t3 = pow((t2-x1),(double)(l1-rho));
   		        t4 = pow((t2-x2),(double)(l1-sigma));

//			cout << rho << "  rho  " << sigma << "  sigma" <<endl;
//			cout << t1 << "  " << "t1" << endl;
//			cout << t2 << "  " << "t2" << endl;
//			cout << t3 << "  " << "t3" << endl;
//			cout << t4 << "  " << "t4" << endl;

                        double Iterm = 0.0;

			for (int tau=0; tau<=npower; tau++){
                                double kappa = (double)(rho+sigma);
                                double kaptau = kappa+(double)tau;
				t5 = kaptau + 1.0;
                                t5 =-0.5*t5;
   				t6=pow(-1.0,(double)(npower-tau))*binomial(npower,tau)*pow(alpsum,t5);
                                t7= pow(-1.0,(kappa))*pow((c+t2),(double)(npower-tau));
//			cout << t7 << "  " << "t7" << tau << "  tau" << endl;

				t7 = t7*DELTA(kaptau,alpsum,c+t2);
                                t8= pow((c-t2),(double)(npower-tau))*DELTA(kaptau,alpsum,c-t2);
				Iterm = Iterm+t6*(t7+t8);
//			cout << t5 << "  " << "t5" << tau << "  tau" <<endl;
//			cout << t6 << "  " << "t6  " << alpsum << "  alpsum" <<  endl;
//			cout << t7 << "  " << "t7" << tau << "  tau" << endl;
//			cout << t8 << "  " << "t8" << tau << "  tau" <<endl;
//			cout << Iterm << "  " << "Iterm" << tau << "  tau" <<endl;
			}

			Iterm = 0.5*Iterm;
      			Total = Total + t1*t3*t4*Iterm;
		}
	}

	return Total;
}

static double DELTA(double k, double a, double c){

//        cout << a << "  a  " << c << "  c  " <<endl;
	double t1,t2,t3,t4,t5,t6;

 	t1 = 0.5*(k+1.0);
	t2 = a*c*c;
	t3 = gsl_sf_gamma(t1);
	t4 = gsl_sf_gamma_inc_P(t1,t2);
//         cout << t4 << "  t4" <<endl;
        t4 = t4*t3;
//        cout << "in delta" << endl;
//	cout << t1 <<  "  "  << t2 << "  " << t3 << "  " << t4 << endl;
	if (c<0){
		 t5=pow(-1.0,k+1.0);
	}
	else{
		t5=1.0;
	}

	t6=    t3-t4*t5;
//        cout << t5 << "t5 in delta" <<endl;
//        cout << t6 << "t6 in delta" <<endl;
        return t6;
}


static double B_term(int i1, int i2, int r1, int r2, int u, int l1, int l2,
	      int l3, int l4, double Px, double Ax, double Bx,
	      double Qx, double Cx, double Dx, double gamma1,
	      double gamma2, double delta){
  /* THO eq. 2.22 */
  return fB(i1,l1,l2,Px,Ax,Bx,r1,gamma1)
    *pow(-1,(double)(i2))*fB(i2,l3,l4,Qx,Cx,Dx,r2,gamma2)
    *pow(-1,(double)(u))*fact_ratio2(i1+i2-2*(r1+r2),u)
    *pow(Qx-Px,(double)(i1+i2-2*(r1+r2)-2*u))
    /pow(delta,(double)(i1+i2-2*(r1+r2)-u));
}


static double kinetic(double alpha1, int l1, int m1, int n1,
	       double xa, double ya, double za,
	       double alpha2, int l2, int m2, int n2,
	       double xb, double yb, double zb){

  double term0,term1,term2;
  term0 = alpha2*(2*(l2+m2+n2)+3)*
    overlap(alpha1,l1,m1,n1,xa,ya,za,
		   alpha2,l2,m2,n2,xb,yb,zb);
  term1 = -2*pow(alpha2,(double)(2))*
    (overlap(alpha1,l1,m1,n1,xa,ya,za,
		    alpha2,l2+2,m2,n2,xb,yb,zb)
     + overlap(alpha1,l1,m1,n1,xa,ya,za,
		      alpha2,l2,m2+2,n2,xb,yb,zb)
     + overlap(alpha1,l1,m1,n1,xa,ya,za,
		      alpha2,l2,m2,n2+2,xb,yb,zb));
  term2 = -0.5*(l2*(l2-1)*overlap(alpha1,l1,m1,n1,xa,ya,za,
					 alpha2,l2-2,m2,n2,xb,yb,zb) +
		m2*(m2-1)*overlap(alpha1,l1,m1,n1,xa,ya,za,
					 alpha2,l2,m2-2,n2,xb,yb,zb) +
		n2*(n2-1)*overlap(alpha1,l1,m1,n1,xa,ya,za,
					 alpha2,l2,m2,n2-2,xb,yb,zb));
  return term0+term1+term2;
}

static double overlap(double alpha1, int l1, int m1, int n1,
		      double xa, double ya, double za,
		      double alpha2, int l2, int m2, int n2,
		      double xb, double yb, double zb){
  /*Taken from THO eq. 2.12*/
  double rab2,gamma,xp,yp,zp,pre,wx,wy,wz;

  rab2 = dist2(xa,ya,za,xb,yb,zb);
  gamma = alpha1+alpha2;
   xp = product_center_1D(alpha1,xa,alpha2,xb);
  yp = product_center_1D(alpha1,ya,alpha2,yb);
  zp = product_center_1D(alpha1,za,alpha2,zb);

  pre = pow(M_PI/gamma,1.5)*exp(-alpha1*alpha2*rab2/gamma);

  wx = overlap_1D(l1,l2,xp-xa,xp-xb,gamma);
  wy = overlap_1D(m1,m2,yp-ya,yp-yb,gamma);
  wz = overlap_1D(n1,n2,zp-za,zp-zb,gamma);
  return pre*wx*wy*wz;
}

static double overlap_1D(int l1, int l2, double PAx,
			 double PBx, double gamma){
  /*Taken from THO eq. 2.12*/
  int i;
  double sum;
  sum = 0.;
  for (i=0; i<(1+floor(0.5*(l1+l2))); i++)
    sum += binomial_prefactor(2*i,l1,l2,PAx,PBx)* 
      fact2(2*i-1)/pow(2*gamma,(double)(i));
  return sum;
}
    
static double nuclear_attraction(double x1, double y1, double z1, double norm1,
				 int l1, int m1, int n1, double alpha1,
				 double x2, double y2, double z2, double norm2,
				 int l2, int m2, int n2, double alpha2,
				 double x3, double y3, double z3){
  int I,J,K;
  double gamma,xp,yp,zp,sum,rab2,rcp2;
  double *Ax,*Ay,*Az;

  gamma = alpha1+alpha2;

  xp = product_center_1D(alpha1,x1,alpha2,x2);
  yp = product_center_1D(alpha1,y1,alpha2,y2);
  zp = product_center_1D(alpha1,z1,alpha2,z2);

  rab2 = dist2(x1,y1,z1,x2,y2,z2);
  rcp2 = dist2(x3,y3,z3,xp,yp,zp);

  Ax = A_array(l1,l2,xp-x1,xp-x2,xp-x3,gamma);
  Ay = A_array(m1,m2,yp-y1,yp-y2,yp-y3,gamma);
  Az = A_array(n1,n2,zp-z1,zp-z2,zp-z3,gamma);

  sum = 0.;
  for (I=0; I<l1+l2+1; I++)
    for (J=0; J<m1+m2+1; J++)
      for (K=0; K<n1+n2+1; K++)
	sum += Ax[I]*Ay[J]*Az[K]*Fgamma(I+J+K,rcp2*gamma);

  free(Ax);
  free(Ay);
  free(Az);
  return -norm1*norm2*
    2*M_PI/gamma*exp(-alpha1*alpha2*rab2/gamma)*sum;
}
    
static double A_term(int i, int r, int u, int l1, int l2,
		     double PAx, double PBx, double CPx, double gamma){
  /* THO eq. 2.18 */
  return pow(-1,(double)(i))*binomial_prefactor(i,l1,l2,PAx,PBx)*
    pow(-1,(double)(u))*fact(i)*pow(CPx,(double)(i-2*r-2*u))*
    pow(0.25/gamma,(double)(r+u))/fact(r)/fact(u)/fact(i-2*r-2*u);
}

static double *A_array(int l1, int l2, double PA, double PB,
		double CP, double g){
  /* THO eq. 2.18 and 3.1 */
  int Imax,i,r,u,I;
  double *A;

  Imax = l1+l2+1;
  A = (double *)malloc(Imax*sizeof(double));
  for (i=0; i<Imax; i++) A[i] = 0.;
  for (i=0; i<Imax; i++)
    for (r=0; r<floor((double)(i/2))+1;r++)
      for (u=0; u<floor((i-2*r)/2.)+1; u++){
	I = i-2*r-u;
	A[I] += A_term(i,r,u,l1,l2,PA,PB,CP,g);
      }
  return A;
}


static int fact(int n){
  if (n <= 1) return 1;
  return n*fact(n-1);
}

static int fact2(int n){ /* double factorial function = 1*3*5*...*n */
  if (n <= 1) return 1;
  return n*fact2(n-2);
}

static double dist2(double x1, double y1, double z1,
		    double x2, double y2, double z2){
  return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
}

//object rewrite
static double dist2(Len& a, Len& b){
  return (a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z);
}

static double dist(double x1, double y1, double z1,
		   double x2, double y2, double z2){
  return sqrt(dist2(x1,y1,z1,x2,y2,z2));
}

//object rewrite
static double dist(Len& a, Len& b){
  return sqrt(dist2(a, b));
}

static double binomial_prefactor(int s, int ia, int ib, double xpa, double xpb){
  int t;
  double sum=0.;
  for (t=0; t<s+1; t++)
    if ((s-ia <= t) && (t <= ib)) 
      sum += binomial(ia,s-t)*binomial(ib,t)*pow(xpa,(double)(ia-s+t))*pow(xpb,(double)(ib-t));
  return sum;
} 

static int binomial(int a, int b){return fact(a)/(fact(b)*fact(a-b));}

static double Fgamma(double m, double x){
  double val;
  if (fabs(x) < SMALL) x = SMALL;
  val = gamm_inc(m+0.5,x);
  /* if (val < SMALL) return 0.; */ /* Gives a bug for D orbitals. */
  return 0.5*pow(x,-m-0.5)*val; 
}

static double gamm_inc(double a, double x){ /* Taken from NR routine gammap */
  double gamser,gammcf,gln;
  
  assert (x >= 0.);
  assert (a > 0.);
  if (x < (a+1.0)) {
    gser(&gamser,a,x,&gln);
    return exp(gln)*gamser;
  } else {
    gcf(&gammcf,a,x,&gln);
    return exp(gln)*(1.0-gammcf);
  }
}
 
static void gser(double *gamser, double a, double x, double *gln){
  int n;
  double sum,del,ap;

  *gln=lgamma(a);
  if (x <= 0.0) {
    assert(x>=0.);
    *gamser=0.0;
    return;
  } else {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=ITMAX;n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) {
	*gamser=sum*exp(-x+a*log(x)-(*gln));
	return;
      }
    }
    printf("a too large, ITMAX too small in routine gser");
    return;
  }
}
 
static void gcf(double *gammcf, double a, double x, double *gln){
  int i;
  double an,b,c,d,del,h;
  
  *gln=lgamma(a);
  b=x+1.0-a;
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;
  for (i=1;i<=ITMAX;i++) {
    an = -i*(i-a);
    b += 2.0;
    d=an*d+b;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=b+an/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;
  }
  assert(i<=ITMAX);
  *gammcf=exp(-x+a*log(x)-(*gln))*h;
}

static int ijkl2intindex(int i, int j, int k, int l){
  int tmp,ij,kl;
  if (i<j){
    tmp = i;
    i = j;
    j = tmp;
  }
  if (k<l){
    tmp = k;
    k = l;
    l = tmp;
  }
  ij = i*(i+1)/2+j;
  kl = k*(k+1)/2+l;
  if (ij<kl){
    tmp = ij;
    ij = kl;
    kl = tmp;
  }
  return ij*(ij+1)/2+kl;
}

static int fact_ratio2(int a, int b){ return fact(a)/fact(b)/fact(a-2*b); }

//object rewrite
static double product_center_1D_all(Len& a, Len& b, Len& c, Len& d){
  return (a.norms[a.i]*a.x + b.norms[b.i]*b.x)/(a.norms[a.i] + b.norms[b.i]);
}

static double three_center_1D(double xi, int ai, double alphai,
			      double xj, int aj, double alphaj,
			      double xk, int ak, double alphak){

  double gamma, dx, px, xpi,xpj,xpk,intgl;
  int q,r,s,n;
  
  gamma = alphai+alphaj+alphak;
  dx = exp(-alphai*alphaj*pow(xi-xj,2.)/gamma) *
    exp(-alphai*alphak*pow(xi-xk,2.)/gamma) *
    exp(-alphaj*alphak*pow(xj-xk,2.)/gamma);
  px = (alphai*xi+alphaj*xj+alphak*xk)/gamma;
    
  xpi = px-xi;
  xpj = px-xj;
  xpk = px-xk;
  intgl = 0;
  for (q=0; q<ai+1; q++){
    for (r=0; r<aj+1; r++){
      for (s=0; s<ak+1; s++){
	n = (q+r+s)/2;
	if ((q+r+s)%2 == 0) {
	  intgl += binomial(ai,q)*binomial(aj,r)*binomial(ak,s)*
	    pow(xpi,(double)(ai-q))*pow(xpj,(double)(aj-r))*pow(xpk,(double)(ak-s))*
	    fact2(2*n-1)/pow(2*gamma,(double)(n))*sqrt(M_PI/gamma);
	}
      }
    }
  }
  return dx*intgl;
}

