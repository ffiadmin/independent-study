      SUBROUTINE qromb(lpow,alpha,ss)
      INTEGER JMAX,JMAXP,K,KM,lpow,it
      REAL*8 a,b,ss,EPS,alpha
      EXTERNAL func
      PARAMETER (a=0.0d+00, b=80.0d+00,
     + EPS=1.d-14, JMAX=50, JMAXP=JMAX+1, K=14, KM=K-1)
CU    USES polint,trapzd
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)
      h(2)=1.d+00
      do 11 j=1,JMAX
        call trapzd(a,b,s(j),j,it,lpow,alpha)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.d+00,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25d+00*h(j)
11    continue
      pause 'too many steps in qromb'
      END


      SUBROUTINE trapzd(a,b,s,n,it,lpow,alpha)
      INTEGER n,lpow
      REAL*8 a,b,s,alpha
      EXTERNAL func
      INTEGER it,j
      REAL*8 del,sum,tnm,x
      if (n.eq.1) then
        s=0.5d+00*(b-a)*(func(a,lpow,alpha)+
     +                   func(b,lpow,alpha))
      it=1
      else
        tnm=dfloat(it)
        del=(b-a)/tnm
        x=a+0.5d+00*del
        sum=0.0d+00
        do 11 j=1,it
          sum=sum+func(x,lpow,alpha)
          x=x+del
11      continue
        s=0.5d+00*(s+(b-a)*sum/tnm)
        it=2*it
        
      endif
      return
      END


      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=30)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=dabs(x-xa(1))
      do 11 i=1,n
        dift=dabs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.d+00)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END

      REAL*8 FUNCTION func(x,lpow,alpha)
C     integral of FUNC
      REAL*8 x,alpha
      INTEGER lpow


      func=7.5d+00*(x**4)*dexp(-alpha*x*x)*dexp(-x)
      END

