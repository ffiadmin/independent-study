      SUBROUTINE wqromb(c,lpow,alpha,ss)
      INTEGER JMAX,JMAXP,K,KM,lpow,it
      REAL*8 a,b,ss,EPS,alpha,c
      EXTERNAL wfunc
      PARAMETER (b=200.0d+00,
     + EPS=1.d-14, JMAX=150, JMAXP=JMAX+1, K=14, KM=K-1)
CU    USES polint,trapzd
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)
      a = c
      h(2)=1.d+00
      do 11 j=1,JMAX
        call wtrpzd(c,a,b,s(j),j,it,lpow,alpha)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.d+00,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25d+00*h(j)
11    continue
      pause 'too many steps in qromb'
      END


      SUBROUTINE wtrpzd(c,a,b,s,n,it,lpow,alpha)
      INTEGER n,lpow
      REAL*8 a,b,c,s,alpha
      EXTERNAL wfunc
      INTEGER it,j
      REAL*8 del,sum,tnm,x
      if (n.eq.1) then
        s=0.5d+00*(b-a)*(wfunc(c,a,lpow,alpha)+
     +                   wfunc(c,b,lpow,alpha))
      it=1
      else
        tnm=dfloat(it)
        del=(b-a)/tnm
        x=a+0.5d+00*del
        sum=0.0d+00
        do 11 j=1,it
          sum=sum+wfunc(c,x,lpow,alpha)
          x=x+del
11      continue
        s=0.5d+00*(s+(b-a)*sum/tnm)
        it=2*it
        
      endif
      return
      END


      REAL*8 FUNCTION wfunc(c,x,lpow,alpha)
C     integral of FUNC
      REAL*8 x,alpha,c
      INTEGER lpow


      wfunc=(x-c)**lpow*dexp(-alpha*x*x)
      END

