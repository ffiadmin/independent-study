c     **** main routine for complex coordinate rotation
c     **** uses integrals fromPRA vol 20 page 814 1979
c     **** uses standard intgral formulas for gto's
      
      EXTERNAL gtest
      EXTERNAL fact2
      EXTERNAL fct
      EXTERNAL OL1D  
      EXTERNAL PBINM 
      EXTERNAL BIN 

      integer maxbas,maxnuc
      parameter (maxstp=500)
      parameter (maxbas=150)
      parameter (maxbs3=450)
      parameter (maxnuc=20)

      double precision alpha(maxbas),
     + x(maxbas),y(maxbas),z(maxbas),
     + nchg(maxnuc),nx(maxnuc),ny(maxnuc),nz(maxnuc),
     + repul(maxnuc),range(maxnuc),repnrg,
     + PI,t1,t2,t3,t4,norm(maxbas),min,max,fact2

      double precision ovrlp,r1
      double precision kenergy

      double precision     smat(maxbas,maxbas),
     + rr1,rr2,alpreal,theta,alpstrt,
     + thstart,thstep,svec(maxbas),
     + aux(maxbs3),rwork(maxbs3),alpstp,
     + dermag,dermin,ang(maxstp),angmin


      double complex scale,vint(maxbas,maxbas),
     + tint(maxbas,maxbas),xmatx(maxbas,maxbas),
     + xmat(maxbas,maxbas),core(maxbas,maxbas),
     + xtmp(maxbas,maxbas),XX1,XX2,eigval(maxbas),
     + VL(maxbas,maxbas),VR(maxbas,maxbas),
     + work(maxbs3),venergy,c6,eta(maxstp),
     + val(maxstp),deriv(maxstp),etamin,valmin

      integer l(maxbas), m(maxbas), n(maxbas),
     + nbasis,nnuc,nthet,naux,info,kk,nalp

      character*80 title

      PI=dacos(-1.0d+00)
c      write(6,*)PI

      scale = 1.0d+00
      XX1=(1.0d+00,0.0d+00)
      XX2=(0.0d+00,0.0d+00)

c     read input
 
      read(5,*)title
      write(6,*)title

      read(5,*)nbasis
      write(6,*) 'there are  ',nbasis,'  basis functions'

c     note: all coordinates are read in in angstroms and then converted to au before printing

      do 10 i = 1, nbasis
           read(5,*)l(i),m(i),n(i),x(i),y(i),z(i),alpha(i)
           x(i) = x(i)/0.529177249
           y(i) = y(i)/0.529177249
           z(i) = z(i)/0.529177249
           write(6,*)l(i),m(i),n(i),x(i),y(i),z(i),alpha(i)
10    continue

      read(5,*)nnuc
      write(6,*)'there are  ',nnuc,'  nuclei.'


      do 20 i = 1,nnuc
           read(5,*)nx(i),ny(i),nz(i),nchg(i),repul(i),range(i)
           nx(i) = nx(i)/0.529177249
           ny(i) = ny(i)/0.529177249
           nz(i) = nz(i)/0.529177249
           write(6,*)nx(i),ny(i),nz(i),nchg(i),repul(i),range(i)
20    continue

      read(5,*)nthet,nalp,alpstrt,alpstp,thstart,thstep

      write(6,*)'the initial values for r and theta are'
      write(6,*)alpstrt,thstart

      write(6,*)'the program will take ',nthet,' steps'
      write(6,*)'each step in theta will be ',thstep,' radians'
      write(6,*)'the program will take ',nalp,' steps'
      write(6,*)'each step in r will be ',alpstp,' au'


      read(5,*)min,max
      write(6,*) 'the range for the real part of the energy is' 
      write(6,*)'from ',min,' to ',max

c     calculate normalization constants

      do 30 i = 1,nbasis
           t1 = 2**(l(i)+m(i)+n(i))
           t2 = alpha(i)**(2*l(i)+2*m(i)+2*n(i)+3)
           t2 = t2**0.25
           t3 = (2.0d+00/PI)**0.75
           t4 = fact2(2*l(i)-1)*fact2(2*m(i)-1)*fact2(2*n(i)-1)
           t4= dsqrt(t4)
           norm(i) = t1*t2*t3/t4
c          write(6,*)norm(i)
30    continue

      do 11 ir=1,nbasis
           do 12 jc=1,nbasis
                call overlap(alpha(ir),l(ir)
     +          ,m(ir),n(ir),x(ir),y(ir),z(ir),
     +          alpha(jc),l(jc),m(jc),n(jc),
     +          x(jc),y(jc),z(jc),ovrlp)

                r1 = norm(ir)*norm(jc)*ovrlp
                smat(ir,jc)=r1                   


                call kinetic
     +          (alpha(ir),l(ir),m(ir),n(ir),x(ir),y(ir),z(ir),
     +          alpha(jc),l(jc),m(jc),n(jc),x(jc),y(jc),z(jc),
     +           kenergy)

                r1 =norm(ir)*norm(jc)*kenergy
                tint(ir,jc)=dcmplx(r1,0.0d+00)
12         continue
11     continue

      info=0
      naux=maxbs3
      call dsyev('V','U',nbasis,smat,maxbas,svec,aux,naux,info)

      do 98 iii = 1,nbasis
c     write(6,*)svec(iii)
98    continue

      alpreal=alpstrt

      do 888 istep=1,nalp 

      theta = thstart

      do 99 nscan=1,nthet
           rr1=alpreal*cos(theta)
           rr2=-alpreal*sin(theta)

           scale = dcmplx(rr1,rr2)
c          write(6,*)'scale ',scale


      do 40 ir=1,nbasis
           do 50 jc=1,nbasis


                vint(ir,jc)=(0.0d+00,0.0d+00)

              do 60 kk=1,nnuc     
                    call attrct(
     +  x(ir),y(ir),z(ir),norm(ir),l(ir),m(ir),n(ir),alpha(ir),
     +  x(jc),y(jc),z(jc),norm(jc),l(jc),m(jc),n(jc),alpha(jc),
     +  nx(kk),ny(kk),nz(kk),venergy,scale)

                    c6 = dcmplx(nchg(kk),0.0d+00)
                    vint(ir,jc)=vint(ir,jc)+c6*venergy

                    call excld(
     +  x(ir),y(ir),z(ir),norm(ir),l(ir),m(ir),n(ir),alpha(ir),
     +  x(jc),y(jc),z(jc),norm(jc),l(jc),m(jc),n(jc),alpha(jc),
     +  nx(kk),ny(kk),nz(kk),repul(kk),range(kk),repnrg)        

                    c6 = dcmplx(repnrg,0.0d+00)
                    vint(ir,jc)=vint(ir,jc)+c6

c       write(6,*) 'row, col, V ',ir,jc,vint(ir,jc),venergy,nchg(kk)


60              continue



                core(ir,jc)=scale*scale*tint(ir,jc)+vint(ir,jc)

c               write(6,*) 'row, col, S ',ir,jc,smat(ir,jc)
c               write(6,*) 'row, col, T ',ir,jc,tint(ir,jc)
c               write(6,*)'core  ',ir,jc,core(ir,jc)

50         continue
40    continue

C     Having built the ie matricies find the eigenvalues/vects of Smat


      do 70 ir = 1,nbasis
          do 80 jc = 1,nbasis
              r1 = smat(ir,jc)/dsqrt(svec(jc))
              xmat(ir,jc)=dcmplx(r1,0.0d+00)
              xmatx(jc,ir)=xmat(ir,jc)
c         write(6,*)xmat(ir,jc),xmatx(ir,jc)
80        continue
70    continue

      call zgemm('N','N',nbasis,nbasis,nbasis,XX1,core,maxbas,
     +   xmat,maxbas,XX2,xtmp,maxbas)

      call zgemm('N','N',nbasis,nbasis,nbasis,XX1,xmatx,maxbas,
     +   xtmp,maxbas,XX2,core,maxbas)

      call zgeev('V','V',nbasis,core,maxbas,eigval,VL,maxbas,
     +   VR,maxbas,work,maxbs3,rwork,info)

      do 86 iii = 1,nbasis
            eigval(iii)=(27.2114d+00,0.0d+00)*eigval(iii)
            ctest=dreal(eigval(iii))
            if ((ctest.ge.min).and.(ctest.le.max)) then
            write(6,*)eigval(iii)
            eta(nscan)=scale
            val(nscan)=eigval(iii)
            ang(nscan)=theta
            endif
86    continue

c     call zgemm('N','N',nbasis,nbasis,nbasis,XX1,VR,maxbas,
c    +   xmatx,maxbas,XX2,xtmp,maxbas)

c     call zgemm('N','N',nbasis,nbasis,nbasis,XX1,xmat,maxbas,
c    +   xtmp,maxbas,XX2,VR,maxbas)

c     do 88 jc =1,nbasis
c        do 89 ir = 1,nbasis
c           write(6,*)ir,jc,VR(ir,jc)
c89       continue
c88    continue
      
      
      theta=theta+thstep

99    continue

      dermin=1.0d+00

      do 999 i=1,nthet-1
          deriv(i)=(val(i+1)-val(i))/
     +    (eta(i+1)-eta(i))
          deriv(i)=deriv(i)
          dermag=dreal(deriv(i)*dconjg(deriv(i)))
          if (dermag.lt.dermin) then
              dermin=dermag
              etamin=eta(i)
              valmin=val(i)
              angmin=ang(i)
          endif
999   continue

c     write(6,*)alpreal,angmin,etamin,valmin,dermin

      alpreal=alpreal+alpstp

888   continue          
      stop
      end



      double precision FUNCTION fact2(n)
C     double factorial function
      INTEGER n


      fact2=1.0d+00

      if ((n.eq.0).or.(n.eq.1)) fact2=1.0d+00
      if (n.eq.3) fact2=3.0d+00
      if (n.eq.5) fact2=15.0d+00
      if (n.eq.7) fact2=15.0d+00*7.0d+00
      if (n.eq.9) fact2=15.0d+00**7.0d+00*9.0d+00
c     write(6,*)n,fact2

      END

      SUBROUTINE overlap(a1,l1,m1,n1,x1,y1,z1,
     +    a2,l2,m2,n2,x2,y2,z2,ovrlp)

      double precision a1,a2,ovrlp,x1,x2,y1,y2,z1,z2
      integer l1,l2,m1,m2,n1,n2

      double precision rab2,gamma,xp,yp,zp,pre,wx,wy,wz,
     +PI,OL1D
    

      PI=dacos(-1.0d+00)

      rab2 = (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2
      gamma = a1+a2

      xp = (a1*x1+a2*x2)/(a1+a2)
      yp = (a1*y1+a2*y2)/(a1+a2)
      zp = (a1*z1+a2*z2)/(a1+a2)

      pre = (PI/gamma)**1.5d+00*dexp(-a1*a2*rab2/gamma)

c     write(6,*)'xp,yp,zp,pre'
c     write(6,*)xp,yp,zp,pre

      wx = OL1D(l1,l2,xp-x1,xp-x2,gamma)
      wy = OL1D(m1,m2,yp-y1,yp-y2,gamma)
      wz = OL1D(n1,n2,zp-z1,zp-z2,gamma)

c     write(6,*)'wx,wy,wz'     
c     write(6,*)l1,l2,wx,wy,wz     

      ovrlp = pre*wx*wy*wz

      return
      end

      double precision FUNCTION OL1D(l1,l2,PAx,PAb,gamma)

      integer l1,l2,max
      double precision PAx, PAb, gamma, sum,fact2,
     +PBINM

c     write(6,*)'l1,l2,PAx,PAb,gamma'
c     write(6,*)l1,l2,PAx,PAb,gamma
      max = int(0.5*(l1+l2))
c     write(6,*)'max  ',max, l1,l2
      sum=0.0d+00

      do 10, i = 0,max

           sum = sum + PBINM(2*i,l1,l2,PAx,PAb)*
     +     fact2(2*i-1)/(2.0d+00*gamma)**i

10    continue
c          write(6,*)'sum ' ,sum

      OL1D = sum
      end

      double precision FUNCTION PBINM(s,ia,ib,xpa,xpb)    


      integer s,ia,ib
      double precision xpa,xpb,sum,BIN

c     write(6,*)'xpa,xpb',xpa,xpb
      sum = 0.0d+00

      do 10 i =0,s
            if(((s-ia).le.i).and.(i.le.ib)) then
                sum = sum + BIN(ia,s-i)*BIN(ib,i)*
     +          xpa**(ia-s+i)*xpb**(ib-i)
            endif

10    continue

      PBINM = sum
c     write(6,*)'PBINM',PBINM

      end

      double precision FUNCTION  BIN(a,b)

      integer a,b,c,t1,t2,t3,fct

      t1 = fct(a)
      t2=fct(b)
      t3=fct(a-b)

c     write(6,*)'t1,t2,t3,',t1,t2,t3

      c = fct(a)/(fct(b)*fct(a-b))

c     write(6,*)'a,b,c ',a,b,c
      BIN=dfloat(c)
c     write(6,*)'BIN**',a,b,BIN

      end

      integer FUNCTION fct(a)

      integer a,sum

      sum =1

      do 10  i=0,a
          if (i.ge.1)sum = sum*i
10     continue

c     write(6,*)'a,sum',a,sum
      fct = sum

      end


      SUBROUTINE kinetic
     + (a1,l1,m1,n1,x1,y1,z1,a2,l2,m2,n2,x2,y2,z2,kenergy)

      double precision a1,x1,y1,z1,a2,x2,y2,z2,kenergy,
     +                 OL1,OL2,OL3,OL4

      double precision term0,term1,term2

      integer l1,m1,n1,l2,m2,n2

                call overlap(a1,l1,m1,n1,   
     +          x1,y1,z1,a2,l2,m2,n2,x2,y2,z2,OL1)

                call overlap(a1,l1,m1,n1,   
     +          x1,y1,z1,a2,l2+2,m2,n2,x2,y2,z2,OL2)

                call overlap(a1,l1,m1,n1,   
     +          x1,y1,z1,a2,l2,m2+2,n2,x2,y2,z2,OL3)

                call overlap(a1,l1,m1,n1,   
     +          x1,y1,z1,a2,l2,m2,n2+2,x2,y2,z2,OL4)

      term0 = a2*OL1*dfloat(2*l2+2*m2+2*n2+3)

      term1= -2.0d+00*a2*a2*(OL2+OL3+OL4) 

                call overlap(a1,l1,m1,n1,   
     +          x1,y1,z1,a2,l2-2,m2,n2,x2,y2,z2,OL2)

                call overlap(a1,l1,m1,n1,   
     +          x1,y1,z1,a2,l2,m2-2,n2,x2,y2,z2,OL3)

                call overlap(a1,l1,m1,n1,   
     +          x1,y1,z1,a2,l2,m2,n2-2,x2,y2,z2,OL4)

      term2= -0.5*(dfloat(l2*(l2-1))*OL2 +
     +             dfloat(m2*(m2-1))*OL3 +
     +             dfloat(n2*(n2-1))*OL4 )

      kenergy=term0+term1+term2

      return
      end


      SUBROUTINE attrct(
     +           x1,y1,z1,norm1,l1,m1,n1,a1,                             
     +           x2,y2,z2,norm2,l2,m2,n2,a2,                             
     +           nx,ny,nz,venergy,scale)

      double precision x1,y1,z1,norm1,a1
      double precision x2,y2,z2,norm2,a2
      double precision nx,ny,nz
      double precision gamma,xp,yp,zp,rab2
      double precision Ax(20),Ay(20),Az(20) 
      double precision r1,r2,r3,PI,r7

      integer l1,m1,n1,l2,m2,n2

      double complex venergy,sum,rcp2,c7

      double complex c1,c2,c3,gtest,c4,scale,
     + cx,cy,cz

      PI=dacos(-1.0d+00)

      gamma=a1+a2

      xp = (a1*x1+a2*x2)/(a1+a2)
      yp = (a1*y1+a2*y2)/(a1+a2)
      zp = (a1*z1+a2*z2)/(a1+a2)

      cx = dcmplx(nx,0.0d+00)
      cy = dcmplx(ny,0.0d+00)
      cz = dcmplx(nz,0.0d+00)

      cx =dcmplx(xp,0.0d+00)-scale*nx 
      cy =dcmplx(yp,0.0d+00)-scale*ny 
      cz =dcmplx(zp,0.0d+00)-scale*nz 

      rab2 = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
c     rcp2 = cx*dconjg(cx) + cy*dconjg(cy) + cz*dconjg(cz)
      rcp2 = cx*cx + cy*cy + cz*cz

      call Aarray(l1,l2,xp-x1,xp-x2,xp-nx,gamma,Ax)
      call Aarray(m1,m2,yp-y1,yp-y2,yp-ny,gamma,Ay)
      call Aarray(n1,n2,zp-z1,zp-z2,zp-nz,gamma,Az)

c      do 555 i = 1,20
c     write(6,*)'Ax,Ay,Az',Ax(i),Ay(i),Az(i)
c555   continue

      sum=(0.0d+00,0.0d+00)
      c4=(0.0d+00,0.0d+00)

      do 10 i = 0,l1+l2
           do 20 j=0,m1+m2
                do 30 k = 0,n1+n2
                     r1 = Ax(i+1)*Ay(j+1)*Az(k+1)
                     r2=dfloat(i+j+k)+0.5d+00
                     c1=dcmplx(r2,0.0d+00)
                     c2 = rcp2*dcmplx(gamma,0.0d+00)
                     r3 = dsqrt(dreal(c2 * dconjg(c2)))
                     if(dabs(r3).lt.(1.0d-17))c3=(1.0d-17,0.0)
                     c3=gtest(r2,c2)        
                     sum =sum+dcmplx(r1,0.0d+00)*c3                       
c     write(6,*)Ax(i+1),Ay(j+1),Az(k+1)
c     write(6,*)'i,j,k,sum ',i,j,k,sum,c3,r4

30              continue
20         continue
10    continue


      r7 = -norm1*norm2*2.0d+00*PI/gamma
      r7 = dexp(-a1*a2*rab2/gamma)*r7
      c7 = dcmplx(r7,0.0d+00)
c     write(6,*)'venergy ',venergy
      venergy = scale*c7*sum                           

c     write(6,*)'venergy ',venergy

      return
      end

      SUBROUTINE Aarray(l1,l2,PA,PB,CP,g,A)

      double precision PA,PB,CP,g,A(20),t1
      integer l1,l2,Imax,lim1,lim2,II

c     write(6,*)'entering Aarray, l1,l2,PA,PB,CP,g'
c     write(6,*)l1,l2,PA,PB,CP,g

      Imax=l1+l2

      Do 10, i=1,20
          A(i)=0.0d+00
10    continue

      do 20 i=0,Imax
          lim1 = int(dfloat(i/2))

          do 30 j = 0,lim1
              lim2=int(dfloat(i-2*j)/2.0d+00)

              do 40 k = 0,lim2
                  II = i-2*j-k+1
                  call Aterm(i,j,k,l1,l2,PA,PB,CP,g,t1)
c                 write(6,*)'return Aterm i,r,u,l1,l2,t1,II'
c                 write(6,*)i,j,k,l1,l2,t1,II
                  A(II)=A(II)+t1
c                 write(6,*)'A(I),I  ',A(II),II
40            continue
30        continue
20    continue
c     write(6,*)'leaving  Aarray, l1,l2,PA,PB,CP,g,A'
c     write(6,*)l1,l2,PA,PB,CP,g
c     do 133 kk=1,5
c         write(6,*)A(kk)
c133   continue

      return
      end

      SUBROUTINE Aterm(i,j,k,l1,l2,PA,PB,CP,g,t1)

      integer i,j,k,l1,l2,fct
      double precision PA,PB,CP,g,t1,PBINM

      t1 = (-1)**i
      t1 = t1*PBINM(i,l1,l2,PA,PB)
      t1 = t1*(-1)**k
      t1 = t1 * fct(i)
      t1 = t1 * CP**(i-2*j-2*k)
      t1 = t1 * (0.25/g)**(j+k)
      t1 = t1 / (fct(j)*fct(k)*fct(i-2*j-2*k))

c     write(6,*)'i,r,u,t1',i,j,k,t1

      return
      end

      double Complex FUNCTION gtest(a,x)
      double complex x, sum,c1,c2,c3,sum2,xx2,cdig,b
      double precision a,r1,r2,sumnrm
      integer n,q,fct

      sum = (0.0d+00,0.0d+00)
      xx2 = (0.0d+00,0.0d+00)
      sumnrm = dsqrt(dreal(x * dconjg(x)))
c     write(6,*)'sumnrm',sumnrm 
      if(sumnrm.lt.1.0d+00)then
          do 10 n=0,20
              r1=a+dfloat(n)
              c1 = dcmplx(r1,0.0d+00)
              q = fct(n)
              r2=dreal(q)
              c2=dcmplx(r2,0.0d+00)
              c3=-x
              sum = sum + (c3)**n/(c1*c2)
10        continue
          sum=  (0.5d+00,0.0d+00)*sum
c         write(6,*)'int',sum   
          gtest=sum
          return
      else
c         write(6,*)a,x
          b=dcmplx(a,0.0d+00)
          sum = cdig(b,x)
c         write(6,*)'sum',sum
          sum2 = cdig(b,xx2)
c         write(6,*)'sum2',sum2
          sum = (0.5d+00,0.0d+00)*x**(-a)*(sum2-sum)
c         write(6,*)'cdig',sum   
          gtest = sum
          return
      endif
      return
      end

      double complex function cdig(alpha,x)
c --- Written By Eric Kostlan & Dmitry Gokhman
c --- March 1986
      double complex alpha,x,cdh
      double complex re,one,p,q
      double precision xlim,zero,dnrm
      data re,one/0.36787944117144232,1./
      data xlim,zero/1.,0./ibuf/34/
c --- If x is near the negative real axis, then shift to x=1.
      if(dnrm(x).lt.xlim.or.dreal(x).lt.zero.and.
     + dabs(dimag(x)).lt.xlim)then
      cdig=re/cdh(alpha,one)
      ilim=dreal(x/re)
      do 1 i=0,ibuf-ilim
      call term(alpha,x,i,p,q)
      cdig=cdig+p*q
1     continue
      else

      cdig=cdexp(-x+alpha*cdlog(x))/cdh(alpha,x)
      endif
      return
      end
c
      subroutine term(alpha,x,i,p,q)
c --- Calculate p*q = -1**i(1-x**(alpha+i))/(alpha+i)i! carefully.
      double complex alpha,x,p,q,ci,alphai
      double complex zero,one,two,cdlx
      double precision tol,xlim,dnrm
      data zero,one,two/0.,1.,2./tol/3.d-7/xlim/39./
      if(i.eq.0)q=one
      ci=i
      alphai=alpha+ci
      if(x.eq.zero)then
      p=one/alphai
      if(i.ne.0)q=-q/ci
      return
      endif
      cdlx=cdlog(x)
c --- If (1-x**alphai)=-x**alphai on the computer,
c --- then change the inductive scheme to avoid overflow.
      if(dreal(alphai*cdlx).gt.xlim.and.i.ne.0)then
      p=p*(alphai-one)/alphai
      q=-q*x/ci
      return
      endif
      if(dnrm(alphai).gt.tol)then
      p=(one-x**alphai)/alphai
      else
      p=-cdlx*(one+cdlx*alphai/two)
      endif
      if(i.ne.0)q=-q/ci
      return
      end
c
      double complex function cdh(alpha,x)
c --- Written By Eric Kostlan & Dmitry Gokhman
c --- March 1986
      double complex alpha,x,cdhs
      double complex one,term,sum,cn,alpha1
      double precision buf
      data one/1./buf/0./
c --- If Re(alpha-x) is too big, shift alpha.
      n=dreal(alpha-x)-buf
      if(n.gt.0)then
      cn=n+1
      alpha1=alpha-cn
      term=one/x
      sum=term
      do 1 i=1,n
      cn=n-i+1
      term=term*(alpha1+cn)/x
      sum=term+sum

1     continue
      sum=sum+term*alpha1/cdhs(alpha1,x)
      cdh=one/sum
      else
      cdh=cdhs(alpha,x)
      endif
      return
      end
c
      double complex function cdhs(alpha,x)
c --- Written By Eric Kostlan & Dmitry Gokhman
c --- March 1986
      double complex zero,half,one,alpha,x
      double complex p0,q0,p1,q1,r0,r1,ci,factor
      double precision tol1,tol2,error,dnrm
      data zero,half,one/0.,0.5,1./
      data tol1,tol2,error/1.d10,1.d-10,5.d-18/ilim/100000/
      q0=one
      q1=one
      p0=x
      p1=x+one-alpha
      do 1 i=1,ilim
      ci=i
      if(p0.ne.zero.and.q0.ne.zero.and.q1.ne.zero)then
      r0=p0/q0
      r1=p1/q1
      if(dnrm(r0-r1).le.dnrm(r1)*error)then
      cdhs=r1
      return
      endif
c --------- Occasionally renormalize the sequences to avoid over(under)flow.
      if(dnrm(p0).gt.tol1.or.dnrm(p0).lt.tol2.or.
     * dnrm(q0).gt.tol1.or.dnrm(q0).lt.tol2)then
      factor=p0*q0
      p0=p0/factor
      q0=q0/factor
      p1=p1/factor
      q1=q1/factor
      endif
      endif
      p0=x*p1+ci*p0
      q0=x*q1+ci*q0
      p1=p0+(ci+one-alpha)*p1
      q1=q0+(ci+one-alpha)*q1
1     continue
c --- If the peripheral routines are written correctly,
c --- the following four statements should never be executed.
      cdhs=half*(r0+r1)
      return
      end
c
      double precision function dnrm(z)

      double complex z
      dnrm=dabs(dreal(z))+dabs(dimag(z))
      return
      end


      subroutine excld(
     +  x1,y1,z1,norm1,l1,m1,n1,a1,
     +  x2,y2,z2,norm2,l2,m2,n2,a2,
     +  nx,ny,nz,repul,range,repnrg)        

      double precision x1,y1,z1,norm1,a1,
     + x2,y2,z2,norm2,a2,
     + nx,ny,nz,repul,range,repnrg,ovrlp  

      integer l1,l2,l3,m1,m2,m3,n1,n2,n3

C     NOTE THAT THIS CODE ASSUME THAT ALL BASIS FUNCTIONS ARE CENTERED AT
C     THE ORIGIN!!

C     This allows us to assume that the gaussians from the two basis functions combine simply to 
C     give a new gaussian at the origin with exponent gamma

      double precision xa,ya,za,gamma

      xa = 0.0d+00
      ya = 0.0d+00
      za = 0.0d+00

      l3 = l1 + l2
      m3 = m1 + m2
      n3 = n1 + n2

      gamma = a1 + a2
c     write(6,*) norm1,norm2,repul,ovrlp

      call overlap(gamma,l3
     + ,m3,n3,xa,ya,za,
     +  range,0,0,0,            
     +  nx,ny,nz,ovrlp)

      repnrg = norm1*norm2*repul*ovrlp
      return
      end

